"""
Estimates every effect on the discovery half only, and writes the predictions
down so they can be tested once on the held-out half.

Run from the project root with the virtual environment active:

    python code/discovery.py

Why this exists. Everything before this module was exploratory. Many
comparisons were made across two drugs, six variant classes, four lineages and
fourteen sites, with no correction for multiple testing, and the striking ones
were reported. That is how findings that fail to replicate are produced, and
three times already a headline dissolved on closer inspection. Those results are
hypotheses.

This module reads only the discovery half and never touches the held-out half
except to count how many samples it contains, which is not outcome data. It
writes docs/PRE_REGISTRATION.md containing numbered, quantitative predictions.
Commit that file before running code/validate.py. The commit timestamp is what
makes the ordering verifiable by someone who was not present, which is the
difference between claiming a pre-specified analysis and being able to show one.

The split. CRyPTIC-v3.0, which CRyPTIC identify as a validation set, holds 9,090
sequenced samples but only 10 with a UKMYC MIC and none carrying an Rv0678
variant, so it cannot test an MIC-level claim. The split used instead is
CRyPTIC-v1.0 against CRyPTIC-v2.0. v1.0 is the frozen pre-2020 collection handed
to WHO to build the first catalogue; v2.0 is everything added afterwards,
including roughly 1,100 NICD samples enriched for bedaquiline resistance. The
halves are therefore separated in time and partly in geography, which makes this
closer to external validation than to cross-validation. That is a stronger test,
and it also means a failure to replicate could reflect a population difference
rather than a false finding. Both must be stated in any write-up.

Cohort decisions, each following from code/audit_cohort.py rather than taste:

  High-quality phenotypes only. Every effect survived or strengthened when
  restricted to measurements where the independent reading methods agreed, so
  the weaker measurements are not carrying anything and excluding them is free.

  One isolate per patient. 5.0% of the cohort belongs to a repeated patient.
  The effects barely moved, but two isolates from one person are one
  observation and the analysis should say so.

  mmpL5 intact. Samples with a disrupted efflux pump are excluded from the
  variant groups, because the pump's absence changes what a repressor variant
  can do.

Outputs:
  outputs/discovery_report.txt
  docs/PRE_REGISTRATION.md
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import cohort  # noqa: E402
from cluster_adjust import add_clusters  # noqa: E402
from mic_model import concentration_series, fit_censored_normal, mic_bounds  # noqa: E402

REPORT = Path("outputs/discovery_report.txt")
PREREG = Path("docs/PRE_REGISTRATION.md")
DISCOVERY = "CRyPTIC-v1.0"
VALIDATION = "CRyPTIC-v2.0"
BOOTSTRAPS = 300
SEED = 20260101

_lines = []


def say(text=""):
    print(text)
    _lines.append(text)


def write_report():
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(_lines) + "\n")


def patient_of(identifier):
    import re

    text = str(identifier)
    site = text.split(".")[1] if "." in text else "?"
    match = re.search(r"\.subj\.(.*?)\.lab\.", text)
    return f"{site}|{match.group(1) if match else text}"


def load_cohort():
    """The full cohort with clusters, patients, MIC intervals and dataset labels."""
    import numpy as np
    import pandas as pd

    df = add_clusters(cohort.assemble(), cohort.load_mutations())
    df["PATIENT"] = [patient_of(i) for i in df.index]

    wgs = pd.read_parquet(cohort.DATA / "WGS_SAMPLES.parquet")
    if "UNIQUEID" not in wgs.columns:
        wgs = wgs.reset_index()
    dataset_column = next(c for c in wgs.columns if c.lower() == "dataset")
    df["DATASET"] = wgs.drop_duplicates("UNIQUEID").set_index("UNIQUEID")[dataset_column].reindex(df.index)

    series, _, _ = concentration_series()
    for drug in cohort.DRUGS:
        bounds = [
            mic_bounds(mic, series.get((design, drug)))
            if series.get((design, drug)) else None
            for mic, design in zip(df[f"MIC_{drug}"], df[f"PLATEDESIGN_{drug}"])
        ]
        df[f"lower_{drug}"] = [b[0] if b else np.nan for b in bounds]
        df[f"upper_{drug}"] = [b[1] if b else np.nan for b in bounds]
    return df


def restrict(df, drug, dataset):
    """Apply the cohort decisions: dataset, high quality, one isolate per patient."""
    subset = df[(df.DATASET == dataset) & (df[f"PHENOTYPE_QUALITY_{drug}"] == "HIGH")]
    return subset.sort_index().groupby("PATIENT", observed=True).head(1)


def binary_effect(frame, groups, drug, rng):
    """Odds ratio against reference, with an interval from resampling clusters."""
    import numpy as np
    from scipy.stats import fisher_exact

    exposed = frame[frame.GROUP.isin(groups)]
    reference = frame[frame.GROUP == "reference"]
    if len(exposed) < 5 or not len(reference):
        return None

    a, n1 = int(exposed[f"resistant_{drug}"].sum()), len(exposed)
    b, n2 = int(reference[f"resistant_{drug}"].sum()), len(reference)

    # With no resistant sample in either group there is no odds ratio to
    # estimate. A continuity correction would return a number, and that number
    # would describe the correction rather than the data.
    if a + b == 0:
        return {
            "exposed": n1, "resistant": 0, "clusters": exposed.CLUSTER.nunique(),
            "reference": n2, "reference_resistant": 0,
            "odds_ratio": None, "p": None, "ci_low": None, "ci_high": None,
            "note": "no resistant samples in either group",
        }

    point, p_value = fisher_exact([[a, n1 - a], [b, n2 - b]])

    exposed_groups = [g.index.to_numpy() for _, g in exposed.groupby("CLUSTER", observed=True)]
    exposed_outcome = exposed[f"resistant_{drug}"]
    reference_outcome = reference[f"resistant_{drug}"].to_numpy()

    estimates = []
    for _ in range(BOOTSTRAPS):
        chosen = rng.integers(0, len(exposed_groups), len(exposed_groups))
        indices = np.concatenate([exposed_groups[i] for i in chosen])
        ea = int(exposed_outcome.loc[indices].sum())
        en = len(indices)
        picked = rng.integers(0, n2, n2)
        rb = int(reference_outcome[picked].sum())
        # Haldane correction keeps a zero cell from producing an infinite odds ratio
        estimates.append(((ea + 0.5) * (n2 - rb + 0.5)) / ((en - ea + 0.5) * (rb + 0.5)))
    estimates = np.array(estimates)

    return {
        "exposed": n1, "resistant": a, "clusters": exposed.CLUSTER.nunique(),
        "reference": n2, "reference_resistant": b,
        "odds_ratio": float(point), "p": float(p_value),
        "ci_low": float(np.percentile(estimates, 2.5)),
        "ci_high": float(np.percentile(estimates, 97.5)),
        "draws": estimates,
    }


def mic_shift(frame, groups, drug, rng):
    """
    Shift in fitted mean log2 MIC against reference, in doublings.

    The interval is on the shift itself, resampling both groups within each
    draw, so it carries the reference group's uncertainty as well. Taking two
    separately fitted intervals and comparing them by eye does not.
    """
    import numpy as np

    exposed = frame[frame.GROUP.isin(groups)].dropna(subset=[f"lower_{drug}"])
    reference = frame[frame.GROUP == "reference"].dropna(subset=[f"lower_{drug}"])
    if len(exposed) < 8 or len(reference) < 50:
        return None

    exposed_fit = fit_censored_normal(exposed[f"lower_{drug}"], exposed[f"upper_{drug}"])
    reference_fit = fit_censored_normal(reference[f"lower_{drug}"], reference[f"upper_{drug}"])
    if not exposed_fit or not reference_fit:
        return None

    exposed_lower = exposed[f"lower_{drug}"].to_numpy()
    exposed_upper = exposed[f"upper_{drug}"].to_numpy()
    positions = {name: np.flatnonzero(exposed.CLUSTER.to_numpy() == name)
                 for name in exposed.CLUSTER.unique()}
    blocks = list(positions.values())
    reference_lower = reference[f"lower_{drug}"].to_numpy()
    reference_upper = reference[f"upper_{drug}"].to_numpy()

    shifts = []
    for _ in range(BOOTSTRAPS):
        chosen = rng.integers(0, len(blocks), len(blocks))
        index = np.concatenate([blocks[i] for i in chosen])
        exposed_draw = fit_censored_normal(exposed_lower[index], exposed_upper[index])
        picked = rng.integers(0, len(reference_lower), len(reference_lower))
        reference_draw = fit_censored_normal(reference_lower[picked], reference_upper[picked])
        if exposed_draw and reference_draw:
            shifts.append(exposed_draw["mu"] - reference_draw["mu"])
    if len(shifts) < BOOTSTRAPS // 4:
        return None
    shifts = np.array(shifts)

    return {
        "exposed": len(exposed), "clusters": exposed.CLUSTER.nunique(),
        "reference_mean": reference_fit["mu"], "exposed_mean": exposed_fit["mu"],
        "exposed_sd": exposed_fit["sigma"],
        "shift": exposed_fit["mu"] - reference_fit["mu"],
        "ci_low": float(np.percentile(shifts, 2.5)),
        "ci_high": float(np.percentile(shifts, 97.5)),
        "above_zero": float((shifts > 0).mean()),
        "draws": shifts,
    }


def censor(values, ladders, rng):
    """Round simulated log2 MICs onto the plate ladders, as a plate would."""
    import numpy as np

    lower, upper = [], []
    for value in values:
        ladder = ladders[rng.integers(0, len(ladders))]
        steps = np.log2(np.asarray(ladder))
        placed = False
        for index, step in enumerate(steps):
            if value <= step:
                lower.append(-np.inf if index == 0 else float(steps[index - 1]))
                upper.append(float(step))
                placed = True
                break
        if not placed:
            lower.append(float(steps[-1]))
            upper.append(np.inf)
    return np.array(lower), np.array(upper)


def simulate_shifts(exposed_n, reference_n, reference_mu, reference_sd,
                    exposed_sd, true_shift, ladders, rng, draws=200):
    """
    Shift estimates from simulated datasets of the held-out size.

    true_shift may be a single value or an array of candidate true values. A
    single zero measures the estimator's standard error under the null. Passing
    the discovery bootstrap draws makes the result a predictive distribution
    carrying both the discovery estimate's uncertainty and the held-out half's
    sampling error, rather than treating the discovery estimate as the truth.
    """
    import numpy as np

    candidates = np.atleast_1d(np.asarray(true_shift, dtype=float))
    estimates = []
    for _ in range(draws):
        drawn_shift = float(candidates[rng.integers(0, len(candidates))])
        exposed_values = rng.normal(reference_mu + drawn_shift, exposed_sd, exposed_n)
        reference_values = rng.normal(reference_mu, reference_sd, reference_n)
        exposed_lower, exposed_upper = censor(exposed_values, ladders, rng)
        reference_lower, reference_upper = censor(reference_values, ladders, rng)
        exposed_fit = fit_censored_normal(exposed_lower, exposed_upper)
        reference_fit = fit_censored_normal(reference_lower, reference_upper)
        if exposed_fit and reference_fit:
            estimates.append(exposed_fit["mu"] - reference_fit["mu"])
    return np.array(estimates)


def simulate_odds_ratios(exposed_n, reference_n, reference_rate, true_odds_ratio,
                         rng, draws=2000):
    """
    Odds ratios from simulated datasets of the held-out size.

    true_odds_ratio may be a single value or an array of candidate true values.
    Passing the discovery bootstrap draws makes the result a predictive
    distribution carrying both the discovery estimate's uncertainty and the
    held-out half's sampling error, rather than treating the discovery estimate
    as if it were the truth.
    """
    import numpy as np

    candidates = np.atleast_1d(np.asarray(true_odds_ratio, dtype=float))
    chosen = candidates[rng.integers(0, len(candidates), draws)]
    odds = reference_rate / (1 - reference_rate) * chosen
    rate = odds / (1 + odds)
    a = rng.binomial(exposed_n, rate)
    b = rng.binomial(reference_n, reference_rate, draws)
    # Haldane correction so a zero cell does not produce an infinite value
    return ((a + 0.5) * (reference_n - b + 0.5)) / ((exposed_n - a + 0.5) * (b + 0.5))


def predictive_interval(simulated, discovery_draws=None, validation_spread=None):
    """
    A 95% interval from simulated estimates.

    Empirical percentiles are unstable at a few hundred draws and biased
    inward, so the interval is built from the mean and standard deviation of
    the simulated estimates, which the simulations show to be close to normal.

    Where the simulated spread understates the combination of discovery and
    validation error, which happens because the estimator's spread depends on
    where the distribution sits relative to the plate, the spread is raised to
    the correct combined value rather than left short. An interval that is too
    narrow converts ordinary sampling variation into recorded failures.
    """
    import numpy as np

    values = np.asarray(simulated, dtype=float)
    centre = float(values.mean())
    spread = float(values.std(ddof=1))
    if discovery_draws is not None and validation_spread is not None:
        combined = float(np.hypot(np.asarray(discovery_draws).std(ddof=1),
                                  validation_spread))
        spread = max(spread, combined)
    return centre - 1.96 * spread, centre + 1.96 * spread


def detectable_odds_ratio(exposed_n, reference_n, reference_rate, rng, power=0.8):
    """
    The smallest odds ratio the held-out half could detect, found by simulation.

    Uses the discovery half's reference rate rather than the held-out half's, so
    no outcome data from the test set enters the calculation.
    """
    import numpy as np
    from scipy.stats import fisher_exact

    if exposed_n < 3 or reference_rate <= 0:
        return None
    for odds_ratio in np.arange(1.2, 60.0, 0.2):
        odds = reference_rate / (1 - reference_rate) * odds_ratio
        rate = odds / (1 + odds)
        significant = 0
        trials = 200
        for _ in range(trials):
            a = int(rng.binomial(exposed_n, rate))
            b = int(rng.binomial(reference_n, reference_rate))
            _, p = fisher_exact([[a, exposed_n - a], [b, reference_n - b]])
            significant += p < 0.05
        if significant / trials >= power:
            return float(odds_ratio)
    return None


def main():
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(SEED)

    say("=" * 72)
    say("Discovery half: estimation and pre-specification")
    say("=" * 72)

    df = load_cohort()
    plate_ladders, _, _ = concentration_series()
    lof_groups = [f"Rv0678 {c}" for c in cohort.LOF_CLASSES]
    comparisons = [
        ("Rv0678 loss of function", lof_groups),
        ("Rv0678 substitution", ["Rv0678 substitution"]),
        ("Rv0678 promoter", ["Rv0678 promoter"]),
        ("pepQ solo", ["pepQ solo"]),
    ]

    say(f"\nFull cohort: {len(df):,}")
    say(f"  {DISCOVERY}: {int((df.DATASET == DISCOVERY).sum()):,}")
    say(f"  {VALIDATION}: {int((df.DATASET == VALIDATION).sum()):,}")

    results, prereg_rows = {}, []
    for drug in cohort.DRUGS:
        discovery = restrict(df, drug, DISCOVERY)
        intact = discovery[(~discovery.mmpL5_LOF) | discovery.GROUP.eq("reference")]
        held_out = restrict(df, drug, VALIDATION)
        held_out = held_out[(~held_out.mmpL5_LOF) | held_out.GROUP.eq("reference")]

        say(f"\n{'=' * 72}")
        say(f"{drug}")
        say("=" * 72)
        say(f"\n  discovery cohort after filters: {len(intact):,}"
            f"   held-out cohort: {len(held_out):,}")

        reference = intact[intact.GROUP == "reference"]
        reference_rate = float(reference[f"resistant_{drug}"].mean())
        say(f"  reference group: {len(reference):,} samples, "
            f"{100 * reference_rate:.2f}% resistant")

        say("\n  Binary effects, interval from resampling clusters:")
        binary_rows = []
        for label, groups in comparisons:
            result = binary_effect(intact, groups, drug, rng)
            if not result:
                continue
            results[(drug, label, "binary")] = result
            undefined = result["odds_ratio"] is None
            binary_rows.append({
                "comparison": label,
                "n": result["exposed"], "clusters": result["clusters"],
                "resistant": result["resistant"],
                "odds ratio": "undefined" if undefined else round(result["odds_ratio"], 1),
                "95% CI": (result.get("note", "") if undefined
                           else f"{result['ci_low']:.1f} to {result['ci_high']:.1f}"),
                "p": "" if undefined else f"{result['p']:.2g}",
            })
        say(pd.DataFrame(binary_rows).to_string(index=False))

        say("\n  MIC shifts in doublings, interval on the shift itself:")
        shift_rows = []
        for label, groups in comparisons:
            result = mic_shift(intact, groups, drug, rng)
            if not result:
                continue
            results[(drug, label, "shift")] = result
            shift_rows.append({
                "comparison": label,
                "n": result["exposed"], "clusters": result["clusters"],
                "reference mean": round(result["reference_mean"], 2),
                "group mean": round(result["exposed_mean"], 2),
                "group sd": round(result["exposed_sd"], 2),
                "shift": round(result["shift"], 2),
                "95% CI": f"{result['ci_low']:.2f} to {result['ci_high']:.2f}",
                "fold": round(2 ** result["shift"], 1),
            })
        say(pd.DataFrame(shift_rows).to_string(index=False))

        say("\n  What the held-out half can detect, at 80% power:")
        say("  Binary power by simulation; MIC power from the estimator's standard")
        say("  error under the null, needing shift over standard error above 2.80.")

        ladders = [series for (design, other), series in plate_ladders.items()
                   if other == drug]
        held_reference = len(held_out[held_out.GROUP == "reference"])
        reference_fit = fit_censored_normal(
            reference[f"lower_{drug}"].dropna(),
            reference[f"upper_{drug}"].dropna(),
        )

        power_rows = []
        for label, groups in comparisons:
            exposed = held_out[held_out.GROUP.isin(groups)]
            exposed_n, exposed_clusters = len(exposed), exposed.CLUSTER.nunique()
            observed = results.get((drug, label, "binary"))
            shift_result = results.get((drug, label, "shift"))

            detectable_or = detectable_odds_ratio(
                exposed_n, held_reference, reference_rate, rng)
            binary_ok = bool(
                observed and observed["odds_ratio"] is not None and detectable_or
                and observed["ci_low"] >= detectable_or
            )

            detectable_shift, shift_ok = None, False
            if shift_result and reference_fit and ladders and exposed_clusters >= 3:
                null = simulate_shifts(
                    exposed_clusters, held_reference, reference_fit["mu"],
                    reference_fit["sigma"], shift_result["exposed_sd"], 0.0,
                    ladders, rng)
                if len(null) > 50:
                    detectable_shift = float(2.80 * null.std(ddof=1))
                    shift_ok = abs(shift_result["ci_low"]) >= detectable_shift and \
                        shift_result["ci_low"] * shift_result["shift"] > 0

            power_rows.append({
                "comparison": label,
                "held-out n": exposed_n,
                "held-out clusters": exposed_clusters,
                "min detectable OR": round(detectable_or, 1) if detectable_or else "none",
                "discovery OR low": (round(observed["ci_low"], 1)
                                     if observed and observed["odds_ratio"] is not None else ""),
                "binary testable": "yes" if binary_ok else "no",
                "min detectable shift": (round(detectable_shift, 2)
                                         if detectable_shift else "none"),
                "discovery shift low": (round(shift_result["ci_low"], 2)
                                        if shift_result else ""),
                "MIC testable": "yes" if shift_ok else "no",
            })

            if not (binary_ok or shift_ok):
                continue

            entry = {
                "drug": drug, "comparison": label,
                "held_out_n": exposed_n, "held_out_clusters": exposed_clusters,
                "binary_ok": binary_ok, "shift_ok": shift_ok,
                "odds_ratio": observed["odds_ratio"] if observed else None,
                "detectable_or": detectable_or,
                "shift": shift_result["shift"] if shift_result else None,
                "detectable_shift": detectable_shift,
            }
            if binary_ok:
                simulated = simulate_odds_ratios(
                    exposed_n, held_reference, reference_rate,
                    observed["draws"], rng, draws=4000)
                # odds ratios are skewed, so the interval is built on the log scale
                logged = np.log(simulated)
                centre, spread = float(logged.mean()), float(logged.std(ddof=1))
                entry["or_prediction"] = (float(np.exp(centre - 1.96 * spread)),
                                          float(np.exp(centre + 1.96 * spread)))
            if shift_ok:
                simulated = simulate_shifts(
                    exposed_clusters, held_reference, reference_fit["mu"],
                    reference_fit["sigma"], shift_result["exposed_sd"],
                    shift_result["draws"], ladders, rng, draws=400)
                entry["shift_prediction"] = predictive_interval(
                    simulated,
                    discovery_draws=shift_result["draws"],
                    validation_spread=detectable_shift / 2.80,
                )
            prereg_rows.append(entry)

        say(pd.DataFrame(power_rows).to_string(index=False))

    write_prereg(prereg_rows, df)
    write_report()
    say(f"\nReport written to {REPORT}")
    say(f"Pre-registration written to {PREREG}")
    say("\nCommit PRE_REGISTRATION.md and push it before running code/validate.py.")
    say("The commit timestamp is what makes the ordering verifiable.")


def write_prereg(rows, df):
    """Write the numbered predictions as a committable document."""
    lines = [
        "# Pre-registration: predictions to be tested on the held-out half",
        "",
        "Researcher: Amina Baktiyarova, Independent Researcher,",
        "ORCID 0009-0007-6265-6493",
        "",
        "Every estimate below comes from the CRyPTIC-v1.0 half of the CRyPTIC",
        "Consortium Dataset v3.4.0, Zenodo version DOI 10.5281/zenodo.15680920.",
        "",
        "## Prior exposure, disclosed",
        "",
        "This is a replication in a later subset with prior exposure, not a",
        "clean hold-out, and the difference is stated here rather than left for",
        "a reader to discover.",
        "",
        "Exploratory analyses run before this split was defined used the full",
        "cohort of 15,158 samples, which includes all 2,761 CRyPTIC-v2.0",
        "samples. Those analyses covered variant class effects, stratification",
        "by site and lineage, correction for clonal clustering, and interval-",
        "censored MIC estimation. A cohort audit additionally reported",
        "CRyPTIC-v2.0 reference resistance rates directly.",
        "",
        "What follows from that. The choice of which comparisons to carry",
        "forward was informed by results computed on combined data, so the",
        "selection is not independent of the held-out half. The quantitative",
        "predictions below were computed on CRyPTIC-v1.0 alone, which the code",
        "in code/discovery.py makes checkable, and the decision criteria and",
        "prediction intervals are fixed before code/validate.py is run.",
        "",
        "So this test can establish that an effect estimated in the earlier",
        "collection holds in a later and partly different one, under criteria",
        "fixed in advance. It cannot establish that the comparisons were chosen",
        "without any knowledge of the held-out samples. Any write-up must say",
        "so.",
        "",
        "## The split",
        "",
        "CRyPTIC-v3.0, which CRyPTIC designate as a validation set, holds 9,090",
        "sequenced samples but only 10 with a UKMYC MIC and none carrying an",
        "Rv0678 variant, so it cannot test a claim about MICs. The split used is",
        "CRyPTIC-v1.0 against CRyPTIC-v2.0. v1.0 is the frozen pre-2020",
        "collection handed to FIND and Seq&Treat to build the first WHO",
        "catalogue. v2.0 is everything added after that freeze, including",
        "approximately 1,100 samples from NICD enriched for bedaquiline",
        "resistance.",
        "",
        "The halves are separated in time and partly in geography, so this is",
        "closer to external validation than to cross-validation. That makes it a",
        "stronger test than a random split, and it also means a failure to",
        "replicate may reflect a difference between the populations rather than",
        "a false finding in the first. Any write-up must state both.",
        "",
        "This split was chosen on scientific grounds, before the power analysis",
        "showed how few comparisons it could test. It is not revised now that",
        "the answer is known, because choosing a split after seeing which one",
        "yields more testable predictions is the specific practice that",
        "pre-specification exists to prevent.",
        "",
        "## Cohort definition",
        "",
        "Applied identically to both halves:",
        "",
        "- Samples with both a UKMYC5 or UKMYC6 MIC and a genome.",
        "- Phenotype quality HIGH only, meaning at least two independent reading",
        "  methods agreed.",
        "- One isolate per patient, the first by identifier.",
        "- Samples with a disrupted mmpL5 excluded from the variant groups.",
        "- Variant groups defined from mutation strings alone. No WHO catalogue",
        "  content, and neither the EFFECTS nor the PREDICTIONS table, is used.",
        "",
        "## What counts as a successful prediction",
        "",
        "A held-out estimate carries its own sampling error, and the held-out",
        "groups are much smaller than the discovery groups. So a prediction is",
        "not that the held-out estimate lands inside the discovery confidence",
        "interval: it would miss far more often than 5% of the time even if the",
        "effect were exactly as estimated, and counting that as a failure to",
        "replicate would be wrong.",
        "",
        "Each prediction interval below is instead a prediction interval for the",
        "held-out estimate itself. It was built by simulating datasets of the",
        "held-out size, with the held-out cluster structure, censored onto the",
        "same plate dilution ladders, and re-estimated. The interval is the",
        "central 95% of those simulated estimates, taken as the mean plus and",
        "minus 1.96 standard deviations rather than as empirical percentiles,",
        "which are unstable and biased inward at a few hundred draws.",
        "",
        "The true effect used in each simulated dataset is itself drawn from the",
        "discovery estimate's bootstrap distribution rather than fixed at the",
        "discovery point estimate, so the interval carries both the uncertainty",
        "in the discovery estimate and the sampling error of the held-out half.",
        "Fixing it at the point estimate would make the interval too narrow and",
        "would turn ordinary sampling variation into recorded failures.",
        "",
        "If the effect is real and the same size in both halves, the held-out",
        "estimate should land inside the interval about 95% of the time.",
        "",
        "A prediction is supported if the held-out estimate falls inside the",
        "stated prediction interval and is in the stated direction. Failing",
        "either way is recorded as a failure.",
        "",
        "Comparisons that cannot be tested are reported as untested estimates",
        "with intervals, not dropped and not described as validated. The",
        "discovery report in outputs/discovery_report.txt lists every comparison",
        "and states which of them the held-out half is powered to test.",
        "",
        "## What is registered, and what is not",
        "",
        "Only comparisons the held-out half is powered to test are registered. A",
        "comparison qualifies when the lower bound of the discovery interval,",
        "not the point estimate, exceeds the smallest effect detectable at 80%",
        "power. Using the point estimate would register comparisons that could",
        "fail for lack of samples rather than lack of effect, leaving the result",
        "uninterpretable.",
        "",
        "Binary power was found by simulation. MIC power comes from the",
        "estimator's standard error under the null, measured by simulating",
        "datasets of the held-out size, with a shift over standard error of 2.80",
        "required for 80% power at a two-sided 5% level.",
        "",
        "The two tests are registered separately, because the MIC test uses the",
        "measurement itself and the binary test uses it only through a",
        "threshold, so one can be powered where the other is not.",
        "",
        "## Predictions",
        "",
    ]

    if not rows:
        lines += [
            "No comparison in the discovery half reaches an effect size the",
            "held-out half is powered to detect at 80%, under either test. No",
            "prediction is registered and the work remains exploratory.",
            "",
        ]
    else:
        for number, row in enumerate(rows, start=1):
            lines.append(f"### P{number}. {row['comparison']}, {row['drug']}")
            lines.append("")
            lines.append(
                f"Held-out samples: {row['held_out_n']}, in "
                f"{row['held_out_clusters']} clusters."
            )
            lines.append("")
            if row.get("binary_ok"):
                low, high = row["or_prediction"]
                lines.append(
                    f"**P{number}a, resistance.** The held-out odds ratio against the "
                    f"reference group will fall between {low:.1f} and {high:.1f}. "
                    f"Discovery estimate {row['odds_ratio']:.1f}; smallest odds ratio "
                    f"the held-out half can detect, {row['detectable_or']:.1f}."
                )
                lines.append("")
            if row.get("shift_ok"):
                low, high = row["shift_prediction"]
                direction = "exceed" if row["shift"] > 0 else "fall below"
                lines.append(
                    f"**P{number}b, MIC.** The held-out fitted mean log2 MIC will "
                    f"{direction} the reference group's, by between {low:.2f} and "
                    f"{high:.2f} doublings. Discovery estimate {row['shift']:.2f}; "
                    f"smallest shift the held-out half can detect, "
                    f"{row['detectable_shift']:.2f} doublings."
                )
                lines.append("")

    lines += [
        "## Analysis to be run",
        "",
        "code/validate.py, once, on CRyPTIC-v2.0. No estimate in this document",
        "will be revised afterwards. If a prediction fails, that is reported as a",
        "failure rather than explained away or replaced with a different",
        "comparison.",
        "",
        "## Data",
        "",
        "The CRyPTIC Consortium Dataset, version v3.4.0, Zenodo version DOI",
        "10.5281/zenodo.15680920, CC BY 4.0. See docs/PROVENANCE.md for the file",
        "inventory and checksums.",
        "",
    ]
    PREREG.parent.mkdir(parents=True, exist_ok=True)
    PREREG.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
