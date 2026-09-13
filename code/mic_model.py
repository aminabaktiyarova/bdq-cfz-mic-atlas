"""
Estimates MIC distributions properly, treating every measurement as the
interval it actually is.

Run from the project root with the virtual environment active:

    python code/mic_model.py

Why medians are not good enough. A broth microdilution plate tests a doubling
series of concentrations. When a plate reports an MIC of 0.25, it means growth
was inhibited at 0.25 and not at 0.12, so the true MIC lies somewhere in
(0.12, 0.25]. Every value on the plate is an interval, not a point. At the ends
it is worse: a reported <=0.008 means the true MIC is somewhere below the lowest
well and could be anything, and a reported >1 means it is above the highest.

Taking a median of those numbers treats interval-censored data as if it were
measured. Where a group's median lands on the plate floor, as it did for the
clofazimine reference group with 5,787 of 14,187 values left-censored, the
median is not an estimate of anything; it is a statement about the plate.

What this does instead. Within a group, log2 MIC is modelled as normally
distributed, and the parameters are estimated by maximum likelihood over the
censoring intervals:

    contribution of one observation = Phi((upper - mu)/sigma) - Phi((lower - mu)/sigma)

with lower = -infinity for a left-censored value and upper = +infinity for a
right-censored one. This is the approach the CRyPTIC ECOFF paper uses to define
wild-type distributions, reimplemented here as reusable code rather than
being reconstructed from a paper.

The estimate that matters is the shift: the difference in fitted mean log2 MIC
between a variant group and the reference group, in doublings. That is the
quantity a binary catalogue cannot express and the reason this resource exists.

Confidence intervals come from resampling clusters, not isolates, so a variant
carried by one clonal outbreak forty times contributes one draw rather than
forty. Everything learned in code/cluster_adjust.py applies here too.

Outputs:
  outputs/mic_estimates.csv     fitted distributions per group, per drug
  outputs/mic_model_report.txt  the readable summary
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import cohort  # noqa: E402
from cluster_adjust import add_clusters  # noqa: E402

REPORT = Path("outputs/mic_model_report.txt")
ESTIMATES = Path("outputs/mic_estimates.csv")
BOOTSTRAPS = 400
SEED = 20260101
MIN_GROUP = 8

_lines = []


def say(text=""):
    print(text)
    _lines.append(text)


def write_report():
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(_lines) + "\n")


def die(message):
    say(f"\nFAILED: {message}")
    write_report()
    sys.exit(1)


def parse_concentration(value):
    """
    Read one CONC entry as a number.

    PLATE_LAYOUT stores CONC as text rather than the float the schema document
    describes, and the entries carry censoring operators: the lowest well reads
    "<=0.008" and the off-plate bin above the highest well reads ">1". The
    number in both cases is a tested concentration, the lowest and the highest
    respectively, so stripping the operator and deduplicating recovers the
    tested ladder exactly. Returns (concentration, operator or None).
    """
    if value is None:
        return None, None
    text = str(value).strip()
    if text in ("", "nan", "None"):
        return None, None
    operator = None
    for candidate in ("<=", ">=", "<", ">"):
        if text.startswith(candidate):
            operator = candidate
            text = text[len(candidate):].strip()
            break
    try:
        return float(text), operator
    except ValueError:
        return None, None


def concentration_series():
    """
    The tested concentrations for each plate design and drug, from PLATE_LAYOUT.

    Returns {(plate design, drug): sorted list of concentrations}. The series is
    validated as a doubling ladder, because the censoring interval below each
    reported value is assumed to be half of it.
    """
    import numpy as np
    import pandas as pd

    layout = pd.read_parquet(cohort.DATA / "PLATE_LAYOUT.parquet")
    if "DRUG" not in layout.columns:
        layout = layout.reset_index()

    series, irregular, unreadable = {}, [], []
    for (design, drug), chunk in layout.groupby(["PLATEDESIGN", "DRUG"], observed=True):
        concentrations, operators = set(), set()
        for raw in chunk.CONC.dropna().unique():
            value, operator = parse_concentration(raw)
            if value is None or value <= 0:
                unreadable.append((design, drug, raw))
                continue
            concentrations.add(value)
            if operator:
                operators.add(operator)
        concentrations = sorted(concentrations)
        if len(concentrations) < 2:
            continue
        ratios = np.array(concentrations[1:]) / np.array(concentrations[:-1])
        # CRyPTIC label concentrations as rounded values, so a true doubling
        # series reads as ratios between about 1.87 and 2.09 rather than exactly
        # 2. A skipped dilution would give a ratio near 4, which is what this
        # check is for.
        if ((ratios < 1.7) | (ratios > 2.4)).any():
            irregular.append((design, drug, [round(float(r), 3) for r in ratios]))
        series[(design, drug)] = concentrations
    return series, irregular, unreadable


def mic_bounds(mic_text, concentrations):
    """
    Convert one reported MIC into the log2 interval it represents.

    Returns (lower, upper) in log2 units, with -inf or +inf at a censored end,
    or None if the value cannot be placed on the tested series.
    """
    import numpy as np

    text = str(mic_text).strip()
    if text in ("", "nan", "None"):
        return None

    if text.startswith("<="):
        value = float(text[2:])
        return (-np.inf, np.log2(value))
    if text.startswith(">="):
        value = float(text[2:])
        return (np.log2(value) - 1.0, np.inf)
    if text.startswith(">"):
        value = float(text[1:])
        return (np.log2(value), np.inf)

    try:
        value = float(text)
    except ValueError:
        return None

    position = None
    for index, concentration in enumerate(concentrations):
        if abs(concentration - value) / value < 0.05:
            position = index
            break
    if position is None:
        return None
    if position == 0:
        return (-np.inf, np.log2(value))
    return (np.log2(concentrations[position - 1]), np.log2(value))


def fit_censored_normal(lower, upper):
    """
    Maximum likelihood mean and standard deviation of a normal distribution
    observed only through censoring intervals.
    """
    import numpy as np
    from scipy.optimize import minimize
    from scipy.stats import norm

    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)

    finite = np.concatenate([lower[np.isfinite(lower)], upper[np.isfinite(upper)]])
    if not len(finite):
        return None
    start = np.array([float(np.mean(finite)), np.log(max(float(np.std(finite)), 0.5))])

    def negative_log_likelihood(parameters):
        mu, log_sigma = parameters
        sigma = np.exp(log_sigma)
        probability = norm.cdf((upper - mu) / sigma) - norm.cdf((lower - mu) / sigma)
        return -np.sum(np.log(np.clip(probability, 1e-12, None)))

    result = minimize(
        negative_log_likelihood, start, method="L-BFGS-B",
        bounds=[(-25, 25), (np.log(0.05), np.log(20))],
    )
    if not result.success:
        return None
    return {"mu": float(result.x[0]), "sigma": float(np.exp(result.x[1]))}


def bootstrap_mu(frame, drug, rng, draws=BOOTSTRAPS):
    """Resample clusters with replacement and refit, returning the spread of mu."""
    import numpy as np

    clusters = frame.CLUSTER.unique()
    by_cluster = {name: chunk for name, chunk in frame.groupby("CLUSTER", observed=True)}
    estimates = []
    for _ in range(draws):
        chosen = rng.choice(clusters, size=len(clusters), replace=True)
        import pandas as pd

        sample = pd.concat([by_cluster[name] for name in chosen])
        fit = fit_censored_normal(sample[f"lower_{drug}"], sample[f"upper_{drug}"])
        if fit:
            estimates.append(fit["mu"])
    if len(estimates) < draws // 4:
        return None
    values = np.array(estimates)
    return {
        "low": float(np.percentile(values, 2.5)),
        "high": float(np.percentile(values, 97.5)),
        "draws": len(values),
    }


def main():
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(SEED)

    say("=" * 72)
    say("Interval-censored MIC estimation")
    say("=" * 72)

    series, irregular, unreadable = concentration_series()
    say(f"\nPlate design and drug combinations in PLATE_LAYOUT: {len(series)}")
    if unreadable:
        say(f"\n  WARNING: {len(unreadable)} CONC entries could not be read as a")
        say("  concentration. Examples:")
        for design, drug, raw in unreadable[:8]:
            say(f"    {design} {drug}: {raw!r}")
    if irregular:
        say("\n  WARNING: a dilution appears to be missing from these series. The")
        say("  interval below a reported value is taken as the previous tested")
        say("  concentration, so a gap makes those bounds wrong:")
        for design, drug, ratios in irregular[:10]:
            say(f"    {design} {drug}: ratios {ratios}")
    for drug in cohort.DRUGS:
        for design in ("UKMYC5", "UKMYC6"):
            if (design, drug) in series:
                concentrations = series[(design, drug)]
                say(f"  {design} {drug}: {concentrations[0]} to {concentrations[-1]}, "
                    f"{len(concentrations)} concentrations")

    df = add_clusters(cohort.assemble(), cohort.load_mutations())
    say(f"\nSamples: {len(df):,}")

    # ------------------------------------------------- build the intervals
    for drug in cohort.DRUGS:
        lowers, uppers, unplaced = [], [], 0
        for mic_text, design in zip(df[f"MIC_{drug}"], df[f"PLATEDESIGN_{drug}"]):
            concentrations = series.get((design, drug))
            bounds = mic_bounds(mic_text, concentrations) if concentrations else None
            if bounds is None:
                unplaced += 1
                lowers.append(np.nan)
                uppers.append(np.nan)
            else:
                lowers.append(bounds[0])
                uppers.append(bounds[1])
        df[f"lower_{drug}"] = lowers
        df[f"upper_{drug}"] = uppers
        placed = int(np.isfinite(df[f"lower_{drug}"]).sum() + np.isinf(df[f"lower_{drug}"]).sum())
        say(f"\n{drug}: {placed:,} of {len(df):,} MICs placed on the tested series, "
            f"{unplaced:,} could not be placed")
        if unplaced > len(df) * 0.02:
            examples = df.loc[df[f"lower_{drug}"].isna(), f"MIC_{drug}"].dropna().unique()[:10]
            die(f"too many unplaceable {drug} MICs. Examples: {list(examples)}")

    # ------------------------------------------------------ fit each group
    say("\n" + "=" * 72)
    say("Fitted log2 MIC distributions, and the shift from reference")
    say("=" * 72)

    records = []
    for drug in cohort.DRUGS:
        reference = df[df.GROUP == "reference"].dropna(subset=[f"lower_{drug}"])
        reference_fit = fit_censored_normal(reference[f"lower_{drug}"], reference[f"upper_{drug}"])
        if reference_fit is None:
            die(f"the reference distribution for {drug} did not fit")

        say(f"\n{drug}")
        say(f"  reference: n = {len(reference):,}, fitted mean log2 MIC "
            f"{reference_fit['mu']:.2f}, standard deviation {reference_fit['sigma']:.2f}")
        say(f"  that mean corresponds to an MIC of {2**reference_fit['mu']:.4f}")
        naive = reference[f"LOG2MIC_{drug}"].median()
        say(f"  the naive median of the same data is {naive:.2f}, a difference of "
            f"{abs(reference_fit['mu'] - naive):.2f} doublings")

        rows = []
        for group in cohort.GROUP_ORDER:
            if group == "reference":
                continue
            sub = df[(df.GROUP == group)].dropna(subset=[f"lower_{drug}"])
            if len(sub) < MIN_GROUP:
                continue
            fit = fit_censored_normal(sub[f"lower_{drug}"], sub[f"upper_{drug}"])
            if fit is None:
                continue
            interval = bootstrap_mu(sub, drug, rng)
            shift = fit["mu"] - reference_fit["mu"]
            rows.append({
                "group": group,
                "isolates": len(sub),
                "clusters": sub.CLUSTER.nunique(),
                "fitted mean log2": round(fit["mu"], 2),
                "sd": round(fit["sigma"], 2),
                "shift (doublings)": round(shift, 2),
                "fold change": round(2 ** shift, 1),
                "95% CI on mean": (f"{interval['low']:.2f} to {interval['high']:.2f}"
                                   if interval else "did not converge"),
            })
            records.append({"drug": drug, "group": group, **rows[-1]})
        say("")
        say(pd.DataFrame(rows).to_string(index=False))

    # ------------------------------------------- by lineage, the key contrast
    say("\n" + "=" * 72)
    say("The same shift, within each lineage")
    say("=" * 72)

    lof_groups = [f"Rv0678 {c}" for c in cohort.LOF_CLASSES]
    intact = df[(~df.mmpL5_LOF) | df.GROUP.eq("reference")]

    for drug in cohort.DRUGS:
        say(f"\n{drug}:")
        rows = []
        for lineage in ["lineage1", "lineage2", "lineage3", "lineage4"]:
            in_lineage = intact[intact.LINEAGE == lineage]
            reference = in_lineage[in_lineage.GROUP == "reference"].dropna(subset=[f"lower_{drug}"])
            if len(reference) < 50:
                continue
            reference_fit = fit_censored_normal(reference[f"lower_{drug}"], reference[f"upper_{drug}"])
            if reference_fit is None:
                continue
            for label, groups in [("loss of function", lof_groups),
                                  ("substitution", ["Rv0678 substitution"])]:
                sub = in_lineage[in_lineage.GROUP.isin(groups)].dropna(subset=[f"lower_{drug}"])
                if len(sub) < MIN_GROUP:
                    continue
                fit = fit_censored_normal(sub[f"lower_{drug}"], sub[f"upper_{drug}"])
                if fit is None:
                    continue
                interval = bootstrap_mu(sub, drug, rng, draws=200)
                shift = fit["mu"] - reference_fit["mu"]
                rows.append({
                    "lineage": lineage,
                    "group": label,
                    "isolates": len(sub),
                    "clusters": sub.CLUSTER.nunique(),
                    "reference mean": round(reference_fit["mu"], 2),
                    "group mean": round(fit["mu"], 2),
                    "shift": round(shift, 2),
                    "fold": round(2 ** shift, 1),
                    "95% CI on group mean": (f"{interval['low']:.2f} to {interval['high']:.2f}"
                                             if interval else "n/a"),
                })
        say(pd.DataFrame(rows).to_string(index=False))

    say("\n  A shift is in doublings of MIC. A shift of 1.0 means the fitted mean MIC")
    say("  is twice the reference. Intervals come from resampling clusters, so a")
    say("  variant carried by one outbreak is not counted as many observations.")

    ESTIMATES.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(ESTIMATES, index=False)
    say(f"\nEstimates written to {ESTIMATES}")
    write_report()
    say(f"Report written to {REPORT}")


if __name__ == "__main__":
    main()
