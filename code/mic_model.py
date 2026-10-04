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


def stream(name):
    """A generator seeded by name.

    Each estimate draws from its own stream, so its interval does not depend on
    which other estimates were computed before it and can be reproduced on its
    own.
    """
    import numpy as np

    return np.random.default_rng([SEED] + [ord(letter) for letter in name])


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


def mic_recorded(value):
    """
    True when a reading was recorded for this sample.

    A plate that produced no MIC stores a null, and the null arrives here as
    None, as a float nan, or as a pandas missing value depending on the pandas
    version reading the parquet file. A sample with no MIC is a missing
    measurement, not a value that failed to match the tested series.
    """
    import pandas as pd

    if value is None or pd.isna(value):
        return False
    return str(value).strip() not in ("", "nan", "None", "<NA>")


def place_mics(mic_texts, designs, series, drug):
    """
    Convert each reported MIC into its censoring interval.

    Returns (lowers, uppers, absent, off_series). A sample with no MIC recorded,
    or on a plate design with no tested series, is counted in absent and carries
    nan bounds. A value that is present and does not match the tested series is
    collected in off_series, which means the series does not describe the data.
    """
    import numpy as np

    lowers, uppers, off_series = [], [], []
    absent = 0
    for mic_text, design in zip(mic_texts, designs):
        concentrations = series.get((design, drug))
        bounds = mic_bounds(mic_text, concentrations) if concentrations else None
        if bounds is not None:
            lowers.append(bounds[0])
            uppers.append(bounds[1])
            continue
        lowers.append(np.nan)
        uppers.append(np.nan)
        if concentrations and mic_recorded(mic_text):
            off_series.append(str(mic_text))
        else:
            absent += 1
    return lowers, uppers, absent, off_series


def log_interval_mass(lower, upper):
    """
    Log of the standard normal probability mass on each interval [lower, upper].

    Written with logcdf, logsf and log1p rather than as the difference of two
    cdfs. Far out in a tail the two cdfs agree to most of their digits and
    subtracting them destroys the rest: the bedaquiline reference group holds
    intervals carrying a mass of 3e-09, where the difference loses enough
    precision to corrupt the gradient and stop the optimiser short of the
    maximum. Each interval is evaluated in the tail it lies in, so the only
    subtraction is of a number below one from one, inside log1p.
    """
    import numpy as np
    from scipy.stats import norm

    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    result = np.empty(len(lower))

    left = np.isneginf(lower)
    right = np.isposinf(upper)
    result[left] = norm.logcdf(upper[left])
    result[right] = norm.logsf(lower[right])

    middle = ~left & ~right
    low, high = lower[middle], upper[middle]
    mass = np.empty(len(low))
    in_upper_tail = low > 0
    tail_low = norm.logsf(low[in_upper_tail])
    tail_high = norm.logsf(high[in_upper_tail])
    mass[in_upper_tail] = tail_low + np.log1p(-np.exp(tail_high - tail_low))
    body_low = norm.logcdf(low[~in_upper_tail])
    body_high = norm.logcdf(high[~in_upper_tail])
    mass[~in_upper_tail] = body_high + np.log1p(-np.exp(body_low - body_high))
    result[middle] = mass
    return result


def log_density_ratios(lower, upper, log_mass):
    """
    The standard normal density at each interval bound, divided by the mass on
    that interval, returned as a pair of arrays.

    Both ratios are formed in log space for the same reason log_interval_mass
    is: where the mass is 1e-09 the density and the mass are both tiny, their
    quotient is of order one, and forming it directly loses the digits that
    carry it. An infinite bound contributes a density of zero.
    """
    import numpy as np
    from scipy.stats import norm

    ratios = []
    for bound in (np.asarray(lower, dtype=float), np.asarray(upper, dtype=float)):
        finite = np.isfinite(bound)
        ratio = np.zeros(len(bound))
        ratio[finite] = np.exp(norm.logpdf(bound[finite]) - log_mass[finite])
        ratios.append(ratio)
    return ratios[0], ratios[1]


def fit_censored_normal(lower, upper):
    """
    Maximum likelihood mean and standard deviation of a normal distribution
    observed only through censoring intervals.
    """
    import numpy as np
    from scipy.optimize import minimize

    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)

    finite = np.concatenate([lower[np.isfinite(lower)], upper[np.isfinite(upper)]])
    if not len(finite):
        return None
    start = np.array([float(np.mean(finite)), np.log(max(float(np.std(finite)), 0.5))])

    def negative_log_likelihood(parameters):
        mu, log_sigma = parameters
        sigma = np.exp(log_sigma)
        return -np.sum(log_interval_mass((lower - mu) / sigma, (upper - mu) / sigma))

    result = minimize(
        negative_log_likelihood, start, method="L-BFGS-B",
        bounds=[(-25, 25), (np.log(0.05), np.log(20))],
    )
    if not result.success:
        return None
    return {"mu": float(result.x[0]), "sigma": float(np.exp(result.x[1]))}


def fit_censored_linear(lower, upper, design):
    """
    Maximum likelihood fit of a normal whose mean is a linear function of the
    columns of `design`, observed only through censoring intervals.

    The same likelihood as fit_censored_normal, with the mean replaced by
    design @ beta, and with the gradient supplied in closed form. A group
    indicator beside a set of site indicators gives that group's shift with
    site held constant, which a single group mean cannot do when the group
    sits mostly at one site.

    Returns {"beta": array, "sigma": float}, or None where the fit does not
    converge. beta[0] is the intercept where the first column is ones.
    """
    import numpy as np
    from scipy.optimize import minimize

    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    design = np.asarray(design, dtype=float)
    if design.ndim != 2 or len(design) != len(lower):
        raise ValueError("design must be one row per observation")

    finite = np.concatenate([lower[np.isfinite(lower)], upper[np.isfinite(upper)]])
    if not len(finite):
        return None
    start = np.zeros(design.shape[1] + 1)
    start[0] = float(finite.mean())
    start[-1] = np.log(max(float(finite.std()), 0.5))

    def negative_log_likelihood(parameters):
        """The objective and its gradient, which share the interval mass."""
        mu = design @ parameters[:-1]
        sigma = np.exp(parameters[-1])
        low = (lower - mu) / sigma
        high = (upper - mu) / sigma
        log_mass = log_interval_mass(low, high)
        at_low, at_high = log_density_ratios(low, high, log_mass)
        # An infinite bound carries a density ratio of zero, and is replaced
        # by zero before the multiplication rather than after it, so that the
        # product is never an infinity against a zero.
        finite_low = np.where(np.isfinite(low), low, 0.0)
        finite_high = np.where(np.isfinite(high), high, 0.0)
        gradient = np.empty(len(parameters))
        gradient[:-1] = design.T @ ((at_high - at_low) / sigma)
        gradient[-1] = np.sum(finite_high * at_high - finite_low * at_low)
        return -np.sum(log_mass), gradient

    bounds = [(-25, 25)] + [(-30, 30)] * (design.shape[1] - 1) + [
        (np.log(0.05), np.log(20))]
    # The default stopping rule halts about 2e-04 of log likelihood short of
    # the maximum here, which is enough to move a shift in its third decimal.
    # With the gradient in closed form the tighter rule costs 0.03 s.
    result = minimize(negative_log_likelihood, start, method="L-BFGS-B",
                      jac=True, bounds=bounds,
                      options={"ftol": 1e-14, "gtol": 1e-9, "maxiter": 2000})
    if not result.success:
        return None
    return {"beta": result.x[:-1], "sigma": float(np.exp(result.x[-1]))}


PLANTED = [(-5.0, 1.0), (-6.5, 1.0), (-7.5, 1.2), (-2.0, 1.5), (-4.0, 0.8)]


def simulate_reports(mu, sd, concentrations, rng, draws=4000):
    """Draw log2 MICs from a known distribution and report them as a plate would.

    A plate reports the lowest tested concentration that inhibited growth, so a
    draw at or below the first well comes back censored at that well and one
    above the last comes back censored above it.
    """
    import numpy as np

    steps = np.log2(np.asarray(concentrations))
    reported = []
    for value in rng.normal(mu, sd, draws):
        index = int(np.searchsorted(steps, value))
        if index == 0:
            reported.append(f"<={concentrations[0]}")
        elif index >= len(steps):
            reported.append(f">{concentrations[-1]}")
        else:
            reported.append(f"{concentrations[index]}")
    return reported


def validate_estimator(concentrations, rng, draws=4000):
    """Recover each planted distribution from reports censored onto a real ladder.

    The naive median of the reported values is carried beside the fit, because
    the gap between the two is what the estimator exists to remove.
    """
    import numpy as np

    rows = []
    for mu, sd in PLANTED:
        reported = simulate_reports(mu, sd, concentrations, rng, draws)
        bounds = [mic_bounds(value, concentrations) for value in reported]
        fit = fit_censored_normal([bound[0] for bound in bounds],
                                  [bound[1] for bound in bounds])
        naive = float(np.median([np.log2(float(value.lstrip("<").lstrip("=").lstrip(">")))
                                 for value in reported]))
        left = sum(1 for value in reported if value.startswith("<="))
        rows.append({
            "true mean": mu,
            "fitted mean": round(fit["mu"], 2) if fit else None,
            "true sd": sd,
            "fitted sd": round(fit["sigma"], 2) if fit else None,
            "naive median": round(naive, 2),
            "left-censored %": round(100 * left / draws, 1),
        })
    return rows


SITE_MINIMUM = 100  # a fit on fewer placeable intervals is not a site effect


def excluded_sensitivity(frame, drug, minimum=SITE_MINIMUM):
    """How far the fitted reference mean could move if the excluded rows differed.

    A row carrying no MIC is excluded from every fit, so it can only act through
    the value it would have had. The movement is the excluded fraction times the
    gap between that value and the fitted mean, and the widest gap between a
    site's fitted mean and the whole group's is what a gap of that size looks
    like in this collection.
    """
    reference = frame[frame.GROUP.eq("reference")]
    fitted = reference.dropna(subset=[f"lower_{drug}"])
    whole = fit_censored_normal(fitted[f"lower_{drug}"], fitted[f"upper_{drug}"])
    if whole is None:
        return None
    widest, where = 0.0, None
    for site, chunk in fitted.groupby("SITEID", observed=True):
        if len(chunk) < minimum:
            continue
        fit = fit_censored_normal(chunk[f"lower_{drug}"], chunk[f"upper_{drug}"])
        if fit and abs(fit["mu"] - whole["mu"]) > widest:
            widest, where = abs(fit["mu"] - whole["mu"]), site
    excluded = len(reference) - len(fitted)
    share = excluded / len(reference)
    return {"drug": drug, "reference": len(reference), "fitted": len(fitted),
            "excluded": excluded, "fitted mean": round(whole["mu"], 3),
            "widest site gap": round(widest, 3), "at site": where,
            "movement at that gap": round(share * widest, 3),
            "movement at 3 doublings": round(share * 3.0, 3)}


def tiny_intervals(frame, drug, floor=1e-08):
    """Observed intervals whose probability mass under the fitted model is tiny.

    Below this floor the mass cannot be formed as a difference of two cumulative
    distribution functions without losing every digit of it, which is why the
    likelihood is evaluated in log space.
    """
    import numpy as np

    fitted = frame[frame.GROUP.eq("reference")].dropna(subset=[f"lower_{drug}"])
    fit = fit_censored_normal(fitted[f"lower_{drug}"], fitted[f"upper_{drug}"])
    if fit is None:
        return None
    low = (fitted[f"lower_{drug}"].to_numpy() - fit["mu"]) / fit["sigma"]
    high = (fitted[f"upper_{drug}"].to_numpy() - fit["mu"]) / fit["sigma"]
    mass = np.exp(log_interval_mass(low, high))
    return {"drug": drug, "intervals": len(mass), "floor": floor,
            "below the floor": int((mass < floor).sum()),
            "smallest": float(mass.min())}


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

    say("\n" + "=" * 72)
    say("Validation against planted parameters")
    say("=" * 72)
    say("")
    say("Each row simulates from a known distribution, censors the draws onto the")
    say("real UKMYC6 bedaquiline ladder, and refits. The naive median beside the")
    say("fit is the median of the reported values.")
    say("")
    say(pd.DataFrame(validate_estimator(series[("UKMYC6", "BDQ")],
                                        stream("planted parameters"))
                     ).to_string(index=False))

    df = add_clusters(cohort.assemble(), cohort.load_mutations())
    say(f"\nSamples: {len(df):,}")

    # ------------------------------------------------- build the intervals
    for drug in cohort.DRUGS:
        lowers, uppers, absent, off_series = place_mics(
            df[f"MIC_{drug}"], df[f"PLATEDESIGN_{drug}"], series, drug)
        df[f"lower_{drug}"] = lowers
        df[f"upper_{drug}"] = uppers
        placed = int(np.isfinite(df[f"lower_{drug}"]).sum()
                     + np.isinf(df[f"lower_{drug}"]).sum())
        if placed + absent + len(off_series) != len(df):
            die(f"{drug}: placed, absent and off-series counts do not sum to "
                f"{len(df):,}")
        say(f"\n{drug}: {placed:,} of {len(df):,} MICs placed on the tested series. "
            f"{absent:,} samples carry no MIC and are excluded from every fit.")
        if off_series:
            die(f"{len(off_series)} {drug} MICs carry a value that is not on the "
                f"tested series, so the series does not describe the data. "
                f"Examples: {sorted(set(off_series))[:10]}")

    # ------------------------------- what the excluded rows and the tails do
    say("\n" + "=" * 72)
    say("The excluded rows, and where the interval mass underflows")
    say("=" * 72)
    say("")
    say("An excluded row acts only through the value it would have had, so the")
    say("movement it can cause is its share of the group times the gap between")
    say("that value and the fitted mean. The widest gap between a site's fitted")
    say(f"mean and the whole group's is shown beside it, over sites carrying at")
    say(f"least {SITE_MINIMUM} placeable intervals.")
    say("")
    say(pd.DataFrame([excluded_sensitivity(df, drug) for drug in cohort.DRUGS]
                     ).to_string(index=False))
    say("")
    say(pd.DataFrame([tiny_intervals(df, drug) for drug in cohort.DRUGS]
                     ).to_string(index=False))

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

        rows, not_converged = [], []
        for group in cohort.GROUP_ORDER:
            if group == "reference":
                continue
            sub = df[(df.GROUP == group)].dropna(subset=[f"lower_{drug}"])
            if len(sub) < MIN_GROUP:
                continue
            fit = fit_censored_normal(sub[f"lower_{drug}"], sub[f"upper_{drug}"])
            if fit is None:
                not_converged.append(f"{group} (n = {len(sub):,})")
                continue
            interval = bootstrap_mu(sub, drug, stream(f"{drug} {group}"))
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
        if not_converged:
            say("  Absent from the table above because the fit did not converge: "
                + "; ".join(not_converged))

    # ------------------------------------------- by lineage, the key contrast
    say("\n" + "=" * 72)
    say("The same shift, within each lineage")
    say("=" * 72)

    lof_groups = [f"Rv0678 {c}" for c in cohort.LOF_CLASSES]
    intact = df[(~df.mmpL5_LOF) | df.GROUP.eq("reference")]

    for drug in cohort.DRUGS:
        say(f"\n{drug}:")
        rows, not_converged = [], []
        for lineage in ["lineage1", "lineage2", "lineage3", "lineage4"]:
            in_lineage = intact[intact.LINEAGE == lineage]
            reference = in_lineage[in_lineage.GROUP == "reference"].dropna(subset=[f"lower_{drug}"])
            if len(reference) < 50:
                continue
            reference_fit = fit_censored_normal(reference[f"lower_{drug}"], reference[f"upper_{drug}"])
            if reference_fit is None:
                not_converged.append(f"{lineage} reference (n = {len(reference):,})")
                continue
            for label, groups in [("loss of function", lof_groups),
                                  ("substitution", ["Rv0678 substitution"])]:
                sub = in_lineage[in_lineage.GROUP.isin(groups)].dropna(subset=[f"lower_{drug}"])
                if len(sub) < MIN_GROUP:
                    continue
                fit = fit_censored_normal(sub[f"lower_{drug}"], sub[f"upper_{drug}"])
                if fit is None:
                    not_converged.append(f"{lineage} {label} (n = {len(sub):,})")
                    continue
                interval = bootstrap_mu(
                    sub, drug, stream(f"{drug} {lineage} {label}"), draws=200)
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
        if not_converged:
            say("  Absent from the table above because the fit did not converge: "
                + "; ".join(not_converged))

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
