"""
Interval-censored estimation of MIC distributions, and the cut-off that
separates a wild-type population from everything above it.

A broth microdilution plate tests a doubling series of concentrations. An MIC
reported as 0.25 means growth was inhibited at 0.25 and not at 0.12, so the
true value lies in (0.12, 0.25]. Every reading is an interval. At the ends it
is worse: a reading of <=0.008 places the true MIC somewhere below the lowest
well and >1 places it above the highest. Treating those numbers as
measurements biases every summary of them, and where a group's median lands on
the plate floor the median reports the lowest tested concentration however far
below it the isolates' MICs lie.

This module models log2 MIC as normally distributed within a group and
estimates the parameters by maximum likelihood over the censoring intervals.
The contribution of one observation is

    Phi((upper - mu)/sigma) - Phi((lower - mu)/sigma)

with lower = -infinity for a left-censored reading and upper = +infinity for a
right-censored one, evaluated in log space throughout for the reason given in
log_interval_mass.

Nothing here reads a particular dataset. The inputs are reported MIC strings, a
tested concentration series, and optionally a cluster label per observation.
"""

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

__all__ = [
    "parse_concentration", "check_doubling_series", "mic_bounds",
    "mic_recorded", "place_mics", "log_interval_mass", "log_density_ratios",
    "fit_censored_normal", "fit_censored_linear", "shift",
    "resample_clusters", "resample_shift", "named_generator",
    "withhold_reason", "round_up_to_series", "ecoff", "simulate_reports",
]

# A labeled doubling series is rounded, so consecutive ratios read between
# about 1.87 and 2.09. A skipped dilution reads near 4.
RATIO_LOW = 1.7
RATIO_HIGH = 2.4

# A reported concentration is matched to a tested one within this relative
# tolerance, which absorbs the rounding in the labels without admitting a
# neighboring well half a dilution away.
MATCH_TOLERANCE = 0.05

# Box bounds on the fitted mean and on the log of the fitted standard
# deviation. A fit that ends on one of them found no interior maximum.
MU_BOUNDS = (-25.0, 25.0)
LOG_SIGMA_BOUNDS = (float(np.log(0.05)), float(np.log(20)))

# Distinct intervals of the series the placed readings must occupy before a
# normal fit is guaranteed a finite maximum. With one or two, the likelihood
# can rise without limit toward a zero standard deviation or an unbounded mean.
MINIMUM_INTERVALS = 3

# Doublings outside the tested series a fitted mean may lie before it is an
# extrapolation past every measurement.
EXTRAPOLATION_MARGIN = 1.0


def parse_concentration(value):
    """
    Read one concentration label as a number and its censoring operator.

    Plate layouts commonly store concentrations as text carrying the operator
    that marks the ends of the series, so "<=0.008" is the lowest tested
    concentration and ">1" the highest. The number is a tested concentration in
    both cases. Returns (concentration, operator or None), or (None, None) for
    anything that is not a concentration.
    """
    if value is None:
        return None, None
    text = str(value).strip()
    if text in ("", "nan", "None", "<NA>"):
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


def check_doubling_series(concentrations):
    """
    The consecutive ratios of a sorted concentration series, and whether every
    one of them is a doubling within the rounding the labels carry.

    Returns (ratios, regular). The censoring interval below a reported value is
    taken as the well beneath it, so a series with a skipped dilution describes
    intervals that are twice as wide as the code assumes and the caller has to
    know.
    """
    concentrations = sorted(float(c) for c in concentrations)
    if len(concentrations) < 2:
        return [], True
    ratios = (np.array(concentrations[1:]) / np.array(concentrations[:-1]))
    regular = bool(((ratios >= RATIO_LOW) & (ratios <= RATIO_HIGH)).all())
    return [float(r) for r in ratios], regular


def mic_recorded(value):
    """
    True when a reading was recorded, false when the plate produced none.

    A missing reading arrives as None, as a float nan, or as a library-specific
    missing value depending on what read the table, and all of them mean the
    same thing: no measurement, as against a measurement that does not match
    the tested series.
    """
    if value is None:
        return False
    try:
        if isinstance(value, float) and np.isnan(value):
            return False
    except TypeError:
        pass
    text = str(value).strip()
    if text in ("", "nan", "NaN", "None", "<NA>", "NaT", "NA"):
        return False
    return True


def mic_bounds(mic_text, concentrations):
    """
    Convert one reported MIC into the log2 interval it represents.

    Returns (lower, upper) in log2 units, with -inf or +inf at a censored end,
    or None where the reading is absent or does not match the tested series.
    An interior reading of c is the interval (the well below c, c]; a reading
    of <=c is everything at or below c; a reading of >c is everything above c.
    """
    if not mic_recorded(mic_text):
        return None
    text = str(mic_text).strip()

    if text.startswith("<="):
        try:
            return (-np.inf, float(np.log2(float(text[2:]))))
        except ValueError:
            return None
    if text.startswith(">="):
        try:
            return (float(np.log2(float(text[2:]))) - 1.0, np.inf)
        except ValueError:
            return None
    if text.startswith(">"):
        try:
            return (float(np.log2(float(text[1:]))), np.inf)
        except ValueError:
            return None

    try:
        value = float(text)
    except ValueError:
        return None
    if value <= 0:
        return None

    series = sorted(float(c) for c in concentrations)
    position = None
    for index, concentration in enumerate(series):
        if abs(concentration - value) / value < MATCH_TOLERANCE:
            position = index
            break
    if position is None:
        return None
    if position == 0:
        return (-np.inf, float(np.log2(value)))
    return (float(np.log2(series[position - 1])), float(np.log2(value)))


def place_mics(mic_texts, concentrations):
    """
    Convert a sequence of reported MICs into censoring intervals.

    Returns (lower, upper, absent, off_series). lower and upper are arrays
    carrying nan where no interval could be formed. absent counts the readings
    that were never recorded. off_series collects the readings that are present
    and do not match the tested series, which means the series does not
    describe the data; the caller has to act on that condition.
    """
    lower, upper, off_series, absent = [], [], [], 0
    for text in mic_texts:
        bounds = mic_bounds(text, concentrations)
        if bounds is not None:
            lower.append(bounds[0])
            upper.append(bounds[1])
            continue
        lower.append(np.nan)
        upper.append(np.nan)
        if mic_recorded(text):
            off_series.append(str(text))
        else:
            absent += 1
    return np.array(lower), np.array(upper), absent, off_series


def log_interval_mass(lower, upper):
    """
    Log of the standard normal probability mass on each interval.

    Written with logcdf, logsf and log1p. Far out in a tail two cdfs agree to
    most of their digits, so their difference loses the rest, which corrupts
    the gradient and stops an optimizer short of the maximum. Each interval is
    evaluated in the tail it lies in, so the only subtraction is of a number
    below one from one, inside log1p.
    """
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
    The standard normal density at each interval bound divided by the mass on
    that interval, as a pair of arrays.

    Both ratios are formed in log space for the same reason log_interval_mass
    is: where the mass is 1e-09 the density and the mass are both tiny and
    their quotient is of order one, so forming it directly loses the digits
    that carry it. An infinite bound contributes a density of zero.
    """
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

    Returns {"mu": float, "sigma": float} in log2 units, or None where the fit
    does not converge or no interval carries a finite bound.
    """
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    keep = ~np.isnan(lower) & ~np.isnan(upper)
    lower, upper = lower[keep], upper[keep]

    finite = np.concatenate([lower[np.isfinite(lower)], upper[np.isfinite(upper)]])
    if not len(finite):
        return None
    start = np.array([float(np.mean(finite)), np.log(max(float(np.std(finite)), 0.5))])

    def negative_log_likelihood(parameters):
        mu, log_sigma = parameters
        sigma = np.exp(log_sigma)
        return -np.sum(log_interval_mass((lower - mu) / sigma, (upper - mu) / sigma))

    result = minimize(negative_log_likelihood, start, method="L-BFGS-B",
                      bounds=[MU_BOUNDS, LOG_SIGMA_BOUNDS])
    if not result.success:
        return None
    return {"mu": float(result.x[0]), "sigma": float(np.exp(result.x[1]))}


def fit_censored_linear(lower, upper, design):
    """
    Maximum likelihood fit of a normal whose mean is a linear function of the
    columns of design, observed only through censoring intervals.

    The same likelihood as fit_censored_normal with the mean replaced by
    design @ beta, and with the gradient supplied in closed form. A group
    indicator beside a set of site indicators gives that group's effect with
    site held constant, which a single group mean cannot do when the group
    sits mostly at one site.

    Returns {"beta": array, "sigma": float}, or None where the fit does not
    converge. beta[0] is the intercept where the first column is ones.
    """
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
        # An infinite bound carries a density ratio of zero, and is replaced by
        # zero before the multiplication, so the product is never an infinity
        # against a zero.
        finite_low = np.where(np.isfinite(low), low, 0.0)
        finite_high = np.where(np.isfinite(high), high, 0.0)
        gradient = np.empty(len(parameters))
        gradient[:-1] = design.T @ ((at_high - at_low) / sigma)
        gradient[-1] = np.sum(finite_high * at_high - finite_low * at_low)
        return -np.sum(log_mass), gradient

    bounds = ([(-25, 25)] + [(-30, 30)] * (design.shape[1] - 1)
              + [(np.log(0.05), np.log(20))])
    # The default stopping rule halts about 2e-04 of log likelihood short of
    # the maximum, which is enough to move an effect in its third decimal.
    # With the gradient in closed form the tighter rule costs about 0.03 s.
    result = minimize(negative_log_likelihood, start, method="L-BFGS-B",
                      jac=True, bounds=bounds,
                      options={"ftol": 1e-14, "gtol": 1e-9, "maxiter": 2000})
    if not result.success:
        return None
    return {"beta": result.x[:-1], "sigma": float(np.exp(result.x[-1]))}


def shift(group_lower, group_upper, reference_lower, reference_upper):
    """
    The difference in fitted mean log2 MIC between a group and a reference
    group, in doublings. Returns None where either fit fails.

    A shift of 1.0 means the group's fitted mean MIC is twice the reference's.
    """
    a = fit_censored_normal(group_lower, group_upper)
    b = fit_censored_normal(reference_lower, reference_upper)
    if a is None or b is None:
        return None
    return a["mu"] - b["mu"]


def resample_clusters(lower, upper, clusters, rng, statistic=None, draws=400,
                      interval=(2.5, 97.5)):
    """
    Resample clusters with replacement and refit, returning the spread of the
    statistic over the draws.

    Isolates sharing a cluster are not independent observations: a variant
    carried by one clonal outbreak forty times is one event, and resampling
    isolates would treat it as forty. Clusters are drawn with replacement and
    every member of a drawn cluster enters the sample.

    statistic takes (lower, upper) and returns a float; it defaults to the
    fitted mean. Returns {"low", "high", "draws"}, or None where fewer than a
    quarter of the draws converged.
    """
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    clusters = np.asarray(clusters)
    if not (len(lower) == len(upper) == len(clusters)):
        raise ValueError("lower, upper and clusters must be the same length")
    if statistic is None:
        def statistic(low, high):
            fit = fit_censored_normal(low, high)
            return None if fit is None else fit["mu"]

    names = np.unique(clusters)
    members = {name: np.flatnonzero(clusters == name) for name in names}
    estimates = []
    for _ in range(draws):
        chosen = rng.choice(names, size=len(names), replace=True)
        index = np.concatenate([members[name] for name in chosen])
        value = statistic(lower[index], upper[index])
        if value is not None and np.isfinite(value):
            estimates.append(float(value))
    if len(estimates) < draws // 4:
        return None
    values = np.array(estimates)
    return {"low": float(np.percentile(values, interval[0])),
            "high": float(np.percentile(values, interval[1])),
            "draws": len(values)}


def resample_shift(group_lower, group_upper, group_clusters,
                   reference_lower, reference_upper, reference_clusters,
                   rng, draws=400, interval=(2.5, 97.5)):
    """
    The spread of the shift between a group and a reference group when
    clusters are resampled within each.

    Each draw resamples the group's clusters with replacement and, separately,
    the reference group's, refits both, and takes the difference of the two
    fitted means, so the interval carries the sampling error of both means.

    Returns {"low", "high", "draws"}, or None where fewer than a quarter of
    the draws converged in both groups.
    """
    arms = []
    for lower, upper, clusters in (
            (group_lower, group_upper, group_clusters),
            (reference_lower, reference_upper, reference_clusters)):
        lower = np.asarray(lower, dtype=float)
        upper = np.asarray(upper, dtype=float)
        clusters = np.asarray(clusters)
        if not (len(lower) == len(upper) == len(clusters)):
            raise ValueError("lower, upper and clusters must be the same length")
        if not len(clusters):
            raise ValueError("each group needs at least one observation")
        members = [np.flatnonzero(clusters == name) for name in np.unique(clusters)]
        arms.append((lower, upper, members))

    estimates = []
    for _ in range(draws):
        means = []
        for lower, upper, members in arms:
            chosen = rng.choice(len(members), size=len(members), replace=True)
            index = np.concatenate([members[i] for i in chosen])
            fit = fit_censored_normal(lower[index], upper[index])
            means.append(np.nan if fit is None else fit["mu"])
        difference = means[0] - means[1]
        if np.isfinite(difference):
            estimates.append(float(difference))
    if len(estimates) < draws // 4:
        return None
    values = np.array(estimates)
    return {"low": float(np.percentile(values, interval[0])),
            "high": float(np.percentile(values, interval[1])),
            "draws": len(values)}


def named_generator(seed, name):
    """
    A random generator seeded by a base seed and the name of the quantity it
    serves.

    Each interval draws from its own stream, so it does not depend on which
    other intervals were computed before it and can be reproduced on its own.
    The name enters as its sequence of code points.
    """
    name = str(name)
    if not name:
        raise ValueError("a stream needs a non-empty name")
    if int(seed) < 0:
        raise ValueError("the seed must be a non-negative integer")
    return np.random.default_rng([int(seed)] + [ord(letter) for letter in name])


def withhold_reason(fit, lower, upper, concentrations,
                    margin=EXTRAPOLATION_MARGIN):
    """
    Why a fitted normal does not support an estimate, or None when it does.

    Checked in this order: no reading was placed; the placed readings occupy
    fewer than MINIMUM_INTERVALS distinct intervals, where the likelihood can
    rise without limit toward a zero standard deviation or an unbounded mean;
    the fit did not converge; a parameter ended on the optimizer's box bound,
    so no interior maximum was found; the fitted mean lies more than `margin`
    doublings outside the tested series, so the estimate extrapolates past
    every measurement.
    """
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    keep = ~np.isnan(lower) & ~np.isnan(upper)
    if not keep.any():
        return "no reading was placed"
    occupied = len(set(zip(lower[keep].tolist(), upper[keep].tolist())))
    if occupied < MINIMUM_INTERVALS:
        return (f"the readings occupy too few intervals of the series "
                f"({occupied} of the {MINIMUM_INTERVALS} a finite maximum needs)")
    if fit is None:
        return "the fit did not converge"
    if (np.isclose(fit["mu"], MU_BOUNDS, rtol=0, atol=1e-6).any()
            or np.isclose(np.log(fit["sigma"]), LOG_SIGMA_BOUNDS,
                          rtol=0, atol=1e-6).any()):
        return "a parameter ended on the optimizer's bound"
    steps = np.log2(sorted(float(c) for c in concentrations))
    if fit["mu"] < steps[0] - margin or fit["mu"] > steps[-1] + margin:
        return (f"the fitted mean lies outside the tested series by more than "
                f"{margin:g} on the log2 scale")
    return None


def round_up_to_series(value, concentrations):
    """
    The lowest tested concentration at or above a log2 value, or None when the
    value exceeds the highest tested concentration.
    """
    for concentration in sorted(float(c) for c in concentrations):
        if np.log2(concentration) >= value - 1e-12:
            return concentration
    return None


def ecoff(fit, concentrations, coverage):
    """
    The cut-off separating a fitted wild-type distribution from everything
    above it, as a tested concentration.

    An epidemiological cut-off is the upper end of the wild-type distribution,
    which is log-normal in MIC, so the cut-off follows from the fit: take the
    quantile of the fitted normal at `coverage`, then round up to the next
    tested concentration, because a cut-off that falls between two wells cannot
    be applied to a plate that does not test it.

    `coverage` has no default and must be supplied. The normative document for
    the convention in a given setting is the EUCAST standard operating
    procedure for MIC distributions and cut-off setting; the percentage it
    prescribes was not retrievable from the sources reachable here, so this
    code does not choose one on the caller's behalf.

    Returns {"quantile": float in log2 units, "concentration": float or None,
    "above the series": bool}. concentration is None when the quantile exceeds
    the highest tested concentration, which means the plate cannot express the
    cut-off and a wider range is needed.
    """
    if fit is None:
        return None
    if not 0 < coverage < 1:
        raise ValueError("coverage must lie strictly between 0 and 1")
    quantile = float(norm.ppf(coverage, loc=fit["mu"], scale=fit["sigma"]))
    concentration = round_up_to_series(quantile, concentrations)
    return {"quantile": quantile, "concentration": concentration,
            "above the series": concentration is None}


def simulate_reports(mu, sd, concentrations, rng, draws=4000):
    """
    Draw log2 MICs from a known distribution and report them as a plate would.

    A plate reports the lowest tested concentration that inhibited growth, so a
    draw at or below the first well comes back censored at that well and one
    above the last comes back censored above it. Used to check that an
    estimator recovers parameters it was given.
    """
    series = sorted(float(c) for c in concentrations)
    steps = np.log2(np.asarray(series))
    reported = []
    for value in rng.normal(mu, sd, draws):
        index = int(np.searchsorted(steps, value))
        if index == 0:
            reported.append(f"<={series[0]}")
        elif index >= len(steps):
            reported.append(f">{series[-1]}")
        else:
            reported.append(f"{series[index]}")
    return reported
