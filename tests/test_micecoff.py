"""
Tests for the micecoff package, which carries the interval-censored estimator
as reusable code.

The package reads no dataset, so every test here builds its own data with the
answer planted in it. Where a test checks an estimator, it checks that the
estimator recovers a planted value.
"""

import numpy as np
import pytest
from scipy.stats import norm

from micecoff import core

# The UKMYC6 bedaquiline ladder, used because its labels are rounded, as a
# real plate's are, so consecutive ratios are not exactly 2.
LADDER = [0.008, 0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1.0]


def objective_and_gradient(lower, upper, design, parameters):
    """
    The module's own objective and closed-form gradient, re-formed here so a
    test can evaluate them at a chosen point. The formula is checked against a
    difference quotient below, which is what licenses the stationarity test
    that uses it.
    """
    mu = design @ parameters[:-1]
    sigma = np.exp(parameters[-1])
    low = (lower - mu) / sigma
    high = (upper - mu) / sigma
    log_mass = core.log_interval_mass(low, high)
    at_low, at_high = core.log_density_ratios(low, high, log_mass)
    finite_low = np.where(np.isfinite(low), low, 0.0)
    finite_high = np.where(np.isfinite(high), high, 0.0)
    gradient = np.empty(len(parameters))
    gradient[:-1] = design.T @ ((at_high - at_low) / sigma)
    gradient[-1] = np.sum(finite_high * at_high - finite_low * at_low)
    return -np.sum(log_mass), gradient


def seeded(offset=0):
    return np.random.default_rng(20260101 + offset)


# parse_concentration


@pytest.mark.parametrize("text,value,operator", [
    ("0.25", 0.25, None),
    ("<=0.008", 0.008, "<="),
    (">1", 1.0, ">"),
    (">=2", 2.0, ">="),
    ("  0.5  ", 0.5, None),
])
def test_a_concentration_label_yields_its_number_and_operator(text, value, operator):
    assert core.parse_concentration(text) == (value, operator)


@pytest.mark.parametrize("text", [None, "", "nan", "None", "<NA>", "POS", "n/a"])
def test_a_label_that_is_not_a_concentration_yields_nothing(text):
    assert core.parse_concentration(text) == (None, None)


# check_doubling_series


def test_a_rounded_doubling_series_is_regular():
    ratios, regular = core.check_doubling_series(LADDER)
    assert regular
    assert len(ratios) == len(LADDER) - 1
    assert not all(abs(r - 2.0) < 1e-9 for r in ratios), (
        "the ladder under test must be a rounded one, or this proves nothing")


def test_a_skipped_dilution_is_not_a_doubling_series():
    _, regular = core.check_doubling_series([0.008, 0.015, 0.06, 0.12])
    assert not regular


def test_a_series_of_one_concentration_has_no_ratios():
    assert core.check_doubling_series([0.25]) == ([], True)


# mic_recorded


@pytest.mark.parametrize("value", [None, float("nan"), "", "nan", "NaN",
                                   "None", "<NA>", "NaT", "NA"])
def test_an_absent_reading_is_not_recorded(value):
    assert not core.mic_recorded(value)


@pytest.mark.parametrize("value", ["0.25", "<=0.008", ">1", 0.25])
def test_a_present_reading_is_recorded(value):
    assert core.mic_recorded(value)


# mic_bounds


def test_an_interior_reading_is_the_interval_below_it():
    assert core.mic_bounds("0.25", LADDER) == (np.log2(0.12), np.log2(0.25))


def test_a_reading_at_the_lowest_well_is_left_censored():
    lower, upper = core.mic_bounds("0.008", LADDER)
    assert lower == -np.inf and upper == np.log2(0.008)


def test_an_explicitly_left_censored_reading_is_left_censored():
    lower, upper = core.mic_bounds("<=0.008", LADDER)
    assert lower == -np.inf and upper == np.log2(0.008)


def test_a_right_censored_reading_is_everything_above_the_highest_well():
    lower, upper = core.mic_bounds(">1", LADDER)
    assert lower == np.log2(1.0) and upper == np.inf


def test_an_at_or_above_reading_opens_at_the_well_below():
    lower, upper = core.mic_bounds(">=1", LADDER)
    assert lower == pytest.approx(np.log2(1.0) - 1.0) and upper == np.inf


def test_a_reading_off_the_tested_series_is_refused():
    assert core.mic_bounds("0.2", LADDER) is None


def test_a_reading_within_the_label_rounding_is_matched():
    """A series labeled 0.12 against a reading of 0.125 is one well."""
    assert core.mic_bounds("0.125", LADDER) is not None


def test_an_absent_reading_yields_no_interval():
    assert core.mic_bounds(None, LADDER) is None
    assert core.mic_bounds(float("nan"), LADDER) is None


# place_mics


def test_absent_readings_and_off_series_readings_are_counted_apart():
    lower, upper, absent, off = core.place_mics(
        ["0.25", None, "0.2", "<=0.008", ""], LADDER)
    assert absent == 2
    assert off == ["0.2"]
    assert np.isnan(lower[1]) and np.isnan(lower[2])
    assert lower[0] == np.log2(0.12)
    assert lower[3] == -np.inf


# log_interval_mass


def test_the_log_space_mass_matches_a_direct_difference_through_the_body():
    low = np.array([-2.0, -1.0, 0.0, 0.5, 1.5])
    high = low + 1.0
    direct = norm.cdf(high) - norm.cdf(low)
    assert np.allclose(np.exp(core.log_interval_mass(low, high)), direct, rtol=1e-12)


def test_the_mass_stays_finite_where_a_direct_difference_underflows():
    """At 40 standard deviations the two cdfs are equal in double precision."""
    low = np.array([40.0])
    high = np.array([41.0])
    assert norm.cdf(high[0]) - norm.cdf(low[0]) == 0.0
    mass = core.log_interval_mass(low, high)
    assert np.isfinite(mass[0]) and mass[0] < -300


def test_a_censored_bound_uses_the_tail_directly():
    assert core.log_interval_mass([-np.inf], [1.5])[0] == pytest.approx(
        norm.logcdf(1.5))
    assert core.log_interval_mass([1.5], [np.inf])[0] == pytest.approx(
        norm.logsf(1.5))


# fit_censored_normal


@pytest.mark.parametrize("mu,sd", [(-5.0, 1.0), (-6.5, 1.0), (-2.0, 1.5)])
def test_the_estimator_recovers_planted_parameters(mu, sd):
    reported = core.simulate_reports(mu, sd, LADDER, seeded(), draws=4000)
    lower, upper, _, _ = core.place_mics(reported, LADDER)
    fit = core.fit_censored_normal(lower, upper)
    assert fit["mu"] == pytest.approx(mu, abs=0.08)
    assert fit["sigma"] == pytest.approx(sd, abs=0.08)


def test_the_estimator_beats_the_median_under_heavy_left_censoring():
    """At the plate floor the median reports the well, not the distribution."""
    reported = core.simulate_reports(-7.5, 1.2, LADDER, seeded(1), draws=4000)
    lower, upper, _, _ = core.place_mics(reported, LADDER)
    fit = core.fit_censored_normal(lower, upper)
    censored = sum(1 for r in reported if r.startswith("<="))
    naive = float(np.median([np.log2(float(r.lstrip("<=>"))) for r in reported]))
    assert censored / len(reported) > 0.5
    assert abs(fit["mu"] - (-7.5)) < abs(naive - (-7.5))


def test_a_fit_with_no_finite_bound_returns_nothing():
    assert core.fit_censored_normal([np.nan], [np.nan]) is None


# fit_censored_linear


def test_the_linear_fit_separates_a_group_effect_from_a_site_effect():
    """
    Carriers concentrated at a site whose baseline is high. The group mean
    confounds the two; the design with a site indicator does not.
    """
    rng = seeded(2)
    site_offset, group_effect = 1.5, 1.0
    rows = []
    for _ in range(600):                      # site A, no carriers
        rows.append((0, 0, rng.normal(-5.0, 1.0)))
    for _ in range(540):                      # site B, baseline raised
        rows.append((1, 0, rng.normal(-5.0 + site_offset, 1.0)))
    for _ in range(60):                       # carriers, almost all at site B
        rows.append((1, 1, rng.normal(-5.0 + site_offset + group_effect, 1.0)))
    site = np.array([r[0] for r in rows], dtype=float)
    group = np.array([r[1] for r in rows], dtype=float)
    latent = np.array([r[2] for r in rows])

    reported = []
    steps = np.log2(np.asarray(LADDER))
    for value in latent:
        index = int(np.searchsorted(steps, value))
        reported.append(f"<={LADDER[0]}" if index == 0 else
                        f">{LADDER[-1]}" if index >= len(steps) else f"{LADDER[index]}")
    lower, upper, _, _ = core.place_mics(reported, LADDER)

    plain = core.fit_censored_linear(lower, upper,
                                     np.column_stack([np.ones(len(rows)), group]))
    adjusted = core.fit_censored_linear(
        lower, upper, np.column_stack([np.ones(len(rows)), group, site]))
    assert plain["beta"][1] > group_effect + 0.4, "the confounded fit must be inflated"
    assert adjusted["beta"][1] == pytest.approx(group_effect, abs=0.35)
    assert adjusted["beta"][2] == pytest.approx(site_offset, abs=0.25)


def test_the_closed_form_gradient_matches_a_difference_quotient():
    rng = seeded(3)
    reported = core.simulate_reports(-5.0, 1.0, LADDER, rng, draws=500)
    lower, upper, _, _ = core.place_mics(reported, LADDER)
    design = np.column_stack([np.ones(len(lower)), rng.normal(size=len(lower))])

    point = np.array([-4.8, 0.3, np.log(1.1)])
    _, analytic = objective_and_gradient(lower, upper, design, point)
    step = 1e-6
    for index in range(len(point)):
        shifted = point.copy()
        shifted[index] += step
        forward, _ = objective_and_gradient(lower, upper, design, shifted)
        shifted[index] -= 2 * step
        backward, _ = objective_and_gradient(lower, upper, design, shifted)
        numeric = (forward - backward) / (2 * step)
        assert analytic[index] == pytest.approx(numeric, rel=1e-4)


def test_the_linear_fit_returns_a_stationary_point():
    """
    The optimizer's default stopping rule halts short of the maximum, which
    moves an effect in its third decimal, so the module asks for a tighter
    one. At the returned solution the gradient must be small. Measured on this
    data the tight rule leaves a largest component near 6e-06 and the default
    near 2e-02, so the threshold below sits orders of magnitude from both.
    """
    rng = seeded(11)
    reported = core.simulate_reports(-5.0, 1.0, LADDER, rng, draws=1200)
    lower, upper, _, _ = core.place_mics(reported, LADDER)
    group = (rng.random(len(lower)) < 0.2).astype(float)
    site = (rng.random(len(lower)) < 0.5).astype(float)
    design = np.column_stack([np.ones(len(lower)), group, site])

    fit = core.fit_censored_linear(lower, upper, design)
    parameters = np.append(fit["beta"], np.log(fit["sigma"]))
    _, gradient = objective_and_gradient(lower, upper, design, parameters)
    assert np.max(np.abs(gradient)) < 1e-3


def test_the_fit_ignores_rows_carrying_no_interval():
    """
    place_mics leaves nan bounds where no reading was recorded. Those rows are
    missing measurements and must not reach the likelihood, nor turn the fit
    into a nan.
    """
    reported = core.simulate_reports(-5.0, 1.0, LADDER, seeded(9), draws=1500)
    clean_lower, clean_upper, _, _ = core.place_mics(reported, LADDER)
    padded_lower, padded_upper, absent, _ = core.place_mics(
        list(reported) + [None] * 300, LADDER)
    assert absent == 300
    clean = core.fit_censored_normal(clean_lower, clean_upper)
    padded = core.fit_censored_normal(padded_lower, padded_upper)
    assert np.isfinite(padded["mu"]) and np.isfinite(padded["sigma"])
    assert padded["mu"] == pytest.approx(clean["mu"], abs=1e-9)
    assert padded["sigma"] == pytest.approx(clean["sigma"], abs=1e-9)


# shift


def test_the_shift_recovers_a_planted_difference():
    reference = core.simulate_reports(-5.0, 1.0, LADDER, seeded(4), draws=3000)
    group = core.simulate_reports(-5.0 + 1.8, 1.0, LADDER, seeded(5), draws=3000)
    ref_lower, ref_upper, _, _ = core.place_mics(reference, LADDER)
    grp_lower, grp_upper, _, _ = core.place_mics(group, LADDER)
    assert core.shift(grp_lower, grp_upper, ref_lower, ref_upper) == pytest.approx(
        1.8, abs=0.1)


# resample_clusters


def test_clustered_data_gives_a_wider_interval_than_the_same_values_alone():
    """
    One dataset, two cluster labellings. Forty isolates in five clonal groups
    carry five independent observations; the same forty relabeled carry forty.
    """
    rng = seeded(6)
    values, clusters = [], []
    for group in range(5):
        center = rng.normal(-5.0, 1.0)
        for _ in range(8):
            values.append(center + rng.normal(0, 0.05))
            clusters.append(f"clone{group}")
    steps = np.log2(np.asarray(LADDER))
    reported = [f"<={LADDER[0]}" if int(np.searchsorted(steps, v)) == 0 else
                f">{LADDER[-1]}" if int(np.searchsorted(steps, v)) >= len(steps) else
                f"{LADDER[int(np.searchsorted(steps, v))]}" for v in values]
    lower, upper, _, _ = core.place_mics(reported, LADDER)

    clustered = core.resample_clusters(lower, upper, clusters,
                                       seeded(7), draws=200)
    independent = core.resample_clusters(lower, upper, [str(i) for i in range(40)],
                                         seeded(7), draws=200)
    assert (clustered["high"] - clustered["low"]) > 1.5 * (
        independent["high"] - independent["low"])


def test_resampling_refuses_mismatched_lengths():
    with pytest.raises(ValueError):
        core.resample_clusters([0.0], [1.0], ["a", "b"], seeded())


# ecoff


def test_the_cut_off_recovers_a_planted_wild_type_boundary():
    """
    A wild-type population whose 99th percentile falls between two wells must
    give the well above it, because a plate cannot express a cut-off it does
    not test.
    """
    fit = {"mu": np.log2(0.015), "sigma": 0.5}
    expected = norm.ppf(0.99, loc=fit["mu"], scale=fit["sigma"])
    assert np.log2(0.03) < expected <= np.log2(0.06), (
        "the planted quantile must fall between two wells")
    result = core.ecoff(fit, LADDER, 0.99)
    assert result["quantile"] == pytest.approx(expected)
    assert result["concentration"] == 0.06
    assert not result["above the series"]


def test_a_cut_off_above_the_tested_range_reports_no_concentration():
    result = core.ecoff({"mu": np.log2(1.0), "sigma": 1.0}, LADDER, 0.99)
    assert result["concentration"] is None
    assert result["above the series"]


def test_wider_coverage_never_lowers_the_cut_off():
    fit = {"mu": np.log2(0.015), "sigma": 0.8}
    quantiles = [core.ecoff(fit, LADDER, c)["quantile"]
                 for c in (0.90, 0.95, 0.99, 0.999)]
    assert quantiles == sorted(quantiles)


@pytest.mark.parametrize("coverage", [0.0, 1.0, -0.1, 1.5])
def test_a_coverage_outside_the_unit_interval_is_refused(coverage):
    with pytest.raises(ValueError):
        core.ecoff({"mu": -5.0, "sigma": 1.0}, LADDER, coverage)


def test_no_fit_gives_no_cut_off():
    assert core.ecoff(None, LADDER, 0.99) is None


# simulate_reports


def test_the_simulator_censors_at_both_ends_of_the_plate():
    reported = core.simulate_reports(-5.0, 4.0, LADDER, seeded(8), draws=2000)
    assert any(r.startswith("<=") for r in reported)
    assert any(r.startswith(">") for r in reported)
    assert set(reported) <= ({f"<={LADDER[0]}", f">{LADDER[-1]}"}
                             | {str(c) for c in LADDER[1:]})


# resample_shift


def placed(mu, sd, n, rng):
    lower, upper, _, _ = core.place_mics(
        core.simulate_reports(mu, sd, LADDER, rng, draws=n), LADDER)
    return lower, upper


def test_the_two_sample_interval_brackets_the_shift_at_the_predicted_width():
    """
    Groups of 300 with a planted shift of 1.5. A normal approximation puts the
    interval's width near 3.92 * sqrt(2 / 300), about 0.32; the check allows
    a factor of two either way. Coverage of the truth across datasets was
    measured by simulation.
    """
    rng = seeded(12)
    group = placed(-3.5, 1.0, 300, rng)
    reference = placed(-5.0, 1.0, 300, rng)
    result = core.resample_shift(group[0], group[1], np.arange(300),
                                 reference[0], reference[1], np.arange(300),
                                 seeded(13), draws=200)
    point = core.shift(group[0], group[1], reference[0], reference[1])
    assert point == pytest.approx(1.5, abs=0.15)
    assert result["low"] < point < result["high"]
    predicted = 3.92 * (2 / 300) ** 0.5
    assert 0.5 * predicted < result["high"] - result["low"] < 2.0 * predicted
    assert result["draws"] == 200


def test_the_two_sample_interval_carries_the_reference_error():
    """
    Equal groups of sixty. Resampling the reference group as well as the group
    roughly doubles the variance of the shift, so the interval is near 1.4
    times the width of the group's own interval, which holds the reference
    mean fixed.
    """
    rng = seeded(14)
    group = placed(-4.0, 1.0, 60, rng)
    reference = placed(-5.0, 1.0, 60, rng)
    both = core.resample_shift(group[0], group[1], np.arange(60),
                               reference[0], reference[1], np.arange(60),
                               seeded(15), draws=300)
    alone = core.resample_clusters(group[0], group[1], np.arange(60),
                                   seeded(15), draws=300)
    ratio = (both["high"] - both["low"]) / (alone["high"] - alone["low"])
    assert ratio > 1.2


def test_the_two_sample_interval_resamples_clusters_in_both_groups():
    """Eight clones of five carry eight observations per group, not forty."""
    rng = seeded(16)
    arms = []
    for mu in (-4.0, -5.0):
        values = np.repeat(rng.normal(mu, 1.0, 8), 5) + rng.normal(0, 0.05, 40)
        steps = np.log2(np.asarray(LADDER))
        reported = []
        for value in values:
            index = int(np.searchsorted(steps, value))
            reported.append(f"<={LADDER[0]}" if index == 0 else
                            f">{LADDER[-1]}" if index >= len(steps) else
                            f"{LADDER[index]}")
        lower, upper, _, _ = core.place_mics(reported, LADDER)
        arms.append((lower, upper))
    clones = np.repeat(np.arange(8), 5)
    rows = np.arange(40)

    def width(group_clusters, reference_clusters):
        result = core.resample_shift(arms[0][0], arms[0][1], group_clusters,
                                     arms[1][0], arms[1][1], reference_clusters,
                                     seeded(17), draws=300)
        return result["high"] - result["low"]

    independent = width(rows, rows)
    assert width(clones, rows) > 1.3 * independent
    assert width(rows, clones) > 1.3 * independent


def test_two_sample_resampling_refuses_mismatched_lengths():
    with pytest.raises(ValueError, match="the same length"):
        core.resample_shift([0.0], [1.0], ["a", "b"], [0.0], [1.0], ["a"], seeded())
    with pytest.raises(ValueError, match="the same length"):
        core.resample_shift([0.0], [1.0], ["a"], [0.0], [1.0], ["a", "b"], seeded())


def test_two_sample_resampling_returns_nothing_when_the_draws_fail():
    """A group holding no interval cannot be fitted in any draw."""
    reference = placed(-5.0, 1.0, 50, seeded(19))
    assert core.resample_shift([np.nan] * 5, [np.nan] * 5, np.arange(5),
                               reference[0], reference[1], np.arange(50),
                               seeded(20), draws=100) is None


def test_two_sample_resampling_refuses_an_empty_group():
    with pytest.raises(ValueError, match="at least one observation"):
        core.resample_shift([], [], [], [0.0], [1.0], ["a"], seeded())
    with pytest.raises(ValueError, match="at least one observation"):
        core.resample_shift([0.0], [1.0], ["a"], [], [], [], seeded())


# named_generator


def test_a_named_stream_is_the_seed_followed_by_the_code_points_of_the_name():
    expected = np.random.default_rng([20260101] + [ord(c) for c in "fit|wild"])
    stream = core.named_generator(20260101, "fit|wild")
    assert np.array_equal(stream.random(5), expected.random(5))


def test_streams_with_different_names_or_seeds_differ():
    first = core.named_generator(20260101, "fit|a").random(5)
    assert not np.array_equal(first, core.named_generator(20260101, "fit|b").random(5))
    assert not np.array_equal(first, core.named_generator(20260102, "fit|a").random(5))


def test_a_stream_needs_a_name_and_a_non_negative_seed():
    with pytest.raises(ValueError, match="non-empty name"):
        core.named_generator(1, "")
    with pytest.raises(ValueError, match="the seed must be a non-negative integer"):
        core.named_generator(-1, "fit|a")


# withhold_reason


def bounds_of(readings):
    lower, upper, _, _ = core.place_mics(readings, LADDER)
    return lower, upper


def test_a_well_supported_fit_is_not_withheld():
    lower, upper = placed(-5.0, 1.0, 500, seeded(18))
    fit = core.fit_censored_normal(lower, upper)
    assert core.withhold_reason(fit, lower, upper, LADDER) is None


def test_a_group_with_no_placed_reading_is_withheld():
    lower, upper = bounds_of([None, ""])
    assert core.withhold_reason(None, lower, upper, LADDER) == "no reading was placed"


def test_readings_in_two_intervals_are_withheld():
    """
    Half at 0.06 and half at 0.12: a vanishing spread centered on the boundary
    gives each interval exactly half its mass, so the likelihood has no finite
    maximum and the optimizer returns whatever spread it stopped at.
    """
    lower, upper = bounds_of(["0.06"] * 20 + ["0.12"] * 20)
    fit = core.fit_censored_normal(lower, upper)
    reason = core.withhold_reason(fit, lower, upper, LADDER)
    assert reason is not None and "too few intervals" in reason


def test_readings_in_three_intervals_support_a_fit():
    lower, upper = bounds_of(["0.06"] * 10 + ["0.12"] * 20 + ["0.25"] * 10)
    fit = core.fit_censored_normal(lower, upper)
    assert core.withhold_reason(fit, lower, upper, LADDER) is None


def test_every_reading_censored_at_the_floor_is_withheld():
    lower, upper = bounds_of(["<=0.008"] * 30)
    fit = core.fit_censored_normal(lower, upper)
    reason = core.withhold_reason(fit, lower, upper, LADDER)
    assert reason is not None and "too few intervals" in reason


def test_a_fit_that_did_not_converge_is_withheld():
    lower, upper = bounds_of(["0.06", "0.12", "0.25"])
    assert core.withhold_reason(None, lower, upper, LADDER) == (
        "the fit did not converge")


@pytest.mark.parametrize("fit", [
    {"mu": core.MU_BOUNDS[0], "sigma": 1.0},
    {"mu": core.MU_BOUNDS[1], "sigma": 1.0},
    {"mu": -4.0, "sigma": float(np.exp(core.LOG_SIGMA_BOUNDS[0]))},
    {"mu": -4.0, "sigma": float(np.exp(core.LOG_SIGMA_BOUNDS[1]))},
])
def test_a_fit_on_the_optimizer_bound_is_withheld(fit):
    lower, upper = bounds_of(["0.06", "0.12", "0.25"])
    assert core.withhold_reason(fit, lower, upper, LADDER) == (
        "a parameter ended on the optimizer's bound")


@pytest.mark.parametrize("readings,parameter", [
    (["<=0.008"] * 200 + ["0.03", ">1"], "mu"),
    (["<=0.008"] * 100 + [">1"] * 100 + ["0.06"], "sigma"),
])
def test_a_fit_the_data_drive_onto_the_optimizer_bound_is_withheld(readings,
                                                                    parameter):
    """
    Both datasets occupy three intervals, so each likelihood has a finite
    maximum, and in both it lies outside the box the optimizer searches: the
    first below a mean of -25, the second above a standard deviation of 20.
    The check must read the same bounds the optimizer stopped on.
    """
    lower, upper = bounds_of(readings)
    fit = core.fit_censored_normal(lower, upper)
    if parameter == "mu":
        assert fit["mu"] == core.MU_BOUNDS[0]
    else:
        assert np.log(fit["sigma"]) == pytest.approx(core.LOG_SIGMA_BOUNDS[1])
    assert core.withhold_reason(fit, lower, upper, LADDER) == (
        "a parameter ended on the optimizer's bound")


@pytest.mark.parametrize("mu,withheld", [
    (np.log2(0.008) - 1.01, True),
    (np.log2(0.008) - 0.99, False),
    (np.log2(1.0) + 0.99, False),
    (np.log2(1.0) + 1.01, True),
])
def test_a_mean_beyond_the_margin_outside_the_series_is_withheld(mu, withheld):
    lower, upper = bounds_of(["0.06", "0.12", "0.25"])
    reason = core.withhold_reason({"mu": mu, "sigma": 1.0}, lower, upper, LADDER)
    assert (reason is not None) is withheld
    if withheld:
        assert "outside the tested series" in reason


# round_up_to_series


def test_a_value_rounds_up_to_the_next_tested_concentration():
    assert core.round_up_to_series(np.log2(0.04), LADDER) == 0.06
    assert core.round_up_to_series(np.log2(0.06), LADDER) == 0.06
    assert core.round_up_to_series(-20.0, LADDER) == 0.008
    assert core.round_up_to_series(0.5, LADDER) is None
