"""
Tests for the interval-censored MIC estimation.

The point of this module is that a plate reading is an interval, not a
measurement, and that taking a median of interval-censored values is biased. The
tests assert both: that the intervals are constructed correctly, and that the
estimator recovers parameters a median cannot.
"""

import numpy as np
import pytest

from mic_model import (concentration_series, fit_censored_linear,
                       fit_censored_normal, log_interval_mass, mic_bounds,
                       mic_recorded, parse_concentration, place_mics)

LADDER = [0.008, 0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1.0]


def test_concentration_text_is_parsed_including_operators():
    """PLATE_LAYOUT stores CONC as text, not the float the schema describes."""
    assert parse_concentration("0.25") == (0.25, None)
    assert parse_concentration("<=0.008") == (0.008, "<=")
    assert parse_concentration(">1.0") == (1.0, ">")
    assert parse_concentration("0")[0] == 0.0
    assert parse_concentration("not a number") == (None, None)


def test_a_reported_value_is_the_interval_below_it():
    """An MIC of 0.25 means growth stopped at 0.25 and not at 0.12."""
    lower, upper = mic_bounds("0.25", LADDER)
    assert lower == pytest.approx(np.log2(0.12), abs=0.01)
    assert upper == pytest.approx(np.log2(0.25), abs=0.01)


def test_lowest_well_is_left_censored():
    lower, upper = mic_bounds("<=0.008", LADDER)
    assert lower == -np.inf
    assert upper == pytest.approx(np.log2(0.008), abs=0.01)


def test_above_the_top_well_is_right_censored():
    lower, upper = mic_bounds(">1.0", LADDER)
    assert lower == pytest.approx(np.log2(1.0), abs=0.01)
    assert upper == np.inf


def test_a_value_not_on_the_ladder_is_refused_rather_than_guessed():
    assert mic_bounds("0.37", LADDER) is None
    assert mic_bounds("", LADDER) is None


def test_ladders_are_read_from_plate_layout(data_dir):
    series, irregular, unreadable = concentration_series()
    assert ("UKMYC6", "BDQ") in series
    assert series[("UKMYC6", "BDQ")][0] == pytest.approx(0.008)
    assert not irregular, "the fixture ladders are clean doubling series"


@pytest.mark.parametrize("true_mu, true_sd", [
    (-5.0, 1.0),
    (-6.5, 1.0),   # heavily left-censored
    (-7.5, 1.2),   # mostly below the plate
    (-2.0, 1.5),
])
def test_estimator_recovers_known_parameters(true_mu, true_sd):
    rng = np.random.default_rng(99)
    steps = np.log2(np.asarray(LADDER))

    def report(value):
        for index, step in enumerate(steps):
            if value <= step:
                return f"<={LADDER[0]}" if index == 0 else f"{LADDER[index]}"
        return f">{LADDER[-1]}"

    draws = rng.normal(true_mu, true_sd, 4000)
    bounds = [mic_bounds(report(x), LADDER) for x in draws]
    fit = fit_censored_normal([b[0] for b in bounds], [b[1] for b in bounds])
    assert fit["mu"] == pytest.approx(true_mu, abs=0.15)
    assert fit["sigma"] == pytest.approx(true_sd, abs=0.15)


def test_median_is_biased_where_the_estimator_is_not():
    """
    The reason this estimator exists. Under heavy left-censoring the median of
    reported values sticks near the plate floor while the fit follows the truth.
    """
    rng = np.random.default_rng(7)
    true_mu = -7.5
    steps = np.log2(np.asarray(LADDER))

    def report(value):
        for index, step in enumerate(steps):
            if value <= step:
                return f"<={LADDER[0]}" if index == 0 else f"{LADDER[index]}"
        return f">{LADDER[-1]}"

    draws = rng.normal(true_mu, 1.2, 4000)
    reported = [report(x) for x in draws]
    bounds = [mic_bounds(value, LADDER) for value in reported]
    fit = fit_censored_normal([b[0] for b in bounds], [b[1] for b in bounds])
    naive = float(np.median([np.log2(float(value.replace("<=", "").replace(">", "")))
                             for value in reported]))

    assert abs(fit["mu"] - true_mu) < 0.2
    assert abs(naive - true_mu) > 0.4
    assert naive > fit["mu"], "the median is biased upward, toward the plate floor"


def test_a_missing_mic_is_not_a_reading():
    """A plate that produced no MIC stores a null, and the null arrives as None,
    as a float nan or as a pandas missing value depending on the pandas version
    reading the parquet file."""
    import pandas as pd

    for missing in (None, np.nan, pd.NA, pd.NaT, "", "nan", "None", "<NA>"):
        assert not mic_recorded(missing)
    for reading in ("0.25", "<=0.008", ">1", " 0.5 "):
        assert mic_recorded(reading)


def test_place_mics_separates_a_missing_mic_from_a_value_off_the_series():
    """A sample with no MIC is a missing measurement. A value that is present
    and does not match the tested series means the series is wrong, which is a
    different fault and must not be hidden among the missing ones."""
    series = {("UKMYC6", "BDQ"): LADDER}
    texts = ["0.25", None, "<=0.008", "0.2", ">1", np.nan]
    designs = ["UKMYC6"] * 6
    lowers, uppers, absent, off_series = place_mics(texts, designs, series, "BDQ")

    assert absent == 2
    assert off_series == ["0.2"]
    placed = sum(1 for value in lowers if not np.isnan(value))
    assert placed == 3
    assert placed + absent + len(off_series) == len(texts)
    assert np.isnan(lowers[1]) and np.isnan(uppers[1])
    assert np.isnan(lowers[3]) and np.isnan(uppers[3])


def test_place_mics_counts_a_design_with_no_tested_series_as_absent():
    """The cohort holds samples with no plate design for a drug. They have no
    MIC either, so they are missing measurements rather than bad values."""
    lowers, uppers, absent, off_series = place_mics(
        ["0.25"], [None], {("UKMYC6", "BDQ"): LADDER}, "BDQ")
    assert absent == 1
    assert off_series == []
    assert np.isnan(lowers[0])


def test_the_interval_mass_survives_the_far_tail():
    """Subtracting two cdfs underflows to zero far out in a tail, which makes
    the log likelihood infinite there and corrupts the gradient. The log-space
    form returns a finite value for the same interval."""
    from scipy.stats import norm

    low, high = 30.0, 31.0
    assert norm.cdf(high) - norm.cdf(low) == 0.0
    mass = log_interval_mass([low], [high])
    assert np.isfinite(mass[0])
    assert mass[0] == pytest.approx(np.log(norm.sf(low) - norm.sf(high)), rel=1e-9)


def test_the_interval_mass_matches_the_direct_difference_where_that_is_safe():
    """In the body of the distribution the two forms agree, so the log-space
    form is the same likelihood and not a different model."""
    from scipy.stats import norm

    low = np.array([-np.inf, -2.0, -0.5, 1.0, 2.0])
    high = np.array([-2.0, -0.5, 1.0, 2.0, np.inf])
    direct = np.log(norm.cdf(high) - norm.cdf(low))
    assert log_interval_mass(low, high) == pytest.approx(direct, rel=1e-12)


def test_the_estimator_recovers_parameters_with_far_tail_observations():
    """A group holding a few observations far from its centre is where the
    difference of cdfs loses precision. The planted parameters must still come
    back."""
    rng = np.random.default_rng(29)
    ladder = np.log2(np.array(LADDER))
    values = np.concatenate([rng.normal(-5.3, 1.1, 4000), rng.normal(1.5, 0.3, 20)])
    lower, upper = [], []
    for value in values:
        if value <= ladder[0]:
            lower.append(-np.inf)
            upper.append(ladder[0])
        elif value > ladder[-1]:
            lower.append(ladder[-1])
            upper.append(np.inf)
        else:
            index = int(np.searchsorted(ladder, value))
            lower.append(ladder[index - 1])
            upper.append(ladder[index])
    fit = fit_censored_normal(lower, upper)
    assert fit["mu"] == pytest.approx(-5.3, abs=0.12)
    assert fit["sigma"] == pytest.approx(1.1, abs=0.15)


def test_a_censored_bound_uses_the_tail_directly():
    """A left-censored observation contributes the lower tail and a
    right-censored one the upper tail, with no interval arithmetic."""
    from scipy.stats import norm

    assert log_interval_mass([-np.inf], [-1.5])[0] == pytest.approx(norm.logcdf(-1.5))
    assert log_interval_mass([2.5], [np.inf])[0] == pytest.approx(norm.logsf(2.5))


SITE_EFFECT = {0: 0.0, 1: 1.5, 2: -0.5}


def confounded_cohort(seed, shift=1.0, carriers_at_high_site=50):
    """A reference group across three sites of different baseline, and a
    carrier group sitting mostly at the high one, which is the shape the real
    cohort has for minor alleles."""
    rng = np.random.default_rng(seed)
    ladder = np.log2(np.array(LADDER))
    site = np.concatenate([rng.integers(0, 3, 3000),
                           np.array([1] * carriers_at_high_site + [0] * 5 + [2] * 5)])
    carrier = np.concatenate([np.zeros(3000),
                              np.ones(carriers_at_high_site + 10)])
    centre = np.array([-5.3 + SITE_EFFECT[s] + shift * c
                       for s, c in zip(site, carrier)])
    lower, upper = [], []
    for value in rng.normal(centre, 1.1):
        if value <= ladder[0]:
            lower.append(-np.inf)
            upper.append(ladder[0])
        elif value > ladder[-1]:
            lower.append(ladder[-1])
            upper.append(np.inf)
        else:
            index = int(np.searchsorted(ladder, value))
            lower.append(ladder[index - 1])
            upper.append(ladder[index])
    design = np.column_stack([np.ones(len(site)), carrier,
                              (site == 1).astype(float), (site == 2).astype(float)])
    return np.array(lower), np.array(upper), design, carrier


def test_the_linear_fit_recovers_planted_coefficients():
    lower, upper, design, _ = confounded_cohort(17)
    fit = fit_censored_linear(lower, upper, design)
    assert fit["beta"][0] == pytest.approx(-5.3, abs=0.2)
    assert fit["beta"][1] == pytest.approx(1.0, abs=0.35)
    assert fit["beta"][2] == pytest.approx(SITE_EFFECT[1], abs=0.2)
    assert fit["beta"][3] == pytest.approx(SITE_EFFECT[2], abs=0.2)
    assert fit["sigma"] == pytest.approx(1.1, abs=0.15)


def test_holding_site_constant_undoes_a_group_concentrated_at_one_site():
    """The difference of group means is biased upward where the carriers sit at
    a site whose baseline is high. The fit with site in the design is not."""
    lower, upper, design, carrier = confounded_cohort(17)
    naive = (fit_censored_normal(lower[carrier == 1], upper[carrier == 1])["mu"]
             - fit_censored_normal(lower[carrier == 0], upper[carrier == 0])["mu"])
    adjusted = fit_censored_linear(lower, upper, design)["beta"][1]
    assert naive > 1.6
    assert abs(adjusted - 1.0) < abs(naive - 1.0)


def test_a_design_of_ones_reproduces_the_single_mean():
    """With nothing to hold constant, the linear fit is the plain fit."""
    lower, upper, _, _ = confounded_cohort(5)
    plain = fit_censored_normal(lower, upper)
    linear = fit_censored_linear(lower, upper, np.ones((len(lower), 1)))
    assert linear["beta"][0] == pytest.approx(plain["mu"], abs=1e-3)
    assert linear["sigma"] == pytest.approx(plain["sigma"], abs=1e-3)


def test_a_design_of_the_wrong_length_is_refused():
    """Numpy would raise on the shape mismatch anyway, further in and with a
    message about broadcasting, so the message is part of what is asserted."""
    lower, upper, design, _ = confounded_cohort(5)
    with pytest.raises(ValueError, match="one row per observation"):
        fit_censored_linear(lower, upper, design[:-1])


def linear_objective(lower, upper, design):
    """The objective fit_censored_linear minimises, and its gradient, as two
    separate callables for checking one against the other."""
    from mic_model import log_density_ratios

    def value(parameters):
        mu = design @ parameters[:-1]
        sigma = np.exp(parameters[-1])
        return -np.sum(log_interval_mass((lower - mu) / sigma,
                                         (upper - mu) / sigma))

    def gradient(parameters):
        mu = design @ parameters[:-1]
        sigma = np.exp(parameters[-1])
        low = (lower - mu) / sigma
        high = (upper - mu) / sigma
        at_low, at_high = log_density_ratios(low, high,
                                             log_interval_mass(low, high))
        finite_low = np.where(np.isfinite(low), low, 0.0)
        finite_high = np.where(np.isfinite(high), high, 0.0)
        return np.concatenate([design.T @ ((at_high - at_low) / sigma),
                               [np.sum(finite_high * at_high
                                       - finite_low * at_low)]])

    return value, gradient


def test_the_closed_form_gradient_matches_a_numerical_one():
    """The gradient is supplied to the optimiser in closed form. A sign or a
    missing term would send it to a different place, so it is checked against a
    difference quotient away from the solution as well as at it."""
    from scipy.optimize import approx_fprime

    lower, upper, design, _ = confounded_cohort(17)
    value, gradient = linear_objective(lower, upper, design)
    rng = np.random.default_rng(5)
    truth = np.array([-5.3, 1.0, SITE_EFFECT[1], SITE_EFFECT[2], np.log(1.1)])
    for _ in range(10):
        point = truth + rng.normal(0, 0.4, len(truth))
        analytic = gradient(point)
        numeric = approx_fprime(point, value, 1e-6)
        assert np.abs(analytic - numeric).max() < 1e-3 * max(
            1.0, float(np.abs(numeric).max()))


def test_the_gradient_survives_observations_censored_at_both_ends():
    """An infinite bound carries no density. Multiplied against it before it is
    discarded it would carry a nan instead, and one nan ruins the sum."""
    lower = np.array([-np.inf, -np.inf, 2.0, -1.0, np.inf * 0 + 3.0])
    upper = np.array([-2.0, 1.0, np.inf, 0.0, np.inf])
    design = np.ones((5, 1))
    _, gradient = linear_objective(lower, upper, design)
    for parameter in (0.0, -3.0, 4.0):
        assert np.isfinite(gradient(np.array([parameter, np.log(1.2)]))).all()
    fit = fit_censored_linear(lower, upper, design)
    assert fit is not None
    assert np.isfinite(fit["beta"]).all()


def test_the_fit_stops_at_a_stationary_point():
    """A looser stopping rule leaves the gradient far from zero and the estimate
    short of the maximum by enough to move a shift in its third decimal. The
    solution is checked against a second optimiser started from the same place.
    """
    from scipy.optimize import minimize

    lower, upper, design, _ = confounded_cohort(17)
    value, gradient = linear_objective(lower, upper, design)
    fit = fit_censored_linear(lower, upper, design)
    solution = np.concatenate([fit["beta"], [np.log(fit["sigma"])]])
    assert np.abs(gradient(solution)).max() < 1e-3

    second = minimize(value, solution, method="Nelder-Mead",
                      options={"xatol": 1e-9, "fatol": 1e-11, "maxfev": 20000})
    assert value(solution) <= second.fun + 1e-6
    assert np.abs(second.x - solution).max() < 1e-3


def test_the_density_ratio_holds_where_the_mass_underflows():
    """Forty standard deviations out, the mass is exp(-805) and the density is
    exp(-801). Both underflow to zero as plain numbers and their quotient
    becomes a nan, while the quotient itself is about forty.

    For an observation censored below, the derivative of the negative log
    likelihood with respect to the mean is that quotient, and it approaches the
    distance to the bound.
    """
    for distance in (8.0, 40.0, 200.0):
        lower = np.array([-np.inf])
        upper = np.array([-distance])
        _, gradient = linear_objective(lower, upper, np.ones((1, 1)))
        slope = gradient(np.array([0.0, 0.0]))[0]
        assert np.isfinite(slope)
        assert slope == pytest.approx(distance, rel=0.02)
