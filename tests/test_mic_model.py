"""
Tests for the interval-censored MIC estimation.

The point of this module is that a plate reading is an interval, not a
measurement, and that taking a median of interval-censored values is biased. The
tests assert both: that the intervals are constructed correctly, and that the
estimator recovers parameters a median cannot.
"""

import numpy as np
import pytest

from mic_model import (concentration_series, fit_censored_normal, mic_bounds,
                       parse_concentration)

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
