"""
Tests for the interval-censored MIC estimation.

The point of this module is that a plate reading is an interval, not a
measurement, and that taking a median of interval-censored values is biased. The
tests assert both: that the intervals are constructed correctly, and that the
estimator recovers parameters a median cannot.
"""

import numpy as np
import pandas as pd
import pytest

from mic_model import (PLANTED, bootstrap_mu, concentration_series,
                       excluded_sensitivity, fit_censored_linear,
                       fit_censored_normal, log_interval_mass, mic_bounds,
                       mic_recorded, parse_concentration, place_mics,
                       simulate_reports, stream, tiny_intervals,
                       validate_estimator)

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


def placed_cohort(seed, samples=300, clone=0, clone_shift=2.0, grouped=True):
    """A cohort on the ladder. The first clone samples sit clone_shift doublings
    above the rest, and share one cluster when grouped and sit in clusters of
    their own when not, so two frames from one seed hold the same MICs and
    differ in nothing but how many independent units the resampling has."""
    rng = np.random.default_rng(seed)
    ladder = np.log2(np.array(LADDER))
    centre = np.where(np.arange(samples) < clone, -5.3 + clone_shift, -5.3)
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
    names = [f"c{index}" for index in range(samples)]
    if grouped:
        names[:clone] = ["outbreak"] * clone
    return pd.DataFrame({"CLUSTER": names, "lower_BDQ": lower, "upper_BDQ": upper})


def test_each_interval_draws_from_its_own_stream():
    """An interval computed after another estimate must match the same interval
    computed on its own, or a figure cannot be reproduced without rerunning
    everything that preceded it."""
    frame = placed_cohort(5)
    alone = bootstrap_mu(frame, "BDQ", stream("BDQ Rv0678 promoter"), draws=20)

    spent = stream("BDQ Rv0678 promoter")
    spent.integers(0, 100, 5000)
    assert bootstrap_mu(frame, "BDQ", spent, draws=20) != alone

    other = stream("BDQ pepQ solo")
    other.integers(0, 100, 5000)
    again = bootstrap_mu(frame, "BDQ", stream("BDQ Rv0678 promoter"), draws=20)
    assert again == alone
    assert stream("BDQ Rv0678 promoter").integers(0, 1 << 30, 3).tolist() != \
        stream("CFZ Rv0678 promoter").integers(0, 1 << 30, 3).tolist()


def test_the_interval_resamples_clusters_rather_than_isolates():
    """A group whose isolates come from one outbreak carries less information
    than the same count of unrelated isolates, which the interval has to show."""
    spread = {}
    for grouped in (False, True):
        interval = bootstrap_mu(placed_cohort(9, clone=100, grouped=grouped),
                                "BDQ", stream("outbreak"), draws=150)
        spread[grouped] = interval["high"] - interval["low"]
    assert spread[True] > 2 * spread[False]


def test_every_interval_in_the_report_is_seeded_by_its_own_name(data_dir, tmp_path,
                                                                monkeypatch):
    """Each interval has to be handed the stream named for the quantity it
    belongs to. With one generator threaded through the run, an interval depends
    on the estimates computed before it."""
    import cohort
    import mic_model

    drawn = []

    def record(frame, drug, rng, draws=None):
        drawn.append(int(rng.integers(0, 1 << 30)))
        return {"low": -5.0, "high": -4.0, "draws": 1}

    monkeypatch.setattr(mic_model, "REPORT", tmp_path / "report.txt")
    monkeypatch.setattr(mic_model, "ESTIMATES", tmp_path / "estimates.csv")
    monkeypatch.setattr(mic_model, "bootstrap_mu", record)
    mic_model.main()

    names = [f"{drug} {group}" for drug in cohort.DRUGS
             for group in cohort.GROUP_ORDER]
    names += [f"{drug} {lineage} {label}" for drug in cohort.DRUGS
              for lineage in ["lineage1", "lineage2", "lineage3", "lineage4"]
              for label in ["loss of function", "substitution"]]
    named = {int(mic_model.stream(name).integers(0, 1 << 30)) for name in names}

    assert drawn, "the run computed no interval"
    assert set(drawn) <= named, "an interval was handed a stream of no quantity"
    assert len(set(drawn)) == len(drawn), "two intervals drew from one stream"


def test_the_planted_validation_recovers_every_planted_distribution():
    """The table in the results document is this function's output, so what it
    claims has to hold for every row of it."""
    rows = validate_estimator(LADDER, stream("planted parameters"))
    assert len(rows) == len(PLANTED)
    for row in rows:
        assert row["fitted mean"] == pytest.approx(row["true mean"], abs=0.15)
        assert row["fitted sd"] == pytest.approx(row["true sd"], abs=0.15)
        assert row["naive median"] > row["fitted mean"], \
            "the median of reported values is biased toward the plate floor"
    heaviest = max(rows, key=lambda row: row["left-censored %"])
    assert heaviest["left-censored %"] > 50
    assert heaviest["naive median"] - heaviest["fitted mean"] > 0.4


def test_a_draw_below_the_plate_is_reported_as_censored():
    """A plate cannot report a value it never tested, so a draw under the first
    well comes back censored at that well and one over the last above it."""
    rng = np.random.default_rng(3)
    low = simulate_reports(-20.0, 0.01, LADDER, rng, draws=20)
    high = simulate_reports(20.0, 0.01, LADDER, rng, draws=20)
    assert set(low) == {f"<={LADDER[0]}"}
    assert set(high) == {f">{LADDER[-1]}"}

    middle = simulate_reports(np.log2(LADDER[3]) - 0.5, 0.01, LADDER, rng, draws=20)
    assert set(middle) == {f"{LADDER[3]}"}, \
        "a draw between two wells is reported at the one that inhibited growth"


def _reference_frame(plan, excluded=0, seed=4, sd=0.5):
    """A reference group of (site, isolates, mean) rows placed on the ladder,
    plus `excluded` rows carrying no MIC at all."""
    rng = np.random.default_rng(seed)
    rows = []
    for site, isolates, mu in plan:
        for value in simulate_reports(mu, sd, LADDER, rng, draws=isolates):
            lower, upper = mic_bounds(value, LADDER)
            rows.append({"GROUP": "reference", "SITEID": site,
                         "lower_BDQ": lower, "upper_BDQ": upper})
    for index in range(excluded):
        rows.append({"GROUP": "reference", "SITEID": "99",
                     "lower_BDQ": np.nan, "upper_BDQ": np.nan})
    return pd.DataFrame(rows)


def test_the_excluded_rows_move_a_mean_by_their_share_of_it():
    """An excluded row acts only through the value it would have had, so the
    movement is its share of the group times the gap."""
    frame = _reference_frame([("A", 400, -5.0), ("B", 400, -5.0), ("C", 400, -4.0),
                              ("D", 400, -7.5)], excluded=100)
    result = excluded_sensitivity(frame, "BDQ", minimum=100)

    assert (result["reference"], result["fitted"], result["excluded"]) == (1700, 1600, 100)
    assert result["at site"] == "D", \
        "the widest gap is the one furthest from the group mean in either direction"
    assert result["widest site gap"] > 0
    share = 100 / 1700
    assert result["movement at that gap"] == pytest.approx(
        round(share * result["widest site gap"], 3), abs=0.001)
    assert result["movement at 3 doublings"] == pytest.approx(round(share * 3, 3),
                                                              abs=0.001)


def test_a_site_too_small_to_fit_is_left_out_of_the_widest_gap():
    """A fit on a handful of intervals is not a site effect, however far it lands
    from the group."""
    frame = _reference_frame([("A", 400, -5.0), ("B", 400, -5.0), ("C", 400, -4.0),
                              ("tiny", 20, 0.0)])
    result = excluded_sensitivity(frame, "BDQ", minimum=100)
    assert result["at site"] == "C"


def test_an_interval_in_the_far_tail_is_where_the_mass_underflows():
    """The mass of an interval far from the fitted mean cannot be formed as a
    difference of two cumulative functions, which is what the log-space
    evaluation exists for."""
    central = _reference_frame([("A", 600, -5.0)], sd=0.3)
    assert tiny_intervals(central, "BDQ", floor=1e-08)["below the floor"] == 0

    doctored = pd.concat([central, _reference_frame([("A", 3, 3.0)], sd=0.3)],
                         ignore_index=True)
    result = tiny_intervals(doctored, "BDQ", floor=1e-08)
    assert result["below the floor"] >= 3
    assert 0 < result["smallest"] < 1e-08, \
        "a mass this small is a number, and a difference of two cumulative "\
        "functions would return zero for it"


# The estimator lives in the micecoff package. These tests hold mic_model to it.


MOVED = ["parse_concentration", "mic_bounds", "mic_recorded", "log_interval_mass",
         "log_density_ratios", "fit_censored_normal", "fit_censored_linear",
         "simulate_reports"]


def test_the_estimator_this_module_exports_is_the_packages():
    """
    The analysis modules import these names from mic_model. Each must be the
    package's own function, so the pipeline and the command line tool run one
    estimator and a fix to it reaches both.
    """
    import mic_model
    from micecoff import core

    for name in MOVED:
        assert getattr(mic_model, name) is getattr(core, name), name


def test_the_streams_follow_the_documented_seed():
    """
    Every published interval was drawn from a generator seeded by 20260101 and
    the code points of the quantity's name. A change to that scheme moves every
    interval in the released tables, so it is held here to the numbers numpy
    gives for the scheme written out.
    """
    import mic_model

    assert mic_model.SEED == 20260101
    for name in ("BDQ reference", "CFZ Rv0678 promoter", "BDQ lineage4 substitution"):
        expected = np.random.default_rng([20260101] + [ord(c) for c in name])
        assert np.array_equal(stream(name).random(4), expected.random(4))


def test_a_skipped_dilution_and_an_unreadable_well_are_reported(tmp_path,
                                                                 monkeypatch):
    """
    A layout whose bedaquiline series skips 0.5 and carries a positive-control
    well recorded as 0, beside a clean clofazimine series. The skipped dilution
    is reported with its ratios rounded to three decimals, where 0.25 over 0.12
    is 2.0833, the control well as unreadable, and both series are still
    returned.
    """
    import cohort

    layout = pd.DataFrame({
        "PLATEDESIGN": ["UKMYC6"] * 9,
        "DRUG": ["BDQ"] * 6 + ["CFZ"] * 3,
        "CONC": ["<=0.12", "0.25", "1", "2", ">2", "0",
                 "<=0.03", "0.06", ">0.06"],
    })
    layout.to_parquet(tmp_path / "PLATE_LAYOUT.parquet")
    monkeypatch.setattr(cohort, "DATA", tmp_path)

    series, irregular, unreadable = concentration_series()
    assert series[("UKMYC6", "BDQ")] == [0.12, 0.25, 1.0, 2.0]
    assert series[("UKMYC6", "CFZ")] == [0.03, 0.06]
    assert irregular == [("UKMYC6", "BDQ", [2.083, 4.0, 2.0])]
    assert unreadable == [("UKMYC6", "BDQ", "0")]


def test_the_module_run_as_a_script_imports_the_repositorys_estimator(tmp_path):
    """
    python code/mic_model.py puts code/ on the path and nothing else. The module
    then has to reach the micecoff package in this repository. A decoy package
    of the same name on PYTHONPATH stands in for a copy installed elsewhere,
    which must not take its place.
    """
    import os
    import pathlib
    import subprocess
    import sys

    root = pathlib.Path(__file__).resolve().parents[1]
    decoy = tmp_path / "decoy" / "micecoff"
    decoy.mkdir(parents=True)
    (decoy / "__init__.py").write_text("")
    program = (
        "import sys; sys.path.insert(0, sys.argv[1]); import mic_model, micecoff; "
        "print(micecoff.__file__)"
    )
    environment = dict(os.environ, PYTHONPATH=str(tmp_path / "decoy"),
                       PYTHONDONTWRITEBYTECODE="1")
    result = subprocess.run([sys.executable, "-c", program, str(root / "code")],
                            cwd=tmp_path, env=environment, capture_output=True,
                            text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    loaded = pathlib.Path(result.stdout.strip()).resolve()
    assert loaded == root / "micecoff" / "__init__.py"
