"""
Tests for the statistical corrections.

Each estimator is checked against data with a known answer planted in it: a
known odds ratio, a known amount of heterogeneity, a known clonal artefact.
"""

import numpy as np
import pandas as pd
import pytest

import analyse_groups
import cluster_adjust
import discovery


def test_mantel_haenszel_recovers_a_known_odds_ratio():
    """Three sites with very different baselines and one shared genotype effect."""
    rng = np.random.default_rng(7)
    true_odds_ratio = 5.0
    rows = []
    for site, baseline in [("01", 0.02), ("02", 0.20), ("03", 0.08)]:
        for group, count in [("reference", 400), ("exposed", 60)]:
            odds = baseline / (1 - baseline)
            if group == "exposed":
                odds *= true_odds_ratio
            probability = odds / (1 + odds)
            for _ in range(count):
                rows.append({"SITEID": site, "GROUP": group,
                             "resistant_BDQ": rng.random() < probability})
    frame = pd.DataFrame(rows)
    frame["EXPOSED"] = frame.GROUP.eq("exposed")
    estimate = cluster_adjust.mh_odds_ratio(frame, "BDQ")
    assert estimate == pytest.approx(true_odds_ratio, rel=0.4)


def test_homogeneity_test_detects_an_effect_that_differs_by_site():
    rng = np.random.default_rng(3)
    analyse_groups._lines.clear()
    rows = []
    for site, odds_ratio in [("01", 1.0), ("02", 12.0), ("03", 1.0)]:
        for group, count in [("reference", 400), ("exposed", 80)]:
            odds = 0.05 / 0.95
            if group == "exposed":
                odds *= odds_ratio
            probability = odds / (1 + odds)
            for _ in range(count):
                rows.append({"SITEID": site, "GROUP": group,
                             "resistant_BDQ": rng.random() < probability})
    analyse_groups.stratified(pd.DataFrame(rows), "exposed", "reference", "BDQ", "SITEID")
    report = "\n".join(analyse_groups._lines)
    assert "differs across strata" in report


def test_collapsing_clusters_removes_an_effect_driven_by_one_clone():
    """
    The exposed group is one clone of 40, all resistant, plus 12 diverse
    isolates with no real effect. The truth is no effect.
    """
    rng = np.random.default_rng(5)
    rows = []
    for index in range(600):
        rows.append({"CLUSTER": f"ref{index}", "EXPOSED": False, "SITEID": "10",
                     "resistant_BDQ": rng.random() < 0.05})
    for _ in range(40):
        rows.append({"CLUSTER": "clone", "EXPOSED": True, "SITEID": "10",
                     "resistant_BDQ": True})
    for index in range(12):
        rows.append({"CLUSTER": f"solo{index}", "EXPOSED": True, "SITEID": "10",
                     "resistant_BDQ": rng.random() < 0.05})
    frame = pd.DataFrame(rows)

    naive = cluster_adjust.mh_odds_ratio(frame, "BDQ")
    collapsed = cluster_adjust.collapsed_draws(frame, "BDQ", rng)
    assert naive > 20, "the uncorrected estimate should be badly inflated"
    assert collapsed["median"] < 5, "collapsing should remove most of the inflation"


def test_prediction_interval_covers_a_future_estimate():
    """
    The interval must carry both the discovery estimate's uncertainty and the
    held-out half's sampling error. Fixing it at the discovery point estimate
    makes it too narrow and turns sampling variation into recorded failures.
    """
    rng = np.random.default_rng(33)
    ladders = [[0.008, 0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1.0]]
    true_shift, reference_mu, reference_sd, exposed_sd = 2.4, -5.3, 1.1, 1.4
    discovery_n, held_out_n, reference_n = 42, 18, 400

    null = discovery.simulate_shifts(held_out_n, reference_n, reference_mu,
                                     reference_sd, exposed_sd, 0.0, ladders, rng,
                                     draws=60)
    covered = trials = 0
    for _ in range(25):
        estimate = discovery.simulate_shifts(discovery_n, reference_n, reference_mu,
                                             reference_sd, exposed_sd, true_shift,
                                             ladders, rng, draws=1)
        bootstrap = discovery.simulate_shifts(discovery_n, reference_n, reference_mu,
                                              reference_sd, exposed_sd,
                                              float(estimate[0]), ladders, rng, draws=40)
        predicted = discovery.simulate_shifts(held_out_n, reference_n, reference_mu,
                                              reference_sd, exposed_sd, bootstrap,
                                              ladders, rng, draws=60)
        low, high = discovery.predictive_interval(
            predicted, discovery_draws=bootstrap,
            validation_spread=float(null.std(ddof=1)))
        future = discovery.simulate_shifts(held_out_n, reference_n, reference_mu,
                                           reference_sd, exposed_sd, true_shift,
                                           ladders, rng, draws=1)
        trials += 1
        covered += low <= float(future[0]) <= high
    assert covered / trials >= 0.85, \
        f"coverage {covered}/{trials} is well below the nominal 95%"


def test_predictive_interval_is_never_narrower_than_the_combined_error():
    """
    The interval must be at least as wide as the discovery and held-out errors
    combined. The simulated spread can understate that, because the estimator's
    spread depends on where the distribution sits relative to the plate, and a
    narrow interval turns ordinary sampling variation into recorded failures.
    """
    rng = np.random.default_rng(1)
    # simulated estimates deliberately too tightly spread
    simulated = rng.normal(2.4, 0.10, 500)
    discovery_draws = rng.normal(2.4, 0.25, 500)
    validation_spread = 0.34

    low, high = discovery.predictive_interval(simulated)
    naive_width = high - low

    low, high = discovery.predictive_interval(
        simulated, discovery_draws=discovery_draws,
        validation_spread=validation_spread)
    corrected_width = high - low

    combined = float(np.hypot(discovery_draws.std(ddof=1), validation_spread))
    assert corrected_width == pytest.approx(2 * 1.96 * combined, rel=0.02)
    assert corrected_width > naive_width * 2, \
        "the floor on the width was not applied"


def test_predictive_interval_keeps_the_wider_of_the_two():
    """Where the simulation is already wide enough, the floor must not shrink it."""
    rng = np.random.default_rng(2)
    simulated = rng.normal(2.4, 0.80, 500)
    discovery_draws = rng.normal(2.4, 0.10, 500)

    low, high = discovery.predictive_interval(
        simulated, discovery_draws=discovery_draws, validation_spread=0.10)
    assert (high - low) == pytest.approx(2 * 1.96 * simulated.std(ddof=1), rel=0.02)


def test_power_threshold_delivers_the_power_it_claims():
    from scipy.stats import fisher_exact

    rng = np.random.default_rng(11)
    exposed_n, reference_n, rate = 40, 2000, 0.02
    threshold = discovery.detectable_odds_ratio(exposed_n, reference_n, rate, rng)
    assert threshold is not None

    odds = rate / (1 - rate) * threshold
    probability = odds / (1 + odds)
    significant = 0
    for _ in range(300):
        a = int(rng.binomial(exposed_n, probability))
        b = int(rng.binomial(reference_n, rate))
        significant += fisher_exact([[a, exposed_n - a], [b, reference_n - b]])[1] < 0.05
    assert 0.7 <= significant / 300 <= 0.92, "the claimed 80% power is not delivered"
