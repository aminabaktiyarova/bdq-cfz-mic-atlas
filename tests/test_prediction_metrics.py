"""
Tests for the genotype rules evaluated as tests for resistance.

The arithmetic is checked against a confusion matrix built by hand, the rules
are checked against the cohort definitions, and the intervals are checked to
resample clusters rather than isolates.
"""

import numpy as np
import pandas as pd
import pytest

import cohort
import prediction_metrics as metrics_module
from cluster_adjust import add_clusters


def test_the_metrics_follow_from_the_cells():
    counts = {"true_positive": 30, "false_positive": 10,
              "false_negative": 20, "true_negative": 40}
    values = metrics_module.metrics(counts)
    assert values["sensitivity"] == pytest.approx(30 / 50)
    assert values["specificity"] == pytest.approx(40 / 50)
    assert values["ppv"] == pytest.approx(30 / 40)
    assert values["npv"] == pytest.approx(40 / 60)


def test_a_missing_denominator_gives_no_ratio():
    """A rule that calls nothing has no positive predictive value, which is not
    the same as one of zero."""
    values = metrics_module.metrics({"true_positive": 0, "false_positive": 0,
                                     "false_negative": 5, "true_negative": 95})
    assert values["ppv"] is None
    assert values["sensitivity"] == pytest.approx(0.0)
    assert values["specificity"] == pytest.approx(1.0)


def test_the_cells_are_counted_from_the_two_columns():
    """The four counts differ from one another, so swapping any two is
    visible."""
    called = pd.Series([True] * 5 + [False] * 5)
    resistant = pd.Series([True, True, True, False, False,
                           True, False, False, False, False])
    assert metrics_module.confusion(called, resistant) == {
        "true_positive": 3, "false_positive": 2,
        "false_negative": 1, "true_negative": 4}


def test_a_perfect_rule_and_a_useless_one():
    resistant = pd.Series([True] * 20 + [False] * 80)
    perfect = metrics_module.metrics(metrics_module.confusion(resistant, resistant))
    assert perfect["sensitivity"] == 1.0 and perfect["specificity"] == 1.0
    inverted = metrics_module.metrics(metrics_module.confusion(~resistant, resistant))
    assert inverted["sensitivity"] == 0.0 and inverted["specificity"] == 0.0


def test_the_rules_agree_with_the_cohort_definitions(data_dir):
    mutations = cohort.load_mutations()
    status = cohort.build_status(mutations)
    flags = metrics_module.rule_flags(mutations, status.index)

    assert (flags["Rv0678 loss"] <= flags["Rv0678 any"]).all()
    assert (flags["Rv0678 any"] <= flags["major variant"]).all()
    assert (flags["major variant"] <= flags["major or minor"]).all()

    solo = status.index[status.IS_SOLO]
    assert flags.loc[solo, "major variant"].all()
    reference = status.index[status.IS_REFERENCE]
    assert not flags.loc[reference, "major variant"].any()
    assert not flags.loc[reference, "major or minor"].any()


def test_an_unreadable_gene_is_not_a_called_variant(data_dir, expected):
    """A null call is not a detected variant, so no rule calls it."""
    mutations = cohort.load_mutations()
    status = cohort.build_status(mutations)
    flags = metrics_module.rule_flags(mutations, status.index)
    target = mutations[mutations.GENE.isin(cohort.BDQ_GENES)]
    nulls = set(target[target.IS_NULL_CALL].UNIQUEID)
    called = set(target[target.IS_MINOR | target.IS_HET_CALL | target.REAL_MAJOR].UNIQUEID)
    only_null = sorted(nulls - called)
    assert len(only_null) == expected["n_null_call"]
    assert not flags.loc[only_null].any().any()


def test_admitting_minor_alleles_cannot_lower_sensitivity(data_dir):
    """The fourth rule is the first plus more isolates, so it calls at least as
    many resistant isolates."""
    mutations = cohort.load_mutations()
    status = cohort.build_status(mutations)
    frame = add_clusters(cohort.assemble(status), mutations)
    flags = metrics_module.rule_flags(mutations, frame.index)
    rows = metrics_module.evaluate(frame, flags, "BDQ")
    by_rule = {row["rule"]: row for row in rows}
    assert by_rule["major or minor"]["true_positive"] >= \
        by_rule["major variant"]["true_positive"]
    assert by_rule["major or minor"]["true_negative"] <= \
        by_rule["major variant"]["true_negative"]


def test_the_interval_resamples_clusters():
    """One dataset, two labellings: ten clonal groups against a hundred
    independent isolates."""
    rng = np.random.default_rng(5)
    rows = []
    for cluster in range(10):
        resistant = rng.random() < 0.5
        for member in range(10):
            rows.append({"CLUSTER": f"c{cluster}", "called": True,
                         "resistant": resistant})
    frame = pd.DataFrame(rows)
    clustered = metrics_module.cluster_interval(
        frame, "called", "resistant", np.random.default_rng(2), draws=200)
    split = frame.assign(CLUSTER=[str(index) for index in range(len(frame))])
    independent = metrics_module.cluster_interval(
        split, "called", "resistant", np.random.default_rng(2), draws=200)
    clustered_width = clustered["ppv"][1] - clustered["ppv"][0]
    independent_width = independent["ppv"][1] - independent["ppv"][0]
    assert clustered_width > 1.5 * independent_width


def test_only_a_rung_on_every_design_is_a_usable_cut_off():
    """0.008 is tested on UKMYC6 and not on UKMYC5, and 2.0 the other way round.
    Neither can serve as a cut-off, because a censored reading on the design
    that does not reach it sits on neither side."""
    series = {("UKMYC6", "BDQ"): [0.008, 0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1.0],
              ("UKMYC5", "BDQ"): [0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1.0, 2.0]}
    shared = metrics_module.shared_thresholds(series, "BDQ")
    assert shared == [0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1.0]
    assert 0.008 not in shared
    assert 2.0 not in shared


def test_a_single_design_contributes_its_whole_ladder():
    series = {("UKMYC5", "PAS"): [0.5, 1.0, 2.0, 4.0]}
    assert metrics_module.shared_thresholds(series, "PAS") == [0.5, 1.0, 2.0, 4.0]
    assert metrics_module.shared_thresholds(series, "BDQ") == []


def test_rounded_concentrations_still_match():
    """CRyPTIC label concentrations as rounded values, so the same rung can read
    0.12 on one design and 0.125 on another."""
    series = {("A", "BDQ"): [0.06, 0.12, 0.25], ("B", "BDQ"): [0.0625, 0.125, 0.25]}
    assert metrics_module.shared_thresholds(series, "BDQ") == [0.06, 0.12, 0.25]


def test_an_interval_is_placed_either_side_of_a_cut_off():
    """A reported MIC of 0.5 is the interval (0.25, 0.5], which lies above a
    cut-off of 0.25 and at or below one of 0.5."""
    import numpy as np

    lower = np.array([np.log2(0.25), np.log2(0.25)])
    upper = np.array([np.log2(0.5), np.log2(0.5)])
    for cut, expected in ((np.log2(0.25), True), (np.log2(0.5), False)):
        above, determinate = metrics_module.above_threshold(lower, upper, cut)
        assert determinate.all()
        assert (above == expected).all()


def test_a_censored_reading_the_plate_cannot_place_is_marked_indeterminate():
    import numpy as np

    # left-censored at 0.015, right-censored at 1.0
    lower = np.array([-np.inf, np.log2(1.0)])
    upper = np.array([np.log2(0.015), np.inf])
    inside = metrics_module.above_threshold(lower, upper, np.log2(0.06))
    assert list(inside[0]) == [False, True]
    assert inside[1].all()

    below_the_floor = metrics_module.above_threshold(lower, upper, np.log2(0.008))
    assert below_the_floor[1][0] is np.False_ or not below_the_floor[1][0]
    above_the_ceiling = metrics_module.above_threshold(lower, upper, np.log2(2.0))
    assert not above_the_ceiling[1][1]


def test_raising_the_cut_off_cannot_raise_the_count_above_it(data_dir):
    """The isolates above a cut-off are nested as the cut-off rises, so the
    count is monotone. A sweep that is not monotone has placed something on the
    wrong side."""
    from mic_model import concentration_series, place_mics

    mutations = cohort.load_mutations()
    frame = add_clusters(cohort.assemble(cohort.build_status(mutations)), mutations)
    series, _, _ = concentration_series()
    for drug in cohort.DRUGS:
        lower, upper, _, off = place_mics(
            frame[f"MIC_{drug}"], frame[f"PLATEDESIGN_{drug}"], series, drug)
        assert not off
        frame[f"lower_{drug}"] = lower
        frame[f"upper_{drug}"] = upper
    flags = metrics_module.rule_flags(mutations, frame.index)

    for drug in cohort.DRUGS:
        thresholds = metrics_module.shared_thresholds(series, drug)
        assert len(thresholds) >= 3
        rows = pd.DataFrame(metrics_module.sweep(
            frame, flags, drug, thresholds, draws=5))
        counts = rows[rows.rule.eq("major variant")].sort_values("threshold_mg_L")
        assert counts.above_threshold.is_monotonic_decreasing
        assert (counts.indeterminate == 0).all()


def test_the_ecoff_row_of_the_sweep_is_the_headline_row(data_dir):
    """Two tables report the same quantity, so they must agree on it, cells and
    interval alike."""
    from mic_model import concentration_series, place_mics

    mutations = cohort.load_mutations()
    frame = add_clusters(cohort.assemble(cohort.build_status(mutations)), mutations)
    series, _, _ = concentration_series()
    for drug in cohort.DRUGS:
        lower, upper, _, _ = place_mics(
            frame[f"MIC_{drug}"], frame[f"PLATEDESIGN_{drug}"], series, drug)
        frame[f"lower_{drug}"] = lower
        frame[f"upper_{drug}"] = upper
    flags = metrics_module.rule_flags(mutations, frame.index)

    for drug in cohort.DRUGS:
        headline = pd.DataFrame(metrics_module.evaluate(frame, flags, drug))
        swept = pd.DataFrame(metrics_module.sweep(
            frame, flags, drug, [metrics_module.ECOFF[drug]]))
        assert swept.is_ecoff.all()
        for rule in metrics_module.RULES:
            a = headline[headline.rule.eq(rule)].iloc[0]
            b = swept[swept.rule.eq(rule)].iloc[0]
            for field in ("true_positive", "false_positive", "false_negative",
                          "true_negative", "sensitivity", "specificity", "ppv",
                          "sensitivity_low", "sensitivity_high", "ppv_low"):
                assert a[field] == b[field], (drug, rule, field)


def test_each_interval_draws_from_its_own_stream(data_dir):
    first = metrics_module.stream("BDQ major variant 0.25").integers(0, 1 << 30, 4)
    again = metrics_module.stream("BDQ major variant 0.25").integers(0, 1 << 30, 4)
    other = metrics_module.stream("BDQ major variant 0.5").integers(0, 1 << 30, 4)
    assert list(first) == list(again)
    assert list(first) != list(other)


def test_a_cut_off_outside_the_shared_rungs_drops_the_readings_it_cannot_place(data_dir):
    """0.008 is a UKMYC6 rung and below the UKMYC5 floor, so a UKMYC5 reading of
    <=0.015 sits on neither side of it. The sweep never uses such a cut-off, and
    the filter that would drop those readings has to work if one is passed."""
    from mic_model import concentration_series, place_mics

    mutations = cohort.load_mutations()
    frame = add_clusters(cohort.assemble(cohort.build_status(mutations)), mutations)
    series, _, _ = concentration_series()
    lower, upper, _, _ = place_mics(
        frame.MIC_BDQ, frame.PLATEDESIGN_BDQ, series, "BDQ")
    frame["lower_BDQ"] = lower
    frame["upper_BDQ"] = upper
    flags = metrics_module.rule_flags(mutations, frame.index)
    with_mic = int(frame.MIC_BDQ.notna().sum())

    shared = pd.DataFrame(metrics_module.sweep(
        frame, flags, "BDQ", [0.015], draws=5))
    assert (shared.indeterminate == 0).all()
    assert (shared.isolates == with_mic).all()

    outside = pd.DataFrame(metrics_module.sweep(
        frame, flags, "BDQ", [0.008], draws=5))
    assert (outside.indeterminate > 0).all()
    assert (outside.isolates == with_mic - outside.indeterminate).all()
    assert (outside.isolates + outside.indeterminate == with_mic).all()
    cells = outside[["true_positive", "false_positive",
                     "false_negative", "true_negative"]].sum(axis=1)
    assert (cells == outside.isolates).all()
