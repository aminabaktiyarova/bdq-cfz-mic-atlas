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
    rows = metrics_module.evaluate(frame, flags, "BDQ", np.random.default_rng(1))
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
