"""
Tests for the wider gene set analysis.

Every frame here is built by hand with a planted answer, so the suite needs no
downloaded data. The estimators are checked against values worked out
independently: the q values against the published Benjamini-Hochberg example,
the stratified odds ratio against a planted effect, and the interval against
the same counts relabelled from independent isolates to a few clones.
"""

import pandas as pd
import pytest

import wider_gene_set as wider


def reference_frame(records):
    """A reference-group frame carrying the columns the module reads."""
    frame = pd.DataFrame(records).set_index("UNIQUEID")
    frame["IS_REFERENCE"] = True
    frame["MIC_BDQ"] = "0.25"
    frame["LINEAGE"] = frame.SUBLINEAGE.str.split(".").str[0]
    return frame


def mutation_frame(records):
    return pd.DataFrame(
        records, columns=["UNIQUEID", "GENE", "MUTATION", "REAL_MAJOR"])


def test_the_q_values_match_the_published_example():
    p = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216]
    expected = [0.01, 0.04, 0.084, 0.084, 0.084, 0.10, 0.10571428,
                0.216, 0.216, 0.216]
    assert wider.benjamini_hochberg(p) == pytest.approx(expected)


def test_the_q_values_come_back_in_the_order_they_were_given():
    """A q attached to the wrong p would misreport every row of the table."""
    p = [0.04, 0.001, 0.2]
    straight = wider.benjamini_hochberg(p)
    shuffled = wider.benjamini_hochberg([p[2], p[0], p[1]])
    assert shuffled == pytest.approx([straight[2], straight[0], straight[1]])


def test_the_largest_q_equals_the_largest_p():
    """The largest p is adjusted by n over n, and the running minimum cannot
    raise anything above it. This is what bounds every q at one, so no clamp
    is applied and this is the property that has to hold."""
    for values in ([0.9, 0.95, 0.99], [0.2, 0.4, 0.6, 0.8], [0.01, 0.5]):
        q = wider.benjamini_hochberg(values)
        assert max(q) == pytest.approx(max(values))
        assert max(q) <= 1.0


def test_q_values_do_not_decrease_as_p_increases():
    values = [0.001, 0.02, 0.03, 0.2, 0.5, 0.9]
    q = wider.benjamini_hochberg(values)
    assert list(q) == sorted(q), "a larger p cannot carry a smaller q"


def test_the_odds_ratio_from_four_cells():
    assert wider.odds_ratio([10, 10, 10, 10]) == pytest.approx(1.0)
    assert wider.odds_ratio([30, 10, 10, 30]) == pytest.approx(9.0)


def test_an_empty_margin_has_no_odds_ratio():
    """A ratio with no denominator is absent, not zero or one."""
    assert wider.odds_ratio([5, 0, 3, 10]) is None
    assert wider.odds_ratio([5, 3, 0, 10]) is None
    assert wider.odds_ratio([0, 5, 3, 0]) is None


def test_no_resistant_carrier_gives_an_odds_ratio_of_zero():
    """Zero is a fact about the cells: no carrier is resistant."""
    assert wider.odds_ratio([0, 20, 5, 100]) == pytest.approx(0.0)


def test_carriers_come_only_from_major_allele_rows():
    mutations = mutation_frame([
        ("s1", "fbiC", "W678G", True),
        ("s2", "fbiC", "W678G", False),
        ("s3", "lpqB", "A10V", True),
    ])
    carriers = wider.carriers_by_gene(mutations)
    assert carriers["fbiC"] == {"s1"}
    assert carriers["lpqB"] == {"s3"}


def test_carriers_of_one_variant_at_one_site_are_one_cluster():
    frame = reference_frame([
        {"UNIQUEID": "a", "SITEID": "10", "SUBLINEAGE": "lineage2.2.1"},
        {"UNIQUEID": "b", "SITEID": "10", "SUBLINEAGE": "lineage2.2.1"},
        {"UNIQUEID": "c", "SITEID": "11", "SUBLINEAGE": "lineage2.2.1"},
        {"UNIQUEID": "d", "SITEID": "10", "SUBLINEAGE": "lineage2.2.1"},
    ])
    mutations = mutation_frame([
        ("a", "fbiC", "W678G", True),
        ("b", "fbiC", "W678G", True),
        ("c", "fbiC", "W678G", True),
        ("d", "fbiC", "S197R", True),
    ])
    clusters = wider.gene_clusters(frame, mutations, "fbiC")
    assert clusters["a"] == clusters["b"], "same site, sublineage and variant"
    assert clusters["a"] != clusters["c"], "a different site is a different event"
    assert clusters["a"] != clusters["d"], "a different variant is a different event"


def test_every_non_carrier_is_its_own_cluster():
    """Isolates carrying nothing share no mutation event, so collapsing them
    would shrink the comparison group for no reason."""
    frame = reference_frame([
        {"UNIQUEID": f"s{i}", "SITEID": "10", "SUBLINEAGE": "lineage2.2.1"}
        for i in range(5)])
    clusters = wider.gene_clusters(frame, mutation_frame([]), "fbiC")
    assert clusters.nunique() == 5


def planted_strata(odds_by_stratum, baselines, carriers=200, others=200):
    """A frame with a known odds ratio in each stratum."""
    records, rows = [], 0
    for stratum, (ratio, baseline) in enumerate(zip(odds_by_stratum, baselines)):
        odds = baseline / (1 - baseline)
        risk = (odds * ratio) / (1 + odds * ratio)
        for carrier, count, chance in ((True, carriers, risk),
                                       (False, others, baseline)):
            resistant = int(round(count * chance))
            for index in range(count):
                records.append({"UNIQUEID": f"s{rows + index}",
                                "SITEID": f"{stratum:02d}",
                                "carrier": carrier,
                                "resistant_BDQ": index < resistant})
            rows += count
    frame = pd.DataFrame(records).set_index("UNIQUEID")
    return frame


def test_the_stratified_estimate_recovers_a_planted_odds_ratio():
    frame = planted_strata([5.0, 5.0, 5.0], [0.05, 0.15, 0.30])
    result = wider.stratified(frame, "carrier", "resistant_BDQ", "SITEID")
    assert result["strata"] == 3
    assert result["or"] == pytest.approx(5.0, rel=0.15)
    assert result["low"] < 5.0 < result["high"]
    assert result["homogeneity_p"] > 0.05, "a common effect is homogeneous"


def test_the_homogeneity_test_fires_when_the_strata_disagree():
    frame = planted_strata([1.0, 12.0, 1.0], [0.10, 0.10, 0.10])
    result = wider.stratified(frame, "carrier", "resistant_BDQ", "SITEID")
    assert result["homogeneity_p"] < 0.05


def test_a_stratum_with_an_empty_margin_is_dropped():
    """A site where nobody carries the variant says nothing about the
    association and must not contribute a zero margin."""
    frame = planted_strata([5.0, 5.0], [0.10, 0.10])
    frame.loc[frame.SITEID.eq("01"), "carrier"] = False
    result = wider.stratified(frame, "carrier", "resistant_BDQ", "SITEID")
    assert result["strata"] == 1
    assert result["or"] is None


def clustered_frame(cluster_count):
    """The same counts spread over a given number of clusters.

    The clusters are contiguous blocks, so a clone is mostly resistant or
    mostly susceptible rather than a sample of the whole. Spreading the
    resistant isolates evenly across clusters would make them interchangeable
    and hide the effect the test is looking for.
    """
    records = []
    block = 200 // cluster_count
    for index in range(400):
        carrier = index < 200
        resistant = (index < 60) if carrier else (index < 240)
        records.append({
            "UNIQUEID": f"s{index}",
            "carrier": carrier,
            "resistant_BDQ": resistant,
            "CLUSTER": f"c{index // block}" if carrier else f"own{index}",
        })
    return pd.DataFrame(records).set_index("UNIQUEID")


def test_the_interval_widens_when_the_carriers_are_clonal():
    """Two hundred isolates in four clones are four events, not two hundred,
    and the interval has to say so."""
    independent = wider.resampled_interval(
        clustered_frame(200), "carrier", "resistant_BDQ",
        wider.stream("independent"), draws=200)
    clonal = wider.resampled_interval(
        clustered_frame(4), "carrier", "resistant_BDQ",
        wider.stream("clonal"), draws=200)
    assert independent is not None and clonal is not None
    assert (clonal[1] - clonal[0]) > 1.5 * (independent[1] - independent[0])


def test_an_interval_is_reproducible_from_its_own_name():
    frame = clustered_frame(20)
    first = wider.resampled_interval(frame, "carrier", "resistant_BDQ",
                                     wider.stream("fbiC BDQ"), draws=100)
    second = wider.resampled_interval(frame, "carrier", "resistant_BDQ",
                                      wider.stream("fbiC BDQ"), draws=100)
    assert first == second


def test_one_gene_against_one_drug_end_to_end():
    """A planted association on a frame small enough to count by hand."""
    records = []
    for index in range(100):
        records.append({"UNIQUEID": f"s{index}", "SITEID": "10",
                        "SUBLINEAGE": f"lineage4.{index % 5}"})
    frame = reference_frame(records)
    frame["resistant_BDQ"] = [index < 20 for index in range(100)]
    # Every resistant isolate carries the variant and no susceptible one does.
    mutations = mutation_frame([(f"s{index}", "fbiC", f"V{index}A", True)
                                for index in range(20)])
    carriers = wider.carriers_by_gene(mutations)
    result = wider.test_one(frame, mutations, carriers, "fbiC", "BDQ")
    assert result["isolates"] == 100
    assert result["resistant"] == 20
    assert result["carriers_resistant"] == 20
    assert result["carriers_susceptible"] == 0
    assert result["odds_ratio"] is None, "an empty cell leaves no finite estimate"
    assert result["p_value"] < 1e-10
    assert result["resistant_carrier_groups"] == 5


def test_no_resistant_carrier_anywhere_has_no_stratified_estimate():
    """With no resistant carrier in any stratum the pooled odds ratio is zero,
    whose logarithm and variance are undefined. Reporting it as absent is the
    honest answer, and it keeps the estimator from warning on every run."""
    frame = planted_strata([5.0, 5.0], [0.10, 0.10])
    frame.loc[frame.carrier, "resistant_BDQ"] = False
    result = wider.stratified(frame, "carrier", "resistant_BDQ", "SITEID")
    assert result["or"] is None
    assert result["homogeneity_p"] is None
    assert result["strata"] == 2, "the strata are counted even where the estimate is not"


def test_a_stratified_estimate_is_absent_rather_than_zero():
    """An estimate of zero and an absent estimate are different things, and a
    row that reports zero would be read as a measured protective effect."""
    frame = planted_strata([5.0, 5.0], [0.10, 0.10])
    frame.loc[frame.carrier, "resistant_BDQ"] = False
    assert wider.stratified(frame, "carrier", "resistant_BDQ", "SITEID")["or"] != 0
