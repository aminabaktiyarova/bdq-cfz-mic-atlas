"""
Tests for the unexplained-resistance tiers.

Two kinds of check. The first builds a cohort frame in which every isolate's
tier is fixed by construction, so the counts have a known answer. The second
runs the tiers through the synthetic dataset and the real cohort definitions in
code/cohort.py, so the module is checked against the definitions it depends on
rather than against a hand-made frame alone.
"""

import pandas as pd
import pytest

import cohort
import unexplained

# Tier counts planted in the frame built by planted_cohort(), counting only the
# resistant isolates.
PLANTED = {
    "attributable": 5,
    "carrier, no demonstrated effect": 3,
    "unexplained": 7,
    "indeterminate": 2,
}

# Site and sublineage combinations spanned by the planted unexplained isolates.
PLANTED_EVENTS = {"10 | lineage2.2.1": 4, "06 | lineage4.10": 2, "02 | lineage3.1": 1}

# Isolates marked LOW quality, one in each tier.
PLANTED_LOW_QUALITY = 1


def planted_cohort():
    """A cohort frame whose tier for every isolate is fixed by construction."""
    rows, mutations = [], []

    def add(unique_id, group, labels, resistant,
            quality="HIGH", site="10", sublineage="lineage2.2.1"):
        rows.append({
            "UNIQUEID": unique_id, "GROUP": group, "resistant_BDQ": resistant,
            "PHENOTYPE_QUALITY_BDQ": quality, "SITEID": site,
            "SUBLINEAGE": sublineage,
        })
        for gene, variant_class in labels:
            mutations.append({"UNIQUEID": unique_id, "GENE": gene,
                              "CLASS": variant_class, "REAL_MAJOR": True})

    add("a1", "Rv0678 frameshift", [("Rv0678", "frameshift")], True)
    add("a2", "Rv0678 substitution", [("Rv0678", "substitution")], True)
    add("a3", "pepQ solo", [("pepQ", "substitution")], True)
    add("a4", "Rv0678 stop codon", [("Rv0678", "stop codon")], True, quality="LOW")
    add("a5", "multiple variants",
        [("Rv0678", "substitution"), ("atpE", "substitution")], True)
    add("a6", "Rv0678 frameshift", [("Rv0678", "frameshift")], False)

    add("c1", "Rv0678 promoter", [("Rv0678", "promoter")], True)
    add("c2", "atpE solo", [("atpE", "substitution")], True)
    add("c3", "multiple variants",
        [("Rv0678", "promoter"), ("atpE", "substitution")], True, quality="LOW")
    add("c4", "Rv0678 in-frame indel", [("Rv0678", "in-frame indel")], False)

    for index in range(4):
        add(f"u{index}", "reference", [], True)
    for index in range(2):
        add(f"v{index}", "reference", [], True, site="06", sublineage="lineage4.10")
    add("w0", "reference", [], True, quality="LOW", site="02", sublineage="lineage3.1")
    add("s1", "reference", [], False)
    add("s2", "reference", [], False, site="06", sublineage="lineage4.10")

    add("i1", "uncertain", [], True)
    add("i2", "uncertain", [("Rv0678", "substitution")], True, quality="LOW")
    add("i3", "uncertain", [], False)

    frame = pd.DataFrame(rows).set_index("UNIQUEID")
    frame["TIER"] = unexplained.assign_tiers(frame, labels_of(mutations))
    return frame


def labels_of(mutations):
    return unexplained.carried_labels(
        pd.DataFrame(mutations, columns=["UNIQUEID", "GENE", "CLASS", "REAL_MAJOR"]))


def test_planted_tier_counts_are_recovered():
    table = unexplained.tier_counts(planted_cohort(), "BDQ")
    assert dict(zip(table.tier, table.isolates)) == PLANTED
    assert int(table.resistant_total.iloc[0]) == sum(PLANTED.values())


def test_susceptible_carriers_are_not_counted():
    """The tiers describe the resistant isolates. A susceptible carrier of the
    same variant belongs to no tier count."""
    frame = planted_cohort()
    assert frame.TIER.value_counts().sum() == len(frame)
    table = unexplained.tier_counts(frame, "BDQ")
    assert table.isolates.sum() == int(frame.resistant_BDQ.sum())
    assert table.isolates.sum() < len(frame)


def test_high_quality_restriction_drops_only_the_low_quality_isolates():
    high = unexplained.tier_counts(planted_cohort(), "BDQ", high_quality=True)
    expected = {tier: count - PLANTED_LOW_QUALITY for tier, count in PLANTED.items()}
    assert dict(zip(high.tier, high.isolates)) == expected


def test_an_uncertain_sample_is_never_counted_as_unexplained():
    """A null, het or minor call is not evidence that no variant is present, so
    such a sample cannot be counted as carrying no genotype."""
    frame = planted_cohort()
    uncertain = frame[frame.GROUP.eq("uncertain")]
    assert len(uncertain) == 3
    assert (uncertain.TIER == "indeterminate").all()


def test_a_class_with_no_demonstrated_effect_is_not_attributable():
    frame = planted_cohort()
    assert frame.loc["c1", "TIER"] == "carrier, no demonstrated effect"
    assert frame.loc["c2", "TIER"] == "carrier, no demonstrated effect"
    assert frame.loc["c3", "TIER"] == "carrier, no demonstrated effect"


def test_one_effect_class_makes_a_multiple_variant_sample_attributable():
    """a5 carries an Rv0678 substitution beside an atpE variant."""
    assert planted_cohort().loc["a5", "TIER"] == "attributable"


def test_carried_labels_separate_the_rv0678_classes():
    labels = labels_of([
        {"UNIQUEID": "x", "GENE": "Rv0678", "CLASS": "promoter", "REAL_MAJOR": True},
        {"UNIQUEID": "x", "GENE": "pepQ", "CLASS": "substitution", "REAL_MAJOR": True},
        {"UNIQUEID": "y", "GENE": "Rv0678", "CLASS": "frameshift", "REAL_MAJOR": True},
        {"UNIQUEID": "z", "GENE": "Rv0678", "CLASS": "substitution", "REAL_MAJOR": False},
    ])
    assert labels["x"] == frozenset({"Rv0678 promoter", "pepQ"})
    assert labels["y"] == frozenset({"Rv0678 frameshift"})
    assert "z" not in labels.index


def test_the_guard_fires_when_a_reference_sample_carries_a_variant():
    frame = pd.DataFrame([{"UNIQUEID": "x", "GROUP": "reference"}]).set_index("UNIQUEID")
    labels = labels_of([{"UNIQUEID": "x", "GENE": "Rv0678",
                         "CLASS": "frameshift", "REAL_MAJOR": True}])
    with pytest.raises(ValueError):
        unexplained.assign_tiers(frame, labels)


def test_the_guard_fires_when_a_carrier_group_has_no_variant():
    frame = pd.DataFrame(
        [{"UNIQUEID": "x", "GROUP": "Rv0678 substitution"}]).set_index("UNIQUEID")
    with pytest.raises(ValueError):
        unexplained.assign_tiers(frame, labels_of([]))


def test_wilson_intervals_bracket_their_own_fractions():
    for row in unexplained.tier_counts(planted_cohort(), "BDQ").itertuples():
        assert row.wilson_low <= row.percent <= row.wilson_high


def test_coarse_events_count_site_and_sublineage_combinations():
    isolates, counts = unexplained.coarse_events(planted_cohort(), "BDQ")
    assert isolates == PLANTED["unexplained"]
    assert dict(counts) == PLANTED_EVENTS


def test_tiers_partition_the_synthetic_cohort(data_dir):
    """Run against the cohort definitions rather than a hand-made frame."""
    mutations = cohort.load_mutations()
    frame = cohort.assemble(cohort.build_status(mutations))
    tiers = unexplained.assign_tiers(frame, unexplained.carried_labels(mutations))
    assert tiers.notna().all()
    assert set(tiers.unique()) <= set(unexplained.TIER_ORDER)
    assert int((tiers == "unexplained").sum()) == int(frame.GROUP.eq("reference").sum())
    assert int((tiers == "indeterminate").sum()) == int(frame.GROUP.eq("uncertain").sum())


def test_mmpl5_variants_do_not_make_a_sample_attributable(data_dir):
    """Every sample in the fixture carries an mmpL5 variant, as almost every
    real sample does. None of them becomes a carrier through mmpL5 alone."""
    mutations = cohort.load_mutations()
    frame = cohort.assemble(cohort.build_status(mutations))
    tiers = unexplained.assign_tiers(frame, unexplained.carried_labels(mutations))
    reference = frame.index[frame.GROUP.eq("reference")]
    assert len(reference) > 0
    assert (tiers.loc[reference] == "unexplained").all()


def planted_structural_case():
    """A cohort where the gene set decides the answer.

    p1 carries a variant in the modifier gene alone, p2 is uncallable in Rv0678
    while carrying a pepQ variant, p3 carries an Rv0678 variant, and p4 carries
    one but is susceptible.
    """
    frame = pd.DataFrame([
        {"UNIQUEID": "p1", "resistant_BDQ": True, "PHENOTYPE_QUALITY_BDQ": "HIGH"},
        {"UNIQUEID": "p2", "resistant_BDQ": True, "PHENOTYPE_QUALITY_BDQ": "HIGH"},
        {"UNIQUEID": "p3", "resistant_BDQ": True, "PHENOTYPE_QUALITY_BDQ": "LOW"},
        {"UNIQUEID": "p4", "resistant_BDQ": False, "PHENOTYPE_QUALITY_BDQ": "HIGH"},
    ]).set_index("UNIQUEID")
    mutations = pd.DataFrame([
        {"UNIQUEID": "p1", "GENE": "mmpL5", "REAL_MAJOR": True, "UNCERTAIN": False},
        {"UNIQUEID": "p2", "GENE": "Rv0678", "REAL_MAJOR": False, "UNCERTAIN": True},
        {"UNIQUEID": "p2", "GENE": "pepQ", "REAL_MAJOR": True, "UNCERTAIN": False},
        {"UNIQUEID": "p3", "GENE": "Rv0678", "REAL_MAJOR": True, "UNCERTAIN": False},
        {"UNIQUEID": "p4", "GENE": "Rv0678", "REAL_MAJOR": True, "UNCERTAIN": False},
    ])
    return frame, mutations


def test_structural_counts_partition_the_resistant_isolates():
    frame, mutations = planted_structural_case()
    counts = unexplained.structural_counts(frame, mutations, cohort.BDQ_GENES, "BDQ")
    assert counts["resistant_total"] == 3
    assert counts["unexplained"] + counts["indeterminate"] + counts["carrier"] == 3
    assert (counts["unexplained"], counts["indeterminate"], counts["carrier"]) == (1, 1, 1)


def test_admitting_the_modifier_gene_removes_the_unexplained_isolate():
    frame, mutations = planted_structural_case()
    with_modifier = unexplained.structural_counts(
        frame, mutations, cohort.ALL_GENES, "BDQ")
    assert with_modifier["unexplained"] == 0
    assert with_modifier["carrier"] == 2


def test_an_uncallable_gene_outranks_a_carried_variant():
    """p2 carries a pepQ variant and is uncallable in Rv0678. It is counted as
    indeterminate, matching the precedence in assign_tiers."""
    frame, mutations = planted_structural_case()
    counts = unexplained.structural_counts(frame, mutations, cohort.BDQ_GENES, "BDQ")
    assert counts["indeterminate"] == 1
    assert counts["carrier"] == 1


def test_structural_counts_respect_the_quality_filter():
    frame, mutations = planted_structural_case()
    high = unexplained.structural_counts(
        frame, mutations, cohort.BDQ_GENES, "BDQ", high_quality=True)
    assert high["resistant_total"] == 2
    assert high["carrier"] == 0


def test_the_three_gene_set_reproduces_the_tier_counts(data_dir):
    mutations = cohort.load_mutations()
    frame = cohort.assemble(cohort.build_status(mutations))
    frame["TIER"] = unexplained.assign_tiers(
        frame, unexplained.carried_labels(mutations))
    for drug in cohort.DRUGS:
        tiers = unexplained.tier_counts(frame, drug).set_index("tier").isolates
        counts = unexplained.structural_counts(
            frame, mutations, cohort.BDQ_GENES, drug)
        assert counts["unexplained"] == tiers["unexplained"]
        assert counts["indeterminate"] == tiers["indeterminate"]


def test_every_gene_set_is_reported_for_both_drugs(data_dir):
    mutations = cohort.load_mutations()
    frame = cohort.assemble(cohort.build_status(mutations))
    table = unexplained.gene_set_table(frame, mutations)
    assert len(table) == len(unexplained.GENE_SETS) * len(cohort.DRUGS)
    assert set(table.gene_set) == {label for label, _ in unexplained.GENE_SETS}


def test_coarse_events_order_ties_by_key():
    """value_counts leaves the order of equal counts to the implementation. The
    isolates below are built so that first appearance and key order disagree."""
    rows = []
    for unique_id, site, sublineage in [
        ("t1", "10", "b"), ("t2", "02", "a"),
        ("t3", "06", "c"), ("t4", "06", "c"),
    ]:
        rows.append({"UNIQUEID": unique_id, "GROUP": "reference",
                     "resistant_BDQ": True, "PHENOTYPE_QUALITY_BDQ": "HIGH",
                     "SITEID": site, "SUBLINEAGE": sublineage})
    frame = pd.DataFrame(rows).set_index("UNIQUEID")
    frame["TIER"] = unexplained.assign_tiers(frame, labels_of([]))
    isolates, counts = unexplained.coarse_events(frame, "BDQ")
    assert isolates == 4
    assert list(counts.index) == ["06 | c", "02 | a", "10 | b"]
    assert list(counts) == [2, 1, 1]
