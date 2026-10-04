"""
Tests for the cohort definitions.

These assert the classifications the whole analysis rests on: which samples
count as reference, which as solo, and which are excluded as uncertain. The
counts come from the synthetic dataset, where the answer is known by
construction.
"""

import pandas as pd

import cohort
from garc import GENE_CODES, gene_codes_protein, parse_frame


def test_target_gene_mutations_all_parse(data_dir):
    mutations = cohort.load_mutations()
    assert mutations.PARSED.all()


def test_synonymous_variants_do_not_make_a_sample_non_reference(data_dir, expected):
    status = cohort.build_status()
    reference = int(status.IS_REFERENCE.sum())
    assert reference == expected["n_reference"] + expected["n_synonymous_only"] \
        + expected["n_serial_patients"]


def test_null_het_and_minor_calls_are_all_excluded_as_uncertain(data_dir, expected):
    status = cohort.build_status()
    uncertain = int(status.any_uncertain.sum())
    assert uncertain == (expected["n_null_call"] + expected["n_het_call"]
                         + expected["n_minor_indel"])


def test_minor_indels_are_excluded_like_het_calls(data_dir, expected):
    """A sub-population carrying an indel is not an absence of one."""
    status = cohort.build_status()
    assert int(status.IS_REFERENCE.sum()) + int(status.any_uncertain.sum()) \
        + int(status.IS_SOLO.sum()) == len(status)
    assert not status.loc[status.index[status["uncertain_Rv0678"]], "IS_REFERENCE"].any()


def test_solo_counts_match_the_planted_groups(data_dir, expected):
    status = cohort.build_status()
    counts = status.GROUP.value_counts()
    assert counts["Rv0678 frameshift"] == expected["n_frameshift"]
    assert counts["Rv0678 substitution"] == expected["n_substitution"]
    assert counts["pepQ solo"] == expected["n_pepq"]


def test_mmpl5_is_not_part_of_the_reference_criterion(data_dir):
    """
    Every sample in the fixture carries an mmpL5 variant, as almost every real
    sample does. If mmpL5 entered the reference criterion the reference group
    would be empty.
    """
    mutations = cohort.load_mutations()
    carriers = mutations[mutations.GENE.eq("mmpL5") & mutations.REAL_MAJOR].UNIQUEID.nunique()
    status = cohort.build_status()
    assert carriers == len(status)
    assert int(status.IS_REFERENCE.sum()) > 0


def test_assemble_joins_phenotypes_and_lineage(data_dir):
    joined = cohort.assemble()
    for drug in cohort.DRUGS:
        assert f"MIC_{drug}" in joined.columns
        assert f"resistant_{drug}" in joined.columns
    assert joined.LINEAGE.notna().all()
    assert joined.index.is_unique


def test_censoring_flags_are_read_from_the_mic_string(data_dir):
    joined = cohort.assemble()
    text = joined.MIC_BDQ.astype(str)
    assert (joined.censored_left_BDQ == text.str.startswith("<=")).all()
    assert (joined.censored_right_BDQ == text.str.startswith(">")).all()


def _classified(rows, gene="Rv0678"):
    """Classify a small mutation table of (sample, mutation, IS_MINOR) rows."""
    frame = pd.DataFrame([{"UNIQUEID": sample, "GENE": gene, "MUTATION": mutation,
                           "CODES_PROTEIN": not mutation.startswith("-"),
                           "IS_MINOR": is_minor}
                          for sample, mutation, is_minor in rows])
    frame[GENE_CODES] = gene_codes_protein(frame)
    return cohort.classify(parse_frame(frame)).set_index(["UNIQUEID", "MUTATION"])


def test_a_deletion_written_twice_counts_as_one_variant():
    """CRyPTIC reports a large deletion both as a fraction of the gene and as the
    sequence removed. Section 8.9 of the project record has the thirteen samples
    in the real table."""
    one = "one.sample"
    classified = _classified([(one, "del_0.81", False), (one, "-3_del_438", False)])
    assert classified.loc[(one, "-3_del_438")].DOUBLE_REPORTED
    assert not classified.loc[(one, "-3_del_438")].REAL_MAJOR
    assert classified.loc[(one, "del_0.81")].REAL_MAJOR
    assert classified.loc[(one, "del_0.81")].CLASS == "gene deletion"
    assert int(classified.REAL_MAJOR.sum()) == 1


def test_a_smaller_deletion_beside_the_one_written_twice_is_kept():
    one = "one.sample"
    classified = _classified([(one, "del_0.81", False), (one, "-3_del_438", False),
                              (one, "300_del_g", False)])
    assert classified.loc[(one, "-3_del_438")].DOUBLE_REPORTED
    assert not classified.loc[(one, "300_del_g")].DOUBLE_REPORTED
    assert int(classified.REAL_MAJOR.sum()) == 2


def test_a_minor_deletion_is_not_taken_for_the_second_report():
    one = "one.sample"
    classified = _classified([(one, "del_0.81", False), (one, "-3_del_438", True)])
    assert not classified.loc[(one, "-3_del_438")].DOUBLE_REPORTED, \
        "a minor call is a detected sub-population, not the same event again"


def test_a_deletion_on_its_own_is_left_alone():
    classified = _classified([("one.sample", "300_del_g", False)])
    assert not classified.DOUBLE_REPORTED.any()
    assert classified.loc[("one.sample", "300_del_g")].REAL_MAJOR


def test_the_collapse_survives_arrays_pandas_will_not_let_us_write_to():
    """From pandas 3 an array handed out by a Series or an Index is read-only,
    and a mask combined in place raises there. Every mask in the collapse is
    combined into a new object instead, which this holds it to."""
    import numpy as np

    def read_only(result):
        array = np.array(result, dtype=bool)
        array.flags.writeable = False
        return array

    isin = pd.MultiIndex.isin
    to_numpy = pd.Series.to_numpy

    def frozen_isin(self, values, level=None):
        return read_only(isin(self, values, level=level))

    def frozen_to_numpy(self, *arguments, **keywords):
        array = to_numpy(self, *arguments, **keywords)
        if array.dtype == bool:
            array = read_only(array)
        return array

    one = "one.sample"
    rows = [(one, "del_0.81", False), (one, "-3_del_438", False)]
    try:
        pd.MultiIndex.isin = frozen_isin
        pd.Series.to_numpy = frozen_to_numpy
        classified = _classified(rows)
    finally:
        pd.MultiIndex.isin = isin
        pd.Series.to_numpy = to_numpy

    assert classified.loc[(one, "-3_del_438")].DOUBLE_REPORTED
    assert int(classified.REAL_MAJOR.sum()) == 1


def test_a_deletion_in_another_sample_is_left_alone():
    """The second report is the same event in the same gene of the same sample,
    so a deletion elsewhere in the table is untouched by it."""
    classified = _classified([("carrier", "del_0.81", False),
                              ("carrier", "-3_del_438", False),
                              ("other", "-3_del_438", False)])
    assert classified.loc[("carrier", "-3_del_438")].DOUBLE_REPORTED
    assert not classified.loc[("other", "-3_del_438")].DOUBLE_REPORTED
    assert classified.loc[("other", "-3_del_438")].REAL_MAJOR


def test_equally_frequent_values_are_listed_in_order_of_value():
    series = pd.Series(["d", "c", "c", "a", "a", "b"])
    counts = cohort.ranked_counts(series)
    assert list(counts.index) == ["a", "c", "b", "d"]
    assert list(counts) == [2, 2, 1, 1]


def test_a_cut_listing_takes_the_same_tied_values_whatever_the_input_order():
    values = ["x"] * 3 + ["y"] * 2 + ["z"] * 2 + ["w"] * 2
    first = cohort.ranked_counts(pd.Series(values), 3)
    second = cohort.ranked_counts(pd.Series(list(reversed(values))), 3)
    assert list(first.index) == ["x", "w", "y"]
    assert list(first.index) == list(second.index)
    assert list(first) == list(second) == [3, 2, 2]


def test_a_long_run_of_equally_frequent_values_keeps_the_order_of_value():
    # Long enough that an unstable sort reorders the tied block: numpy's
    # quicksort falls back to insertion sort on short arrays and preserves
    # their order, so a short fixture cannot detect the loss of stability.
    labels = [f"g{index:02d}" for index in range(40)]
    series = pd.Series(["top"] * 3 + [label for label in labels for _ in range(2)])
    counts = cohort.ranked_counts(series)
    assert list(counts.index) == ["top"] + labels
    assert list(counts) == [3] + [2] * len(labels)


def test_an_indel_confined_to_the_promoter_is_classed_as_a_promoter_change():
    one = "one.sample"
    classified = _classified([(one, "-21_ins_ttc", False), (one, "19_del_gtc", False),
                              (one, "c-11a", False)])
    assert classified.loc[(one, "-21_ins_ttc")].CLASS == "promoter"
    assert classified.loc[(one, "19_del_gtc")].CLASS == "in-frame indel"
    assert classified.loc[(one, "c-11a")].CLASS == "promoter"
