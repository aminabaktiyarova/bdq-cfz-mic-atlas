"""
Tests for the cohort definitions.

These assert the classifications the whole analysis rests on: which samples
count as reference, which as solo, and which are excluded as uncertain. The
counts come from the synthetic dataset, where the answer is known by
construction.
"""

import cohort


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
