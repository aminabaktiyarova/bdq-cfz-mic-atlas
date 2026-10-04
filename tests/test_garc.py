"""
Tests for the GARC mutation grammar parser.

The rules being checked come from NOMENCLATURE.md in the piezo repository, which
is the authoritative definition of the grammar the CRyPTIC MUTATION column is
written in. The consequential ones are that X is a null call and Z a het call
rather than residues, and that neither is evidence a variant is present or
absent.
"""

import pandas as pd
import pytest

from garc import GENE_CODES, gene_codes_protein, parse_frame, parse_mutation


@pytest.mark.parametrize("mutation, codes_protein, kind, affects", [
    ("N4T", True, "SNP", "CDS"),
    ("T78I", True, "SNP", "CDS"),
    ("c-11a", True, "SNP", "PROM"),
    ("a-4t", True, "SNP", "PROM"),
    ("g1052t", False, "SNP", "RNA"),
    ("c344t", False, "SNP", "RNA"),
    ("148_ins_g", True, "INDEL", "CDS"),
    ("192_del_4", True, "INDEL", "CDS"),
    ("del_0.81", True, "GENE_DELETION", "GENE"),
    ("141_minorindel", True, "INDEL", "CDS"),
])
def test_every_observed_form_parses(mutation, codes_protein, kind, affects):
    result = parse_mutation(mutation, codes_protein)
    assert result["PARSED"], f"{mutation} did not parse"
    assert result["KIND"] == kind
    assert result["AFFECTS"] == affects


def test_x_is_a_null_call_not_a_residue():
    result = parse_mutation("G65X", True)
    assert result["IS_NULL_CALL"]
    assert not result["IS_REAL_VARIANT"], "a null call is not a variant"
    assert not result["IS_SYNONYMOUS"], "a null call is not a synonymous change"


def test_z_is_a_het_call_not_a_residue():
    result = parse_mutation("C46Z", True)
    assert result["IS_HET_CALL"]
    assert not result["IS_REAL_VARIANT"]
    assert not result["IS_SYNONYMOUS"]


def test_lowercase_x_and_z_behave_the_same_for_nucleotides():
    assert parse_mutation("c1079x", False)["IS_NULL_CALL"]
    assert parse_mutation("g3047z", False)["IS_HET_CALL"]


def test_same_residue_both_ends_is_synonymous():
    result = parse_mutation("A69A", True)
    assert result["IS_SYNONYMOUS"]
    assert not result["IS_REAL_VARIANT"], "a synonymous change is not a variant"


def test_substitution_is_a_real_variant():
    result = parse_mutation("N4T", True)
    assert result["IS_REAL_VARIANT"]
    assert not result["IS_SYNONYMOUS"]


def test_stop_codon_is_flagged():
    assert parse_mutation("Q10!", True)["IS_STOP"]
    assert not parse_mutation("Q10L", True)["IS_STOP"]


@pytest.mark.parametrize("mutation, frameshift", [
    ("148_ins_g", True),      # one base
    ("192_del_4", True),      # four bases
    ("300_ins_atg", False),   # three bases, in frame
    ("300_del_6", False),     # six bases, in frame
])
def test_frameshift_needs_a_length_not_divisible_by_three(mutation, frameshift):
    assert parse_mutation(mutation, True)["IS_FRAMESHIFT"] is frameshift


def test_frameshift_only_applies_where_there_is_a_reading_frame():
    assert parse_mutation("-12_ins_at", True)["IS_FRAMESHIFT"] is False, \
        "a promoter has no reading frame to shift"
    assert parse_mutation("2814_ins_g", False)["IS_FRAMESHIFT"] is False, \
        "rRNA has no reading frame to shift"


@pytest.mark.parametrize("mutation, affects, frameshift", [
    ("-8_del_at", "PROM", False),           # two bases, both upstream
    ("-3_del_cttgtgagc", "CDS", False),     # nine bases, six of them coding
    ("-3_del_cttgtgag", "CDS", True),       # eight bases, five of them coding
    ("-1_del_at", "CDS", True),             # one base upstream, one coding
    ("-1_del_atg", "CDS", True),            # three bases, two of them coding
    ("-3_ins_cttgtgag", "PROM", False),     # inserted upstream of the gene
])
def test_a_deletion_from_the_promoter_can_reach_the_gene(mutation, affects, frameshift):
    result = parse_mutation(mutation, True)
    assert result["AFFECTS"] == affects
    assert result["IS_FRAMESHIFT"] is frameshift


def test_a_deletion_reaching_rrna_shifts_no_frame():
    result = parse_mutation("-3_del_cttgtgag", False)
    assert result["AFFECTS"] == "RNA"
    assert result["IS_FRAMESHIFT"] is False, "rRNA has no reading frame to shift"


def test_the_position_of_an_upstream_deletion_is_still_its_first_base():
    result = parse_mutation("-3_del_cttgtgag", True)
    assert result["POSITION"] == -3
    assert result["INDEL_SIZE"] == 8


def _rows(pairs):
    """A mutation frame of (gene, mutation, CODES_PROTEIN) rows."""
    return pd.DataFrame([{"GENE": gene, "MUTATION": mutation,
                          "CODES_PROTEIN": codes_protein}
                         for gene, mutation, codes_protein in pairs])


def test_the_coding_flag_is_read_over_the_gene_not_the_row():
    frame = _rows([("Rv0678", "N4T", True),
                   ("Rv0678", "c-11a", False),
                   ("rrl", "g2814c", False)])
    coding = gene_codes_protein(frame)
    assert list(coding) == [True, True, False], \
        "CODES_PROTEIN is False on every promoter row, so the gene is only "\
        "coding if some mutation in it is placed in the coding sequence"


def test_a_missing_coding_flag_is_not_evidence_of_coding():
    frame = _rows([("Rv0678", "N4T", None)])
    assert list(gene_codes_protein(frame)) == [False]


def test_a_deletion_reaching_the_gene_is_parsed_against_the_gene():
    frame = _rows([("Rv0678", "N4T", True), ("Rv0678", "-3_del_cttgtgag", False)])
    frame[GENE_CODES] = gene_codes_protein(frame)
    parsed = parse_frame(frame)
    deletion = parsed[parsed.MUTATION.eq("-3_del_cttgtgag")].iloc[0]
    assert deletion.AFFECTS == "CDS"
    assert deletion.IS_FRAMESHIFT


def test_a_deletion_reaching_rrna_is_parsed_against_that_gene():
    frame = _rows([("rrl", "g2814c", False), ("rrl", "-3_del_cttgtgag", False)])
    frame[GENE_CODES] = gene_codes_protein(frame)
    parsed = parse_frame(frame)
    deletion = parsed[parsed.MUTATION.eq("-3_del_cttgtgag")].iloc[0]
    assert deletion.AFFECTS == "RNA"
    assert deletion.IS_FRAMESHIFT is not None, "an unknown length is not a known one"
    assert not deletion.IS_FRAMESHIFT


def test_rrna_change_is_never_synonymous():
    result = parse_mutation("g1052g", False)
    assert not result["IS_SYNONYMOUS"], "rRNA has no codon for a change to be silent in"


def test_gene_deletion_carries_the_deleted_fraction():
    assert parse_mutation("del_0.81", True)["DELETED_FRACTION"] == pytest.approx(0.81)
    assert parse_mutation("del_1.0", True)["DELETED_FRACTION"] == pytest.approx(1.0)


def test_minor_indel_length_is_unknown_not_zero():
    result = parse_mutation("141_minorindel", True)
    assert result["SIZE_RESOLVED"] is False
    assert result["IS_FRAMESHIFT"] is None, \
        "an unknown length must not be recorded as not a frameshift"


def test_unrecognised_string_is_reported_rather_than_guessed():
    result = parse_mutation("something_unexpected", True)
    assert result["PARSED"] is False
    assert result["IS_REAL_VARIANT"] is False
