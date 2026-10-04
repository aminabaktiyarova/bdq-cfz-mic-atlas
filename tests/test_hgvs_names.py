"""
Tests for the GARC to HGVS translation.

The translation is a mapping between two naming conventions, so the property
that matters is that it is exact in both directions. Every residue pair in the
genetic code is round tripped, rather than a few examples being spot checked.
"""

import itertools

import pandas as pd
import pytest

import hgvs_names as hgvs


def test_every_residue_pair_survives_the_round_trip():
    """441 pairs, every residue against every residue, the stop codon
    included. A single wrong entry in the three-letter map shows up here."""
    residues = sorted(hgvs.THREE_LETTER)
    assert len(residues) == 21
    for reference, alternate in itertools.product(residues, repeat=2):
        mutation = f"{reference}42{alternate}"
        assert hgvs.round_trip_holds(mutation), mutation


def test_the_names_are_the_ones_the_conventions_specify():
    assert hgvs.to_hgvs("S2F")[0] == "p.Ser2Phe"
    assert hgvs.to_hgvs("A36V")[0] == "p.Ala36Val"
    assert hgvs.to_hgvs("R38!")[0] == "p.Arg38Ter", "the stop codon is Ter"
    assert hgvs.to_hgvs("c-11a")[0] == "c.-11C>A"
    assert hgvs.to_hgvs("g-17t")[0] == "c.-17G>T"


def test_a_nucleotide_substitution_is_uppercased_and_keeps_its_direction():
    """HGVS writes the reference first and the alternate after the arrow.
    Reversing them names a different variant."""
    name, _ = hgvs.to_hgvs("a100g")
    assert name == "c.100A>G"
    assert name != "c.100G>A"
    assert hgvs.to_garc(name) == "a100g"


def test_a_gene_that_does_not_code_protein_takes_the_other_prefix():
    assert hgvs.to_hgvs("a100g", codes_protein=False)[0] == "n.100A>G"
    name, reason = hgvs.to_hgvs("S2F", codes_protein=False)
    assert name is None and "does not code protein" in reason


def test_a_null_or_het_call_is_refused_rather_than_named():
    """Neither is a residue, so neither has a name. The reason says which
    case it is: a position that could not be read is a different thing from
    a letter outside the genetic code, and the column is read by people."""
    for mutation in ("G65X", "G65Z"):
        name, reason = hgvs.to_hgvs(mutation)
        assert name is None, mutation
        assert reason == "a null or het call is not a residue", mutation
    for mutation in ("X65G", "Z65G"):
        name, reason = hgvs.to_hgvs(mutation)
        assert name is None, mutation
        assert reason == "a null or het call is not a residue", mutation


def test_a_letter_outside_the_genetic_code_says_so_instead():
    """B, J, O and U are not residues and are not null or het calls."""
    name, reason = hgvs.to_hgvs("B65G")
    assert name is None
    assert reason == "residue outside the genetic code"


def test_a_lowercase_null_or_het_call_matches_no_form_at_all():
    """Written lowercase it is neither a protein nor a nucleotide change."""
    for mutation in ("g65x", "g65z"):
        name, reason = hgvs.to_hgvs(mutation)
        assert name is None
        assert reason == "not a form this translation covers"


def test_an_indel_carries_no_name_and_says_why():
    for mutation in ("141_ins_c", "16_del_g", "-3_del_cttgtgag"):
        name, reason = hgvs.to_hgvs(mutation)
        assert name is None
        assert "protein consequence" in reason


def test_a_gene_deletion_carries_no_name():
    for mutation in ("del_1.0", "del_0.55", "del_minorindel"):
        name, reason = hgvs.to_hgvs(mutation)
        assert name is None
        assert "no HGVS variant name" in reason


def test_a_negative_position_is_kept_on_both_sides():
    """A negative position is a base upstream of the start codon in both
    conventions, so it must survive unchanged."""
    assert hgvs.to_hgvs("c-11a")[0] == "c.-11C>A"
    assert hgvs.to_garc("c.-11C>A") == "c-11a"


def test_the_table_records_a_reason_for_every_unnamed_variant():
    table = hgvs.name_table([("Rv0678", "S2F"), ("Rv0678", "141_ins_c"),
                             ("Rv0678", "del_1.0"), ("pepQ", "c-11a")])
    assert len(table) == 4
    named = table.hgvs.notna()
    assert int(named.sum()) == 2
    assert table[~named].reason.notna().all()
    assert table[named].reason.isna().all()


def test_a_name_that_does_not_translate_back_fails_the_run(monkeypatch):
    """The guard is what makes the table trustworthy, so it has to fire."""
    monkeypatch.setattr(hgvs, "to_garc", lambda name, codes_protein=True: "wrong")
    with pytest.raises(ValueError, match="does not translate back"):
        hgvs.name_table([("Rv0678", "S2F")])


def test_the_two_maps_are_inverses_of_one_another():
    assert set(hgvs.ONE_LETTER.values()) == set(hgvs.THREE_LETTER)
    for one, three in hgvs.THREE_LETTER.items():
        assert hgvs.ONE_LETTER[three] == one
