"""
Tests for the comparison against the WHO catalogue.

The catalogue the module reads in anger is restricted-licence and is not in
this repository, so every test here builds its own catalogue in piezo's format
with grades planted in it, and its own mutation frames. The suite therefore
runs on a machine that has never obtained a catalogue or downloaded the CRyPTIC
data.
"""

from pathlib import Path

import pandas as pd
import pytest

import benchmark_catalogue as benchmark

HEADER = ("GENBANK_REFERENCE,CATALOGUE_NAME,CATALOGUE_VERSION,CATALOGUE_GRAMMAR,"
          "PREDICTION_VALUES,DRUG,MUTATION,PREDICTION,SOURCE,EVIDENCE,OTHER")

# gene, mutation, drug, grade. geneT carries no clofazimine row, which is how a
# drug the catalogue does not associate with a gene reaches the module.
PLANTED = [
    ("geneR", "S2F", "BDQ", "R"), ("geneR", "S2F", "CFZ", "R"),
    ("geneR", "A10V", "BDQ", "U"), ("geneR", "A10V", "CFZ", "U"),
    ("geneR", "T30A", "BDQ", "S"), ("geneR", "T30A", "CFZ", "S"),
    ("geneP", "G50D", "BDQ", "F"), ("geneP", "G50D", "CFZ", "F"),
    ("geneT", "D28G", "BDQ", "R"),
    ("geneM", "I948V", "BDQ", "S"), ("geneM", "I948V", "CFZ", "S"),
    ("geneM", "100_del_gg", "BDQ", "S"), ("geneM", "100_del_gg", "CFZ", "S"),
    ("geneR", "del_0.55", "BDQ", "U"), ("geneR", "del_0.55", "CFZ", "U"),
    ("geneR", "200_del_gtca", "BDQ", "R"), ("geneR", "200_del_gtca", "CFZ", "R"),
]


@pytest.fixture()
def catalogue_file(tmp_path):
    rows = [HEADER]
    for gene, mutation, drug, grade in PLANTED:
        rows.append(f"NC_000962.3,WHO-UCN-GTB-PCI-2023.5,2.0,GARC1,RFUS,"
                    f"{drug},{gene}@{mutation},{grade},{{}},{{}},{{}}")
    path = tmp_path / "planted_GARC1_RFUS.csv"
    path.write_text("\n".join(rows) + "\n")
    return path


@pytest.fixture()
def catalogue(catalogue_file):
    return benchmark.load_catalogue(path=catalogue_file)


def frame(rows):
    """A mutation frame carrying only the columns the module reads."""
    return pd.DataFrame(
        rows,
        columns=["UNIQUEID", "GENE", "MUTATION", "CLASS",
                 "IS_REAL_VARIANT", "IS_MINOR"])


def status_for(index, lof=()):
    """The per-sample status the module reads. The column is named for the
    real modifier gene because that is the interface, not fixture data."""
    return pd.DataFrame({"mmpL5_LOF": [uid in set(lof) for uid in index]},
                        index=pd.Index(index, name="UNIQUEID"))


def call(mutations, catalogue, epistasis=False, lof=()):
    index = pd.Index(sorted(set(mutations.UNIQUEID)), name="UNIQUEID")
    grades = benchmark.grade_variants(mutations, catalogue)
    return benchmark.isolate_calls(mutations, grades, index,
                                   status_for(index, lof), epistasis=epistasis)


def test_the_planted_grades_come_back(catalogue):
    mutations = frame([
        ("s1", "geneR", "S2F", "substitution", True, False),
        ("s2", "geneR", "A10V", "substitution", True, False),
        ("s3", "geneR", "T30A", "substitution", True, False),
    ])
    grades = benchmark.grade_variants(mutations, catalogue).set_index("MUTATION")
    assert grades.loc["S2F", "BDQ"] == "R"
    assert grades.loc["A10V", "BDQ"] == "U"
    assert grades.loc["T30A", "BDQ"] == "S"


def test_a_drug_the_catalogue_does_not_list_for_the_gene(catalogue):
    """geneT has a bedaquiline row and no clofazimine row. The absence of an
    association is not a grade, and it is not evidence of resistance."""
    mutations = frame([("s1", "geneT", "D28G", "substitution", True, False)])
    grades = benchmark.grade_variants(mutations, catalogue)
    assert grades.loc[0, "BDQ"] == "R"
    assert grades.loc[0, "CFZ"] == "not listed"
    calls = call(mutations, catalogue)
    assert calls.loc["s1", "BDQ"] == "R"
    assert calls.loc["s1", "CFZ"] == "S"


def test_the_strongest_grade_an_isolate_carries_wins(catalogue):
    mutations = frame([
        ("resistant", "geneR", "S2F", "substitution", True, False),
        ("resistant", "geneR", "T30A", "substitution", True, False),
        ("unknown", "geneR", "A10V", "substitution", True, False),
        ("unknown", "geneR", "T30A", "substitution", True, False),
        ("failed", "geneP", "G50D", "substitution", True, False),
        ("failed", "geneR", "T30A", "substitution", True, False),
        ("susceptible", "geneR", "T30A", "substitution", True, False),
    ])
    calls = call(mutations, catalogue)
    assert calls.loc["resistant", "BDQ"] == "R"
    assert calls.loc["unknown", "BDQ"] == "U"
    assert calls.loc["failed", "BDQ"] == "F"
    assert calls.loc["susceptible", "BDQ"] == "S"


def test_unknown_outranks_failed(catalogue):
    """A position the catalogue grades U and one it grades F are both short of
    a call, and U is the one reported."""
    mutations = frame([
        ("s1", "geneR", "A10V", "substitution", True, False),
        ("s1", "geneP", "G50D", "substitution", True, False),
    ])
    assert call(mutations, catalogue).loc["s1", "BDQ"] == "U"


def test_an_isolate_carrying_no_graded_variant_is_susceptible(catalogue):
    mutations = frame([("s1", "geneR", "T30A", "substitution", True, False)])
    index = pd.Index(["s1", "s2"], name="UNIQUEID")
    grades = benchmark.grade_variants(mutations, catalogue)
    calls = benchmark.isolate_calls(mutations, grades, index, status_for(index))
    assert calls.loc["s2", "BDQ"] == "S"


def test_minor_alleles_and_unreal_variants_are_not_graded(catalogue):
    """A minor allele is reported Susceptible by CRyPTIC's own application of
    this catalogue, and a synonymous change is not a variant. Neither reaches
    a grade, so neither can make an isolate resistant."""
    mutations = frame([
        ("minor", "geneR", "S2F", "substitution", True, True),
        ("synonymous", "geneR", "S2F", "substitution", False, False),
    ])
    assert not benchmark.graded_rows(mutations).any()
    calls = call(mutations, catalogue)
    assert set(calls["BDQ"]) == {"S"}


def test_the_duplicate_spelling_of_a_deletion_is_graded(catalogue):
    """One deletion is written twice in v3.4.0, as del_<fraction> and as the
    sequence removed, and the two do not grade alike. cohort.REAL_MAJOR drops
    the second so one event counts once; dropping it here would lose the grade
    the catalogue assigns to the event."""
    mutations = frame([
        ("s1", "geneR", "del_0.55", "gene deletion", True, False),
        ("s1", "geneR", "200_del_gtca", "frameshift", True, False),
    ])
    assert benchmark.graded_rows(mutations).sum() == 2
    assert call(mutations, catalogue).loc["s1", "BDQ"] == "R"


def test_the_epistasis_rule_calls_the_isolate_susceptible(catalogue):
    mutations = frame([
        ("disrupted", "geneR", "S2F", "substitution", True, False),
        ("disrupted", "geneM", "100_del_gg", "frameshift", True, False),
        ("intact", "geneR", "S2F", "substitution", True, False),
        ("intact", "geneM", "I948V", "substitution", True, False),
    ])
    plain = call(mutations, catalogue)
    assert plain.loc["disrupted", "BDQ"] == "R"
    assert plain.loc["intact", "BDQ"] == "R"

    overridden = call(mutations, catalogue, epistasis=True, lof=["disrupted"])
    assert overridden.loc["disrupted", "BDQ"] == "S"
    assert overridden.loc["disrupted", "CFZ"] == "S"
    assert overridden.loc["intact", "BDQ"] == "R"


def test_the_catalogue_identity_is_checked(tmp_path):
    rows = [HEADER]
    # piezo requires every declared prediction value to occur, so the wrong
    # name is the only thing this catalogue gets wrong.
    for mutation, grade in (("S2F", "R"), ("A10V", "U"),
                            ("T30A", "S"), ("G50D", "F")):
        rows.append("NC_000962.3,SOMETHING-ELSE,2.0,GARC1,RFUS,BDQ,"
                    f"geneR@{mutation},{grade},{{}},{{}},{{}}")
    path = tmp_path / "wrong_GARC1_RFUS.csv"
    path.write_text("\n".join(rows) + "\n")
    with pytest.raises(ValueError, match="catalogue name"):
        benchmark.load_catalogue(path=path)


def test_a_missing_catalogue_says_where_it_comes_from(tmp_path):
    with pytest.raises(FileNotFoundError, match="PROVENANCE"):
        benchmark.load_catalogue(path=tmp_path / "absent.csv")


def test_every_output_stays_inside_quarantine():
    benchmark.check_output_paths()
    for name in ("REPORT", "METRICS", "GRADES"):
        assert Path("quarantine") in getattr(benchmark, name).parents


def test_an_output_outside_quarantine_is_refused(monkeypatch):
    monkeypatch.setattr(benchmark, "REPORT", Path("outputs/benchmark_report.txt"))
    with pytest.raises(ValueError, match="outside"):
        benchmark.check_output_paths()
