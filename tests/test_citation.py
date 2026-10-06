"""
Tests holding the citation metadata to the repository.

CITATION.cff, for the software, and docs/ATLAS_CITATION.cff, for the atlas
table, are what a reader copies into a reference list, so their author, ORCID,
affiliation, license, repository and the dataset they cite must say what the
rest of the repository says. The repository pins no YAML parser, so the file's key
lines are read directly; the file was validated against the CFF 1.2.0 schema
separately.
"""

import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
CFF = ROOT / "CITATION.cff"
ATLAS = ROOT / "docs" / "ATLAS_CITATION.cff"
ORCID = "0009-0007-6265-6493"


def values(key, text=None):
    """Every value written for a key, at any indentation, quotes removed."""
    text = CFF.read_text() if text is None else text
    return [value.strip().strip('"') for value in
            re.findall(rf"^\s*(?:- )?{re.escape(key)}:\s*(.+)$", text, re.M)]


def first_author(path=CFF):
    block = path.read_text().split("\nauthors:\n")[1].split("\n  - ")[0]
    return {key: values(key, block)[0]
            for key in ("family-names", "given-names", "orcid", "affiliation")}


def test_the_file_declares_cff_1_2_0_and_its_required_keys():
    text = CFF.read_text()
    assert values("cff-version") == ["1.2.0"]
    for key in ("message", "title", "authors"):
        assert re.search(rf"^{key}:", text, re.M), key
    assert values("type")[0] == "software"


def test_the_author_carries_the_identity_the_repository_states():
    author = first_author()
    assert (author["given-names"], author["family-names"]) == ("Amina", "Baktiyarova")
    assert author["orcid"] == f"https://orcid.org/{ORCID}"
    assert author["affiliation"] == "Independent Researcher"
    readme = (ROOT / "README.md").read_text()
    assert f"Amina Baktiyarova, Independent Researcher\nORCID {ORCID}" in readme


def test_no_affiliation_other_than_independent_researcher_is_stated():
    assert set(values("affiliation")) == {"Independent Researcher"}


def test_the_license_and_repository_are_the_packages():
    tomllib = pytest.importorskip("tomllib")
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    assert values("license")[0] == project["license"] == "MIT"
    assert (ROOT / "LICENSE").read_text().startswith("MIT License\n")
    assert values("repository-code") == [project["urls"]["Repository"]]
    assert project["authors"][0]["name"] == "Amina Baktiyarova"


def test_the_cited_dataset_is_the_recorded_release():
    provenance = (ROOT / "docs" / "PROVENANCE.md").read_text()
    version = re.search(r"^Version: (v[\d.]+)$", provenance, re.M).group(1)
    doi = re.search(r"^Version DOI \(cite this in every output\): (\S+)$",
                    provenance, re.M).group(1)
    reference = CFF.read_text().split("\nreferences:\n")[1]
    assert values("version", reference) == [version]
    assert values("doi", reference) == [doi]
    assert values("license", reference) == ["CC-BY-4.0"]
    assert f"version {version}, Zenodo version DOI\n{doi}" in (
        ROOT / "README.md").read_text()


def test_the_cited_dataset_is_typed_in_the_reference_vocabulary():
    """
    CFF 1.2.0 types the work a file describes as software or dataset, and a
    reference from a separate vocabulary in which a dataset is data.
    """
    reference = CFF.read_text().split("\nreferences:\n")[1]
    assert values("type", reference) == ["data"]


def test_the_entry_claims_no_release():
    """No release has been made, so no version or release date is cited."""
    software = CFF.read_text().split("\nreferences:\n")[0]
    assert values("version", software) == []
    assert values("date-released") == []


# The atlas table


def atlas_section(path):
    text = path.read_text()
    return text.split("## outputs/atlas_evidence.csv")[1].split("\n---")[0]


def test_the_atlas_entry_declares_cff_1_2_0_as_a_dataset():
    text = ATLAS.read_text()
    assert values("cff-version", text) == ["1.2.0"]
    for key in ("message", "title", "authors"):
        assert re.search(rf"^{key}:", text, re.M), key
    assert values("type", text.split("\nreferences:\n")[0]) == ["dataset"]


def test_the_atlas_entry_carries_the_same_author_identity():
    assert first_author(ATLAS) == first_author(CFF)
    assert set(values("affiliation", ATLAS.read_text())) == {"Independent Researcher"}


def test_the_atlas_entry_is_licensed_as_the_derived_data_are():
    assert values("license", ATLAS.read_text().split("\nreferences:\n")[0]) == [
        "CC-BY-4.0"]
    assert "derived data under CC BY 4.0" in (ROOT / "README.md").read_text()
    licenses = " ".join((ROOT / "LICENSES.md").read_text().split())
    assert "Derived data tables released by this project are licensed under" in licenses
    assert "Attribution 4.0 International (CC BY 4.0)" in licenses


def test_the_atlas_entry_names_a_table_the_pipeline_writes_and_documents():
    title = values("title", ATLAS.read_text())[0]
    assert "(outputs/atlas_evidence.csv)" in title
    readme = (ROOT / "README.md").read_text()
    assert "| `atlas_evidence.csv` | `code/build_atlas.py` |" in readme
    assert "Written by `code/build_atlas.py`. One row per distinct mutation." in (
        atlas_section(ROOT / "docs" / "DATA_DICTIONARY.md"))


def test_the_atlas_abstract_says_what_the_code_and_dictionary_say():
    """
    The genes named are the ones the cohort definition reads, the row is the
    dictionary's, and the prevalence caveat the dictionary carries is kept.
    """
    import cohort

    abstract = " ".join(ATLAS.read_text().split("\nabstract: >-\n")[1]
                        .split("\ntype:")[0].split())
    genes = cohort.BDQ_GENES
    assert f"One row per distinct mutation in {genes[0]}, {genes[1]} and {genes[2]}" \
        in abstract
    assert cohort.MODIFIER_GENE not in abstract
    dictionary = " ".join(atlas_section(ROOT / "docs" / "DATA_DICTIONARY.md").split())
    caveat = ("Resistance in this collection is concentrated at one site by design, "
              "so no count in this table is a prevalence estimate.")
    assert caveat in dictionary
    assert caveat.replace("this table", "the table") in abstract
    assert "docs/DATA_DICTIONARY.md defines every column." in abstract


def test_the_atlas_entry_cites_the_source_release_and_the_software():
    references = ATLAS.read_text().split("\nreferences:\n")[1]
    software, source = CFF.read_text().split("\nreferences:\n")
    data, cited = references.split("\n  - type: software\n")
    assert values("type", data)[0] == "data"
    assert values("doi", data) == values("doi", source)
    assert values("version", data) == values("version", source)
    assert values("title", cited) == values("title", software)
    assert values("repository-code", cited) == values("repository-code", software)
    assert values("license", cited) == values("license", software)
    assert first_author(ATLAS) == {key: values(key, cited)[0] for key in (
        "family-names", "given-names", "orcid", "affiliation")}


def test_the_atlas_entry_claims_no_deposit():
    """No identifier, version or release date exists until the table is deposited."""
    entry = ATLAS.read_text().split("\nreferences:\n")[0]
    for key in ("doi", "identifiers", "version", "date-released"):
        assert values(key, entry) == [], key


def test_the_readme_and_the_software_entry_point_at_the_atlas_entry():
    assert "docs/ATLAS_CITATION.cff" in values("message")[0]
    section = (ROOT / "README.md").read_text().split("\n## Citing\n")[1]
    section = " ".join(section.split("\n## ")[0].split())
    assert "`CITATION.cff`" in section and "`docs/ATLAS_CITATION.cff`" in section
    assert "`outputs/atlas_evidence.csv`" in section
    reference = CFF.read_text().split("\nreferences:\n")[1]
    assert (f"v3.4.0, version DOI {values('doi', reference)[0]}" in section)
