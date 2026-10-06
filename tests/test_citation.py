"""
Tests holding the citation metadata to the repository.

CITATION.cff, for the software, and docs/DATA_CITATION.cff, for the released
data tables, are what a reader copies into a reference list, so their author,
ORCID, affiliation, license, repository and the dataset they cite must say what
the rest of the repository says. The repository pins no YAML parser, so the
file's key lines are read directly; the file was validated against the CFF
1.2.0 schema separately.
"""

import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
CFF = ROOT / "CITATION.cff"
DATA = ROOT / "docs" / "DATA_CITATION.cff"
ORCID = "0009-0007-6265-6493"
SOFTWARE_DOI = "10.5281/zenodo.23196270"
SOFTWARE_CONCEPT_DOI = "10.5281/zenodo.23196269"
DATA_DOI = "10.5281/zenodo.23194535"
DATA_CONCEPT_DOI = "10.5281/zenodo.23194534"
NUMBER_WORDS = ("zero", "one", "two", "three", "four", "five", "six", "seven",
                "eight", "nine", "ten", "eleven", "twelve", "thirteen",
                "fourteen", "fifteen", "sixteen", "seventeen", "eighteen",
                "nineteen", "twenty")


def values(key, text=None):
    """Every value written for a key, at any indentation, quotes removed."""
    text = CFF.read_text() if text is None else text
    return [value.strip().strip('"') for value in
            re.findall(rf"^\s*(?:- )?{re.escape(key)}:\s*(.+)$", text, re.M)]


def first_author(path=CFF):
    block = path.read_text().split("\nauthors:\n")[1].split("\n  - ")[0]
    return {key: values(key, block)[0]
            for key in ("family-names", "given-names", "orcid", "affiliation")}


def citing_section():
    section = (ROOT / "README.md").read_text().split("\n## Citing\n")[1]
    return " ".join(section.split("\n## ")[0].split())


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


def test_the_entry_cites_the_archived_release():
    """
    The entry carries the version DOI of the archived release, the version
    micecoff reports, the date it was archived, and the concept DOI that
    resolves to the latest release. The README carries both DOIs under the
    same two names.
    """
    import micecoff

    head = CFF.read_text().split("\nreferences:\n")[0]
    identifiers = head.split("\nidentifiers:\n")[1].split("\nauthors:\n")[0]
    assert values("doi", head) == [SOFTWARE_DOI]
    assert values("version", head) == [micecoff.__version__]
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", values("date-released", head)[0])
    assert values("type", identifiers) == ["doi"]
    assert values("value", identifiers) == [SOFTWARE_CONCEPT_DOI]
    section = citing_section()
    assert f"version DOI {SOFTWARE_DOI}" in section
    assert f"concept DOI {SOFTWARE_CONCEPT_DOI}" in section


# The released data tables


def atlas_section(path):
    text = path.read_text()
    return text.split("## outputs/atlas_evidence.csv")[1].split("\n---")[0]


def data_head():
    """Everything the entry says about the tables, before its references."""
    return DATA.read_text().split("\nreferences:\n")[0]


def data_abstract():
    return " ".join(DATA.read_text().split("\nabstract: >-\n")[1]
                    .split("\ntype:")[0].split())


def defined_tables():
    """Every table docs/DATA_DICTIONARY.md gives a column list for."""
    dictionary = (ROOT / "docs" / "DATA_DICTIONARY.md").read_text()
    return re.findall(r"^## outputs/(\S+\.csv)$", dictionary, re.M)


def test_the_data_entry_declares_cff_1_2_0_as_a_dataset():
    text = DATA.read_text()
    assert values("cff-version", text) == ["1.2.0"]
    for key in ("message", "title", "authors"):
        assert re.search(rf"^{key}:", text, re.M), key
    assert re.search(r"^type: dataset$", data_head(), re.M)


def test_the_data_entry_carries_the_same_author_identity():
    assert first_author(DATA) == first_author(CFF)
    affiliations = set(values("affiliation", DATA.read_text()))
    assert affiliations == {"Independent Researcher"}


def test_the_data_entry_is_licensed_as_the_derived_data_are():
    assert values("license", data_head()) == ["CC-BY-4.0"]
    assert "derived data under CC BY 4.0" in (ROOT / "README.md").read_text()
    licenses = " ".join((ROOT / "LICENSES.md").read_text().split())
    assert "Derived data tables released by this project are licensed under" in licenses
    assert "Attribution 4.0 International (CC BY 4.0)" in licenses


def test_the_data_entry_titles_the_drugs_the_genes_and_the_source_release():
    import cohort

    title = values("title", DATA.read_text())[0]
    for token in ("Bedaquiline", "clofazimine", "Mycobacterium tuberculosis",
                  *cohort.BDQ_GENES, "v3.4.0"):
        assert token in title, token


def test_the_data_entry_counts_the_tables_the_dictionary_defines():
    """
    One table is deposited for every column list in the dictionary, so the
    count the abstract opens with is read back from the dictionary, and every
    table it counts is one the README records a writer for.
    """
    tables = defined_tables()
    assert data_abstract().split()[0].lower() == NUMBER_WORDS[len(tables)]
    readme = (ROOT / "README.md").read_text()
    for table in tables:
        assert f"`{table}`" in readme, table


def test_the_data_abstract_says_what_the_code_and_dictionary_say():
    """
    The genes named in the per-variant row are the ones the cohort definition
    reads, the row is the dictionary's, and the prevalence caveat the
    dictionary carries is kept.
    """
    import cohort

    abstract = data_abstract()
    first, second, third = cohort.BDQ_GENES
    phrase = f"row per distinct mutation in {first}, {second} and {third}"
    assert phrase in abstract
    sentence = phrase + abstract.split(phrase)[1].split(". ")[0]
    assert cohort.MODIFIER_GENE.lower() not in sentence.lower()
    dictionary = " ".join(atlas_section(ROOT / "docs" / "DATA_DICTIONARY.md").split())
    assert ("Written by `code/build_atlas.py`. One row per distinct mutation."
            in dictionary)
    caveat = ("Resistance in this collection is concentrated at one site by design, "
              "so no count in this table is a prevalence estimate.")
    assert caveat in dictionary
    assert caveat.replace("this table", "these tables") in abstract
    assert "docs/DATA_DICTIONARY.md defines every column." in abstract


def test_the_data_entry_cites_the_source_release_and_the_software():
    references = DATA.read_text().split("\nreferences:\n")[1]
    software, source = CFF.read_text().split("\nreferences:\n")
    data, cited = references.split("\n  - type: software\n")
    assert values("type", data)[0] == "data"
    assert values("doi", data) == values("doi", source)
    assert values("version", data) == values("version", source)
    assert values("title", cited) == values("title", software)
    assert values("repository-code", cited) == values("repository-code", software)
    assert values("license", cited) == values("license", software)
    assert first_author(DATA) == {key: values(key, cited)[0] for key in (
        "family-names", "given-names", "orcid", "affiliation")}


def test_the_data_entry_cites_the_deposit():
    """
    The entry carries the version DOI of the deposit, the version and release
    date it was published under, and the concept DOI that resolves to the
    latest version. The README carries both DOIs under the same two names.
    """
    head = data_head()
    identifiers = head.split("\nidentifiers:\n")[1].split("\nauthors:\n")[0]
    assert values("doi", head) == [DATA_DOI]
    assert values("version", head) == ["1.0.0"]
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", values("date-released", head)[0])
    assert values("type", identifiers) == ["doi"]
    assert values("value", identifiers) == [DATA_CONCEPT_DOI]
    section = citing_section()
    assert f"version DOI {DATA_DOI}" in section
    assert f"concept DOI {DATA_CONCEPT_DOI}" in section


def test_the_readme_and_the_software_entry_point_at_the_data_entry():
    assert "docs/DATA_CITATION.cff" in values("message")[0]
    section = citing_section()
    assert "`CITATION.cff`" in section
    assert "`docs/DATA_CITATION.cff`" in section
    reference = CFF.read_text().split("\nreferences:\n")[1]
    assert f"v3.4.0, version DOI {values('doi', reference)[0]}" in section
