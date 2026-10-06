"""
Tests that the repository's documentation matches what the code produces.

A column released without an entry in the data dictionary ships undocumented,
and an entry for a column that does not exist misleads. The atlas table is
checked against a run on the synthetic dataset, in tests/test_atlas.py. The
remaining tables are checked here against whatever the modules last wrote, and
a table not yet written is skipped, so the suite collects the same tests and
passes in a fresh clone that has never downloaded the CRyPTIC data.
"""

import csv
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DICTIONARY = ROOT / "docs" / "DATA_DICTIONARY.md"


def documented_columns(table_name):
    text = DICTIONARY.read_text()
    parts = text.split(f"## outputs/{table_name}")
    if len(parts) != 2:
        return None
    body = parts[1].split("\n---")[0]
    return re.findall(r"^\| `([A-Za-z0-9_ ()%]+)` \|", body, re.M)


def written_tables():
    return sorted((ROOT / "outputs").glob("*.csv"))


def listed_tables():
    """The tables the README's outputs section says the modules write."""
    readme = (ROOT / "README.md").read_text()
    section = readme.split("\n## Outputs\n")[1].split("\n## ")[0]
    return sorted(set(re.findall(r"`([\w.]+\.csv)`", section)))


def test_the_dictionary_exists_and_names_its_tables():
    text = DICTIONARY.read_text()
    assert "# Data dictionary" in text
    assert re.findall(r"^## outputs/[\w.]+\.csv$", text, re.M)


@pytest.mark.parametrize("name", listed_tables())
def test_every_listed_table_matches_its_dictionary_entry(name):
    """
    Parametrized over the tables the README lists, so the tests collected do
    not depend on what is on disk. The dictionary entry is required whether or
    not the table has been written; its columns are compared once it has.
    """
    documented = documented_columns(name)
    assert documented is not None, f"{name} has no section in the data dictionary"
    path = ROOT / "outputs" / name
    if not path.exists():
        pytest.skip(f"{name} has not been written; run the module that writes it")
    assert next(csv.reader(path.open())) == documented


def test_every_written_table_is_listed_in_the_readme():
    """A table on disk that the README does not list escapes the check above."""
    assert {path.name for path in written_tables()} <= set(listed_tables())


def runnable_modules():
    """Modules with a main block, which are the ones a reader runs."""
    return sorted(path.name for path in (ROOT / "code").glob("*.py")
                  if "__main__" in path.read_text())


def run_list():
    """The fenced block under the heading that tells a reader what to run."""
    readme = (ROOT / "README.md").read_text()
    section = readme.split("## Running the analysis")[1]
    return section.split("```")[1]


def named_in(section):
    return set(re.findall(r"code/([\w.]+\.py)", section))


def test_the_readme_accounts_for_every_runnable_module():
    assert set(runnable_modules()) <= named_in((ROOT / "README.md").read_text())


def test_the_run_list_holds_every_analysis_module():
    """A module added without a line in that block would be undiscoverable.
    Naming it elsewhere in the README does not put it in a reader's hands.
    The setup modules are exempt: they have their own section."""
    readme = (ROOT / "README.md").read_text()
    setup = readme.split("## Reproducing the environment")[1].split("\n## ")[0]
    expected = set(runnable_modules()) - named_in(setup)
    assert expected
    assert expected <= named_in(run_list())


def test_the_readme_names_no_module_that_is_absent():
    readme = (ROOT / "README.md").read_text()
    named = set(re.findall(r"code/([\w.]+\.py)", readme))
    on_disk = {path.name for path in (ROOT / "code").glob("*.py")}
    assert named <= on_disk


def test_the_readme_points_at_the_data_dictionary():
    readme = (ROOT / "README.md").read_text()
    assert "docs/DATA_DICTIONARY.md" in readme


def documented_values(table_name, column):
    """The values a dictionary entry names in backticks, for a column whose
    entry enumerates them."""
    text = DICTIONARY.read_text()
    body = text.split(f"## outputs/{table_name}")[1].split("\n---")[0]
    row = [line for line in body.splitlines()
           if line.startswith(f"| `{column}` |")]
    return set(re.findall(r"`([^`]+)`", row[0])) - {column} if row else set()


def test_every_kind_of_estimate_is_named_in_the_dictionary():
    """The estimates table stacks several kinds of estimate in one column. A
    kind added without an entry ships unexplained."""
    path = ROOT / "outputs" / "heteroresistance_estimates.csv"
    if not path.exists():
        pytest.skip("the estimates table has not been written")
    rows = list(csv.DictReader(path.open()))
    written = {row["estimate"] for row in rows}
    assert written
    assert written <= documented_values("heteroresistance_estimates.csv",
                                        "estimate")


def collected_test_count():
    """
    The number of tests pytest collects, read from pytest's own collector
    rather than counted by hand. Collection imports the test modules and runs
    no test body, so this does not recurse.
    """
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "--collect-only", "-q",
         "-p", "no:cacheprovider", str(ROOT / "tests")],
        capture_output=True, text=True, cwd=ROOT, check=False)
    assert result.returncode == 0, result.stderr or result.stdout
    per_file = re.findall(r"^tests/\S+\.py: (\d+)$", result.stdout, re.M)
    assert per_file, result.stdout
    return sum(int(count) for count in per_file)


def test_the_readme_states_the_current_test_count():
    """
    A suite size written into prose drifts as tests are added, and a reader has
    no way to tell. The count is compared against the collector.
    """
    readme = (ROOT / "README.md").read_text()
    stated = re.search(r"The suite is ([\d,]+) tests", readme)
    assert stated, "the README does not state a test count"
    assert int(stated.group(1).replace(",", "")) == collected_test_count()


def test_the_readme_names_every_document():
    """A document nobody is pointed at is a document nobody reads."""
    readme = (ROOT / "README.md").read_text()
    missing = [path.name for path in sorted((ROOT / "docs").glob("*.md"))
               if path.name not in readme]
    assert not missing


def test_the_results_document_cites_only_tables_the_pipeline_writes():
    """
    Every cited table is one the README lists as written by a module, which
    holds in a fresh clone. Where tables have been written, every cited one is
    among them.
    """
    results = ROOT / "docs" / "RESULTS.md"
    if not results.exists():
        pytest.skip("no results document")
    cited = set(re.findall(r"outputs/([\w.]+\.csv)", results.read_text()))
    assert cited
    assert cited <= set(listed_tables())
    present = {path.name for path in written_tables()}
    if present:
        assert cited <= present


def test_the_results_document_names_no_quarantined_table():
    """The licence boundary has to hold in the released text, not only in the
    code that produced it."""
    results = ROOT / "docs" / "RESULTS.md"
    if not results.exists():
        pytest.skip("no results document")
    paragraphs = results.read_text().split("\n\n")
    for table in ("EFFECTS", "PREDICTIONS", "ANTIBIOGRAM"):
        for paragraph in paragraphs:
            if table in paragraph:
                assert "excluded from" in paragraph, paragraph
