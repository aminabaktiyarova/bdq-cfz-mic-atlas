"""
Tests that micecoff/README.md describes what the package does.

The README is the package's documentation and its data dictionary: it defines
every column the command writes, every option and its default, the reasons an
estimate is withheld, the names the intervals are seeded by, and the constants
the method uses. Each of those is read back from the README here and compared
with the code, so a change to either side that the other does not follow fails
the suite.
"""

import argparse
import csv
import io
import math
import os
import pathlib
import re
import shlex
import subprocess
import sys

import numpy as np
import pytest

from micecoff import cli, core

ROOT = pathlib.Path(__file__).resolve().parents[1]
README = ROOT / "micecoff" / "README.md"
LADDER = [0.008, 0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1.0]
COLUMNS = {"fit": cli.FIT_COLUMNS, "shift": cli.SHIFT_COLUMNS,
           "ecoff": cli.ECOFF_COLUMNS}


def text():
    return README.read_text()


def flat(value):
    """Prose with its line breaks and runs of spaces collapsed."""
    return " ".join(value.split())


def section(heading):
    """The body under a level-two heading, up to the next one."""
    parts = text().split(f"\n## {heading}\n")
    assert len(parts) == 2, f"no single section headed {heading!r}"
    return parts[1].split("\n## ")[0]


def subsection(heading):
    parts = section("Output").split(f"\n### {heading}\n")
    assert len(parts) == 2, f"no single subsection headed {heading!r}"
    return parts[1].split("\n### ")[0]


def table(body):
    """The rows of a Markdown table whose first cell is a backticked name."""
    rows = []
    for line in body.splitlines():
        match = re.match(r"^\| `([^`]+)` \|(.*)\|$", line)
        if match:
            rows.append([match.group(1)] + [cell.strip()
                                             for cell in match.group(2).split("|")])
    return rows


def fenced(body):
    return re.findall(r"```\n(.*?)```", body, re.S)


def planted_input(path):
    """A file with the README's columns and labels, large enough to estimate."""
    rng = np.random.default_rng(20260101)
    lines = ["MIC,group,cluster"]
    for label, mu, n in (("wild type", -5.0, 200), ("carrier", -3.5, 120)):
        for index, reading in enumerate(core.simulate_reports(mu, 1.0, LADDER, rng,
                                                              draws=n)):
            lines.append(f"{reading},{label},{label[:4]}{index % 20}")
    path.write_text("\n".join(lines) + "\n")
    return path


def run(capsys, argv):
    assert cli.main(argv) == 0
    return list(csv.DictReader(io.StringIO(capsys.readouterr().out)))


# Columns


@pytest.mark.parametrize("command", sorted(COLUMNS))
def test_each_subcommand_documents_its_columns_in_order(command):
    documented = [row[0] for row in table(subsection(command))]
    assert documented == COLUMNS[command]


@pytest.mark.parametrize("command", sorted(COLUMNS))
def test_every_column_entry_states_type_definition_and_when_empty(command):
    for row in table(subsection(command)):
        assert len(row) == 4, row
        assert all(row), row


SERIES = "0.008,0.015,0.03,0.06,0.12,0.25,0.5,1"


@pytest.mark.parametrize("command,extra", [
    ("fit", []), ("shift", ["--reference", "wild type"]),
    ("ecoff", ["--coverage", "0.95"])])
def test_the_empty_cells_are_the_ones_the_readme_says(tmp_path, capsys, command,
                                                       extra):
    """
    One group reported and one withheld. A column documented as never empty
    is filled in both rows, a column documented as empty when withheld is
    empty only in the withheld row, and withheld itself is empty only where
    the estimate is reported.
    """
    path = planted_input(tmp_path / "mics.csv")
    with open(path, "a") as handle:
        handle.writelines("<=0.008,floor,f\n" for _ in range(30))
    rows = {row["group"]: row
            for row in run(capsys, [command, str(path), "--series", SERIES,
                                    "--group-column", "group", "--draws", "100",
                                    *extra])}
    reported, withheld = rows["carrier"], rows["floor"]
    assert reported["withheld"] == "" and withheld["withheld"] != ""
    for name, _, _, empty_when in table(subsection(command)):
        if empty_when == "never":
            assert reported[name] != "" and withheld[name] != "", name
        elif empty_when == "withheld":
            assert reported[name] != "" and withheld[name] == "", name
        elif empty_when == "the estimate is reported":
            assert reported[name] == "" and withheld[name] != "", name
        else:
            pytest.fail(f"{name}: no rule for {empty_when!r}")


# Options


def subparsers():
    parser = cli.build_parser()
    choices = next(action for action in parser._actions
                   if isinstance(action, argparse._SubParsersAction)).choices
    return parser, choices


def option_name(action):
    if not action.option_strings:
        return action.dest
    return max(action.option_strings, key=len)


def default_cell(action):
    if action.required:
        return "required"
    if isinstance(action, argparse._StoreTrueAction):
        return "off"
    if action.default in (None, argparse.SUPPRESS):
        return "none"
    return f"`{action.default}`"


def code_options():
    """Each option with the subcommands it belongs to and its expected cells."""
    parser, choices = subparsers()
    found = {}
    for command, sub in choices.items():
        for action in sub._actions:
            if isinstance(action, argparse._HelpAction):
                continue
            found.setdefault(option_name(action), []).append((command, action))
    expected = {}
    for name, entries in found.items():
        commands = {command for command, _ in entries}
        required = sorted(command for command, action in entries if action.required)
        optional = {default_cell(action) for _, action in entries
                    if not action.required}
        if commands == set(choices):
            where = "all" + (", required by " + ", ".join(
                f"`{command}`" for command in required)
                if required and optional else "")
        else:
            where = ", ".join(f"`{command}`" for command in sorted(commands))
        if not optional:
            optional = {"required"}
        assert len(optional) == 1, (name, optional)
        expected[name] = (where, optional.pop())
    for action in parser._actions:
        if action.option_strings and not isinstance(action, argparse._HelpAction):
            expected[option_name(action)] = ("none", "none")
    return expected


def test_every_option_is_documented_with_its_subcommands_and_default():
    documented = {row[0]: (row[1], row[2]) for row in table(section("Commands"))}
    assert documented == code_options()


def test_the_bounds_the_options_table_states_are_enforced():
    """
    Each "at least n" in the table is the smallest value the parser takes, the
    seed is refused below zero, and the coverage is refused at 0 and 1.
    """
    documented = {row[0]: row[3] for row in table(section("Commands"))}
    common = ["fit", "x.csv", "--series", SERIES]

    def accepted(option, value):
        try:
            cli.build_parser().parse_args(common + [option, value])
        except SystemExit:
            return False
        return True

    checked = 0
    for option, meaning in documented.items():
        bound = re.search(r"at least (\d+)", meaning)
        if bound:
            least = int(bound.group(1))
            assert accepted(option, str(least)) and not accepted(option, str(least - 1))
            checked += 1
    assert checked == 2
    assert "non-negative" in documented["--seed"]
    assert accepted("--seed", "0") and not accepted("--seed", "-1")
    assert "strictly between 0 and 1" in documented["--coverage"]
    coverage = ["ecoff", "x.csv", "--series", SERIES, "--coverage"]
    for value, ok in (("0", False), ("1", False), ("0.5", True)):
        try:
            cli.build_parser().parse_args(coverage + [value])
            assert ok
        except SystemExit:
            assert not ok


# Examples


def test_the_readme_commands_run_and_write_the_documented_columns(tmp_path, capsys,
                                                                  monkeypatch):
    blocks = fenced(section("Commands"))
    commands = [line for block in blocks for line in block.splitlines()
                if line.startswith("micecoff ")]
    assert sorted(shlex.split(line)[1] for line in commands) == sorted(COLUMNS)
    planted_input(tmp_path / "mics.csv")
    monkeypatch.chdir(tmp_path)
    for line in commands:
        argv = shlex.split(line)[1:]
        rows = run(capsys, argv)
        assert list(rows[0]) == COLUMNS[argv[0]]
        assert all(row["withheld"] == "" for row in rows), line


def test_the_input_excerpt_places_on_the_example_series(tmp_path, capsys):
    excerpt = fenced(section("Input"))[0]
    path = tmp_path / "mics.csv"
    path.write_text(excerpt)
    series = shlex.split(fenced(section("Commands"))[0].splitlines()[0])
    series = series[series.index("--series") + 1]
    rows = run(capsys, ["fit", str(path), "--series", series, "--group-column",
                        "group", "--draws", "100"])
    assert sum(int(row["placed"]) for row in rows) == len(excerpt.splitlines()) - 1
    assert all(row["off_series"] == "0" for row in rows)


# Withholding


def documented_reasons():
    body = section("Withheld estimates")
    reasons = re.findall(r"^- `([^`]+)`", body, re.M)
    assert reasons
    return reasons


def pattern(reason):
    return re.compile(re.escape(reason).replace(re.escape("<n>"), r"\d+"))


def produced_reasons():
    """Every reason the code writes, produced by the case that writes it."""

    def bounds(readings):
        lower, upper, _, _ = core.place_mics(readings, LADDER)
        return lower, upper

    three = bounds(["0.06", "0.12", "0.25"])
    out = [
        core.withhold_reason(None, *bounds([None]), LADDER),
        core.withhold_reason(None, *bounds(["0.06", "0.12"]), LADDER),
        core.withhold_reason(None, *three, LADDER),
        core.withhold_reason({"mu": core.MU_BOUNDS[0], "sigma": 1.0}, *three, LADDER),
        core.withhold_reason({"mu": -9.0, "sigma": 1.0}, *three, LADDER),
        cli.DRAWS_FAILED,
    ]
    lower, upper = bounds(["0.06"] * 10 + ["0.12"] * 10 + ["0.25"] * 10)
    data = {"lower": lower, "upper": upper, "series": LADDER,
            "clusters": np.array(["a", "b", "c"] * 10)}
    _, units = cli.assess(data, np.arange(30), 5)
    out.append(units)
    out.append("reference: " + units)
    return out


def test_every_reason_the_code_writes_is_documented_and_no_other():
    documented = documented_reasons()
    prefix = [reason for reason in documented if reason.endswith(": ")]
    assert prefix == ["reference: "]
    patterns = [pattern(reason) for reason in documented if reason not in prefix]
    matched = set()
    for reason in produced_reasons():
        assert reason is not None
        if reason.startswith(prefix[0]):
            reason = reason[len(prefix[0]):]
            matched.add(prefix[0])
        hits = [p.pattern for p in patterns if p.fullmatch(reason)]
        assert len(hits) == 1, reason
        matched.add(hits[0])
    assert matched == {p.pattern for p in patterns} | set(prefix)


def test_the_reasons_are_checked_in_the_order_the_readme_lists_them():
    """Where two reasons apply, the code reports the one listed first."""
    order = [pattern(reason) for reason in documented_reasons()]

    def rank(reason):
        return next(i for i, p in enumerate(order) if p.fullmatch(reason))

    def bounds(readings):
        lower, upper, _, _ = core.place_mics(readings, LADDER)
        return lower, upper

    two = bounds(["0.06", "0.12"])
    three = bounds(["0.06", "0.12", "0.25"])
    assert rank(core.withhold_reason(None, *bounds([None]), LADDER)) == 0
    assert rank(core.withhold_reason(None, *two, LADDER)) == 1
    on_bound_and_outside = {"mu": core.MU_BOUNDS[0], "sigma": 1.0}
    assert rank(core.withhold_reason(on_bound_and_outside, *three, LADDER)) == 3
    lower, upper = bounds(["<=0.008"] * 28 + ["0.015", "0.03"])
    data = {"lower": lower, "upper": upper, "series": LADDER,
            "clusters": np.array(["a", "b", "c"] * 10)}
    fit, reason = cli.assess(data, np.arange(30), 5)
    assert fit["mu"] < math.log2(LADDER[0]) - core.EXTRAPOLATION_MARGIN
    assert rank(reason) == 4


# Seeding


def test_the_documented_stream_names_reproduce_the_intervals(tmp_path, capsys):
    names = re.findall(r"`((?:fit|shift|ecoff)\|[^`]+)`", flat(section("Intervals")))
    assert sorted(name.split("|")[0] for name in names) == ["ecoff", "fit", "shift"]
    stream = {name.split("|")[0]: name for name in names}
    path = planted_input(tmp_path / "mics.csv")
    rows = [line.split(",") for line in path.read_text().splitlines()[1:]]
    bounds = {}
    for label in ("wild type", "carrier"):
        readings = [r[0] for r in rows if r[1] == label]
        lower, upper, _, _ = core.place_mics(readings, LADDER)
        bounds[label] = (lower, upper, np.arange(len(readings)))
    seed = 20260101
    common = [str(path), "--series", SERIES, "--group-column", "group",
              "--draws", "100", "--seed", str(seed)]

    def named(kind, group, reference=None):
        name = stream[kind].replace("<group>", group)
        return core.named_generator(seed, name.replace("<reference>", reference or ""))

    fit = {row["group"]: row for row in run(capsys, ["fit", *common])}
    mean = core.resample_clusters(*bounds["carrier"], named("fit", "carrier"),
                                  draws=100)
    assert fit["carrier"]["mean_low"] == f"{mean['low']:.4f}"

    shift = run(capsys, ["shift", *common, "--reference", "wild type"])
    interval = core.resample_shift(*bounds["carrier"], *bounds["wild type"],
                                   named("shift", "carrier", "wild type"), draws=100)
    assert shift[0]["shift_low"] == f"{interval['low']:.4f}"

    ecoff = {row["group"]: row for row in run(capsys, ["ecoff", *common,
                                                       "--coverage", "0.95"])}
    quantile = core.resample_clusters(*bounds["carrier"], named("ecoff", "carrier"),
                                      draws=100, statistic=cli.quantile_of(0.95))
    assert ecoff["carrier"]["quantile_low"] == f"{quantile['low']:.4f}"


# Constants and metadata


def test_the_constants_the_readme_states_are_the_codes():
    prose = flat(text())
    assert f"within {core.MATCH_TOLERANCE:.0%}" in prose
    low, high = core.MU_BOUNDS
    sd_low, sd_high = (math.exp(value) for value in core.LOG_SIGMA_BOUNDS)
    assert (f"means from {low:g} to {high:g} and standard deviations from "
            f"{sd_low:.2g} to {sd_high:.2g}") in prose
    assert (f"from the {cli.INTERVAL[0]:g}th to the {cli.INTERVAL[1]:g}th "
            "percentile") in prose
    assert f"{cli.INTERVAL[1] - cli.INTERVAL[0]:g}% interval" in prose
    assert "four decimal places" in prose and cli.log2_text(1.234567) == "1.2346"
    assert ("four significant figures" in prose
            and cli.concentration_text(0.0123456) == "0.01235")


def project():
    tomllib = pytest.importorskip("tomllib")
    return tomllib.loads((ROOT / "pyproject.toml").read_text())


def test_the_package_metadata_carries_this_readme():
    assert project()["project"]["readme"] == "micecoff/README.md"
    assert README.is_file()


def test_the_readme_states_the_declared_requirements():
    configuration = project()["project"]
    python = re.fullmatch(r">=(\d+\.\d+)", configuration["requires-python"]).group(1)
    floors = dict(re.fullmatch(r"([a-z]+)>=([\d.]+)", item).groups()
                  for item in configuration["dependencies"])
    assert sorted(floors) == ["numpy", "scipy"]
    assert (f"requires Python {python} or later, numpy {floors['numpy']} or later "
            f"and scipy {floors['scipy']} or later") in flat(text())


def run_module(*argv, cwd):
    environment = dict(os.environ, PYTHONPATH=str(ROOT), PYTHONDONTWRITEBYTECODE="1")
    return subprocess.run([sys.executable, "-m", "micecoff", *argv], cwd=cwd,
                          env=environment, capture_output=True, text=True,
                          timeout=120)


def test_the_exit_status_is_the_one_the_readme_states(tmp_path):
    prose = flat(section("Output"))
    assert ("The exit status is 0 when the run completes, including a run in which "
            "estimates are withheld, 1 when the input stops the run, and 2 when "
            "the arguments are invalid") in prose
    path = planted_input(tmp_path / "mics.csv")
    with open(path, "a") as handle:
        handle.writelines("<=0.008,floor,f\n" for _ in range(30))
    completed = run_module("fit", str(path), "--series", SERIES, "--group-column",
                           "group", "--draws", "100", cwd=tmp_path)
    assert completed.returncode == 0 and "the readings occupy" in completed.stdout
    stopped = run_module("fit", str(path), "--series", "0.5,1", cwd=tmp_path)
    assert stopped.returncode == 1
    invalid = run_module("fit", str(path), "--series", SERIES, "--draws", "1",
                         cwd=tmp_path)
    assert invalid.returncode == 2


def test_every_link_resolves_outside_the_repository():
    """
    The README is the package's description wherever the package is installed,
    where a path relative to this repository leads nowhere.
    """
    targets = re.findall(r"\]\(([^)]+)\)", text())
    urls = re.findall(r"\bhttps?://\S+", text())
    assert urls
    assert all(target.startswith("https://") for target in targets)
    assert all(url.startswith("https://") for url in urls)


def test_the_repository_readme_points_at_the_package():
    readme = (ROOT / "README.md").read_text()
    assert "micecoff/README.md" in readme
    assert "pip install ." in readme
    assert "micecoff/README.md" in (ROOT / "docs" / "DATA_DICTIONARY.md").read_text()
