"""
Tests for the micecoff command line interface and its packaging.

Every input is written by the test with its parameters planted, so a test
checks that the interface recovers something known. The interface is called
through main(argv), except where the test concerns how the package runs or
installs.
"""

import ast
import csv
import io
import os
import pathlib
import re
import subprocess
import sys

import numpy as np
import pytest
from scipy.stats import norm

from micecoff import __version__, cli, core

ROOT = pathlib.Path(__file__).resolve().parents[1]
LADDER = [0.008, 0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1.0]
SERIES = "0.008,0.015,0.03,0.06,0.12,0.25,0.5,1"


def reports(mu, sd, n, seed):
    return core.simulate_reports(mu, sd, LADDER, np.random.default_rng(seed),
                                 draws=n)


def write_csv(path, rows, header=("MIC", "group", "cluster"), encoding="utf-8"):
    with open(path, "w", newline="", encoding=encoding) as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return str(path)


def planted(path, groups, seed=20260101, clones=None):
    """
    One row per simulated report. groups maps a label to (mu, sd, n); clones
    maps a label to the number of cluster labels its rows cycle through.
    """
    rows = []
    for offset, (name, (mu, sd, n)) in enumerate(groups.items()):
        size = (clones or {}).get(name, n)
        for index, text in enumerate(reports(mu, sd, n, seed + offset)):
            rows.append((text, name, f"{name}-{index % size}"))
    return write_csv(path, rows)


def run(capsys, *argv):
    assert cli.main([str(a) for a in argv]) == 0
    captured = capsys.readouterr()
    return list(csv.DictReader(io.StringIO(captured.out))), captured.err


def stopped(capsys, *argv):
    with pytest.raises(SystemExit) as raised:
        cli.main([str(a) for a in argv])
    capsys.readouterr()
    return str(raised.value.code)


def refused(capsys, *argv):
    """argparse rejects the arguments with its usage exit code."""
    with pytest.raises(SystemExit) as raised:
        cli.main([str(a) for a in argv])
    capsys.readouterr()
    return raised.value.code


def by_group(rows):
    return {row["group"]: row for row in rows}


def number(text):
    return float(text)


def brackets(row, estimate, low, high, expected_width):
    """
    The interval holds its own estimate and is within a factor of two of the
    width a normal approximation predicts. Whether an interval covers the
    truth at its nominal rate is a property of many datasets and was measured
    by simulation; one fixed dataset can only check that the interval belongs
    to its estimate and is of the right size.
    """
    value, lower, upper = number(row[estimate]), number(row[low]), number(row[high])
    assert lower < value < upper
    assert 0.5 * expected_width < upper - lower < 2.0 * expected_width


# The interface: column order is held here, so a rename is a visible change.


def test_each_subcommand_writes_its_documented_columns(tmp_path, capsys):
    path = planted(tmp_path / "in.csv", {"a": (-5.0, 1.0, 200), "b": (-4.0, 1.0, 200)})
    expected = {
        "fit": ["group", "rows", "absent", "off_series", "placed",
                "left_censored", "right_censored", "resampling_unit",
                "clusters", "mean_log2", "mean_low", "mean_high", "sd_log2",
                "sd_low", "sd_high", "geometric_mean", "geometric_mean_low",
                "geometric_mean_high", "draws", "withheld"],
        "shift": ["group", "reference", "rows", "absent", "off_series",
                  "placed", "reference_placed", "resampling_unit", "clusters",
                  "reference_clusters", "mean_log2", "reference_mean_log2",
                  "shift", "shift_low", "shift_high", "draws", "withheld"],
        "ecoff": ["group", "rows", "absent", "off_series", "placed",
                  "resampling_unit", "clusters", "mean_log2", "sd_log2",
                  "coverage", "quantile_log2", "quantile_low",
                  "quantile_high", "ecoff", "ecoff_low", "ecoff_high",
                  "draws", "withheld"],
    }
    extra = {"fit": [], "shift": ["--reference", "a"],
             "ecoff": ["--coverage", "0.99"]}
    for command, columns in expected.items():
        cli.main([command, path, "--series", SERIES, "--group-column", "group",
                  "--draws", "100"] + extra[command])
        header = capsys.readouterr().out.splitlines()[0]
        assert header.split(",") == columns


# fit


def test_fit_recovers_planted_parameters_for_each_group(tmp_path, capsys):
    """Planted in reverse order of name, so the output order is the sort's."""
    truth = {"wild": (-5.3, 1.0, 800), "carrier": (-3.5, 1.2, 500)}
    path = planted(tmp_path / "in.csv", truth)
    rows, _ = run(capsys, "fit", path, "--series", SERIES, "--group-column",
                  "group", "--draws", "200")
    assert [row["group"] for row in rows] == ["carrier", "wild"]
    for name, (mu, sd, n) in truth.items():
        row = by_group(rows)[name]
        assert row["withheld"] == ""
        assert int(row["placed"]) == n
        assert number(row["mean_log2"]) == pytest.approx(mu, abs=0.12)
        brackets(row, "mean_log2", "mean_low", "mean_high", 3.92 * sd / n ** 0.5)
        assert number(row["sd_log2"]) == pytest.approx(sd, abs=0.12)
        brackets(row, "sd_log2", "sd_low", "sd_high", 3.92 * sd / (2 * n) ** 0.5)
        assert number(row["geometric_mean"]) == pytest.approx(
            2 ** number(row["mean_log2"]), rel=1e-3)
        assert number(row["geometric_mean_low"]) == pytest.approx(
            2 ** number(row["mean_low"]), rel=1e-3)
        assert int(row["draws"]) == 200


def test_fit_matches_the_estimator_it_wraps(tmp_path, capsys):
    path = planted(tmp_path / "in.csv", {"a": (-5.0, 1.0, 300)})
    rows, _ = run(capsys, "fit", path, "--series", SERIES, "--draws", "100")
    lower, upper, _, _ = core.place_mics(reports(-5.0, 1.0, 300, 20260101), LADDER)
    fit = core.fit_censored_normal(lower, upper)
    assert rows[0]["mean_log2"] == f"{fit['mu']:.4f}"
    assert rows[0]["sd_log2"] == f"{fit['sigma']:.4f}"


def test_fit_without_a_group_column_reports_one_group(tmp_path, capsys):
    path = planted(tmp_path / "in.csv", {"a": (-5.0, 1.0, 100), "b": (-4.0, 1.0, 100)})
    rows, _ = run(capsys, "fit", path, "--series", SERIES, "--draws", "100")
    assert len(rows) == 1
    assert rows[0]["group"] == "all" and rows[0]["placed"] == "200"


def test_censored_counts_follow_the_readings(tmp_path, capsys):
    texts = reports(-3.5, 2.0, 400, 3)
    path = write_csv(tmp_path / "in.csv", [(t, "a", "c") for t in texts])
    rows, _ = run(capsys, "fit", path, "--series", SERIES, "--draws", "100")
    assert int(rows[0]["left_censored"]) == sum(t.startswith("<=") for t in texts)
    assert int(rows[0]["right_censored"]) == sum(t.startswith(">") for t in texts)
    assert int(rows[0]["left_censored"]) > 0 and int(rows[0]["right_censored"]) > 0


def test_absent_readings_are_counted_and_leave_the_estimate_unchanged(tmp_path,
                                                                       capsys):
    texts = reports(-5.0, 1.0, 200, 4)
    clean = write_csv(tmp_path / "clean.csv", [(t, "a", "c") for t in texts])
    padded_rows = []
    for index, text in enumerate(texts):
        padded_rows.append((text, "a", "c"))
        if index % 4 == 0:
            padded_rows.append(("", "a", "c"))
    padded = write_csv(tmp_path / "padded.csv", padded_rows)
    plain, _ = run(capsys, "fit", clean, "--series", SERIES, "--draws", "100")
    gapped, err = run(capsys, "fit", padded, "--series", SERIES, "--draws", "100")
    assert gapped[0]["absent"] == "50" and gapped[0]["rows"] == "250"
    assert "50 rows carry no reading" in err
    for key in plain[0]:
        if key not in ("rows", "absent"):
            assert gapped[0][key] == plain[0][key], key


def test_every_row_is_absent_off_series_or_placed(tmp_path, capsys):
    rows = [("0.06", "a", "c")] * 30 + [("0.12", "a", "c")] * 30
    rows += [("0.25", "a", "c")] * 30 + [("", "a", "c")] * 4 + [("0.2", "a", "c")] * 3
    rows += [("<=0.015", "a", "c")] * 2
    path = write_csv(tmp_path / "in.csv", rows)
    out, err = run(capsys, "fit", path, "--series", SERIES, "--draws", "100",
                   "--allow-off-series")
    row = out[0]
    assert (row["rows"], row["absent"], row["off_series"], row["placed"]) == (
        "99", "4", "5", "90")
    assert row["left_censored"] == "0"
    assert "5 readings off the series were dropped" in err


# Withholding


def test_a_withheld_group_keeps_its_counts_and_states_its_reason(tmp_path, capsys):
    rows = [("<=0.008", "floor", "c")] * 40
    rows += [(t, "wild", "c") for t in reports(-5.0, 1.0, 200, 5)]
    path = write_csv(tmp_path / "in.csv", rows)
    out, _ = run(capsys, "fit", path, "--series", SERIES, "--group-column", "group",
                 "--draws", "100")
    floor = by_group(out)["floor"]
    assert floor["placed"] == "40" and floor["left_censored"] == "40"
    assert "too few intervals" in floor["withheld"]
    assert floor["mean_log2"] == floor["mean_low"] == floor["draws"] == ""
    assert by_group(out)["wild"]["withheld"] == ""


def test_a_group_with_fewer_clusters_than_the_minimum_is_withheld(tmp_path, capsys):
    path = planted(tmp_path / "in.csv", {"a": (-5.0, 1.0, 120)}, clones={"a": 4})
    out, _ = run(capsys, "fit", path, "--series", SERIES, "--cluster-column",
                 "cluster", "--draws", "100")
    assert out[0]["clusters"] == "4"
    assert out[0]["withheld"] == (
        "too few resampling units (4 of the 5 an interval needs)")
    out, _ = run(capsys, "fit", path, "--series", SERIES, "--cluster-column",
                 "cluster", "--draws", "100", "--min-clusters", "4")
    assert out[0]["withheld"] == "" and out[0]["mean_log2"] != ""


def test_a_mean_extrapolated_below_the_plate_is_withheld(tmp_path, capsys):
    texts = reports(-9.5, 1.5, 3000, 6)
    lower, upper, _, _ = core.place_mics(texts, LADDER)
    occupied = len(set(zip(lower.tolist(), upper.tolist())))
    assert occupied >= core.MINIMUM_INTERVALS, "the planted data must pass the "\
        "interval rule so that the extrapolation rule is the one tested"
    path = write_csv(tmp_path / "in.csv", [(t, "a", "c") for t in texts])
    out, _ = run(capsys, "fit", path, "--series", SERIES, "--draws", "100")
    assert "outside the tested series" in out[0]["withheld"]


# shift


def test_shift_recovers_a_planted_shift(tmp_path, capsys):
    path = planted(tmp_path / "in.csv", {"carrier": (-3.2, 1.0, 300),
                                         "wild": (-5.0, 1.0, 900)})
    out, _ = run(capsys, "shift", path, "--series", SERIES, "--group-column",
                 "group", "--reference", "wild", "--draws", "200")
    assert [row["group"] for row in out] == ["carrier"]
    row = out[0]
    assert row["reference"] == "wild" and row["reference_placed"] == "900"
    assert number(row["shift"]) == pytest.approx(1.8, abs=0.15)
    brackets(row, "shift", "shift_low", "shift_high", 3.92 * (1 / 300 + 1 / 900) ** 0.5)
    assert number(row["shift"]) == pytest.approx(
        number(row["mean_log2"]) - number(row["reference_mean_log2"]), abs=2e-4)


def test_the_shift_interval_carries_the_reference_groups_error(tmp_path, capsys):
    """
    Equal groups of sixty. The shift's interval must be wider than the
    group's own interval on its mean, near 1.4 times, because the reference
    mean is estimated from as few observations as the group's.
    """
    path = planted(tmp_path / "in.csv", {"carrier": (-4.0, 1.0, 60),
                                         "wild": (-5.0, 1.0, 60)})
    shift, _ = run(capsys, "shift", path, "--series", SERIES, "--group-column",
                   "group", "--reference", "wild", "--draws", "300")
    fit, _ = run(capsys, "fit", path, "--series", SERIES, "--group-column",
                 "group", "--draws", "300")
    shift_width = number(shift[0]["shift_high"]) - number(shift[0]["shift_low"])
    carrier = by_group(fit)["carrier"]
    mean_width = number(carrier["mean_high"]) - number(carrier["mean_low"])
    assert shift_width > 1.2 * mean_width


def test_a_cluster_column_widens_the_shift_interval(tmp_path, capsys):
    rng = np.random.default_rng(7)
    steps = np.log2(np.asarray(LADDER))
    rows = []
    for clone in range(8):
        center = rng.normal(-3.5, 1.0)
        for member in range(6):
            value = center + rng.normal(0, 0.05)
            index = int(np.searchsorted(steps, value))
            text = (f"<={LADDER[0]}" if index == 0 else f">{LADDER[-1]}"
                    if index >= len(steps) else f"{LADDER[index]}")
            rows.append((text, "carrier", f"clone{clone}"))
    rows += [(t, "wild", f"w{i}") for i, t in enumerate(reports(-5.0, 1.0, 600, 8))]
    path = write_csv(tmp_path / "in.csv", rows)
    common = ["shift", path, "--series", SERIES, "--group-column", "group",
              "--reference", "wild", "--draws", "300"]
    by_row, _ = run(capsys, *common)
    by_cluster, _ = run(capsys, *common, "--cluster-column", "cluster")
    assert by_row[0]["resampling_unit"] == "row" and by_row[0]["clusters"] == "48"
    assert by_cluster[0]["resampling_unit"] == "cluster"
    assert by_cluster[0]["clusters"] == "8"
    assert by_cluster[0]["shift"] == by_row[0]["shift"]

    def width(row):
        return number(row["shift_high"]) - number(row["shift_low"])

    assert width(by_cluster[0]) > 1.5 * width(by_row[0])


def test_a_withheld_reference_withholds_every_shift(tmp_path, capsys):
    rows = [("<=0.008", "wild", "c")] * 50
    rows += [(t, "carrier", "c") for t in reports(-4.0, 1.0, 100, 9)]
    path = write_csv(tmp_path / "in.csv", rows)
    out, _ = run(capsys, "shift", path, "--series", SERIES, "--group-column",
                 "group", "--reference", "wild", "--draws", "100")
    assert out[0]["withheld"].startswith("reference: ")
    assert out[0]["shift"] == ""


def test_shift_refuses_a_reference_the_file_does_not_hold(tmp_path, capsys):
    path = planted(tmp_path / "in.csv", {"a": (-5.0, 1.0, 50), "b": (-4.0, 1.0, 50)})
    message = stopped(capsys, "shift", path, "--series", SERIES, "--group-column",
                      "group", "--reference", "wild")
    assert "no group 'wild'" in message


def test_shift_refuses_a_file_holding_only_the_reference(tmp_path, capsys):
    path = planted(tmp_path / "in.csv", {"wild": (-5.0, 1.0, 50)})
    message = stopped(capsys, "shift", path, "--series", SERIES, "--group-column",
                      "group", "--reference", "wild")
    assert "holds only the reference group" in message


def test_shift_requires_a_group_column(tmp_path, capsys):
    path = planted(tmp_path / "in.csv", {"wild": (-5.0, 1.0, 50)})
    assert refused(capsys, "shift", path, "--series", SERIES,
                   "--reference", "wild") == 2


# ecoff


def test_ecoff_recovers_a_planted_wild_type_boundary(tmp_path, capsys):
    """
    A wild type at 0.015 with a spread of half a doubling has its 99th
    percentile between the 0.03 and 0.06 wells, so the cut-off is 0.06.
    """
    mu, sd = np.log2(0.015), 0.5
    planted_quantile = norm.ppf(0.99, loc=mu, scale=sd)
    assert np.log2(0.03) < planted_quantile <= np.log2(0.06)
    path = planted(tmp_path / "in.csv", {"wild": (mu, sd, 3000)})
    out, _ = run(capsys, "ecoff", path, "--series", SERIES, "--coverage", "0.99",
                 "--draws", "200")
    row = out[0]
    assert row["coverage"] == "0.99"
    assert row["ecoff"] == "0.06"
    assert number(row["quantile_log2"]) == pytest.approx(planted_quantile, abs=0.15)
    # The standard error of a normal quantile at z is sd * sqrt((1 + z^2/2) / n).
    z = norm.ppf(0.99)
    brackets(row, "quantile_log2", "quantile_low", "quantile_high",
             3.92 * sd * ((1 + z * z / 2) / 3000) ** 0.5)
    assert row["ecoff_low"] in ("0.03", "0.06")
    assert row["ecoff_high"] in ("0.06", "0.12")


def test_ecoff_bounds_round_up_from_the_quantile_bounds(tmp_path, capsys):
    """
    The planted 95th percentile sits on the 0.12 well, so the interval on the
    quantile straddles it and its two ends round to different wells, which is
    the case where a bound read from the wrong end is visible.
    """
    mu = np.log2(0.12) - norm.ppf(0.95)
    path = planted(tmp_path / "in.csv", {"wild": (mu, 1.0, 200)})
    out, _ = run(capsys, "ecoff", path, "--series", SERIES, "--coverage", "0.95",
                 "--draws", "200")
    row = out[0]
    assert row["ecoff_low"] != row["ecoff_high"]
    for quantile, well in (("quantile_log2", "ecoff"), ("quantile_low", "ecoff_low"),
                           ("quantile_high", "ecoff_high")):
        expected = core.round_up_to_series(number(row[quantile]), LADDER)
        assert number(row[well]) == expected


def test_a_cut_off_above_the_plate_is_written_as_above_the_highest_well(tmp_path,
                                                                        capsys):
    path = planted(tmp_path / "in.csv", {"wild": (np.log2(0.5), 1.0, 1000)})
    out, _ = run(capsys, "ecoff", path, "--series", SERIES, "--coverage", "0.99",
                 "--draws", "100")
    assert out[0]["ecoff"] == ">1"
    assert number(out[0]["quantile_log2"]) > 0


@pytest.mark.parametrize("coverage", ["0", "1", "1.5", "-0.1", "high"])
def test_ecoff_refuses_a_coverage_outside_the_unit_interval(tmp_path, capsys,
                                                            coverage):
    path = planted(tmp_path / "in.csv", {"wild": (-5.0, 1.0, 50)})
    assert refused(capsys, "ecoff", path, "--series", SERIES,
                   "--coverage", coverage) == 2


def test_ecoff_has_no_default_coverage(tmp_path, capsys):
    path = planted(tmp_path / "in.csv", {"wild": (-5.0, 1.0, 50)})
    assert refused(capsys, "ecoff", path, "--series", SERIES) == 2


# Seeding and determinism


def test_a_groups_interval_does_not_depend_on_the_other_groups(tmp_path, capsys):
    groups = {"a": (-5.0, 1.0, 150), "b": (-4.5, 1.0, 150), "c": (-4.0, 1.0, 150)}
    full = planted(tmp_path / "full.csv", groups)
    rows = list(csv.reader(open(full)))
    write_csv(tmp_path / "without_b.csv", [r for r in rows[1:] if r[1] != "b"])
    for command, extra in (("fit", []), ("ecoff", ["--coverage", "0.99"])):
        everything, _ = run(capsys, command, full, "--series", SERIES,
                            "--group-column", "group", "--draws", "100", *extra)
        fewer, _ = run(capsys, command, str(tmp_path / "without_b.csv"),
                       "--series", SERIES, "--group-column", "group",
                       "--draws", "100", *extra)
        assert by_group(fewer)["c"] == by_group(everything)["c"]


def test_each_interval_is_reproducible_from_its_stream_name(tmp_path, capsys):
    """
    A group's interval draws from the generator named by the subcommand and
    the group, so it can be recomputed alone from the seed and that name.
    """
    texts = {"a": reports(-5.0, 1.0, 150, 21), "b": reports(-4.0, 1.0, 150, 22)}
    rows = [(t, name, "c") for name, ts in texts.items() for t in ts]
    path = write_csv(tmp_path / "in.csv", rows)
    common = [path, "--series", SERIES, "--group-column", "group", "--draws", "100"]
    bounds = {name: core.place_mics(ts, LADDER)[:2] for name, ts in texts.items()}
    rows_a = np.arange(150)

    fit, _ = run(capsys, "fit", *common)
    mean = core.resample_clusters(*bounds["b"], rows_a,
                                  core.named_generator(20260101, "fit|b"), draws=100)
    sd = core.resample_clusters(*bounds["b"], rows_a,
                                core.named_generator(20260101, "fit|b"), draws=100,
                                statistic=cli.sd_of)
    row = by_group(fit)["b"]
    assert row["mean_low"] == f"{mean['low']:.4f}"
    assert row["mean_high"] == f"{mean['high']:.4f}"
    assert row["sd_low"] == f"{sd['low']:.4f}" and row["sd_high"] == f"{sd['high']:.4f}"

    shift, _ = run(capsys, "shift", *common, "--reference", "a")
    interval = core.resample_shift(*bounds["b"], rows_a, *bounds["a"], rows_a,
                                   core.named_generator(20260101, "shift|b|a"),
                                   draws=100)
    assert shift[0]["shift_low"] == f"{interval['low']:.4f}"
    assert shift[0]["shift_high"] == f"{interval['high']:.4f}"

    ecoff, _ = run(capsys, "ecoff", *common, "--coverage", "0.95")
    quantile = core.resample_clusters(*bounds["b"], rows_a,
                                      core.named_generator(20260101, "ecoff|b"),
                                      draws=100, statistic=cli.quantile_of(0.95))
    assert by_group(ecoff)["b"]["quantile_low"] == f"{quantile['low']:.4f}"


@pytest.mark.parametrize("command,extra,target", [
    ("fit", [], "resample_clusters"),
    ("shift", ["--reference", "a"], "resample_shift"),
    ("ecoff", ["--coverage", "0.99"], "resample_clusters"),
])
def test_an_interval_that_could_not_be_formed_withholds_the_estimate(
        tmp_path, capsys, monkeypatch, command, extra, target):
    path = planted(tmp_path / "in.csv", {"a": (-5.0, 1.0, 100), "b": (-4.0, 1.0, 100)})
    monkeypatch.setattr(core, target, lambda *args, **kwargs: None)
    out, _ = run(capsys, command, path, "--series", SERIES, "--group-column",
                 "group", "--draws", "100", *extra)
    for row in out:
        assert row["withheld"] == cli.DRAWS_FAILED
        assert row["draws"] == ""


def test_two_runs_write_identical_output(tmp_path, capsys):
    path = planted(tmp_path / "in.csv", {"a": (-5.0, 1.0, 200), "b": (-4.0, 1.0, 200)})
    argv = ["shift", path, "--series", SERIES, "--group-column", "group",
            "--reference", "a", "--draws", "100"]
    cli.main(argv)
    first = capsys.readouterr().out
    cli.main(argv)
    assert capsys.readouterr().out == first


def test_the_seed_moves_the_interval_and_leaves_the_estimate(tmp_path, capsys):
    path = planted(tmp_path / "in.csv", {"a": (-5.0, 1.0, 200)})
    one, _ = run(capsys, "fit", path, "--series", SERIES, "--draws", "100")
    two, _ = run(capsys, "fit", path, "--series", SERIES, "--draws", "100",
                 "--seed", "7")
    assert one[0]["mean_log2"] == two[0]["mean_log2"]
    assert one[0]["mean_low"] != two[0]["mean_low"]


# Readings that do not fit the series


def test_a_reading_off_the_series_stops_the_run(tmp_path, capsys):
    rows = [("0.06", "a", "c"), ("0.2", "a", "c"), ("0.12", "a", "c")]
    message = stopped(capsys, "fit", write_csv(tmp_path / "in.csv", rows),
                      "--series", SERIES)
    assert "1 readings do not match the tested series" in message
    assert "line 3: '0.2'" in message


@pytest.mark.parametrize("reading", ["<=0.015", ">2", ">=0.5", "<0.008"])
def test_a_censored_reading_from_another_plate_design_stops_the_run(tmp_path,
                                                                     capsys,
                                                                     reading):
    """
    A plate testing 0.008 to 1 reports <=0.008 and >1 at its ends. <=0.015 and
    >2 came from a plate testing a different series.
    """
    rows = [("0.06", "a", "c"), (reading, "a", "c"), ("<=0.008", "a", "c"),
            (">1", "a", "c")]
    message = stopped(capsys, "fit", write_csv(tmp_path / "in.csv", rows),
                      "--series", SERIES)
    assert "1 readings do not match the tested series" in message
    assert f"line 3: {reading!r}" in message


def test_censored_readings_at_the_ends_of_the_series_are_placed(tmp_path, capsys):
    rows = [("<=0.008", "a", "c")] * 20 + [("0.06", "a", "c")] * 20
    rows += [(">1", "a", "c")] * 20 + [(">=1", "a", "c")] * 5
    out, _ = run(capsys, "fit", write_csv(tmp_path / "in.csv", rows),
                 "--series", SERIES, "--draws", "100")
    assert out[0]["placed"] == "65" and out[0]["off_series"] == "0"


# Malformed input


def test_a_missing_mic_column_names_the_columns_present(tmp_path, capsys):
    path = write_csv(tmp_path / "in.csv", [("0.06", "a", "c")],
                     header=("reading", "group", "cluster"))
    message = stopped(capsys, "fit", path, "--series", SERIES)
    assert "no column 'MIC'" in message and "'reading'" in message


def test_the_mic_column_can_be_named(tmp_path, capsys):
    rows = [(t, "a", "c") for t in reports(-5.0, 1.0, 100, 10)]
    path = write_csv(tmp_path / "in.csv", rows, header=("BDQ_MIC", "group", "c"))
    out, _ = run(capsys, "fit", path, "--series", SERIES, "--mic-column",
                 "BDQ_MIC", "--draws", "100")
    assert out[0]["placed"] == "100"


@pytest.mark.parametrize("label", ["", "  "])
def test_a_row_with_no_group_label_stops_the_run(tmp_path, capsys, label):
    rows = [("0.06", "a", "c"), ("0.12", label, "c")]
    message = stopped(capsys, "fit", write_csv(tmp_path / "in.csv", rows),
                      "--series", SERIES, "--group-column", "group")
    assert "1 rows carry no group label" in message and "lines 3" in message


def test_a_row_with_no_cluster_label_stops_the_run(tmp_path, capsys):
    rows = [("0.06", "a", "c"), ("0.12", "a", "")]
    message = stopped(capsys, "fit", write_csv(tmp_path / "in.csv", rows),
                      "--series", SERIES, "--cluster-column", "cluster")
    assert "1 rows carry no cluster label" in message


def test_labels_are_kept_exactly_as_written(tmp_path, capsys):
    rows = [(t, "A", "c") for t in reports(-5.0, 1.0, 60, 11)]
    rows += [(t, " A", "c") for t in reports(-4.0, 1.0, 60, 12)]
    out, _ = run(capsys, "fit", write_csv(tmp_path / "in.csv", rows), "--series",
                 SERIES, "--group-column", "group", "--draws", "100")
    assert sorted(row["group"] for row in out) == [" A", "A"]


def test_a_repeated_column_name_stops_the_run(tmp_path, capsys):
    path = write_csv(tmp_path / "in.csv", [("0.06", "a", "b")],
                     header=("MIC", "group", "MIC"))
    assert "the header repeats 'MIC'" in stopped(capsys, "fit", path,
                                                 "--series", SERIES)


def test_a_row_of_the_wrong_length_stops_the_run(tmp_path, capsys):
    path = tmp_path / "in.csv"
    path.write_text("MIC,group,cluster\n0.06,a,c\n0.12,a\n")
    assert "line 3 holds 2 fields where the header holds 3" in stopped(
        capsys, "fit", str(path), "--series", SERIES)


def test_an_empty_file_and_a_header_alone_stop_the_run(tmp_path, capsys):
    empty = tmp_path / "empty.csv"
    empty.write_text("")
    assert "the input is empty" in stopped(capsys, "fit", str(empty),
                                           "--series", SERIES)
    header = tmp_path / "header.csv"
    header.write_text("MIC,group,cluster\n")
    assert "a header and no rows" in stopped(capsys, "fit", str(header),
                                             "--series", SERIES)


def test_an_unreadable_path_stops_the_run(tmp_path, capsys):
    assert "cannot read" in stopped(capsys, "fit", str(tmp_path / "absent.csv"),
                                    "--series", SERIES)


def test_a_byte_order_mark_and_blank_lines_are_tolerated(tmp_path, capsys):
    rows = [(t, "a", "c") for t in reports(-5.0, 1.0, 100, 13)]
    path = write_csv(tmp_path / "in.csv", rows, encoding="utf-8-sig")
    with open(path, "a", encoding="utf-8") as handle:
        handle.write("\n\n")
    out, _ = run(capsys, "fit", path, "--series", SERIES, "--draws", "100")
    assert out[0]["placed"] == "100" and out[0]["rows"] == "100"


def test_a_byte_order_mark_on_standard_input_is_tolerated(capsys, monkeypatch):
    text = "\ufeffMIC,group,cluster\n" + "".join(
        f"{t},a,c\n" for t in reports(-5.0, 1.0, 100, 15))
    monkeypatch.setattr(sys, "stdin", io.StringIO(text))
    out, _ = run(capsys, "fit", "-", "--series", SERIES, "--draws", "100")
    assert out[0]["placed"] == "100"


def test_input_can_come_from_standard_input(tmp_path, capsys, monkeypatch):
    rows = [(t, "a", "c") for t in reports(-5.0, 1.0, 100, 14)]
    path = write_csv(tmp_path / "in.csv", rows)
    from_file, _ = run(capsys, "fit", path, "--series", SERIES, "--draws", "100")
    monkeypatch.setattr(sys, "stdin", io.StringIO(open(path).read()))
    from_stdin, _ = run(capsys, "fit", "-", "--series", SERIES, "--draws", "100")
    assert from_stdin == from_file


# The series


@pytest.mark.parametrize("series,message", [
    ("<=0.008,0.015,0.03", "with no operator"),
    ("0.008,0.015,0.015,0.03", "the series repeats 0.015"),
    ("0,0.5,1", "not a concentration: '0'"),
    ("-0.5,0.5,1", "not a concentration: '-0.5'"),
    ("0.5,abc,1", "not a concentration: 'abc'"),
    ("0.5,inf", "not a concentration: 'inf'"),
    ("0.5", "at least two concentrations"),
    ("0.12,0.125,0.25", "too close"),
])
def test_a_malformed_series_stops_the_run(tmp_path, capsys, series, message):
    """
    The series is passed as --series=value, which every argparse version reads
    as the option's value even where the value begins with a minus sign.
    """
    path = write_csv(tmp_path / "in.csv", [("0.5", "a", "c")])
    assert message in stopped(capsys, "fit", path, f"--series={series}")


def test_a_series_with_a_skipped_dilution_warns_and_spans_the_gap(tmp_path, capsys):
    rows = [("0.008", "a", "c")] * 10 + [("0.015", "a", "c")] * 10
    rows += [("0.06", "a", "c")] * 10 + [("0.12", "a", "c")] * 10
    path = write_csv(tmp_path / "in.csv", rows)
    out, err = run(capsys, "fit", path, "--series", "0.008,0.015,0.06,0.12",
                   "--draws", "100")
    assert "not a doubling series" in err
    lower, upper, _, _ = core.place_mics(["0.06"], [0.008, 0.015, 0.06, 0.12])
    assert lower[0] == np.log2(0.015)
    assert out[0]["placed"] == "40"


def test_the_series_may_be_given_in_any_order(tmp_path, capsys):
    path = planted(tmp_path / "in.csv", {"a": (-5.0, 1.0, 100)})
    forward, _ = run(capsys, "fit", path, "--series", SERIES, "--draws", "100")
    backward, _ = run(capsys, "fit", path, "--series",
                      ",".join(reversed(SERIES.split(","))), "--draws", "100")
    assert forward == backward


@pytest.mark.parametrize("option,value", [("--draws", "99"), ("--draws", "x"),
                                          ("--seed", "-1"),
                                          ("--min-clusters", "1")])
def test_resampling_options_are_bounded(tmp_path, capsys, option, value):
    path = planted(tmp_path / "in.csv", {"a": (-5.0, 1.0, 50)})
    assert refused(capsys, "fit", path, "--series", SERIES, option, value) == 2


# Running and packaging


def run_module(*argv, cwd):
    environment = dict(os.environ, PYTHONPATH=str(ROOT), PYTHONDONTWRITEBYTECODE="1")
    return subprocess.run([sys.executable, "-m", "micecoff", *argv], cwd=cwd,
                          env=environment, capture_output=True, text=True,
                          timeout=120)


def test_the_package_runs_as_a_module(tmp_path):
    version = run_module("--version", cwd=tmp_path)
    assert version.returncode == 0
    assert version.stdout.strip() == f"micecoff {__version__}"
    path = planted(tmp_path / "in.csv", {"a": (-5.0, 1.0, 100)})
    result = run_module("fit", path, "--series", SERIES, "--draws", "100",
                        cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines()[0].startswith("group,rows,")


def test_a_stopped_run_exits_nonzero_with_its_message_on_stderr(tmp_path):
    path = write_csv(tmp_path / "in.csv", [("0.2", "a", "c")])
    result = run_module("fit", path, "--series", SERIES, cwd=tmp_path)
    assert result.returncode == 1
    assert result.stdout == ""
    assert result.stderr.startswith("micecoff: 1 readings do not match")


def project():
    tomllib = pytest.importorskip("tomllib")
    return tomllib.loads((ROOT / "pyproject.toml").read_text())


def test_the_console_script_names_the_interface_entry_point():
    scripts = project()["project"]["scripts"]
    assert scripts == {"micecoff": "micecoff.cli:main"}
    module, _, function = scripts["micecoff"].partition(":")
    assert module == cli.__name__ and getattr(cli, function) is cli.main


def test_the_version_is_read_from_the_package():
    configuration = project()
    assert "version" not in configuration["project"]
    assert "version" in configuration["project"]["dynamic"]
    assert configuration["tool"]["setuptools"]["dynamic"]["version"] == {
        "attr": "micecoff.__version__"}
    assert re.fullmatch(r"\d+\.\d+\.\d+", __version__)


def test_a_license_expression_is_built_by_a_setuptools_that_reads_one():
    """
    A license given as an SPDX expression string, with license-files, is a
    project table setuptools 76.1.0 refuses and 77.0.1 builds, measured on
    this project with the floor removed.
    """
    configuration = project()
    assert isinstance(configuration["project"]["license"], str)
    floors = [re.fullmatch(r"setuptools>=(\d+)(?:\.\d+)*", item)
              for item in configuration["build-system"]["requires"]]
    floor = [int(match.group(1)) for match in floors if match]
    assert floor and floor[0] >= 77
    for name in configuration["project"]["license-files"]:
        assert (ROOT / name).is_file()


def declared():
    names = {}
    for item in project()["project"]["dependencies"]:
        match = re.fullmatch(r"([A-Za-z0-9_.-]+)>=([\d.]+)", item)
        assert match, f"a dependency must be a plain lower bound: {item}"
        names[match.group(1).lower()] = match.group(2)
    return names


def version_tuple(text):
    return tuple(int(part) for part in text.split("."))


def test_every_third_party_import_in_the_package_is_declared():
    imported = set()
    for path in (ROOT / "micecoff").glob("*.py"):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Import):
                imported.update(alias.name.split(".")[0] for alias in node.names)
            if isinstance(node, ast.ImportFrom) and node.level == 0:
                imported.add(node.module.split(".")[0])
    third_party = imported - set(sys.stdlib_module_names) - {"micecoff"}
    assert third_party == set(declared())


def test_the_pinned_environment_satisfies_every_declared_floor():
    pins = {}
    for line in (ROOT / "requirements.txt").read_text().splitlines():
        match = re.fullmatch(r"([A-Za-z0-9_.-]+)==([\d.]+)", line.strip())
        if match:
            pins[match.group(1).lower()] = match.group(2)
    for name, floor in declared().items():
        assert name in pins, f"{name} is declared and not pinned"
        assert version_tuple(pins[name]) >= version_tuple(floor)


def test_the_build_includes_the_package_alone():
    assert project()["tool"]["setuptools"]["packages"] == ["micecoff"]


def test_the_source_distribution_leaves_out_the_atlas_readme_and_tests():
    """
    setuptools adds README.md and tests/test*.py to a source distribution by
    default. The manifest takes both out, measured by listing a built one.
    """
    commands = [line.split() for line in
                (ROOT / "MANIFEST.in").read_text().splitlines()
                if line.strip() and not line.startswith("#")]
    assert ["exclude", "README.md"] in commands
    assert ["prune", "tests"] in commands
