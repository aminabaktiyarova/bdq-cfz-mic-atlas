"""
Command line interface to the interval-censored MIC estimator.

Three subcommands read a CSV holding one reported MIC per row and write a CSV
to standard output:

  fit     the fitted log2 MIC distribution of each group
  shift   each group's shift against a reference group, in doublings
  ecoff   the cut-off at a chosen quantile of a fitted wild-type
          distribution, as a tested concentration

The tested concentration series is supplied with --series. A well that no
isolate reported still bounds the interval below the well above it, so the
series cannot be read off the data.

Every estimate carries a 95% interval from resampling with
replacement: clusters when --cluster-column names one, single rows otherwise,
which treats the rows as independent. Each group's interval draws from its own
generator, seeded by --seed and the group's name, so a group's interval does
not depend on which other groups the file holds.

An estimate is withheld, with its reason in the withheld column and its counts
kept, when the group holds fewer resampling units than --min-clusters, when
core.withhold_reason finds the fit unsupported, or when fewer than a quarter
of the resampling draws converge.

The ecoff subcommand fits one normal distribution to every placed reading in a
group, so each group must hold a wild-type population, defined for example by
genotype. Non-wild-type isolates in a group widen the fit and raise the
cut-off.

Input the tool cannot interpret stops it with a message: a reading that does
not match the series, a censored reading whose bound is not the end of the
series, a malformed series, a repeated column name, a row whose field count
differs from the header's, and a row with no group or cluster label.
--allow-off-series drops the readings in the first two classes, counts them
per group, and reports the total.
"""

import argparse
import csv
import sys

import numpy as np
from scipy.stats import norm

from . import __version__, core

INTERVAL = (2.5, 97.5)
MINIMUM_DRAWS = 100


def stop(message):
    raise SystemExit(f"micecoff: {message}")


def note(message):
    print(f"micecoff: {message}", file=sys.stderr)


# Reading the input


def read_rows(path):
    """The header, the rows and each row's line number, from a path or stdin."""
    if path == "-":
        handle = sys.stdin
    else:
        try:
            handle = open(path, newline="", encoding="utf-8")
        except OSError as error:
            stop(f"cannot read {path}: {error.strerror}")
    try:
        reader = csv.reader(handle)
        header = next(reader, None)
        rows, lines = [], []
        for row in reader:
            if row:
                rows.append(row)
                lines.append(reader.line_num)
    finally:
        if handle is not sys.stdin:
            handle.close()

    if not header:
        stop("the input is empty")
    # A byte order mark, which spreadsheet software writes at the start of a
    # CSV, reaches the first column name from a file and from standard input.
    if header[0].startswith("\ufeff"):
        header[0] = header[0][1:]
    repeated = sorted({name for name in header if header.count(name) > 1})
    if repeated:
        stop(f"the header repeats {', '.join(map(repr, repeated))}")
    if not rows:
        stop("the input holds a header and no rows")
    for row, line in zip(rows, lines):
        if len(row) != len(header):
            stop(f"line {line} holds {len(row)} fields where the header "
                 f"holds {len(header)}")
    return header, rows, lines


def column(header, rows, name, what):
    if name not in header:
        stop(f"no column {name!r} for the {what}; the header holds "
             f"{', '.join(map(repr, header))}")
    index = header.index(name)
    return [row[index] for row in rows]


def labels(header, rows, lines, name, what):
    """A label column in which every row must carry a label."""
    values = column(header, rows, name, what)
    empty = [line for value, line in zip(values, lines) if not value.strip()]
    if empty:
        shown = ", ".join(str(line) for line in empty[:5])
        stop(f"{len(empty)} rows carry no {what} label in column {name!r} "
             f"(lines {shown})")
    return values


def parse_series(text):
    """A comma-separated series of tested concentrations, sorted and checked."""
    values = []
    for item in text.split(","):
        item = item.strip()
        value, operator = core.parse_concentration(item)
        if operator is not None:
            stop(f"the series takes tested concentrations with no operator; "
                 f"got {item!r}")
        if value is None or not np.isfinite(value) or value <= 0:
            stop(f"not a concentration: {item!r}")
        values.append(value)
    if len(values) < 2:
        stop("the series needs at least two concentrations")
    repeated = sorted({value for value in values if values.count(value) > 1})
    if repeated:
        stop(f"the series repeats {', '.join(f'{v:g}' for v in repeated)}")
    series = sorted(values)

    # A reading matches a tested concentration within core.MATCH_TOLERANCE, so
    # two concentrations closer than this ratio could both match one reading.
    separable = (1 + core.MATCH_TOLERANCE) / (1 - core.MATCH_TOLERANCE)
    ratios, regular = core.check_doubling_series(series)
    for ratio, below, above in zip(ratios, series, series[1:]):
        if ratio < separable:
            stop(f"the series holds {below:g} and {above:g}, which are too "
                 "close for a reading to be matched to one of them")
    if not regular:
        note("warning: the series is not a doubling series, so the interval "
             "below a reading spans the gap to the next tested concentration; "
             f"ratios {', '.join(f'{r:.3f}' for r in ratios)}")
    return series


def matches(value, concentration):
    return abs(concentration - value) / value < core.MATCH_TOLERANCE


def place(readings, lines, series, allow_off_series):
    """
    Censoring intervals for every reading, with absent readings and off-series
    readings marked. A reading enters a fit only where it is neither.

    A censored reading carries its own bound, and on a plate testing this
    series that bound is the end of the series: <= at the lowest
    concentration, > or >= at the highest. A censored reading at any other
    concentration came from a plate testing a different series, so it counts
    as off the series.
    """
    lower, upper, _, _ = core.place_mics(readings, series)
    absent = np.array([not core.mic_recorded(text) for text in readings])
    off = ~absent & np.isnan(lower)
    for index, text in enumerate(readings):
        if absent[index] or off[index]:
            continue
        value, operator = core.parse_concentration(text)
        if operator == "<=" and not matches(value, series[0]):
            off[index] = True
        if operator in (">", ">=") and not matches(value, series[-1]):
            off[index] = True

    if off.any() and not allow_off_series:
        shown = "; ".join(f"line {lines[i]}: {readings[i]!r}"
                          for i in np.flatnonzero(off)[:5])
        stop(f"{int(off.sum())} readings do not match the tested series "
             f"({shown}). The series does not describe this data. Correct the "
             "series, split the input by plate design, or pass "
             "--allow-off-series to drop them.")
    return lower, upper, absent, off


def load(args):
    """Everything a subcommand needs from the input file."""
    header, rows, lines = read_rows(args.input)
    readings = column(header, rows, args.mic_column, "MIC")
    groups = (labels(header, rows, lines, args.group_column, "group")
              if args.group_column else ["all"] * len(rows))
    clusters = (labels(header, rows, lines, args.cluster_column, "cluster")
                if args.cluster_column else None)
    series = parse_series(args.series)
    lower, upper, absent, off = place(readings, lines, series,
                                      args.allow_off_series)
    groups = np.array(groups)
    indices = {name: np.flatnonzero(groups == name)
               for name in sorted(set(groups.tolist()))}
    if int(absent.sum()):
        note(f"note: {int(absent.sum())} rows carry no reading and enter no fit")
    if int(off.sum()):
        note(f"note: {int(off.sum())} readings off the series were dropped "
             "under --allow-off-series")
    return {"series": series, "lower": lower, "upper": upper,
            "absent": absent, "off": off, "usable": ~absent & ~off,
            "groups": indices,
            "clusters": None if clusters is None else np.array(clusters)}


# Estimating


def counts(data, index):
    placed = index[data["usable"][index]]
    return placed, {
        "rows": len(index),
        "absent": int(data["absent"][index].sum()),
        "off_series": int(data["off"][index].sum()),
        "placed": len(placed),
    }


def units(data, placed):
    """The resampling unit of each placed row: its cluster, or the row."""
    if data["clusters"] is None:
        return placed.copy()
    return data["clusters"][placed]


def unit_name(data):
    return "row" if data["clusters"] is None else "cluster"


def assess(data, placed, minimum):
    """The fit of one group's placed readings, and why it is withheld if it is."""
    lower, upper = data["lower"][placed], data["upper"][placed]
    fit = core.fit_censored_normal(lower, upper)
    reason = core.withhold_reason(fit, lower, upper, data["series"])
    available = len(set(units(data, placed).tolist()))
    if reason is None and available < minimum:
        reason = (f"too few resampling units ({available} of the {minimum} "
                  "an interval needs)")
    return fit, reason


def spread(data, placed, rng, draws, statistic=None):
    return core.resample_clusters(data["lower"][placed], data["upper"][placed],
                                  units(data, placed), rng, statistic=statistic,
                                  draws=draws, interval=INTERVAL)


def sd_of(lower, upper):
    fit = core.fit_censored_normal(lower, upper)
    return None if fit is None else fit["sigma"]


def quantile_of(coverage):
    def statistic(lower, upper):
        fit = core.fit_censored_normal(lower, upper)
        if fit is None:
            return None
        return float(norm.ppf(coverage, loc=fit["mu"], scale=fit["sigma"]))
    return statistic


DRAWS_FAILED = "fewer than a quarter of the resampling draws converged"


def log2_text(value):
    return f"{value:.4f}"


def concentration_text(value):
    return f"{value:.4g}"


def well_text(log2_value, series):
    """The tested concentration a cut-off rounds up to, or >c above the series."""
    concentration = core.round_up_to_series(log2_value, series)
    return f">{series[-1]:g}" if concentration is None else f"{concentration:g}"


FIT_COLUMNS = ["group", "rows", "absent", "off_series", "placed",
               "left_censored", "right_censored", "resampling_unit", "clusters",
               "mean_log2", "mean_low", "mean_high", "sd_log2", "sd_low",
               "sd_high", "geometric_mean", "geometric_mean_low",
               "geometric_mean_high", "draws", "withheld"]

SHIFT_COLUMNS = ["group", "reference", "rows", "absent", "off_series", "placed",
                 "reference_placed", "resampling_unit", "clusters",
                 "reference_clusters", "mean_log2", "reference_mean_log2",
                 "shift", "shift_low", "shift_high", "draws", "withheld"]

ECOFF_COLUMNS = ["group", "rows", "absent", "off_series", "placed",
                 "resampling_unit", "clusters", "mean_log2", "sd_log2",
                 "coverage", "quantile_log2", "quantile_low", "quantile_high",
                 "ecoff", "ecoff_low", "ecoff_high", "draws", "withheld"]


def blank(columns):
    return {name: "" for name in columns}


def command_fit(args):
    data = load(args)
    out = []
    for name, index in data["groups"].items():
        placed, row = counts(data, index)
        result = blank(FIT_COLUMNS)
        result.update(row, group=name, resampling_unit=unit_name(data),
                      clusters=len(set(units(data, placed).tolist())),
                      left_censored=int(np.isneginf(data["lower"][placed]).sum()),
                      right_censored=int(np.isposinf(data["upper"][placed]).sum()))
        fit, reason = assess(data, placed, args.min_clusters)
        if reason is None:
            # One stream name for both passes, so the mean and the standard
            # deviation are read from the same resampled datasets.
            stream = f"fit|{name}"
            mean = spread(data, placed, core.named_generator(args.seed, stream),
                          args.draws)
            sd = spread(data, placed, core.named_generator(args.seed, stream),
                        args.draws, statistic=sd_of)
            if mean is None or sd is None:
                reason = DRAWS_FAILED
        if reason is None:
            result.update(
                mean_log2=log2_text(fit["mu"]),
                mean_low=log2_text(mean["low"]), mean_high=log2_text(mean["high"]),
                sd_log2=log2_text(fit["sigma"]),
                sd_low=log2_text(sd["low"]), sd_high=log2_text(sd["high"]),
                geometric_mean=concentration_text(2 ** fit["mu"]),
                geometric_mean_low=concentration_text(2 ** mean["low"]),
                geometric_mean_high=concentration_text(2 ** mean["high"]),
                draws=mean["draws"])
        else:
            result["withheld"] = reason
        out.append(result)
    write(out, FIT_COLUMNS)


def command_shift(args):
    data = load(args)
    if args.reference not in data["groups"]:
        stop(f"no group {args.reference!r} in column {args.group_column!r}; "
             f"the column holds {', '.join(map(repr, data['groups']))}")
    if len(data["groups"]) < 2:
        stop(f"column {args.group_column!r} holds only the reference group "
             f"{args.reference!r}")

    reference, _ = counts(data, data["groups"][args.reference])
    reference_fit, reference_reason = assess(data, reference, args.min_clusters)

    out = []
    for name, index in data["groups"].items():
        if name == args.reference:
            continue
        placed, row = counts(data, index)
        result = blank(SHIFT_COLUMNS)
        result.update(row, group=name, reference=args.reference,
                      reference_placed=len(reference),
                      resampling_unit=unit_name(data),
                      clusters=len(set(units(data, placed).tolist())),
                      reference_clusters=len(set(units(data, reference).tolist())))
        fit, reason = assess(data, placed, args.min_clusters)
        if reference_reason is not None:
            reason = f"reference: {reference_reason}"
        if reason is None:
            interval = core.resample_shift(
                data["lower"][placed], data["upper"][placed], units(data, placed),
                data["lower"][reference], data["upper"][reference],
                units(data, reference),
                core.named_generator(args.seed, f"shift|{name}|{args.reference}"),
                draws=args.draws, interval=INTERVAL)
            if interval is None:
                reason = DRAWS_FAILED
        if reason is None:
            result.update(
                mean_log2=log2_text(fit["mu"]),
                reference_mean_log2=log2_text(reference_fit["mu"]),
                shift=log2_text(fit["mu"] - reference_fit["mu"]),
                shift_low=log2_text(interval["low"]),
                shift_high=log2_text(interval["high"]),
                draws=interval["draws"])
        else:
            result["withheld"] = reason
        out.append(result)
    write(out, SHIFT_COLUMNS)


def command_ecoff(args):
    data = load(args)
    series = data["series"]
    out = []
    for name, index in data["groups"].items():
        placed, row = counts(data, index)
        result = blank(ECOFF_COLUMNS)
        result.update(row, group=name, resampling_unit=unit_name(data),
                      clusters=len(set(units(data, placed).tolist())),
                      coverage=f"{args.coverage:g}")
        fit, reason = assess(data, placed, args.min_clusters)
        if reason is None:
            interval = spread(data, placed,
                              core.named_generator(args.seed, f"ecoff|{name}"),
                              args.draws, statistic=quantile_of(args.coverage))
            if interval is None:
                reason = DRAWS_FAILED
        if reason is None:
            cut = core.ecoff(fit, series, args.coverage)
            result.update(
                mean_log2=log2_text(fit["mu"]), sd_log2=log2_text(fit["sigma"]),
                quantile_log2=log2_text(cut["quantile"]),
                quantile_low=log2_text(interval["low"]),
                quantile_high=log2_text(interval["high"]),
                ecoff=well_text(cut["quantile"], series),
                ecoff_low=well_text(interval["low"], series),
                ecoff_high=well_text(interval["high"], series),
                draws=interval["draws"])
        else:
            result["withheld"] = reason
        out.append(result)
    write(out, ECOFF_COLUMNS)


def write(rows, columns):
    writer = csv.DictWriter(sys.stdout, fieldnames=columns, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)


# Arguments


def bounded_int(minimum, what):
    def convert(text):
        try:
            value = int(text)
        except ValueError:
            raise argparse.ArgumentTypeError(f"{what} must be an integer")
        if value < minimum:
            raise argparse.ArgumentTypeError(f"{what} must be at least {minimum}")
        return value
    return convert


def proportion(text):
    try:
        value = float(text)
    except ValueError:
        raise argparse.ArgumentTypeError("coverage must be a number")
    if not 0 < value < 1:
        raise argparse.ArgumentTypeError(
            "coverage must lie strictly between 0 and 1")
    return value


def build_parser():
    parser = argparse.ArgumentParser(
        prog="micecoff",
        description="Interval-censored MIC estimation and cut-off derivation.")
    parser.add_argument("--version", action="version",
                        version=f"micecoff {__version__}")
    subcommands = parser.add_subparsers(dest="command", required=True)

    def shared(sub, group_required=False):
        sub.add_argument("input", help="CSV path, or - for standard input")
        sub.add_argument("--series", required=True,
                         help="tested concentrations, comma separated, for "
                              "example 0.008,0.015,0.03,0.06,0.12")
        sub.add_argument("--mic-column", default="MIC",
                         help="column holding the reported MIC (default MIC)")
        sub.add_argument("--group-column", required=group_required,
                         help="column holding the group label")
        sub.add_argument("--cluster-column",
                         help="column holding a cluster label; intervals then "
                              "resample clusters, and single rows otherwise")
        sub.add_argument("--draws", type=bounded_int(MINIMUM_DRAWS, "draws"),
                         default=400, help="resampling draws (default 400)")
        sub.add_argument("--seed", type=bounded_int(0, "the seed"),
                         default=20260101,
                         help="base seed for the resampling (default 20260101)")
        sub.add_argument("--min-clusters", type=bounded_int(2, "min-clusters"),
                         default=5,
                         help="resampling units a group needs before an "
                              "estimate is reported (default 5)")
        sub.add_argument("--allow-off-series", action="store_true",
                         help="drop readings that do not match the series, "
                              "and count them, instead of stopping")

    fit = subcommands.add_parser(
        "fit", help="the fitted log2 MIC distribution of each group")
    shared(fit)
    fit.set_defaults(handler=command_fit)

    shift = subcommands.add_parser(
        "shift", help="each group's shift against a reference group, in doublings")
    shared(shift, group_required=True)
    shift.add_argument("--reference", required=True,
                       help="the label of the reference group")
    shift.set_defaults(handler=command_shift)

    ecoff = subcommands.add_parser(
        "ecoff", help="the cut-off at a quantile of a fitted wild-type "
                      "distribution")
    shared(ecoff)
    ecoff.add_argument("--coverage", type=proportion, required=True,
                       help="the quantile of the fitted wild-type distribution "
                            "the cut-off must reach, strictly between 0 and 1. "
                            "There is no default: the convention is set by the "
                            "EUCAST standard operating procedure for MIC "
                            "distributions and cut-off setting")
    ecoff.set_defaults(handler=command_ecoff)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.handler(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
