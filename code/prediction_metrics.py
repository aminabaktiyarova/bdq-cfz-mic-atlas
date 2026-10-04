"""
What a genotype rule would achieve as a test for resistance at the ECOFF.

Run from the project root with the virtual environment active:

    python code/prediction_metrics.py

The project reports odds ratios and MIC shifts. Those say how far a variant
moves the MIC; they do not say how much of the resistance in a collection a rule
would catch, or what it would call resistant that is not. This module reports
sensitivity, specificity and the predictive values for four rules, with
intervals that resample clusters.

The rules, each applied per drug to the samples carrying a MIC for it:

  major variant       a real major-allele variant in Rv0678, pepQ or atpE
  Rv0678 any          a real major-allele variant in Rv0678, any class
  Rv0678 loss         an Rv0678 frameshift, stop codon or gene deletion
  major or minor      the major-variant rule, or a detected minor allele in any
                      of the three genes

A rule calls an isolate resistant or not. A null call is not a detected variant,
so a sample whose gene could not be read is called not-resistant by every rule,
which is how a catalogue applied to a real sequence behaves. The fourth rule
exists because Section 8.8 finds detected minor alleles associated with MIC
shifts as large as major alleles of the same class, and the question is what
admitting them would cost.

Predictive values depend on how much resistance the collection holds. Section
5.3 records that resistance here is concentrated at one site by design, so the
positive predictive value below describes this collection and transfers to no
other. Sensitivity and specificity are less exposed to that, and they still
describe this collection's mix of variants.

Intervals resample clusters, as everywhere else in the project, because isolates
sharing a site, a sublineage and a mutation are one event observed many times.

Outputs:
  outputs/prediction_metrics.txt
  outputs/prediction_metrics.csv
  outputs/prediction_thresholds.csv

The ECOFF is one cut-off among the tested concentrations, and a rule that looks
weak against it may look different against another. The same four rules are
therefore also evaluated at every concentration shared by both plate designs,
which is where each isolate's position relative to the cut-off is determinate.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import cohort  # noqa: E402
from cluster_adjust import add_clusters  # noqa: E402
from mic_model import concentration_series, place_mics  # noqa: E402

REPORT = Path("outputs/prediction_metrics.txt")
TABLE = Path("outputs/prediction_metrics.csv")
THRESHOLDS = Path("outputs/prediction_thresholds.csv")

BOOTSTRAPS = 400
SEED = 20260101
RULES = ["major variant", "Rv0678 any", "Rv0678 loss", "major or minor"]
# The concentration the ECOFF sits at, in mg/L, for both plate designs. An
# isolate is resistant when its MIC exceeds it, which on a doubling ladder means
# reaching the next rung up.
ECOFF = {"BDQ": 0.25, "CFZ": 0.25}

_lines = []


def say(text=""):
    print(text)
    _lines.append(text)


def stream(name):
    """A generator seeded by name.

    Each interval draws from its own stream, so a figure does not depend on
    which figures were computed before it. The ECOFF row of the sweep and the
    row for the same rule in the headline table are keyed identically, so the
    two tables carry one interval for one quantity rather than two.
    """
    import numpy as np

    return np.random.default_rng([SEED] + [ord(letter) for letter in name])


def rule_flags(mutations, index):
    """One boolean column per rule, for every sample in the index."""
    import pandas as pd

    target = mutations[mutations.GENE.isin(cohort.BDQ_GENES)]
    major = target[target.REAL_MAJOR]
    rv = major[major.GENE.eq("Rv0678")]
    minor = target[target.IS_HET_CALL | target.IS_MINOR]

    flags = pd.DataFrame(index=index)
    flags["major variant"] = index.isin(set(major.UNIQUEID))
    flags["Rv0678 any"] = index.isin(set(rv.UNIQUEID))
    flags["Rv0678 loss"] = index.isin(set(rv[rv.IS_LOF].UNIQUEID))
    flags["major or minor"] = flags["major variant"] | index.isin(set(minor.UNIQUEID))
    return flags


def confusion(called, resistant):
    """Counts of the four cells, as plain integers."""
    return {
        "true_positive": int((called & resistant).sum()),
        "false_positive": int((called & ~resistant).sum()),
        "false_negative": int((~called & resistant).sum()),
        "true_negative": int((~called & ~resistant).sum()),
    }


def metrics(counts):
    """Sensitivity, specificity and the predictive values. A ratio with no
    denominator is None rather than zero."""
    def ratio(numerator, denominator):
        return numerator / denominator if denominator else None

    positives = counts["true_positive"] + counts["false_negative"]
    negatives = counts["false_positive"] + counts["true_negative"]
    called = counts["true_positive"] + counts["false_positive"]
    not_called = counts["false_negative"] + counts["true_negative"]
    return {
        "sensitivity": ratio(counts["true_positive"], positives),
        "specificity": ratio(counts["true_negative"], negatives),
        "ppv": ratio(counts["true_positive"], called),
        "npv": ratio(counts["true_negative"], not_called),
    }


def cluster_cells(frame, called_column, resistant_column):
    """The four confusion cells for each cluster, as one row per cluster.

    A cluster enters a resample whole, so its isolates only ever contribute
    together. Holding the counts per cluster makes a resample a sum over rows
    rather than a rebuild of the cohort.
    """
    import numpy as np
    import pandas as pd

    called = frame[called_column].to_numpy(dtype=bool)
    resistant = frame[resistant_column].to_numpy(dtype=bool)
    codes, names = pd.factorize(frame.CLUSTER)
    size = len(names)
    columns = [called & resistant, called & ~resistant,
               ~called & resistant, ~called & ~resistant]
    return np.column_stack([
        np.bincount(codes, weights=weights, minlength=size) for weights in columns])


def cluster_interval(frame, called_column, resistant_column, rng, draws=BOOTSTRAPS):
    """Percentile intervals for each metric, resampling clusters."""
    import numpy as np

    cells = cluster_cells(frame, called_column, resistant_column)
    size = len(cells)
    collected = {key: [] for key in ("sensitivity", "specificity", "ppv", "npv")}
    for _ in range(draws):
        totals = cells[rng.integers(0, size, size)].sum(axis=0)
        values = metrics(dict(zip(
            ("true_positive", "false_positive", "false_negative", "true_negative"),
            totals)))
        for key, value in values.items():
            if value is not None:
                collected[key].append(value)
    intervals = {}
    for key, values in collected.items():
        if len(values) < draws // 4:
            intervals[key] = None
        else:
            intervals[key] = (float(np.percentile(values, 2.5)),
                              float(np.percentile(values, 97.5)))
    return intervals


def evaluate(frame, flags, drug):
    """Every rule against one drug at the ECOFF, on the samples carrying a MIC
    for it."""
    subset = frame[frame[f"MIC_{drug}"].notna()].copy()
    for rule in RULES:
        subset[rule] = flags[rule].reindex(subset.index).fillna(False)
    resistant = f"resistant_{drug}"

    records = []
    for rule in RULES:
        counts = confusion(subset[rule], subset[resistant])
        values = metrics(counts)
        intervals = cluster_interval(
            subset, rule, resistant, stream(f"{drug} {rule} {ECOFF[drug]:g}"))
        record = {"drug": drug, "rule": rule, "isolates": len(subset),
                  "clusters": subset.CLUSTER.nunique(),
                  "resistant": int(subset[resistant].sum()), **counts}
        for key, value in values.items():
            record[key] = round(value, 4) if value is not None else None
            interval = intervals[key]
            record[f"{key}_low"] = round(interval[0], 4) if interval else None
            record[f"{key}_high"] = round(interval[1], 4) if interval else None
        records.append(record)
    return records


def shared_thresholds(series, drug):
    """The tested concentrations present on every plate design for this drug.

    A cut-off is usable only where every isolate's position relative to it is
    determinate. Below the highest of the designs' lowest rungs, a left-censored
    reading on the design with the lower floor sits on neither side of the
    cut-off; above the lowest of the designs' highest rungs, a right-censored
    reading sits on neither side either. A concentration that is a rung on every
    design is inside both bounds by construction, so the shared rungs are exactly
    the usable cut-offs.

    Concentrations are matched within 5%, as mic_bounds matches them, because
    CRyPTIC label them as rounded values.
    """
    ladders = [concentrations for (_, code), concentrations in series.items()
               if code == drug]
    if len(ladders) < 2:
        return list(ladders[0]) if ladders else []
    shared = []
    for candidate in ladders[0]:
        if all(any(abs(rung - candidate) / candidate < 0.05 for rung in ladder)
               for ladder in ladders[1:]):
            shared.append(candidate)
    return shared


def above_threshold(lower, upper, cut):
    """Whether each interval lies above a cut-off, and whether that is
    determinate.

    The MIC is known only as an interval. It lies above the cut-off when its
    lower bound reaches it, and at or below when its upper bound does not exceed
    it. Anything else is a reading the plate cannot place either side.
    """
    import numpy as np

    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    above = lower >= cut
    at_or_below = upper <= cut
    return above, above | at_or_below


def sweep(frame, flags, drug, thresholds, draws=BOOTSTRAPS):
    """Every rule against every cut-off, rather than against the ECOFF alone."""
    import numpy as np

    subset = frame[frame[f"MIC_{drug}"].notna()].copy()
    for rule in RULES:
        subset[rule] = flags[rule].reindex(subset.index).fillna(False)

    records = []
    for concentration in thresholds:
        cut = float(np.log2(concentration))
        above, determinate = above_threshold(
            subset[f"lower_{drug}"], subset[f"upper_{drug}"], cut)
        usable = subset[determinate].copy()
        usable["above"] = above[determinate]
        for rule in RULES:
            counts = confusion(usable[rule], usable["above"])
            values = metrics(counts)
            intervals = cluster_interval(
                usable, rule, "above",
                stream(f"{drug} {rule} {concentration:g}"), draws=draws)
            record = {"drug": drug, "rule": rule,
                      "threshold_mg_L": concentration,
                      "threshold_log2": round(cut, 3),
                      "is_ecoff": bool(
                          abs(concentration - ECOFF[drug]) / ECOFF[drug] < 0.05),
                      "isolates": len(usable),
                      "indeterminate": int((~determinate).sum()),
                      "clusters": usable.CLUSTER.nunique(),
                      "above_threshold": int(usable["above"].sum()), **counts}
            for key, value in values.items():
                record[key] = round(value, 4) if value is not None else None
                interval = intervals[key]
                record[f"{key}_low"] = round(interval[0], 4) if interval else None
                record[f"{key}_high"] = round(interval[1], 4) if interval else None
            records.append(record)
    return records


def main():
    import numpy as np
    import pandas as pd

    say("=" * 72)
    say("Genotype rules as tests for resistance at the ECOFF")
    say("=" * 72)

    mutations = cohort.load_mutations()
    status = cohort.build_status(mutations)
    frame = add_clusters(cohort.assemble(status), mutations)
    flags = rule_flags(mutations, frame.index)

    series, _, _ = concentration_series()
    for drug in cohort.DRUGS:
        lower, upper, _, off_series = place_mics(
            frame[f"MIC_{drug}"], frame[f"PLATEDESIGN_{drug}"], series, drug)
        if off_series:
            raise ValueError(f"{drug} MICs off the tested series: {off_series[:5]}")
        frame[f"lower_{drug}"] = lower
        frame[f"upper_{drug}"] = upper

    say("\nIsolates each rule calls, across the whole cohort:")
    for rule in RULES:
        say(f"  {rule:16s} {int(flags[rule].sum()):,} of {len(frame):,}")

    records = []
    for drug in cohort.DRUGS:
        records.extend(evaluate(frame, flags, drug))
    table = pd.DataFrame(records)

    for drug in cohort.DRUGS:
        block = table[table.drug.eq(drug)]
        say(f"\n{'-' * 72}")
        say(f"{drug}: {int(block.isolates.iloc[0]):,} isolates in "
            f"{int(block.clusters.iloc[0]):,} clusters, "
            f"{int(block.resistant.iloc[0]):,} resistant")
        say("-" * 72)
        say(f"\n  {'rule':16s} {'sens':>18s} {'spec':>18s} {'PPV':>18s}")
        for row in block.itertuples():
            cells = []
            for key in ("sensitivity", "specificity", "ppv"):
                value = getattr(row, key)
                low = getattr(row, f"{key}_low")
                high = getattr(row, f"{key}_high")
                cells.append(f"{100 * value:5.1f} ({100 * low:4.1f} to {100 * high:4.1f})"
                             if value is not None and low is not None else "n/a")
            say(f"  {row.rule:16s} " + " ".join(f"{cell:>18s}" for cell in cells))
        say("\n  cells, as true positive, false positive, false negative, true negative:")
        for row in block.itertuples():
            say(f"  {row.rule:16s} {row.true_positive:5d} {row.false_positive:6d} "
                f"{row.false_negative:5d} {row.true_negative:7d}")

    say("\nThe predictive values describe this collection, whose resistance is")
    say("concentrated at one site by design, and transfer to no other.")

    say("\n" + "=" * 72)
    say("The same rules against every usable cut-off")
    say("=" * 72)
    say("")
    say("A cut-off is usable where every isolate's position relative to it is")
    say("determinate, which is at the concentrations both plate designs tested.")
    say("The ECOFF row of each table reproduces the figures above.")

    sweeps = []
    for drug in cohort.DRUGS:
        thresholds = shared_thresholds(series, drug)
        say(f"\n{drug}: {len(thresholds)} usable cut-offs, "
            f"{thresholds[0]:g} to {thresholds[-1]:g} mg/L")
        sweeps.extend(sweep(frame, flags, drug, thresholds))
    sweep_table = pd.DataFrame(sweeps)

    for drug in cohort.DRUGS:
        for rule in RULES:
            block = sweep_table[sweep_table.drug.eq(drug) & sweep_table.rule.eq(rule)]
            say(f"\n  {drug}, {rule}")
            say(f"    {'mg/L':>7s} {'above':>7s} {'sens':>18s} {'spec':>18s} "
                f"{'PPV':>18s}")
            for row in block.itertuples():
                cells = []
                for key in ("sensitivity", "specificity", "ppv"):
                    value = getattr(row, key)
                    low = getattr(row, f"{key}_low")
                    high = getattr(row, f"{key}_high")
                    cells.append(
                        f"{100 * value:5.1f} ({100 * low:4.1f} to {100 * high:4.1f})"
                        if value is not None and low is not None else "n/a")
                marker = " <- ECOFF" if row.is_ecoff else ""
                say(f"    {row.threshold_mg_L:7g} {row.above_threshold:7d} "
                    + " ".join(f"{cell:>18s}" for cell in cells) + marker)

    say("")
    say("Sensitivity rises and the positive predictive value falls as the cut-off")
    say("drops, because a lower cut-off counts more isolates as resistant and the")
    say("rules call a fixed set. Where a rule's specificity stays high across the")
    say("range, what it calls is not following the cut-off at all.")

    TABLE.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(TABLE, index=False)
    sweep_table.to_csv(THRESHOLDS, index=False)
    REPORT.write_text("\n".join(_lines) + "\n")
    print(f"\nTable written to {TABLE}")
    print(f"Cut-off sweep written to {THRESHOLDS} ({len(sweep_table)} rows)")
    print(f"Report written to {REPORT}")


if __name__ == "__main__":
    main()
