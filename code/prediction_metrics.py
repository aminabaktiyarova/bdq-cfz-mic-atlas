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
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import cohort  # noqa: E402
from cluster_adjust import add_clusters  # noqa: E402

REPORT = Path("outputs/prediction_metrics.txt")
TABLE = Path("outputs/prediction_metrics.csv")

BOOTSTRAPS = 400
SEED = 20260101
RULES = ["major variant", "Rv0678 any", "Rv0678 loss", "major or minor"]

_lines = []


def say(text=""):
    print(text)
    _lines.append(text)


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


def evaluate(frame, flags, drug, rng):
    """Every rule against one drug, on the samples carrying a MIC for it."""
    subset = frame[frame[f"MIC_{drug}"].notna()].copy()
    for rule in RULES:
        subset[rule] = flags[rule].reindex(subset.index).fillna(False)
    resistant = f"resistant_{drug}"

    records = []
    for rule in RULES:
        counts = confusion(subset[rule], subset[resistant])
        values = metrics(counts)
        intervals = cluster_interval(subset, rule, resistant, rng)
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

    say("\nIsolates each rule calls, across the whole cohort:")
    for rule in RULES:
        say(f"  {rule:16s} {int(flags[rule].sum()):,} of {len(frame):,}")

    rng = np.random.default_rng(SEED)
    records = []
    for drug in cohort.DRUGS:
        records.extend(evaluate(frame, flags, drug, rng))
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

    TABLE.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(TABLE, index=False)
    REPORT.write_text("\n".join(_lines) + "\n")
    print(f"\nTable written to {TABLE}")
    print(f"Report written to {REPORT}")


if __name__ == "__main__":
    main()
