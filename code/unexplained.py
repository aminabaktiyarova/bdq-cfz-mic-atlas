"""
Genotypically unexplained bedaquiline and clofazimine resistance.

Run from the project root with the virtual environment active:

    python code/unexplained.py

A resistant isolate is counted as genotypically unexplained when it carries no
real major-allele variant in Rv0678, pepQ or atpE, which is the reference group
defined in code/cohort.py. Every classification here comes from mutation
strings. No catalogue grading is read, and neither the EFFECTS table nor the
PREDICTIONS table is opened, so these counts contain only CC BY 4.0
CRyPTIC-derived content and can be released beside the rest of the derived data.

Four tiers partition the resistant isolates.

  attributable. The isolate carries at least one variant in a class whose MIC
  shift was estimated on this cohort with a bootstrap interval above zero:
  Rv0678 loss of function, Rv0678 substitution, and pepQ.

  carrier, no demonstrated effect. The isolate carries a variant only in classes
  with no estimated upward shift here: Rv0678 promoter variants, Rv0678 in-frame
  indels, and atpE. Rv0678 promoter variants shift the bedaquiline MIC downward
  by 1.63 doublings with an interval excluding zero, and atpE and the in-frame
  indels hold too few isolates for an estimate. The 1.63 figure is the
  discovery-half estimate produced by code/discovery.py.

  unexplained. No real major-allele variant in any of the three genes.

  indeterminate. One of the three genes carries a null call, a het call or a
  minor allele, so the isolate can be asserted neither to carry a variant nor to
  lack one.

The boundary between the first two tiers uses effect estimates computed on this
same cohort, so that split is descriptive and is not independent evidence. The
unexplained tier is structural and uses no effect estimate.

mmpL5 enters no tier. It carries a real major-allele variant in 53,361 of 54,057
genomes, so admitting it would make almost every isolate a carrier. Its role is
a covariate, as recorded in code/cohort.py.

These are not prevalence figures. Bedaquiline and clofazimine resistance in this
collection is concentrated at one site by design, so the resistant denominator
is not a population sample. The Wilson intervals describe sampling error within
this collection alone.

Site and sublineage are reported beside the unexplained tier because its
isolates share no defining mutation and cannot be clustered on one. The count of
distinct site and sublineage combinations gives a coarse floor on the number of
independent events.

The unexplained count depends on which genes are admitted, so the report gives
it under three gene sets: Rv0678 alone, the three genes used throughout, and
those three with mmpL5 added. The middle set is the definition used everywhere
else here, and the report checks that it returns the same counts as the tiers.

Outputs:
  outputs/unexplained_report.txt
  outputs/unexplained_counts.csv
  outputs/unexplained_gene_sets.csv
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import cohort  # noqa: E402

REPORT = Path("outputs/unexplained_report.txt")
TABLE = Path("outputs/unexplained_counts.csv")
GENE_SET_TABLE = Path("outputs/unexplained_gene_sets.csv")

# Variant classes with an upward MIC shift estimated on this cohort, bootstrap
# interval above zero.
EFFECT_LABELS = frozenset({
    "Rv0678 frameshift",
    "Rv0678 stop codon",
    "Rv0678 gene deletion",
    "Rv0678 substitution",
    "pepQ",
})

# Variant classes carrying no estimated upward shift on this cohort.
NO_EFFECT_LABELS = frozenset({
    "Rv0678 promoter",
    "Rv0678 in-frame indel",
    "atpE",
})

TIER_ORDER = [
    "attributable",
    "carrier, no demonstrated effect",
    "unexplained",
    "indeterminate",
]

# Gene sets the unexplained count is reported under. The second is the
# definition used everywhere else in the project.
GENE_SETS = [
    ("Rv0678 alone", ["Rv0678"]),
    ("Rv0678, pepQ, atpE", cohort.BDQ_GENES),
    ("those three and mmpL5", cohort.ALL_GENES),
]

_lines = []


def say(text=""):
    print(text)
    _lines.append(text)


def carried_labels(mutations):
    """Per sample, the set of variant-class labels it carries in the three genes.

    Rv0678 is labelled by class because its classes do not behave alike. pepQ
    and atpE are labelled by gene, matching the grouping in code/cohort.py.
    """
    import numpy as np
    import pandas as pd

    real = mutations[mutations.REAL_MAJOR & mutations.GENE.isin(cohort.BDQ_GENES)]
    gene = real.GENE.astype(str)
    label = np.where(gene.eq("Rv0678"), "Rv0678 " + real.CLASS.astype(str), gene)
    held = pd.DataFrame({"UNIQUEID": real.UNIQUEID.astype(str), "LABEL": label})
    return held.groupby("UNIQUEID").LABEL.agg(frozenset)


def assign_tiers(frame, labels):
    """Tier every sample in the cohort frame. Raises if a label set and a group
    label disagree about whether the sample carries a variant."""
    import pandas as pd

    held = labels.reindex(frame.index.astype(str))
    tiers = []
    for group, carried in zip(frame.GROUP, held):
        carried = frozenset() if not isinstance(carried, frozenset) else carried
        if group == "reference" and carried:
            raise ValueError("a reference sample carries a variant label")
        if group not in ("reference", "uncertain") and not carried:
            raise ValueError(f"group {group!r} carries no variant label")
        if group == "uncertain":
            tiers.append("indeterminate")
        elif not carried:
            tiers.append("unexplained")
        elif carried & EFFECT_LABELS:
            tiers.append("attributable")
        else:
            tiers.append("carrier, no demonstrated effect")
    return pd.Series(tiers, index=frame.index, name="TIER")


def tier_counts(frame, drug, high_quality=False):
    """Tier counts among the resistant isolates for one drug, with Wilson
    intervals on each fraction."""
    import pandas as pd
    from statsmodels.stats.proportion import proportion_confint

    subset = frame
    if high_quality:
        subset = subset[subset[f"PHENOTYPE_QUALITY_{drug}"].eq("HIGH")]
    resistant = subset[subset[f"resistant_{drug}"]]
    total = len(resistant)

    rows = []
    for tier in TIER_ORDER:
        count = int((resistant.TIER == tier).sum())
        if total:
            low, high = proportion_confint(count, total, method="wilson")
        else:
            low = high = float("nan")
        rows.append({
            "drug": drug,
            "quality": "HIGH" if high_quality else "all",
            "tier": tier,
            "isolates": count,
            "resistant_total": total,
            "percent": round(100 * count / total, 1) if total else float("nan"),
            "wilson_low": round(100 * low, 1) if total else float("nan"),
            "wilson_high": round(100 * high, 1) if total else float("nan"),
        })
    return pd.DataFrame(rows)


def coarse_events(frame, drug, tier="unexplained"):
    """Site and sublineage combinations spanned by one tier's resistant isolates.

    value_counts leaves the order of equal counts to the implementation, and
    that order differs between pandas versions, so ties are ordered by key and
    the report is reproducible across environments.
    """
    subset = frame[frame[f"resistant_{drug}"] & frame.TIER.eq(tier)]
    combinations = subset.SITEID.astype(str) + " | " + subset.SUBLINEAGE.astype(str)
    counts = combinations.value_counts().sort_index().sort_values(
        ascending=False, kind="stable")
    return len(subset), counts


def structural_counts(frame, mutations, genes, drug, high_quality=False):
    """Unexplained, indeterminate and carrier counts among the resistant
    isolates when only the named genes are admitted.

    This reads the mutation frame directly instead of the group labels, so it
    can be run under a gene set other than the one the cohort is defined on. An
    uncallable gene takes precedence over a carried variant, as it does in
    assign_tiers.
    """
    subset = frame
    if high_quality:
        subset = subset[subset[f"PHENOTYPE_QUALITY_{drug}"].eq("HIGH")]
    resistant = subset.index[subset[f"resistant_{drug}"]].astype(str)

    rows = mutations[mutations.GENE.isin(genes)]
    carriers = set(rows[rows.REAL_MAJOR].UNIQUEID.astype(str))
    uncallable = set(rows[rows.UNCERTAIN].UNIQUEID.astype(str))

    indeterminate = sum(1 for sample in resistant if sample in uncallable)
    unexplained = sum(1 for sample in resistant
                      if sample not in uncallable and sample not in carriers)
    return {
        "drug": drug,
        "quality": "HIGH" if high_quality else "all",
        "genes": " ".join(genes),
        "resistant_total": len(resistant),
        "unexplained": unexplained,
        "percent": round(100 * unexplained / len(resistant), 1) if len(resistant) else float("nan"),
        "indeterminate": indeterminate,
        "carrier": len(resistant) - indeterminate - unexplained,
    }


def gene_set_table(frame, mutations):
    """The unexplained count under each gene set, for both drugs."""
    import pandas as pd

    rows = []
    for label, genes in GENE_SETS:
        for drug in cohort.DRUGS:
            row = structural_counts(frame, mutations, genes, drug)
            row["gene_set"] = label
            rows.append(row)
    # gene_set names the set that genes lists, so the two sit together.
    return pd.DataFrame(rows)[["drug", "quality", "genes", "gene_set",
                               "resistant_total", "unexplained", "percent",
                               "indeterminate", "carrier"]]


def report_drug(frame, drug):
    import pandas as pd

    tables = [tier_counts(frame, drug), tier_counts(frame, drug, high_quality=True)]
    for table in tables:
        quality = table.quality.iloc[0]
        total = int(table.resistant_total.iloc[0])
        say(f"\n  {drug}, {quality} phenotype quality, {total:,} resistant isolates")
        say(f"    {'tier':32s} {'n':>6s} {'%':>7s}   95% Wilson")
        for row in table.itertuples():
            say(f"    {row.tier:32s} {row.isolates:6,d} {row.percent:6.1f}%"
                f"   {row.wilson_low:5.1f} to {row.wilson_high:5.1f}")

    isolates, counts = coarse_events(frame, drug)
    say(f"\n  Unexplained resistant isolates: {isolates:,} across {len(counts):,} "
        f"site and sublineage combinations")
    if len(counts):
        say(f"    largest combination: {counts.index[0]}, {counts.iloc[0]:,} isolates")
        say("    combinations holding four or more isolates:")
        for key, value in counts[counts >= 4].items():
            say(f"      {key:44s} {value:,}")
    return pd.concat(tables, ignore_index=True)


def main():
    import pandas as pd

    say("=" * 72)
    say("Genotypically unexplained bedaquiline and clofazimine resistance")
    say("=" * 72)
    say("\nClassification from mutation strings only. No catalogue grading is read.")

    mutations = cohort.load_mutations()
    status = cohort.build_status(mutations)
    frame = cohort.assemble(status)
    frame["TIER"] = assign_tiers(frame, carried_labels(mutations))

    say(f"\nCohort: {len(frame):,} samples with a genome and a UKMYC MIC")
    say("\nTier sizes across the whole cohort:")
    sizes = frame.TIER.value_counts()
    for tier in TIER_ORDER:
        say(f"  {tier:32s} {int(sizes.get(tier, 0)):,}")

    tables = []
    for drug in cohort.DRUGS:
        say(f"\n{'-' * 72}")
        say(f"{drug}")
        say("-" * 72)
        tables.append(report_drug(frame, drug))

    say(f"\n{'-' * 72}")
    say("Sensitivity of the unexplained count to the gene set")
    say("-" * 72)
    sets = gene_set_table(frame, mutations)
    say(f"\n  {'gene set':24s} {'drug':5s} {'resistant':>10s} "
        f"{'unexplained':>12s} {'%':>7s} {'indeterminate':>14s}")
    for row in sets.itertuples():
        say(f"  {row.gene_set:24s} {row.drug:5s} {row.resistant_total:10,d} "
            f"{row.unexplained:12,d} {row.percent:6.1f}% {row.indeterminate:14,d}")

    # The gene set the cohort is defined on must return the tier counts.
    table = pd.concat(tables, ignore_index=True)
    for drug in cohort.DRUGS:
        tier_row = table[table.drug.eq(drug) & table.quality.eq("all")]
        check = sets[sets.drug.eq(drug) & sets.gene_set.eq("Rv0678, pepQ, atpE")].iloc[0]
        for tier, column in (("unexplained", "unexplained"),
                             ("indeterminate", "indeterminate")):
            counted = int(tier_row[tier_row.tier.eq(tier)].isolates.iloc[0])
            if counted != int(check[column]):
                raise ValueError(
                    f"{drug} {tier}: tiers give {counted}, gene set gives {check[column]}")
    say("\n  The three-gene row reproduces the tier counts for both drugs.")

    TABLE.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(TABLE, index=False)
    sets.to_csv(GENE_SET_TABLE, index=False)
    REPORT.write_text("\n".join(_lines) + "\n")
    print(f"\nTables written to {TABLE} and {GENE_SET_TABLE}")
    print(f"Report written to {REPORT}")


if __name__ == "__main__":
    main()
