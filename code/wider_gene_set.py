"""
Whether a wider gene set accounts for the resistance the three genes do not.

Run from the project root with the virtual environment active:

    python code/wider_gene_set.py

code/unexplained.py reports that 80 of 171 bedaquiline-resistant and 527 of 683
clofazimine-resistant isolates carry no real major-allele variant in Rv0678,
pepQ or atpE, and records as a limit that the fraction describes that gene set
rather than the genome. The obvious objection is that the residue is an artifact
of stopping at three genes. This module tests it against ten further genes.

The gene set, Rv1979c, Rv2983, fbiA, fbiB, fbiC, fgd1, lpqB, mmpS5, mtrA and
mtrB, is the union of the genes TB-Profiler's database associates with
bedaquiline and with clofazimine, less the three already in scope and the
modifier gene. mmpR5 is excluded because it is that database's alias for
Rv0678 and the CRyPTIC tables carry no gene of that name. Only the gene names
are taken from it. No grade, confidence or other catalogue content is read
here, so every figure this module writes is licence-clean and belongs in
outputs/ beside the rest.

What is compared. The reference group code/cohort.py defines is the unexplained tier
for a resistant isolate, so the comparison is made inside it: among isolates
carrying no real major-allele variant in the three genes and nothing uncertain
in them, does carrying a variant in a candidate gene track resistance. An
isolate counts as a carrier when it holds a real major-allele variant in that
gene, by the same definition the rest of the project uses.

A carrier count on its own answers nothing. mtrB carries a variant in 91 per
cent of the cohort, so admitting it would empty the unexplained tier the way
admitting mmpL5 does in code/unexplained.py, and name no cause. The question is
association, not coverage.

Three things are reported against each crude odds ratio, because this
collection makes a crude figure easy to believe and hard to trust:

  Multiplicity. Twenty tests are run, ten genes on each of two drugs. Every
  p value carries a Benjamini-Hochberg q over that family of twenty, and the
  family is fixed before the tests rather than chosen after them.

  Stratification. Mantel-Haenszel estimates holding site constant and holding
  lineage constant, each with a homogeneity test. Where homogeneity rejects,
  the pooled figure summarises strata that disagree and is reported with that
  attached rather than read as one quantity.

  Clonality. Isolates sharing a site, a sublineage and a mutation are one event
  observed many times, so the interval resamples clusters keyed that way, and
  the report names how many site and sublineage groups the resistant carriers
  fall into, how large the biggest is, and what share the commonest single
  variant holds.

Outputs:
  outputs/wider_gene_set.csv
  outputs/wider_gene_set_report.txt
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import cohort  # noqa: E402
from prediction_metrics import cluster_cells, stream  # noqa: E402

REPORT = Path("outputs/wider_gene_set_report.txt")
TABLE = Path("outputs/wider_gene_set.csv")

# The union of what TB-Profiler's database lists for the two drugs, less the
# three genes in scope, the modifier gene, and the Rv0678 alias.
CANDIDATE_GENES = ["Rv1979c", "Rv2983", "fbiA", "fbiB", "fbiC", "fgd1",
                   "lpqB", "mmpS5", "mtrA", "mtrB"]

BOOTSTRAPS = 400

_lines = []


def say(text=""):
    print(text)
    _lines.append(text)


def carriers_by_gene(mutations):
    """The isolates carrying a real major-allele variant in each gene."""
    major = mutations[mutations.REAL_MAJOR]
    return {gene: set(block.UNIQUEID)
            for gene, block in major.groupby("GENE", observed=True)}


def gene_clusters(frame, mutations, gene):
    """A cluster key for one candidate gene.

    Carriers cluster on site, sublineage and the variant they carry, which is
    the same rule add_clusters applies to the three genes in scope. Every
    non-carrier is its own cluster: they share no mutation event.
    """
    import pandas as pd

    defining = (mutations[mutations.REAL_MAJOR & mutations.GENE.eq(gene)]
                .groupby("UNIQUEID", observed=True).MUTATION.first())
    carried = defining.reindex(frame.index)
    grouped = (frame.SITEID.astype(str) + " | " + frame.SUBLINEAGE.astype(str)
               + " | " + carried.astype(str))
    own = pd.Series(frame.index, index=frame.index).astype(str)
    return grouped.where(carried.notna(), own)


def odds_ratio(cells):
    """The odds ratio from the four cells, or None where a margin is empty."""
    a, b, c, d = cells
    if b == 0 or c == 0 or (a == 0 and d == 0):
        return None
    return float((a * d) / (b * c))


def resampled_interval(frame, carrier_column, resistant_column, rng,
                       draws=BOOTSTRAPS):
    """A percentile interval for the odds ratio, resampling clusters."""
    import numpy as np

    cells = cluster_cells(frame, carrier_column, resistant_column)
    size = len(cells)
    values = []
    for _ in range(draws):
        totals = cells[rng.integers(0, size, size)].sum(axis=0)
        value = odds_ratio(totals)
        if value is not None:
            values.append(value)
    if len(values) < draws // 4:
        return None
    return float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))


def stratified(frame, carrier_column, resistant_column, by):
    """Mantel-Haenszel across the strata of one column, with a homogeneity test.

    A stratum in which every isolate is a carrier, or none is, or every isolate
    is resistant, or none is, carries no information about the association and
    is dropped rather than contributing a zero margin.
    """
    import numpy as np
    from statsmodels.stats.contingency_tables import StratifiedTable

    tables = []
    for _, block in frame.groupby(by, observed=True):
        carrier = block[carrier_column].to_numpy(dtype=bool)
        resistant = block[resistant_column].to_numpy(dtype=bool)
        table = np.array(
            [[(carrier & resistant).sum(), (carrier & ~resistant).sum()],
             [(~carrier & resistant).sum(), (~carrier & ~resistant).sum()]],
            dtype=float)
        if table.sum(axis=0).min() > 0 and table.sum(axis=1).min() > 0:
            tables.append(table)
    if len(tables) < 2:
        return {"or": None, "low": None, "high": None,
                "strata": len(tables), "homogeneity_p": None}
    # The Mantel-Haenszel estimate is the sum of a*d/n over the sum of b*c/n.
    # Where no carrier anywhere is resistant the numerator is zero and the
    # estimate is zero, whose logarithm and variance are undefined, and the
    # homogeneity test then divides by zero. Both are reported as absent rather
    # than computed, which is what an estimate with no information is.
    stacked = np.dstack(tables)
    totals = stacked.sum(axis=2)
    numerator = sum(t[0, 0] * t[1, 1] / t.sum() for t in tables)
    denominator = sum(t[0, 1] * t[1, 0] / t.sum() for t in tables)
    if numerator == 0 or denominator == 0 or totals.min() == 0:
        return {"or": None, "low": None, "high": None,
                "strata": len(tables), "homogeneity_p": None}
    stratified_table = StratifiedTable(stacked)
    low, high = stratified_table.oddsratio_pooled_confint()
    return {"or": float(stratified_table.oddsratio_pooled),
            "low": float(low), "high": float(high), "strata": len(tables),
            "homogeneity_p": float(stratified_table.test_equal_odds().pvalue)}


def benjamini_hochberg(p_values):
    """q values over the whole family, in the order the p values were given."""
    import numpy as np

    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranks = np.arange(1, len(p) + 1)
    adjusted = p[order] * len(p) / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    # No clamp at 1 is needed. After the running minimum from the largest p
    # downwards, every value is at most the largest p times n over n, which is
    # the largest p itself, so no q can exceed it.
    q = np.empty_like(adjusted)
    q[order] = adjusted
    return q


def test_one(frame, mutations, carriers, gene, drug):
    """One gene against one drug, inside the reference group."""
    import math

    from scipy.stats import fisher_exact

    subset = frame[frame.IS_REFERENCE & frame[f"MIC_{drug}"].notna()].copy()
    subset["carrier"] = subset.index.isin(carriers.get(gene, set()))
    subset["CLUSTER"] = gene_clusters(subset, mutations, gene)
    resistant = f"resistant_{drug}"

    carrier = subset["carrier"].to_numpy(dtype=bool)
    is_resistant = subset[resistant].to_numpy(dtype=bool)
    cells = [int((carrier & is_resistant).sum()), int((carrier & ~is_resistant).sum()),
             int((~carrier & is_resistant).sum()), int((~carrier & ~is_resistant).sum())]
    table = [[cells[0], cells[2]], [cells[1], cells[3]]]
    estimate, p_value = fisher_exact(table)

    interval = resampled_interval(subset, "carrier", resistant,
                                  stream(f"wider {gene} {drug}"))
    site = stratified(subset, "carrier", resistant, "SITEID")
    lineage = stratified(subset, "carrier", resistant, "LINEAGE")

    resistant_carriers = subset.index[carrier & is_resistant]
    groups = (subset.loc[resistant_carriers, "SITEID"].astype(str) + " | "
              + subset.loc[resistant_carriers, "SUBLINEAGE"].astype(str))
    held = mutations[mutations.REAL_MAJOR & mutations.GENE.eq(gene)
                     & mutations.UNIQUEID.isin(set(resistant_carriers))]
    commonest = cohort.ranked_counts(held.MUTATION, 1)

    return {
        "drug": drug, "gene": gene,
        "isolates": len(subset), "resistant": int(is_resistant.sum()),
        "carriers": int(carrier.sum()),
        "carriers_resistant": cells[0], "carriers_susceptible": cells[1],
        "clusters": int(subset.CLUSTER.nunique()),
        "odds_ratio": round(float(estimate), 4) if math.isfinite(estimate) else None,
        "p_value": float(p_value),
        "or_low": round(interval[0], 4) if interval else None,
        "or_high": round(interval[1], 4) if interval else None,
        "site_or": round(site["or"], 4) if site["or"] is not None else None,
        "site_low": round(site["low"], 4) if site["low"] else None,
        "site_high": round(site["high"], 4) if site["high"] else None,
        "site_strata": site["strata"],
        "site_homogeneity_p": site["homogeneity_p"],
        "lineage_or": round(lineage["or"], 4) if lineage["or"] is not None else None,
        "lineage_low": round(lineage["low"], 4) if lineage["low"] else None,
        "lineage_high": round(lineage["high"], 4) if lineage["high"] else None,
        "lineage_strata": lineage["strata"],
        "lineage_homogeneity_p": lineage["homogeneity_p"],
        "resistant_carrier_groups": int(groups.nunique()),
        "largest_group": int(groups.value_counts().iloc[0]) if len(groups) else 0,
        "commonest_variant": commonest.index[0] if len(commonest) else None,
        "commonest_variant_count": int(commonest.iloc[0]) if len(commonest) else 0,
    }


def main():
    import pandas as pd

    say("=" * 76)
    say("Does a wider gene set account for the resistance the three genes do not")
    say("=" * 76)

    core = cohort.load_mutations()
    status = cohort.build_status(core)
    frame = cohort.assemble(status)
    wide = cohort.load_mutations(genes=CANDIDATE_GENES)
    carriers = carriers_by_gene(wide)

    say(f"\nCohort {len(frame):,} isolates, reference group "
        f"{int(frame.IS_REFERENCE.sum()):,}.")
    say("\nCarrier prevalence across the cohort:")
    say(f"  {'gene':9s} {'carriers':>9s} {'per cent':>9s}")
    for gene in CANDIDATE_GENES:
        n = len(carriers.get(gene, set()) & set(frame.index))
        say(f"  {gene:9s} {n:>9,} {100 * n / len(frame):>8.1f}%")

    records = [test_one(frame, wide, carriers, gene, drug)
               for drug in cohort.DRUGS for gene in CANDIDATE_GENES]
    table = pd.DataFrame(records)
    table["q_value"] = benjamini_hochberg(table.p_value)
    # The q belongs beside the p it adjusts rather than at the end of the row.
    order = list(table.columns)
    order.insert(order.index("p_value") + 1, order.pop(order.index("q_value")))
    table = table[order]

    say("\nTwenty tests, ten genes on each of two drugs. "
        "Benjamini-Hochberg q is over that family.")
    for drug in cohort.DRUGS:
        block = table[table.drug.eq(drug)].sort_values("p_value")
        first = block.iloc[0]
        say(f"\n{'-' * 76}")
        say(f"{drug}: {int(first.isolates):,} reference isolates, "
            f"{int(first.resistant):,} resistant")
        say("-" * 76)
        say(f"  {'gene':9s} {'R carr':>7s} {'S carr':>7s} {'OR':>7s} "
            f"{'95% resampling clusters':>26s} {'p':>9s} {'q':>9s}")
        for row in block.itertuples():
            interval = ("n/a" if pd.isna(row.or_low)
                        else f"{row.or_low:.2f} to {row.or_high:.2f}")
            estimate = ("n/a" if pd.isna(row.odds_ratio)
                        else f"{row.odds_ratio:.2f}")
            say(f"  {row.gene:9s} {row.carriers_resistant:>7,} "
                f"{row.carriers_susceptible:>7,} {estimate:>7s} {interval:>26s} "
                f"{row.p_value:>9.2g} {row.q_value:>9.2g}")

    say(f"\n{'=' * 76}")
    say("Every gene whose crude q falls below 0.05, held constant by site and by lineage")
    say("=" * 76)
    flagged = table[table.q_value < 0.05]
    if flagged.empty:
        say("\nNone.")
    for row in flagged.itertuples():
        say(f"\n{row.gene} against {row.drug}: crude OR "
            f"{row.odds_ratio:.2f}, p = {row.p_value:.2g}, q = {row.q_value:.2g}")
        for label, estimate, low, high, strata, homogeneity in (
                ("site", row.site_or, row.site_low, row.site_high,
                 row.site_strata, row.site_homogeneity_p),
                ("lineage", row.lineage_or, row.lineage_low, row.lineage_high,
                 row.lineage_strata, row.lineage_homogeneity_p)):
            if estimate is None:
                say(f"  {label:8s} too few usable strata")
                continue
            verdict = ("the strata disagree, so this pools incompatible things"
                       if homogeneity < 0.05 else "the strata agree")
            say(f"  {label:8s} MH OR {estimate:.2f} ({low:.2f} to {high:.2f}), "
                f"{strata} strata, homogeneity p = {homogeneity:.2g}, {verdict}")
        say(f"  {row.carriers_resistant} resistant carriers in "
            f"{row.resistant_carrier_groups} site and sublineage groups, "
            f"largest {row.largest_group}; commonest variant "
            f"{row.commonest_variant} in {row.commonest_variant_count}")

    survives = flagged[
        flagged.site_or.notna() & flagged.site_homogeneity_p.ge(0.05)
        & flagged.site_low.gt(1.0) & flagged.or_low.gt(1.0)]
    say(f"\n{'=' * 76}")
    say("Result")
    say("=" * 76)
    if survives.empty:
        say("\nNo gene in this set shows an association with unexplained resistance")
        say("that survives the multiplicity correction, an interval that resamples")
        say("clusters, and site stratification with its strata in agreement. Every")
        say("crude association above either fails one of those or rests on a few")
        say("clonal groups. The unexplained tier that code/unexplained.py counts")
        say("is not accounted for by these ten genes, which leaves that figure a")
        say("statement about a named gene set rather than about the genome.")
    else:
        for row in survives.itertuples():
            say(f"\n{row.gene} against {row.drug} survives: site-adjusted OR "
                f"{row.site_or:.2f} ({row.site_low:.2f} to {row.site_high:.2f}), "
                f"q = {row.q_value:.2g}, {row.resistant_carrier_groups} groups")

    TABLE.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(TABLE, index=False)
    REPORT.write_text("\n".join(_lines) + "\n")
    print(f"\nWritten: {TABLE}, {REPORT}")


if __name__ == "__main__":
    main()
