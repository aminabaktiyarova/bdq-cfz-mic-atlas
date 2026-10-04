"""
Splits the cohort by variant class and tests the key comparisons with site and
lineage held constant.

Run from the project root with the virtual environment active:

    python code/analyse_groups.py

Three things motivate this module.

Loss of function was pooling gene deletions, frameshifts and stop codons, and
those did not look alike in the first pass. They are separated here.

The samples excluded as uncertain were carrying far more resistance than the
reference group, so they are profiled rather than discarded silently.

Every earlier comparison was unadjusted. Resistance in this dataset is heavily
concentrated at one contributing site by design, and lineage is unevenly
distributed across sites, so an unadjusted odds ratio confounds genotype with
where the sample came from. The Mantel-Haenszel estimate here holds the
stratifying variable constant. The accompanying homogeneity test asks whether
one pooled estimate is defensible at all: if the odds ratio differs across
strata, the pooled figure is a summary of things that should not be summarised.

Output: outputs/group_analysis_report.txt
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import cohort  # noqa: E402

REPORT = Path("outputs/group_analysis_report.txt")
MIN_STRATUM = 20  # a stratum smaller than this carries no usable information

# The lineages carrying enough Rv0678 loss-of-function isolates for a stratum to
# say anything. lineage1 carries two, which is why a homogeneity test over every
# lineage answers a different question from one over these three.
PRINCIPAL_LINEAGES = ("lineage2", "lineage3", "lineage4")

_lines = []


def say(text=""):
    print(text)
    _lines.append(text)


def write_report():
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(_lines) + "\n")


def fisher(a, n1, b, n2):
    """Two-sided Fisher exact test on counts a of n1 against b of n2."""
    from scipy.stats import fisher_exact

    odds, p = fisher_exact([[a, n1 - a], [b, n2 - b]])
    return odds, p


def group_rows(df, drug):
    """Every group except the reference group, against the reference group.

    Each group in the cohort appears, including the ones that are not a single
    variant class: a sample carrying several variants and a sample whose calls
    are uncertain both have a resistance rate, and a table that leaves them out
    cannot be read against the cohort.
    """
    reference = df[df.GROUP.eq("reference")]
    b, n2 = int(reference[f"resistant_{drug}"].sum()), len(reference)
    rows = []
    for group in df.GROUP.dropna().unique():
        if group == "reference":
            continue
        sub = df[df.GROUP.eq(group)]
        a, n1 = int(sub[f"resistant_{drug}"].sum()), len(sub)
        if not n1:
            continue
        odds, p = fisher(a, n1, b, n2)
        rows.append({
            "group": group, "n": n1, "resistant": a,
            "% R": round(100 * a / n1, 1),
            "odds ratio": round(odds, 1) if a else 0.0,
            "p": f"{p:.2g}",
        })
    return sorted(rows, key=lambda row: -row["n"])


def direct(df, group_a, group_b, drug):
    """group_a against group_b, rather than each against the reference group."""
    first, second = df[df.GROUP.eq(group_a)], df[df.GROUP.eq(group_b)]
    a, n1 = int(first[f"resistant_{drug}"].sum()), len(first)
    b, n2 = int(second[f"resistant_{drug}"].sum()), len(second)
    odds, p = fisher(a, n1, b, n2)
    return {"drug": drug, "group": group_a, "n": n1, "resistant": a,
            "against": group_b, "against n": n2, "against resistant": b,
            "odds ratio": round(odds, 2), "p": f"{p:.3g}"}


def stratified(df, group_a, group_b, drug, by, only=None):
    """
    Mantel-Haenszel odds ratio for group_a against group_b, holding `by` constant.

    Strata contributing no information (too few samples, or no exposed samples)
    are dropped and reported, because including them adds noise without adding
    evidence. `only` restricts the estimate to the strata it names, which is a
    different question from an estimate over every stratum and is reported as
    such.
    """
    import numpy as np
    import pandas as pd
    from statsmodels.stats.contingency_tables import StratifiedTable

    subset = df[df.GROUP.isin([group_a, group_b])]
    if only is not None:
        subset = subset[subset[by].isin(only)]
        say(f"  strata restricted to {', '.join(map(str, only))}")
    rows, tables, dropped = [], [], []

    for stratum, chunk in subset.groupby(by, observed=True):
        exposed = chunk[chunk.GROUP == group_a]
        unexposed = chunk[chunk.GROUP == group_b]
        a, n1 = int(exposed[f"resistant_{drug}"].sum()), len(exposed)
        b, n2 = int(unexposed[f"resistant_{drug}"].sum()), len(unexposed)
        rows.append({
            "stratum": stratum, "exposed_n": n1, "exposed_R": a,
            "reference_n": n2, "reference_R": b,
        })
        # A stratum with no exposed samples, too few reference samples, or no
        # resistant samples at all carries no information about an odds ratio.
        # Leaving the last kind in leaves the homogeneity test undefined.
        if n1 == 0 or n2 < MIN_STRATUM or (a + b) == 0:
            dropped.append(stratum)
            continue
        tables.append(np.array([[a, n1 - a], [b, n2 - b]]))

    table = pd.DataFrame(rows).sort_values("exposed_n", ascending=False)
    say(table.to_string(index=False))

    if dropped:
        say(f"\n  strata carrying no information, dropped: {', '.join(map(str, dropped))}")
    if len(tables) < 2:
        say("  fewer than two informative strata; no adjusted estimate possible")
        return

    def fmt(p):
        """Very small p-values underflow to zero; report them as a bound."""
        return "< 1e-300" if p == 0 else f"{p:.3g}"

    st = StratifiedTable(tables)
    say(f"\n  Mantel-Haenszel odds ratio, {by} held constant: {st.oddsratio_pooled:.2f}")
    low, high = st.oddsratio_pooled_confint()
    say(f"  95% confidence interval: {low:.2f} to {high:.2f}")
    say(f"  test of no association: p = {fmt(st.test_null_odds().pvalue)}")
    homogeneity = st.test_equal_odds()
    say(f"  test of equal odds across strata: p = {fmt(homogeneity.pvalue)}")
    if homogeneity.pvalue < 0.05:
        say("  the odds ratio differs across strata, so the pooled figure should not")
        say("  be reported on its own; the per-stratum table above is the result")


def main():
    import pandas as pd

    say("=" * 72)
    say("Group analysis: variant class, site and lineage")
    say("=" * 72)

    df = cohort.assemble()
    say(f"\nSamples with genotype status and a UKMYC MIC: {len(df):,}")

    # ------------------------------------------------------- 1. refined groups
    say("\n1. Refined groups, with loss of function split by class")
    for drug in cohort.DRUGS:
        say(f"\n{drug}:")
        rows = []
        for group in cohort.GROUP_ORDER:
            sub = df[df.GROUP == group]
            if not len(sub):
                continue
            resistant = int(sub[f"resistant_{drug}"].sum())
            rows.append({
                "group": group,
                "n": len(sub),
                "resistant": resistant,
                "% R": round(100 * resistant / len(sub), 1),
                "median log2MIC": round(sub[f"LOG2MIC_{drug}"].median(), 2),
                "left-cens": int(sub[f"censored_left_{drug}"].sum()),
                "right-cens": int(sub[f"censored_right_{drug}"].sum()),
            })
        say(pd.DataFrame(rows).to_string(index=False))

    # ---------------------------------------------- 2. group against reference
    say("\n2. Each group against the reference group, unadjusted")
    reference = df[df.GROUP == "reference"]
    for drug in cohort.DRUGS:
        say(f"\n{drug}:  reference is {int(reference[f'resistant_{drug}'].sum())}"
            f" of {len(reference):,} resistant")
        say(pd.DataFrame(group_rows(df, drug)).to_string(index=False))

    # ------------------------------------------------- 3. mmpL5 within classes
    say("\n3. The mmpL5 covariate, within each Rv0678 class")
    for drug in cohort.DRUGS:
        say(f"\n{drug}:")
        rows = []
        for group in cohort.GROUP_ORDER:
            if not group.startswith("Rv0678"):
                continue
            sub = df[df.GROUP == group]
            if len(sub) < 5:
                continue
            for disrupted in (False, True):
                part = sub[sub.mmpL5_LOF == disrupted]
                if not len(part):
                    continue
                rows.append({
                    "group": group,
                    "mmpL5 disrupted": disrupted,
                    "n": len(part),
                    "resistant": int(part[f"resistant_{drug}"].sum()),
                    "% R": round(100 * part[f"resistant_{drug}"].mean(), 1),
                    "median log2MIC": round(part[f"LOG2MIC_{drug}"].median(), 2),
                })
        say(pd.DataFrame(rows).to_string(index=False))

    say("\n  mmpL5 disruption in the reference group, as a control:")
    rows = []
    for disrupted in (False, True):
        part = reference[reference.mmpL5_LOF == disrupted]
        rows.append({
            "mmpL5 disrupted": disrupted, "n": len(part),
            "BDQ % R": round(100 * part.resistant_BDQ.mean(), 2),
            "CFZ % R": round(100 * part.resistant_CFZ.mean(), 2),
            "median log2 BDQ": round(part.LOG2MIC_BDQ.median(), 2),
            "median log2 CFZ": round(part.LOG2MIC_CFZ.median(), 2),
        })
    say(pd.DataFrame(rows).to_string(index=False))
    say("\n  (if losing mmpL5 lowers MICs even without an Rv0678 variant, the effect")
    say("   is about the pump itself, not only about the epistasis)")

    # ------------------------------------------------- 4. gene deletion detail
    say("\n4. Every Rv0678 gene deletion, with provenance")
    mutations = cohort.load_mutations(genes=["Rv0678"])
    deletions = mutations[mutations.KIND == "GENE_DELETION"][
        ["UNIQUEID", "DELETED_FRACTION", "IS_MINOR"]
    ].set_index("UNIQUEID")
    sites = cohort.load_sites()

    genomes = pd.read_parquet(cohort.DATA / "GENOMES.parquet").reset_index()
    meta = genomes.set_index("UNIQUEID")[["LINEAGE", "SUBLINEAGE", "TB_COVERAGE", "TB_DEPTH"]]
    detail = deletions.join(meta)
    detail["SITEID"] = [str(u).split(".")[1] if "." in str(u) else "?" for u in detail.index]
    if sites is not None:
        detail["COUNTRY"] = detail.SITEID.map(sites.COUNTRY)
    detail["has_MIC"] = detail.index.isin(df.index)
    for drug in cohort.DRUGS:
        detail[f"MIC_{drug}"] = df[f"MIC_{drug}"].reindex(detail.index)
    say(detail.sort_values(["SITEID", "DELETED_FRACTION"]).to_string())

    say("\n  Deletion calls by site and deleted fraction:")
    say(pd.crosstab(detail.SITEID, detail.DELETED_FRACTION).to_string())
    say("\n  Coverage and depth of the deleted samples against everything else:")
    others = meta.loc[~meta.index.isin(detail.index)]
    say(f"    deleted samples:  median coverage {detail.TB_COVERAGE.median():.2f},"
        f" median depth {detail.TB_DEPTH.median():.1f}")
    say(f"    all others:       median coverage {others.TB_COVERAGE.median():.2f},"
        f" median depth {others.TB_DEPTH.median():.1f}")
    say("\n  (repeated identical fractions at one site, in one lineage, would point to")
    say("   a clonal cluster or a calling artefact rather than independent events)")

    # ----------------------------------------------------- 5. uncertain group
    say("\n5. The uncertain group, profiled rather than discarded")
    uncertain = df[df.GROUP == "uncertain"]
    say(f"\n  samples: {len(uncertain):,}")
    for gene in cohort.BDQ_GENES:
        flagged = uncertain[uncertain[f"uncertain_{gene}"]]
        if not len(flagged):
            continue
        say(f"\n  uncertain at {gene}: {len(flagged):,}")
        for drug in cohort.DRUGS:
            resistant = int(flagged[f"resistant_{drug}"].sum())
            b, n2 = int(reference[f"resistant_{drug}"].sum()), len(reference)
            odds, p = fisher(resistant, len(flagged), b, n2)
            say(f"    {drug}: {resistant} of {len(flagged):,} resistant "
                f"({100*resistant/len(flagged):.1f}%), odds ratio vs reference "
                f"{odds:.1f}, p = {p:.2g}")

    say("\n  Samples carrying two or more variants:")
    multiple = df[df.GROUP == "multiple variants"]
    say(f"    n = {len(multiple):,}")
    for drug in cohort.DRUGS:
        if len(multiple):
            say(f"    {drug}: {int(multiple[f'resistant_{drug}'].sum())} resistant "
                f"({100*multiple[f'resistant_{drug}'].mean():.1f}%)")

    # ------------------------------------------------ 6. stratified estimates
    say("\n6. Key comparisons with site held constant")
    lof_groups = [f"Rv0678 {c}" for c in cohort.LOF_CLASSES]
    df = df.copy()
    df["GROUP_POOLED"] = df.GROUP.where(~df.GROUP.isin(lof_groups), "Rv0678 loss of function")

    intact = df[(~df.mmpL5_LOF) | (df.GROUP == "reference")].copy()
    intact["GROUP"] = intact.GROUP_POOLED

    for drug in cohort.DRUGS:
        say(f"\n{drug}: Rv0678 loss of function against reference, mmpL5 intact, by site")
        stratified(intact, "Rv0678 loss of function", "reference", drug, "SITEID")

    for drug in cohort.DRUGS:
        say(f"\n{drug}: Rv0678 substitution against reference, mmpL5 intact, by site")
        stratified(intact, "Rv0678 substitution", "reference", drug, "SITEID")

    say("\n7. The same comparisons with lineage held constant")
    for drug in cohort.DRUGS:
        say(f"\n{drug}: Rv0678 loss of function against reference, mmpL5 intact, by lineage")
        stratified(intact, "Rv0678 loss of function", "reference", drug, "LINEAGE")

    say("\n8. The lineage homogeneity test over the principal lineages")
    for drug in cohort.DRUGS:
        say(f"\n{drug}: Rv0678 loss of function against reference, mmpL5 intact, by lineage")
        stratified(intact, "Rv0678 loss of function", "reference", drug, "LINEAGE",
                   only=PRINCIPAL_LINEAGES)

    say("\n9. Loss of function against substitution, mmpL5 intact")
    say("")
    say(pd.DataFrame([direct(intact, "Rv0678 loss of function", "Rv0678 substitution",
                             drug) for drug in cohort.DRUGS]).to_string(index=False))

    write_report()
    say(f"\nWritten to {REPORT}")


if __name__ == "__main__":
    main()
