"""
Recomputes the key comparisons treating clonally related isolates as what they
are: repeated observations of one event rather than independent evidence.

Run from the project root with the virtual environment active:

    python code/cluster_adjust.py

Why this exists. Every count reported so far has treated each isolate as an
independent observation. Tuberculosis spreads clonally, so a single successful
strain can appear dozens of times in a collection, and forty-one isolates
sharing one mutation, one sublineage and one site are one event observed
forty-one times. Treating them as forty-one independent events inflates the
apparent evidence without adding any.

A cluster here is defined as the combination of site, sublineage, and the
genotype-defining mutation, with reference samples carrying "none" for the
mutation. This is a deliberately conservative proxy. A proper transmission
analysis would cluster on genomic distance, which needs the full variant data
and more compute than this project has. Two genuinely unrelated patients at one
site sharing a sublineage and a common mutation will be merged here, which
understates the evidence rather than overstating it. The direction of that error
is the acceptable one.

Three estimates are produced for each comparison, and they should be read
together rather than one being chosen:

  Isolate level. What has been reported so far. Correct only if isolates are
  independent, which they are not.

  Cluster-robust. A logistic model fitted to every isolate, with standard errors
  that allow observations within a cluster to be correlated. The point estimate
  uses all the data; the uncertainty acknowledges the clustering. This is the
  primary estimate.

  One isolate per cluster. Each cluster reduced to a single randomly chosen
  representative, repeated many times, reporting the spread of results across
  draws. This discards real information and is the most conservative view. If an
  effect survives here it is not an artefact of clonal expansion.

Output: outputs/cluster_report.txt
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import cohort  # noqa: E402

REPORT = Path("outputs/cluster_report.txt")
DRAWS = 200
SEED = 20260101
MIN_STRATUM = 20
MAIN_LINEAGES = ["lineage1", "lineage2", "lineage3", "lineage4"]

_lines = []


def say(text=""):
    print(text)
    _lines.append(text)


def write_report():
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(_lines) + "\n")


def add_clusters(df, mutations):
    """Attach a cluster key: site, sublineage, and the defining mutation."""
    import pandas as pd

    defining = mutations[mutations.REAL_MAJOR & mutations.GENE.isin(cohort.BDQ_GENES)]
    defining = defining.groupby("UNIQUEID", observed=True).MUTATION.first()

    df = df.copy()
    df["defining_mutation"] = defining.reindex(df.index).fillna("none")
    carrier = df.defining_mutation.ne("none")
    grouped = (
        df.SITEID.astype(str)
        + " | " + df.SUBLINEAGE.astype(str)
        + " | " + df.defining_mutation.astype(str)
    )
    # Carriers cluster on site, sublineage and mutation. Every non-carrier is
    # its own cluster: they share no mutation event, so there is nothing to
    # collapse, and collapsing them would shrink the reference group for no
    # reason.
    df["CLUSTER"] = grouped.where(carrier, pd.Series(df.index, index=df.index).astype(str))
    return df


def mh_odds_ratio(frame, drug, by="SITEID"):
    """Mantel-Haenszel odds ratio for EXPOSED against reference. None if unusable."""
    import numpy as np
    from statsmodels.stats.contingency_tables import StratifiedTable

    tables = []
    for _, chunk in frame.groupby(by, observed=True):
        exposed = chunk[chunk.EXPOSED]
        reference = chunk[~chunk.EXPOSED]
        a, n1 = int(exposed[f"resistant_{drug}"].sum()), len(exposed)
        b, n2 = int(reference[f"resistant_{drug}"].sum()), len(reference)
        if n1 == 0 or n2 < MIN_STRATUM or (a + b) == 0:
            continue
        tables.append(np.array([[a, n1 - a], [b, n2 - b]]))
    if len(tables) < 1:
        return None
    if len(tables) == 1:
        table = tables[0]
        if (table == 0).any():
            table = table + 0.5
        return float((table[0, 0] * table[1, 1]) / (table[0, 1] * table[1, 0]))
    return float(StratifiedTable(tables).oddsratio_pooled)


def robust_logit(frame, drug):
    """Logistic fit with standard errors clustered on CLUSTER."""
    import numpy as np
    import statsmodels.formula.api as smf

    data = frame.copy()
    data["y"] = data[f"resistant_{drug}"].astype(int)
    data["exposed"] = data.EXPOSED.astype(int)
    data["site"] = data.SITEID.astype(str)

    terms = ["exposed"]
    if data.site.nunique() > 1:
        counts = data.site.value_counts()
        data.loc[data.site.isin(counts[counts < 30].index), "site"] = "other"
        if data.site.nunique() > 1:
            terms.append("C(site)")
    try:
        fit = smf.logit("y ~ " + " + ".join(terms), data=data).fit(
            disp=False, maxiter=200, cov_type="cluster",
            cov_kwds={"groups": data.CLUSTER},
        )
    except Exception:  # noqa: BLE001
        return None
    if "exposed" not in fit.params.index:
        return None
    low, high = fit.conf_int().loc["exposed"]
    return {
        "odds_ratio": float(np.exp(fit.params["exposed"])),
        "ci_low": float(np.exp(low)),
        "ci_high": float(np.exp(high)),
        "p": float(fit.pvalues["exposed"]),
        "clusters": int(data.CLUSTER.nunique()),
    }


def collapsed_draws(frame, drug, rng):
    """Reduce each carrier cluster to one isolate, repeatedly, and re-estimate."""
    import numpy as np
    import pandas as pd

    carriers = frame[frame.EXPOSED]
    reference = frame[~frame.EXPOSED]
    estimates = []
    for _ in range(DRAWS):
        picked = carriers.groupby("CLUSTER", observed=True).sample(
            n=1, random_state=int(rng.integers(1 << 31))
        )
        estimate = mh_odds_ratio(pd.concat([picked, reference]), drug)
        if estimate is not None and np.isfinite(estimate):
            estimates.append(estimate)
    if len(estimates) < DRAWS // 4:
        return None
    values = np.array(estimates)
    return {
        "median": float(np.median(values)),
        "low": float(np.percentile(values, 2.5)),
        "high": float(np.percentile(values, 97.5)),
        "above_one": float((values > 1).mean()),
        "draws": len(values),
    }


def compare(frame, drug, label, rng):
    import pandas as pd

    exposed = frame[frame.EXPOSED]
    reference = frame[~frame.EXPOSED]
    if not len(exposed) or not len(reference):
        say(f"  {label}: no data")
        return

    from scipy.stats import fisher_exact

    a, n1 = int(exposed[f"resistant_{drug}"].sum()), len(exposed)
    b, n2 = int(reference[f"resistant_{drug}"].sum()), len(reference)
    crude, p_crude = fisher_exact([[a, n1 - a], [b, n2 - b]])

    say(f"\n  {label}")
    say(f"    isolates:  {n1} exposed ({a} resistant), {n2:,} reference ({b} resistant)")
    say(f"    clusters:  {exposed.CLUSTER.nunique()} exposed, {reference.CLUSTER.nunique():,} reference")
    say(f"    largest exposed cluster: {exposed.CLUSTER.value_counts().iloc[0]} isolates")
    say(f"    isolate-level odds ratio: {crude:.1f}  (p = {p_crude:.2g})")

    robust = robust_logit(frame, drug)
    if robust:
        say(f"    cluster-robust:           {robust['odds_ratio']:.1f}  "
            f"({robust['ci_low']:.1f} to {robust['ci_high']:.1f}), p = {robust['p']:.2g}")
    else:
        say("    cluster-robust:           did not fit")

    collapsed = collapsed_draws(frame, drug, rng)
    if collapsed:
        say(f"    one per cluster:          {collapsed['median']:.1f}  "
            f"({collapsed['low']:.1f} to {collapsed['high']:.1f} across {collapsed['draws']} draws), "
            f"above 1 in {100*collapsed['above_one']:.0f}% of draws")
    else:
        say("    one per cluster:          too little data after collapsing")


def main():
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(SEED)

    say("=" * 72)
    say("Clonal clustering: effects recomputed with honest sample sizes")
    say("=" * 72)

    df = add_clusters(cohort.assemble(), cohort.load_mutations())
    say(f"\nSamples: {len(df):,}   distinct clusters: {df.CLUSTER.nunique():,}")

    # ---------------------------------------------------- 1. cluster structure
    say("\n1. How much clustering is there, by group")
    rows = []
    for group in cohort.GROUP_ORDER:
        sub = df[df.GROUP == group]
        if not len(sub):
            continue
        clusters = sub.CLUSTER.nunique()
        rows.append({
            "group": group,
            "isolates": len(sub),
            "clusters": clusters,
            "isolates per cluster": round(len(sub) / clusters, 2),
            "largest cluster": int(sub.CLUSTER.value_counts().iloc[0]),
        })
    say(pd.DataFrame(rows).to_string(index=False))

    say("\n  Largest clusters carrying an Rv0678 variant:")
    carriers = df[df.GROUP.str.startswith("Rv0678")]
    top = carriers.CLUSTER.value_counts().head(10)
    detail = pd.DataFrame({"isolates": top})
    detail["BDQ resistant"] = [int(carriers[carriers.CLUSTER == c].resistant_BDQ.sum()) for c in top.index]
    detail["CFZ resistant"] = [int(carriers[carriers.CLUSTER == c].resistant_CFZ.sum()) for c in top.index]
    say(detail.to_string())

    # ------------------------------------------------- 2. the surviving claim
    say("\n2. The claim that survived site adjustment: does the clofazimine")
    say("   response to losing Rv0678 differ between lineage2 and lineage4?")

    lof_groups = [f"Rv0678 {c}" for c in cohort.LOF_CLASSES]
    intact = df[(~df.mmpL5_LOF) | df.GROUP.eq("reference")].copy()

    for drug in cohort.DRUGS:
        say(f"\n{drug}: Rv0678 loss of function against reference")
        for lineage in MAIN_LINEAGES:
            frame = intact[
                intact.LINEAGE.eq(lineage) & intact.GROUP.isin(lof_groups + ["reference"])
            ].copy()
            frame["EXPOSED"] = frame.GROUP.isin(lof_groups)
            if int(frame.EXPOSED.sum()) < 5:
                continue
            compare(frame, drug, lineage, rng)

    # --------------------------------------------------- 3. substitutions
    say("\n3. The same for Rv0678 substitutions")
    for drug in cohort.DRUGS:
        say(f"\n{drug}: Rv0678 substitution against reference")
        for lineage in MAIN_LINEAGES:
            frame = intact[
                intact.LINEAGE.eq(lineage) & intact.GROUP.isin(["Rv0678 substitution", "reference"])
            ].copy()
            frame["EXPOSED"] = frame.GROUP.eq("Rv0678 substitution")
            if int(frame.EXPOSED.sum()) < 5:
                continue
            compare(frame, drug, lineage, rng)

    # ------------------------------------------------------- 4. pooled effect
    say("\n4. Pooled across lineages, for reference")
    for drug in cohort.DRUGS:
        say(f"\n{drug}:")
        for label, groups in [
            ("Rv0678 loss of function", lof_groups),
            ("Rv0678 substitution", ["Rv0678 substitution"]),
        ]:
            frame = intact[intact.GROUP.isin(groups + ["reference"])].copy()
            frame["EXPOSED"] = frame.GROUP.isin(groups)
            compare(frame, drug, label, rng)

    say("\n\nHow to read these three numbers together.")
    say("")
    say("Cluster-robust standard errors correct the uncertainty, not the point")
    say("estimate. A comparison dominated by one large clone will keep an inflated")
    say("odds ratio here while its confidence interval widens to admit the truth. So")
    say("a wide cluster-robust interval is itself the warning, and the point estimate")
    say("beside it should not be quoted on its own.")
    say("")
    say("The one-per-cluster estimate is the one that moves when a clone is driving")
    say("the result, because it gives that clone a single vote. Where it sits close")
    say("to the isolate-level figure, clonal expansion is not doing the work. Where it")
    say("collapses toward 1 while the isolate-level figure stays large, the effect was")
    say("one event counted many times.")
    say("")
    say("Look at the largest exposed cluster and the cluster count beside every")
    say("comparison. A group of sixty isolates in six clusters is six observations.")

    write_report()
    say(f"\nWritten to {REPORT}")


if __name__ == "__main__":
    main()
