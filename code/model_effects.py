"""
Separates the lineage effect from site confounding, and tests whether the
mmpL5-disrupted samples are independent events or one clone counted many times.

Run from the project root with the virtual environment active:

    python code/model_effects.py

Three questions, in order of how much they can invalidate what came before.

First, clonality. Forty-three samples carry both a frameshift Rv0678 and a
disrupted mmpL5, and none of them is resistant. If those forty-three share one
sublineage, one site, and the same two mutations, they are one clonal outbreak
observed forty-three times, and the effective sample size is one. The epistasis
would still be real but the evidence for it would be far weaker than the count
suggests. This is checked by counting distinct combinations rather than
distinct samples.

Second, disentangling lineage from site. The Rv0678 effect looks far stronger in
lineage2 than lineage4, but lineage2 concentrates at the contributing site that
was deliberately enriched for bedaquiline resistance. Computing a site-adjusted
estimate separately within each lineage answers whether the lineage difference
survives once site is held constant, without relying on a model.

Third, a logistic model carrying genotype, mmpL5 status, lineage and site
together, with and without a genotype-by-lineage interaction, tested against
each other. Sparse cells make separation a real risk here, so the fit is
checked for it and the result is reported as secondary to the stratified
estimates rather than in place of them.

Output: outputs/model_report.txt
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import cohort  # noqa: E402

REPORT = Path("outputs/model_report.txt")
MIN_STRATUM = 20
MAIN_LINEAGES = ["lineage1", "lineage2", "lineage3", "lineage4"]

_lines = []


def say(text=""):
    print(text)
    _lines.append(text)


def write_report():
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(_lines) + "\n")


def mh_estimate(frame, drug, by):
    """Mantel-Haenszel odds ratio with `by` held constant. Returns None if unusable."""
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
    if len(tables) < 2:
        return None
    st = StratifiedTable(tables)
    low, high = st.oddsratio_pooled_confint()
    return {
        "strata": len(tables),
        "odds_ratio": st.oddsratio_pooled,
        "ci_low": low,
        "ci_high": high,
        "p_association": st.test_null_odds().pvalue,
        "p_homogeneity": st.test_equal_odds().pvalue,
    }


def main():
    import numpy as np
    import pandas as pd

    say("=" * 72)
    say("Lineage against site, and the clonality question")
    say("=" * 72)

    df = cohort.assemble()
    mutations = cohort.load_mutations()
    sites = cohort.load_sites()
    say(f"\nSamples with genotype status and a MIC: {len(df):,}")

    # ------------------------------------------------------- 1. clonality
    say("\n" + "=" * 72)
    say("1. Are the mmpL5-disrupted samples independent events?")
    say("=" * 72)

    frameshift = df[df.GROUP == "Rv0678 frameshift"]
    disrupted = frameshift[frameshift.mmpL5_LOF]
    say(f"\nRv0678 frameshift samples: {len(frameshift):,}")
    say(f"  of which mmpL5 also disrupted: {len(disrupted):,}")

    def variant_of(ids, gene, lof_only):
        sub = mutations[
            mutations.UNIQUEID.isin(ids) & mutations.GENE.eq(gene) & mutations.REAL_MAJOR
        ]
        if lof_only:
            sub = sub[sub.IS_LOF]
        return sub.groupby("UNIQUEID", observed=True).MUTATION.first()

    for label, ids in [("mmpL5 disrupted", disrupted.index), ("mmpL5 intact", frameshift.index[~frameshift.mmpL5_LOF])]:
        if not len(ids):
            continue
        detail = pd.DataFrame({
            "Rv0678": variant_of(ids, "Rv0678", lof_only=True),
            "mmpL5": variant_of(ids, "mmpL5", lof_only=True),
        }).reindex(ids)
        detail["SUBLINEAGE"] = df.SUBLINEAGE.reindex(ids)
        detail["SITEID"] = df.SITEID.reindex(ids)
        if sites is not None:
            detail["COUNTRY"] = detail.SITEID.map(sites.COUNTRY)

        say(f"\n{label}: {len(detail):,} samples")
        say(f"  distinct Rv0678 mutations:    {detail.Rv0678.nunique()}")
        if detail.mmpL5.notna().any():
            say(f"  distinct mmpL5 mutations:     {detail.mmpL5.nunique()}")
        say(f"  distinct sublineages:         {detail.SUBLINEAGE.nunique()}")
        say(f"  distinct sites:               {detail.SITEID.nunique()}")
        combination = detail[["Rv0678", "mmpL5", "SUBLINEAGE"]].apply(
            lambda row: " | ".join(str(value) for value in row), axis=1
        )
        say(f"  distinct mutation-and-sublineage combinations: {combination.nunique()}")
        say("\n  most frequent combinations:")
        say(cohort.ranked_counts(combination, 8).to_string())
        say("\n  by site:")
        by_site = detail.SITEID.value_counts()
        if sites is not None:
            by_site.index = [f"{s} ({sites.COUNTRY.get(s, '?')})" for s in by_site.index]
        say(by_site.to_string())

    say("\n  If one combination accounts for most of the mmpL5-disrupted samples, the")
    say("  effective sample size is the number of independent events, not the number")
    say("  of isolates, and the epistasis evidence must be stated on that basis.")

    # ------------------------------- 2. lineage effect with site held constant
    say("\n" + "=" * 72)
    say("2. The lineage effect, with site held constant within each lineage")
    say("=" * 72)

    intact = df[~df.mmpL5_LOF | df.GROUP.eq("reference")].copy()
    lof_groups = [f"Rv0678 {c}" for c in cohort.LOF_CLASSES]

    for comparison, groups in [
        ("Rv0678 loss of function", lof_groups),
        ("Rv0678 substitution", ["Rv0678 substitution"]),
    ]:
        say(f"\n{comparison} against reference:")
        for drug in cohort.DRUGS:
            say(f"\n  {drug}:")
            rows = []
            for lineage in MAIN_LINEAGES:
                chunk = intact[intact.LINEAGE.eq(lineage) & intact.GROUP.isin(groups + ["reference"])].copy()
                chunk["EXPOSED"] = chunk.GROUP.isin(groups)
                exposed_n = int(chunk.EXPOSED.sum())
                if exposed_n < 5:
                    rows.append({"lineage": lineage, "exposed n": exposed_n,
                                 "site-adjusted OR": "too few", "95% CI": "", "p": ""})
                    continue
                result = mh_estimate(chunk, drug, "SITEID")
                if result is None:
                    crude_a = int(chunk[chunk.EXPOSED][f"resistant_{drug}"].sum())
                    crude_b = int(chunk[~chunk.EXPOSED][f"resistant_{drug}"].sum())
                    rows.append({
                        "lineage": lineage, "exposed n": exposed_n,
                        "site-adjusted OR": "one site only",
                        "95% CI": f"{crude_a} exposed R, {crude_b} reference R", "p": "",
                    })
                    continue
                rows.append({
                    "lineage": lineage,
                    "exposed n": exposed_n,
                    "site-adjusted OR": round(result["odds_ratio"], 1),
                    "95% CI": f"{result['ci_low']:.1f} to {result['ci_high']:.1f}",
                    "p": f"{result['p_association']:.2g}",
                })
            say(pd.DataFrame(rows).to_string(index=False))

    say("\n  Non-overlapping confidence intervals between lineages mean the difference")
    say("  is not explained by site. Overlapping ones mean it may be.")

    # --------------------------------------------------- 3. logistic model
    say("\n" + "=" * 72)
    say("3. Logistic model with genotype, mmpL5, lineage and site together")
    say("=" * 72)

    import statsmodels.formula.api as smf

    model_data = df[
        df.GROUP.isin(lof_groups + ["Rv0678 substitution", "reference"])
        & df.LINEAGE.isin(MAIN_LINEAGES)
    ].copy()
    model_data["genotype"] = np.where(
        model_data.GROUP.isin(lof_groups), "loss_of_function",
        np.where(model_data.GROUP.eq("Rv0678 substitution"), "substitution", "reference"))
    # Unused categories become all-zero dummy columns, which make the design
    # matrix singular, so they are removed after the reference level is fixed.
    model_data["genotype"] = pd.Categorical(
        model_data.genotype, categories=["reference", "substitution", "loss_of_function"]
    ).remove_unused_categories()
    model_data["lineage"] = pd.Categorical(
        model_data.LINEAGE, categories=["lineage4", "lineage1", "lineage2", "lineage3"]
    ).remove_unused_categories()
    model_data["mmpl5_lof"] = model_data.mmpL5_LOF.astype(int)
    model_data["site"] = model_data.SITEID.astype(str)

    # Sites too small to estimate their own effect are pooled.
    site_counts = model_data.site.value_counts()
    small = site_counts[site_counts < 50].index
    model_data.loc[model_data.site.isin(small), "site"] = "other"
    say(f"\nModel rows: {len(model_data):,}   sites pooled as 'other': {len(small)}")

    for drug in cohort.DRUGS:
        say(f"\n--- {drug} ---")
        model_data["y"] = model_data[f"resistant_{drug}"].astype(int)
        say(f"  resistant: {int(model_data.y.sum()):,} of {len(model_data):,}")

        # A term with no variation cannot be estimated and makes the design
        # matrix singular, so terms are included only if they actually vary.
        terms = ["genotype"]
        skipped = []
        for term, series in [("mmpl5_lof", model_data.mmpl5_lof),
                             ("lineage", model_data.lineage),
                             ("C(site)", model_data.site)]:
            if series.nunique(dropna=True) > 1:
                terms.append(term)
            else:
                skipped.append(term)
        if skipped:
            say(f"  terms dropped for having no variation: {', '.join(skipped)}")

        base = "y ~ " + " + ".join(terms)
        full = base + " + genotype:lineage" if "lineage" in terms else None
        try:
            fit_base = smf.logit(base, data=model_data).fit(disp=False, maxiter=200)
            fit_full = (smf.logit(full, data=model_data).fit(disp=False, maxiter=200)
                        if full else None)
        except Exception as error:  # noqa: BLE001
            say(f"  model did not fit: {error}")
            say("  the stratified estimates in section 2 stand on their own and do not")
            say("  depend on this model fitting")
            continue

        extreme = fit_base.params[fit_base.bse > 5]
        if len(extreme):
            say("\n  WARNING: separation. These terms have standard errors above 5, meaning")
            say("  some cell is all-resistant or all-susceptible and the estimate is not")
            say("  trustworthy. Read the stratified results in section 2 instead:")
            for name in extreme.index:
                say(f"    {name}")

        say("\n  Genotype effects, site and lineage held constant:")
        rows = []
        for name in fit_base.params.index:
            if not name.startswith(("genotype", "mmpl5")):
                continue
            coefficient = fit_base.params[name]
            low, high = fit_base.conf_int().loc[name]
            rows.append({
                "term": name,
                "odds ratio": round(float(np.exp(coefficient)), 2),
                "95% CI": f"{np.exp(low):.2f} to {np.exp(high):.2f}",
                "p": f"{fit_base.pvalues[name]:.2g}",
            })
        say(pd.DataFrame(rows).to_string(index=False))

        if fit_full is None:
            say("\n  lineage does not vary here, so the interaction cannot be tested")
            continue

        from scipy.stats import chi2

        statistic = 2 * (fit_full.llf - fit_base.llf)
        degrees = int(fit_full.df_model - fit_base.df_model)
        p_interaction = chi2.sf(statistic, degrees) if degrees > 0 else float("nan")
        say("\n  Does the genotype effect differ by lineage?")
        say(f"    likelihood ratio statistic {statistic:.2f} on {degrees} degrees of freedom,"
            f" p = {p_interaction:.3g}")
        if p_interaction < 0.05:
            say("    the effect differs by lineage, so no single odds ratio describes it")
        else:
            say("    no evidence the effect differs by lineage once site is accounted for")

    write_report()
    say(f"\nWritten to {REPORT}")


if __name__ == "__main__":
    main()
