"""
Audits the foundations the analysis rests on, before the atlas is built on them.

Run from the project root with the virtual environment active:

    python code/audit_cohort.py

Four questions, each capable of changing a conclusion rather than refining one.

1. Join integrity. 6,525 samples have a UKMYC MIC and no genome. That has been
   treated as sequencing that was never done. The CRyPTIC release notes record
   two historic failures where identifiers in the genetics tables would not
   match the phenotype tables, one involving lab identifiers containing "/" and
   "." that were rewritten in one table and not the other. If some of those
   6,525 are unmatched identifiers rather than absent genomes, the cohort is
   losing real data and losing it non-randomly.

2. Patient-level replication. Sample identifiers carry a subject field and an
   isolate number, so serial isolates from one patient are present. Two isolates
   from one patient are not two observations. This is separate from the clonal
   clustering already corrected for: that concerned one strain spreading between
   patients, this concerns one patient sampled repeatedly, and correcting for
   one does not correct for the other.

3. Phenotype quality. PHENOTYPE_QUALITY is LOW for a substantial share of rows,
   meaning the independent reading methods disagreed and the value reported is
   one reader's judgement. Every effect so far has been computed on all rows
   regardless. If an effect depends on the disputed measurements, that needs to
   be known.

4. Discovery and validation. WGS_SAMPLES carries a DATASET column separating
   CRyPTIC-v1.0, the collection handed to WHO to build the first catalogue, from
   CRyPTIC-v3.0, which CRyPTIC identify as a genuine validation set. Everything
   so far has been exploratory, with many comparisons made and the interesting
   ones reported, which is how findings that do not replicate are produced. This
   reports whether the split is large enough to fix hypotheses on one half and
   test them once on the other.

Output: outputs/audit_report.txt
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import cohort  # noqa: E402
from mic_model import parse_concentration  # noqa: E402
from cluster_adjust import add_clusters  # noqa: E402

REPORT = Path("outputs/audit_report.txt")

_lines = []


def say(text=""):
    print(text)
    _lines.append(text)


def write_report():
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(_lines) + "\n")


def normalise(identifier):
    """Collapse an identifier to a form insensitive to separator rewriting."""
    return re.sub(r"[^a-z0-9]+", "_", str(identifier).lower())


def patient_of(identifier):
    """
    The patient an isolate came from: site plus subject field.

    Identifiers look like site.06.subj.06TB_1276.lab.06MIL2678.iso.2. Subject
    fields are not assumed unique across sites, so the site is included.
    """
    text = str(identifier)
    site = text.split(".")[1] if "." in text else "?"
    match = re.search(r"\.subj\.(.*?)\.lab\.", text)
    subject = match.group(1) if match else text
    return f"{site}|{subject}"


def fisher(a, n1, b, n2):
    from scipy.stats import fisher_exact

    return fisher_exact([[a, n1 - a], [b, n2 - b]])


def pooled_odds_ratio(tables):
    """Mantel-Haenszel odds ratio over 2x2 tables, with its interval.

    Returns None where fewer than two strata carry information, because a pooled
    estimate over one stratum is that stratum.
    """
    import numpy as np
    from statsmodels.stats.contingency_tables import StratifiedTable

    usable = [np.array(table) for table in tables
              if table[0][0] + table[1][0] > 0 and min(sum(row) for row in table) > 0]
    if len(usable) < 2:
        return None
    pooled = StratifiedTable(usable)
    low, high = pooled.oddsratio_pooled_confint()
    return {"odds_ratio": float(pooled.oddsratio_pooled), "low": float(low),
            "high": float(high), "strata": len(usable)}


SCOPE_DRUGS = ("BDQ", "CFZ", "DLM", "LZD")


def censoring_profile(phenotypes, drug, keep=None):
    """How much of a drug's MIC range the plates can measure.

    A left-censored reading says the MIC is below the lowest tested well and a
    right-censored one that it is above the highest, so a drug read mostly at one
    end of its ladder carries little quantitative information. `keep` restricts
    the count to a set of samples.
    """
    rows = phenotypes[phenotypes.DRUG.eq(drug) & phenotypes.MIC.notna()]
    if keep is not None:
        rows = rows[rows.UNIQUEID.isin(keep)]
    if not len(rows):
        return None
    text = rows.MIC.astype(str)
    left, right = text.str.startswith("<="), text.str.startswith(">")
    resistant = rows.BINARY_PHENOTYPE.eq("R")
    return {
        "drug": drug, "MICs": len(rows),
        "left-censored %": round(100 * left.mean(), 1),
        "right-censored %": round(100 * right.mean(), 1),
        "resistant": int(resistant.sum()),
        "resistant at the ceiling %": (round(100 * right[resistant].mean(), 1)
                                       if int(resistant.sum()) else None),
    }


def resistant_wells(layout, drug):
    """Per plate design, the tested concentrations a resistant MIC can land on.

    PLATE_LAYOUT labels every well S or R, so the wells labelled R are the
    measurable range above the breakpoint. A drug with one such well can only
    report resistance as a single value or as off-scale.
    """
    rows = layout[layout.DRUG.eq(drug)]
    counts = {}
    for design, chunk in rows.groupby("PLATEDESIGN", observed=True):
        wells = {}
        for concentration, label in zip(chunk.CONC, chunk.BINARY_PHENOTYPE):
            value, operator = parse_concentration(concentration)
            if value is None or value <= 0:
                continue
            # A ">x" row is the bin above the highest well and repeats that
            # well's concentration, so the concentration counts towards the
            # ladder while the label belongs to the bin and not to the well.
            if operator in (">", ">=") and value in wells:
                continue
            wells[value] = label
        resistant = sorted(value for value, label in wells.items() if label == "R")
        susceptible = sorted(value for value, label in wells.items() if label == "S")
        counts[design] = {
            "wells": len(wells),
            "resistant wells": len(resistant),
            "highest susceptible": susceptible[-1] if susceptible else None,
        }
    return counts


def mic_is_missing(df, drug):
    """Which rows record a plate design, and which of those record no MIC.

    A sample with no row for the drug at all is not a missing measurement of it,
    so the denominator is the rows that carry a design.
    """
    recorded = df[f"PLATEDESIGN_{drug}"].notna()
    return recorded, recorded & df[f"MIC_{drug}"].isna()


def main():
    import numpy as np
    import pandas as pd
    from scipy.stats import chi2_contingency

    say("=" * 72)
    say("Cohort audit")
    say("=" * 72)

    phenotypes = pd.read_parquet(cohort.DATA / "UKMYC_PHENOTYPES.parquet").reset_index()
    genomes = pd.read_parquet(cohort.DATA / "GENOMES.parquet").reset_index()

    # ------------------------------------------------------- 1. join integrity
    say("\n" + "=" * 72)
    say("1. Join integrity: are the unmatched samples absent, or unmatched?")
    say("=" * 72)

    mic_ids = set(phenotypes.UNIQUEID.unique())
    genome_ids = set(genomes.UNIQUEID.unique())
    unmatched = mic_ids - genome_ids
    say(f"\n  samples with a UKMYC MIC:     {len(mic_ids):,}")
    say(f"  samples in GENOMES:           {len(genome_ids):,}")
    say(f"  matched exactly:              {len(mic_ids & genome_ids):,}")
    say(f"  MIC samples with no genome:   {len(unmatched):,}")

    from collections import defaultdict

    mic_normal, genome_normal = defaultdict(set), defaultdict(set)
    for identifier in mic_ids:
        mic_normal[normalise(identifier)].add(identifier)
    for identifier in genome_ids:
        genome_normal[normalise(identifier)].add(identifier)
    recovered = set(mic_normal) & set(genome_normal)
    extra = sorted(
        {i for key in recovered for i in mic_normal[key]} - (mic_ids & genome_ids)
    )

    say(f"\n  matched after collapsing separators: {len(recovered):,}")
    say(f"  additional samples recovered:        {len(extra):,}")

    # Normalising is only safe if it does not merge genuinely distinct samples.
    mic_collisions = len(mic_ids) - len({normalise(i) for i in mic_ids})
    genome_collisions = len(genome_ids) - len({normalise(i) for i in genome_ids})
    say(f"\n  distinct identifiers lost to collision if normalised:"
        f" {mic_collisions} in phenotypes, {genome_collisions} in genomes")
    if mic_collisions or genome_collisions:
        say("  normalising the join would merge distinct samples and must not be done")
    if extra:
        say("\n  WARNING: these matched only after normalising the identifier, so the")
        say("  exact join is dropping real data. Examples of the two spellings:")
        for identifier in extra[:10]:
            partners = sorted(genome_normal[normalise(identifier)])
            say(f"    phenotypes: {identifier}")
            say(f"    genomes:    {', '.join(partners)}")
    else:
        say("  no samples recovered by normalising, so the unmatched samples are")
        say("  genuinely absent from GENOMES rather than spelled differently")

    # ------------------------------------------- 2. patient-level replication
    say("\n" + "=" * 72)
    say("2. Patient-level replication: how many isolates per patient?")
    say("=" * 72)

    df = add_clusters(cohort.assemble(), cohort.load_mutations())
    df["PATIENT"] = [patient_of(i) for i in df.index]

    isolates, patients = len(df), df.PATIENT.nunique()
    say(f"\n  isolates in the analysis cohort: {isolates:,}")
    say(f"  distinct patients:               {patients:,}")
    say(f"  isolates per patient:            {isolates / patients:.3f}")

    counts = df.PATIENT.value_counts()
    repeated = counts[counts > 1]
    say(f"  patients contributing more than one isolate: {len(repeated):,}")
    say(f"  isolates belonging to such patients:         {int(repeated.sum()):,}"
        f" ({100 * repeated.sum() / isolates:.1f}% of the cohort)")
    if len(repeated):
        say("\n  distribution of isolates per patient:")
        say(counts.value_counts().sort_index().to_string())

    say("\n  By group:")
    rows = []
    for group in cohort.GROUP_ORDER:
        sub = df[df.GROUP == group]
        if not len(sub):
            continue
        rows.append({
            "group": group,
            "isolates": len(sub),
            "patients": sub.PATIENT.nunique(),
            "clusters": sub.CLUSTER.nunique(),
            "isolates per patient": round(len(sub) / sub.PATIENT.nunique(), 2),
        })
    say(pd.DataFrame(rows).to_string(index=False))
    say("\n  Clusters and patients measure different things. A cluster is one strain")
    say("  spreading between people; a repeated patient is one person sampled twice.")
    say("  An effective sample size has to account for both.")

    # --------------------------------- effect sensitivity to one isolate per patient
    lof_groups = [f"Rv0678 {c}" for c in cohort.LOF_CLASSES]
    intact = df[(~df.mmpL5_LOF) | df.GROUP.eq("reference")]

    say("\n  Key effects with one isolate per patient, keeping the first of each:")
    single = intact.sort_index().groupby("PATIENT", observed=True).head(1)
    for drug in cohort.DRUGS:
        say(f"\n  {drug}:")
        rows = []
        for label, groups in [("Rv0678 loss of function", lof_groups),
                              ("Rv0678 substitution", ["Rv0678 substitution"])]:
            for frame, name in [(intact, "all isolates"), (single, "one per patient")]:
                exposed = frame[frame.GROUP.isin(groups)]
                reference = frame[frame.GROUP == "reference"]
                if not len(exposed):
                    continue
                a, n1 = int(exposed[f"resistant_{drug}"].sum()), len(exposed)
                b, n2 = int(reference[f"resistant_{drug}"].sum()), len(reference)
                odds, p = fisher(a, n1, b, n2)
                rows.append({
                    "comparison": label, "cohort": name,
                    "exposed n": n1, "resistant": a,
                    "reference n": n2, "odds ratio": round(odds, 1),
                    "p": f"{p:.2g}",
                })
        say(pd.DataFrame(rows).to_string(index=False))

    # ------------------------------------------------------ 3. phenotype quality
    say("\n" + "=" * 72)
    say("3. Phenotype quality: do the effects depend on disputed measurements?")
    say("=" * 72)

    for drug in cohort.DRUGS:
        quality = df[f"PHENOTYPE_QUALITY_{drug}"]
        say(f"\n  {drug} quality in the cohort:")
        say(quality.value_counts(dropna=False).to_string())

        high = intact[intact[f"PHENOTYPE_QUALITY_{drug}"] == "HIGH"]
        say(f"\n  {drug}, effects restricted to HIGH quality only:")
        rows = []
        for label, groups in [("Rv0678 loss of function", lof_groups),
                              ("Rv0678 substitution", ["Rv0678 substitution"])]:
            for frame, name in [(intact, "all quality"), (high, "HIGH only")]:
                exposed = frame[frame.GROUP.isin(groups)]
                reference = frame[frame.GROUP == "reference"]
                if not len(exposed) or not len(reference):
                    continue
                a, n1 = int(exposed[f"resistant_{drug}"].sum()), len(exposed)
                b, n2 = int(reference[f"resistant_{drug}"].sum()), len(reference)
                odds, p = fisher(a, n1, b, n2)
                rows.append({
                    "comparison": label, "cohort": name,
                    "exposed n": n1, "resistant": a,
                    "reference n": n2, "odds ratio": round(odds, 1),
                    "p": f"{p:.2g}",
                })
        say(pd.DataFrame(rows).to_string(index=False))

    # ------------------------------------------- 4. discovery and validation
    say("\n" + "=" * 72)
    say("4. Discovery and validation: is the held-out set usable?")
    say("=" * 72)

    path = cohort.DATA / "WGS_SAMPLES.parquet"
    if not path.is_file():
        say(f"\n  {path} not found, so the split cannot be assessed")
        write_report()
        return

    try:
        wgs = pd.read_parquet(path)
    except Exception as error:  # noqa: BLE001
        size_mb = path.stat().st_size / 1_000_000
        say(f"\n  {path} could not be read: {error}")
        say(f"  the file on disk is {size_mb:.2f} MB; the Zenodo record lists 9.4 MB")
        say("  A truncated or failed download produces exactly this. Re-download the")
        say("  file and verify its MD5 against the Zenodo record before rerunning.")
        write_report()
        return
    if "UNIQUEID" not in wgs.columns:
        wgs = wgs.reset_index()
    say(f"\n  WGS_SAMPLES columns: {list(wgs.columns)}")

    dataset_column = next(
        (c for c in wgs.columns if c.upper() == "DATASET"),
        None,
    )
    if dataset_column is None:
        say("\n  no DATASET column found. The split cannot be made from this table.")
        say("  Candidate columns carrying few distinct values:")
        for column in wgs.columns:
            distinct = wgs[column].nunique(dropna=True)
            if 1 < distinct <= 10:
                say(f"    {column}: {sorted(wgs[column].dropna().unique())[:10]}")
        write_report()
        return

    say(f"\n  DATASET across all WGS samples:")
    say(wgs[dataset_column].value_counts(dropna=False).to_string())

    assignment = wgs.drop_duplicates("UNIQUEID").set_index("UNIQUEID")[dataset_column]
    df["DATASET"] = assignment.reindex(df.index)
    say(f"\n  DATASET within the analysis cohort ({len(df):,} samples):")
    say(df.DATASET.value_counts(dropna=False).to_string())

    say("\n  Group sizes by dataset:")
    say(pd.crosstab(df.GROUP, df.DATASET.fillna("unassigned")).to_string())

    say("\n  Resistance by dataset, reference group only, as a sanity check that the")
    say("  halves are comparable:")
    reference = df[df.GROUP == "reference"]
    rows = []
    for name, chunk in reference.groupby(df.DATASET.fillna("unassigned"), observed=True):
        rows.append({
            "dataset": name, "n": len(chunk),
            "BDQ % R": round(100 * chunk.resistant_BDQ.mean(), 2),
            "CFZ % R": round(100 * chunk.resistant_CFZ.mean(), 2),
        })
    say(pd.DataFrame(rows).to_string(index=False))

    say("\n  Whether a held-out test is possible depends on the exposed groups having")
    say("  enough samples in the smaller half. The crosstab above is the answer.")

    # ------------------------------------------- 5. rows that carry no MIC
    say("\n" + "=" * 72)
    say("5. Rows carrying no MIC, and whether the loss is random")
    say("=" * 72)

    worst_site = None
    for drug in cohort.DRUGS:
        recorded, missing = mic_is_missing(df, drug)
        say(f"\n  {drug}: {int(recorded.sum()):,} rows carry a plate design, "
            f"{int(missing.sum())} of them no MIC. "
            f"{int((~recorded).sum())} samples carry no {drug} row at all.")

        with_mic = df[df[f"MIC_{drug}"].notna()]
        left = int(with_mic[f"censored_left_{drug}"].sum())
        right = int(with_mic[f"censored_right_{drug}"].sum())
        say(f"    of the {len(with_mic):,} rows carrying a MIC, {left:,} are "
            f"left-censored, {100 * left / len(with_mic):.2f}%, and {right} are "
            f"right-censored")

        by_site = pd.DataFrame({
            "rows": recorded.groupby(df.SITEID).sum(),
            "missing": missing.groupby(df.SITEID).sum(),
        })
        by_site["percent"] = (100 * by_site.missing / by_site.rows).round(2)
        say("")
        say(by_site.sort_values("percent", ascending=False).head(4).to_string())
        worst_site = by_site.percent.idxmax() if worst_site is None else worst_site
        here = df.SITEID.eq(worst_site)
        a, n1 = int(missing[here].sum()), int(recorded[here].sum())
        b, n2 = int(missing[~here].sum()), int(recorded[~here].sum())
        odds, p = fisher(a, n1, b, n2)
        say(f"    site {worst_site}: {a} of {n1:,}, {100 * a / n1:.2f}%, against "
            f"{b} of {n2:,} elsewhere, {100 * b / n2:.2f}%, odds ratio {odds:.1f}, "
            f"p = {p:.2g}")

        say("\n    By plate design:")
        tables = []
        for design, chunk in df[recorded].groupby(df[f"PLATEDESIGN_{drug}"][recorded],
                                                  observed=True):
            lost = int(missing[chunk.index].sum())
            say(f"      {design}: {lost} of {len(chunk):,}, "
                f"{100 * lost / len(chunk):.2f}%")
        for site, chunk in df[recorded].groupby(df.SITEID[recorded], observed=True):
            designs = chunk[f"PLATEDESIGN_{drug}"].unique()
            if not {"UKMYC5", "UKMYC6"} <= set(designs):
                continue
            counts = []
            for design in ("UKMYC5", "UKMYC6"):
                part = chunk[chunk[f"PLATEDESIGN_{drug}"].eq(design)]
                lost = int(missing[part.index].sum())
                counts.append([lost, len(part) - lost])
            tables.append(counts)
        pooled = pooled_odds_ratio(tables)
        if pooled:
            say(f"      UKMYC5 against UKMYC6, SITEID held constant across "
                f"{pooled['strata']} sites that ran both: "
                f"{pooled['odds_ratio']:.2f} ({pooled['low']:.2f} to "
                f"{pooled['high']:.2f})")

        say("\n    By genotype status:")
        status = pd.Series("carrier", index=df.index)
        status[df.GROUP.eq("reference")] = "reference"
        status[df.GROUP.eq("uncertain")] = "uncertain"
        for label, excluded in (("all sites", df.SITEID.ne(df.SITEID)),
                                (f"without site {worst_site}", here)):
            keep = recorded & ~excluded
            table = []
            for name in ("uncertain", "reference", "carrier"):
                part = keep & status.eq(name)
                lost = int(missing[part].sum())
                table.append([lost, int(part.sum()) - lost])
                say(f"      {label}, {name}: {lost} of {int(part.sum()):,}, "
                    f"{100 * lost / max(int(part.sum()), 1):.2f}%")
            statistic, p, _, _ = chi2_contingency(np.array(table))
            say(f"      {label}: chi-square {statistic:.1f}, p = {p:.3g}")

    # ------------------------------------------------ 6. why these two drugs
    say("\n" + "=" * 72)
    say("6. What the plates can measure, for the four drugs in this axis")
    say("=" * 72)
    say("")
    rows = [censoring_profile(phenotypes, drug, keep=set(df.index))
            for drug in SCOPE_DRUGS]
    say(pd.DataFrame([row for row in rows if row]).to_string(index=False))
    say("")
    layout = pd.read_parquet(cohort.DATA / "PLATE_LAYOUT.parquet").reset_index()
    for drug in SCOPE_DRUGS:
        for design, counts in resistant_wells(layout, drug).items():
            say(f"  {drug} {design}: {counts['wells']} tested concentrations, "
                f"{counts['resistant wells']} of them above the breakpoint, "
                f"highest susceptible well {counts['highest susceptible']}")

    write_report()
    say(f"\nWritten to {REPORT}")


if __name__ == "__main__":
    main()
