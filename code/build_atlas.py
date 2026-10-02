"""
Per-variant evidence for the bedaquiline and clofazimine atlas.

Run from the project root with the virtual environment active:

    python code/build_atlas.py

One row per distinct mutation, counting the evidence that supports it. A
catalogue grades a mutation from the isolates carrying it, and isolate counts
overstate the evidence whenever a strain has spread: Rv0678 192_ins_g is carried
by 52 isolates in this cohort, 41 of them one outbreak at one site in Peru, so
the independent evidence is a fraction of the isolate count. Every row here
carries its cluster count beside its isolate count, and the cluster count is
what the estimation layer uses.

A cluster is the combination of site, sublineage and the mutation itself, the
same definition as code/cluster_adjust.py. Two unrelated patients at one site
sharing a sublineage and a mutation are merged by that rule, so it understates
independence rather than overstating it.

Only solo samples enter a variant's row: exactly one real major-allele variant
across Rv0678, pepQ and atpE, and nothing uncertain in any of them. A sample
carrying two variants cannot attribute its MIC to either, and a sample whose
gene could not be called cannot be asserted to carry the variant or to lack it.

mmpL5 is never a subject. It carries a real major-allele variant in 53,361 of
54,057 genomes, so every sample solo on an Rv0678 variant also carries one, and
admitting mmpL5 as a subject credits the modifier with the repressor's effect.
It is carried as a per-row covariate instead, as the count of the variant's
isolates whose mmpL5 is disrupted.

A shift is estimated only where the independent evidence supports one. The
threshold is five clusters. Resampling clusters with replacement draws from as
many distinct values as there are clusters, so at three clusters the interval is
built from ten distinct resamples and describes the resampling rather than the
uncertainty. Every variant below the threshold keeps its counts and carries no
estimate, with the reason recorded in the row.

The shift is the variant's fitted mean log2 MIC minus the reference group's,
where the reference group is the 14,187 samples of Section 5.7 that carry no
real major-allele variant in the three genes. The interval covers the variant's
own resampling only. The reference group holds about 14,000 measurements and its
mean carries a standard error near 0.01 doublings, which is an order of
magnitude below the narrowest variant interval here.

Outputs:
  outputs/atlas_evidence.csv  one row per variant
  outputs/atlas_report.txt    the readable summary
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import cohort  # noqa: E402
from cluster_adjust import add_clusters  # noqa: E402
from mic_model import (bootstrap_mu, concentration_series,  # noqa: E402
                       fit_censored_normal, place_mics)

TABLE = Path("outputs/atlas_evidence.csv")
REPORT = Path("outputs/atlas_report.txt")

# Clusters a variant needs before a shift is estimated for it.
MIN_CLUSTERS = 5
# How far a fitted mean may fall outside the tested concentration range, in
# doublings, before the estimate is withheld as extrapolation.
MAX_EXTRAPOLATION = 1.0
BOOTSTRAPS = 400
SEED = 20260101

# Columns printed in the report. The table written to disk carries every column.
INFLATION_COLUMNS = ["gene", "mutation", "class", "isolates", "clusters",
                     "isolates_per_cluster", "largest_cluster", "sites",
                     "BDQ_resistant", "CFZ_resistant"]
EVIDENCE_COLUMNS = ["gene", "mutation", "class", "isolates", "clusters", "sites",
                    "sublineages", "BDQ_isolates", "BDQ_resistant",
                    "CFZ_isolates", "CFZ_resistant"]
SHIFT_COLUMNS = ["gene", "mutation", "class", "isolates", "clusters",
                 "BDQ_shift", "CFZ_shift"]

_lines = []


def say(text=""):
    print(text)
    _lines.append(text)


def subject_variants(mutations, status):
    """The single variant carried by each solo sample, one row per sample.

    Raises if the modifier gene reaches the subject column, which is the fault
    that once credited mmpL5 with the repressor's effect.
    """
    import pandas as pd

    solo = status.index[status.IS_SOLO]
    carried = mutations[
        mutations.REAL_MAJOR
        & mutations.GENE.isin(cohort.BDQ_GENES)
        & mutations.UNIQUEID.isin(solo)
    ]
    per_sample = carried.groupby("UNIQUEID", observed=True).agg(
        GENE=("GENE", "first"), MUTATION=("MUTATION", "first"),
        CLASS=("CLASS", "first"), variants=("REAL_MAJOR", "size"))
    if (per_sample.variants != 1).any():
        raise ValueError("a solo sample carries more than one subject variant")
    if cohort.MODIFIER_GENE in set(per_sample.GENE.astype(str)):
        raise ValueError(f"{cohort.MODIFIER_GENE} reached the subject column")
    return per_sample.drop(columns="variants")


def evidence_table(frame, subjects):
    """One row per variant: how many isolates, and how much of that is independent."""
    import pandas as pd

    rows = frame.join(subjects, how="inner")
    records = []
    for (gene, mutation), chunk in rows.groupby(["GENE", "MUTATION"], observed=True):
        record = {
            "gene": gene,
            "mutation": mutation,
            "class": chunk["CLASS"].iloc[0],
            "isolates": len(chunk),
            "clusters": chunk.CLUSTER.nunique(),
            "sites": chunk.SITEID.nunique(),
            "sublineages": chunk.SUBLINEAGE.nunique(),
            "largest_cluster": int(chunk.CLUSTER.value_counts().max()),
            "mmpL5_disrupted": int(chunk.mmpL5_LOF.sum()),
        }
        for drug in cohort.DRUGS:
            with_mic = chunk[chunk[f"MIC_{drug}"].notna()]
            record[f"{drug}_isolates"] = len(with_mic)
            record[f"{drug}_resistant"] = int(with_mic[f"resistant_{drug}"].sum())
            record[f"{drug}_left_censored"] = int(with_mic[f"censored_left_{drug}"].sum())
            record[f"{drug}_right_censored"] = int(with_mic[f"censored_right_{drug}"].sum())
        records.append(record)
    table = pd.DataFrame(records)
    return table.sort_values(["clusters", "isolates"], ascending=False, kind="stable")


def attach_intervals(frame, series):
    """Add the censoring bounds the fits read, one pair of columns per drug."""
    frame = frame.copy()
    for drug in cohort.DRUGS:
        lowers, uppers, _, off_series = place_mics(
            frame[f"MIC_{drug}"], frame[f"PLATEDESIGN_{drug}"], series, drug)
        if off_series:
            raise ValueError(
                f"{len(off_series)} {drug} MICs are not on the tested series: "
                f"{sorted(set(off_series))[:5]}")
        frame[f"lower_{drug}"] = lowers
        frame[f"upper_{drug}"] = uppers
    return frame


def tested_range(rows, series, drug):
    """The lowest and highest tested concentration, in log2, across the plate
    designs these isolates were measured on. None if no design is known."""
    import numpy as np

    floors, ceilings = [], []
    for design in rows[f"PLATEDESIGN_{drug}"].dropna().unique():
        concentrations = series.get((design, drug))
        if concentrations:
            floors.append(np.log2(min(concentrations)))
            ceilings.append(np.log2(max(concentrations)))
    if not floors:
        return None
    return min(floors), max(ceilings)


def estimate_shifts(frame, subjects, table, series, rng=None):
    """Add a fitted mean, a shift and an interval to every variant the evidence
    supports, and a reason to every variant it does not."""
    import numpy as np
    import pandas as pd

    rng = rng if rng is not None else np.random.default_rng(SEED)
    carriers = frame.join(subjects, how="inner")
    table = table.copy()

    for drug in cohort.DRUGS:
        reference = frame[frame.GROUP.eq("reference")].dropna(subset=[f"lower_{drug}"])
        reference_fit = fit_censored_normal(
            reference[f"lower_{drug}"], reference[f"upper_{drug}"])
        if reference_fit is None:
            raise ValueError(f"the reference distribution for {drug} did not fit")

        means, shifts, lows, highs, reasons = [], [], [], [], []
        for row in table.itertuples():
            rows = carriers[
                carriers.MUTATION.eq(row.mutation) & carriers.GENE.eq(row.gene)
            ].dropna(subset=[f"lower_{drug}"])
            clusters = rows.CLUSTER.nunique()
            fit = (fit_censored_normal(rows[f"lower_{drug}"], rows[f"upper_{drug}"])
                   if clusters >= MIN_CLUSTERS else None)
            limits = tested_range(rows, series, drug)
            if clusters < MIN_CLUSTERS:
                reason = f"{clusters} clusters, below {MIN_CLUSTERS}"
            elif fit is None:
                reason = "the fit did not converge"
            elif limits and fit["mu"] < limits[0] - MAX_EXTRAPOLATION:
                reason = (f"fitted mean {fit['mu']:.2f} lies more than "
                          f"{MAX_EXTRAPOLATION:g} doubling below the lowest tested "
                          f"concentration, {limits[0]:.2f}")
                fit = None
            elif limits and fit["mu"] > limits[1] + MAX_EXTRAPOLATION:
                reason = (f"fitted mean {fit['mu']:.2f} lies more than "
                          f"{MAX_EXTRAPOLATION:g} doubling above the highest tested "
                          f"concentration, {limits[1]:.2f}")
                fit = None
            else:
                reason = ""
            interval = (bootstrap_mu(rows, drug, rng, draws=BOOTSTRAPS)
                        if fit is not None else None)
            if fit is not None and interval is None:
                reason = "the resampled fits did not converge"
            means.append(fit["mu"] if fit else np.nan)
            shifts.append(fit["mu"] - reference_fit["mu"] if fit else np.nan)
            lows.append(interval["low"] - reference_fit["mu"] if interval else np.nan)
            highs.append(interval["high"] - reference_fit["mu"] if interval else np.nan)
            reasons.append(reason)

        table[f"{drug}_reference_mean"] = round(reference_fit["mu"], 3)
        table[f"{drug}_mean"] = np.round(means, 3)
        table[f"{drug}_shift"] = np.round(shifts, 3)
        table[f"{drug}_shift_low"] = np.round(lows, 3)
        table[f"{drug}_shift_high"] = np.round(highs, 3)
        table[f"{drug}_not_estimated"] = reasons
    return table


def main():
    import pandas as pd

    say("=" * 72)
    say("Per-variant evidence")
    say("=" * 72)

    mutations = cohort.load_mutations()
    status = cohort.build_status(mutations)
    series, _, _ = concentration_series()
    frame = attach_intervals(add_clusters(cohort.assemble(status), mutations), series)
    subjects = subject_variants(mutations, status)

    table = estimate_shifts(frame, subjects, evidence_table(frame, subjects), series)
    say(f"\nSolo samples with a MIC: {int(table.isolates.sum()):,}")
    say(f"Distinct variants: {len(table):,}")
    say(f"Variants carried by one isolate: {int((table.isolates == 1).sum()):,}")
    say(f"Variants resting on one cluster: {int((table.clusters == 1).sum()):,}")

    say("\nBy gene:")
    by_gene = table.groupby("gene", observed=True).agg(
        variants=("mutation", "size"), isolates=("isolates", "sum"),
        clusters=("clusters", "sum"))
    say(by_gene.to_string())

    say("\nBy variant class:")
    by_class = table.groupby("class", observed=True).agg(
        variants=("mutation", "size"), isolates=("isolates", "sum"),
        clusters=("clusters", "sum"))
    say(by_class.to_string())

    say("\nWhere isolate counts most overstate the independent evidence:")
    inflated = table[table.isolates >= 5].copy()
    inflated["isolates_per_cluster"] = (inflated.isolates / inflated.clusters).round(2)
    say(inflated.sort_values("isolates_per_cluster", ascending=False, kind="stable")
        .head(12)[INFLATION_COLUMNS].to_string(index=False))

    say("\nVariants resting on the most independent evidence:")
    say(table.head(12)[EVIDENCE_COLUMNS].to_string(index=False))

    say(f"\nShifts, for variants resting on {MIN_CLUSTERS} clusters or more.")
    say("A shift is in doublings of MIC against the reference group, whose fitted")
    say(f"mean log2 MIC is {table.BDQ_reference_mean.iloc[0]} for BDQ and "
        f"{table.CFZ_reference_mean.iloc[0]} for CFZ.")
    estimated = table[table.BDQ_shift.notna() | table.CFZ_shift.notna()]
    if len(estimated):
        shown = estimated[SHIFT_COLUMNS].copy()
        for drug in cohort.DRUGS:
            shown[f"{drug} 95%"] = [
                f"{low:.2f} to {high:.2f}" if low == low else "n/a"
                for low, high in zip(estimated[f"{drug}_shift_low"],
                                     estimated[f"{drug}_shift_high"])]
        say("")
        say(shown.to_string(index=False))
    else:
        say("\n  No variant reaches the threshold.")

    declined = table[table.BDQ_shift.isna() & table.CFZ_shift.isna()]
    say(f"\nVariants carrying no estimate: {len(declined):,} of {len(table):,}, "
        f"of which {int((declined.clusters == 1).sum()):,} rest on one cluster.")
    say("Their counts stand; only the shift is withheld.")

    TABLE.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(TABLE, index=False)
    REPORT.write_text("\n".join(_lines) + "\n")
    print(f"\nTable written to {TABLE} ({len(table):,} variants)")
    print(f"Report written to {REPORT}")


if __name__ == "__main__":
    main()
