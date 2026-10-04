"""
The samples excluded from the core analysis as uncertain, analysed on their own.

Run from the project root with the virtual environment active:

    python code/heteroresistance.py

The cohort definitions in code/cohort.py exclude a sample from both the
reference group and the solo groups when a target gene carries a null call, a
het call or a minor allele, because none of the three can be asserted to carry
the variant or to lack it. That exclusion is correct for attribution and it
pools two different things. A null call is a position that could not be read,
which is missing data. A het call or a minor indel is a variant detected in a
minority of reads, which is a measurement rather than its absence.

This module separates them and estimates what the second group shows. Three
analyses:

  The split. Every excluded sample classified by whether its uncertainty is a
  null call, a detected minor allele, or both, per gene, with resistance rates
  beside each.

  The resolved classes. MINOR_MUTATION carries the resolved form of a minor
  allele, so `141_minorindel` reads as `141_ins_c` and `C46Z` as `C46G`. The
  resolved form is parsed by the project's own parser and classified by the
  same rules as a major allele, and the MIC shift against the reference group
  is estimated per class with intervals from resampling clusters.

  The read fraction. FRS gives the fraction of reads supporting the minor
  allele. The mean log2 MIC is fitted as a line in FRS with a common standard
  deviation, over the same censoring intervals, so the slope is in doublings of
  MIC per unit of read fraction.

  More than one allele. A sample carrying two minor alleles at the subject gene
  cannot attribute its MIC to either, so eligible_samples excludes it. Nothing
  about such a sample is wild type and no major allele explains it, so the group
  is measured on its own: what it holds, how much resistance it carries, and its
  MIC shift from the same joint fit.

  Site held constant. The minor-allele groups are small and unevenly spread
  across the collection sites, and the sites differ in mean MIC, so a group
  mean confounds the two. Every group is refitted jointly with the reference
  group in one censored regression, once without site indicators and once with
  them. The two share a reference mean and a standard deviation, so what
  separates them is the site adjustment alone.

Eligibility for the second and third analyses is strict, because the point is
attribution: no real major-allele variant in Rv0678, pepQ or atpE, no null call
in any of them, no minor allele outside Rv0678, and exactly one minor allele at
Rv0678. A sample carrying two minor alleles at Rv0678 cannot attribute its MIC
to either.

Everything here is exploratory. The group was chosen for analysis after the
record showed it carrying more resistance than the reference group and no
comparison is pre-registered. These are estimates with intervals, not tested
hypotheses.

Outputs:
  outputs/heteroresistance_report.txt
  outputs/heteroresistance_estimates.csv
  outputs/multi_allele_counts.csv

The adjusted fit carries the module's running cost. It refits a regression of
sixteen columns over fourteen thousand rows once per bootstrap draw, per drug.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import cohort  # noqa: E402
from garc import parse_frame  # noqa: E402
from cluster_adjust import add_clusters  # noqa: E402
from mic_model import (concentration_series, fit_censored_linear,  # noqa: E402
                       fit_censored_normal, log_interval_mass, place_mics)

REPORT = Path("outputs/heteroresistance_report.txt")
ESTIMATES = Path("outputs/heteroresistance_estimates.csv")
MULTI = Path("outputs/multi_allele_counts.csv")

SUBJECT_GENE = "Rv0678"
BOOTSTRAPS = 400
SEED = 20260101
# Isolates a group needs before its shift or slope is estimated.
MIN_ISOLATES = 12
# Isolates a site needs before it gets its own indicator in the adjusted fit.
# Smaller sites are pooled, so that a handful of rows cannot claim a column.
SITE_MINIMUM = 30

# The columns the parser reads. The resolved form replaces MUTATION, so the
# frame must be rebuilt from these rather than reusing an already parsed one.
RAW_COLUMNS = [
    "UNIQUEID", "GENE", "MUTATION", "GENE_POSITION", "REF", "ALT",
    "NUCLEOTIDE_NUMBER", "NUCLEOTIDE_INDEX", "CODES_PROTEIN", "INDEL_LENGTH",
    "INDEL_NUCLEOTIDES", "AMINO_ACID_NUMBER", "AMINO_ACID_SEQUENCE",
    "NUMBER_NUCLEOTIDE_CHANGES", "IS_NULL", "IS_MINOR", "MINOR_MUTATION",
    "MINOR_READS", "COVERAGE", "FRS", cohort.GENE_CODES,
]

_lines = []


def say(text=""):
    print(text)
    _lines.append(text)


def stream(name):
    """A generator seeded by name.

    Each estimate draws from its own stream, so its interval does not depend on
    which other estimates were computed before it and can be reproduced on its
    own.
    """
    import numpy as np

    return np.random.default_rng([SEED] + [ord(letter) for letter in name])


def call_types(mutations):
    """Per sample, whether its uncertainty is a null call, a detected minor
    allele, or both, and which genes it touches."""
    import numpy as np
    import pandas as pd

    rows = mutations[mutations.GENE.isin(cohort.BDQ_GENES) & mutations.UNCERTAIN]
    per_sample = rows.groupby("UNIQUEID", observed=True).agg(
        null=("IS_NULL_CALL", "any"),
        detected=("IS_MINOR", "any"),
        genes=("GENE", lambda values: ",".join(sorted(set(map(str, values))))))
    per_sample["detected"] = per_sample.detected | rows.groupby(
        "UNIQUEID", observed=True).IS_HET_CALL.any()
    per_sample["TYPE"] = np.select(
        [per_sample.detected & ~per_sample.null,
         per_sample.null & ~per_sample.detected],
        ["minor allele only", "null only"], "both")
    return per_sample


def eligible_samples(mutations):
    """Samples whose only finding in the three genes is one minor allele at the
    subject gene."""
    target = mutations[mutations.GENE.isin(cohort.BDQ_GENES)]
    major = set(target[target.REAL_MAJOR].UNIQUEID)
    nulls = set(target[target.IS_NULL_CALL].UNIQUEID)
    minor = target[target.IS_HET_CALL | target.IS_MINOR]
    elsewhere = set(minor[~minor.GENE.eq(SUBJECT_GENE)].UNIQUEID)
    subject = minor[minor.GENE.eq(SUBJECT_GENE)]
    counts = subject.groupby("UNIQUEID", observed=True).size()
    keep = {unique_id for unique_id, count in counts.items()
            if count == 1 and unique_id not in major and unique_id not in nulls
            and unique_id not in elsewhere}
    return keep, int((counts > 1).sum())


def resolve(mutations, keep):
    """The resolved form of each eligible sample's minor allele, parsed and
    classified by the same rules as a major allele."""
    import numpy as np

    minor = mutations[mutations.GENE.eq(SUBJECT_GENE)
                      & (mutations.IS_HET_CALL | mutations.IS_MINOR)
                      & mutations.UNIQUEID.isin(keep)
                      & mutations.MINOR_MUTATION.notna()]
    raw = minor[[column for column in RAW_COLUMNS if column in minor.columns]].copy()
    raw["MUTATION"] = raw.MINOR_MUTATION.astype(str)
    parsed = cohort.classify(parse_frame(raw))
    unparsed = int((~parsed.PARSED).sum())
    parsed = parsed[parsed.PARSED].set_index("UNIQUEID")
    parsed["MINOR_GROUP"] = np.where(
        parsed.CLASS.isin(cohort.LOF_CLASSES), "loss of function",
        np.where(parsed.AFFECTS.eq("PROM"), "promoter", "substitution"))
    return parsed, unparsed


def multi_allele_samples(mutations):
    """Samples whose only finding across the three genes is more than one minor
    allele at the subject gene.

    eligible_samples admits exactly one, because a sample carrying two cannot
    attribute its MIC to either. The group is still worth measuring: nothing
    about it is wild type, and no major allele explains it.
    """
    target = mutations[mutations.GENE.isin(cohort.BDQ_GENES)]
    major = set(target[target.REAL_MAJOR].UNIQUEID)
    nulls = set(target[target.IS_NULL_CALL].UNIQUEID)
    minor = target[target.IS_HET_CALL | target.IS_MINOR]
    elsewhere = set(minor[~minor.GENE.eq(SUBJECT_GENE)].UNIQUEID)
    subject = minor[minor.GENE.eq(SUBJECT_GENE)]
    counts = subject.groupby("UNIQUEID", observed=True).size()
    keep = {unique_id for unique_id, count in counts.items()
            if count > 1 and unique_id not in major and unique_id not in nulls
            and unique_id not in elsewhere}
    return keep, int((counts > 1).sum())


def multi_allele_profile(mutations, keep):
    """Per sample, how many minor alleles it carries at the subject gene, how
    many distinct positions they sit at, how many are loss of function, and
    what their read fractions sum to."""
    import pandas as pd

    columns = ["alleles", "positions", "loss_of_function", "read_fraction_sum",
               "read_fraction_low", "read_fraction_high"]
    minor = mutations[mutations.GENE.eq(SUBJECT_GENE)
                      & (mutations.IS_HET_CALL | mutations.IS_MINOR)
                      & mutations.UNIQUEID.isin(keep)
                      & mutations.MINOR_MUTATION.notna()]
    if minor.empty:
        # A cohort with no multi-allele carrier is a cohort, not an error, and
        # the parser returns none of the columns to classify when given no rows.
        empty = pd.DataFrame(columns=columns, index=pd.Index([], name="UNIQUEID"))
        return empty, 0
    raw = minor[[column for column in RAW_COLUMNS if column in minor.columns]].copy()
    raw["MUTATION"] = raw.MINOR_MUTATION.astype(str)
    parsed = cohort.classify(parse_frame(raw))
    unparsed = int((~parsed.PARSED).sum())
    parsed = parsed[parsed.PARSED]
    parsed = parsed.assign(LOF=parsed.CLASS.isin(cohort.LOF_CLASSES))
    profile = parsed.groupby("UNIQUEID", observed=True).agg(
        alleles=("MUTATION", "size"),
        positions=("POSITION", "nunique"),
        loss_of_function=("LOF", "sum"),
        read_fraction_sum=("FRS", "sum"),
        read_fraction_low=("FRS", "min"),
        read_fraction_high=("FRS", "max"))
    return profile, unparsed


def multi_allele_rows(frame, profile, drug):
    """What the group holds, and how its resistance compares with the reference
    group. Counts only; the MIC shift comes from the joint fit."""
    carriers = frame[frame.MULTI_ALLELES.notna()].dropna(subset=[f"lower_{drug}"])
    reference = frame[frame.GROUP.eq("reference")].dropna(subset=[f"lower_{drug}"])
    rows = []
    for label, subset in (("all", carriers),
                          ("two alleles", carriers[carriers.MULTI_ALLELES.eq(2)]),
                          ("three or more", carriers[carriers.MULTI_ALLELES.ge(3)])):
        resistant = int(subset[f"resistant_{drug}"].sum())
        rows.append({
            "drug": drug, "subset": label, "isolates": len(subset),
            "clusters": subset.CLUSTER.nunique(),
            "resistant": resistant,
            "percent": round(100 * resistant / len(subset), 1) if len(subset) else None,
            "reference_percent": round(
                100 * reference[f"resistant_{drug}"].mean(), 2)})
    return rows


def fit_slope(lower, upper, fraction):
    """Mean log2 MIC as a line in the read fraction, with a common standard
    deviation, over the censoring intervals. The slope is in doublings of MIC
    per unit of read fraction."""
    import numpy as np
    from scipy.optimize import minimize

    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    fraction = np.asarray(fraction, dtype=float)
    finite = np.concatenate([lower[np.isfinite(lower)], upper[np.isfinite(upper)]])
    if not len(finite):
        return None
    start = np.array([float(finite.mean()), 0.0,
                      np.log(max(float(finite.std()), 0.5))])

    def negative_log_likelihood(parameters):
        mu = parameters[0] + parameters[1] * fraction
        sigma = np.exp(parameters[2])
        return -np.sum(log_interval_mass((lower - mu) / sigma, (upper - mu) / sigma))

    result = minimize(
        negative_log_likelihood, start, method="L-BFGS-B",
        bounds=[(-25, 25), (-30, 30), (np.log(0.05), np.log(20))])
    if not result.success:
        return None
    return {"intercept": float(result.x[0]), "slope": float(result.x[1]),
            "sigma": float(np.exp(result.x[2]))}


def resample_clusters(frame, rng, draws, statistic):
    """Percentile interval for a statistic, resampling clusters with
    replacement."""
    import numpy as np
    import pandas as pd

    groups = {name: chunk for name, chunk in frame.groupby("CLUSTER", observed=True)}
    names = list(groups)
    values = []
    for _ in range(draws):
        picked = rng.choice(names, size=len(names), replace=True)
        value = statistic(pd.concat([groups[name] for name in picked]))
        if value is not None:
            values.append(value)
    if len(values) < draws // 4:
        return None
    return float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))


def prepare(mutations, status):
    """The cohort frame with censoring bounds, the resolved minor allele, its
    read fraction, and a cluster key."""
    frame = add_clusters(cohort.assemble(status), mutations)
    series, _, _ = concentration_series()
    for drug in cohort.DRUGS:
        lower, upper, _, off_series = place_mics(
            frame[f"MIC_{drug}"], frame[f"PLATEDESIGN_{drug}"], series, drug)
        if off_series:
            raise ValueError(f"{drug} MICs off the tested series: {off_series[:5]}")
        frame[f"lower_{drug}"] = lower
        frame[f"upper_{drug}"] = upper

    keep, multiple = eligible_samples(mutations)
    resolved, unparsed = resolve(mutations, keep)
    several, _ = multi_allele_samples(mutations)
    profile, multi_unparsed = multi_allele_profile(mutations, several)
    frame["MULTI_ALLELES"] = profile.alleles.reindex(frame.index)
    frame["MULTI_FRS_SUM"] = profile.read_fraction_sum.reindex(frame.index)
    frame["MINOR_GROUP"] = resolved.MINOR_GROUP.reindex(frame.index)
    frame["MINOR_FORM"] = resolved.MUTATION.reindex(frame.index)
    frame["MINOR_FRS"] = resolved.FRS.reindex(frame.index)
    # add_clusters keys a major-allele carrier on site, sublineage and its
    # defining mutation. A minor-allele carrier has no defining major mutation,
    # so it is keyed on the resolved form the same way. Everything else keeps
    # the cluster of its own index.
    resolved_key = (frame.SITEID.astype(str) + " | " + frame.SUBLINEAGE.astype(str)
                    + " | " + frame.MINOR_FORM.astype(str))
    frame["CLUSTER"] = resolved_key.where(frame.MINOR_FORM.notna(), frame.CLUSTER)
    return frame, {"eligible": len(keep), "multiple": multiple, "unparsed": unparsed,
                   "several": len(several), "several_unparsed": multi_unparsed,
                   "profile": profile}


def shift_rows(frame, drug, rng):
    """MIC shift against the reference group for every minor-allele class, and
    for the major-allele classes on the same cohort and filters."""
    reference = frame[frame.GROUP.eq("reference")].dropna(subset=[f"lower_{drug}"])
    reference_fit = fit_censored_normal(
        reference[f"lower_{drug}"], reference[f"upper_{drug}"])
    if reference_fit is None:
        raise ValueError(f"the reference distribution for {drug} did not fit")

    major = {f"major {label}": groups for label, groups in (
        ("loss of function", [f"Rv0678 {c}" for c in cohort.LOF_CLASSES]),
        ("substitution", ["Rv0678 substitution"]))}
    rows = []
    for label in ("loss of function", "substitution", "promoter"):
        subset = frame[frame.MINOR_GROUP.eq(label)].dropna(subset=[f"lower_{drug}"])
        rows.append((f"minor {label}", subset))
    for label, groups in major.items():
        rows.append((label, frame[frame.GROUP.isin(groups)]
                     .dropna(subset=[f"lower_{drug}"])))

    records = []
    for label, subset in rows:
        record = {"drug": drug, "group": label, "isolates": len(subset),
                  "clusters": subset.CLUSTER.nunique(),
                  "resistant": int(subset[f"resistant_{drug}"].sum()),
                  "reference_mean": round(reference_fit["mu"], 3)}
        fit = (fit_censored_normal(subset[f"lower_{drug}"], subset[f"upper_{drug}"])
               if len(subset) >= MIN_ISOLATES else None)
        if fit is None:
            record.update(shift=None, shift_low=None, shift_high=None,
                          withheld=f"{len(subset)} isolates, below {MIN_ISOLATES}")
        else:
            interval = resample_clusters(
                subset, rng, BOOTSTRAPS,
                lambda sample, d=drug: (
                    lambda f: f["mu"] if f else None)(
                    fit_censored_normal(sample[f"lower_{d}"], sample[f"upper_{d}"])))
            record.update(
                shift=round(fit["mu"] - reference_fit["mu"], 3),
                shift_low=(round(interval[0] - reference_fit["mu"], 3)
                           if interval else None),
                shift_high=(round(interval[1] - reference_fit["mu"], 3)
                            if interval else None),
                withheld="")
        records.append(record)
    return records


def group_masks(frame):
    """The carrier groups the joint fit estimates, as boolean masks.

    The reference group is the baseline and carries no column. Minor and major
    carriers are disjoint by construction, because eligible_samples admits no
    sample carrying a real major variant in the three genes.
    """
    masks = {}
    for label in ("loss of function", "substitution", "promoter"):
        masks[f"minor {label}"] = frame.MINOR_GROUP.eq(label)
    masks["major loss of function"] = frame.GROUP.isin(
        [f"{SUBJECT_GENE} {class_}" for class_ in cohort.LOF_CLASSES])
    masks["major substitution"] = frame.GROUP.eq(f"{SUBJECT_GENE} substitution")
    masks["minor, two or more alleles"] = frame.MULTI_ALLELES.notna()
    return masks


def site_levels(frame):
    """Site labels for the design, with sites below SITE_MINIMUM pooled."""
    sites = frame.SITEID.astype(str)
    counts = sites.value_counts()
    small = set(counts[counts < SITE_MINIMUM].index)
    return sites.where(~sites.isin(small), "other")


def design(frame, masks, sites=None):
    """Intercept, one column per group, and one indicator per site beyond the
    first. Returns the matrix and its column names."""
    import numpy as np

    columns = [np.ones(len(frame))]
    names = ["intercept"]
    for label, mask in masks.items():
        columns.append(mask.to_numpy(dtype=float))
        names.append(label)
    if sites is not None:
        for level in sorted(sites.unique())[1:]:
            columns.append(sites.eq(level).to_numpy(dtype=float))
            names.append(f"site {level}")
    return np.column_stack(columns), names


def cluster_positions(clusters):
    """Row positions of each cluster, for resampling arrays rather than
    concatenating frames."""
    import numpy as np
    import pandas as pd

    codes = pd.Categorical(clusters).codes
    order = np.argsort(codes, kind="stable")
    bounds = np.searchsorted(codes[order], np.arange(codes.max() + 2))
    return [order[start:stop] for start, stop in zip(bounds[:-1], bounds[1:])]


def adjusted_rows(frame, drug, rng, draws=None):
    """Every group's MIC shift against the reference group, fitted jointly with
    it in one censored regression, once plain and once with site held constant.

    The two fits share a reference mean and a standard deviation and differ
    only in the site indicators, so the distance between them is the site
    adjustment. Both differ from shift_rows, which fits a separate standard
    deviation for each group against a reference fitted on its own.
    """
    import numpy as np
    import pandas as pd

    draws = BOOTSTRAPS if draws is None else draws
    masks = group_masks(frame)
    carried = frame.GROUP.eq("reference")
    for mask in masks.values():
        carried = carried | mask
    available = frame[carried & frame[f"lower_{drug}"].notna()]
    masks = {label: mask.loc[available.index] for label, mask in masks.items()}
    stacked = np.column_stack([mask.to_numpy() for mask in masks.values()])
    if stacked.sum(axis=1).max(initial=0) > 1:
        raise ValueError("a sample falls in more than one group")

    sizes = {label: int(mask.sum()) for label, mask in masks.items()}
    withheld = [label for label in masks if sizes[label] < MIN_ISOLATES]
    # A group too small to carry a column is dropped from the fit rather than
    # left in it, because a carrier without a column of its own would be
    # counted as a reference isolate and pull the reference mean toward itself.
    keep = ~np.logical_or.reduce(
        [masks[label].to_numpy() for label in withheld], axis=0
    ) if withheld else np.ones(len(available), dtype=bool)
    rows = available[keep]
    fitted = {label: mask[keep] for label, mask in masks.items()
              if label not in withheld}
    lower = rows[f"lower_{drug}"].to_numpy(dtype=float)
    upper = rows[f"upper_{drug}"].to_numpy(dtype=float)
    sites = site_levels(rows)
    matrices = {"joint": design(rows, fitted),
                "site-adjusted": design(rows, fitted, sites)}

    point = {}
    for kind, (matrix, names) in matrices.items():
        fit = fit_censored_linear(lower, upper, matrix)
        if fit is None:
            raise ValueError(f"the {kind} fit for {drug} did not converge")
        point[kind] = dict(zip(names, fit["beta"]))

    groups = cluster_positions(rows.CLUSTER)
    membership = {label: mask.to_numpy() for label, mask in fitted.items()}
    drawn = {kind: {label: [] for label in fitted} for kind in matrices}
    for _ in range(draws):
        picked = rng.integers(0, len(groups), size=len(groups))
        index = np.concatenate([groups[position] for position in picked])
        # A draw that misses every isolate of a group leaves that group's
        # coefficient unidentified, and the optimiser returns its starting
        # value. The draw carries no information about that group, so it is
        # dropped for that group rather than entering the interval as a zero.
        present = [label for label in fitted if membership[label][index].any()]
        for kind, (matrix, names) in matrices.items():
            fit = fit_censored_linear(lower[index], upper[index], matrix[index])
            if fit is None:
                continue
            beta = dict(zip(names, fit["beta"]))
            for label in present:
                drawn[kind][label].append(beta[label])

    records = []
    for kind in matrices:
        for label, mask in masks.items():
            subset = available[mask]
            record = {"drug": drug, "estimate": f"{kind} shift", "group": label,
                      "isolates": len(subset),
                      "clusters": subset.CLUSTER.nunique(),
                      "resistant": int(subset[f"resistant_{drug}"].sum()),
                      "reference_mean": round(point[kind]["intercept"], 3)}
            values = drawn[kind].get(label, [])
            if label in withheld:
                record.update(shift=None, shift_low=None, shift_high=None,
                              withheld=f"{sizes[label]} isolates, "
                                       f"below {MIN_ISOLATES}")
            elif len(values) < draws // 4:
                record.update(shift=round(point[kind][label], 3),
                              shift_low=None, shift_high=None,
                              withheld=f"{len(values)} of {draws} draws fitted")
            else:
                record.update(
                    shift=round(point[kind][label], 3),
                    shift_low=round(float(np.percentile(values, 2.5)), 3),
                    shift_high=round(float(np.percentile(values, 97.5)), 3),
                    withheld="")
            records.append(record)
    return records


def slope_rows(frame, drug, rng):
    """The read-fraction slope, for the minor-allele group and its classes."""
    carriers = frame[frame.MINOR_GROUP.notna() & frame.MINOR_FRS.notna()]
    records = []
    for label in ("all", "loss of function", "substitution", "promoter"):
        subset = carriers if label == "all" else carriers[carriers.MINOR_GROUP.eq(label)]
        subset = subset.dropna(subset=[f"lower_{drug}"])
        record = {"drug": drug, "group": f"minor {label}", "isolates": len(subset),
                  "clusters": subset.CLUSTER.nunique()}
        fit = (fit_slope(subset[f"lower_{drug}"], subset[f"upper_{drug}"],
                         subset.MINOR_FRS) if len(subset) >= MIN_ISOLATES else None)
        if fit is None:
            record.update(slope=None, slope_low=None, slope_high=None,
                          withheld=f"{len(subset)} isolates, below {MIN_ISOLATES}")
        else:
            interval = resample_clusters(
                subset, rng, BOOTSTRAPS,
                lambda sample, d=drug: (
                    lambda f: f["slope"] if f else None)(
                    fit_slope(sample[f"lower_{d}"], sample[f"upper_{d}"],
                              sample.MINOR_FRS)))
            record.update(slope=round(fit["slope"], 3),
                          slope_low=round(interval[0], 3) if interval else None,
                          slope_high=round(interval[1], 3) if interval else None,
                          withheld="")
        records.append(record)
    return records


def main():
    import numpy as np
    import pandas as pd

    say("=" * 72)
    say("The uncertain group, split by call type")
    say("=" * 72)

    mutations = cohort.load_mutations()
    status = cohort.build_status(mutations)
    frame, counts = prepare(mutations, status)

    types = call_types(mutations)
    excluded = frame[frame.GROUP.eq("uncertain")].join(types)
    say(f"\nExcluded as uncertain, with a MIC: {len(excluded):,}")
    say(excluded.TYPE.value_counts().to_string())

    reference = frame[frame.GROUP.eq("reference")]
    for drug in cohort.DRUGS:
        rate = 100 * reference[f"resistant_{drug}"].mean()
        say(f"\n{drug} resistance, reference {rate:.2f}%:")
        for kind, chunk in excluded.groupby("TYPE", observed=True):
            with_mic = chunk[chunk[f"MIC_{drug}"].notna()]
            resistant = int(with_mic[f"resistant_{drug}"].sum())
            say(f"  {kind:18s} {resistant:3d} of {len(with_mic):3d}"
                + (f"  ({100 * resistant / len(with_mic):5.1f}%)" if len(with_mic) else ""))
        for gene in cohort.BDQ_GENES:
            touches = excluded.genes.fillna("").str.contains(gene)
            for kind in ("minor allele only", "null only"):
                chunk = excluded[touches & excluded.TYPE.eq(kind)]
                with_mic = chunk[chunk[f"MIC_{drug}"].notna()]
                if len(with_mic):
                    resistant = int(with_mic[f"resistant_{drug}"].sum())
                    say(f"    {gene:7s} {kind:18s} {resistant:3d} of {len(with_mic):3d}"
                        f"  ({100 * resistant / len(with_mic):5.1f}%)")

    say("\n" + "=" * 72)
    say(f"Resolved {SUBJECT_GENE} minor alleles")
    say("=" * 72)
    say(f"\nEligible isolates: {counts['eligible']:,}. Excluded for carrying more "
        f"than one minor allele at {SUBJECT_GENE}: {counts['multiple']:,}. "
        f"Resolved forms that did not parse: {counts['unparsed']:,}.")
    carriers = frame[frame.MINOR_GROUP.notna()]
    say(f"Of the eligible, {len(carriers):,} carry a MIC.")
    say("")
    say(carriers.groupby("MINOR_GROUP", observed=True).agg(
        isolates=("MINOR_FRS", "size"), clusters=("CLUSTER", "nunique"),
        median_read_fraction=("MINOR_FRS",
                              lambda s: round(float(s.median()), 3))).to_string())

    shifts, slopes, adjusted = [], [], []
    for drug in cohort.DRUGS:
        shifts.extend(shift_rows(frame, drug, stream(f"shift {drug}")))
        slopes.extend(slope_rows(frame, drug, stream(f"slope {drug}")))
        adjusted.extend(adjusted_rows(frame, drug, stream(f"adjusted {drug}")))

    shift_table = pd.DataFrame(shifts)
    say("\nMIC shift against the reference group, in doublings. The major rows use")
    say("the same cohort, filters, estimator and resampling as the minor rows.")
    say("")
    say(shift_table.to_string(index=False))

    slope_table = pd.DataFrame(slopes)
    say("\nDoublings of MIC per unit of read fraction.")
    say("")
    say(slope_table.to_string(index=False))

    adjusted_table = pd.DataFrame(adjusted)
    say("\n" + "=" * 72)
    say(f"Samples carrying more than one minor allele at {SUBJECT_GENE}")
    say("=" * 72)
    say("")
    say(f"Across every genome, {counts['multiple']:,} samples carry more than one "
        f"minor allele at {SUBJECT_GENE}.")
    say(f"{counts['several']:,} of those carry nothing else in the three genes, and "
        f"{int(frame.MULTI_ALLELES.notna().sum()):,} of")
    say("those carry a MIC. Resolved forms that did not parse: "
        f"{counts['several_unparsed']:,}.")
    profile = counts["profile"].reindex(
        frame.index[frame.MULTI_ALLELES.notna()]).dropna(how="all")
    say("")
    say("Per sample, among those carrying a MIC:")
    say("")
    say(profile.describe().round(3).to_string())
    say("")
    say("Alleles per sample: "
        + ", ".join(f"{int(k)} in {v}" for k, v in
                    profile.alleles.value_counts().sort_index().items()) + ".")
    say(f"Samples whose alleles all sit at distinct positions: "
        f"{int((profile.positions == profile.alleles).sum())} of {len(profile)}.")
    say(f"Samples carrying at least one loss of function: "
        f"{int((profile.loss_of_function > 0).sum())}; carrying nothing else: "
        f"{int((profile.loss_of_function == profile.alleles).sum())}.")
    say("")
    say("The read fractions are per-position ratios, so summing them across")
    say("positions is not bounded by one and the sum is not a wild-type")
    say(f"complement: {int((profile.read_fraction_sum > 1.0).sum())} of "
        f"{len(profile)} samples sum above 1.0, to a maximum of "
        f"{profile.read_fraction_sum.max():.3f}.")
    pairs = profile[profile.alleles == 2]
    if len(pairs):
        difference = (pairs.read_fraction_high - pairs.read_fraction_low)
        say(f"Among the {len(pairs)} two-allele samples the two fractions differ by a")
        say(f"median of {difference.median():.3f}, and by more than 0.3 in "
            f"{int((difference > 0.3).sum())} of them, so the")
        say("alleles are mostly not at equal shares.")

    multi = []
    for drug in cohort.DRUGS:
        multi.extend(multi_allele_rows(frame, counts["profile"], drug))
    multi_table = pd.DataFrame(multi)
    say("")
    say("Resistance against the reference group, in percent:")
    say("")
    say(multi_table.to_string(index=False))

    say("\n" + "=" * 72)
    say("The same shifts with site held constant")
    say("=" * 72)
    say("")
    say("Every group refitted jointly with the reference group in one censored")
    say("regression, once plain and once with an indicator per site. The two share")
    say("a reference mean and a standard deviation, so the distance between them is")
    say("the site adjustment. The plain rows differ from the table above in fitting")
    say("one standard deviation across all groups rather than one per group.")
    say("")
    say(adjusted_table.to_string(index=False))
    say("")
    say("Read the intercept as the reference mean at the baseline site, which is the")
    say("first site label in sorted order, rather than across the collection.")

    ESTIMATES.parent.mkdir(parents=True, exist_ok=True)
    # The estimates answer different questions about the same groups, so they
    # are stacked with a column naming which, rather than joined.
    combined = pd.concat([shift_table.assign(estimate="shift"),
                          slope_table.assign(estimate="slope"),
                          adjusted_table], ignore_index=True)
    MULTI.parent.mkdir(parents=True, exist_ok=True)
    multi_table.to_csv(MULTI, index=False)
    columns = ["drug", "estimate", "group", "isolates", "clusters", "resistant",
               "reference_mean", "shift", "shift_low", "shift_high",
               "slope", "slope_low", "slope_high", "withheld"]
    combined.reindex(columns=columns).to_csv(ESTIMATES, index=False)
    REPORT.write_text("\n".join(_lines) + "\n")
    print(f"\nEstimates written to {ESTIMATES}")
    print(f"Multi-allele counts written to {MULTI}")
    print(f"Report written to {REPORT}")


if __name__ == "__main__":
    main()
