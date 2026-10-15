"""
Rv0678 frameshifts with the contributing site held constant.

Run from the project root with the virtual environment active:

    python code/frameshift_site.py

Why this exists. The Rv0678 frameshift class carries the widest spread of any
class in this cohort, a fitted standard deviation of 2.81 doublings on
bedaquiline against 1.09 in the reference group, and one variant accounts for
most of it: 192_ins_g is carried by 51 solo isolates with a MIC, 41 of them at
one site in one sublineage, and its fitted mean falls more than three doublings
below the lowest tested concentration, which is why the atlas withholds its
shift.

Two published explanations for a frameshift that reads as susceptible are in
print. One is biological: a frameshift that conserves the reading sequence
through the open frame or a forward alternative frame leaves a functional
product. The other is a measurement artefact: the isolates behind that
observation came disproportionately from one laboratory group whose MICs sit
systematically low. The second is testable here, because every isolate carries
the site that measured it. Holding the site constant takes the site means out
of the comparison, and what remains of the class effect and of its spread is
what the genotype carries.

What is fitted. Every estimate is one censored linear regression over the
reference group and the carriers together, with an intercept, one indicator per
carrier group, and, in the adjusted fit, one indicator per site beyond the
first. The two fits of a pair run on the same rows under one standard
deviation and differ only in the site indicators, so the distance between them
is the site adjustment alone. Three group structures are fitted:

  class            the Rv0678 frameshift class as one group.
  class less one   the same class with 192_ins_g's carriers removed from the
                   rows rather than folded into the reference group.
  variants         141_ins_c, 138_ins_g, 192_ins_g and the remaining
                   frameshifts, each with a column of its own.

A group resting on fewer than five clusters is left out of the design, and its
rows out of the fit, as the atlas leaves such a variant out of its shifts, and
the table carries the reason in place of the estimate. The draws column counts
the resampled fits an interval rests on: a draw that misses every isolate of a
group is dropped for that group rather than entering its interval as a zero.

A joint fit shares one standard deviation across the reference group and the
carriers, so its shifts differ from the per-class table in outputs/
mic_estimates.csv, which fits each class separately against a separately
fitted reference. The comparison this module supports is between its own plain
and adjusted columns.

Spread is reported separately, because a shared standard deviation cannot
carry it: the class is fitted on its own rows, once with an intercept alone and
once with an intercept and site indicators, and the second standard deviation
is what is left of the spread after the site means are removed.

The table also carries each site's own mean in the reference group, fitted with
an indicator per site. A claim that one laboratory group reads low is a claim
about those means, so they are reported beside the shifts that rest on them.

What a mean of 192_ins_g rests on. 39 of its 51 bedaquiline readings and 40 of
its 51 clofazimine readings sit at or below the lowest tested concentration, so
its mean and its standard deviation are identified jointly by the few readings
above the floor, and the atlas withholds the estimate on that ground. Those
readings are also a single site's: all 39 of the bedaquiline ones, and 39 of
the 40 clofazimine ones, come from the site that carries 41 of the 51 isolates,
where 39 of the 41 sit at the floor, while the 10 isolates at the four other
sites were every one measured on the plate. So the site level that carries the
variant is the level that carries the censoring, and its mean follows the
standard deviation the other levels fix. The censored counts are in the table
beside every estimate, and the findings this module supports are the class's
shift and the class's spread, not a mean for that variant.

Outputs:
  outputs/frameshift_site.csv          one row per drug, estimate and group
  outputs/frameshift_site_report.txt   the readable summary
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import cohort  # noqa: E402
from cluster_adjust import add_clusters  # noqa: E402
from build_atlas import (MIN_CLUSTERS, attach_intervals,  # noqa: E402
                         subject_variants, tested_range)
from heteroresistance import SITE_MINIMUM, cluster_positions, design  # noqa: E402
from mic_model import (concentration_series, fit_censored_linear,  # noqa: E402
                       fit_censored_normal, stream)

TABLE = Path("outputs/frameshift_site.csv")
REPORT = Path("outputs/frameshift_site_report.txt")

SUBJECT_CLASS = "Rv0678 frameshift"
WITHHELD = "192_ins_g"
NAMED = ["141_ins_c", "138_ins_g", "192_ins_g"]
BOOTSTRAPS = 400
# How far a fitted mean may fall outside the tested concentration range, in
# doublings, before the atlas withholds the estimate as extrapolation. The same
# rule is applied here to the adjusted mean, so that the two are comparable.
MAX_EXTRAPOLATION = 1.0
# Isolates a site needs before it gets its own indicator. A joint fit carries
# the reference group's fourteen thousand rows and takes the project's
# threshold; a fit on one class's rows alone carries a hundred and fifty, where
# the same threshold would pool every site into one level and leave nothing for
# the site indicators to remove.
SPREAD_SITE_MINIMUM = 5

_lines = []


def say(text=""):
    print(text)
    _lines.append(text)


def write_report():
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(_lines) + "\n")


def die(message):
    say(f"\nFAILED: {message}")
    write_report()
    sys.exit(1)


def outside_range(mean, limits):
    """The note the atlas prints when a fitted mean leaves the tested range."""
    if not limits:
        return ""
    if mean < limits[0] - MAX_EXTRAPOLATION:
        return (f"fitted mean {mean:.2f} lies more than {MAX_EXTRAPOLATION:g} "
                f"doubling below the lowest tested concentration, {limits[0]:.2f}")
    if mean > limits[1] + MAX_EXTRAPOLATION:
        return (f"fitted mean {mean:.2f} lies more than {MAX_EXTRAPOLATION:g} "
                f"doubling above the highest tested concentration, {limits[1]:.2f}")
    return ""


def pooled_sites(frame, minimum):
    """Site labels for a design, with sites below the minimum pooled.

    The project's rule, with the threshold exposed, because the two fits here
    differ by two orders of magnitude in how many rows they carry.
    """
    sites = frame.SITEID.astype(str)
    counts = sites.value_counts()
    small = set(counts[counts < minimum].index)
    return sites.where(~sites.isin(small), "other")


def build_frame():
    """The atlas frame, with each solo sample's own variant joined to it."""
    mutations = cohort.load_mutations()
    status = cohort.build_status(mutations)
    series, irregular, unreadable = concentration_series()
    # The layout carries a few entries that are not a concentration, across
    # plate designs and drugs this module does not read. mic_model reports them
    # and goes on, and the series it builds is the one every estimate here
    # rests on, so a drug of its own is checked rather than the whole table.
    missing = [drug for drug in cohort.DRUGS
               if not any(key[1] == drug for key in series)]
    if missing:
        die(f"no tested concentrations were read for {', '.join(missing)}")
    frame = attach_intervals(add_clusters(cohort.assemble(status), mutations), series)
    subjects = subject_variants(mutations, status)
    return frame, frame.join(subjects, how="inner"), series


def group_masks(carriers, structure):
    """The carrier columns of one design, as masks over the carrier rows."""
    frameshift = carriers.GROUP.eq(SUBJECT_CLASS)
    if structure == "class":
        return {SUBJECT_CLASS: frameshift}
    if structure == "class less one":
        return {f"{SUBJECT_CLASS} without {WITHHELD}":
                frameshift & carriers.MUTATION.ne(WITHHELD)}
    if structure == "variants":
        masks = {f"Rv0678 {name}": frameshift & carriers.MUTATION.eq(name)
                 for name in NAMED}
        masks[f"{SUBJECT_CLASS}, the other variants"] = (
            frameshift & ~carriers.MUTATION.isin(NAMED))
        return masks
    raise ValueError(f"unknown structure {structure!r}")


def groups_present(sampled, names, labels):
    """The groups one resampled design carries at least one isolate of.

    A draw that misses every isolate of a group leaves that group's
    coefficient unidentified, and the optimizer returns its starting value,
    which is zero. Recording it would pull the group's interval onto zero on
    evidence that is not in the draw, so the draw is dropped for that group
    and kept for the rest.
    """
    return [label for label in labels if sampled[:, names.index(label)].any()]


def counts(rows, drug):
    """What the evidence behind one group is, beside its estimate."""
    return {"isolates": len(rows),
            "clusters": int(rows.CLUSTER.nunique()),
            "sites": int(rows.SITEID.nunique()),
            "largest_cluster": int(rows.CLUSTER.value_counts().max()) if len(rows) else 0,
            "resistant": int(rows[f"resistant_{drug}"].sum()),
            "left_censored": int(rows[f"censored_left_{drug}"].sum()),
            "right_censored": int(rows[f"censored_right_{drug}"].sum())}


def shift_rows(frame, carriers, structure, drug, series):
    """One group's shift against the reference group, plain and with site held
    constant, from one pair of fits over the same rows."""
    import numpy as np
    import pandas as pd

    masks = group_masks(carriers, structure)
    measured = carriers[f"lower_{drug}"].notna()
    counted = {label: carriers[mask & measured] for label, mask in masks.items()}
    # The atlas's own rule, applied to the shifts: a group resting on fewer
    # than MIN_CLUSTERS clusters is left out of the design, and its rows out of
    # the fit, where they would move the shared standard deviation without
    # being estimable themselves. The table carries a row for it either way,
    # with the reason in place of the estimate.
    reasons = {}
    for label, rows in counted.items():
        clusters = rows.CLUSTER.nunique()
        if not len(rows):
            reasons[label] = "no carriers with a measured MIC"
        elif clusters < MIN_CLUSTERS:
            reasons[label] = f"{clusters} clusters, below {MIN_CLUSTERS}"
    labels = [label for label in counted if label not in reasons]

    reference = frame[frame.GROUP.eq("reference")].dropna(subset=[f"lower_{drug}"])
    rows = pd.concat([reference] + [counted[label].reindex(columns=reference.columns)
                                    for label in labels])
    fitted = {label: pd.Series(rows.index.isin(counted[label].index), index=rows.index)
              for label in labels}

    lower = rows[f"lower_{drug}"].to_numpy(dtype=float)
    upper = rows[f"upper_{drug}"].to_numpy(dtype=float)
    sites = pooled_sites(rows, SITE_MINIMUM)
    designs = {"plain shift": design(rows, fitted),
               "site-adjusted shift": design(rows, fitted, sites)}

    point = {}
    for kind, (matrix, names) in designs.items():
        fit = fit_censored_linear(lower, upper, matrix)
        if fit is None:
            die(f"the {kind} fit for {drug} on the {structure} design did not converge")
        point[kind] = dict(zip(names, fit["beta"]), sigma=fit["sigma"])

    groups = cluster_positions(rows.CLUSTER)
    drawn = {kind: {label: [] for label in fitted} for kind in designs}
    for kind, (matrix, names) in designs.items():
        rng = stream(f"frameshift site|{drug}|{structure}|{kind}")
        for _ in range(BOOTSTRAPS):
            picked = rng.integers(0, len(groups), size=len(groups))
            index = np.concatenate([groups[position] for position in picked])
            sampled = matrix[index]
            present = groups_present(sampled, names, fitted)
            fit = fit_censored_linear(lower[index], upper[index], sampled)
            if fit is None:
                continue
            beta = dict(zip(names, fit["beta"]))
            for label in present:
                drawn[kind][label].append(beta[label])

    records = []
    for kind in designs:
        for label in masks:
            subset = counted[label]
            limits = tested_range(subset, series, drug)
            record = {"drug": drug, "estimate": kind, "structure": structure,
                      "group": label, **counts(subset, drug),
                      "reference_mean": round(point[kind]["intercept"], 3),
                      "shift": None, "shift_low": None, "shift_high": None,
                      "draws": None, "fitted_mean": None,
                      "sigma": round(point[kind]["sigma"], 3),
                      "tested_floor": round(limits[0], 3) if limits else None,
                      "tested_ceiling": round(limits[1], 3) if limits else None,
                      "note": reasons.get(label, "")}
            if label in fitted:
                mean = point[kind]["intercept"] + point[kind][label]
                values = drawn[kind][label]
                record.update(shift=round(point[kind][label], 3),
                              fitted_mean=round(mean, 3), draws=len(values),
                              note=outside_range(mean, limits))
                if len(values) >= BOOTSTRAPS // 4:
                    record.update(
                        shift_low=round(float(np.percentile(values, 2.5)), 3),
                        shift_high=round(float(np.percentile(values, 97.5)), 3))
                else:
                    short = f"{len(values)} of {BOOTSTRAPS} draws fitted"
                    record["note"] = (f"{short}; {record['note']}"
                                      if record["note"] else short)
            records.append(record)
    return records


def spread_rows(frame, carriers, drug, series):
    """The class's own standard deviation, before and after the site means."""
    import numpy as np
    import pandas as pd

    records = []
    sets = {SUBJECT_CLASS: carriers[carriers.GROUP.eq(SUBJECT_CLASS)],
            f"{SUBJECT_CLASS} without {WITHHELD}":
                carriers[carriers.GROUP.eq(SUBJECT_CLASS)
                         & carriers.MUTATION.ne(WITHHELD)],
            f"Rv0678 {WITHHELD}":
                carriers[carriers.MUTATION.eq(WITHHELD)],
            "reference": frame[frame.GROUP.eq("reference")]}
    for label, rows in sets.items():
        rows = rows.dropna(subset=[f"lower_{drug}"])
        if not len(rows):
            records.append({
                "drug": drug, "estimate": "plain spread", "structure": "spread",
                "group": label, **counts(rows, drug),
                "reference_mean": None, "shift": None, "shift_low": None,
                "shift_high": None, "draws": None, "fitted_mean": None,
                "sigma": None, "tested_floor": None, "tested_ceiling": None,
                "note": "no carriers with a measured MIC"})
            continue
        lower = rows[f"lower_{drug}"].to_numpy(dtype=float)
        upper = rows[f"upper_{drug}"].to_numpy(dtype=float)
        plain = fit_censored_normal(lower, upper)
        sites = pooled_sites(rows, SPREAD_SITE_MINIMUM)
        levels = sorted(sites.unique())
        matrix = design(rows, {}, sites)[0] if len(levels) > 1 else None
        within = (fit_censored_linear(lower, upper, matrix)
                  if matrix is not None else None)
        limits = tested_range(rows, series, drug)
        fits = []
        if plain is not None:
            fits.append(("plain spread", plain["sigma"], plain["mu"], ""))
        if within is not None:
            # Each row's fitted value under this design is its own site's mean,
            # so the average over rows is the isolate-weighted average of the
            # site means, which is what the plain fit's mean is comparable
            # with. The intercept alone would be one site's mean.
            mean = float((matrix @ within["beta"]).mean())
            fits.append(("within-site spread", within["sigma"], mean,
                         f"{len(levels)} site levels"))
        for kind, sigma, mean, remark in fits:
            note = "; ".join(part for part in (remark, outside_range(mean, limits))
                             if part)
            records.append({
                "drug": drug, "estimate": kind, "structure": "spread",
                "group": label, **counts(rows, drug),
                "reference_mean": None, "shift": None, "shift_low": None,
                "shift_high": None, "draws": None,
                "fitted_mean": round(mean, 3),
                "sigma": round(sigma, 3),
                "tested_floor": round(limits[0], 3) if limits else None,
                "tested_ceiling": round(limits[1], 3) if limits else None,
                "note": note})
    return records


def site_rows(frame, series, drug):
    """Each site's own mean in the reference group.

    The measurement explanation is a claim about a laboratory group reading
    low, so a reader needs to know where the site that carries most of a
    variant's isolates sits when only reference isolates are measured there.
    The reference group alone is fitted, with an indicator per site and the
    same pooling the joint fits use, because carriers would bring their
    genotype into the site mean.
    """
    rows = frame[frame.GROUP.eq("reference")].dropna(subset=[f"lower_{drug}"])
    lower = rows[f"lower_{drug}"].to_numpy(dtype=float)
    upper = rows[f"upper_{drug}"].to_numpy(dtype=float)
    sites = pooled_sites(rows, SITE_MINIMUM)
    matrix, names = design(rows, {}, sites)
    fit = fit_censored_linear(lower, upper, matrix)
    if fit is None:
        die(f"the reference fit by site for {drug} did not converge")
    beta = dict(zip(names, fit["beta"]))
    levels = sorted(sites.unique())

    records = []
    for level in levels:
        part = rows[sites.eq(level)]
        limits = tested_range(part, series, drug)
        shift = beta.get(f"site {level}", 0.0)
        records.append({
            "drug": drug, "estimate": "reference site mean",
            "structure": "site", "group": f"site {level}",
            **counts(part, drug),
            "reference_mean": round(beta["intercept"], 3),
            "shift": round(shift, 3), "shift_low": None, "shift_high": None,
            "draws": None,
            "fitted_mean": round(beta["intercept"] + shift, 3),
            "sigma": round(fit["sigma"], 3),
            "tested_floor": round(limits[0], 3) if limits else None,
            "tested_ceiling": round(limits[1], 3) if limits else None,
            "note": "the baseline level" if level == levels[0] else ""})
    return records


COLUMNS = ["drug", "estimate", "structure", "group", "isolates", "clusters",
           "sites", "largest_cluster", "resistant", "left_censored",
           "right_censored", "reference_mean", "shift",
           "shift_low", "shift_high", "draws", "fitted_mean", "sigma",
           "tested_floor", "tested_ceiling", "note"]


def interval(row):
    if row.shift_low != row.shift_low or row.shift_low is None:
        return "n/a"
    return f"{row.shift_low:.2f} to {row.shift_high:.2f}"


def main():
    import pandas as pd

    say("=" * 72)
    say("Rv0678 frameshifts with the contributing site held constant")
    say("=" * 72)

    frame, carriers, series = build_frame()
    frameshift = carriers[carriers.GROUP.eq(SUBJECT_CLASS)]
    withheld = frameshift[frameshift.MUTATION.eq(WITHHELD)]
    say(f"\nRv0678 frameshift carriers: {len(frameshift):,} isolates in "
        f"{frameshift.CLUSTER.nunique()} clusters at {frameshift.SITEID.nunique()} sites")
    say(f"{WITHHELD}: {len(withheld):,} isolates, "
        f"{int(withheld.SITEID.value_counts().max())} of them at one site")
    say(f"\n{WITHHELD} by site, beside the readings at or below the lowest")
    say("tested concentration:")
    say(withheld.groupby(withheld.SITEID.astype(str)).agg(
        isolates=("SITEID", "size"),
        **{f"{drug} censored": (f"censored_left_{drug}", "sum")
           for drug in cohort.DRUGS}).to_string())
    say(f"\nReference group: {int(frame.GROUP.eq('reference').sum()):,} isolates")

    records = []
    for drug in cohort.DRUGS:
        for structure in ("class", "class less one", "variants"):
            say(f"\nFitting {drug}, the {structure} design, "
                f"{BOOTSTRAPS} draws per interval.")
            records += shift_rows(frame, carriers, structure, drug, series)
        records += spread_rows(frame, carriers, drug, series)
        records += site_rows(frame, series, drug)
    table = pd.DataFrame(records)[COLUMNS]

    for drug in cohort.DRUGS:
        say(f"\n{'=' * 72}\n{drug}\n{'=' * 72}")
        shifts = table[table.drug.eq(drug) & table.estimate.str.endswith("shift")].copy()
        shifts["95%"] = [interval(row) for row in shifts.itertuples()]
        say("\nShifts in doublings against the reference group, one pair of fits")
        say("per design, the pair differing only in the site indicators. The")
        say("censored column counts readings at or below the lowest tested")
        say("concentration.")
        say("")
        say(shifts[["estimate", "group", "isolates", "clusters", "sites",
                    "largest_cluster", "left_censored", "shift", "95%", "sigma"]]
            .rename(columns={"left_censored": "censored"}).to_string(index=False))
        remarks = shifts[shifts.note.ne("") & shifts.note.notna()]
        if len(remarks):
            say("\nWhere an estimate carries a note:")
            for row in remarks.itertuples():
                say(f"  {row.estimate}, {row.group}: {row.note}")
        thin = shifts[shifts.draws.notna() & shifts.draws.lt(BOOTSTRAPS)]
        if len(thin):
            say("\nIntervals resting on fewer than every draw, the rest having")
            say("missed every isolate of the group:")
            for row in thin.itertuples():
                say(f"  {row.estimate}, {row.group}: {int(row.draws)} of {BOOTSTRAPS}")

        places = table[table.drug.eq(drug) & table.structure.eq("site")]
        say("\nEach site's mean in the reference group, from one fit of that")
        say("group with an indicator per site. A claim that a laboratory reads")
        say("low is a claim about these.")
        say("")
        say(places[["group", "isolates", "left_censored", "shift",
                    "fitted_mean", "note"]]
            .rename(columns={"left_censored": "censored"}).to_string(index=False))

        spreads = table[table.drug.eq(drug) & table.estimate.str.endswith("spread")]
        say("\nStandard deviation of each group fitted on its own rows, before")
        say("and after the site means are removed. The mean of a within-site fit")
        say("is the isolate-weighted average of its site means.")
        say("")
        say(spreads[["estimate", "group", "isolates", "sites", "left_censored",
                     "sigma", "fitted_mean", "tested_floor"]]
            .rename(columns={"left_censored": "censored"}).to_string(index=False))
        outside = spreads[spreads.note.str.contains("tested concentration")]
        if len(outside):
            say("\nWhere a fitted mean lies outside the tested range:")
            for row in outside.itertuples():
                say(f"  {row.estimate}, {row.group}: {row.note}")

    TABLE.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(TABLE, index=False)
    write_report()
    print(f"\nTable written to {TABLE} ({len(table)} rows)")
    print(f"Report written to {REPORT}")


if __name__ == "__main__":
    main()
