"""
Tests for the Rv0678 frameshift fits with the contributing site held constant.

The point of the module is that a class whose carriers are concentrated at one
site cannot be read off a plain fit, so the tests plant a site effect and a
genotype shift separately and assert which fit recovers which. A cohort where
the carriers are spread evenly across the sites carries the same genotype shift
and no confound, and there the two fits must agree; that control is what makes
the gap between them in the confounded cohort the confound rather than the
model.

The withholding rules, the censored counts and the number of draws an interval
rests on are asserted too, because the supplementary table is read by people
who cannot rerun the fits.
"""

import numpy as np
import pandas as pd
import pytest

import build_atlas
import cohort
import frameshift_site as fs
from conftest import LADDERS, report_mic
from heteroresistance import SITE_MINIMUM

DESIGN = "UKMYC6"
DRUG = "BDQ"
SERIES = dict(LADDERS)
# The baseline site's mean and the within-site spread the cohorts are drawn
# with. The baseline sits well inside the bedaquiline ladder, which runs from
# log2(0.008) = -6.97 to log2(1.0) = 0, so a planted site effect of -2 and a
# planted genotype shift of +2 both stay on the plate.
BASELINE = -4.0
SPREAD = 1.0
SITE_EFFECTS = {"01": 0.0, "02": 0.0, "03": -2.0}
ECOFF = np.log2(0.25)


def frame_from(blocks, seed=11):
    """A frame carrying the columns the module reads, from planted truth.

    Each block names a group, a mutation, a site, how many isolates, the mean
    they are drawn around and how many clusters they share. Reference isolates
    are each their own cluster, as the project's cluster key makes them. MICs
    are reported onto a real plate ladder and placed by the pipeline's own
    functions, so the intervals here are the intervals the fits see.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for block in blocks:
        for member in range(block["count"]):
            name = f"s{len(rows)}"
            cluster = (f"{block['site']}|{block['mutation']}|"
                       f"{member % block['clusters']}"
                       if block.get("clusters") else name)
            row = {"UNIQUEID": name, "GROUP": block["group"],
                   "MUTATION": block["mutation"], "SITEID": block["site"],
                   "CLUSTER": cluster}
            value = rng.normal(block["mean"], block.get("sd", SPREAD))
            # One latent value reported on both ladders, so one frame can be
            # read for either drug. Each fit reads one drug's columns.
            for drug in cohort.DRUGS:
                row[f"PLATEDESIGN_{drug}"] = DESIGN
                row[f"MIC_{drug}"] = report_mic(value, LADDERS[(DESIGN, drug)])
            rows.append(row)
    frame = build_atlas.attach_intervals(
        pd.DataFrame(rows).set_index("UNIQUEID"), SERIES)
    for drug in cohort.DRUGS:
        text = frame[f"MIC_{drug}"].astype(str)
        frame[f"censored_left_{drug}"] = text.str.startswith("<=")
        frame[f"censored_right_{drug}"] = text.str.startswith(">")
        # The flag only feeds the counts beside each estimate; no fit reads it.
        frame[f"resistant_{drug}"] = frame[f"lower_{drug}"] >= ECOFF
    return frame


def reference_blocks(per_site=200):
    return [{"group": "reference", "mutation": None, "site": site,
             "count": per_site, "mean": BASELINE + effect, "clusters": None}
            for site, effect in SITE_EFFECTS.items()]


def carrier_block(mutation, site, count, shift, clusters):
    return {"group": fs.SUBJECT_CLASS, "mutation": mutation, "site": site,
            "count": count, "mean": BASELINE + SITE_EFFECTS[site] + shift,
            "clusters": clusters}


def spread_evenly(mutation, total, shift, clusters, sites=("01", "02", "03")):
    return [carrier_block(mutation, site, total // len(sites), shift, clusters)
            for site in sites]


def confounded(withheld_at_low=90, withheld_elsewhere=15, withheld_shift=0.0,
               seed=11):
    """A cohort shaped like the real one: three frameshift variants carrying a
    genotype shift of two doublings, spread across the sites, and a fourth
    whose carriers sit almost entirely at the site that reads two doublings
    low and which carries no shift of its own.
    """
    blocks = (reference_blocks()
              + spread_evenly("141_ins_c", 90, 2.0, 10)
              + spread_evenly("138_ins_g", 60, 2.0, 5)
              + spread_evenly("200_ins_t", 60, 2.0, 7)
              + [carrier_block(fs.WITHHELD, "03", withheld_at_low,
                               withheld_shift, 8),
                 carrier_block(fs.WITHHELD, "01", withheld_elsewhere,
                               withheld_shift, 1)])
    frame = frame_from(blocks, seed=seed)
    return frame, frame[frame.MUTATION.notna()]


def class_spread(frame):
    """The spread records for the class, keyed by estimate."""
    carriers = frame[frame.MUTATION.notna()]
    return {record["estimate"]: record
            for record in fs.spread_rows(frame, carriers, DRUG, SERIES)
            if record["group"] == fs.SUBJECT_CLASS}


def shifts(frame, carriers, structure, draws=60, monkeypatch=None):
    """The module's own records for one structure, keyed by estimate and group."""
    if monkeypatch is not None:
        monkeypatch.setattr(fs, "BOOTSTRAPS", draws)
    else:
        fs.BOOTSTRAPS = draws
    records = fs.shift_rows(frame, carriers, structure, DRUG, SERIES)
    return {(record["estimate"], record["group"]): record for record in records}


@pytest.fixture(autouse=True)
def draws_restored():
    """Every test sets the draw count; none of them leaves it changed."""
    original = fs.BOOTSTRAPS
    yield
    fs.BOOTSTRAPS = original


def test_the_site_adjusted_fit_recovers_the_shift_the_plain_fit_misses():
    """The withheld variant carries no shift of its own and sits at the site
    that reads low, so the plain fit reads it as susceptible and the adjusted
    fit returns it to the reference."""
    frame, carriers = confounded()
    table = shifts(frame, carriers, "variants")
    plain = table[("plain shift", f"Rv0678 {fs.WITHHELD}")]
    adjusted = table[("site-adjusted shift", f"Rv0678 {fs.WITHHELD}")]
    assert plain["shift"] < -0.6
    assert abs(adjusted["shift"]) < 0.45
    # The two intervals do not meet, so the adjustment moves the estimate by
    # more than the resampling carries.
    assert plain["shift_high"] < adjusted["shift_low"]
    for row in (plain, adjusted):
        assert row["shift_low"] < row["shift"] < row["shift_high"]


def test_both_fits_recover_a_shift_carried_evenly_across_the_sites():
    """Carriers spread evenly carry no confound, so the plain fit is right and
    the adjustment must leave it where it is."""
    frame = frame_from(reference_blocks()
                       + spread_evenly("141_ins_c", 300, 2.0, 10))
    carriers = frame[frame.MUTATION.notna()]
    table = shifts(frame, carriers, "variants", draws=25)
    for kind in ("plain shift", "site-adjusted shift"):
        row = table[(kind, "Rv0678 141_ins_c")]
        assert row["shift"] == pytest.approx(2.0, abs=0.3)
        assert row["shift_low"] < row["shift"] < row["shift_high"]


def test_the_confound_and_not_the_model_is_what_moves_the_estimate():
    """Spreading the same carriers evenly across the sites leaves the two fits
    together, so the gap between them in the confounded cohort is the
    confound."""
    frame, carriers = confounded(withheld_at_low=35, withheld_elsewhere=70)
    table = shifts(frame, carriers, "variants")
    gap = abs(table[("plain shift", f"Rv0678 {fs.WITHHELD}")]["shift"]
              - table[("site-adjusted shift", f"Rv0678 {fs.WITHHELD}")]["shift"])
    assert gap < 0.3


def test_the_class_shift_is_pulled_down_by_the_variant_at_the_low_site():
    """The class holds three variants carrying two doublings and one carrying
    none, so neither fit reports two. The adjusted fit reports more of it,
    because the plain fit also carries the low site the fourth sits at."""
    frame, carriers = confounded()
    table = shifts(frame, carriers, "class")
    plain = table[("plain shift", fs.SUBJECT_CLASS)]["shift"]
    adjusted = table[("site-adjusted shift", fs.SUBJECT_CLASS)]["shift"]
    assert 0.0 < plain < 2.0
    assert adjusted > plain + 0.2

    less = shifts(frame, carriers, "class less one")
    label = f"{fs.SUBJECT_CLASS} without {fs.WITHHELD}"
    for kind in ("plain shift", "site-adjusted shift"):
        assert less[(kind, label)]["shift"] == pytest.approx(2.0, abs=0.4)


def test_the_class_less_one_design_drops_the_carriers_from_the_rows():
    """Folding them into the reference instead would move the reference mean,
    which is what every shift is measured against."""
    frame, carriers = confounded()
    label = f"{fs.SUBJECT_CLASS} without {fs.WITHHELD}"
    kept = shifts(frame, carriers, "class less one", draws=5)[
        ("plain shift", label)]

    without = frame[frame.MUTATION.ne(fs.WITHHELD) | frame.MUTATION.isna()]
    same = shifts(without, without[without.MUTATION.notna()],
                  "class less one", draws=5)[("plain shift", label)]
    assert kept["reference_mean"] == same["reference_mean"]
    assert kept["isolates"] == same["isolates"] == 210


def test_the_pair_of_fits_runs_on_the_same_rows():
    """The two estimates differ only in the site indicators, so every count
    beside them, and the rows behind them, are the same."""
    frame, carriers = confounded()
    table = shifts(frame, carriers, "variants", draws=5)
    for group in {group for _, group in table}:
        plain = table[("plain shift", group)]
        adjusted = table[("site-adjusted shift", group)]
        for field in ("isolates", "clusters", "sites", "largest_cluster",
                      "resistant", "left_censored", "right_censored",
                      "tested_floor", "tested_ceiling"):
            assert plain[field] == adjusted[field]
    # The site indicators are what the adjusted fit adds, so its shared
    # standard deviation is the smaller of the two.
    assert (table[("site-adjusted shift", f"Rv0678 {fs.WITHHELD}")]["sigma"]
            < table[("plain shift", f"Rv0678 {fs.WITHHELD}")]["sigma"])


def test_the_variant_columns_partition_the_class():
    frame, carriers = confounded()
    masks = fs.group_masks(carriers, "variants")
    assert list(masks) == [f"Rv0678 {name}" for name in fs.NAMED] + [
        f"{fs.SUBJECT_CLASS}, the other variants"]
    stacked = pd.DataFrame(masks)
    assert stacked.sum(axis=1).max() == 1
    assert (stacked.any(axis=1) == carriers.GROUP.eq(fs.SUBJECT_CLASS)).all()

    one = fs.group_masks(carriers, "class less one")
    label = f"{fs.SUBJECT_CLASS} without {fs.WITHHELD}"
    assert carriers[one[label]].MUTATION.eq(fs.WITHHELD).sum() == 0
    assert one[label].sum() == 210

    with pytest.raises(ValueError, match="unknown structure"):
        fs.group_masks(carriers, "class less two")


def test_a_group_below_the_cluster_minimum_is_withheld_and_left_out():
    """The atlas's own rule. The rows go too, where they would move the shared
    standard deviation without being estimable themselves."""
    blocks = (reference_blocks()
              + spread_evenly("141_ins_c", 90, 2.0, 10)
              + [carrier_block(fs.WITHHELD, "03", 20, 0.0,
                               build_atlas.MIN_CLUSTERS - 1)])
    frame = frame_from(blocks)
    carriers = frame[frame.MUTATION.notna()]
    table = shifts(frame, carriers, "variants", draws=5)

    small = table[("plain shift", f"Rv0678 {fs.WITHHELD}")]
    assert small["note"] == f"4 clusters, below {build_atlas.MIN_CLUSTERS}"
    assert small["shift"] is None and small["shift_low"] is None
    assert small["draws"] is None
    assert small["isolates"] == 20

    without = frame[frame.MUTATION.ne(fs.WITHHELD) | frame.MUTATION.isna()]
    same = shifts(without, without[without.MUTATION.notna()], "variants",
                  draws=5)
    assert (table[("plain shift", "Rv0678 141_ins_c")]["reference_mean"]
            == same[("plain shift", "Rv0678 141_ins_c")]["reference_mean"])


def test_a_named_variant_with_no_carriers_carries_the_reason():
    blocks = reference_blocks() + spread_evenly("141_ins_c", 90, 2.0, 10)
    frame = frame_from(blocks)
    carriers = frame[frame.MUTATION.notna()]
    table = shifts(frame, carriers, "variants", draws=5)
    absent = table[("plain shift", f"Rv0678 {fs.WITHHELD}")]
    assert absent["note"] == "no carriers with a measured MIC"
    assert absent["isolates"] == 0
    assert absent["shift"] is None
    assert absent["draws"] is None
    assert table[("plain shift", "Rv0678 141_ins_c")]["shift"] is not None


def test_a_draw_that_misses_a_group_carries_nothing_about_it():
    """The rule is read off the resampled design: a group whose column came
    back all zeros is not in that draw, and its coefficient there is the
    optimizer's starting value rather than an estimate."""
    sampled = np.array([[1.0, 0.0, 1.0], [1.0, 0.0, 0.0]])
    names = ["intercept", "Rv0678 192_ins_g", "Rv0678 141_ins_c"]
    labels = names[1:]
    assert fs.groups_present(sampled, names, labels) == ["Rv0678 141_ins_c"]
    assert fs.groups_present(np.ones((2, 3)), names, labels) == labels


def test_a_dropped_draw_does_not_enter_the_interval_as_a_zero(monkeypatch):
    """The atlas's cluster minimum keeps a group this thin out of the fits, so
    the minimum is lowered here to reach the case: one clone, which rather
    more than a third of the draws miss."""
    monkeypatch.setattr(fs, "MIN_CLUSTERS", 1)
    blocks = (reference_blocks(per_site=60)
              + [carrier_block("141_ins_c", "01", 20, 2.0, 1)])
    frame = frame_from(blocks)
    carriers = frame[frame.MUTATION.notna()]
    table = shifts(frame, carriers, "class", draws=80, monkeypatch=monkeypatch)
    row = table[("plain shift", fs.SUBJECT_CLASS)]
    assert 20 < row["draws"] < 80
    assert row["shift_low"] > 0.5
    assert row["shift_low"] < row["shift"] < row["shift_high"]


def test_an_interval_resting_on_too_few_fits_is_withheld(monkeypatch):
    """The two point fits run before any draw, so failing every later fit
    leaves the shift standing and takes the interval away."""
    frame, carriers = confounded()
    real = fs.fit_censored_linear
    seen = {"calls": 0}

    def flaky(lower, upper, matrix):
        seen["calls"] += 1
        return real(lower, upper, matrix) if seen["calls"] <= 2 else None

    monkeypatch.setattr(fs, "fit_censored_linear", flaky)
    table = shifts(frame, carriers, "class", draws=12, monkeypatch=monkeypatch)
    row = table[("plain shift", fs.SUBJECT_CLASS)]
    assert row["shift"] is not None
    assert row["shift_low"] is None and row["shift_high"] is None
    assert row["note"] == "0 of 12 draws fitted"


def test_the_estimates_are_reproducible_and_draw_from_their_own_streams():
    frame, carriers = confounded()
    first = shifts(frame, carriers, "class", draws=20)
    second = shifts(frame, carriers, "class", draws=20)
    assert first == second

    keyed = fs.stream(f"frameshift site|{DRUG}|class|plain shift")
    same = fs.stream(f"frameshift site|{DRUG}|class|plain shift")
    other = fs.stream(f"frameshift site|{DRUG}|class|site-adjusted shift")
    draw = keyed.integers(0, 1000, size=5)
    assert (draw == same.integers(0, 1000, size=5)).all()
    assert (draw != other.integers(0, 1000, size=5)).any()


def test_a_site_below_the_minimum_is_pooled_and_the_spread_rule_keeps_it():
    """A joint fit carries the reference group's thousands of rows and takes
    the project's threshold. One class's rows are two orders of magnitude
    fewer: on the real class, 148 isolates across ten sites, that threshold
    leaves two levels and the site indicators almost nothing to remove."""
    sizes = [41, 30, 20, 18, 12, 9, 7, 5, 4, 2]
    assert sum(sizes) == 148
    class_rows = frame_from([
        {"group": fs.SUBJECT_CLASS, "mutation": "141_ins_c",
         "site": f"{index:02d}", "count": size, "mean": BASELINE,
         "clusters": 1}
        for index, size in enumerate(sizes, start=1)])
    assert set(fs.pooled_sites(class_rows, SITE_MINIMUM)) == {"01", "02",
                                                              "other"}
    assert len(set(fs.pooled_sites(class_rows, fs.SPREAD_SITE_MINIMUM))) == 9

    frame, _ = confounded()
    joint = fs.pooled_sites(frame, SITE_MINIMUM)
    assert set(joint) == set(SITE_EFFECTS)
    thin = frame.copy()
    thin.loc[thin.index[:5], "SITEID"] = "99"
    pooled = fs.pooled_sites(thin, SITE_MINIMUM)
    assert "99" not in set(pooled)
    assert (pooled == "other").sum() == 5


def test_the_within_site_fit_removes_a_planted_site_effect():
    """One variant, the same within-site spread at both sites, two doublings
    between them. The plain fit must carry the gap between the site means and
    the within-site fit must not."""
    # Both sites at -3: the planted gap is zero and both fits must report the
    # planted spread.
    together = (reference_blocks(per_site=200)
                + [carrier_block("141_ins_c", "01", 15, 1.0, 5),
                   carrier_block("141_ins_c", "03", 45, 3.0, 9)])
    rows = class_spread(frame_from(together, seed=4))
    assert rows["plain spread"]["sigma"] == pytest.approx(SPREAD, abs=0.3)
    assert rows["within-site spread"]["sigma"] == pytest.approx(SPREAD, abs=0.3)

    # Site 01 at -5 and site 03 at -3, the same spread inside each.
    apart = (reference_blocks(per_site=200)
             + [carrier_block("141_ins_c", "01", 15, -1.0, 5),
                carrier_block("141_ins_c", "03", 45, 3.0, 9)])
    rows = class_spread(frame_from(apart, seed=4))
    assert rows["within-site spread"]["sigma"] == pytest.approx(SPREAD, abs=0.3)
    assert (rows["plain spread"]["sigma"]
            > rows["within-site spread"]["sigma"] + 0.2)
    assert rows["within-site spread"]["note"].startswith("2 site levels")


def test_the_within_site_mean_is_the_isolate_weighted_average():
    """The intercept alone would be the mean of whichever site sorts first,
    which on an unbalanced group is not the group's mean."""
    blocks = (reference_blocks(per_site=200)
              + [carrier_block("141_ins_c", "01", 15, -1.0, 5),
                 carrier_block("141_ins_c", "03", 45, 3.0, 9)])
    rows = class_spread(frame_from(blocks, seed=4))
    # Site 01 carries 15 isolates around -5 and site 03 carries 45 around -3,
    # so the weighted average sits at -3.5 and the baseline level at -5.
    assert rows["within-site spread"]["fitted_mean"] == pytest.approx(
        -3.5, abs=0.35)
    assert rows["within-site spread"]["fitted_mean"] > -4.2


def test_each_site_mean_is_fitted_on_the_reference_group_alone():
    """A claim that a laboratory reads low is a claim about these means, so the
    carriers stay out of them: their genotype would raise the site they sit
    at."""
    frame, _ = confounded()
    rows = {record["group"]: record
            for record in fs.site_rows(frame, SERIES, DRUG)}
    assert set(rows) == {f"site {site}" for site in SITE_EFFECTS}
    assert rows["site 01"]["note"] == "the baseline level"
    assert rows["site 01"]["shift"] == 0.0
    for site, effect in SITE_EFFECTS.items():
        row = rows[f"site {site}"]
        assert row["isolates"] == 200
        assert row["shift"] == pytest.approx(effect, abs=0.25)
        assert row["fitted_mean"] == pytest.approx(BASELINE + effect, abs=0.25)
        assert set(row) == set(fs.COLUMNS)


def test_the_extrapolation_note_follows_the_atlas_rule():
    """The wording is the atlas's own, because the two tables are read side by
    side, and the boundary is the same: more than one doubling outside."""
    limits = (-6.966, 0.0)
    assert fs.outside_range(-7.9, limits) == ""
    assert fs.outside_range(-6.966 - fs.MAX_EXTRAPOLATION, limits) == ""
    assert fs.outside_range(-10.082, limits) == (
        "fitted mean -10.08 lies more than 1 doubling below the lowest "
        "tested concentration, -6.97")
    assert fs.outside_range(1.5, limits) == (
        "fitted mean 1.50 lies more than 1 doubling above the highest "
        "tested concentration, 0.00")
    assert fs.outside_range(-99.0, None) == ""


def test_every_record_carries_the_table_columns_and_the_counts():
    frame, carriers = confounded()
    records = fs.shift_rows(frame, carriers, "variants", DRUG, SERIES)
    records += fs.spread_rows(frame, carriers, DRUG, SERIES)
    for record in records:
        assert set(record) == set(fs.COLUMNS)
    table = pd.DataFrame(records)[fs.COLUMNS]
    assert len(table) == 8 + 8

    withheld = carriers[carriers.MUTATION.eq(fs.WITHHELD)]
    row = next(record for record in records
               if record["group"] == f"Rv0678 {fs.WITHHELD}")
    assert row["isolates"] == len(withheld) == 105
    assert row["sites"] == 2
    assert row["left_censored"] == int(
        withheld[f"censored_left_{DRUG}"].sum())
    assert row["resistant"] == int(withheld[f"resistant_{DRUG}"].sum())


def test_the_module_runs_end_to_end(data_dir, tmp_path, monkeypatch):
    """Every block of the report is printed for each drug, and the table it
    writes carries the documented columns and every kind of estimate."""
    monkeypatch.setattr(fs, "BOOTSTRAPS", 6)
    monkeypatch.setattr(fs, "TABLE", tmp_path / "frameshift_site.csv")
    monkeypatch.setattr(fs, "REPORT", tmp_path / "frameshift_site_report.txt")
    monkeypatch.setattr(fs, "_lines", [])
    fs.main()

    table = pd.read_csv(fs.TABLE)
    assert list(table.columns) == fs.COLUMNS
    assert set(table.drug) == set(cohort.DRUGS)
    assert set(table.estimate) == {"plain shift", "site-adjusted shift",
                                   "plain spread", "within-site spread",
                                   "reference site mean"}
    assert set(table.structure) == {"class", "class less one", "variants",
                                    "spread", "site"}
    report = fs.REPORT.read_text()
    for heading in ("Shifts in doublings", "Each site's mean",
                    "Standard deviation of each group"):
        assert report.count(heading) == len(cohort.DRUGS), heading


def test_the_module_reads_the_cohort_and_splits_the_planted_clone(data_dir,
                                                                  expected):
    """The fixture cohort carries the real one's shape: a clone of one
    frameshift variant at a single site. It plants no site effect, so the two
    fits must agree, and the clone is one cluster, which is below the atlas's
    minimum."""
    fs.BOOTSTRAPS = 10
    frame, carriers, series = fs.build_frame()
    class_rows = carriers[carriers.GROUP.eq(fs.SUBJECT_CLASS)]
    assert len(class_rows) == expected["n_frameshift"]
    clone = class_rows[class_rows.MUTATION.eq(fs.WITHHELD)]
    assert len(clone) == expected["n_clone"]
    assert clone.CLUSTER.nunique() == 1
    assert set(clone.SITEID) == {"10"}

    table = {(record["estimate"], record["group"]): record
             for record in fs.shift_rows(frame, carriers, "class", DRUG,
                                         series)}
    for kind in ("plain shift", "site-adjusted shift"):
        assert table[(kind, fs.SUBJECT_CLASS)]["shift"] == pytest.approx(
            expected["frameshift_shift"], abs=0.5)

    variants = {(record["estimate"], record["group"]): record
                for record in fs.shift_rows(frame, carriers, "variants", DRUG,
                                            series)}
    assert variants[("plain shift", f"Rv0678 {fs.WITHHELD}")]["note"] == (
        f"1 clusters, below {build_atlas.MIN_CLUSTERS}")
