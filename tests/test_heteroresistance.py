"""
Tests for the analysis of the samples excluded as uncertain.

The point of the module is that a null call and a detected minor allele are
different things, and that a minor allele resolves to a variant that can be
classified like any other. The tests assert both, and assert that the estimates
resample clusters rather than isolates.
"""

import numpy as np
import pandas as pd
import pytest

import cohort
import heteroresistance as hetero


def resolved_mutations(data_dir):
    """The fixture's minor alleles with a resolved form attached, as the real
    table carries one in MINOR_MUTATION."""
    mutations = cohort.load_mutations()
    mutations = mutations.copy()
    resolved = np.where(mutations.MUTATION.eq("C46Z"), "C46G",
                        np.where(mutations.MUTATION.eq("141_minorindel"),
                                 "141_ins_c", None))
    mutations["MINOR_MUTATION"] = resolved
    return mutations


def test_a_null_call_and_a_minor_allele_are_counted_apart(data_dir, expected):
    types = hetero.call_types(cohort.load_mutations())
    counts = types.TYPE.value_counts()
    assert counts["null only"] == expected["n_null_call"]
    assert counts["minor allele only"] == (expected["n_het_call"]
                                           + expected["n_minor_indel"])
    assert "both" not in counts


def test_eligibility_excludes_every_other_finding(data_dir):
    mutations = resolved_mutations(data_dir)
    keep, multiple = hetero.eligible_samples(mutations)
    status = cohort.build_status(cohort.load_mutations())
    assert multiple == 0
    assert keep
    assert not (keep & set(status.index[status.IS_SOLO]))
    assert not (keep & set(status.index[status.IS_REFERENCE]))
    nulls = set(mutations[mutations.IS_NULL_CALL].UNIQUEID)
    assert not (keep & nulls)


def test_a_sample_with_two_minor_alleles_is_excluded(data_dir):
    mutations = resolved_mutations(data_dir)
    first = sorted(hetero.eligible_samples(mutations)[0])[0]
    extra = pd.DataFrame([{
        "UNIQUEID": first, "GENE": "Rv0678", "MUTATION": "T33Z",
        "MINOR_MUTATION": "T33P", "IS_MINOR": True, "IS_NULL": False,
        "IS_HET_CALL": True, "IS_NULL_CALL": False, "REAL_MAJOR": False,
        "UNCERTAIN": True, "FRS": 0.4,
    }])
    doctored = pd.concat([mutations, extra], ignore_index=True)
    keep, multiple = hetero.eligible_samples(doctored)
    assert first not in keep
    assert multiple == 1


def test_a_minor_allele_resolves_to_a_classified_variant(data_dir):
    mutations = resolved_mutations(data_dir)
    keep, _ = hetero.eligible_samples(mutations)
    resolved, unparsed = hetero.resolve(mutations, keep)
    assert unparsed == 0
    assert set(resolved.MUTATION) == {"C46G", "141_ins_c"}
    by_form = resolved.set_index("MUTATION").MINOR_GROUP.to_dict()
    assert by_form["141_ins_c"] == "loss of function"
    assert by_form["C46G"] == "substitution"


def test_a_major_carrier_keeps_its_own_cluster_and_a_minor_carrier_gets_one(data_dir,
                                                                            expected):
    """The clone of 25 shares one cluster through its major variant, and a
    resolved minor allele is keyed on its own form rather than on nothing."""
    mutations = resolved_mutations(data_dir)
    status = cohort.build_status(cohort.load_mutations())
    frame, counts = hetero.prepare(mutations, status)

    clone = frame[frame.GROUP.eq("Rv0678 frameshift")]
    assert len(clone) == expected["n_frameshift"]
    assert clone.CLUSTER.value_counts().max() == expected["n_clone"]

    carriers = frame[frame.MINOR_FORM.notna()]
    assert len(carriers) == counts["eligible"]
    assert carriers.CLUSTER.str.contains("141_ins_c").any()
    assert (carriers.CLUSTER.str.count(r"\|") == 2).all()


def test_the_read_fraction_slope_recovers_a_planted_one():
    ladder = np.log2([0.008, 0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1.0])
    rng = np.random.default_rng(7)
    fraction = rng.uniform(0.05, 0.9, 3000)
    truth = rng.normal(-5.3 + 3.0 * fraction, 1.1)
    lower, upper = [], []
    for value in truth:
        if value <= ladder[0]:
            lower.append(-np.inf)
            upper.append(ladder[0])
        elif value > ladder[-1]:
            lower.append(ladder[-1])
            upper.append(np.inf)
        else:
            index = int(np.searchsorted(ladder, value))
            lower.append(ladder[index - 1])
            upper.append(ladder[index])
    fit = hetero.fit_slope(lower, upper, fraction)
    assert fit["slope"] == pytest.approx(3.0, abs=0.35)
    assert fit["intercept"] == pytest.approx(-5.3, abs=0.25)


def test_the_slope_is_flat_when_the_read_fraction_carries_nothing():
    ladder = np.log2([0.008, 0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1.0])
    rng = np.random.default_rng(11)
    fraction = rng.uniform(0.05, 0.9, 3000)
    truth = rng.normal(-4.0, 1.1, 3000)
    lower, upper = [], []
    for value in truth:
        index = int(np.searchsorted(ladder, value))
        if index == 0:
            lower.append(-np.inf)
            upper.append(ladder[0])
        elif index >= len(ladder):
            lower.append(ladder[-1])
            upper.append(np.inf)
        else:
            lower.append(ladder[index - 1])
            upper.append(ladder[index])
    fit = hetero.fit_slope(lower, upper, fraction)
    assert fit["slope"] == pytest.approx(0.0, abs=0.35)


def test_the_interval_resamples_clusters():
    """Two labellings of one dataset: five clonal groups against forty
    independent isolates."""
    rng = np.random.default_rng(3)
    rows = []
    for cluster in range(5):
        offset = rng.normal(0, 1.2)
        for member in range(8):
            rows.append({"CLUSTER": f"c{cluster}", "value": rng.normal(offset, 0.5)})
    frame = pd.DataFrame(rows)
    statistic = lambda chunk: float(chunk.value.mean())

    clustered = hetero.resample_clusters(frame, np.random.default_rng(1), 300, statistic)
    split = frame.assign(CLUSTER=[str(i) for i in range(len(frame))])
    independent = hetero.resample_clusters(split, np.random.default_rng(1), 300, statistic)
    assert (clustered[1] - clustered[0]) > 1.5 * (independent[1] - independent[0])


def minor_row(unique_id, mutation, resolved, **flags):
    row = {"UNIQUEID": unique_id, "GENE": "Rv0678", "MUTATION": mutation,
           "MINOR_MUTATION": resolved, "IS_MINOR": True, "IS_NULL": False,
           "IS_HET_CALL": False, "IS_NULL_CALL": False, "REAL_MAJOR": False,
           "UNCERTAIN": True, "FRS": 0.4}
    row.update(flags)
    return row


def test_a_het_call_counts_even_where_the_minor_flag_says_otherwise(data_dir):
    """One row in the real table is a het call that IS_MINOR marks false,
    recorded in Section 3.4 of the project record. It is still a detected
    sub-population."""
    mutations = cohort.load_mutations()
    odd = pd.DataFrame([minor_row("odd.sample", "C46Z", "C46G",
                                  IS_MINOR=False, IS_HET_CALL=True)])
    types = hetero.call_types(pd.concat([mutations, odd], ignore_index=True))
    assert types.loc["odd.sample", "TYPE"] == "minor allele only"


def test_a_minor_allele_beside_a_null_call_is_not_eligible(data_dir):
    """A gene that also carries an unreadable position cannot be asserted to
    carry only the minor allele."""
    mutations = resolved_mutations(data_dir)
    sample = "paired.sample"
    rows = pd.DataFrame([
        minor_row(sample, "141_minorindel", "141_ins_c"),
        {"UNIQUEID": sample, "GENE": "pepQ", "MUTATION": "G65X",
         "MINOR_MUTATION": None, "IS_MINOR": False, "IS_NULL": True,
         "IS_HET_CALL": False, "IS_NULL_CALL": True, "REAL_MAJOR": False,
         "UNCERTAIN": True, "FRS": None},
    ])
    doctored = pd.concat([mutations, rows], ignore_index=True)
    keep, _ = hetero.eligible_samples(doctored)
    assert sample not in keep

    without_null = pd.concat([mutations, rows.head(1)], ignore_index=True)
    assert sample in hetero.eligible_samples(without_null)[0]


def confounded_cohort(seed=7, carriers=40, shift=1.0,
                      site_effects=(0.0, 1.5, -0.5), at_high=34):
    """A cohort with a planted group shift and a planted site effect, with the
    carriers concentrated at the site whose mean is highest.

    The group mean then answers a different question from the planted shift,
    because the carriers are drawn from a site that is 1.5 doublings above the
    baseline to begin with. Censoring is the same doubling ladder the plates
    use, so the frame exercises the interval likelihood rather than point MICs.
    """
    rng = np.random.default_rng(seed)
    ladder = np.array([-7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0])
    sites = [f"{index:02d}" for index in range(len(site_effects))]

    def rows(count, site_index, mean):
        values = rng.normal(mean, 1.0, count)
        step = np.searchsorted(ladder, values)
        return pd.DataFrame({
            "SITEID": sites[site_index],
            "lower_BDQ": np.where(step == 0, -np.inf,
                                  ladder[np.clip(step - 1, 0, None)]),
            "upper_BDQ": np.where(step >= len(ladder), np.inf,
                                  ladder[np.clip(step, None, len(ladder) - 1)]),
        })

    reference = pd.concat(
        [rows(1000, index, -5.0 + effect)
         for index, effect in enumerate(site_effects)], ignore_index=True)
    reference["GROUP"] = "reference"
    reference["MINOR_GROUP"] = None

    high = int(np.argmax(site_effects))
    other = (high + 1) % len(site_effects)
    carrier = pd.concat([
        rows(at_high, high, -5.0 + site_effects[high] + shift),
        rows(carriers - at_high, other, -5.0 + site_effects[other] + shift),
    ], ignore_index=True)
    carrier["GROUP"] = "uncertain"
    carrier["MINOR_GROUP"] = "loss of function"

    frame = pd.concat([reference, carrier], ignore_index=True)
    frame["resistant_BDQ"] = frame.lower_BDQ > -3.0
    frame["CLUSTER"] = frame.index.astype(str)
    return frame


def test_the_adjusted_fit_recovers_a_shift_the_group_mean_misses():
    frame = confounded_cohort()
    rng = np.random.default_rng(1)
    rows = pd.DataFrame(hetero.adjusted_rows(frame, "BDQ", rng, draws=40))
    rows = rows[rows.group.eq("minor loss of function")].set_index("estimate")
    assert rows.loc["joint shift", "shift"] > 1.6
    assert abs(rows.loc["site-adjusted shift", "shift"] - 1.0) < 0.35
    assert (rows.loc["site-adjusted shift", "shift_low"]
            < 1.0 < rows.loc["site-adjusted shift", "shift_high"])


def test_the_adjustment_is_what_moves_the_estimate():
    """Spreading the same carriers evenly across the sites leaves the two fits
    together, so the gap between them is the confound and not the model."""
    even = confounded_cohort(at_high=14)
    rng = np.random.default_rng(1)
    rows = pd.DataFrame(hetero.adjusted_rows(even, "BDQ", rng, draws=20))
    rows = rows[rows.group.eq("minor loss of function")].set_index("estimate")
    gap = abs(rows.loc["joint shift", "shift"]
              - rows.loc["site-adjusted shift", "shift"])
    assert gap < 0.2


def test_a_site_below_the_minimum_is_pooled():
    frame = confounded_cohort()
    frame.loc[frame.index[:5], "SITEID"] = "99"
    levels = hetero.site_levels(frame)
    assert "99" not in set(levels)
    assert (levels == "other").sum() == 5


def test_the_design_names_one_column_per_group_and_site_beyond_the_first():
    frame = confounded_cohort()
    masks = {"minor loss of function": frame.MINOR_GROUP.eq("loss of function")}
    plain, names = hetero.design(frame, masks)
    assert names == ["intercept", "minor loss of function"]
    assert plain.shape == (len(frame), 2)

    sites = hetero.site_levels(frame)
    adjusted, adjusted_names = hetero.design(frame, masks, sites)
    assert adjusted_names == names + ["site 01", "site 02"]
    assert adjusted.shape == (len(frame), 4)
    assert (adjusted[:, 0] == 1).all()
    assert adjusted[:, 2:].sum(axis=1).max() == 1


def test_a_group_below_the_minimum_is_withheld_and_left_out_of_the_fit():
    """A group too small to estimate is also too small to leave in the rows,
    where it would be counted as reference and move the reference mean."""
    small = confounded_cohort(carriers=8, at_high=8)
    rows = pd.DataFrame(hetero.adjusted_rows(
        small, "BDQ", np.random.default_rng(1), draws=5))
    group = rows[rows.group.eq("minor loss of function")]
    assert group["shift"].isna().all()
    assert group["withheld"].str.contains("8 isolates, below 12").all()

    none = small[small.MINOR_GROUP.isna()].copy()
    same = pd.DataFrame(hetero.adjusted_rows(
        none, "BDQ", np.random.default_rng(1), draws=5))
    assert (group.reference_mean.to_numpy()
            == same[same.group.eq("minor loss of function")]
            .reference_mean.to_numpy()).all()


def test_a_sample_in_two_groups_is_refused():
    frame = confounded_cohort()
    frame.loc[frame.index[-1], "GROUP"] = "Rv0678 substitution"
    with pytest.raises(ValueError, match="more than one group"):
        hetero.adjusted_rows(frame, "BDQ", np.random.default_rng(1), draws=2)


def test_the_bootstrap_resamples_clusters_not_isolates():
    frame = confounded_cohort()
    frame["CLUSTER"] = "one"
    positions = hetero.cluster_positions(frame.CLUSTER)
    assert len(positions) == 1
    assert len(positions[0]) == len(frame)

    frame["CLUSTER"] = np.where(frame.MINOR_GROUP.notna(), "carriers",
                                frame.index.astype(str))
    positions = hetero.cluster_positions(frame.CLUSTER)
    sizes = sorted(len(block) for block in positions)
    assert sizes[-1] == int(frame.MINOR_GROUP.notna().sum())
    assert sum(sizes) == len(frame)


def test_the_intervals_follow_the_cluster_key():
    """Collapsing the carriers into one clone per site changes both intervals.
    Resampling rows would leave them where the singleton keys put them, because
    the rows themselves are the same.

    It widens the plain interval, where the two clones carry different site
    means and their relative weight moves with the draw, and narrows the
    adjusted one, where the draw can no longer vary the composition inside a
    clone. A draw that misses both clones carries nothing about the group, and
    the lower bound shows whether such a draw was dropped or recorded as zero.
    """
    cloned = confounded_cohort()
    carrier = cloned.MINOR_GROUP.notna()
    cloned["CLUSTER"] = np.where(
        carrier, "clone " + cloned.SITEID.astype(str), cloned.CLUSTER)

    widths, bounds = {}, {}
    for name, frame in (("singleton", confounded_cohort()), ("cloned", cloned)):
        rows = pd.DataFrame(hetero.adjusted_rows(
            frame, "BDQ", np.random.default_rng(1), draws=60))
        rows = rows[rows.group.eq("minor loss of function")].set_index("estimate")
        widths[name] = (rows["shift_high"] - rows["shift_low"]).to_dict()
        bounds[name] = rows["shift_low"].to_dict()

    assert widths["cloned"]["joint shift"] > 1.5 * widths["singleton"]["joint shift"]
    assert (widths["cloned"]["site-adjusted shift"]
            < 0.5 * widths["singleton"]["site-adjusted shift"])
    assert bounds["cloned"]["site-adjusted shift"] > 0.5


def test_each_estimate_draws_from_its_own_stream():
    """An interval computed after another estimate must match the same interval
    computed on its own, or a figure cannot be reproduced without rerunning
    everything that preceded it."""
    frame = confounded_cohort()
    alone = pd.DataFrame(hetero.adjusted_rows(
        frame, "BDQ", hetero.stream("adjusted BDQ"), draws=20))

    spent = hetero.stream("adjusted BDQ")
    spent.integers(0, 100, 5000)
    assert not alone.equals(pd.DataFrame(hetero.adjusted_rows(
        frame, "BDQ", spent, draws=20)))

    after = hetero.stream("shift BDQ")
    after.integers(0, 100, 5000)
    again = pd.DataFrame(hetero.adjusted_rows(
        frame, "BDQ", hetero.stream("adjusted BDQ"), draws=20))
    assert alone.equals(again)
    assert not hetero.stream("adjusted BDQ").integers(0, 1 << 30, 3).tolist() == \
        hetero.stream("adjusted CFZ").integers(0, 1 << 30, 3).tolist()
