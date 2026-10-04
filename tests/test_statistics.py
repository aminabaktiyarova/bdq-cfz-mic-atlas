"""
Tests for the statistical corrections.

Each estimator is checked against data with a known answer planted in it: a
known odds ratio, a known amount of heterogeneity, a known clonal artefact.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import analyse_groups
import audit_cohort
import cluster_adjust
import discovery


def test_mantel_haenszel_recovers_a_known_odds_ratio():
    """Three sites with very different baselines and one shared genotype effect."""
    rng = np.random.default_rng(7)
    true_odds_ratio = 5.0
    rows = []
    for site, baseline in [("01", 0.02), ("02", 0.20), ("03", 0.08)]:
        for group, count in [("reference", 400), ("exposed", 60)]:
            odds = baseline / (1 - baseline)
            if group == "exposed":
                odds *= true_odds_ratio
            probability = odds / (1 + odds)
            for _ in range(count):
                rows.append({"SITEID": site, "GROUP": group,
                             "resistant_BDQ": rng.random() < probability})
    frame = pd.DataFrame(rows)
    frame["EXPOSED"] = frame.GROUP.eq("exposed")
    estimate = cluster_adjust.mh_odds_ratio(frame, "BDQ")
    assert estimate == pytest.approx(true_odds_ratio, rel=0.4)


def test_homogeneity_test_detects_an_effect_that_differs_by_site():
    rng = np.random.default_rng(3)
    analyse_groups._lines.clear()
    rows = []
    for site, odds_ratio in [("01", 1.0), ("02", 12.0), ("03", 1.0)]:
        for group, count in [("reference", 400), ("exposed", 80)]:
            odds = 0.05 / 0.95
            if group == "exposed":
                odds *= odds_ratio
            probability = odds / (1 + odds)
            for _ in range(count):
                rows.append({"SITEID": site, "GROUP": group,
                             "resistant_BDQ": rng.random() < probability})
    analyse_groups.stratified(pd.DataFrame(rows), "exposed", "reference", "BDQ", "SITEID")
    report = "\n".join(analyse_groups._lines)
    assert "differs across strata" in report


def test_collapsing_clusters_removes_an_effect_driven_by_one_clone():
    """
    The exposed group is one clone of 40, all resistant, plus 12 diverse
    isolates with no real effect. The truth is no effect.
    """
    rng = np.random.default_rng(5)
    rows = []
    for index in range(600):
        rows.append({"CLUSTER": f"ref{index}", "EXPOSED": False, "SITEID": "10",
                     "resistant_BDQ": rng.random() < 0.05})
    for _ in range(40):
        rows.append({"CLUSTER": "clone", "EXPOSED": True, "SITEID": "10",
                     "resistant_BDQ": True})
    for index in range(12):
        rows.append({"CLUSTER": f"solo{index}", "EXPOSED": True, "SITEID": "10",
                     "resistant_BDQ": rng.random() < 0.05})
    frame = pd.DataFrame(rows)

    naive = cluster_adjust.mh_odds_ratio(frame, "BDQ")
    collapsed = cluster_adjust.collapsed_draws(frame, "BDQ", rng)
    assert naive > 20, "the uncorrected estimate should be badly inflated"
    assert collapsed["median"] < 5, "collapsing should remove most of the inflation"


def test_prediction_interval_covers_a_future_estimate():
    """
    The interval must carry both the discovery estimate's uncertainty and the
    held-out half's sampling error. Fixing it at the discovery point estimate
    makes it too narrow and turns sampling variation into recorded failures.
    """
    rng = np.random.default_rng(33)
    ladders = [[0.008, 0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1.0]]
    true_shift, reference_mu, reference_sd, exposed_sd = 2.4, -5.3, 1.1, 1.4
    discovery_n, held_out_n, reference_n = 42, 18, 400

    null = discovery.simulate_shifts(held_out_n, reference_n, reference_mu,
                                     reference_sd, exposed_sd, 0.0, ladders, rng,
                                     draws=60)
    covered = trials = 0
    for _ in range(25):
        estimate = discovery.simulate_shifts(discovery_n, reference_n, reference_mu,
                                             reference_sd, exposed_sd, true_shift,
                                             ladders, rng, draws=1)
        bootstrap = discovery.simulate_shifts(discovery_n, reference_n, reference_mu,
                                              reference_sd, exposed_sd,
                                              float(estimate[0]), ladders, rng, draws=40)
        predicted = discovery.simulate_shifts(held_out_n, reference_n, reference_mu,
                                              reference_sd, exposed_sd, bootstrap,
                                              ladders, rng, draws=60)
        low, high = discovery.predictive_interval(
            predicted, discovery_draws=bootstrap,
            validation_spread=float(null.std(ddof=1)))
        future = discovery.simulate_shifts(held_out_n, reference_n, reference_mu,
                                           reference_sd, exposed_sd, true_shift,
                                           ladders, rng, draws=1)
        trials += 1
        covered += low <= float(future[0]) <= high
    assert covered / trials >= 0.85, \
        f"coverage {covered}/{trials} is well below the nominal 95%"


def test_predictive_interval_is_never_narrower_than_the_combined_error():
    """
    The interval must be at least as wide as the discovery and held-out errors
    combined. The simulated spread can understate that, because the estimator's
    spread depends on where the distribution sits relative to the plate, and a
    narrow interval turns ordinary sampling variation into recorded failures.
    """
    rng = np.random.default_rng(1)
    # simulated estimates deliberately too tightly spread
    simulated = rng.normal(2.4, 0.10, 500)
    discovery_draws = rng.normal(2.4, 0.25, 500)
    validation_spread = 0.34

    low, high = discovery.predictive_interval(simulated)
    naive_width = high - low

    low, high = discovery.predictive_interval(
        simulated, discovery_draws=discovery_draws,
        validation_spread=validation_spread)
    corrected_width = high - low

    combined = float(np.hypot(discovery_draws.std(ddof=1), validation_spread))
    assert corrected_width == pytest.approx(2 * 1.96 * combined, rel=0.02)
    assert corrected_width > naive_width * 2, \
        "the floor on the width was not applied"


def test_predictive_interval_keeps_the_wider_of_the_two():
    """Where the simulation is already wide enough, the floor must not shrink it."""
    rng = np.random.default_rng(2)
    simulated = rng.normal(2.4, 0.80, 500)
    discovery_draws = rng.normal(2.4, 0.10, 500)

    low, high = discovery.predictive_interval(
        simulated, discovery_draws=discovery_draws, validation_spread=0.10)
    assert (high - low) == pytest.approx(2 * 1.96 * simulated.std(ddof=1), rel=0.02)


def test_power_threshold_delivers_the_power_it_claims():
    from scipy.stats import fisher_exact

    rng = np.random.default_rng(11)
    exposed_n, reference_n, rate = 40, 2000, 0.02
    threshold = discovery.detectable_odds_ratio(exposed_n, reference_n, rate, rng)
    assert threshold is not None

    odds = rate / (1 - rate) * threshold
    probability = odds / (1 + odds)
    significant = 0
    for _ in range(300):
        a = int(rng.binomial(exposed_n, probability))
        b = int(rng.binomial(reference_n, rate))
        significant += fisher_exact([[a, exposed_n - a], [b, reference_n - b]])[1] < 0.05
    assert 0.7 <= significant / 300 <= 0.92, "the claimed 80% power is not delivered"


def _resistance_frame(counts):
    """A cohort frame from {group: (isolates, resistant)}."""
    rows = []
    for group, (isolates, resistant) in counts.items():
        for index in range(isolates):
            rows.append({"GROUP": group, "resistant_BDQ": index < resistant})
    return pd.DataFrame(rows)


def test_every_group_but_the_reference_is_compared_with_it():
    """A group that is not a single variant class still has a resistance rate,
    and a table that leaves it out cannot be read against the cohort."""
    frame = _resistance_frame({"reference": (1000, 10),
                               "Rv0678 substitution": (200, 20),
                               "multiple variants": (14, 3),
                               "uncertain": (289, 35)})
    rows = analyse_groups.group_rows(frame, "BDQ")
    assert [row["group"] for row in rows] == ["uncertain", "Rv0678 substitution",
                                              "multiple variants"]
    assert [row["n"] for row in rows] == [289, 200, 14]
    assert rows[2]["resistant"] == 3


def test_a_direct_comparison_recovers_a_planted_odds_ratio():
    """60 of 200 against 20 of 200 is an odds ratio of 3.857."""
    frame = _resistance_frame({"loss of function": (200, 60),
                               "substitution": (200, 20),
                               "reference": (1000, 10)})
    row = analyse_groups.direct(frame, "loss of function", "substitution", "BDQ")
    assert row["n"] == 200 and row["against n"] == 200
    assert row["odds ratio"] == pytest.approx((60 * 180) / (140 * 20), abs=0.01)
    assert float(row["p"]) < 1e-06


def test_restricting_the_strata_asks_a_different_question():
    """Three strata share an effect and a fourth reverses it. Over all four the
    homogeneity test has to fire; over the three it has to stay quiet."""
    rng = np.random.default_rng(11)
    rows = []
    for site, odds_ratio in [("01", 6.0), ("02", 6.0), ("03", 6.0), ("04", 0.15)]:
        for group, count in [("reference", 500), ("exposed", 120)]:
            odds = 0.08 / 0.92
            if group == "exposed":
                odds *= odds_ratio
            probability = odds / (1 + odds)
            for _ in range(count):
                rows.append({"SITEID": site, "GROUP": group,
                             "resistant_BDQ": rng.random() < probability})
    frame = pd.DataFrame(rows)

    analyse_groups._lines.clear()
    analyse_groups.stratified(frame, "exposed", "reference", "BDQ", "SITEID")
    assert "differs across strata" in "\n".join(analyse_groups._lines)

    analyse_groups._lines.clear()
    analyse_groups.stratified(frame, "exposed", "reference", "BDQ", "SITEID",
                              only=("01", "02", "03"))
    restricted = "\n".join(analyse_groups._lines)
    assert "strata restricted to 01, 02, 03" in restricted
    assert "04" not in restricted
    assert "differs across strata" not in restricted


def _clonal_frame(seed=5, references=600, clone=40, solo=12, solo_rate=0.05,
                  clone_rate=1.0):
    """A reference group, one clone of carriers, and diverse carriers.

    Every diverse carrier sits in a cluster of its own, so a collapsed draw
    always takes it and only the clone's representative varies. clone_rate below
    1 makes the clone's members differ, which is what lets the collapsed median
    move from one seed to the next.
    """
    rng = np.random.default_rng(seed)
    rows = [{"CLUSTER": f"ref{index}", "GROUP": "reference", "LINEAGE": "lineage2",
             "SITEID": "10", "resistant_BDQ": rng.random() < 0.05}
            for index in range(references)]
    rows += [{"CLUSTER": "clone", "GROUP": "Rv0678 frameshift", "LINEAGE": "lineage2",
              "SITEID": "10", "resistant_BDQ": rng.random() < clone_rate}
             for _ in range(clone)]
    rows += [{"CLUSTER": f"solo{index}", "GROUP": "Rv0678 frameshift",
              "LINEAGE": "lineage4", "SITEID": "10",
              "resistant_BDQ": rng.random() < solo_rate} for index in range(solo)]
    return pd.DataFrame(rows)


def test_the_comparison_frame_marks_the_carriers_and_keeps_the_reference():
    frame = _clonal_frame()
    built = cluster_adjust.exposed_frame(frame, ["Rv0678 frameshift"])
    assert len(built) == len(frame)
    assert int(built.EXPOSED.sum()) == 52
    assert set(built[~built.EXPOSED].GROUP) == {"reference"}

    one_lineage = cluster_adjust.exposed_frame(frame, ["Rv0678 frameshift"], "lineage4")
    assert set(one_lineage.LINEAGE) == {"lineage4"}
    assert int(one_lineage.EXPOSED.sum()) == 12


def test_the_sweep_reports_the_range_over_its_own_seeds(monkeypatch):
    """The sweep has to recompute under a different stream each time, so its
    range must be the range of the five estimates those streams give."""
    monkeypatch.setattr(cluster_adjust, "DRAWS", 40)
    monkeypatch.setattr(cluster_adjust, "SWEEP_SEEDS", 5)
    frame = cluster_adjust.exposed_frame(
        _clonal_frame(solo=30, solo_rate=0.3, clone_rate=0.5), ["Rv0678 frameshift"])

    swept = cluster_adjust.seed_sweep(frame, "BDQ", "planted")
    alone = [cluster_adjust.collapsed_draws(
        frame, "BDQ", cluster_adjust.stream(f"planted BDQ sweep {index}"))
        for index in range(5)]
    medians = [result["median"] for result in alone]

    assert swept["seeds"] == 5
    assert swept["median low"] == pytest.approx(round(min(medians), 2))
    assert swept["median high"] == pytest.approx(round(max(medians), 2))
    assert swept["median low"] < swept["median high"], \
        "five seeds that all gave one median would mean one stream, not five"


def test_a_committed_pre_registration_is_left_alone(tmp_path, monkeypatch):
    """Rewriting it would put a prediction made after the result was known where
    one made before it used to be."""
    registered = tmp_path / "PRE_REGISTRATION.md"
    registered.write_text("the registered predictions\n")
    monkeypatch.setattr(discovery, "PREREG", registered)
    monkeypatch.setattr(discovery, "prereg_is_registered", lambda: True)
    discovery._lines.clear()

    assert discovery.write_prereg([]) is False
    assert registered.read_text() == "the registered predictions\n"
    assert "is in git, so it is not rewritten" in "\n".join(discovery._lines)


def test_a_pre_registration_not_yet_committed_is_written(tmp_path, monkeypatch):
    target = tmp_path / "PRE_REGISTRATION.md"
    monkeypatch.setattr(discovery, "PREREG", target)
    monkeypatch.setattr(discovery, "prereg_is_registered", lambda: False)
    discovery._lines.clear()

    assert discovery.write_prereg([]) is True
    written = target.read_text()
    assert "Pre-registration" in written
    assert "No comparison in the discovery half reaches an effect size" in written, \
        "with no powered comparison the document has to say so rather than be empty"


def test_a_pre_registration_counts_as_registered_only_once_git_holds_it(tmp_path,
                                                                        monkeypatch):
    """A file on disk is not a registered document. Only git fixes the moment the
    predictions were made, so the guard has to ask git rather than the filesystem."""
    import subprocess

    def run(*arguments):
        subprocess.run(["git", *arguments], capture_output=True, text=True, check=True)

    repository = tmp_path / "repository"
    (repository / "docs").mkdir(parents=True)
    relative = Path("docs/PRE_REGISTRATION.md")
    (repository / relative).write_text("draft predictions\n")
    monkeypatch.chdir(repository)
    monkeypatch.setattr(discovery, "PREREG", relative)

    run("init", "-q")
    assert discovery.prereg_is_registered() is False, \
        "a file git does not track is a draft, whatever it says"

    run("add", str(relative))
    assert discovery.prereg_is_registered() is True, \
        "a staged document is no more ours to rewrite than a committed one"

    run("-c", "user.email=t@example.invalid", "-c", "user.name=t",
        "commit", "-q", "-m", "register the predictions")
    assert discovery.prereg_is_registered() is True


def test_a_missing_mic_is_counted_only_where_the_drug_was_tested():
    """A sample with no row for the drug is not a missing measurement of it, so
    the denominator is the rows that carry a plate design."""
    frame = pd.DataFrame({
        "PLATEDESIGN_BDQ": ["UKMYC6", "UKMYC6", None, None],
        "MIC_BDQ": ["0.25", None, None, "0.5"],
    })
    recorded, missing = audit_cohort.mic_is_missing(frame, "BDQ")
    assert list(recorded) == [True, True, False, False]
    assert list(missing) == [False, True, False, False]
    assert int(recorded.sum()) == 2 and int(missing.sum()) == 1


def test_the_pooled_odds_ratio_recovers_a_planted_ratio():
    """Three strata with different baselines and one shared odds ratio of 4.75."""
    tables = [[[20, 80], [5, 95]], [[40, 60], [12, 88]], [[8, 92], [2, 98]]]
    pooled = audit_cohort.pooled_odds_ratio(tables)
    assert pooled["strata"] == 3
    assert pooled["odds_ratio"] == pytest.approx(4.0, abs=1.0)
    assert pooled["low"] < pooled["odds_ratio"] < pooled["high"]


def test_a_stratum_with_nothing_in_the_numerator_is_left_out():
    """A site that lost no row says nothing about why rows are lost, and a pooled
    estimate over one stratum is that stratum."""
    informative = [[20, 80], [5, 95]]
    empty = [[0, 100], [0, 100]]
    assert audit_cohort.pooled_odds_ratio([informative, empty, informative])["strata"] == 2
    assert audit_cohort.pooled_odds_ratio([informative, empty]) is None
    assert audit_cohort.pooled_odds_ratio([empty, empty]) is None


def _layout_rows(design, pairs):
    """PLATE_LAYOUT rows for one design: (CONC text, S or R)."""
    return [{"PLATEDESIGN": design, "DRUG": "LZD", "CONC": conc,
             "BINARY_PHENOTYPE": label} for conc, label in pairs]


def test_the_lowest_well_counts_and_the_off_scale_bin_does_not():
    """PLATE_LAYOUT writes the lowest well as <=x and adds a >x row for the bin
    above the highest well. The first is a tested concentration and the second is
    not."""
    layout = pd.DataFrame(_layout_rows("UKMYC5", [
        ("<=0.06", "S"), ("0.12", "S"), ("0.25", "S"), ("0.5", "S"), ("1.0", "S"),
        ("2.0", "R"), (">2", "R"),
    ]))
    counts = audit_cohort.resistant_wells(layout, "LZD")["UKMYC5"]
    assert counts["wells"] == 6, "six wells were tested, and >2 is not one of them"
    assert counts["resistant wells"] == 1
    assert counts["highest susceptible"] == 1.0


def test_the_off_scale_bin_does_not_relabel_the_highest_well():
    """The bin above the highest well repeats that well's concentration. Reading
    the bin's label onto the well would move the breakpoint down a dilution."""
    layout = pd.DataFrame(_layout_rows("UKMYC6", [
        ("<=0.06", "S"), ("0.12", "S"), ("1.0", "S"), (">1", "R"),
    ]))
    counts = audit_cohort.resistant_wells(layout, "LZD")["UKMYC6"]
    assert counts["wells"] == 3
    assert counts["resistant wells"] == 0
    assert counts["highest susceptible"] == 1.0


def test_the_censoring_profile_reads_the_operators_and_the_resistant_rows():
    frame = pd.DataFrame({
        "UNIQUEID": ["a", "b", "c", "d", "e"],
        "DRUG": ["DLM"] * 5,
        "MIC": ["<=0.015", "<=0.015", "0.12", ">0.5", ">0.5"],
        "BINARY_PHENOTYPE": ["S", "S", "S", "R", "R"],
    })
    profile = audit_cohort.censoring_profile(frame, "DLM")
    assert profile["MICs"] == 5
    assert profile["left-censored %"] == 40.0
    assert profile["right-censored %"] == 40.0
    assert profile["resistant"] == 2
    assert profile["resistant at the ceiling %"] == 100.0, \
        "the share is of the resistant rows, not of every row"

    restricted = audit_cohort.censoring_profile(frame, "DLM", keep={"a", "b", "c"})
    assert restricted["MICs"] == 3 and restricted["resistant"] == 0
    assert restricted["resistant at the ceiling %"] is None
