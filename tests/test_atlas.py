"""
Tests for the per-variant evidence table.

The fixture plants one clonal outbreak inside the frameshift group: 25 isolates
sharing one mutation, one site and one sublineage. A catalogue counting isolates
would read that as 25 independent observations. The table must report it as one
cluster.
"""

import pytest

import build_atlas
import cohort
from cluster_adjust import add_clusters


@pytest.fixture()
def table(data_dir):
    mutations = cohort.load_mutations()
    status = cohort.build_status(mutations)
    frame = add_clusters(cohort.assemble(status), mutations)
    subjects = build_atlas.subject_variants(mutations, status)
    return build_atlas.evidence_table(frame, subjects)


def test_the_planted_clone_is_one_cluster(table, expected):
    """25 isolates carrying 192_ins_g at one site in one sublineage."""
    row = table[table.mutation.eq("192_ins_g")].iloc[0]
    assert row.isolates == expected["n_clone"]
    assert row.clusters == 1
    assert row.largest_cluster == expected["n_clone"]


def test_every_planted_variant_appears_once(table, expected):
    assert table.mutation.is_unique
    frameshifts = table[table["class"].eq("frameshift")]
    assert int(frameshifts.isolates.sum()) == expected["n_frameshift"]
    substitutions = table[table.gene.eq("Rv0678") & table["class"].eq("substitution")]
    assert int(substitutions.isolates.sum()) == expected["n_substitution"]
    assert int(table[table.gene.eq("pepQ")].isolates.sum()) == expected["n_pepq"]


def test_the_modifier_gene_is_never_a_subject(data_dir):
    """Every sample in the fixture carries an mmpL5 variant, as almost every
    real sample does. None of them may reach the subject column."""
    mutations = cohort.load_mutations()
    status = cohort.build_status(mutations)
    subjects = build_atlas.subject_variants(mutations, status)
    assert cohort.MODIFIER_GENE not in set(subjects.GENE.astype(str))
    assert set(subjects.GENE.astype(str)) <= set(cohort.BDQ_GENES)


def test_a_sample_carrying_two_variants_is_not_counted(data_dir):
    """Only solo samples enter a variant's row, so the isolate totals cannot
    exceed the solo count."""
    mutations = cohort.load_mutations()
    status = cohort.build_status(mutations)
    frame = add_clusters(cohort.assemble(status), mutations)
    subjects = build_atlas.subject_variants(mutations, status)
    table = build_atlas.evidence_table(frame, subjects)
    solo_with_mic = int(frame.index.isin(status.index[status.IS_SOLO]).sum())
    assert int(table.isolates.sum()) == solo_with_mic


def test_uncertain_samples_are_absent(data_dir, expected):
    """A null, het or minor call means the sample cannot be asserted to carry
    the variant, so it supports no variant's row."""
    mutations = cohort.load_mutations()
    status = cohort.build_status(mutations)
    subjects = build_atlas.subject_variants(mutations, status)
    uncertain = status.index[status.any_uncertain]
    assert len(uncertain) == (expected["n_null_call"] + expected["n_het_call"]
                              + expected["n_minor_indel"])
    assert not subjects.index.isin(uncertain).any()


def test_the_drug_columns_count_only_isolates_with_that_mic(table):
    for drug in cohort.DRUGS:
        assert (table[f"{drug}_isolates"] <= table.isolates).all()
        assert (table[f"{drug}_resistant"] <= table[f"{drug}_isolates"]).all()


def test_the_table_is_ordered_by_independent_evidence(table):
    clusters = list(table.clusters)
    assert clusters == sorted(clusters, reverse=True)


def test_the_guard_fires_if_the_modifier_gene_reaches_the_subject_column(data_dir):
    """The filter excludes mmpL5. The guard is what catches a filter that stops
    excluding it, which is the fault that once gave mmpL5 the repressor's
    effect."""
    import pandas as pd

    mutations = cohort.load_mutations()
    status = cohort.build_status(mutations)
    solo = status.index[status.IS_SOLO][0]
    doctored = pd.concat([
        mutations,
        pd.DataFrame([{
            "UNIQUEID": solo, "GENE": cohort.MODIFIER_GENE, "MUTATION": "D767N",
            "CLASS": "substitution", "REAL_MAJOR": True, "IS_MINOR": False,
        }])], ignore_index=True)
    original = cohort.BDQ_GENES
    try:
        cohort.BDQ_GENES = original + [cohort.MODIFIER_GENE]
        with pytest.raises(ValueError):
            build_atlas.subject_variants(doctored, status)
    finally:
        cohort.BDQ_GENES = original


def planted_frame():
    """Two variants with known evidence. A carries three isolates in two
    clusters, one of them with no bedaquiline MIC. B carries one isolate."""
    import pandas as pd

    rows = [
        ("s1", "A", "10", "lineage2.2.1", "0.25", True),
        ("s2", "A", "10", "lineage2.2.1", "0.5", True),
        ("s3", "A", "06", "lineage4.10", None, False),
        ("s4", "B", "02", "lineage3", "0.12", False),
    ]
    frame = pd.DataFrame([{
        "UNIQUEID": unique_id,
        "CLUSTER": f"{site} | {sublineage} | {mutation}",
        "SITEID": site, "SUBLINEAGE": sublineage, "mmpL5_LOF": False,
        "MIC_BDQ": mic, "resistant_BDQ": resistant,
        "censored_left_BDQ": False, "censored_right_BDQ": False,
        "MIC_CFZ": "0.12", "resistant_CFZ": False,
        "censored_left_CFZ": False, "censored_right_CFZ": False,
    } for unique_id, mutation, site, sublineage, mic, resistant in rows]
    ).set_index("UNIQUEID")
    subjects = pd.DataFrame([
        {"UNIQUEID": unique_id, "GENE": "Rv0678", "MUTATION": mutation,
         "CLASS": "substitution"}
        for unique_id, mutation, _, _, _, _ in rows]).set_index("UNIQUEID")
    return frame, subjects


def test_the_drug_columns_count_only_isolates_carrying_that_mic():
    frame, subjects = planted_frame()
    table = build_atlas.evidence_table(frame, subjects).set_index("mutation")
    assert table.loc["A", "isolates"] == 3
    assert table.loc["A", "BDQ_isolates"] == 2
    assert table.loc["A", "BDQ_resistant"] == 2
    assert table.loc["A", "CFZ_isolates"] == 3
    assert table.loc["B", "isolates"] == 1


def test_a_variant_spanning_two_clusters_is_ordered_above_a_single_one():
    frame, subjects = planted_frame()
    table = build_atlas.evidence_table(frame, subjects)
    assert list(table.mutation) == ["A", "B"]
    assert list(table.clusters) == [2, 1]
    assert list(table.largest_cluster) == [2, 1]


def test_a_sample_carrying_two_variants_never_reaches_a_variant_row(data_dir):
    """The real cohort holds 88 such samples. They are excluded by the solo
    restriction, before the guard against a miscounted subject can fire."""
    import pandas as pd

    mutations = cohort.load_mutations()
    status = cohort.build_status(mutations)
    solo = status.index[status.IS_SOLO][0]
    second = pd.DataFrame([{
        "UNIQUEID": solo, "GENE": "pepQ", "MUTATION": "P60L",
        "CLASS": "substitution", "REAL_MAJOR": True, "IS_MINOR": False,
    }])
    doctored = pd.concat([mutations, second], ignore_index=True)
    status = status.copy()
    status.loc[solo, "IS_SOLO"] = False
    subjects = build_atlas.subject_variants(doctored, status)
    assert solo not in subjects.index


def test_the_module_runs_end_to_end(data_dir, tmp_path, monkeypatch):
    """Every column the report prints must exist on the frame it prints from."""
    monkeypatch.setattr(build_atlas, "TABLE", tmp_path / "atlas_evidence.csv")
    monkeypatch.setattr(build_atlas, "REPORT", tmp_path / "atlas_report.txt")
    build_atlas._lines.clear()
    build_atlas.main()
    assert (tmp_path / "atlas_evidence.csv").is_file()
    report = (tmp_path / "atlas_report.txt").read_text()
    assert "Per-variant evidence" in report
    for column in build_atlas.EVIDENCE_COLUMNS + build_atlas.INFLATION_COLUMNS:
        if column != "isolates_per_cluster":
            assert column in report or column.replace("_", " ") in report


LADDERS = {("UKMYC6", "BDQ"): [0.008, 0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1.0],
           ("UKMYC6", "CFZ"): [0.03, 0.06, 0.12, 0.25, 0.5, 1.0, 2.0]}
REFERENCE_MU = -5.3
REFERENCE_SD = 1.1


def censored_cohort(variant_shift, n_clusters, per_cluster=5, n_reference=800,
                    seed=5, cluster_sd=0.0):
    """A reference group and one variant carrying a planted shift, with MICs
    censored onto a real plate ladder.

    cluster_sd plants a shared offset per cluster, which is what a clone is:
    every isolate in it repeats one event rather than adding an independent
    observation.
    """
    import numpy as np
    import pandas as pd
    from conftest import report_mic

    rng = np.random.default_rng(seed)
    rows = []

    def add(unique_id, group, cluster, true_log2):
        row = {"UNIQUEID": unique_id, "GROUP": group, "CLUSTER": cluster,
               "SITEID": "10", "SUBLINEAGE": "lineage2.2.1", "mmpL5_LOF": False}
        for drug in cohort.DRUGS:
            row[f"PLATEDESIGN_{drug}"] = "UKMYC6"
            row[f"MIC_{drug}"] = report_mic(true_log2, LADDERS[("UKMYC6", drug)])
        rows.append(row)

    for index in range(n_reference):
        add(f"r{index}", "reference", f"r{index}",
            rng.normal(REFERENCE_MU, REFERENCE_SD))
    for cluster in range(n_clusters):
        offset = rng.normal(0.0, cluster_sd) if cluster_sd else 0.0
        for member in range(per_cluster):
            add(f"v{cluster}_{member}", "Rv0678 substitution", f"c{cluster}",
                rng.normal(REFERENCE_MU + variant_shift + offset, REFERENCE_SD))

    frame = pd.DataFrame(rows).set_index("UNIQUEID")
    carriers = [index for index in frame.index if index.startswith("v")]
    subjects = pd.DataFrame([
        {"UNIQUEID": index, "GENE": "Rv0678", "MUTATION": "V1A",
         "CLASS": "substitution"} for index in carriers]).set_index("UNIQUEID")
    table = pd.DataFrame([{"gene": "Rv0678", "mutation": "V1A"}])
    return build_atlas.attach_intervals(frame, LADDERS), subjects, table


def estimate(variant_shift, n_clusters, **kwargs):
    import numpy as np

    frame, subjects, table = censored_cohort(variant_shift, n_clusters, **kwargs)
    return build_atlas.estimate_shifts(
        frame, subjects, table, LADDERS, rng=np.random.default_rng(3)).iloc[0]


def test_a_planted_variant_shift_is_recovered():
    """Eight clusters of five isolates, planted 1.8 doublings above reference."""
    row = estimate(1.8, 8)
    assert row.BDQ_not_estimated == ""
    assert row.BDQ_shift == pytest.approx(1.8, abs=0.4)
    assert row.BDQ_shift_low > 0
    assert row.BDQ_shift_low < row.BDQ_shift < row.BDQ_shift_high


def test_a_variant_below_the_cluster_threshold_carries_no_estimate():
    row = estimate(1.8, build_atlas.MIN_CLUSTERS - 1)
    import math
    assert math.isnan(row.BDQ_shift)
    assert math.isnan(row.BDQ_shift_low)
    assert f"below {build_atlas.MIN_CLUSTERS}" in row.BDQ_not_estimated


def test_a_fitted_mean_far_below_the_tested_range_carries_no_estimate():
    """Every isolate at the plate floor leaves the mean unidentified, and the
    fit answers with a number far below anything measured."""
    import math

    row = estimate(-8.0, 6)
    assert math.isnan(row.BDQ_shift)
    assert "lowest tested concentration" in row.BDQ_not_estimated


def test_the_interval_widens_when_the_same_isolates_share_fewer_clusters():
    """Forty isolates in eight clusters carry more independent evidence than
    forty in five, and the interval must say so."""
    spread = estimate(1.8, 8, per_cluster=5, cluster_sd=0.8)
    clumped = estimate(1.8, 5, per_cluster=8, cluster_sd=0.8)
    assert (clumped.BDQ_shift_high - clumped.BDQ_shift_low) > (
        spread.BDQ_shift_high - spread.BDQ_shift_low)


def test_the_interval_is_resampled_over_clusters_and_not_over_isolates():
    """One dataset, two labellings. Five clonal groups resampled as five
    clusters must give a wider interval than the same isolates resampled as
    forty independent ones."""
    import numpy as np

    frame, subjects, table = censored_cohort(1.8, 5, per_cluster=8, cluster_sd=1.2)
    clustered = build_atlas.estimate_shifts(
        frame, subjects, table, LADDERS, rng=np.random.default_rng(3)).iloc[0]

    split = frame.copy()
    split["CLUSTER"] = list(split.index)
    independent = build_atlas.estimate_shifts(
        split, subjects, table, LADDERS, rng=np.random.default_rng(3)).iloc[0]

    assert clustered.BDQ_shift == pytest.approx(independent.BDQ_shift, abs=1e-9)
    clustered_width = clustered.BDQ_shift_high - clustered.BDQ_shift_low
    independent_width = independent.BDQ_shift_high - independent.BDQ_shift_low
    assert clustered_width > 1.5 * independent_width


def test_the_shift_is_taken_against_the_reference_group_alone():
    """A large, strongly shifted variant group would drag a whole-cohort mean
    with it. The reference mean must not move."""
    row = estimate(4.0, 20, per_cluster=10, n_reference=800)
    assert row.BDQ_reference_mean == pytest.approx(REFERENCE_MU, abs=0.2)
    assert row.BDQ_shift == pytest.approx(4.0, abs=0.4)


def test_the_data_dictionary_documents_every_column_in_order(data_dir, tmp_path,
                                                             monkeypatch):
    """A column added to the table without an entry in docs/DATA_DICTIONARY.md
    would ship undocumented."""
    import csv
    import re
    from pathlib import Path

    written = tmp_path / "atlas_evidence.csv"
    monkeypatch.setattr(build_atlas, "TABLE", written)
    monkeypatch.setattr(build_atlas, "REPORT", tmp_path / "atlas_report.txt")
    build_atlas._lines.clear()
    build_atlas.main()

    columns = next(csv.reader(written.open()))
    dictionary = (Path(__file__).resolve().parents[1]
                  / "docs" / "DATA_DICTIONARY.md").read_text()
    section = dictionary.split("## outputs/atlas_evidence.csv")[1].split("\n---")[0]
    documented = re.findall(r"^\| `([A-Za-z0-9_]+)` \|", section, re.M)
    assert documented, "no columns documented for the atlas table"
    assert columns == documented
