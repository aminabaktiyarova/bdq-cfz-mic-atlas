"""
Tests for sample identifier handling.

Identifiers look like site.06.subj.06TB_1276.lab.06MIL2678.iso.2. Two things
have to be read out of them correctly: which patient an isolate came from, so
serial isolates are not counted as independent observations, and whether two
spellings refer to the same sample, since the CRyPTIC release notes record
historic failures where separators were rewritten in one table and not another.
"""

import audit_cohort
import discovery


def test_patient_is_the_site_and_subject():
    assert audit_cohort.patient_of(
        "site.06.subj.06TB_1276.lab.06MIL2678.iso.2") == "06|06TB_1276"


def test_serial_isolates_resolve_to_one_patient():
    first = audit_cohort.patient_of("site.06.subj.06TB_1276.lab.06MIL2678.iso.1")
    second = audit_cohort.patient_of("site.06.subj.06TB_1276.lab.06MIL2678.iso.2")
    assert first == second


def test_subject_fields_are_not_assumed_unique_across_sites():
    left = audit_cohort.patient_of("site.06.subj.X.lab.A.iso.1")
    right = audit_cohort.patient_of("site.10.subj.X.lab.A.iso.1")
    assert left != right


def test_ena_identifiers_parse():
    assert audit_cohort.patient_of(
        "site.ENA.subj.SRR6824540.lab.1.iso.1") == "ENA|SRR6824540"


def test_subject_fields_containing_separators_parse():
    assert audit_cohort.patient_of("site.07.subj.A/B.1.lab.X.Y.iso.1") == "07|A/B.1"


def test_discovery_and_audit_agree_on_the_patient():
    identifier = "site.06.subj.06TB_1276.lab.06MIL2678.iso.2"
    assert discovery.patient_of(identifier) == audit_cohort.patient_of(identifier)


def test_normalising_reconciles_the_two_historic_spellings():
    genetics = "site.07.subj.X.lab.06MIL/2678.1.iso.1"
    phenotypes = "site.07.subj.X.lab.06MIL_2678_1.iso.1"
    assert audit_cohort.normalise(genetics) == audit_cohort.normalise(phenotypes)


def test_normalising_does_not_merge_different_isolates():
    first = "site.06.subj.06TB_1276.lab.X.iso.1"
    second = "site.06.subj.06TB_1276.lab.X.iso.2"
    assert audit_cohort.normalise(first) != audit_cohort.normalise(second)
