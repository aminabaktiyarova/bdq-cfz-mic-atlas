"""
Shared fixtures for the test suite.

The tests run against a synthetic dataset built to the same schema as the
CRyPTIC tables, with known ground truth planted in it. Nothing here touches the
real data, so the suite runs in seconds on a machine that has never downloaded
the 1.1 GB MUTATIONS table, and a failure points at the code rather than at the
data having changed.

Ground truth planted in the fixture is recorded in the EXPECTED dictionary and
asserted by the tests. Where a statistical estimator is being checked, the
fixture generates data from a distribution with known parameters and the test
asserts the estimator recovers them.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))

# Concentrations tested on each plate design, matching the real ladders.
LADDERS = {
    ("UKMYC6", "BDQ"): [0.008, 0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1.0],
    ("UKMYC5", "BDQ"): [0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1.0, 2.0],
    ("UKMYC6", "CFZ"): [0.03, 0.06, 0.12, 0.25, 0.5, 1.0, 2.0],
    ("UKMYC5", "CFZ"): [0.06, 0.12, 0.25, 0.5, 1.0, 2.0, 4.0],
}

# The truth the fixture is built to contain.
EXPECTED = {
    "reference_mu": -5.3,
    "reference_sd": 1.1,
    "frameshift_shift": 2.4,
    "substitution_shift": 1.2,
    "n_reference": 600,
    "n_frameshift": 40,
    "n_substitution": 40,
    "n_synonymous_only": 15,   # carry a synonymous variant only, so still reference
    "n_null_call": 12,         # uncallable
    "n_het_call": 9,           # uncallable
    "n_minor_indel": 7,        # uncallable
    "n_pepq": 20,
    "n_clone": 25,             # one clonal outbreak inside the frameshift group
    "n_serial_patients": 10,   # patients contributing a second isolate
}


def _mutation_row(unique_id, gene, mutation, codes_protein=True,
                  is_minor=False, is_null=False):
    return {
        "UNIQUEID": unique_id, "GENE": gene, "MUTATION": mutation,
        "GENE_POSITION": 1.0, "REF": "a", "ALT": "c",
        "NUCLEOTIDE_NUMBER": 1.0, "NUCLEOTIDE_INDEX": 1.0,
        "CODES_PROTEIN": codes_protein, "INDEL_LENGTH": None,
        "INDEL_NUCLEOTIDES": None, "AMINO_ACID_NUMBER": 1.0,
        "AMINO_ACID_SEQUENCE": None, "NUMBER_NUCLEOTIDE_CHANGES": 1.0,
        "IS_NULL": is_null, "IS_MINOR": is_minor, "MINOR_MUTATION": None,
        "MINOR_READS": None, "COVERAGE": 50.0, "FRS": None,
    }


def report_mic(true_log2, ladder):
    """Round a true log2 MIC onto a plate ladder, as a plate reading would."""
    steps = np.log2(np.asarray(ladder))
    for index, step in enumerate(steps):
        if true_log2 <= step:
            return f"<={ladder[0]}" if index == 0 else f"{ladder[index]}"
    return f">{ladder[-1]}"


@pytest.fixture(scope="session")
def dataset(tmp_path_factory):
    """Build the synthetic dataset once and return the directory holding it."""
    rng = np.random.default_rng(1234)
    directory = tmp_path_factory.mktemp("cryptic")

    mutations, samples = [], []
    counter = [0]

    def new_sample(kind, site, lineage, sublineage, shift, patient=None):
        counter[0] += 1
        index = counter[0]
        subject = patient if patient else f"S{index}"
        isolate = 2 if patient else 1
        unique_id = f"site.{site}.subj.{subject}.lab.L{index}.iso.{isolate}"
        samples.append({
            "UNIQUEID": unique_id, "KIND": kind, "SITEID": site,
            "LINEAGE": lineage, "SUBLINEAGE": sublineage, "SHIFT": shift,
            "SUBJECT": subject,
        })
        mutations.append(_mutation_row(unique_id, "mmpL5", "D767N"))
        return unique_id, subject

    sites = ["02", "06", "10"]
    lineages = [("lineage2", "lineage2.2.1"), ("lineage4", "lineage4.10")]

    for i in range(EXPECTED["n_reference"]):
        new_sample("reference", sites[i % 3], *lineages[i % 2], 0.0)

    for i in range(EXPECTED["n_synonymous_only"]):
        unique_id, _ = new_sample("reference", sites[i % 3], *lineages[i % 2], 0.0)
        mutations.append(_mutation_row(unique_id, "Rv0678", "A69A"))

    # Frameshift group: a clonal outbreak plus genuinely diverse isolates.
    for i in range(EXPECTED["n_clone"]):
        unique_id, _ = new_sample("frameshift", "10", "lineage2", "lineage2.2.1",
                                  EXPECTED["frameshift_shift"])
        mutations.append(_mutation_row(unique_id, "Rv0678", "192_ins_g"))
    for i in range(EXPECTED["n_frameshift"] - EXPECTED["n_clone"]):
        unique_id, _ = new_sample("frameshift", sites[i % 3], *lineages[i % 2],
                                  EXPECTED["frameshift_shift"])
        mutations.append(_mutation_row(unique_id, "Rv0678", f"{100 + i}_ins_c"))

    for i in range(EXPECTED["n_substitution"]):
        unique_id, _ = new_sample("substitution", sites[i % 3], *lineages[i % 2],
                                  EXPECTED["substitution_shift"])
        mutations.append(_mutation_row(unique_id, "Rv0678", f"N{4 + i}T"))

    for i in range(EXPECTED["n_pepq"]):
        unique_id, _ = new_sample("pepq", sites[i % 3], *lineages[i % 2], 0.6)
        mutations.append(_mutation_row(unique_id, "pepQ", f"P{60 + i}L"))

    for i in range(EXPECTED["n_null_call"]):
        unique_id, _ = new_sample("uncertain", sites[i % 3], *lineages[i % 2], 0.0)
        mutations.append(_mutation_row(unique_id, "Rv0678", "G65X", is_null=True))
    for i in range(EXPECTED["n_het_call"]):
        unique_id, _ = new_sample("uncertain", sites[i % 3], *lineages[i % 2], 0.0)
        mutations.append(_mutation_row(unique_id, "Rv0678", "C46Z", is_minor=True))
    for i in range(EXPECTED["n_minor_indel"]):
        unique_id, _ = new_sample("uncertain", sites[i % 3], *lineages[i % 2], 0.0)
        mutations.append(_mutation_row(unique_id, "Rv0678", "141_minorindel",
                                       is_minor=True))

    # Serial isolates: a second isolate from an existing patient.
    for i in range(EXPECTED["n_serial_patients"]):
        source = samples[i]
        new_sample("reference", source["SITEID"], source["LINEAGE"],
                   source["SUBLINEAGE"], 0.0, patient=source["SUBJECT"])

    frame = pd.DataFrame(samples)

    # Phenotypes, generated from the planted distributions and censored onto plates.
    phenotypes = []
    for _, row in frame.iterrows():
        design = "UKMYC6" if rng.random() < 0.7 else "UKMYC5"
        for drug in ("BDQ", "CFZ"):
            ladder = LADDERS[(design, drug)]
            true_log2 = rng.normal(
                EXPECTED["reference_mu"] + row.SHIFT, EXPECTED["reference_sd"])
            mic = report_mic(true_log2, ladder)
            numeric = float(mic.replace("<=", "").replace(">", ""))
            phenotypes.append({
                "UNIQUEID": row.UNIQUEID, "DRUG": drug, "PLATEDESIGN": design,
                "BELONGS_GPI": True, "SITEID": row.SITEID, "DILUTION": 1.0,
                "PHENOTYPE_QUALITY": "HIGH" if rng.random() < 0.8 else "LOW",
                "READINGDAY": "14", "PRIMARY_DILUTION": 1.0, "PRIMARY_METHOD": "VZ",
                "AMYGDA_DILUTION": None, "BASHTHEBUG_DILUTION": None,
                "TMAS_DILUTION": None, "PHENOTYPE_DESCRIPTION": "VZ,TM AGREE",
                "BASHTHEBUG_NUMBER_CLASSIFICATIONS": None,
                "MIC": mic, "LOG2MIC": float(np.log2(numeric)),
                "BINARY_PHENOTYPE": "R" if true_log2 > -1.0 else "S",
            })

    pd.DataFrame(mutations).set_index(["UNIQUEID", "GENE", "MUTATION"]).to_parquet(
        directory / "MUTATIONS.parquet")
    pd.DataFrame({
        "UNIQUEID": frame.UNIQUEID, "SPECIES": "M. tuberculosis", "N_LINEAGES": 1,
        "LINEAGE": frame.LINEAGE, "SUBLINEAGE": frame.SUBLINEAGE,
        "MYCOBACTERIAL_READS": 1, "TB_READS": 1, "TB_COVERAGE": 99.0,
        "TB_DEPTH": 90.0, "ANTIBIOGRAM": "S" * 14, "PIPELINE_BUILD": "x",
    }).set_index("UNIQUEID").to_parquet(directory / "GENOMES.parquet")
    pd.DataFrame(phenotypes).set_index(["UNIQUEID", "DRUG"]).to_parquet(
        directory / "UKMYC_PHENOTYPES.parquet")

    # PLATE_LAYOUT stores CONC as text carrying censoring operators, as the real
    # file does despite the schema document describing a float.
    layout = []
    for (design, drug), ladder in LADDERS.items():
        for position, concentration in enumerate(ladder, start=1):
            layout.append({
                "PLATEDESIGN": design, "DRUG": drug, "DILUTION": position,
                "CONC": f"<={concentration}" if position == 1 else str(concentration),
                "ROW": 1, "COL": position, "BINARY_PHENOTYPE": "S",
            })
        layout.append({
            "PLATEDESIGN": design, "DRUG": drug, "DILUTION": len(ladder) + 1,
            "CONC": f">{ladder[-1]}", "ROW": None, "COL": None,
            "BINARY_PHENOTYPE": "R",
        })
    pd.DataFrame(layout).set_index(["PLATEDESIGN", "DRUG", "DILUTION"]).to_parquet(
        directory / "PLATE_LAYOUT.parquet")

    pd.DataFrame({
        "SITEID": sites, "COUNTRY": ["China", "Italy", "South Africa"],
        "DESCRIPTION": ["a", "b", "c"],
    }).to_csv(directory / "SITES.csv.gz", index=False)

    assignment = np.where(np.arange(len(frame)) % 5 == 0, "CRyPTIC-v2.0", "CRyPTIC-v1.0")
    pd.DataFrame({
        "UNIQUEID": frame.UNIQUEID, "run_accession": "R", "study_accession": "P",
        "sample_accession": "S", "center_name": "c", "country": "x",
        "location": "y", "first_public": "z", "fastq_ftp": "f", "fastq_md5": "m",
        "fastq_bytes": "b", "dataset": assignment, "status": "complete",
    }).set_index("UNIQUEID").to_parquet(directory / "WGS_SAMPLES.parquet")

    return directory


@pytest.fixture()
def data_dir(dataset, monkeypatch):
    """Point every module at the synthetic dataset instead of the real one."""
    import cohort

    monkeypatch.setattr(cohort, "DATA", dataset)
    return dataset


@pytest.fixture()
def expected():
    return EXPECTED
