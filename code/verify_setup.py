"""
Verifies that the analysis environment is correctly installed and
that the local copy of the CRyPTIC v3.4.0 tables reproduces known figures.

Run from the project root with the virtual environment active:

    python code/verify_setup.py

It checks three things in order and stops at the first failure:

  1. Every required package imports and reports its pinned version.
  2. The expected data files are present where the project expects them.
  3. The four target drugs produce the sample and resistance counts recorded in
     docs/PROVENANCE.md.

Step 3 is the real test. If the counts match, the local data is the same release
that the recorded findings were computed from, and any later number computed
locally can be trusted to sit on the same foundation.
"""

import sys
from pathlib import Path

EXPECTED_VERSIONS = {
    "pandas": "3.0.5",
    "numpy": "2.5.2",
    "pyarrow": "25.0.1",
    "scipy": "1.18.1",
    "statsmodels": "0.15.0",
    "matplotlib": "3.11.1",
}

# CRyPTIC v3.4.0, Zenodo version DOI 10.5281/zenodo.15680920.
# Counts are for samples holding both a UKMYC MIC and a genome.
EXPECTED_MATCHED_SAMPLES = 15160
EXPECTED_COUNTS = {
    # drug: (rows, resistant, resistant at HIGH phenotype quality)
    "BDQ": (15156, 171, 141),
    "CFZ": (15158, 683, 426),
    "DLM": (15158, 229, 152),
    "LZD": (15156, 217, 160),
}

DATA_DIR = Path("data/cryptic-v3.4.0")
REQUIRED_FILES = ["UKMYC_PHENOTYPES.parquet", "GENOMES.parquet"]


def fail(message):
    print(f"\nFAILED: {message}")
    sys.exit(1)


def check_packages():
    print("Step 1: packages")
    print(f"  python {sys.version.split()[0]}")
    if sys.version_info[:2] < (3, 14):
        print(f"  note: expected Python 3.14, found {sys.version_info.major}.{sys.version_info.minor}")
    import importlib

    for name, expected in EXPECTED_VERSIONS.items():
        try:
            module = importlib.import_module(name)
        except ImportError:
            fail(f"{name} is not installed. Run: pip install -r requirements.txt")
        found = getattr(module, "__version__", "unknown")
        flag = "" if found == expected else f"   (expected {expected})"
        print(f"  {name:12s} {found}{flag}")
    print("  all packages import\n")


def check_files():
    print("Step 2: data files")
    if not DATA_DIR.is_dir():
        fail(f"{DATA_DIR} does not exist. Create it and put the CRyPTIC tables inside.")
    for name in REQUIRED_FILES:
        path = DATA_DIR / name
        if not path.is_file():
            fail(f"{path} is missing.")
        size_mb = path.stat().st_size / 1_000_000
        print(f"  {name:28s} {size_mb:8.1f} MB")
    print("  required files present\n")


def check_counts():
    import pandas as pd

    print("Step 3: reproducing recorded counts")
    phenotypes = pd.read_parquet(DATA_DIR / "UKMYC_PHENOTYPES.parquet").reset_index()
    genomes = pd.read_parquet(DATA_DIR / "GENOMES.parquet").reset_index()

    genome_ids = set(genomes.UNIQUEID)
    matched = phenotypes[phenotypes.UNIQUEID.isin(genome_ids)]
    n_matched = matched.UNIQUEID.nunique()

    print(f"  samples with a UKMYC MIC:            {phenotypes.UNIQUEID.nunique():,}")
    print(f"  samples in GENOMES:                  {len(genome_ids):,}")
    print(f"  samples with both:                   {n_matched:,}")
    if n_matched != EXPECTED_MATCHED_SAMPLES:
        fail(
            f"expected {EXPECTED_MATCHED_SAMPLES:,} matched samples, found {n_matched:,}. "
            "The local data is probably a different CRyPTIC version."
        )

    print(f"\n  {'drug':6s} {'rows':>8s} {'resistant':>10s} {'at HIGH':>9s}   result")
    problems = []
    for drug, (exp_rows, exp_r, exp_r_high) in EXPECTED_COUNTS.items():
        subset = matched[matched.DRUG == drug]
        rows = len(subset)
        resistant = int((subset.BINARY_PHENOTYPE == "R").sum())
        high = subset[subset.PHENOTYPE_QUALITY == "HIGH"]
        resistant_high = int((high.BINARY_PHENOTYPE == "R").sum())
        ok = (rows, resistant, resistant_high) == (exp_rows, exp_r, exp_r_high)
        if not ok:
            problems.append(
                f"{drug}: got ({rows}, {resistant}, {resistant_high}), "
                f"expected ({exp_rows}, {exp_r}, {exp_r_high})"
            )
        print(f"  {drug:6s} {rows:8,d} {resistant:10,d} {resistant_high:9,d}   {'match' if ok else 'MISMATCH'}")

    if problems:
        fail("counts do not match:\n    " + "\n    ".join(problems))
    print("\n  all counts match the recorded figures\n")


if __name__ == "__main__":
    check_packages()
    check_files()
    check_counts()
    print("Setup verified. The environment and the local data both check out.")
