"""
Builds the per-sample genotype status table and reports the group summary.

Run from the project root with the virtual environment active:

    python code/build_sample_status.py

All cohort definitions live in code/cohort.py so that there is one place where
"solo", "reference", "uncertain" and the variant classes are decided. This
script assembles the table, writes it, and prints a summary. The stratified
analysis is in code/analyse_groups.py.

Outputs:
  outputs/sample_status.parquet     one row per sample with a MIC
  outputs/sample_status_report.txt  the readable summary
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import cohort  # noqa: E402

OUT_TABLE = Path("outputs/sample_status.parquet")
REPORT = Path("outputs/sample_status_report.txt")

_lines = []


def say(text=""):
    print(text)
    _lines.append(text)


def main():
    import pandas as pd

    say("=" * 72)
    say("Per-sample genotype status")
    say("=" * 72)

    mutations = cohort.load_mutations()
    say(f"\nTarget-gene mutation rows: {len(mutations):,}")

    say("\nReal major-allele variants by gene and class:")
    say(pd.crosstab(
        mutations[mutations.REAL_MAJOR].GENE,
        mutations[mutations.REAL_MAJOR].CLASS,
    ).to_string())

    status = cohort.build_status(mutations)

    say("\nSample classification across all genomes:")
    say(f"  total samples:                             {len(status):,}")
    say(f"  reference (no variant, nothing uncertain): {int(status.IS_REFERENCE.sum()):,}")
    say(f"  solo (exactly one variant):                {int(status.IS_SOLO.sum()):,}")
    say(f"  two or more variants:                      "
        f"{int((~status.IS_SOLO & ~status.IS_REFERENCE & ~status.any_uncertain).sum()):,}")
    say(f"  excluded as uncertain:                     {int(status.any_uncertain.sum()):,}")

    say("\n  Effect of the uncertainty exclusion on the reference group:")
    say(f"    would be reference ignoring uncertainty: "
        f"{int((status.n_bdq_variants == 0).sum()):,}")
    say(f"    actually reference:                      {int(status.IS_REFERENCE.sum()):,}")
    say(f"    excluded, a gene carries a null, het or minor call: "
        f"{int(((status.n_bdq_variants == 0) & status.any_uncertain).sum()):,}")

    say("\n  What makes a sample uncertain, by gene:")
    for gene in cohort.BDQ_GENES:
        say(f"    {gene:8s} {int(status[f'uncertain_{gene}'].sum()):,}")

    say(f"\n  mmpL5 loss of function, the covariate: {int(status.mmpL5_LOF.sum()):,} samples")

    say("\nGroup sizes across all genomes:")
    counts = status.GROUP.value_counts()
    for group in cohort.GROUP_ORDER:
        if group in counts:
            say(f"  {group:24s} {counts[group]:,}")

    joined = cohort.assemble(status)
    say(f"\nSamples with genotype status and a UKMYC MIC: {len(joined):,}")

    say("\nGroup sizes among those with a MIC:")
    counts = joined.GROUP.value_counts()
    for group in cohort.GROUP_ORDER:
        if group in counts:
            say(f"  {group:24s} {counts[group]:,}")

    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    joined.to_parquet(OUT_TABLE)
    say(f"\nPer-sample table written to {OUT_TABLE} ({len(joined):,} rows)")
    REPORT.write_text("\n".join(_lines) + "\n")
    say(f"Report written to {REPORT}")


if __name__ == "__main__":
    main()
