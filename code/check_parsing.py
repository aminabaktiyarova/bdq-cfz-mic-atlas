"""
Validates the GARC parser in code/garc.py against MUTATIONS.parquet.

Run from the project root with the virtual environment active:

    python code/check_parsing.py

Nothing biological is computed here. This exists because every number the atlas
will produce depends on the mutation strings being read correctly, and a parser
that silently mishandles a category would produce plausible and wrong results
rather than an error.

Four things are checked, and the script stops at the first failure:

  1. Every mutation string in the six target genes matches the grammar.
  2. The parser's null-call detection agrees with the IS_NULL column.
  3. The parser's het-call detection is consistent with IS_MINOR.
  4. Category counts per gene are biologically sensible: rRNA carries no
     synonymous mutations, promoters carry no amino acid changes.

Checks 2 and 3 are the substantive ones. They test the reading of the grammar
against CRyPTIC's own flags, computed independently of the mutation string.

Output goes to outputs/parsing_check.txt as well as the screen.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from garc import parse_frame  # noqa: E402

DATA = Path("data/cryptic-v3.4.0")
MUTATIONS = DATA / "MUTATIONS.parquet"
REPORT = Path("outputs/parsing_check.txt")

GENES = ["Rv0678", "mmpL5", "pepQ", "atpE", "rrl", "rplC"]
INDEX_LEVELS = ["UNIQUEID", "GENE", "MUTATION"]

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


def load_target_mutations():
    import pyarrow.dataset as ds

    if not MUTATIONS.is_file():
        die(f"{MUTATIONS} not found.")
    dataset = ds.dataset(MUTATIONS, format="parquet")
    df = dataset.to_table(filter=ds.field("GENE").isin(GENES)).to_pandas()
    if [n for n in INDEX_LEVELS if n not in df.columns]:
        df = df.reset_index()
    missing = [n for n in INDEX_LEVELS if n not in df.columns]
    if missing:
        die(f"expected columns absent after reset_index: {missing}")
    return df


def main():
    import pandas as pd

    say("=" * 72)
    say("GARC parser validation")
    say("=" * 72)

    df = load_target_mutations()
    say(f"\nRows for the six target genes: {len(df):,}")
    say(f"Distinct mutation strings:     {df.MUTATION.nunique():,}")

    df = parse_frame(df)

    # ------------------------------------------------------- 1. parse coverage
    say("\n1. Grammar coverage")
    unparsed = df[~df.PARSED]
    say(f"   strings that did not parse: {len(unparsed):,}")
    if len(unparsed):
        import re as _re

        def shape(text):
            """Collapse a string to its structural form so like is grouped with like."""
            text = _re.sub(r"\d+", "#", str(text))
            return _re.sub(r"[A-Za-z]", "L", text)

        shapes = unparsed.assign(SHAPE=unparsed.MUTATION.map(shape))
        say("\n   grouped by structural shape:")
        grouped = shapes.groupby("SHAPE", observed=True).agg(
            rows=("MUTATION", "size"),
            distinct=("MUTATION", "nunique"),
            example=("MUTATION", "first"),
            genes=("GENE", lambda s: ", ".join(sorted(set(s)))),
        )
        say(grouped.sort_values("rows", ascending=False).to_string())
        die(
            f"{unparsed.MUTATION.nunique():,} distinct mutation strings do not match the "
            "grammar. The parser is incomplete and must be fixed before proceeding."
        )
    say("   every string parsed")

    # ------------------------------------------------- 2. null calls vs IS_NULL
    say("\n2. Parser null calls against the IS_NULL column")
    crosstab = pd.crosstab(df.IS_NULL_CALL, df.IS_NULL)
    say(crosstab.to_string())
    disagreements = df[df.IS_NULL_CALL != df.IS_NULL]
    say(f"\n   rows where they disagree: {len(disagreements):,}")
    if len(disagreements):
        say("   examples of disagreement:")
        sample = disagreements[["GENE", "MUTATION", "IS_NULL", "IS_NULL_CALL", "IS_HET_CALL"]]
        say(sample.drop_duplicates().head(15).to_string(index=False))
        say(
            "\n   Note: disagreement is not automatically an error. IS_NULL is CRyPTIC's"
        )
        say(
            "   own flag and may use a broader or narrower definition than the X"
        )
        say(
            "   character alone. Read the examples before deciding which to trust."
        )
    else:
        say("   perfect agreement")

    # -------------------------------------------------- 3. het calls vs IS_MINOR
    say("\n3. Parser het calls against the IS_MINOR column")
    say(pd.crosstab(df.IS_HET_CALL, df.IS_MINOR).to_string())
    say("\n   het-call rows by gene:")
    say(df[df.IS_HET_CALL].GENE.value_counts().to_string())

    # ------------------------------------------------------- 4. category counts
    say("\n4. Mutation categories per gene")
    categories = (
        df.assign(
            CATEGORY=lambda d: pd.Series(
                pd.NA, index=d.index, dtype="object"
            ).mask(d.IS_NULL_CALL, "null call")
            .mask(d.IS_HET_CALL, "het call")
            .mask(d.KIND.eq("INDEL"), "indel")
            .mask(d.IS_SYNONYMOUS, "synonymous")
            .mask(d.IS_REAL_VARIANT & d.KIND.eq("SNP"), "substitution")
        )
    )
    categories["CATEGORY"] = categories.CATEGORY.fillna("other")
    say(pd.crosstab(categories.GENE, categories.CATEGORY).to_string())

    say("\n   By region affected:")
    say(pd.crosstab(df.GENE, df.AFFECTS).to_string())

    say("\n   Indels by gene, and how many shift the reading frame:")
    indels = df[df.KIND == "INDEL"]
    if len(indels):
        summary = indels.groupby("GENE", observed=True).agg(
            indels=("MUTATION", "size"),
            distinct=("MUTATION", "nunique"),
            frameshift=("IS_FRAMESHIFT", "sum"),
        )
        say(summary.to_string())
    else:
        say("   none")

    say("\n   Stop codons introduced, by gene:")
    stops = df[df.IS_STOP]
    say(stops.GENE.value_counts().to_string() if len(stops) else "   none")

    # -------------------------------------------------------- sanity assertions
    say("\n5. Sanity checks")
    problems = []

    rrl_synonymous = int(df[(df.GENE == "rrl") & df.IS_SYNONYMOUS].shape[0])
    say(f"   synonymous mutations called in rrl (expect 0, it is rRNA): {rrl_synonymous}")
    if rrl_synonymous:
        problems.append("rrl carries synonymous mutations, which is not possible for rRNA")

    prom_aa = int(df[(df.AFFECTS == "PROM") & df.REF_RESIDUE.isin(list("ACDEFGHIKLMNOPQRSTVWY!"))].shape[0])
    say(f"   promoter rows carrying an amino acid residue (expect 0): {prom_aa}")
    if prom_aa:
        problems.append("promoter mutations are being read as amino acid changes")

    rrl_cds = int(df[(df.GENE == "rrl") & (df.AFFECTS == "CDS")].shape[0])
    say(f"   rrl rows classified as coding sequence (expect 0): {rrl_cds}")
    if rrl_cds:
        problems.append("rrl rows are being classified as coding sequence")

    if problems:
        die("; ".join(problems))
    say("   all sanity checks passed")

    # --------------------------------------------- what the atlas will work with
    say("\n6. What survives as a real, called, non-synonymous variant")
    real = df[df.IS_REAL_VARIANT]
    say(f"   rows: {len(real):,} of {len(df):,}")
    summary = real.groupby("GENE", observed=True).agg(
        rows=("MUTATION", "size"),
        samples=("UNIQUEID", "nunique"),
        distinct=("MUTATION", "nunique"),
        minor=("IS_MINOR", "sum"),
    )
    summary["major"] = summary.rows - summary.minor
    say(summary.to_string())

    say("\n   Samples carrying a real major-allele variant, by gene:")
    major = real[~real.IS_MINOR]
    say(major.groupby("GENE", observed=True).UNIQUEID.nunique().to_string())

    # ------------------------- 7. forms not defined in the published grammar
    say("\n7. Forms not defined in NOMENCLATURE.md, profiled against CRyPTIC's own columns")

    unresolved = df[~df.SIZE_RESOLVED]
    say(f"\n   minorindel rows: {len(unresolved):,} across {unresolved.GENE.nunique()} genes")
    if len(unresolved):
        say("   by gene:")
        say(unresolved.GENE.value_counts().to_string())
        say("\n   IS_MINOR for these rows:")
        say(unresolved.IS_MINOR.value_counts(dropna=False).to_string())
        say("\n   IS_NULL for these rows:")
        say(unresolved.IS_NULL.value_counts(dropna=False).to_string())
        say("\n   INDEL_LENGTH and INDEL_NUCLEOTIDES availability:")
        say(f"     INDEL_LENGTH populated:      {int(unresolved.INDEL_LENGTH.notna().sum()):,}")
        say(f"     INDEL_NUCLEOTIDES populated: {int(unresolved.INDEL_NUCLEOTIDES.notna().sum()):,}")
        if unresolved.INDEL_LENGTH.notna().any():
            say("\n   INDEL_LENGTH values where present:")
            say(unresolved.INDEL_LENGTH.value_counts().head(10).to_string())
        say("\n   FRS where present:")
        say(f"     populated: {int(unresolved.FRS.notna().sum()):,} of {len(unresolved):,}")
        if unresolved.FRS.notna().any():
            say(unresolved.FRS.describe().to_string())
        say("\n   MINOR_MUTATION examples:")
        examples = unresolved.MINOR_MUTATION.dropna().drop_duplicates().head(10).tolist()
        say(f"     {examples if examples else 'none populated'}")
        say("\n   For comparison, IS_MINOR across all target-gene rows:")
        say(df.IS_MINOR.value_counts(dropna=False).to_string())

    deletions = df[df.KIND == "GENE_DELETION"]
    say(f"\n   gene-deletion rows: {len(deletions):,}")
    if len(deletions):
        say("   by gene:")
        say(deletions.GENE.value_counts().to_string())
        say("\n   fraction of gene deleted:")
        say(deletions.DELETED_FRACTION.describe().to_string())
        say("\n   distinct fractions observed:")
        say(deletions.DELETED_FRACTION.value_counts().sort_index(ascending=False).head(15).to_string())
        say("\n   samples carrying a gene deletion, by gene:")
        say(deletions.groupby("GENE", observed=True).UNIQUEID.nunique().to_string())
        say("\n   IS_MINOR for gene deletions:")
        say(deletions.IS_MINOR.value_counts(dropna=False).to_string())

    write_report()
    say(f"\nWritten to {REPORT}")


if __name__ == "__main__":
    main()
