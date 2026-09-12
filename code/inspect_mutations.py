"""
Inspects MUTATIONS.parquet without loading it into memory, and reports what the
atlas needs to know before any ingestion code is written.

Run from the project root with the virtual environment active:

    python code/inspect_mutations.py

MUTATIONS.parquet is around 1.1 GB on disk and 1.3 GB uncompressed, so nothing
here reads the whole file at once. The gene vocabulary is counted by streaming
one column in batches, and only rows belonging to the six target genes are ever
materialised.

Note on structure: UNIQUEID, GENE and MUTATION are stored as pandas index levels
rather than ordinary columns. Arrow-level filters address them as fields, but
anything read into pandas needs the index reset before those names work as
columns. Every read in this project has to do that.

The gene vocabulary is cached to outputs/gene_vocabulary.csv on the first run and
reused afterwards, so reruns skip the streaming pass. Delete that file to force a
fresh count.

Output is printed and written to outputs/mutations_inspection.txt.
"""

import sys
from collections import Counter
from pathlib import Path

DATA = Path("data/cryptic-v3.4.0")
MUTATIONS = DATA / "MUTATIONS.parquet"
REPORT = Path("outputs/mutations_inspection.txt")
VOCAB_CACHE = Path("outputs/gene_vocabulary.csv")

# Target loci. Matched case-insensitively against the file's own vocabulary,
# because the exact spelling CRyPTIC uses is not assumed.
TARGETS = {
    "rv0678": "bedaquiline and clofazimine, efflux pump repressor",
    "mmpl5": "bedaquiline and clofazimine, efflux pump, epistatic modifier",
    "pepq": "bedaquiline and clofazimine, low-level resistance",
    "atpe": "bedaquiline, drug target",
    "rrl": "linezolid, 23S rRNA",
    "rplc": "linezolid, ribosomal protein L3",
}

# Names stored as index levels rather than columns in the CRyPTIC parquet tables.
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


def flatten(df):
    """Move any index levels back into ordinary columns."""
    if [n for n in INDEX_LEVELS if n not in df.columns]:
        df = df.reset_index()
    still_missing = [n for n in INDEX_LEVELS if n not in df.columns]
    if still_missing:
        die(f"expected columns absent after reset_index: {still_missing}")
    return df


def gene_vocabulary(pf, pa, pc):
    """Count rows per gene, streaming one column in batches. Cached."""
    import pandas as pd

    if VOCAB_CACHE.is_file():
        cached = pd.read_csv(VOCAB_CACHE)
        say(f"Gene vocabulary read from cache: {VOCAB_CACHE}")
        say(f"  distinct genes: {len(cached):,}")
        return Counter(dict(zip(cached.GENE, cached.ROWS)))

    say("Streaming the GENE column to build the vocabulary.")
    say("This reads one column in batches and holds only counts in memory.")
    counts = Counter()
    rows_seen = 0
    for batch in pf.iter_batches(batch_size=2_000_000, columns=["GENE"]):
        column = batch.column(0)
        if pa.types.is_dictionary(column.type):
            column = column.dictionary_decode()
        for entry in pc.value_counts(column).to_pylist():
            counts[entry["values"]] += entry["counts"]
        rows_seen += batch.num_rows
        print(f"  ...{rows_seen:,} rows", end="\r", flush=True)
    print(" " * 40, end="\r")
    say(f"  streamed {rows_seen:,} rows")
    say(f"  distinct genes in file: {len(counts):,}")

    VOCAB_CACHE.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(sorted(counts.items()), columns=["GENE", "ROWS"]).to_csv(VOCAB_CACHE, index=False)
    say(f"  cached to {VOCAB_CACHE}")
    return counts


def main():
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    import pyarrow.dataset as ds
    import pandas as pd

    if not MUTATIONS.is_file():
        die(f"{MUTATIONS} not found.")

    say("=" * 72)
    say("MUTATIONS.parquet inspection")
    say("=" * 72)
    say(f"\nFile size on disk: {MUTATIONS.stat().st_size / 1_000_000_000:.2f} GB")

    # ---------------------------------------------------------------- structure
    pf = pq.ParquetFile(MUTATIONS)
    meta = pf.metadata
    say(f"Rows:       {meta.num_rows:,}")
    say(f"Columns:    {meta.num_columns}")
    say(f"Row groups: {meta.num_row_groups}")
    uncompressed = sum(
        meta.row_group(i).column(j).total_uncompressed_size
        for i in range(meta.num_row_groups)
        for j in range(meta.num_columns)
    )
    say(f"Uncompressed size: {uncompressed / 1_000_000_000:.2f} GB")

    say("\nColumns and types:")
    for field in pf.schema_arrow:
        note = "   (index level)" if field.name in INDEX_LEVELS else ""
        say(f"  {field.name:28s} {str(field.type):45s}{note}")

    # ------------------------------------------------------- gene vocabulary
    say("")
    counts = gene_vocabulary(pf, pa, pc)

    # ------------------------------------------------------------ target match
    say("\nMatching target loci against the file's own vocabulary:")
    resolved = {}
    for target, purpose in TARGETS.items():
        hits = [g for g in counts if str(g).lower() == target]
        if not hits:
            hits = [g for g in counts if target in str(g).lower()]
        if len(hits) == 1:
            resolved[target] = hits[0]
            say(f"  {target:8s} -> {hits[0]:10s} {counts[hits[0]]:>10,} rows   ({purpose})")
        elif not hits:
            say(f"  {target:8s} -> NOT FOUND   ({purpose})")
        else:
            say(f"  {target:8s} -> AMBIGUOUS: {hits}   ({purpose})")

    missing = [t for t in TARGETS if t not in resolved]
    if missing:
        die(f"target loci not resolved unambiguously: {missing}. Do not proceed.")

    gene_labels = sorted(resolved.values())

    # -------------------------------------------------- read target rows only
    say(f"\nReading only rows for {gene_labels}.")
    dataset = ds.dataset(MUTATIONS, format="parquet")
    df = flatten(dataset.to_table(filter=ds.field("GENE").isin(gene_labels)).to_pandas())
    say(f"  rows read: {len(df):,}")
    say(f"  memory held: {df.memory_usage(deep=True).sum() / 1_000_000:.0f} MB")

    # ------------------------------------------------------------- per gene
    say("\nPer gene:")
    say(f"  {'gene':10s} {'rows':>9s} {'samples':>9s} {'distinct mut':>13s} {'minor':>9s} {'null':>8s}")
    for gene in gene_labels:
        sub = df[df.GENE == gene]
        say(
            f"  {gene:10s} {len(sub):9,d} {sub.UNIQUEID.nunique():9,d} "
            f"{sub.MUTATION.nunique():13,d} {int(sub.IS_MINOR.sum()):9,d} {int(sub.IS_NULL.sum()):8,d}"
        )

    # -------------------------------------------------- coding vs non-coding
    say("\nCODES_PROTEIN by gene:")
    say(pd.crosstab(df.GENE, df.CODES_PROTEIN).to_string())

    # ------------------------------------------------------ read support (FRS)
    say("\nFRS, the fraction of reads supporting each call:")
    say(f"  null FRS values: {int(df.FRS.isna().sum()):,} of {len(df):,}")
    say(df.FRS.describe(percentiles=[0.01, 0.05, 0.25, 0.5, 0.95]).to_string())
    say("\n  FRS by IS_MINOR:")
    say(df.groupby("IS_MINOR", observed=True).FRS.describe()[["count", "min", "50%", "max"]].to_string())

    say("\n  major-allele rows (IS_MINOR False) per gene, and their FRS floor:")
    for gene in gene_labels:
        sub = df[(df.GENE == gene) & (~df.IS_MINOR)]
        floor = f"{sub.FRS.min():.3f}" if sub.FRS.notna().any() else "all null"
        say(f"    {gene:10s} {len(sub):8,d} rows   min FRS {floor}")

    # ------------------------------------------------------- mutation shapes
    say("\nExample mutation strings per gene (first eight distinct):")
    for gene in gene_labels:
        say(f"  {gene:10s} {df[df.GENE == gene].MUTATION.drop_duplicates().head(8).tolist()}")

    # ------------------------------------------------------------ join check
    say("\nJoin coverage against the phenotype and genome tables:")
    phenotypes = pd.read_parquet(DATA / "UKMYC_PHENOTYPES.parquet").reset_index()
    genomes = pd.read_parquet(DATA / "GENOMES.parquet").reset_index()
    mic_ids = set(phenotypes.UNIQUEID)
    genome_ids = set(genomes.UNIQUEID)
    matched = mic_ids & genome_ids
    mut_ids = set(df.UNIQUEID)

    say(f"  samples with a UKMYC MIC:                 {len(mic_ids):,}")
    say(f"  samples in GENOMES:                       {len(genome_ids):,}")
    say(f"  samples with both:                        {len(matched):,}")
    say(f"  samples with a target-gene mutation:      {len(mut_ids):,}")
    say(f"  of those, also have MIC and genome:       {len(mut_ids & matched):,}")
    say(f"  matched samples with no target mutation:  {len(matched - mut_ids):,}")
    say("\n  (the last line is the wild-type reference group the atlas measures shifts against)")

    write_report()
    say(f"\nWritten to {REPORT}")


if __name__ == "__main__":
    main()
