"""
Canonical cohort definitions for the atlas.

Every analysis imports from here rather than redefining a cohort, so there is
exactly one place where "solo", "reference" and "loss of function" are decided.

Nothing here reads the WHO catalogue, the EFFECTS table or the PREDICTIONS
table. Every classification comes from the mutation strings themselves, so
anything built on this module carries only CC BY 4.0 CRyPTIC-derived content and
can be released.

Definitions:

  Real variant. A called, non-synonymous change: substitution, indel, or gene
  deletion. Null calls (X), het calls (Z) and synonymous substitutions are not
  variants.

  Major allele. IS_MINOR false. Where FRS exists it separates cleanly at 0.90,
  but FRS is absent from most rows, so IS_MINOR is the operative flag.

  Uncertain. The gene carries a null call, a het call, or any minor allele. A
  position that could not be read might carry anything, and a sub-population
  carrying a variant is not an absence of one. Such a sample cannot be asserted
  wild type at that gene, so it is excluded from both the reference and the solo
  groups rather than counted as carrying nothing.

  Variant class. Gene deletion, frameshift, stop codon, in-frame indel,
  promoter, or substitution. Kept separate rather than collapsed into "loss of
  function", because the classes do not behave alike and pooling them hides
  that.

  Solo. Exactly one real major-allele variant across Rv0678, pepQ and atpE, and
  nothing uncertain in any of them. mmpL5 is excluded from this count: it
  carries variants in almost every sample, and its role is to modify the effect
  of an Rv0678 variant rather than to cause resistance. It is carried as a
  covariate instead.

  Reference. No real major-allele variant in Rv0678, pepQ or atpE, and nothing
  uncertain in any of them.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from garc import GENE_CODES, gene_codes_protein, parse_frame  # noqa: E402

DATA = Path("data/cryptic-v3.4.0")

# Genes whose variants define the bedaquiline and clofazimine genotype.
BDQ_GENES = ["Rv0678", "pepQ", "atpE"]
# Carried as a covariate, never part of the solo or reference definition.
MODIFIER_GENE = "mmpL5"
ALL_GENES = BDQ_GENES + [MODIFIER_GENE]

DRUGS = ["BDQ", "CFZ"]
INDEX_LEVELS = ["UNIQUEID", "GENE", "MUTATION"]

# Classes counted as loss of function of the gene product.
LOF_CLASSES = ["gene deletion", "frameshift", "stop codon"]

# Group labels, in the order they should be reported.
GROUP_ORDER = [
    "reference",
    "Rv0678 gene deletion",
    "Rv0678 frameshift",
    "Rv0678 stop codon",
    "Rv0678 in-frame indel",
    "Rv0678 promoter",
    "Rv0678 substitution",
    "pepQ solo",
    "atpE solo",
    "multiple variants",
    "uncertain",
]


def flatten(df):
    """Move parquet index levels back into ordinary columns."""
    if [n for n in INDEX_LEVELS if n not in df.columns]:
        df = df.reset_index()
    return df


def ranked_counts(series, limit=None, dropna=True):
    """Value counts ordered by frequency, ties ordered by value.

    value_counts leaves the order of equally frequent values to the grouping,
    and that order differs between pandas versions, so a listing cut to a fixed
    length can take a different tied value on a different installation. Sorting
    the index and then sorting by count with a stable sort fixes both the order
    of the listing and which values fall inside the cut.
    """
    counts = series.value_counts(dropna=dropna)
    counts = counts.sort_index().sort_values(ascending=False, kind="stable")
    return counts if limit is None else counts.head(limit)


def double_reported_deletion(df):
    """Which rows report a deletion that another row already reports.

    A large deletion is written twice in v3.4.0, once as del_<fraction> and once
    as the sequence removed, so one event in one gene of one sample would count
    as two variants. The row carrying the sequence is marked and the
    del_<fraction> row kept, because that row carries the fraction and classifies
    as a gene deletion. Where a sample and gene hold more than one deletion
    beside a del_<fraction> row, the largest is the one written twice.
    """
    import numpy as np
    import pandas as pd

    keys = ["UNIQUEID", "GENE"]
    columns = keys + ["KIND", "INDEL_TYPE", "INDEL_SIZE", "IS_MINOR"]
    blank = pd.Series(np.zeros(len(df), dtype=bool), index=df.index)
    if [column for column in columns if column not in df.columns]:
        return blank
    work = df.reset_index(drop=True)
    whole = work.loc[work.KIND.eq("GENE_DELETION"), keys]
    if whole.empty:
        return blank
    # Every mask below is combined into a new object rather than updated in
    # place. pandas hands out read-only arrays from version 3, and an in-place
    # update of one raises.
    beside_one = pd.Series(
        pd.MultiIndex.from_frame(work[keys]).isin(
            set(map(tuple, whole.to_numpy()))),
        index=work.index)
    partner = (work.KIND.eq("INDEL") & work.INDEL_TYPE.eq("del")
               & ~work.IS_MINOR.eq(True) & work.INDEL_SIZE.notna()
               & beside_one)
    if not partner.any():
        return blank
    largest = work[partner].groupby(keys, observed=True).INDEL_SIZE.idxmax()
    marked = np.zeros(len(df), dtype=bool)
    marked[largest.to_numpy(dtype=int)] = True
    return pd.Series(marked, index=df.index)


def classify(df):
    """Add CLASS, IS_LOF, REAL_MAJOR and UNCERTAIN to a parsed mutation frame."""
    import numpy as np

    df = df.copy()
    df["DOUBLE_REPORTED"] = double_reported_deletion(df)
    df["REAL_MAJOR"] = df.IS_REAL_VARIANT & ~df.IS_MINOR & ~df.DOUBLE_REPORTED
    df["UNCERTAIN"] = df.IS_NULL_CALL | df.IS_HET_CALL | df.IS_MINOR
    # A promoter has no reading frame, so an indel confined to one is a change
    # to the promoter, and the promoter test comes before the indel test.
    df["CLASS"] = np.where(
        df.KIND.eq("GENE_DELETION"), "gene deletion",
        np.where(df.IS_FRAMESHIFT.eq(True), "frameshift",
        np.where(df.IS_STOP, "stop codon",
        np.where(df.AFFECTS.eq("PROM"), "promoter",
        np.where(df.KIND.eq("INDEL"), "in-frame indel", "substitution")))))
    df["IS_LOF"] = df.REAL_MAJOR & df.CLASS.isin(LOF_CLASSES)
    return df


def load_mutations(genes=None):
    """Read, parse and classify the target-gene mutations. Memory-safe."""
    import pyarrow.dataset as ds

    genes = genes or ALL_GENES
    path = DATA / "MUTATIONS.parquet"
    if not path.is_file():
        raise FileNotFoundError(f"{path} not found")
    dataset = ds.dataset(path, format="parquet")
    df = flatten(dataset.to_table(filter=ds.field("GENE").isin(genes)).to_pandas())
    # The grammar needs to know whether the gene codes protein, which the
    # table's own CODES_PROTEIN column does not say.
    df[GENE_CODES] = gene_codes_protein(df)
    df = parse_frame(df)
    if not df.PARSED.all():
        raise ValueError("some mutation strings did not parse; run code/check_parsing.py")
    return classify(df)


def build_status(mutations=None):
    """One row per genome, carrying genotype status and the group label."""
    import pandas as pd

    if mutations is None:
        mutations = load_mutations()

    per_gene = (
        mutations.groupby(["UNIQUEID", "GENE"], observed=True)
        .agg(
            n_variants=("REAL_MAJOR", "sum"),
            n_lof=("IS_LOF", "sum"),
            n_minor=("IS_MINOR", "sum"),
            uncertain=("UNCERTAIN", "any"),
        )
        .reset_index()
    )

    def spread(column, prefix):
        wide = per_gene.pivot(index="UNIQUEID", columns="GENE", values=column)
        wide = wide.reindex(columns=ALL_GENES)
        wide.columns = [f"{prefix}_{g}" for g in wide.columns]
        return wide

    status = pd.concat(
        [
            spread("n_variants", "variants"),
            spread("n_lof", "lof"),
            spread("uncertain", "uncertain"),
            spread("n_minor", "minor"),
        ],
        axis=1,
    )

    genomes = pd.read_parquet(DATA / "GENOMES.parquet").reset_index()
    status = status.reindex(genomes.UNIQUEID.unique())
    for column in status.columns:
        if column.startswith("uncertain_"):
            status[column] = status[column].fillna(False).astype(bool)
        else:
            status[column] = status[column].fillna(0).astype(int)

    variants = status[[f"variants_{g}" for g in BDQ_GENES]].sum(axis=1)
    uncertain = status[[f"uncertain_{g}" for g in BDQ_GENES]].any(axis=1)

    status["n_bdq_variants"] = variants
    status["any_uncertain"] = uncertain
    status["IS_SOLO"] = (variants == 1) & ~uncertain
    status["IS_REFERENCE"] = (variants == 0) & ~uncertain
    status["mmpL5_LOF"] = status["lof_mmpL5"] > 0

    # For solo samples, which class the single variant belongs to.
    solo_ids = status.index[status.IS_SOLO]
    solo_variants = mutations[
        mutations.REAL_MAJOR & mutations.UNIQUEID.isin(solo_ids)
        & mutations.GENE.isin(BDQ_GENES)
    ]
    labels = solo_variants.set_index("UNIQUEID")[["GENE", "CLASS"]]
    labels = labels[~labels.index.duplicated()]
    status = status.join(labels.rename(columns={"GENE": "solo_gene", "CLASS": "solo_class"}))

    status["GROUP"] = [
        _group_label(row) for row in status.itertuples()
    ]
    return status


def _group_label(row):
    if row.IS_REFERENCE:
        return "reference"
    if row.any_uncertain:
        return "uncertain"
    if not row.IS_SOLO:
        return "multiple variants"
    if row.solo_gene == "Rv0678":
        return f"Rv0678 {row.solo_class}"
    return f"{row.solo_gene} solo"


def assemble(status=None):
    """Join genotype status to the MICs, lineage and site. One row per sample."""
    import pandas as pd

    if status is None:
        status = build_status()

    phenotypes = pd.read_parquet(DATA / "UKMYC_PHENOTYPES.parquet").reset_index()
    phenotypes = phenotypes[phenotypes.DRUG.isin(DRUGS)]
    mic = phenotypes.pivot_table(
        index="UNIQUEID", columns="DRUG",
        values=["MIC", "LOG2MIC", "BINARY_PHENOTYPE", "PHENOTYPE_QUALITY", "PLATEDESIGN"],
        aggfunc="first", observed=True,
    )
    mic.columns = [f"{a}_{b}" for a, b in mic.columns]

    genomes = pd.read_parquet(DATA / "GENOMES.parquet").reset_index()
    lineage = genomes.set_index("UNIQUEID")[["LINEAGE", "SUBLINEAGE"]]
    site = phenotypes.drop_duplicates("UNIQUEID").set_index("UNIQUEID")[["SITEID"]]

    joined = status.join(mic, how="inner").join(lineage).join(site)
    joined["SITEID"] = joined.SITEID.astype(str)

    for drug in DRUGS:
        text = joined[f"MIC_{drug}"].astype(str)
        joined[f"censored_left_{drug}"] = text.str.startswith("<=")
        joined[f"censored_right_{drug}"] = text.str.startswith(">")
        joined[f"resistant_{drug}"] = joined[f"BINARY_PHENOTYPE_{drug}"].eq("R")

    return joined


def load_sites():
    """SITEID to country and institution. Returns None if the file is absent."""
    import pandas as pd

    path = DATA / "SITES.csv.gz"
    if not path.is_file():
        return None
    sites = pd.read_csv(path, dtype={"SITEID": str})
    sites["SITEID"] = sites.SITEID.str.strip()
    return sites.set_index("SITEID")
