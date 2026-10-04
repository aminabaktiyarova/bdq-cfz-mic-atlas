"""
HGVS names for the variants the atlas reports.

Run from the project root with the virtual environment active:

    python code/hgvs_names.py

Variant names in this project are GARC, the grammar the CRyPTIC tables are
written in. TB-Profiler, the WHO catalogue and most of the literature name the
same variants in HGVS. Without a translation the atlas cannot be joined to any
of them by variant name, which the results document records as a limitation.

This module translates. It reads no catalogue: the mapping follows from the two
naming conventions and from nothing else, so what it writes carries no
restricted content and is released beside the rest.

What translates:

  A protein substitution. GARC writes the one-letter residues either side of
  the position, S2F. HGVS writes three-letter codes, p.Ser2Phe. The stop codon
  is ! in GARC and Ter in HGVS.

  A nucleotide substitution. GARC writes lowercase bases, c-11a. HGVS writes
  uppercase with the reference first, c.-11C>A, and a negative position carries
  the same meaning in both, a base upstream of the start codon. The prefix is
  c. for a gene that codes protein and n. for one that does not.

What does not translate, and why. An insertion or deletion is named at a
different layer by the two conventions. GARC gives the nucleotide change,
141_ins_c. The catalogues give the protein consequence, p.Ala110fs. Deriving
one from the other means translating the reference coding sequence to find the
codon at which the frame first breaks, and the result could not be checked
against those catalogues, because their frameshift entries do not record which
nucleotide change produced them. A whole or partial gene deletion has no HGVS
variant name at all; the catalogues record it as a consequence, feature
ablation. Rows of both kinds are written with no name and the reason recorded.

A null call and a het call are not variants and are refused rather than named.

The translation is checked rather than asserted. Every name this module emits
is translated back to GARC and compared against the string it came from, and a
row that does not survive the round trip fails the run. The test suite does the
same over every residue pair in the genetic code.

Outputs:
  outputs/hgvs_names.csv
  outputs/hgvs_names_report.txt
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import cohort  # noqa: E402

REPORT = Path("outputs/hgvs_names_report.txt")
TABLE = Path("outputs/hgvs_names.csv")
ATLAS = Path("outputs/atlas_evidence.csv")

THREE_LETTER = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
    "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val", "!": "Ter",
}
ONE_LETTER = {three: one for one, three in THREE_LETTER.items()}

# Z is a het call and X a null call. Neither is a residue, so neither is named.
NOT_A_RESIDUE = set("ZzXx")

GARC_PROTEIN = re.compile(r"^([A-Z!])(-?\d+)([A-Z!])$")
GARC_NUCLEOTIDE = re.compile(r"^([acgt])(-?\d+)([acgt])$")
GARC_INDEL = re.compile(r"^(-?\d+)_(ins|del)_(.+)$")
GARC_GENE_DELETION = re.compile(r"^del_(\d+(?:\.\d+)?|minorindel)$")

HGVS_PROTEIN = re.compile(r"^p\.([A-Z][a-z]{2})(-?\d+)([A-Z][a-z]{2})$")
HGVS_NUCLEOTIDE = re.compile(r"^[cn]\.(-?\d+)([ACGT])>([ACGT])$")

_lines = []


def say(text=""):
    print(text)
    _lines.append(text)


def to_hgvs(mutation, codes_protein=True):
    """The HGVS name for a GARC mutation, or None and the reason there is none."""
    mutation = str(mutation)

    match = GARC_PROTEIN.match(mutation)
    if match:
        reference, position, alternate = match.groups()
        if not codes_protein:
            return None, "a protein change in a gene that does not code protein"
        if reference in NOT_A_RESIDUE or alternate in NOT_A_RESIDUE:
            return None, "a null or het call is not a residue"
        if reference not in THREE_LETTER or alternate not in THREE_LETTER:
            return None, "residue outside the genetic code"
        return (f"p.{THREE_LETTER[reference]}{position}"
                f"{THREE_LETTER[alternate]}"), None

    match = GARC_NUCLEOTIDE.match(mutation)
    if match:
        reference, position, alternate = match.groups()
        prefix = "c" if codes_protein else "n"
        return f"{prefix}.{position}{reference.upper()}>{alternate.upper()}", None

    if GARC_INDEL.match(mutation):
        return None, ("an indel is named by its nucleotide change here and by "
                      "its protein consequence in the catalogues")
    if GARC_GENE_DELETION.match(mutation):
        return None, "a gene deletion has no HGVS variant name"
    return None, "not a form this translation covers"


def to_garc(name, codes_protein=True):
    """The GARC mutation for an HGVS name, which is the inverse of to_hgvs."""
    match = HGVS_PROTEIN.match(str(name))
    if match:
        reference, position, alternate = match.groups()
        if reference in ONE_LETTER and alternate in ONE_LETTER:
            return f"{ONE_LETTER[reference]}{position}{ONE_LETTER[alternate]}"
        return None
    match = HGVS_NUCLEOTIDE.match(str(name))
    if match:
        position, reference, alternate = match.groups()
        return f"{reference.lower()}{position}{alternate.lower()}"
    return None


def round_trip_holds(mutation, codes_protein=True):
    """Whether a name translates back to the string it came from."""
    name, _ = to_hgvs(mutation, codes_protein)
    return name is not None and to_garc(name, codes_protein) == str(mutation)


def name_table(variants):
    """One row per variant, with its name or the reason it has none.

    variants is an iterable of (gene, mutation) pairs. Every name is checked by
    translating it back, and a failure raises rather than being written.
    """
    import pandas as pd

    records = []
    for gene, mutation in variants:
        name, reason = to_hgvs(mutation, codes_protein=True)
        if name is not None and not round_trip_holds(mutation):
            raise ValueError(
                f"{gene}@{mutation} translated to {name}, which does not "
                "translate back to the string it came from")
        records.append({"gene": gene, "mutation": mutation,
                        "hgvs": name, "reason": reason})
    return pd.DataFrame(records).sort_values(["gene", "mutation"],
                                             ignore_index=True)


def main():
    import pandas as pd

    say("=" * 72)
    say("HGVS names for the variants the atlas reports")
    say("=" * 72)

    if not ATLAS.is_file():
        raise FileNotFoundError(f"{ATLAS} not found; run code/build_atlas.py first")
    atlas = pd.read_csv(ATLAS)
    table = name_table(zip(atlas.gene, atlas.mutation))

    named = table.hgvs.notna()
    say(f"\nVariants in the atlas: {len(table)}")
    say(f"  named:   {int(named.sum())}")
    say(f"  unnamed: {int((~named).sum())}")

    say("\nWhy a variant carries no name:")
    for reason, count in cohort.ranked_counts(table.reason).items():
        say(f"  {count:>4}  {reason}")

    say("\nBy gene:")
    say(f"  {'gene':9s} {'named':>6s} {'unnamed':>8s}")
    for gene, block in table.groupby("gene", observed=True):
        say(f"  {gene:9s} {int(block.hgvs.notna().sum()):>6} "
            f"{int(block.hgvs.isna().sum()):>8}")

    say("\nEvery name translated back to the string it came from: "
        f"{int(named.sum())} of {int(named.sum())}.")

    say("\nExamples:")
    for row in table[named].head(4).itertuples():
        say(f"  {row.gene}@{row.mutation:12s} {row.hgvs}")

    TABLE.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(TABLE, index=False)
    REPORT.write_text("\n".join(_lines) + "\n")
    print(f"\nWritten: {TABLE}, {REPORT}")


if __name__ == "__main__":
    main()
