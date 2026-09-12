"""
Parser for mutation strings in the GARC grammar used by the CRyPTIC tables.

GARC is the grammar the Fowler lab's gnomonicus, gumpy and piezo tools emit, and
the MUTATION column of MUTATIONS.parquet is written in it. The rules below are
taken from NOMENCLATURE.md in the piezo repository (oxfordmmm/piezo), which is
the authoritative definition:

  Amino acids are always uppercase, nucleotides always lowercase.

  Het calls are given by the letter z, hence Z for an amino acid or z for a
  base. Null calls are given by the letter x, hence X for an amino acid or x
  for a base.

  ! is reserved for the STOP codon.

  Position is context dependent. For a coding sequence it is the amino acid
  residue number. For a gene encoding rRNA, and for a promoter, it is the
  nucleotide number. A negative position is within the promoter.

  Insertions and deletions take the form position_ins_bases or position_del_n.
  A frameshift is an insertion or deletion whose length is not divisible by
  three.

  A whole-gene deletion is written del_1.0, and a partial one del_0.<number>,
  where the number is the fraction of the gene deleted. So del_0.81 means 81
  per cent of the gene is absent.

Two further notes.

The form <position>_minorindel appears in the CRyPTIC tables but is not defined
in NOMENCLATURE.md, and its definition was not located in the piezo, gumpy or
gnomonicus source. It is parsed here as an indel of unresolved length at that
position, and IS_FRAMESHIFT is left as None rather than False, because an
unknown length is not the same as a length divisible by three. Callers must test
for True explicitly rather than relying on falsiness. code/check_parsing.py
profiles these rows against IS_MINOR, FRS and the indel columns so the reading
rests on evidence rather than on the token's name.

A fractional gene deletion has no single position and no residues, so POSITION,
REF_RESIDUE and ALT_RESIDUE are None for it and DELETED_FRACTION carries the
figure instead.

The consequence that matters most: X and Z are not residues. A mutation string
ending in X is a position that could not be called, and one ending in Z is a
heterozygous call. Neither is evidence that a variant is present, and neither is
evidence that one is absent. Counting them as substitutions would invent
mutations; ignoring them silently would treat unknown positions as wild type.
Both are wrong, so this parser labels them explicitly and leaves the decision to
the caller.
"""

import re

# Amino acids, including O (pyrrolysine), the stop codon, and the X and Z
# placeholders. Uppercase throughout, per the grammar.
AMINO_ACIDS = "ACDEFGHIKLMNOPQRSTVWXYZ!"
NUCLEOTIDES = "acgtxz"

NULL_CHARS = {"X", "x"}
HET_CHARS = {"Z", "z"}

AA_SNP = re.compile(rf"^([{AMINO_ACIDS}])(-?\d+)([{AMINO_ACIDS}])$")
NT_SNP = re.compile(rf"^([{NUCLEOTIDES}])(-?\d+)([{NUCLEOTIDES}])$")
INDEL = re.compile(rf"^(-?\d+)_(ins|del)_([{NUCLEOTIDES}]+|\d+)$")
MINOR_INDEL = re.compile(r"^(-?\d+)_minorindel$")
GENE_DELETION = re.compile(r"^del_(\d+(?:\.\d+)?)$")

# Fields every parsed mutation carries, so callers can rely on the shape.
FIELDS = (
    "PARSED",        # bool, whether the string matched the grammar at all
    "KIND",          # SNP or INDEL
    "AFFECTS",       # CDS, RNA or PROM
    "REF_RESIDUE",   # reference amino acid or base, None for indels
    "POSITION",      # int, amino acid number for CDS, nucleotide number otherwise
    "ALT_RESIDUE",   # alternate amino acid or base, None for indels
    "INDEL_TYPE",    # ins or del, None for SNPs
    "INDEL_SIZE",    # int number of bases, None for SNPs
    "IS_SYNONYMOUS", # bool, reference and alternate residues identical
    "IS_NULL_CALL",  # bool, X or x present
    "IS_HET_CALL",   # bool, Z or z present
    "IS_STOP",       # bool, alternate residue is the stop codon
    "IS_FRAMESHIFT", # bool or None, indel length not divisible by three, None if length unknown
    "IS_REAL_VARIANT",  # bool, a called, non-synonymous change
    "SIZE_RESOLVED",    # bool, whether the indel length is actually known
    "DELETED_FRACTION", # float, fraction of the gene deleted, gene deletions only
)


def _blank(parsed=False):
    return dict.fromkeys(FIELDS, None) | {
        "PARSED": parsed,
        "IS_SYNONYMOUS": False,
        "IS_NULL_CALL": False,
        "IS_HET_CALL": False,
        "IS_STOP": False,
        "IS_FRAMESHIFT": False,
        "IS_REAL_VARIANT": False,
        "SIZE_RESOLVED": True,
        "DELETED_FRACTION": None,
    }


def parse_mutation(mutation, codes_protein):
    """
    Parse one GARC mutation string.

    codes_protein is the CODES_PROTEIN flag from the same row, used to tell a
    coding sequence from ribosomal RNA. A negative position always means the
    promoter regardless of that flag.

    Returns a dict with the keys in FIELDS. PARSED is False if the string did
    not match the grammar, which callers must treat as an error rather than
    skipping the row.
    """
    if mutation is None:
        return _blank()
    text = str(mutation).strip()

    match = AA_SNP.match(text)
    if match:
        ref, position, alt = match.group(1), int(match.group(2)), match.group(3)
        result = _blank(True)
        result.update(
            KIND="SNP",
            AFFECTS="PROM" if position < 0 else "CDS",
            REF_RESIDUE=ref,
            POSITION=position,
            ALT_RESIDUE=alt,
            IS_NULL_CALL=bool({ref, alt} & NULL_CHARS),
            IS_HET_CALL=bool({ref, alt} & HET_CHARS),
            IS_STOP=alt == "!",
        )
        uncalled = result["IS_NULL_CALL"] or result["IS_HET_CALL"]
        result["IS_SYNONYMOUS"] = (ref == alt) and not uncalled
        result["IS_REAL_VARIANT"] = (ref != alt) and not uncalled
        return result

    match = NT_SNP.match(text)
    if match:
        ref, position, alt = match.group(1), int(match.group(2)), match.group(3)
        result = _blank(True)
        if position < 0:
            affects = "PROM"
        else:
            affects = "CDS" if codes_protein else "RNA"
        result.update(
            KIND="SNP",
            AFFECTS=affects,
            REF_RESIDUE=ref,
            POSITION=position,
            ALT_RESIDUE=alt,
            IS_NULL_CALL=bool({ref, alt} & NULL_CHARS),
            IS_HET_CALL=bool({ref, alt} & HET_CHARS),
        )
        uncalled = result["IS_NULL_CALL"] or result["IS_HET_CALL"]
        # A nucleotide change in rRNA or a promoter is never synonymous: there
        # is no codon for it to be silent in.
        result["IS_SYNONYMOUS"] = (ref == alt) and not uncalled and affects == "CDS"
        result["IS_REAL_VARIANT"] = (ref != alt) and not uncalled
        return result

    match = INDEL.match(text)
    if match:
        position, indel_type, payload = int(match.group(1)), match.group(2), match.group(3)
        size = int(payload) if payload.isdigit() else len(payload)
        affects = "PROM" if position < 0 else ("CDS" if codes_protein else "RNA")
        result = _blank(True)
        result.update(
            KIND="INDEL",
            AFFECTS=affects,
            POSITION=position,
            INDEL_TYPE=indel_type,
            INDEL_SIZE=size,
            # Only a coding sequence has a reading frame to shift. An indel in a
            # promoter or in rRNA changes length without shifting any frame.
            IS_FRAMESHIFT=(affects == "CDS") and (size % 3) != 0,
            IS_REAL_VARIANT=True,
        )
        return result

    match = MINOR_INDEL.match(text)
    if match:
        position = int(match.group(1))
        result = _blank(True)
        result.update(
            KIND="INDEL",
            AFFECTS="PROM" if position < 0 else ("CDS" if codes_protein else "RNA"),
            POSITION=position,
            INDEL_TYPE=None,
            INDEL_SIZE=None,
            SIZE_RESOLVED=False,
            IS_FRAMESHIFT=None,
            IS_REAL_VARIANT=True,
        )
        return result

    match = GENE_DELETION.match(text)
    if match:
        result = _blank(True)
        result.update(
            KIND="GENE_DELETION",
            AFFECTS="GENE",
            DELETED_FRACTION=float(match.group(1)),
            IS_REAL_VARIANT=True,
        )
        return result

    return _blank(False)


def parse_frame(df, mutation_column="MUTATION", codes_column="CODES_PROTEIN"):
    """
    Parse a whole DataFrame of mutations and return it with the parsed fields
    added as columns. Distinct mutation strings are parsed once and reused,
    which matters because the same variant recurs across thousands of samples.
    """
    import pandas as pd

    keys = df[[mutation_column, codes_column]].drop_duplicates()
    parsed = [parse_mutation(m, c) for m, c in zip(keys[mutation_column], keys[codes_column])]
    lookup = pd.concat([keys.reset_index(drop=True), pd.DataFrame(parsed)], axis=1)
    return df.merge(lookup, on=[mutation_column, codes_column], how="left", validate="many_to_one")
