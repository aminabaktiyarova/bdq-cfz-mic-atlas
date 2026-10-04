"""
The genotype layer against the WHO catalogue, second edition.

Run from the project root with the virtual environment active:

    python code/benchmark_catalogue.py

Every other analysis in this project classifies variants from the mutation
strings and reads no catalogue. This module is the one comparison against an
external grading, and it answers two questions on this cohort: how much of the
measured resistance the catalogue's calls recover, and which variants it grades
differently from what the MIC layer measures.

Licence boundary. The catalogue is CC BY-NC-SA 3.0 IGO, whose non-commercial and
share-alike terms are incompatible with this project's CC BY 4.0 release. The
catalogue file is read from quarantine/ and every output of this module is
written to quarantine/ as well, including the aggregate counts. Whether an
aggregate agreement count is a reproduction of the catalogue's variant-to-grade
mapping or a fact about a comparison is unsettled, so nothing here enters
outputs/ and nothing here is released without that decision being taken first.
The module raises if an output path falls outside quarantine/.

Catalogue files, obtained from github.com/oxfordmmm/tuberculosis_amr_catalogues
at commit cf60292082da179f072479270dfa09ba4567b800 and recorded with their
checksums in quarantine/catalogues/PROVENANCE.txt:

  v2.1  primary. The current conversion.
  v2.0  sensitivity. Differs by 171 rows, all additions, none regraded. For the
        two drugs here the additions are 6 rows each, every one graded R and
        every one a minor-allele form of a loss-of-function wildcard in Rv0678
        or pepQ, which v2.0 cannot express.

Both files carry CATALOGUE_NAME WHO-UCN-GTB-PCI-2023.5 and CATALOGUE_VERSION
2.0 in every row, which is the catalogue CRyPTIC applied to these same samples
per the EFFECTS and PREDICTIONS tables.

How an isolate is called. piezo grades a mutation, not a sample, so the
aggregation is this module's own and is stated rather than inherited:

  R if any graded variant is R, else U if any is U, else F if any is F, else S.

A mutation the catalogue does not list for a drug contributes nothing, and an
isolate carrying no graded variant is S. That is the catalogue's own logic: the
absence of an association is not evidence of resistance.

Which rows are graded. A called, non-synonymous change carried in every read:
IS_REAL_VARIANT and not a minor allele. A null call is a position that could not
be read and is not a detected variant, and a minor allele is reported
Susceptible by CRyPTIC's own application of this catalogue. Section 7.8 treats
both the same way, so the rules here and there are comparable.

Rows the cohort marks DOUBLE_REPORTED are graded here, which is the one place
this module departs from cohort.REAL_MAJOR. A large deletion is written twice in
v3.4.0, once as del_<fraction> and once as the sequence removed, and
cohort.classify drops the second so that one event counts as one variant. The
two spellings need not grade alike, because they match different rules: a
fractional deletion matches only a rule written for the fraction it states,
while the sequence removed is read as an indel and can match a frameshift
wildcard. One cohort isolate carries a partial Rv0678 deletion written both
ways and the two spellings reach different rules. Dropping the duplicate would
hide whichever grade the catalogue assigns to the event. Grading both spellings cannot double-count an isolate,
because the call takes the strongest grade the isolate carries.

The epistasis rule is reported both ways. CRyPTIC apply a rule under which a
loss-of-function mutation in mmpL5 overrides any resistance-associated mutation
in Rv0678 and the isolate is called Susceptible. Section 11.1 records that this
project's own epistasis finding was retracted: the 43 samples carrying both are
2 independent events, not 43. Adopting the rule silently would import a claim
this cohort does not support, and ignoring it would misrepresent what the
catalogue does in practice, so both calls are computed and reported side by
side.

For the comparison against the measured phenotype, U is not a call of
resistance. A three-way breakdown of R, U and S against the phenotype is
reported beside the two-way confusion, because an isolate graded U is a
different failure from one graded S.

Intervals resample clusters, through the same functions Section 7.8 uses, so
the catalogue's sensitivity and specificity sit on the same footing as the
project's own rules.

Outputs, all inside quarantine/:
  quarantine/benchmark/benchmark_report.txt
  quarantine/benchmark/benchmark_metrics.csv
  quarantine/benchmark/variant_grades.csv
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import cohort  # noqa: E402
from cluster_adjust import add_clusters  # noqa: E402
from prediction_metrics import (  # noqa: E402
    cluster_interval, confusion, metrics, stream)

QUARANTINE = Path("quarantine")
CATALOGUES = QUARANTINE / "catalogues"
OUTPUT_DIR = QUARANTINE / "benchmark"

REPORT = OUTPUT_DIR / "benchmark_report.txt"
METRICS = OUTPUT_DIR / "benchmark_metrics.csv"
GRADES = OUTPUT_DIR / "variant_grades.csv"

CATALOGUE_FILES = {
    "v2.1": CATALOGUES / "NC_000962.3_WHO-UCN-TB-2023.5_v2.1_GARC1_RFUS.csv",
    "v2.0": CATALOGUES / "NC_000962.3_WHO-UCN-TB-2023.5_v2.0_GARC1_RFUS.csv",
}
PRIMARY = "v2.1"

# Highest first. A grade earlier in this list wins when an isolate carries
# variants of more than one grade.
PRECEDENCE = ["R", "U", "F", "S"]

_lines = []


def say(text=""):
    print(text)
    _lines.append(text)


def check_output_paths():
    """Refuse to write catalogue-derived material outside quarantine/."""
    for path in (REPORT, METRICS, GRADES):
        if QUARANTINE not in path.parents:
            raise ValueError(f"{path} is outside {QUARANTINE}")


def load_catalogue(version=PRIMARY, path=None):
    """The catalogue as piezo reads it, with its identity checked."""
    import piezo

    path = Path(path) if path is not None else CATALOGUE_FILES[version]
    if not path.is_file():
        raise FileNotFoundError(
            f"{path} not found. See quarantine/catalogues/PROVENANCE.txt for "
            "where the file comes from. It is never committed.")
    catalogue = piezo.ResistanceCatalogue(str(path))
    if catalogue.catalogue.name != "WHO-UCN-GTB-PCI-2023.5":
        raise ValueError(f"unexpected catalogue name {catalogue.catalogue.name}")
    if catalogue.catalogue.grammar != "GARC1":
        raise ValueError(f"unexpected grammar {catalogue.catalogue.grammar}")
    return catalogue


def graded_rows(mutations):
    """The rows whose grade the catalogue is asked for.

    cohort.REAL_MAJOR additionally excludes the duplicate spelling of a
    double-reported deletion, which is right for counting variants and wrong
    for grading them. See the module docstring.
    """
    return mutations.IS_REAL_VARIANT & ~mutations.IS_MINOR


def grade_variants(mutations, catalogue):
    """One row per distinct major-allele variant, with its grade for each drug.

    A drug the catalogue does not associate with the variant's gene returns no
    key, which is recorded as "not listed" and counts as S when an isolate is
    called.
    """
    import pandas as pd

    pairs = (mutations.loc[graded_rows(mutations), ["GENE", "MUTATION", "CLASS"]]
             .drop_duplicates(subset=["GENE", "MUTATION"]))
    columns = ["GENE", "MUTATION", "CLASS"] + list(cohort.DRUGS)
    if pairs.empty:
        return pd.DataFrame(columns=columns)

    records = []
    for gene, mutation, variant_class in pairs.itertuples(index=False):
        prediction = catalogue.predict(f"{gene}@{mutation}")
        if not isinstance(prediction, dict):
            raise ValueError(f"{gene}@{mutation} returned {prediction!r}")
        record = {"GENE": gene, "MUTATION": mutation, "CLASS": variant_class}
        for drug in cohort.DRUGS:
            record[drug] = prediction.get(drug, "not listed")
        records.append(record)
    return pd.DataFrame(records, columns=columns).sort_values(
        ["GENE", "MUTATION"], ignore_index=True)


def isolate_calls(mutations, grades, index, status, epistasis=False):
    """The catalogue's call for each isolate and drug.

    With epistasis true, an isolate whose mmpL5 carries a loss-of-function
    variant is called S whatever its other variants grade, which is the rule
    CRyPTIC apply.
    """
    import pandas as pd

    major = mutations[graded_rows(mutations)]
    joined = major.merge(grades, on=["GENE", "MUTATION"], how="left",
                         suffixes=("", "_graded"))
    if joined[cohort.DRUGS[0]].isna().any():
        raise ValueError("a graded row reached no grade")

    calls = pd.DataFrame(index=index)
    for drug in cohort.DRUGS:
        rank = joined[drug].map(
            {grade: position for position, grade in enumerate(PRECEDENCE)})
        # "not listed" contributes nothing, so it takes the lowest rank.
        rank = rank.fillna(len(PRECEDENCE) - 1)
        best = rank.groupby(joined.UNIQUEID, observed=True).min()
        column = best.reindex(index).fillna(len(PRECEDENCE) - 1).astype(int)
        calls[drug] = [PRECEDENCE[position] for position in column]
        if epistasis:
            override = status.mmpL5_LOF.reindex(index).fillna(False).to_numpy()
            calls[drug] = calls[drug].where(~override, "S")
    return calls


def evaluate(frame, calls, drug, label):
    """The catalogue as a test for resistance, on the isolates carrying a MIC."""
    subset = frame[frame[f"MIC_{drug}"].notna()].copy()
    subset["called"] = calls[drug].reindex(subset.index).eq("R")
    resistant = f"resistant_{drug}"

    counts = confusion(subset["called"], subset[resistant])
    values = metrics(counts)
    intervals = cluster_interval(
        subset, "called", resistant, stream(f"benchmark {label} {drug}"))

    record = {"catalogue": label, "drug": drug, "isolates": len(subset),
              "clusters": subset.CLUSTER.nunique(),
              "resistant": int(subset[resistant].sum()), **counts}
    for key, value in values.items():
        record[key] = round(value, 4) if value is not None else None
        interval = intervals[key]
        record[f"{key}_low"] = round(interval[0], 4) if interval else None
        record[f"{key}_high"] = round(interval[1], 4) if interval else None
    return record, subset


def three_way(subset, calls, drug):
    """Counts of each grade against the measured phenotype."""
    import pandas as pd

    grade = calls[drug].reindex(subset.index)
    return pd.crosstab(grade, subset[f"resistant_{drug}"]).reindex(
        index=[g for g in PRECEDENCE if g in set(grade)], fill_value=0)


def main():
    import pandas as pd

    check_output_paths()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    say("=" * 72)
    say("The genotype layer against the WHO catalogue, second edition")
    say("=" * 72)
    say("Catalogue-derived. This report stays in quarantine/.")

    mutations = cohort.load_mutations()
    status = cohort.build_status(mutations)
    frame = add_clusters(cohort.assemble(status), mutations)

    records = []
    grade_tables = {}
    for version in CATALOGUE_FILES:
        catalogue = load_catalogue(version)
        grades = grade_variants(mutations, catalogue)
        grade_tables[version] = grades

        say(f"\n{'-' * 72}")
        say(f"Catalogue {version}: {len(grades):,} distinct major-allele variants graded")
        say("-" * 72)
        for drug in cohort.DRUGS:
            counts = cohort.ranked_counts(grades[drug])
            say(f"  {drug}: " + ", ".join(
                f"{grade} {count:,}" for grade, count in counts.items()))

        for epistasis in (False, True):
            label = f"{version} epistasis" if epistasis else version
            calls = isolate_calls(mutations, grades, frame.index, status,
                                  epistasis=epistasis)
            for drug in cohort.DRUGS:
                record, subset = evaluate(frame, calls, drug, label)
                records.append(record)
                if version == PRIMARY and not epistasis:
                    say(f"\n  {drug}, grade against measured phenotype:")
                    table = three_way(subset, calls, drug)
                    say("    grade  susceptible  resistant")
                    for grade, row in table.iterrows():
                        say(f"    {grade:5s} {row.get(False, 0):12,} {row.get(True, 0):10,}")

    table = pd.DataFrame(records)

    say(f"\n{'=' * 72}")
    say("The catalogue as a test for resistance, R against everything else")
    say("=" * 72)
    say(f"\n  {'catalogue':20s} {'drug':4s} {'sens':>18s} {'spec':>18s} {'PPV':>18s}")
    for row in table.itertuples():
        cells = []
        for key in ("sensitivity", "specificity", "ppv"):
            value = getattr(row, key)
            low = getattr(row, f"{key}_low")
            high = getattr(row, f"{key}_high")
            cells.append(f"{100 * value:5.1f} ({100 * low:4.1f} to {100 * high:4.1f})"
                         if value is not None and low is not None else "n/a")
        say(f"  {row.catalogue:20s} {row.drug:4s} " + " ".join(f"{c:>18s}" for c in cells))

    say("\n  cells, as true positive, false positive, false negative, true negative:")
    for row in table.itertuples():
        say(f"  {row.catalogue:20s} {row.drug:4s} {row.true_positive:5d} "
            f"{row.false_positive:6d} {row.false_negative:5d} {row.true_negative:7d}")

    primary = grade_tables[PRIMARY].copy()
    other = grade_tables["v2.0"]
    merged = primary.merge(other, on=["GENE", "MUTATION"], suffixes=("", "_v2.0"))
    moved = merged[(merged.BDQ.ne(merged["BDQ_v2.0"]))
                   | (merged.CFZ.ne(merged["CFZ_v2.0"]))]
    say(f"\nVariants this cohort carries whose grade differs between the two "
        f"conversions: {len(moved)}")

    primary.to_csv(GRADES, index=False)
    table.to_csv(METRICS, index=False)
    REPORT.write_text("\n".join(_lines) + "\n")
    say(f"\nWritten: {REPORT}, {METRICS}, {GRADES}")
    REPORT.write_text("\n".join(_lines) + "\n")


if __name__ == "__main__":
    main()
