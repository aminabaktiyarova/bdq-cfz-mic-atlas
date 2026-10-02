# bdq-cfz-mic-atlas

An open, versioned atlas of how mutations in the Rv0678 and MmpL5 efflux axis
shift bedaquiline and clofazimine minimum inhibitory concentrations in
*Mycobacterium tuberculosis*, built from the CRyPTIC Consortium's quantitative
phenotype data.

## What this is

Existing resistance catalogues grade a mutation as associated with resistance to
one drug or another. They do not quantify how far a mutation moves the MIC, and
they do not describe the shared-mechanism cross-resistance between bedaquiline
and clofazimine, which arise from the same efflux pump and the same regulatory
gene.

This project builds that quantitative layer: for each mutation and mutation
class in the axis, the distribution of bedaquiline and clofazimine MICs it is
observed with, the size of the shift relative to isolates carrying no graded
mutation, and how the two drugs' responses differ. It also releases the pipeline
that produces it, so the resource can be regenerated against each new CRyPTIC
release rather than becoming a snapshot.

MIC values are interval-censored: an isolate whose growth is inhibited at every
concentration tested has an MIC somewhere below the lowest well, and one growing
at every concentration has an MIC above the highest. These are handled as
censored observations throughout rather than imputed to the boundary value.

## Data

Built on the CRyPTIC Consortium Dataset, version v3.4.0, Zenodo version DOI
10.5281/zenodo.15680920, CC BY 4.0.

The dataset is not redistributed here. `docs/PROVENANCE.md` records the exact
version, licence, file inventory, checksums, schema and the discrepancies found
between the published schema document and the data as shipped, together with the
measured properties of the input tables. Every figure this project produces traces
back to that record.

## Results

`docs/RESULTS.md` records what was computed: the cohort and what it excludes, the
interval-censored estimator and its validations, the variant class effects before
and after holding site constant, the clonal-clustering correction and the
effective sample sizes it leaves, the per-variant layer and how few variants it
reaches, what detected minor alleles carry, the fraction of resistance no variant
in the three genes explains, the four genotype rules read as diagnostic tests, the
one pre-registered comparison and its result, and the limitations that apply to
all of it.

No resistance catalogue takes part in any figure. Every genotype is read from the
mutation strings in `MUTATIONS.parquet`.

## Repository layout

```
code/         analysis and pipeline code (MIT)
data/         CRyPTIC source tables, downloaded by the user, not tracked
quarantine/   restricted-licence and personal-data inputs, not tracked
docs/         provenance record, results, pre-registration, data dictionary
outputs/      derived tables for release (CC BY 4.0)
```

`data/` and `quarantine/` are excluded from version control. The separation is
physical so that restricted inputs cannot reach a released artifact. See
`LICENSES.md` for what is excluded and why.

## Reproducing the environment

```
python3.14 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python code/verify_setup.py
```

`verify_setup.py` checks the installed packages, confirms the CRyPTIC tables are
present, and recomputes known sample and resistance counts against the recorded
figures. It fails loudly if the local data is a different release.

## Running the analysis

Each module reads the CRyPTIC tables directly and writes to `outputs/`. None
reads another module's output, so they can be run in any order, with one
exception: `discovery.py` writes `docs/PRE_REGISTRATION.md`, and `validate.py`
refuses to run until that file is committed and unmodified.

```
python code/inspect_mutations.py    # structure of the MUTATIONS table
python code/check_parsing.py        # the GARC parser against CRyPTIC's own flags
python code/build_sample_status.py  # per-sample genotype status
python code/audit_cohort.py         # join integrity, replication, quality, the split
python code/analyse_groups.py       # variant class effects, site and lineage strata
python code/model_effects.py        # lineage against site, logistic model
python code/cluster_adjust.py       # the same effects corrected for clonal clustering
python code/mic_model.py            # interval-censored MIC distributions
python code/discovery.py            # discovery-half estimates, writes the pre-registration
python code/validate.py             # the held-out test, runs once
python code/build_atlas.py          # per-variant evidence and shifts
python code/unexplained.py          # resistance carrying no variant in the three genes
python code/heteroresistance.py     # the samples excluded as uncertain
python code/prediction_metrics.py   # genotype rules as tests for resistance
```

## Outputs

`outputs/` holds a readable report beside each table. Every column of every
table is defined in `docs/DATA_DICTIONARY.md`.

| Table | Written by |
| --- | --- |
| `atlas_evidence.csv` | `code/build_atlas.py` |
| `mic_estimates.csv` | `code/mic_model.py` |
| `unexplained_counts.csv`, `unexplained_gene_sets.csv` | `code/unexplained.py` |
| `heteroresistance_estimates.csv` | `code/heteroresistance.py` |
| `prediction_metrics.csv` | `code/prediction_metrics.py` |
| `gene_vocabulary.csv` | `code/inspect_mutations.py` |

Outputs are regenerated from the CRyPTIC release and are not tracked.

## Tests

```
pytest
```

The suite runs against a synthetic dataset built to the CRyPTIC schema with
known ground truth planted in it, so it needs no downloaded data and finishes in
about half a minute.

It checks the things the results depend on rather than the things that are easy
to check. That a null call and a het call are not variants and not wild type.
That a minor-allele indel excludes a sample from the reference group, as a het
call does. That a reported MIC is treated as the interval below it rather than
as a measurement, and that the censored estimator recovers parameters a median
cannot. That the Mantel-Haenszel estimate recovers a known odds ratio across
sites of different baseline risk, and that the homogeneity test fires when the
effect differs between them. That collapsing clonal clusters removes an effect
planted so as to be driven entirely by one outbreak. That the prediction
interval is never narrower than the discovery and held-out errors combined. That
a resistant isolate whose genes could not be called is counted neither as
carrying a variant nor as lacking one.

Each test was verified by breaking the code it covers and confirming the test
fails.

## Licensing

Code under MIT, derived data under CC BY 4.0. WHO catalogue content and the
CRyPTIC `EFFECTS`, `PREDICTIONS` and `BASHTHEBUG_CLASSIFICATIONS` tables are
excluded from every released artifact. Full statement in `LICENSES.md`.

## Author

Amina Baktiyarova, Independent Researcher
ORCID 0009-0007-6265-6493

## Status

No versioned release has been made. The pipeline, the derived tables and
`docs/RESULTS.md` are current with each other; the atlas table has no citable
identifier of its own yet.
