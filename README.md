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
between the published schema document and the data as shipped. Every figure this
project produces traces back to that record.

## Repository layout

```
code/         analysis and pipeline code (MIT)
data/         CRyPTIC source tables, downloaded by the user, not tracked
quarantine/   restricted-licence and personal-data inputs, not tracked
docs/         provenance record and methodology notes
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
interval is never narrower than the discovery and held-out errors combined.

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

Under construction. No release has been made and no results are published yet.
