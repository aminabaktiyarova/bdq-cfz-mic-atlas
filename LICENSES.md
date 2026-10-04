# Licensing

Different parts of this project carry different licences. This file states which
applies to what, and why.

## Code

Everything under `code/` is licensed under the MIT License. The full text is in
`LICENSE`.

## Data

Derived data tables released by this project are licensed under Creative Commons
Attribution 4.0 International (CC BY 4.0).

https://creativecommons.org/licenses/by/4.0/

CC0 is not available for these tables. They are derived from the CRyPTIC
Consortium Dataset, which is released under CC BY 4.0, and the attribution
requirement of that licence cannot be waived by a downstream user. Any table
containing CRyPTIC-derived content therefore carries CC BY 4.0.

## Source data

The CRyPTIC Consortium Dataset, version v3.4.0, Zenodo version DOI
10.5281/zenodo.15680920, licensed CC BY 4.0. Curated by Philip Fowler,
University of Oxford, ORCID 0000-0003-0912-4483, on behalf of the CRyPTIC
Consortium.

The dataset is not redistributed in this repository. It is downloaded from
Zenodo by the user. `docs/PROVENANCE.md` records the exact version, licence,
file inventory and checksums that the analysis depends on.

## Excluded inputs

Two categories of input are deliberately kept out of this repository and out of
every released artifact.

**WHO catalogue content.** The WHO Catalogue of mutations in Mycobacterium
tuberculosis complex and their association with drug resistance is licensed
CC BY-NC-SA 3.0 IGO. Its non-commercial and share-alike terms are incompatible
with a CC BY 4.0 release, so catalogue content is never merged into this
project's data. It is used only for benchmarking, and any comparison against it
is pulled at run time rather than stored here.

What this means for the benchmark in `code/benchmark_catalogue.py` is a
division. The per-variant table it produces pairs each variant with the grade
the catalogue assigns it, which is the catalogue's mapping rewritten, so that
table is written to `quarantine/` and is not released. The confusion cells and
the sensitivity, specificity and predictive values derived from them are
measurements of how the catalogue performs on this cohort. No grade for any
variant can be recovered from a confusion matrix over fourteen thousand
isolates, and a rule under which reporting such a measurement required the
catalogue's own licence would make an independent evaluation of that catalogue
unpublishable by anyone. Those measurements are released under CC BY 4.0 with
the catalogue cited. The division is enforced in the code: the per-variant
table must resolve inside `quarantine/`, the released files inside `outputs/`,
and the released report is searched for every variant the module graded before
it is written.

This exclusion extends to the CRyPTIC `EFFECTS` and `PREDICTIONS` tables.
Although CRyPTIC distribute those two tables under CC BY 4.0, their content is
WHO catalogue classifications applied to CRyPTIC samples, and their `EVIDENCE`
field reproduces WHO catalogue text verbatim. They are treated as WHO catalogue
content for licensing purposes.

**Personal data.** The CRyPTIC `BASHTHEBUG_CLASSIFICATIONS` table contains
volunteer usernames, user identifiers and IP addresses. It is not used, not
stored in this repository, and no derived table retains any of those fields.

## Attribution

When reusing material from this project, cite the Zenodo record for the release
used. When reusing the underlying phenotype and genotype data, cite the CRyPTIC
Consortium Dataset version DOI above, as CC BY 4.0 requires.
