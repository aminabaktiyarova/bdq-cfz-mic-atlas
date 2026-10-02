# Data dictionary

Every column of every table this project releases, with the definition behind
it. Tables are written to `outputs/` by the modules named against each one and
are regenerated from the CRyPTIC release recorded in `docs/PROVENANCE.md`.

Derived data is released under CC BY 4.0. No table here carries WHO catalogue
content, CRyPTIC's EFFECTS or PREDICTIONS tables, or any field from
BASHTHEBUG_CLASSIFICATIONS. See `LICENSES.md`.

Concentrations are in mg/L. A mean or a bound written as log2 is the base-2
logarithm of a concentration in mg/L; a shift is in doublings of MIC, so a
shift of 1.0 is a doubling.

---

## outputs/atlas_evidence.csv

Written by `code/build_atlas.py`. One row per distinct mutation.

A row counts only solo isolates: samples carrying exactly one real major-allele
variant across Rv0678, pepQ and atpE, with no null call, het call or minor
allele in any of the three. A sample carrying two variants cannot attribute its
MIC to either, and a sample whose gene could not be read can be asserted neither
to carry the variant nor to lack it. mmpL5 is never the subject of a row.

| Column | Type | Definition |
| --- | --- | --- |
| `gene` | text | Rv0678, pepQ or atpE |
| `mutation` | text | The variant, in the GARC grammar the CRyPTIC tables use |
| `class` | text | gene deletion, frameshift, stop codon, in-frame indel, promoter, or substitution |
| `isolates` | integer | Solo isolates carrying this variant that have both a genome and a UKMYC MIC for at least one of the two drugs |
| `clusters` | integer | Distinct combinations of site, sublineage and this mutation among those isolates. The independent evidence behind the row |
| `sites` | integer | Distinct collection sites among those isolates |
| `sublineages` | integer | Distinct sublineages among those isolates |
| `largest_cluster` | integer | Isolates in the largest single cluster. A number close to `isolates` means the evidence is one outbreak |
| `mmpL5_disrupted` | integer | Isolates whose mmpL5 carries a frameshift, stop codon or gene deletion. Carried as a covariate, never as the subject |
| `BDQ_isolates` | integer | Of `isolates`, those carrying a bedaquiline MIC |
| `BDQ_resistant` | integer | Of `BDQ_isolates`, those above the bedaquiline ECOFF |
| `BDQ_left_censored` | integer | Of `BDQ_isolates`, those reported at or below the lowest tested concentration |
| `BDQ_right_censored` | integer | Of `BDQ_isolates`, those reported above the highest tested concentration |
| `CFZ_isolates` | integer | As `BDQ_isolates`, for clofazimine |
| `CFZ_resistant` | integer | As `BDQ_resistant`, for clofazimine |
| `CFZ_left_censored` | integer | As `BDQ_left_censored`, for clofazimine |
| `CFZ_right_censored` | integer | As `BDQ_right_censored`, for clofazimine |
| `BDQ_reference_mean` | log2 mg/L | Fitted mean log2 bedaquiline MIC of the reference group, identical in every row |
| `BDQ_mean` | log2 mg/L | Fitted mean log2 bedaquiline MIC for this variant. Empty where no estimate is made |
| `BDQ_shift` | doublings | `BDQ_mean` minus `BDQ_reference_mean`. Positive means a higher MIC than the reference group |
| `BDQ_shift_low` | doublings | Lower bound of the 95% interval on the shift, from resampling clusters |
| `BDQ_shift_high` | doublings | Upper bound of the same interval |
| `BDQ_not_estimated` | text | Why the bedaquiline shift is empty. Blank where a shift is given |
| `CFZ_reference_mean` | log2 mg/L | As `BDQ_reference_mean`, for clofazimine |
| `CFZ_mean` | log2 mg/L | As `BDQ_mean`, for clofazimine |
| `CFZ_shift` | doublings | As `BDQ_shift`, for clofazimine |
| `CFZ_shift_low` | doublings | As `BDQ_shift_low`, for clofazimine |
| `CFZ_shift_high` | doublings | As `BDQ_shift_high`, for clofazimine |
| `CFZ_not_estimated` | text | As `BDQ_not_estimated`, for clofazimine |

### How an estimate is withheld

A shift is given only where the evidence supports one, and `*_not_estimated`
records the reason otherwise. Two rules withhold it.

Fewer than five clusters. Resampling k clusters draws from as many distinct
values as there are clusters, so at three clusters the interval describes the
resampling rather than the uncertainty.

A fitted mean more than one doubling outside the tested concentration range.
Where most of a variant's isolates are censored at the plate floor, the fit
answers with a mean far below anything measured, and that number is
extrapolation rather than an estimate.

### Reading a row

`clusters`, not `isolates`, is the evidence. A variant carried by fifty
isolates in two clusters has been observed twice. `largest_cluster` beside
`isolates` shows how much of the count is one event.

A shift is against the reference group of the same drug, which is every sample
with no real major-allele variant in the three genes and nothing uncertain in
them. The interval covers the variant's own resampling; the reference group
holds about fourteen thousand measurements and its mean carries a standard
error near 0.01 doublings.

Resistance in this collection is concentrated at one site by design, so no
count in this table is a prevalence estimate.

---

## outputs/mic_estimates.csv

Written by `code/mic_model.py`. One row per variant class per drug, for classes
reaching the module's minimum group size.

| Column | Type | Definition |
| --- | --- | --- |
| `drug` | text | BDQ or CFZ |
| `group` | text | The variant class, as defined in `code/cohort.py` |
| `isolates` | integer | Isolates in the group carrying a placeable MIC for this drug |
| `clusters` | integer | Distinct combinations of site, sublineage and defining mutation among them |
| `fitted mean log2` | log2 mg/L | Maximum likelihood mean of the group's log2 MIC, over the censoring intervals |
| `sd` | doublings | Maximum likelihood standard deviation of the same distribution |
| `shift (doublings)` | doublings | The group's fitted mean minus the reference group's |
| `fold change` | ratio | Two raised to the shift |
| `95% CI on mean` | text | Interval on the group's fitted mean, from resampling clusters, written as two log2 values |

The interval is on the group mean rather than on the shift. The reference
group's own mean carries a standard error near 0.01 doublings.

---

## outputs/unexplained_counts.csv

Written by `code/unexplained.py`. One row per drug, phenotype-quality subset and
tier. The four tiers partition the resistant isolates.

| Column | Type | Definition |
| --- | --- | --- |
| `drug` | text | BDQ or CFZ |
| `quality` | text | `all` phenotype qualities, or `HIGH` only |
| `tier` | text | attributable; carrier, no demonstrated effect; unexplained; indeterminate |
| `isolates` | integer | Resistant isolates in this tier |
| `resistant_total` | integer | Resistant isolates in the subset, the denominator |
| `percent` | percent | `isolates` over `resistant_total` |
| `wilson_low` | percent | Lower bound of the 95% Wilson interval on that percentage |
| `wilson_high` | percent | Upper bound of the same interval |

The intervals describe sampling error inside this collection. They are not
population estimates.

---

## outputs/unexplained_gene_sets.csv

Written by `code/unexplained.py`. The same count under three gene sets, which
shows how much of the figure is the choice of genes.

| Column | Type | Definition |
| --- | --- | --- |
| `drug` | text | BDQ or CFZ |
| `quality` | text | `all` phenotype qualities, or `HIGH` only |
| `genes` | text | The genes admitted, space separated |
| `gene_set` | text | The same set, named as the report prints it |
| `resistant_total` | integer | Resistant isolates in the subset |
| `unexplained` | integer | Of those, carrying no real major-allele variant in the admitted genes |
| `percent` | percent | `unexplained` over `resistant_total` |
| `indeterminate` | integer | Of those, whose admitted genes carry a null, het or minor call |
| `carrier` | integer | Of those, carrying a real major-allele variant in the admitted genes |

---

## outputs/heteroresistance_estimates.csv

Written by `code/heteroresistance.py`. Four kinds of estimate about the samples
the cohort definitions exclude as uncertain, stacked with a column naming which.

| Column | Type | Definition |
| --- | --- | --- |
| `drug` | text | BDQ or CFZ |
| `estimate` | text | `shift` for a per-group MIC shift against the reference group, `slope` for the read-fraction regression, `joint shift` and `site-adjusted shift` for the two fits that hold every group in one regression with the reference group |
| `group` | text | `minor` groups carry one detected minor allele at Rv0678, resolved to a class; `major` groups carry that class as a major allele, on the same cohort and filters |
| `isolates` | integer | Isolates in the group carrying a placeable MIC |
| `clusters` | integer | Distinct combinations of site, sublineage and the defining or resolved mutation |
| `resistant` | integer | Of those, above the ECOFF. Empty on a slope row |
| `reference_mean` | log2 mg/L | Fitted mean log2 MIC of the reference group. On a `site-adjusted shift` row, that mean at the baseline site, which is the first site label in sorted order. Empty on a slope row |
| `shift` | doublings | The group's fitted mean minus the reference mean. Empty on a slope row or where withheld |
| `shift_low` | doublings | Lower bound of the 95% interval on the shift, from resampling clusters |
| `shift_high` | doublings | Upper bound of the same interval |
| `slope` | doublings per unit | Change in fitted mean log2 MIC per unit of read fraction. Empty on a shift row or where withheld |
| `slope_low` | doublings per unit | Lower bound of the 95% interval on the slope, from resampling clusters |
| `slope_high` | doublings per unit | Upper bound of the same interval |
| `withheld` | text | Why an estimate is empty. Blank where one is given |

The read fraction is FRS, the fraction of reads supporting the minor allele, and
runs from 0.077 to 0.890 in this cohort. Every estimate in this table is
exploratory.

A `shift` row fits each group against a reference group fitted on its own, so
each carries its own standard deviation. A `joint shift` row fits every group
and the reference group in one censored regression under one standard
deviation, and a `site-adjusted shift` row adds an indicator per collection
site to that same regression. The two joint rows are fitted on the same rows,
so the distance between them is the site adjustment. A site with fewer than 30
isolates is pooled into one `other` level, and a group with fewer than 12 is
dropped from the joint fits rather than left in them as reference.

---

## outputs/prediction_metrics.csv

Written by `code/prediction_metrics.py`. One row per rule per drug.

| Column | Type | Definition |
| --- | --- | --- |
| `drug` | text | BDQ or CFZ |
| `rule` | text | major variant; Rv0678 any; Rv0678 loss; major or minor |
| `isolates` | integer | Isolates carrying an MIC for this drug, the denominator |
| `clusters` | integer | Distinct clusters among them |
| `resistant` | integer | Of those, above the ECOFF |
| `true_positive` | integer | Called by the rule and resistant |
| `false_positive` | integer | Called by the rule and not resistant |
| `false_negative` | integer | Not called and resistant |
| `true_negative` | integer | Not called and not resistant |
| `sensitivity` | fraction | True positives over resistant isolates |
| `sensitivity_low` | fraction | Lower bound of the 95% interval, from resampling clusters |
| `sensitivity_high` | fraction | Upper bound of the same interval |
| `specificity` | fraction | True negatives over isolates that are not resistant |
| `specificity_low` | fraction | Lower bound of the 95% interval |
| `specificity_high` | fraction | Upper bound of the same interval |
| `ppv` | fraction | True positives over isolates the rule calls |
| `ppv_low` | fraction | Lower bound of the 95% interval |
| `ppv_high` | fraction | Upper bound of the same interval |
| `npv` | fraction | True negatives over isolates the rule does not call |
| `npv_low` | fraction | Lower bound of the 95% interval |
| `npv_high` | fraction | Upper bound of the same interval |

A rule calls a sample whose gene could not be read as not resistant, because a
null call is not a detected variant. The predictive values depend on how much
resistance the collection holds, which here is concentrated at one site by
design, so they transfer to no other collection.

---

## outputs/gene_vocabulary.csv

Written by `code/inspect_mutations.py`. Every gene appearing in the CRyPTIC
MUTATIONS table with the number of rows it holds. A description of the source
table rather than a result of this project.

| Column | Type | Definition |
| --- | --- | --- |
| `GENE` | text | The gene name as the source table spells it |
| `ROWS` | integer | Rows in MUTATIONS naming that gene, across every sample |
