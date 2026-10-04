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

A large deletion is written twice in the source table, once as the fraction of
the gene absent and once as the sequence removed, and the two rows are one
event. The fraction row is the variant here and the sequence row is dropped, so
a sample whose only finding is such a deletion counts as a solo isolate.
`docs/PROVENANCE.md`, schema question 9, has the counts.

The class of a variant is decided in the order the column lists: a gene deletion
first, then a frameshift, then a stop codon, then a change confined to the
promoter, then an indel that leaves the reading frame intact, and a substitution
otherwise. A promoter has no reading frame, so an insertion or deletion inside
one is classed as a promoter change rather than as an in-frame indel. A deletion
written from a promoter position reaches the coding sequence when it is long
enough, and then the bases it removes from the gene decide whether it shifts the
frame.

| Column | Type | Definition |
| --- | --- | --- |
| `gene` | text | Rv0678, pepQ or atpE |
| `mutation` | text | The variant, in the GARC grammar the CRyPTIC tables use. A gene deletion is written `del_` and the fraction of the gene absent, so `del_0.81` is 81 per cent deleted |
| `class` | text | gene deletion, frameshift, stop codon, promoter, in-frame indel, or substitution, decided in that order |
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

## outputs/multi_allele_counts.csv

Written by `code/heteroresistance.py`. Resistance in the samples carrying more
than one detected minor allele at Rv0678 and nothing else in the three genes.
One row per drug and subset.

| Column | Type | Definition |
| --- | --- | --- |
| `drug` | text | BDQ or CFZ |
| `subset` | text | `all`, `two alleles`, or `three or more`. The last two partition the first |
| `isolates` | integer | Samples in the subset carrying a placeable MIC for this drug |
| `clusters` | integer | Distinct combinations of site, sublineage and the defining mutation. Equal to `isolates` throughout, since these samples carry no defining major-allele mutation to share |
| `resistant` | integer | Of those, above the ECOFF |
| `percent` | percent | `resistant` over `isolates` |
| `reference_percent` | percent | The same rate in the reference group, over the reference samples carrying a placeable MIC for this drug, which is 14,013 for bedaquiline and 14,038 for clofazimine rather than all 14,187 |

A sample carrying two minor alleles at Rv0678 cannot attribute its MIC to either,
which is why `heteroresistance_estimates.csv` excludes it from the per-class
estimates. It is still a sample with no wild-type assertion available at the gene
and no major allele to explain it, so the group is counted here and carries its
own column in the joint and site-adjusted fits, under the group name
`minor, two or more alleles`.

The MIC shift for the group is in `heteroresistance_estimates.csv`, not here.

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

## outputs/prediction_thresholds.csv

Written by `code/prediction_metrics.py`. The same four rules as
`prediction_metrics.csv`, evaluated at every concentration both plate designs
tested rather than at the ECOFF alone. One row per drug, rule and cut-off.

| Column | Type | Definition |
| --- | --- | --- |
| `drug` | text | BDQ or CFZ |
| `rule` | text | The genotype rule, as in `prediction_metrics.csv` |
| `threshold_mg_L` | mg/L | The cut-off, a concentration tested on every plate design for this drug |
| `threshold_log2` | log2 mg/L | The same cut-off in log2 units |
| `is_ecoff` | boolean | True on the cut-off the ECOFF sits at, whose row reproduces `prediction_metrics.csv` |
| `isolates` | integer | Isolates whose position relative to the cut-off is determinate |
| `indeterminate` | integer | Isolates the plate cannot place either side of the cut-off. Zero at every cut-off in this table, by the choice of cut-offs |
| `clusters` | integer | Distinct combinations of site, sublineage and the defining mutation among those isolates |
| `above_threshold` | integer | Of those, whose MIC interval lies above the cut-off |
| `true_positive` | integer | Called by the rule and above the cut-off |
| `false_positive` | integer | Called by the rule and at or below the cut-off |
| `false_negative` | integer | Not called and above the cut-off |
| `true_negative` | integer | Not called and at or below the cut-off |
| `sensitivity` | fraction | True positives over isolates above the cut-off |
| `sensitivity_low` | fraction | Lower bound of the 95% interval, from resampling clusters |
| `sensitivity_high` | fraction | Upper bound of the same interval |
| `specificity` | fraction | True negatives over isolates at or below the cut-off |
| `specificity_low` | fraction | Lower bound of the 95% interval |
| `specificity_high` | fraction | Upper bound of the same interval |
| `ppv` | fraction | True positives over isolates the rule calls |
| `ppv_low` | fraction | Lower bound of the 95% interval |
| `ppv_high` | fraction | Upper bound of the same interval |
| `npv` | fraction | True negatives over isolates the rule does not call |
| `npv_low` | fraction | Lower bound of the 95% interval |
| `npv_high` | fraction | Upper bound of the same interval |

A cut-off is usable only where every isolate's position relative to it is
determinate. Below the highest of the designs' lowest rungs a left-censored
reading sits on neither side, and above the lowest of the designs' highest rungs
a right-censored reading sits on neither side. A concentration that is a rung on
every design is inside both bounds, which is why the cut-offs are the shared
rungs: 0.015 to 1 mg/L for bedaquiline, seven of them, and 0.06 to 2 mg/L for
clofazimine, six.

Position is read from the censoring interval rather than the reported number. An
interval lies above the cut-off when its lower bound reaches it and at or below
when its upper bound does not exceed it. On the real cohort that rule and
CRyPTIC's own `BINARY_PHENOTYPE` agree at the ECOFF on all 14,972 bedaquiline and
14,997 clofazimine isolates, with no disagreement and nothing indeterminate.

Predictive values depend on how much resistance the collection holds, and
resistance here is concentrated at one site by design, so they describe this
collection and transfer to no other.

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

---

## outputs/wider_gene_set.csv

Written by `code/wider_gene_set.py`. One row per candidate gene per drug,
twenty rows. Every row is computed inside the reference group, which is the
unexplained tier of a resistant isolate, so the comparison is resistant against
susceptible among isolates carrying no real major-allele variant in Rv0678,
pepQ or atpE and nothing uncertain in them.

| Column | Type | Definition |
| --- | --- | --- |
| `drug` | text | BDQ or CFZ |
| `gene` | text | The candidate gene |
| `isolates` | integer | Reference-group isolates carrying an MIC for this drug |
| `resistant` | integer | Of those, above the ECOFF |
| `carriers` | integer | Of those, carrying a real major-allele variant in the gene |
| `carriers_resistant` | integer | Carriers that are resistant |
| `carriers_susceptible` | integer | Carriers that are not resistant |
| `clusters` | integer | Distinct clusters, keyed on site, sublineage and the carried variant |
| `odds_ratio` | number | Carrying a variant, resistant against susceptible. Zero where no carrier is resistant. Empty where the estimate is not finite |
| `p_value` | number | Fisher exact, two sided |
| `q_value` | number | Benjamini-Hochberg over the family of twenty tests |
| `or_low` | number | Lower bound of the 95% interval, from resampling clusters |
| `or_high` | number | Upper bound of the same interval |
| `site_or` | number | Mantel-Haenszel odds ratio holding site constant. Empty where fewer than two strata carry information |
| `site_low` | number | Lower bound of its 95% interval |
| `site_high` | number | Upper bound of its 95% interval |
| `site_strata` | integer | Sites contributing a table with no empty margin |
| `site_homogeneity_p` | number | Test of equal odds across those sites. Below 0.05 the pooled estimate summarises strata that disagree |
| `lineage_or` | number | The same estimate holding lineage constant |
| `lineage_low` | number | Lower bound of its 95% interval |
| `lineage_high` | number | Upper bound of its 95% interval |
| `lineage_strata` | integer | Lineages contributing a table with no empty margin |
| `lineage_homogeneity_p` | number | Test of equal odds across those lineages |
| `resistant_carrier_groups` | integer | Site and sublineage groups the resistant carriers fall into, a floor on the number of independent events |
| `largest_group` | integer | Resistant carriers in the largest of those groups |
| `commonest_variant` | text | The variant most of the resistant carriers hold. Empty where none is resistant |
| `commonest_variant_count` | integer | Resistant carriers holding it |

The gene set is the union of the genes TB-Profiler's database associates with
the two drugs, less the three in scope, the modifier gene, and `mmpR5`, which
is that database's alias for Rv0678. Only the gene names are taken from it: no
grade or confidence is read, so this table carries no catalogue content.

An odds ratio here is an association inside one collection whose resistance is
concentrated at one site by design, and the genes are polymorphic, so a crude
figure is read against the three columns beside it: the interval that resamples
clusters, the homogeneity test, and the number of groups the resistant carriers
fall into.

---

## outputs/benchmark_metrics.csv

Written by `code/benchmark_catalogue.py`. One row per catalogue variant per
drug, eight rows: two conversions of the WHO catalogue, second edition, each
with and without the epistasis rule.

These are measurements of how an external catalogue performs on this cohort.
They carry no grade for any variant, and no grade can be recovered from them.
The per-variant table the same module produces is the catalogue's own mapping
and is not released; see `LICENSES.md`.

| Column | Type | Definition |
| --- | --- | --- |
| `catalogue` | text | The conversion, and whether the epistasis rule was applied |
| `drug` | text | BDQ or CFZ |
| `isolates` | integer | Isolates carrying an MIC for this drug, the denominator |
| `clusters` | integer | Distinct clusters among them |
| `resistant` | integer | Of those, above the ECOFF |
| `true_positive` | integer | Graded R by the catalogue and resistant |
| `false_positive` | integer | Graded R and not resistant |
| `false_negative` | integer | Not graded R and resistant |
| `true_negative` | integer | Not graded R and not resistant |
| `sensitivity` | fraction | True positives over resistant isolates |
| `sensitivity_low` | fraction | Lower bound of the 95% interval, from resampling clusters |
| `sensitivity_high` | fraction | Upper bound of the same interval |
| `specificity` | fraction | True negatives over isolates that are not resistant |
| `specificity_low` | fraction | Lower bound of the 95% interval |
| `specificity_high` | fraction | Upper bound of the same interval |
| `ppv` | fraction | True positives over isolates graded R |
| `ppv_low` | fraction | Lower bound of the 95% interval |
| `ppv_high` | fraction | Upper bound of the same interval |
| `npv` | fraction | True negatives over isolates not graded R |
| `npv_low` | fraction | Lower bound of the 95% interval |
| `npv_high` | fraction | Upper bound of the same interval |

An isolate is graded R when any variant it carries is, U when any is U and none
is R, and so on down the precedence R, U, F, S. A grade of U is not a call of
resistance, so these cells count R against everything else. Predictive values
depend on how much resistance the collection holds, which here is concentrated
at one site by design.
