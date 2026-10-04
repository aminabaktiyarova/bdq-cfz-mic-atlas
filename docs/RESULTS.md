# Methods and results

Amina Baktiyarova, Independent Researcher, ORCID 0009-0007-6265-6493

Built on the CRyPTIC Consortium Dataset v3.4.0, Zenodo version DOI
10.5281/zenodo.15680920, CC BY 4.0. `docs/PROVENANCE.md` records the dataset, its
licence, the schema and the properties of the input tables. This document records
what was computed from them: the cohort, the estimator, the effects, and what the
evidence does and does not support. Every column of every released table is
defined in `docs/DATA_DICTIONARY.md`.

## Purpose and scope

The subject is the Rv0678 and MmpL5 efflux axis and its effect on bedaquiline and
clofazimine minimum inhibitory concentrations. The question is quantitative: how
far does a given mutation or mutation class move the MIC, with what evidence
behind the estimate, and how do the two drugs differ. Two drugs are in scope.
Delamanid and linezolid were considered and dropped, for reasons given under
"Drugs considered and dropped".

No resistance catalogue takes part in any definition, cohort, or figure here.
Every genotype in this document is read from the mutation strings in
`MUTATIONS.parquet`. The CRyPTIC `EFFECTS` and `PREDICTIONS` tables, which carry
WHO catalogue classifications applied to CRyPTIC samples, and the
`GENOMES.ANTIBIOGRAM` column, which is derived the same way, are excluded from
every figure in this repository.

One comparison in this work was pre-specified and tested against a held-out half
of the cohort, recorded in `docs/PRE_REGISTRATION.md` and under "The
pre-registered test" below. Everything else is an estimate with an interval, not
a tested hypothesis. "Limitations" states the constraints that apply to all of
it.

## Cohort construction

### Definitions

Every definition below is read from a mutation string.

**Real variant.** A called, non-synonymous change: a substitution, an indel, or a
gene deletion. A null call, a het call and a synonymous substitution are none of
these.

**Major allele.** `IS_MINOR` false.

**Uncertain.** The gene carries a null call, a het call, or any minor allele. Such
a sample cannot be asserted to carry a variant at that gene nor to lack one, so
it enters neither the carrier groups nor the reference group. What these samples
show is estimated separately, under "Minor alleles and null calls".

**Variant class.** Gene deletion, frameshift, stop codon, in-frame indel,
promoter, or substitution. The classes are kept apart rather than pooled, because
they do not behave alike.

**Loss of function.** Frameshift, stop codon, or gene deletion.

**Solo.** Exactly one real major-allele variant across Rv0678, pepQ and atpE,
with nothing uncertain in any of them.

**Reference.** No real major-allele variant in Rv0678, pepQ or atpE, with nothing
uncertain in any of them.

mmpL5 enters neither the solo criterion nor the reference criterion, and is
carried as a covariate instead. It holds a real major-allele variant in 53,361 of
54,057 genomes, so requiring a clean mmpL5 would leave a reference group of 109
rather than 49,971. Its role in this axis is to modify the effect of a repressor
variant rather than to produce resistance on its own.

### Classification

Across all 54,057 genomes:

| Category | Genomes |
| --- | --- |
| Reference | 49,971 |
| Solo, exactly one variant | 2,588 |
| Two or more variants | 75 |
| Excluded as uncertain | 1,423 |

51,305 genomes carry no real major-allele variant in the three genes. 1,334 of
them are excluded because a gene carries a null, het or minor call, leaving the
reference group of 49,971. Uncertainty falls on Rv0678 in 836 genomes, pepQ in
598 and atpE in 61. mmpL5 carries a loss-of-function variant in 622.

15,158 genomes also carry a UKMYC bedaquiline or clofazimine MIC, and those are
the analysis cohort:

| Group | Isolates |
| --- | --- |
| reference | 14,187 |
| uncertain | 289 |
| pepQ solo | 217 |
| Rv0678 substitution | 195 |
| Rv0678 frameshift | 148 |
| Rv0678 promoter | 83 |
| multiple variants | 14 |
| atpE solo | 10 |
| Rv0678 gene deletion | 6 |
| Rv0678 stop codon | 5 |
| Rv0678 in-frame indel | 4 |

Real major-allele variants across the full 54,057, by gene and class:

| Gene | Frameshift | Gene deletion | In-frame indel | Promoter | Stop codon | Substitution |
| --- | --- | --- | --- | --- | --- | --- |
| Rv0678 | 519 | 13 | 11 | 654 | 32 | 756 |
| atpE | 0 | 0 | 0 | 12 | 0 | 54 |
| mmpL5 | 286 | 2 | 42 | 0 | 336 | 96,337 |
| pepQ | 15 | 0 | 2 | 11 | 1 | 754 |

atpE is descriptive rather than inferential in this cohort: 10 solo isolates with
an MIC. The Rv0678 gene-deletion, stop-codon and in-frame-indel groups hold 6, 5
and 4 isolates, and none of the three carries a fitted estimate anywhere in this
work.

A deletion of most of the gene is recorded twice in v3.4.0, once as a fraction of
the gene and once as the sequence removed, and the two rows describe one event.
The counts above take the event once, which is why a sample whose only finding is
such a deletion is a solo sample here rather than one carrying two variants.

### Isolates from one patient

The 15,158 isolates come from 14,708 patients, 1.031 isolates per patient. 309
patients contribute more than one isolate, accounting for 759 isolates, 5.0% of
the cohort: 214 patients with two, 54 with three, 36 with four, 5 with five.

Restricting to one isolate per patient leaves every effect where it was:

| Comparison | Drug | All isolates | One per patient |
| --- | --- | --- | --- |
| Rv0678 loss of function | BDQ | OR 50.9 | OR 46.9 |
| Rv0678 substitution | BDQ | OR 23.6 | OR 21.4 |
| Rv0678 loss of function | CFZ | OR 13.1 | OR 12.5 |
| Rv0678 substitution | CFZ | OR 5.9 | OR 5.7 |

The loss-of-function rows are restricted to isolates whose mmpL5 is intact, 116
of the 159 carriers, because a loss-of-function variant in mmpL5 removes the
efflux pump the repressor controls. The substitution rows need no such
restriction: none of the 195 substitution carriers has an mmpL5 loss-of-function
variant. Every odds ratio here is against the reference group by Fisher's exact
test.

Clonal relatedness, which is a larger effect than patient-level replication, is
treated under "Clonal clustering".

### Phenotype quality

Restricted to MICs read by at least two agreeing methods, three of the four
effects strengthen:

| Comparison | Drug | All quality | HIGH only |
| --- | --- | --- | --- |
| Rv0678 loss of function | BDQ | OR 50.9 (n=116) | OR 74.6 (n=96) |
| Rv0678 substitution | BDQ | OR 23.6 (n=195) | OR 29.3 (n=157) |
| Rv0678 loss of function | CFZ | OR 13.1 (n=116) | OR 11.7 (n=71) |
| Rv0678 substitution | CFZ | OR 5.9 (n=195) | OR 7.6 (n=142) |

The same mmpL5 restriction applies to the loss-of-function rows. No effect rests
on the disputed readings, so the high-quality subset was used as the primary
analysis for the pre-registered work.

### Isolates carrying no MIC

186 bedaquiline and 161 clofazimine rows in the cohort carry no MIC. Each also
carries a null `BINARY_PHENOTYPE`, so none enters the resistant counts, and none
is a value that failed to match the tested concentration series. They are
excluded from every fit, which takes the reference group from 14,187 to 14,013
for bedaquiline and 14,038 for clofazimine.

Counts below are over the rows that carry a plate design, 15,156 for bedaquiline
and 15,158 for clofazimine; two isolates have no bedaquiline row at all.

The loss is concentrated rather than random. Site 11 loses 86 of its 443
bedaquiline rows, 19.41%, and 87 of 443 clofazimine rows, 19.64%, against 98 of
14,713 and 74 of 14,715 elsewhere, 0.67% and 0.50%. Across all 14 drugs, 86
samples are missing every one, which is that block.

Plate design carries a second association. UKMYC5 loses 135 of 6,346 bedaquiline
rows, 2.13%, against 49 of 8,810 on UKMYC6, 0.56%; for clofazimine 131 of 6,348,
2.06%, against 30 of 8,810, 0.34%. Eight sites ran both designs and six of them
lost at least one row, and over those six the Mantel-Haenszel odds ratio for
UKMYC5 against UKMYC6 is 4.89 (3.64 to 6.57) for bedaquiline and 8.34 (5.82 to
11.96) for clofazimine, so the design association is not the site association
restated.

An apparent association between an uncallable gene and a missing MIC does not
survive the same adjustment. Crude, for bedaquiline, 4.15% of the uncertain group
is missing against 1.21% of the reference group and 0.00% of carriers,
chi-square 29.2, p = 4.6e-07. Excluding site 11 gives 0.77%, 0.70% and 0.00%,
chi-square 4.8, p = 0.091, and clofazimine moves from p = 5.35e-05 to p = 0.366.

The exclusion cannot move a fitted mean far, because a missing isolate acts only
through the value it would have had, weighted 174 of 14,187 for bedaquiline and
149 of 14,187 for clofazimine. The widest gap between a site's fitted reference
mean and the whole group's, over the sites carrying at least 100 placeable
intervals, is 0.931 doublings for bedaquiline and 1.096 for clofazimine, both at
site 14. A missing isolate differing from the fitted mean by that much would move
the whole-group mean by 0.011 and 0.012. At an implausible gap of 3 doublings the
movement is 0.037 and 0.032, below the precision the means are reported to.

## The MIC model

### Why a median does not answer the question

A broth microdilution plate tests a doubling series, so every MIC it reports is
an interval rather than a value. An MIC of 0.25 means growth was inhibited at
0.25 and not at 0.12, so the true value lies in (0.12, 0.25]. At the ends the
interval is unbounded: a reported `<=0.008` means the MIC is somewhere below the
lowest well, and a reported `>1` that it is above the highest. Of the isolates
carrying a bedaquiline MIC, 13.89% are left-censored and 27 are right-censored; of
those carrying a clofazimine MIC, 40.43% and 25.

A median of the reported values answers a different question, because it returns
the well where growth stopped rather than the concentration where it would have.
On the reference group the median sits 0.25 doublings above the fitted mean for
bedaquiline and 0.34 above it for clofazimine, and the gap widens with censoring.

### The estimator

Within a group, log2 MIC is taken as normally distributed and the two parameters
are estimated by maximum likelihood over the censoring intervals. One
observation contributes

    Phi((upper - mu) / sigma) - Phi((lower - mu) / sigma)

with `lower` minus infinity on a left-censored value and `upper` plus infinity on
a right-censored one.

That contribution is evaluated in log space throughout, with `logcdf`, `logsf`
and `log1p`, and never as the difference of two cumulative distribution
functions. Far out in a tail the two functions agree to most of their digits and
subtracting them destroys the rest. The bedaquiline reference group holds 17
intervals carrying a probability mass below 1e-08, the smallest 3.6e-09, which is
enough for the subtraction to corrupt the gradient and stop the optimiser short
of the maximum. The clofazimine group holds none below that floor, its smallest
mass 3.8e-07, so the two drugs do not stress the arithmetic equally.

Intervals on a group mean come from resampling clusters with replacement rather
than isolates, for the reason given under "Clonal clustering".

### Validation against planted parameters

Simulating from known distributions and censoring the result onto the real UKMYC6
bedaquiline ladder, 4,000 draws per row:

| True mean | Fitted mean | True sd | Fitted sd | Naive median | Left-censored |
| --- | --- | --- | --- | --- | --- |
| -5.0 | -5.03 | 1.0 | 1.00 | -4.06 | 2.5% |
| -6.5 | -6.50 | 1.0 | 1.01 | -6.06 | 31.9% |
| -7.5 | -7.52 | 1.2 | 1.18 | -6.97 | 67.9% |
| -2.0 | -1.98 | 1.5 | 1.48 | -1.00 | 0.0% |
| -4.0 | -4.01 | 0.8 | 0.79 | -3.06 | 0.0% |

Every mean is recovered within 0.03 doublings and every standard deviation
within 0.02. The naive median is biased upward by half a doubling to a full
doubling throughout, and at 67.9% censoring it sticks near the plate floor. The
table is section 2 of `outputs/mic_model_report.txt`.

### The reference distributions

Bedaquiline: 14,013 isolates, fitted mean log2 MIC -5.31, standard deviation
1.09, corresponding to an MIC of 0.0253 mg/L. Clofazimine: 14,038 isolates,
fitted mean -4.40, standard deviation 1.29, MIC 0.0474 mg/L.

Per-group fitted means, standard deviations, shifts, fold changes and intervals
are in `outputs/mic_estimates.csv`. A shift is the difference of two fitted means,
and its second decimal is not stable across Python stacks: several of the
unrounded differences fall within 0.002 of a two-decimal rounding boundary, so a
regeneration on a different stack can move a last digit by 0.01. Shifts are
quoted here to one decimal for that reason, and the interval rather than the last
digit carries the uncertainty.

### A class is not one distribution

The Rv0678 frameshift group has a fitted standard deviation of 2.81 doublings
against 1.09 in the reference group, and its interval spans four doublings. The
group holds a 41-isolate clonal outbreak whose MICs sit near the plate floor
alongside frameshifts with raised MICs, so it is at least bimodal and the class
mean describes neither mode. The substitution group shows the same at 1.72,
consistent with lineage-marker polymorphisms mixed in with resistance mutations.
This is the argument for a per-variant layer rather than stopping at classes.

### What a binary call cannot express

Within lineage4, clofazimine loss of function gives an odds ratio of 1.1 at
p = 0.94, on 2 resistant isolates in 42. The MIC says otherwise: the reference
mean is -4.54, the loss-of-function mean -3.36, a shift of 1.2 doublings whose
interval, -3.83 to -2.97 on the group mean, does not touch the reference. The
effect is there and lands about 1.4 doublings short of the clofazimine ECOFF at
0.25 mg/L, so almost nothing crosses it.

Read across lineages, the clofazimine response to losing Rv0678 runs from 2.2
doublings in lineage3 and 2.1 in lineage2 down to 1.2 in lineage4, where the
site-adjusted odds ratios of 13.4, 10.2 and 1.1 suggested presence against
absence. Lineage-stratified figures are in `outputs/mic_estimates.csv`, and the
lineage contrast is confounded by site, treated under "Site and lineage".

### Concentrations that are not concentrations

Three `CONC` entries in `PLATE_LAYOUT` cannot be read as concentrations: the
positive-control wells for MYCOTB, UKMYC5 and UKMYC6, all recorded as 0. They are
skipped when the tested series is built. No MIC in the cohort fails to place on
the series it was measured against; the isolates without a usable value are the
ones carrying no MIC at all, counted under "Isolates carrying no MIC".

## Variant class effects

### Against the reference group

14,187 isolates carry no real major-allele variant in the three genes and nothing
uncertain in any of them. 80 of them are bedaquiline resistant, 0.56%, and 527 are
clofazimine resistant, 3.71%. Every odds ratio below is that group's comparison by
Fisher's exact test, on all isolates in the group, with the denominator being the
isolates carrying a MIC for the drug.

| Group | Isolates | BDQ resistant | BDQ OR | BDQ p | CFZ resistant | CFZ OR | CFZ p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Rv0678 frameshift | 148 | 25, 16.9% | 35.8 | 1.3e-27 | 35, 23.6% | 8.0 | 3.7e-18 |
| Rv0678 substitution | 195 | 23, 11.8% | 23.6 | 6.1e-22 | 36, 18.5% | 5.9 | 5.2e-15 |
| Rv0678 stop codon | 5 | 1, 20.0% | 44.1 | 0.028 | 2, 40.0% | 17.3 | 0.013 |
| Rv0678 gene deletion | 6 | 0 | 0.0 | 1 | 2, 33.3% | 13.0 | 0.019 |
| Rv0678 in-frame indel | 4 | 0 | 0.0 | 1 | 0 | 0.0 | 1 |
| Rv0678 promoter | 83 | 1, 1.2% | 2.2 | 0.38 | 2, 2.4% | 0.6 | 0.77 |
| pepQ solo | 217 | 4, 1.8% | 3.3 | 0.038 | 16, 7.4% | 2.1 | 0.010 |
| atpE solo | 10 | 1, 10.0% | 19.6 | 0.056 | 0 | 0.0 | 1 |
| multiple variants | 14 | 1, 7.1% | 13.6 | 0.077 | 3, 21.4% | 7.1 | 0.014 |
| excluded as uncertain | 289 | 35, 12.1% | 24.3 | 3.9e-32 | 60, 20.8% | 6.8 | 2.4e-26 |

Four, five and six isolates support no inference, so the in-frame-indel,
stop-codon and gene-deletion rows are counts rather than estimates. They are
shown because pooling the stop codons and the deletions into a loss-of-function
group is what the next table does, and a reader should see how little they
weigh.

### Promoter variants are not substitutions

Rv0678 promoter variants show no effect on either drug: bedaquiline odds ratio 2.2
at p = 0.38, clofazimine 0.6 at p = 0.77. Pooling them with substitutions, which
an earlier grouping did, dilutes a real effect with a group that has none. On the
MIC scale they move the bedaquiline MIC downward by 1.6 doublings, which is
treated under "The per-variant layer": that shift belongs to one variant rather
than to the class.

### Losing the repressor against changing it

Restricted to isolates whose mmpL5 is intact, against the same reference group:

| Drug | Group | Isolates | Resistant | Odds ratio | p |
| --- | --- | --- | --- | --- | --- |
| BDQ | Rv0678 substitution | 195 | 23, 11.8% | 23.6 | 6.1e-22 |
| BDQ | Rv0678 loss of function | 116 | 26, 22.4% | 50.9 | 5.2e-32 |
| CFZ | Rv0678 substitution | 195 | 36, 18.5% | 5.9 | 5.2e-15 |
| CFZ | Rv0678 loss of function | 116 | 39, 33.6% | 13.1 | 3.7e-26 |

Loss of function against substitution directly, rather than each against the
reference group: odds ratio 2.16 at p = 0.0158 for bedaquiline and 2.24 at
p = 0.0038 for clofazimine. Losing the repressor raises resistance further than
changing it, on both drugs, and clofazimine shows the same ordering on a shallower
gradient. These are unadjusted for site and for clonal relatedness, both of which
are treated below.

### Gene deletions

13 isolates carry an Rv0678 gene deletion and 2 an mmpL5 deletion, all as major
alleles. Deleted fractions run 0.52 to 0.98 with a median of 0.81, and they
cluster by sublineage: 0.81 across six lineage4.10 isolates, 0.82 in two
consecutive accessions of lineage4.3.4.2.1, and 0.52 to 0.55 in three lineage2.2.1
isolates. Identical fractions within a sublineage and different fractions between
them is what clonal inheritance of a real deletion looks like; a coverage artefact
would give ragged fractions. Coverage and depth of the deleted isolates are close
to everything else, median 99.27% against 99.30% and 71.3 against 90.0.

Six of the 13 carry an MIC, and those six are the Rv0678 gene deletion group in
the tables above. None is bedaquiline resistant. Measured against the other 110
isolates of the loss-of-function group, 26 resistant at 23.6%, the probability of
seeing zero resistant among six is 0.20 and Fisher's exact test gives p = 0.34.
Against the frameshift group alone, 25 of 148 at 16.9%, the figures are 0.33 and
p = 0.59. On bedaquiline the deletions are underpowered rather than anomalous.

On clofazimine two of the six are resistant, 33.3% against 33.6% in the other
110, with an odds ratio of 0.99 and p = 1.0. There the deletions sit on the rate
of the group they belong to.

## Site and lineage

Resistance in this collection is concentrated by design. One site, the South
African NICD, supplies 61 of the 171 bedaquiline-resistant and 264 of the 683
clofazimine-resistant isolates, because roughly 1,100 of its samples were added
specifically for being enriched in bedaquiline resistance. Site therefore enters
every stratified estimate below, and no figure anywhere in this work is a
prevalence estimate for any population.

### The most stable effect

Clofazimine substitution against the reference group, site held constant:
Mantel-Haenszel odds ratio 6.54, 95% interval 4.36 to 9.79, across 11 sites, with
a homogeneity p of 0.65. The strata agree, so the pooled figure summarises one
thing rather than several.

### The lineage contrast

Rv0678 loss of function against the reference group, within lineage, unadjusted:

| Lineage | Carriers | BDQ resistant | BDQ OR | BDQ p | CFZ resistant | CFZ OR | CFZ p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| lineage2 | 62 | 16 | 62.3 | 9.1e-21 | 32 | 20.0 | 1.8e-24 |
| lineage3 | 10 | 7 | 324.6 | 1.6e-12 | 4 | 14.0 | 8.5e-04 |
| lineage4 | 42 | 2 | 8.4 | 0.028 | 2 | 1.8 | 0.33 |

The strata disagree across these three lineages, homogeneity p = 0.00065 for
bedaquiline and 0.0014 for clofazimine, so a single pooled odds ratio would
summarise incompatible things. That is the reason nothing in this work reports
one.

Holding site constant inside each lineage:

| Lineage | BDQ, site held constant | CFZ, site held constant |
| --- | --- | --- |
| lineage2 | 59.3 (23.3 to 150.5) | 10.2 (5.8 to 18.0) |
| lineage3 | 255.1 (46.0 to 1414.5) | 13.4 (3.2 to 56.1) |
| lineage4 | 17.6 (3.5 to 88.0) | 1.1 (0.3 to 4.5) |

For clofazimine, lineage4's interval reaches 4.5 and lineage2's starts at 5.8, so
the two do not overlap and that difference is not explained by site. For
bedaquiline every interval overlaps every other, so the apparent lineage
difference does not survive adjustment and is not reported as a finding. On the
MIC scale the surviving clofazimine difference is 1.3 doublings between the two
lineage means, -2.07 against -3.36, which is treated under "What a binary call
cannot express".

## Clonal clustering

Tuberculosis spreads clonally, so one successful strain can appear dozens of times
and be counted as dozens of observations. A cluster here is the combination of
site, sublineage and the genotype-defining mutation, applied to variant carriers;
reference isolates are treated as independent, because they share no mutation
event. This is a conservative proxy. A transmission analysis would cluster on
genomic distance, which needs the full variant table and more compute than this
project has, and two unrelated patients at one site sharing a sublineage and a
common mutation are merged here, which understates the evidence rather than
inflating it.

15,158 isolates fall into 14,758 clusters.

| Group | Isolates | Clusters | Isolates per cluster | Largest cluster |
| --- | --- | --- | --- | --- |
| reference | 14,187 | 14,187 | 1.00 | 1 |
| Rv0678 frameshift | 148 | 60 | 2.47 | 41 |
| Rv0678 promoter | 83 | 16 | 5.19 | 41 |
| Rv0678 substitution | 195 | 111 | 1.76 | 16 |
| pepQ solo | 217 | 87 | 2.49 | 42 |
| excluded as uncertain | 289 | 275 | 1.05 | 15 |
| multiple variants | 14 | 7 | 2.00 | 8 |
| atpE solo | 10 | 10 | 1.00 | 1 |
| Rv0678 gene deletion | 6 | 2 | 3.00 | 5 |
| Rv0678 stop codon | 5 | 5 | 1.00 | 1 |
| Rv0678 in-frame indel | 4 | 4 | 1.00 | 1 |

The largest clusters carrying an Rv0678 variant:

| Cluster | Isolates | BDQ resistant | CFZ resistant |
| --- | --- | --- | --- |
| site 05, lineage4, 192_ins_g | 41 | 0 | 0 |
| site 10, lineage2.2, c-11a | 41 | 0 | 1 |
| site 20, lineage2.2, c-11a | 19 | 1 | 0 |
| site 06, lineage4.4.1.1, M146T | 16 | 0 | 3 |
| site 05, lineage4.1.2.1, L40V | 14 | 0 | 0 |
| site 10, lineage2.2.1, 141_ins_c | 12 | 6 | 7 |
| site 03, lineage2.2.3, E55D | 11 | 0 | 1 |
| site 05, lineage4.1.1.3, R90C | 9 | 0 | 0 |
| site 10, lineage2.2.1, 138_ins_g | 8 | 1 | 7 |
| site 01, lineage2.2, c-11a | 8 | 0 | 0 |

Three of those clusters carry the promoter variant c-11a, 68 isolates between
them with one clofazimine-resistant isolate among them. Counted across the whole
cohort rather than by cluster, c-11a in lineage2.2 is 98 isolates at 7 sites with
7 clofazimine-resistant and none bedaquiline-resistant. L40V gives 14 isolates
with no resistance on either drug, M146T 18 with 3 clofazimine-resistant, E55D 18
with 2. Several of the commonest Rv0678 variants are lineage markers carrying
little or no phenotypic effect, and a catalogue grading on isolate counts weights
them by the size of the outbreak that carries them.

### The same effects three ways

Isolate level treats every isolate as independent. Cluster-robust fits every
isolate with standard errors that allow correlation within a cluster, which
corrects the uncertainty rather than the point estimate. One isolate per cluster
reduces each carrier cluster to a single representative, repeated 200 times, and
its value is that it shows whether an effect survives when a clone gets one vote.

| Drug | Comparison | Isolates | Clusters | Isolate level | Cluster-robust | Draws above 1 |
| --- | --- | --- | --- | --- | --- | --- |
| BDQ | Rv0678 loss of function | 116 | 65 | 50.9 | 61.1 (31.5 to 118.3) | 100% |
| BDQ | Rv0678 substitution | 195 | 111 | 23.6 | 33.0 (18.5 to 58.6) | 100% |
| CFZ | Rv0678 loss of function | 116 | 65 | 13.1 | 7.7 (4.7 to 12.6) | 100% |
| CFZ | Rv0678 substitution | 195 | 111 | 5.9 | 6.5 (4.1 to 10.2) | 100% |

Every effect survives all three treatments. Where the collapsed estimate exceeds
the isolate-level one, the large clusters are the susceptible ones, so clonality
was biasing those effects downward rather than upward.

The one-per-cluster column reports the fraction of draws above 1 rather than a
median, because at 200 draws the median carries Monte Carlo noise. Across five
seeds the bedaquiline loss-of-function median moves between 57.8 and 63.1 and the
clofazimine one between 6.3 and 6.4, while both are above 1 in every draw of every
seed. For the lineage4 clofazimine comparison below, the median holds at 0.77
across all five seeds and the fraction above 1 moves between 6% and 12%. Section 5
of `outputs/cluster_report.txt` carries the sweep.

The lineage contrast under clustering, clofazimine loss of function:

| Lineage | Clusters | Isolate level | Cluster-robust | Draws above 1 |
| --- | --- | --- | --- | --- |
| lineage2 | 31 | 20.0 | 10.8 (6.0 to 19.4) | 100% |
| lineage3 | 8 | 14.0 | 11.7 (4.0 to 34.4) | 100% |
| lineage4 | 24 | 1.8 | 1.1 (0.3 to 4.1) | 9% |

lineage2 is above 1 in every draw and lineage4 in 9% of them, 6% to 12% across
the five seeds, so the clofazimine lineage difference is not an artefact of clonal
expansion.

The effective sample sizes that every interval in this work rests on are 116
isolates in 65 clusters for loss of function and 195 in 111 for substitution.

## The per-variant layer

`outputs/atlas_evidence.csv` carries one row per distinct mutation: its isolate
count, its cluster count, the sites and sublineages it spans, its resistant,
left-censored and right-censored counts on each drug, how many of its carriers
have a disrupted mmpL5, and a fitted shift with an interval where the independent
evidence supports one. Only solo isolates enter a row, and mmpL5 is never the
subject of a row.

217 distinct variants across 668 solo isolates in 295 clusters: 138 Rv0678
variants on 441 isolates in 198 clusters, 69 pepQ variants on 217 isolates in 87
clusters, and 10 atpE variants on 10 isolates in 10 clusters.

152 of the 217 are carried by a single isolate, and 182 rest on a single cluster.
By class: 160 substitutions, 34 frameshifts, 12 promoter variants, 5 stop codons,
4 in-frame indels and 2 gene deletions.

### What is withheld, and why

A variant resting on fewer than five clusters carries no shift. Resampling
clusters draws from as many distinct values as there are clusters, so at three
clusters the resulting interval describes the resampling rather than the
uncertainty in the estimate.

A variant whose fitted mean falls more than one doubling outside the tested
concentration range carries no shift either, because the fit is then extrapolating
past every measurement that produced it.

A withheld row keeps every count and records which rule withheld it, so the
evidence is still readable.

### The four variants that support an estimate

Shifts in doublings against the reference group, to one decimal. The released
table carries three, and the third is not reproducible across Python stacks for
the reason given under "The reference distributions".

| Variant | Class | Isolates | Clusters | Sites | Sublineages | BDQ shift | 95% | CFZ shift | 95% |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Rv0678 141_ins_c | frameshift | 36 | 14 | 5 | 10 | 2.8 | 2.3 to 3.1 | 2.2 | 1.5 to 2.7 |
| Rv0678 c-11a | promoter | 74 | 7 | 7 | 1 | -1.8 | -2.3 to -1.5 | -0.4 | -0.7 to -0.3 |
| Rv0678 138_ins_g | frameshift | 18 | 7 | 5 | 5 | 2.6 | 1.6 to 3.3 | 2.6 | 1.0 to 3.4 |
| pepQ P69L | substitution | 18 | 5 | 5 | 1 | -0.2 | -0.4 to 0.1 | 0.1 | -0.5 to 0.4 |

Four of 217. That is what this dataset supports at variant resolution, and no
amount of further analysis changes it.

### A variant the second rule catches

Rv0678 192_ins_g clears the cluster rule with 51 solo isolates in 9 clusters, 41
of them in one, and fails the range rule. 39 of its 51 bedaquiline MICs are
left-censored, and its fitted mean is -10.08 against a lowest tested concentration
of -6.97, so the fit sits three doublings below anything the plate could measure.
The clofazimine fit lands at -6.99 against a floor of -5.06. The bootstrap interval
on the bedaquiline shift spans roughly -8 to +3 doublings, and its bounds move by
about 0.3 with the resampling seed, which is what a fit resting on censored
observations and 9 clusters produces. Both shifts are withheld and its counts
stand. What this variant is not is a variant with a measured MIC near the
reference: its MICs are mostly below the plate, and how far below cannot be said
from these data.

### Two decompositions a class layer cannot express

The Rv0678 frameshift class shifts the bedaquiline MIC by 0.8 doublings. Remove
192_ins_g and the 97 isolates that remain shift by 2.3, and the clofazimine figure
moves from 0.7 to 1.8. One clone at one site halves the class estimate, and a
catalogue grading on isolate counts has no way to show it.

The Rv0678 promoter class shift of -1.6 doublings is one variant. The 8 promoter
isolates that are not c-11a shift by +0.5 for bedaquiline and +0.9 for
clofazimine. Eight isolates support no claim of their own, so what this
establishes is the negative: the class figure describes c-11a in lineage2.2 and
not promoter variants as a class.

### Verification

Every point estimate in the table was reproduced by a second route that imports
none of this project's modules, takes variant status from CRyPTIC's own `IS_NULL`
and `IS_MINOR` flags, parses the concentration ladders from `PLATE_LAYOUT`, and
writes the censored likelihood out again. The two agree to three decimals on all
four variants and on the withheld one. The class-level shifts recomputed on the
atlas frame return the values in `outputs/mic_estimates.csv` before any variant is
removed. Two consecutive runs of the builder produce a byte-identical table.

## Minor alleles and null calls

The cohort definitions exclude an isolate from both the reference group and the
carrier groups when a target gene carries a null call, a het call or a minor
allele. That exclusion is right for attribution and it pools two different
things. A null call is a position that could not be read, which is missing data.
A het call or a minor indel is a variant detected in a minority of reads, which
is a measurement rather than the absence of one.

289 excluded isolates carry an MIC. 200 carry only a detected minor allele, 88
only null calls, and 1 both.

| Subset | BDQ resistant | CFZ resistant |
| --- | --- | --- |
| reference group | 0.56% | 3.71% |
| minor allele only, any of the three genes | 34 of 199, 17.1% | 56 of 200, 28.0% |
| null calls only, any of the three genes | 1 of 78, 1.3% | 4 of 79, 5.1% |
| minor allele at Rv0678 | 32 of 139, 23.0% | 51 of 140, 36.4% |
| null call at Rv0678 | 1 of 42, 2.4% | 2 of 43, 4.7% |
| minor allele at pepQ | 2 of 62, 3.2% | 5 of 62, 8.1% |
| null call at pepQ | 0 of 46, 0.0% | 2 of 46, 4.3% |

The excess resistance in the excluded group belongs to the detected minor alleles
at Rv0678. The null calls sit beside the reference group, which is what missing
data should do.

### Resolving a minor allele to a variant

`MINOR_MUTATION` carries the resolved form of a detected minor allele, so
`141_minorindel` reads as `141_ins_c` and `C46Z` as `C46G`. It is populated for
958 of the 959 Rv0678 minor-allele rows and 955 of those parse under the same
grammar parser the major alleles use.

Het calls resolve to 346 substitutions, 36 promoter variants and 35 stop codons.
Minor indels resolve to 524 frameshifts, 14 in-frame indels and 3 substitutions.
So 559 of 958 resolve to a loss-of-function class. The commonest resolved forms
are the same variants that dominate the major-allele layer: `141_ins_c` in 139
rows, `192_ins_g` in 66, `138_ins_g` in 58.

Eligibility for an estimate is strict, because the point is attribution: no real
major-allele variant in the three genes, no null call in any of them, no minor
allele outside Rv0678, and exactly one minor allele at Rv0678. That leaves 240
isolates, 68 of them carrying an MIC. 289 isolates are excluded for carrying more
than one minor allele at Rv0678, and one resolved form does not parse.

### What a minor allele is worth, with site held constant

Shifts in doublings against the reference group, to one decimal. The major rows
use the same cohort, the same filters, the same estimator and the same resampling
as the minor rows, so the four are comparable. `outputs/heteroresistance_estimates.csv`
carries three decimals and the plain fits beside these.

| Group | Isolates | Clusters | BDQ shift | 95% | CFZ shift | 95% |
| --- | --- | --- | --- | --- | --- | --- |
| minor loss of function | 40 | 36 | 1.9 | 1.3 to 2.5 | 1.6 | 1.0 to 2.1 |
| minor substitution | 23 | 22 | 1.3 | 0.7 to 1.9 | 1.0 | 0.4 to 1.5 |
| major loss of function | 159 | 67 | 1.0 | -0.2 to 2.2 | 1.0 | 0.3 to 1.7 |
| major substitution | 195 | 111 | 1.3 | 0.9 to 1.7 | 1.1 | 0.8 to 1.4 |

The isolate and cluster counts are the bedaquiline ones. One substitution carrier
has no clofazimine MIC, so the clofazimine figures in that row rest on 194
isolates in 110 clusters.

Site is held constant because 84 of the 139 Rv0678 minor-allele isolates come
from the one site enriched for bedaquiline resistance by design, and the sites
differ in mean MIC. Every group is refitted jointly with the reference group in
one censored regression carrying an indicator per site, so the shifts above and
the plain fits in the released table differ by the site adjustment alone.

Every minor-allele interval excludes zero on both drugs. The minor and major
intervals overlap, so what the data supports is a shift of the same size for a
variant detected in a minority of reads as for the same class carried in every
read. The bedaquiline major loss-of-function interval crosses zero on 159
isolates because those isolates sit in 67 clusters, one of which holds 41 of them,
and resampling clusters gives that outbreak one vote.

CRyPTIC report minor alleles as Susceptible, so this is a quantity the catalogue
layer discards by construction.

### Samples carrying more than one minor allele

A sample with two detected minor alleles at Rv0678 cannot attribute its MIC to
either, which is why the estimates above exclude it. It is still a sample with no
wild-type assertion available at the gene and no major allele to explain it, so
the group is worth measuring on its own terms.

289 samples across all 54,057 genomes carry more than one minor allele at Rv0678.
257 of them carry nothing else in the three genes, and 49 of those carry a
placeable MIC. That last figure is the analysable group.

Those 49 samples fall into 49 clusters, one each. They carry no defining
major-allele mutation to share, so clonal expansion inflates nothing here, which
is not true of any other carrier group in this work. 38 carry two alleles, 9
carry three and 2 carry four. All 111 alleles sit at distinct positions, in every
one of the samples, so no two compete at a locus. 48 of the 49 carry at least one
loss-of-function allele and 19 carry nothing else.

| Drug | Subset | Isolates | Resistant | Reference |
| --- | --- | --- | --- | --- |
| BDQ | all | 49 | 9, 18.4% | 0.57% |
| BDQ | two alleles | 38 | 7, 18.4% | 0.57% |
| BDQ | three or more | 11 | 2, 18.2% | 0.57% |
| CFZ | all | 49 | 19, 38.8% | 3.75% |
| CFZ | two alleles | 38 | 15, 39.5% | 3.75% |
| CFZ | three or more | 11 | 4, 36.4% | 3.75% |

The reference rates are over the reference samples carrying a placeable MIC for
that drug, 14,013 and 14,038, rather than all 14,187. Carrying three or more
alleles is no worse than carrying two on either drug.

With site held constant, in the same joint regression as the groups above, the
shift is 1.8 doublings for bedaquiline and 1.3 for clofazimine, against 1.9 and
1.6 for a single minor loss-of-function allele and 1.0 on both drugs for a
major-allele loss of function. Adding this group to that regression moves every other group's
estimate by at most 0.002 doublings.

### A reading the read fractions do not support

The read fractions of a sample's alleles sum to a median of 0.923, which invites
reading the remainder as a wild-type share: several sub-populations, each
carrying a different variant, with little wild type left. The fractions do not
support it.

12 of the 49 samples sum above 1.0, to a maximum of 1.077, which no partition can
do. Among the 38 two-allele samples the two fractions differ by a median of 0.309
and by more than 0.3 in 20 of them, so the alleles mostly do not sit at
comparable shares. FRS is a ratio of supporting reads to coverage at one
position, and coverage differs between positions, so a sum across positions has
no common denominator and is bounded by nothing. The summed fraction is reported
in `outputs/multi_allele_counts.csv` as a measurement and carries no inference
about population structure.

### A prediction that failed

If a minor allele acts through the sub-population carrying it, the MIC should
rise with the fraction of reads supporting the allele. Fitting the mean log2 MIC
as a line in that fraction, over the same censoring intervals with a common
standard deviation, gives a slope indistinguishable from zero everywhere.

| Drug | Subset | Isolates | Clusters | Slope, doublings per unit | 95% |
| --- | --- | --- | --- | --- | --- |
| BDQ | all | 68 | 63 | 0.0 | -1.7 to 1.8 |
| BDQ | loss of function | 40 | 36 | 0.2 | -2.3 to 2.7 |
| BDQ | substitution | 23 | 22 | -1.2 | -4.1 to 2.0 |
| CFZ | all | 68 | 63 | -0.4 | -2.3 to 1.3 |
| CFZ | loss of function | 40 | 36 | -0.1 | -2.2 to 2.2 |
| CFZ | substitution | 23 | 22 | -0.7 | -2.6 to 2.6 |

The read fractions run from 0.077 to 0.890 with a median of 0.460, so their range
is not the limitation. The sample is. At 40 isolates the interval admits anything
between about -2.3 and +2.7 doublings per unit, so a moderate dose-response is
not excluded; it is not found. The estimator is not the reason: on 3,000 to 4,000
simulated isolates censored onto the real bedaquiline ladder it recovers a planted
slope of 3.000 and a planted slope of zero.

Everything in this section is exploratory. The group was chosen for analysis after
the excluded isolates were seen to carry more resistance than the reference group,
and no comparison here was pre-specified.

## Resistance the genotype does not explain

An isolate is genotypically unexplained when it is resistant and carries no real
major-allele variant in Rv0678, pepQ or atpE. `code/unexplained.py` computes the
counts from mutation strings alone, so they sit beside the rest of the released
data under the same licence.

Four tiers partition the resistant isolates. **Attributable** carries at least one
variant in a class whose shift was estimated on this cohort with an interval above
zero: Rv0678 loss of function, Rv0678 substitution, pepQ. **Carrier, no
demonstrated effect** carries variants only in classes with no estimated upward
shift here: Rv0678 promoter variants, Rv0678 in-frame indels, atpE.
**Unexplained** carries no real major-allele variant in the three genes.
**Indeterminate** has one of the three genes carrying a null call, a het call or a
minor allele, so the isolate can be asserted neither to carry a variant nor to
lack one.

The boundary between the first two tiers rests on effect estimates computed on
this same cohort, so that split is descriptive rather than independent evidence.
The unexplained tier is structural and uses no estimate. mmpL5 enters no tier.

| Drug | Resistant | Attributable | Carrier, no demonstrated effect | Unexplained | Indeterminate |
| --- | --- | --- | --- | --- | --- |
| BDQ | 171 | 54 | 2 | 80 | 35 |
| CFZ | 683 | 94 | 2 | 527 | 60 |

The unexplained fraction, with Wilson intervals:

| Drug | Phenotype quality | Unexplained | Resistant | % | 95% interval |
| --- | --- | --- | --- | --- | --- |
| BDQ | all | 80 | 171 | 46.8 | 39.5 to 54.3 |
| BDQ | HIGH | 54 | 141 | 38.3 | 30.7 to 46.5 |
| CFZ | all | 527 | 683 | 77.2 | 73.9 to 80.2 |
| CFZ | HIGH | 327 | 426 | 76.8 | 72.5 to 80.5 |

The fraction persists at high phenotype quality, so it is not measurement error
alone. The intervals describe sampling error inside this collection and are not
population estimates.

### The figure is about a gene set, not a genome

| Gene set | BDQ unexplained | % | CFZ unexplained | % |
| --- | --- | --- | --- | --- |
| Rv0678 alone | 87 of 171 | 50.9 | 550 of 683 | 80.5 |
| Rv0678, pepQ, atpE | 80 of 171 | 46.8 | 527 of 683 | 77.2 |
| those three and mmpL5 | 0 of 171 | 0.0 | 6 of 683 | 0.9 |

Admitting mmpL5 removes the quantity, because almost every isolate carries a
variant in it. The indeterminate count rises with the gene set as well, 33 to 35
to 41 for bedaquiline and 53 to 60 to 84 for clofazimine, since each gene brings
its own unreadable positions. Any unexplained fraction therefore states something
about a named gene set and nothing about the genome.

Unexplained isolates carry no defining mutation and so cannot be clustered on one.
Grouped by site and sublineage, the 80 bedaquiline-resistant unexplained isolates
span 44 combinations with 8 in the largest, and the 527 clofazimine-resistant ones
span 122 combinations with 62 in the largest. That is a coarse floor on the number
of independent events behind the fraction.

## Genotype rules as tests for resistance

A shift says how far a variant moves the MIC. It does not say how much of a
collection's resistance a rule would catch, or what it would wrongly call
resistant. `outputs/prediction_metrics.csv` reports sensitivity, specificity and
the predictive values for four rules, with intervals from resampling clusters.

A null call is not a detected variant, so an isolate whose gene could not be read
is called not-resistant by every rule. That is how a catalogue applied to a real
sequence behaves, and counting it otherwise would flatter the rule.

| Drug | Rule | Sensitivity | Specificity | PPV |
| --- | --- | --- | --- | --- |
| BDQ | major variant in the three genes | 34.5 (27.0 to 43.4) | 95.7 (94.1 to 96.8) | 8.4 (5.8 to 11.8) |
| BDQ | any Rv0678 major variant | 31.6 (24.2 to 39.8) | 97.2 (96.0 to 98.2) | 11.4 (8.5 to 17.0) |
| BDQ | Rv0678 loss of function | 16.4 (9.6 to 24.7) | 99.1 (98.4 to 99.5) | 16.6 (8.9 to 29.1) |
| BDQ | major variant or detected minor allele | 52.6 (44.6 to 61.4) | 94.6 (93.4 to 95.6) | 10.2 (7.9 to 13.0) |
| CFZ | major variant in the three genes | 15.2 (11.0 to 19.7) | 95.8 (94.6 to 97.0) | 14.8 (10.9 to 18.9) |
| CFZ | any Rv0678 major variant | 12.7 (8.9 to 17.3) | 97.3 (96.2 to 98.2) | 18.4 (12.5 to 25.8) |
| CFZ | Rv0678 loss of function | 6.2 (3.5 to 9.6) | 99.1 (98.5 to 99.5) | 24.9 (13.5 to 41.7) |
| CFZ | major variant or detected minor allele | 22.2 (18.1 to 27.0) | 94.9 (93.5 to 96.0) | 17.2 (13.6 to 20.7) |

Denominators are 14,972 bedaquiline isolates in 14,573 clusters with 171
resistant, and 14,997 clofazimine isolates in 14,598 clusters with 683 resistant:
the cohort less the isolates carrying no MIC.

Admitting detected minor alleles raises bedaquiline sensitivity from 34.5% to
52.6% at a cost of 1.1 points of specificity, 643 false positives becoming 794.
For clofazimine it raises sensitivity from 15.2% to 22.2% at a similar cost. That
is the practical form of the minor-allele shifts above.

The predictive values describe this collection, whose resistance is concentrated at
one site by design, and transfer to no other. Sensitivity and specificity are less
exposed to that, and still describe this collection's particular mix of variants.

The intervals are wide for a quantity computed on fourteen thousand isolates
because the clusters are uneven. The 643 bedaquiline false positives sit in 268
clusters with 60 in the largest, so a Wilson interval on isolates gives
specificity 95.3 to 96.0 where resampling clusters gives 94.3 to 96.8, three and a
half times wider. The cluster interval is the one reported.

The rule counts and the tier counts above differ by a known amount. The
major-variant rule calls 59 of the 171 bedaquiline-resistant isolates where the
tiers count 56 carriers, and the difference is exactly the 3 resistant isolates
that carry both a major-allele variant and an uncertain call: the tiers place them
in the indeterminate tier because uncertainty takes precedence there, while the
rule calls them. For clofazimine it is 104 called against 96 tiered, with 8 in
both states. Every cell of every confusion matrix was reproduced by a route that
reads the source tables directly and uses none of this project's cohort code.

### The same rules at other cut-offs

The ECOFF is one concentration among those the plates tested, and a rule that
looks weak against it can look different against another.
`outputs/prediction_thresholds.csv` carries all four rules at every usable
cut-off: seven for bedaquiline, 0.015 to 1 mg/L, and six for clofazimine, 0.06 to
2 mg/L.

A cut-off is usable only where every isolate's position relative to it is
determinate. Below the highest of the two designs' lowest rungs, a left-censored
reading on the design with the lower floor lies on neither side of it; above the
lowest of the two designs' highest rungs, a right-censored reading lies on
neither side either. A concentration that is a rung on both designs is inside
both bounds, so the shared rungs are exactly the usable cut-offs, and no isolate
is dropped at any of them.

Position is read from the censoring interval rather than from the reported
number: an interval lies above the cut-off when its lower bound reaches it, and
at or below when its upper bound does not exceed it. At the ECOFF that rule and
CRyPTIC's own resistant or susceptible call agree on all 14,972 bedaquiline and
14,997 clofazimine isolates, with no disagreement, which is an independent check
on the censoring placement.

The broadest rule, a major variant or a detected minor allele, across the
bedaquiline range:

| Cut-off, mg/L | Above it | Sensitivity | Specificity | PPV |
| --- | --- | --- | --- | --- |
| 0.015 | 11446 | 6.0 (5.0 to 7.2) | 94.5 (90.9 to 97.3) | 77.9 (67.9 to 88.0) |
| 0.03 | 6529 | 8.5 (7.0 to 10.1) | 96.1 (94.5 to 97.4) | 62.7 (53.7 to 71.9) |
| 0.06 | 1614 | 21.9 (18.8 to 24.7) | 96.0 (94.8 to 97.0) | 40.1 (34.4 to 45.9) |
| 0.12 | 508 | 45.7 (39.4 to 51.1) | 95.5 (94.0 to 96.6) | 26.2 (21.2 to 33.2) |
| 0.25 | 171 | 52.6 (44.6 to 61.4) | 94.6 (93.4 to 95.6) | 10.2 (7.9 to 13.0) |
| 0.5 | 79 | 27.9 (18.0 to 39.8) | 94.2 (92.4 to 95.2) | 2.5 (1.6 to 3.7) |
| 1 | 33 | 15.2 (3.7 to 28.2) | 94.1 (92.7 to 95.2) | 0.6 (0.1 to 1.0) |

And across the clofazimine range:

| Cut-off, mg/L | Above it | Sensitivity | Specificity | PPV |
| --- | --- | --- | --- | --- |
| 0.06 | 5599 | 9.1 (7.6 to 11.1) | 96.0 (94.6 to 97.1) | 58.0 (51.3 to 64.7) |
| 0.12 | 2186 | 13.8 (11.7 to 16.8) | 95.5 (94.2 to 96.4) | 34.2 (29.1 to 39.2) |
| 0.25 | 683 | 22.2 (18.1 to 27.0) | 94.9 (93.5 to 96.0) | 17.2 (13.6 to 20.7) |
| 0.5 | 238 | 31.5 (24.5 to 38.3) | 94.5 (93.3 to 95.7) | 8.5 (6.3 to 11.2) |
| 1 | 81 | 35.8 (22.5 to 48.2) | 94.3 (93.0 to 95.4) | 3.3 (1.9 to 4.8) |
| 2 | 30 | 26.7 (10.7 to 45.0) | 94.2 (92.9 to 95.3) | 0.9 (0.3 to 1.8) |

The ECOFF sits at 0.25 mg/L in both tables. Two things are visible there that the
single-cut-off view hides.

Sensitivity on bedaquiline peaks at the ECOFF, at 52.6%, and falls away on both
sides: 6.0% at 0.015 mg/L, where 11,446 of the 14,972 isolates sit above the
cut-off and the rules call a fixed few hundred, and 15.2% at 1 mg/L, where only
33 isolates remain above it and most are ones no rule calls. Specificity moves
between 94.1% and 96.1% across the whole range, so the isolates these rules call
barely track the cut-off.

Clofazimine peaks elsewhere. Sensitivity climbs from 9.1% at 0.06 mg/L to 35.8%
at 1 mg/L, four doublings above the ECOFF, and falls to 26.7% only at the top of
the range where 30 isolates remain. The genotype accounts for the highest
clofazimine MICs better than for those just over the breakpoint, where the ECOFF
figure of 22.2% is measured. That asymmetry between the two drugs is what a
single cut-off cannot express.

## The pre-registered test

Everything above is an estimate with an interval. One comparison was specified in
advance, committed, and then tested against a held-out half of the cohort.

### The split

`WGS_SAMPLES` carries a dataset label. Within the analysis cohort: CRyPTIC-v1.0
holds 12,387 isolates, CRyPTIC-v2.0 holds 2,761, CRyPTIC-v3.0 holds 10.

CRyPTIC designate v3.0 as a validation set, and it cannot test an MIC-level claim:
only 10 of its 9,090 isolates carry a UKMYC MIC and none carries an Rv0678
variant, because those isolates were sequenced without the broth microdilution
plates. The split used was v1.0 against v2.0. v1.0 is the frozen pre-2020
collection; v2.0 is everything added after that freeze, including roughly 1,100
isolates from one site enriched for bedaquiline resistance. The halves are
separated in time and partly in geography, which is closer to external validation
than to cross-validation.

Reference resistance is comparable across the halves, 0.60% against 0.43% for
bedaquiline and 3.98% against 2.51% for clofazimine. The split was fixed on those
grounds before the power analysis showed how few comparisons it could test, and
was not revised afterwards.

### Prior exposure, disclosed

This is a replication in a later subset with prior exposure, not a clean hold-out.
The exploratory work that preceded the split used all 15,158 isolates, including
every v2.0 isolate: the class effects, the stratification by site and lineage, the
clustering correction and the censored estimation. The comparison that was
registered was chosen knowing the direction and the rough size of the effect in
the whole cohort. What the test establishes is that the magnitude holds in the
later half under criteria fixed beforehand, not that the comparison was chosen
blind.

### What a prediction interval has to be

A prediction is not that the held-out estimate lands inside the discovery
confidence interval. The held-out estimate carries its own sampling error and the
held-out groups are much smaller, so that criterion would miss far more often than
5% of the time even when the effect is exactly as estimated.

Each interval here is a prediction interval for the held-out estimate itself,
built by simulating datasets of the held-out size, with the held-out cluster
structure, censored onto the same plate ladders, and re-estimated. The true effect
in each simulated dataset is drawn from the discovery estimate's bootstrap
distribution rather than fixed at the point estimate, so the interval carries both
sources of error. It is taken as the mean plus and minus 1.96 standard deviations
rather than as empirical percentiles, which are unstable and biased inward at a
few hundred draws.

Over full discovery-then-validate cycles with a known planted truth, the interval
built this way covers a newly simulated estimate in 98% of cycles against a
nominal 95%, slightly conservative. An interval fixed at the discovery point
estimate covers in 84%, and the discovery confidence interval covers in 74%.

### The registered predictions and the result

`docs/PRE_REGISTRATION.md` was committed before the test ran, and `validate.py`
refuses to run if that file is absent, uncommitted or modified.

Rv0678 loss of function against the reference group, bedaquiline, on 23 held-out
isolates in 19 clusters:

| Prediction | Quantity | Predicted | Observed | Verdict |
| --- | --- | --- | --- | --- |
| P1a | Odds ratio | 17.50 to 278.00 | 81.99 | Supported |
| P1b | MIC shift, doublings | 1.60 to 3.22 | 2.55 | Supported |

P1a: 23 held-out isolates in 19 clusters with 7 resistant, against 1,884 reference
isolates with 10 resistant. Discovery estimate 73.0, smallest detectable odds
ratio 27.4, observed 82.0 with its own interval 21.2 to 296.5.

P1b: reference mean -5.17, group mean -2.61, observed shift 2.55 doublings.
Discovery estimate 2.43, smallest detectable shift 1.00 doublings, and the
observed shift carries its own interval of 1.93 to 3.11.

The predicted ranges, the discovery estimates and the smallest detectable effects
are the registered ones, read out of the committed pre-registration. That document
records 22 held-out isolates in 18 clusters, which is what the classification in
use when it was committed gave; the corrected classification puts 23 in 19 in the
same held-out half, which is where the observed figures above come from. Both
observed values fall inside the registered ranges either way.

Two of two supported. Both concern the same comparison on the same 23 isolates, so
this is one replication rather than two independent ones, and the held-out sample
is small enough that its own intervals are wide.

## Limitations

1. Resistance in this collection is concentrated at one site by design, so no
   figure here is a prevalence estimate for any population, and the predictive
   values in particular transfer to no other collection.
2. The validation carries prior exposure, disclosed above and in the
   pre-registration itself.
3. One comparison was testable. Everything else is an estimate with an interval
   rather than a tested hypothesis, and many comparisons were made during the
   exploratory phase without correction for multiple testing.
4. The effective sample sizes after clustering are 63 clusters for loss of
   function and 111 for substitution. Every interval rests on those, not on the
   isolate counts.
5. Where censoring is heavy the estimate leans on the assumption that log2 MIC is
   normally distributed within a group. 41.2% of the clofazimine reference group's
   MICs are left-censored, and 40.4% across the whole cohort.
6. Clustering is approximated by site, sublineage and defining mutation. A
   transmission analysis would cluster on genomic distance, which needs the full
   variant table and more compute than this project has. The approximation merges
   unrelated isolates and so understates the evidence rather than inflating it.
7. The per-variant layer reaches 4 variants of 217. The rest carry counts and no
   shift, 182 of them resting on a single cluster. A quantitative resource built
   on this dataset is a resource about a handful of variants and a large amount of
   carefully bounded ignorance.
8. Variant naming here is not interoperable with the WHO catalogue or
   TB-Profiler, because HGVS normalisation is not implemented, so the per-variant
   table cannot yet be joined to either by variant name.
9. Delamanid and linezolid are out of scope. In this cohort delamanid is 61.8%
   left-censored with 43.7% of its resistant isolates at the plate ceiling, and
   linezolid's plates carry one tested concentration above the breakpoint on
   UKMYC5 and two on UKMYC6, with 32.7% of its resistant isolates off the top of
   the plate. Section 6 of `outputs/audit_report.txt` carries both profiles.
