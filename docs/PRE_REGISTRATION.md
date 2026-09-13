# Pre-registration: predictions to be tested on the held-out half

Researcher: Amina Baktiyarova, Independent Researcher,
ORCID 0009-0007-6265-6493

Every estimate below comes from the CRyPTIC-v1.0 half of the CRyPTIC
Consortium Dataset v3.4.0, Zenodo version DOI 10.5281/zenodo.15680920.

## Prior exposure, disclosed

This is a replication in a later subset with prior exposure, not a
clean hold-out, and the difference is stated here rather than left for
a reader to discover.

Exploratory analyses run before this split was defined used the full
cohort of 15,158 samples, which includes all 2,761 CRyPTIC-v2.0
samples. Those analyses covered variant class effects, stratification
by site and lineage, correction for clonal clustering, and interval-
censored MIC estimation. A cohort audit additionally reported
CRyPTIC-v2.0 reference resistance rates directly.

What follows from that. The choice of which comparisons to carry
forward was informed by results computed on combined data, so the
selection is not independent of the held-out half. The quantitative
predictions below were computed on CRyPTIC-v1.0 alone, which the code
in code/discovery.py makes checkable, and the decision criteria and
prediction intervals are fixed before code/validate.py is run.

So this test can establish that an effect estimated in the earlier
collection holds in a later and partly different one, under criteria
fixed in advance. It cannot establish that the comparisons were chosen
without any knowledge of the held-out samples. Any write-up must say
so.

## The split

CRyPTIC-v3.0, which CRyPTIC designate as a validation set, holds 9,090
sequenced samples but only 10 with a UKMYC MIC and none carrying an
Rv0678 variant, so it cannot test a claim about MICs. The split used is
CRyPTIC-v1.0 against CRyPTIC-v2.0. v1.0 is the frozen pre-2020
collection handed to FIND and Seq&Treat to build the first WHO
catalogue. v2.0 is everything added after that freeze, including
approximately 1,100 samples from NICD enriched for bedaquiline
resistance.

The halves are separated in time and partly in geography, so this is
closer to external validation than to cross-validation. That makes it a
stronger test than a random split, and it also means a failure to
replicate may reflect a difference between the populations rather than
a false finding in the first. Any write-up must state both.

This split was chosen on scientific grounds, before the power analysis
showed how few comparisons it could test. It is not revised now that
the answer is known, because choosing a split after seeing which one
yields more testable predictions is the specific practice that
pre-specification exists to prevent.

## Cohort definition

Applied identically to both halves:

- Samples with both a UKMYC5 or UKMYC6 MIC and a genome.
- Phenotype quality HIGH only, meaning at least two independent reading
  methods agreed.
- One isolate per patient, the first by identifier.
- Samples with a disrupted mmpL5 excluded from the variant groups.
- Variant groups defined from mutation strings alone. No WHO catalogue
  content, and neither the EFFECTS nor the PREDICTIONS table, is used.

## What counts as a successful prediction

A held-out estimate carries its own sampling error, and the held-out
groups are much smaller than the discovery groups. So a prediction is
not that the held-out estimate lands inside the discovery confidence
interval: it would miss far more often than 5% of the time even if the
effect were exactly as estimated, and counting that as a failure to
replicate would be wrong.

Each prediction interval below is instead a prediction interval for the
held-out estimate itself. It was built by simulating datasets of the
held-out size, with the held-out cluster structure, censored onto the
same plate dilution ladders, and re-estimated. The interval is the
central 95% of those simulated estimates, taken as the mean plus and
minus 1.96 standard deviations rather than as empirical percentiles,
which are unstable and biased inward at a few hundred draws.

The true effect used in each simulated dataset is itself drawn from the
discovery estimate's bootstrap distribution rather than fixed at the
discovery point estimate, so the interval carries both the uncertainty
in the discovery estimate and the sampling error of the held-out half.
Fixing it at the point estimate would make the interval too narrow and
would turn ordinary sampling variation into recorded failures.

If the effect is real and the same size in both halves, the held-out
estimate should land inside the interval about 95% of the time.

A prediction is supported if the held-out estimate falls inside the
stated prediction interval and is in the stated direction. Failing
either way is recorded as a failure.

Comparisons that cannot be tested are reported as untested estimates
with intervals, not dropped and not described as validated. The
discovery report in outputs/discovery_report.txt lists every comparison
and states which of them the held-out half is powered to test.

## What is registered, and what is not

Only comparisons the held-out half is powered to test are registered. A
comparison qualifies when the lower bound of the discovery interval,
not the point estimate, exceeds the smallest effect detectable at 80%
power. Using the point estimate would register comparisons that could
fail for lack of samples rather than lack of effect, leaving the result
uninterpretable.

Binary power was found by simulation. MIC power comes from the
estimator's standard error under the null, measured by simulating
datasets of the held-out size, with a shift over standard error of 2.80
required for 80% power at a two-sided 5% level.

The two tests are registered separately, because the MIC test uses the
measurement itself and the binary test uses it only through a
threshold, so one can be powered where the other is not.

## Predictions

### P1. Rv0678 loss of function, BDQ

Held-out samples: 22, in 18 clusters.

**P1a, resistance.** The held-out odds ratio against the reference group will fall between 17.5 and 278.0. Discovery estimate 73.0; smallest odds ratio the held-out half can detect, 27.4.

**P1b, MIC.** The held-out fitted mean log2 MIC will exceed the reference group's, by between 1.60 and 3.22 doublings. Discovery estimate 2.43; smallest shift the held-out half can detect, 1.00 doublings.

## Analysis to be run

code/validate.py, once, on CRyPTIC-v2.0. No estimate in this document
will be revised afterwards. If a prediction fails, that is reported as a
failure rather than explained away or replaced with a different
comparison.

## Data

The CRyPTIC Consortium Dataset, version v3.4.0, Zenodo version DOI
10.5281/zenodo.15680920, CC BY 4.0. See docs/PROVENANCE.md for the file
inventory and checksums.

