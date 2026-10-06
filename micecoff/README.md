# micecoff

Interval-censored estimation of minimum inhibitory concentration (MIC)
distributions from broth microdilution readings. For each group of readings,
micecoff reports the fitted MIC distribution, the shift of the group against a
reference group, and the epidemiological cut-off (ECOFF) at a chosen quantile of
a fitted wild-type distribution. Every estimate carries a 95% interval from
resampling.

## Readings as intervals

A broth microdilution plate tests a series of concentrations. A reported MIC of
0.25 mg/L on a plate that also tests 0.12 mg/L means growth was inhibited at
0.25 and not at 0.12, so the MIC lies in the interval (0.12, 0.25]. A reading of
`<=0.008` at the lowest tested concentration places the MIC anywhere at or below
0.008, and `>1` at the highest places it anywhere above 1. A median of the
reported numbers treats these intervals as points and is biased upward, toward
the concentration where growth stopped, and where many readings sit at the
plate floor the median reports the floor.

micecoff models log2 MIC within a group as normally distributed and estimates
the mean and standard deviation by maximum likelihood over the censoring
intervals. Each reading contributes the probability mass the fitted normal
places on its interval, computed in log space so that an interval far in a tail
keeps its precision.

## Installation

micecoff requires Python 3.10 or later, numpy 1.26 or later and scipy 1.11 or
later. It is not published to a package index. From a clone of the repository
at https://github.com/aminabaktiyarova/bdq-cfz-mic-atlas, the command below
installs it and puts the `micecoff` command on the path, and
`python -m micecoff` runs the same interface.

```
pip install .
```

## Input

The input is a CSV file with a header row and one reading per row, or `-` to
read standard input. Fields are separated by commas, the file is read as UTF-8,
and a byte order mark at its start is ignored. The MIC column holds each reading
as the plate reports it: a tested concentration such as `0.25`, `<=c` where c is
the lowest tested concentration, or `>c` or `>=c` where c is the highest. An
empty cell, `nan`, `NaN`, `NA`, `NaT`, `None` or `<NA>` is a missing reading,
which counts in `absent` and enters no fit. A group column and a cluster column,
when named, must carry a label on every row, and labels are compared exactly as
written, so `wild type` and ` wild type` are two groups. The excerpt below is
from an input file, `mics.csv`, with a group column and a cluster column.

```
MIC,group,cluster
<=0.008,wild type,site01
0.015,wild type,site02
0.03,wild type,site01
0.12,carrier,site03
>1,carrier,site04
```

`--series` lists the concentrations the plate tested, comma separated, in the
units the readings use. It includes every tested concentration, among them any
that no reading reports, because each well bounds the interval of the reading
above it. The series takes plain concentrations: at least two, none repeated,
none at or below zero, and none carrying an operator. A series that is not a
doubling series runs with a warning, and the interval below each reading then
spans the gap to the next tested concentration.

A reading that matches no tested concentration within 5% stops the run. So does
a censored reading whose bound is not the end of the series, such as `<=0.015`
against a series starting at 0.008, which came from a plate testing another
series; readings from plates testing different series go in separate runs, one
series each. `--allow-off-series` drops both kinds of reading, counts them per
group in `off_series` and reports the total on standard error. The run also
stops on a repeated column name, a row whose field count differs from the
header's, an empty file, a header with no rows, and an empty group or cluster
label.

## Commands

`fit` reports each group's fitted distribution, `shift` reports each group's
shift against the reference group, and `ecoff` reports the cut-off at the given
quantile of each group's fitted distribution. Each group is fitted on its own,
so in the `ecoff` output only a wild-type group's row is an ECOFF. The `0.99`
below shows the syntax and is not a recommendation; the section on the cut-off
explains why micecoff sets no default.

```
micecoff fit mics.csv --series 0.008,0.015,0.03,0.06,0.12,0.25,0.5,1 --group-column group --cluster-column cluster
micecoff shift mics.csv --series 0.008,0.015,0.03,0.06,0.12,0.25,0.5,1 --group-column group --cluster-column cluster --reference "wild type"
micecoff ecoff mics.csv --series 0.008,0.015,0.03,0.06,0.12,0.25,0.5,1 --group-column group --cluster-column cluster --coverage 0.99
```

| Option | Subcommands | Default | Meaning |
| --- | --- | --- | --- |
| `input` | all | required | CSV path, or `-` for standard input |
| `--series` | all | required | Tested concentrations, comma separated |
| `--mic-column` | all | `MIC` | Column holding the reading |
| `--group-column` | all, required by `shift` | none | Column holding the group label; without it, every row is one group named `all` |
| `--cluster-column` | all | none | Column holding a cluster label; intervals then resample clusters, and single rows otherwise |
| `--draws` | all | `400` | Resampled datasets per interval, at least 100 |
| `--seed` | all | `20260101` | Base seed for the resampling, a non-negative integer |
| `--min-clusters` | all | `5` | Resampling units a group needs before an estimate is reported, at least 2 |
| `--allow-off-series` | all | off | Drop and count readings that do not match the series |
| `--reference` | `shift` | required | Label of the reference group |
| `--coverage` | `ecoff` | required | Quantile of the fitted distribution the cut-off must reach, strictly between 0 and 1 |
| `--version` | none | none | Print the version and exit |

## Output

The result is a CSV on standard output with one row per group, in sorted order
of label; `shift` writes every group except the reference. Notes and warnings
go to standard error. The exit status is 0 when the run completes, including a
run in which estimates are withheld, 1 when the input stops the run, and 2 when
the arguments are invalid. Values on the log2 scale are written to four decimal
places and concentrations to four significant figures. A withheld estimate
leaves its columns empty and keeps its counts, with the reason in `withheld`.

### fit

| Column | Type | Definition | Empty when |
| --- | --- | --- | --- |
| `group` | text | The group label as written in the input, or `all` without a group column | never |
| `rows` | integer | Rows carrying this label | never |
| `absent` | integer | Rows whose reading is missing | never |
| `off_series` | integer | Rows whose reading does not match the series, dropped under `--allow-off-series` | never |
| `placed` | integer | Rows whose reading became a censoring interval, which is `rows` less `absent` and `off_series` | never |
| `left_censored` | integer | Placed rows whose interval is open below | never |
| `right_censored` | integer | Placed rows whose interval is open above | never |
| `resampling_unit` | text | `cluster` with a cluster column, `row` without | never |
| `clusters` | integer | Distinct resampling units among the placed rows | never |
| `mean_log2` | log2 concentration | Maximum likelihood mean of log2 MIC | withheld |
| `mean_low` | log2 concentration | 2.5th percentile of the fitted mean over the resampled datasets | withheld |
| `mean_high` | log2 concentration | 97.5th percentile of the fitted mean over the resampled datasets | withheld |
| `sd_log2` | doublings | Maximum likelihood standard deviation of log2 MIC | withheld |
| `sd_low` | doublings | 2.5th percentile of the fitted standard deviation over the resampled datasets | withheld |
| `sd_high` | doublings | 97.5th percentile of the fitted standard deviation over the resampled datasets | withheld |
| `geometric_mean` | concentration | 2 raised to `mean_log2`, the geometric mean of the fitted MIC distribution, in the units of the series | withheld |
| `geometric_mean_low` | concentration | 2 raised to `mean_low` | withheld |
| `geometric_mean_high` | concentration | 2 raised to `mean_high` | withheld |
| `draws` | integer | Resampled datasets whose fit converged | withheld |
| `withheld` | text | Why the estimate is withheld | the estimate is reported |

### shift

| Column | Type | Definition | Empty when |
| --- | --- | --- | --- |
| `group` | text | The group label as written in the input | never |
| `reference` | text | The reference group's label | never |
| `rows` | integer | Rows carrying the group's label | never |
| `absent` | integer | The group's rows whose reading is missing | never |
| `off_series` | integer | The group's rows whose reading does not match the series, dropped under `--allow-off-series` | never |
| `placed` | integer | The group's rows whose reading became a censoring interval | never |
| `reference_placed` | integer | The reference group's rows whose reading became a censoring interval | never |
| `resampling_unit` | text | `cluster` with a cluster column, `row` without | never |
| `clusters` | integer | Distinct resampling units among the group's placed rows | never |
| `reference_clusters` | integer | Distinct resampling units among the reference group's placed rows | never |
| `mean_log2` | log2 concentration | Maximum likelihood mean of the group's log2 MIC | withheld |
| `reference_mean_log2` | log2 concentration | Maximum likelihood mean of the reference group's log2 MIC | withheld |
| `shift` | doublings | `mean_log2` less `reference_mean_log2`; 1.0 means the group's fitted mean MIC is twice the reference's | withheld |
| `shift_low` | doublings | 2.5th percentile of the shift over the resampled datasets | withheld |
| `shift_high` | doublings | 97.5th percentile of the shift over the resampled datasets | withheld |
| `draws` | integer | Resampled datasets in which both fits converged | withheld |
| `withheld` | text | Why the estimate is withheld | the estimate is reported |

### ecoff

| Column | Type | Definition | Empty when |
| --- | --- | --- | --- |
| `group` | text | The group label as written in the input, or `all` without a group column | never |
| `rows` | integer | Rows carrying this label | never |
| `absent` | integer | Rows whose reading is missing | never |
| `off_series` | integer | Rows whose reading does not match the series, dropped under `--allow-off-series` | never |
| `placed` | integer | Rows whose reading became a censoring interval | never |
| `resampling_unit` | text | `cluster` with a cluster column, `row` without | never |
| `clusters` | integer | Distinct resampling units among the placed rows | never |
| `mean_log2` | log2 concentration | Maximum likelihood mean of log2 MIC | withheld |
| `sd_log2` | doublings | Maximum likelihood standard deviation of log2 MIC | withheld |
| `coverage` | proportion | The quantile given with `--coverage` | never |
| `quantile_log2` | log2 concentration | The quantile of the fitted normal at `coverage` | withheld |
| `quantile_low` | log2 concentration | 2.5th percentile of that quantile over the resampled datasets | withheld |
| `quantile_high` | log2 concentration | 97.5th percentile of that quantile over the resampled datasets | withheld |
| `ecoff` | concentration | The lowest tested concentration at or above `quantile_log2`, or `>c` when the quantile lies above the highest tested concentration c | withheld |
| `ecoff_low` | concentration | The same rounding applied to `quantile_low` | withheld |
| `ecoff_high` | concentration | The same rounding applied to `quantile_high` | withheld |
| `draws` | integer | Resampled datasets whose fit converged | withheld |
| `withheld` | text | Why the estimate is withheld | the estimate is reported |

## Intervals

Each interval runs from the 2.5th to the 97.5th percentile of the estimate over
`--draws` resampled datasets. A resampled dataset draws, with replacement, as
many units as the group holds: clusters when `--cluster-column` names one, and
single rows otherwise, which treats the rows as independent. Isolates from one
clonal outbreak, or repeated isolates from one patient, are one observation
between them, so data carrying either belong under a cluster column. For
`shift`, each resampled dataset resamples the group's units and, separately,
the reference group's, refits both and takes the difference of the fitted
means, so the interval carries the sampling error of both. The `fit` intervals
on the mean and on the standard deviation are read from the same resampled
datasets.

Each group's interval draws from its own generator, seeded with the base seed
followed by the code points of a name: `fit|<group>`, `shift|<group>|<reference>`
or `ecoff|<group>`, with the labels as written in the input. A group's interval
therefore does not depend on which other groups the file holds, and numpy's
`default_rng([seed] + [ord(c) for c in name])` reproduces it on its own.

## Withheld estimates

An estimate is withheld for one of the reasons below, checked in this order.
The reason is written in `withheld` as shown, with each `<n>` replaced by a
count.

- `no reading was placed`
- `the readings occupy too few intervals of the series (<n> of the 3 a finite maximum needs)`,
  because with readings in one or two distinct intervals the likelihood can
  rise without limit toward a zero standard deviation or an unbounded mean.
- `the fit did not converge`
- `a parameter ended on the optimizer's bound`, where the optimizer searches
  means from -25 to 25 and standard deviations from 0.05 to 20 on the log2
  scale, and a fit ending on one of these limits found no interior maximum.
- `the fitted mean lies outside the tested series by more than 1 on the log2 scale`,
  which is an extrapolation past every measurement.
- `too few resampling units (<n> of the <n> an interval needs)`, the second
  count being `--min-clusters`.
- `fewer than a quarter of the resampling draws converged`
- `reference: ` followed by one of the reasons above, in `shift`, when the
  reference group's own estimate is withheld, which withholds every row.

## The cut-off

`ecoff` takes the quantile of the fitted normal at `--coverage` and rounds it up
to the lowest tested concentration at or above it, since a plate expresses a
cut-off only at a concentration it tests. A quantile above the highest tested
concentration c is written `>c`, and the plate then cannot express the cut-off.
`ecoff_low` and `ecoff_high` apply the same rounding to the ends of the interval
on the quantile, and a pair of different values means the cut-off could move by
at least one dilution.

An ECOFF is the upper end of the wild-type MIC distribution, identified with
statistical methods (Kahlmeter and Turnidge, 2023). The quantile that marks it
is a convention of the procedure a laboratory follows, which for EUCAST is SOP
10.2. micecoff takes the quantile as a required argument and sets no default.

`ecoff` fits one normal distribution to every placed reading in a group, so the
group has to hold the wild-type population alone, defined for example by
genotype. Non-wild-type isolates in a group widen the fit and raise the
cut-off, and micecoff does not separate a mixed distribution into its parts.

## Assumptions and limitations

- log2 MIC is taken as normal within a group. Where most readings are censored,
  the estimate leans on that assumption beyond the measured range.
- Percentile intervals from resampling run narrower than their nominal coverage
  in small groups.
- Rows are taken as independent unless a cluster column groups them.
- One run takes one tested series.

## License and source

MIT. Source: https://github.com/aminabaktiyarova/bdq-cfz-mic-atlas. Author:
Amina Baktiyarova, ORCID 0009-0007-6265-6493.

## References

EUCAST. SOP 10.2, MIC distributions and the setting of epidemiological cut-off
(ECOFF), version 10.2 (2021), re-affirmed 12 May 2025.
https://www.eucast.org/eucastsops

Kahlmeter G, Turnidge J. Clinical Microbiology Reviews 2023, 36(4):e00100-22.
https://doi.org/10.1128/cmr.00100-22
