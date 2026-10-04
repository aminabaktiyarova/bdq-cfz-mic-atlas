# Trajectory B: data provenance record

Researcher: Amina Baktiyarova, Independent Researcher, ORCID 0009-0007-6265-6493

This file records exactly which data version, licence and documentation this project rests
on, what the schema says, and what the input tables contain. Every figure published from
this project traces back to what is written here. What was computed from these inputs is
recorded separately, in RESULTS.md. Facts below were read directly off the Zenodo record, DATA_SCHEMA.pdf and
RELEASE_NOTES.md. Items marked PLACEHOLDER are not yet verified and must not be filled in
with an assumption.

## Dataset

Title: The CRyPTIC Consortium Dataset
Version: v3.4.0
Published: 21 May 2025
Confirmed latest version at time of recording: yes, no newer-version banner on the record
Version DOI (cite this in every output): 10.5281/zenodo.15680920
Collection DOI (all versions): 10.5281/zenodo.15679730
Record URL: https://zenodo.org/records/15680920
Licence: Creative Commons Attribution 4.0 International (CC BY 4.0)

Creators:
- The CRyPTIC Consortium (data collector)
- Fowler, Philip (data curator), ORCID 0000-0003-0912-4483, University of Oxford

Funding acknowledged on the record:
- Wellcome Trust, The CRyPTIC Consortium, 200205/Z/15/Z
- Bill and Melinda Gates Foundation, The CRyPTIC Consortium, OPP1133541

Schema document: "CRyPTIC Release 3.4.0 Schema", PWF, 21 May 2025, DATA_SCHEMA.pdf

## Contents as stated on the record

- 53,897 samples have both WGS and pDST data
- An additional 11,945 samples have pDST data only

The v3.4.0 record does not state how many samples have MICs measured on UKMYC5 or UKMYC6
plates. The v2.1.2 record stated 21,570 for that subset, but v2.1.2 was processed
differently and its figure does not carry over.
Computed from UKMYC_PHENOTYPES: 21,685 samples across 288,904 rows. Per-drug counts are in
"UKMYC_PHENOTYPES: actual shape" below.

## File inventory (19 files, 2.6 GB total)

| File | Size |
| --- | --- |
| BASHTHEBUG.parquet | 3.7 MB |
| BASHTHEBUG_CLASSIFICATIONS.parquet | 162.1 MB |
| COUNTRIES_LOOKUP.csv.gz | 4.6 kB |
| DATA_SCHEMA.pdf | 102.7 kB |
| DRUG_CODES.csv.gz | 385 Bytes |
| DST_MEASUREMENTS.parquet | 2.7 MB |
| DST_SAMPLES.parquet | 820.5 kB |
| EFFECTS.parquet | 5.9 MB |
| GENOMES.parquet | 2.6 MB |
| MUTATIONS.parquet | 1.1 GB |
| PLATE_LAYOUT.parquet | 5.9 kB |
| PREDICTIONS.parquet | 2.4 MB |
| RELEASE_NOTES.md | 18.8 kB |
| SITES.csv.gz | 1.3 kB |
| UKMYC_GROWTH.parquet | 15.1 MB |
| UKMYC_PHENOTYPES.parquet | 1.6 MB |
| UKMYC_PLATES.parquet | 2.4 MB |
| VARIANTS.parquet | 1.3 GB |
| WGS_SAMPLES.parquet | 9.4 MB |

MD5 checksums for every file are published on the Zenodo record. The copies this analysis
reads were verified against them: all nine match, including MUTATIONS.parquet at
1,062,006,199 bytes. No file differs from its published checksum at full size.

| File | Bytes | MD5 |
| --- | --- | --- |
| DATA_SCHEMA.pdf | 102,717 | `8b11bfbb9255da6dfc1b8b7aede767d0` |
| DRUG_CODES.csv.gz | 385 | `923d3a193df21698bd6a00f857ab337e` |
| GENOMES.parquet | 2,604,147 | `ebd82e85f71e36de5da10e776b6afe4e` |
| MUTATIONS.parquet | 1,062,006,199 | `d5feeeae14304006ba67aaaef84cff03` |
| PLATE_LAYOUT.parquet | 5,916 | `cb403c7517ec847467b7980cbc3e5389` |
| RELEASE_NOTES.md | 18,798 | `dfac1d2007ad30bb749f7bb3bcb9645b` |
| SITES.csv.gz | 1,325 | `c24c882c8988b5af9940232ada27fb60` |
| UKMYC_PHENOTYPES.parquet | 1,559,813 | `020b6c0af6c05e19610a59f5ef97b832` |
| WGS_SAMPLES.parquet | 9,410,557 | `ea798f4cfc28525cf394ff9196c93021` |

Four further files downloaded as 92-byte stubs rather than data: DST_MEASUREMENTS.parquet,
DST_SAMPLES.parquet, UKMYC_GROWTH.parquet and UKMYC_PLATES.parquet. All four carry the same
checksum as one another, which is the signature of a server response rather than of
truncated data, and none is read by any module here. They are recorded as absent rather
than treated as present.

Naming mismatch to watch: the schema calls the drug lookup table DRUG_CODE, the file on
Zenodo is DRUG_CODES.csv.gz.

## Differences from v2.1.2 that affect this project

1. Plate images in v3.4.0 were read by TMAS, a convolutional neural network, replacing
   AMyGDA. The release notes state this raised the proportion of MICs held at high
   confidence (at least two of the three measurement methods agreeing) from 79.2% to 88.7%.
   Neither figure reproduces from the table shipped with this version; see the section on
   the high-confidence proportion below.
2. All CRyPTIC 96-well plate data was redownloaded from clires2.org, picking up 425
   additional samples from National University of Singapore. Their FASTQ files are not yet
   processed, so these phenotypes have no matching genetics.
3. Two BashTheBug tables were added, holding 4.75 million volunteer classifications.
4. Genetics were processed through a Mycobacterial pipeline on the EIT Pathogena cloud
   platform. Clockwork remains the variant caller, other components differ.
5. VARIANTS.parquet fell from 2.9 GB to 1.3 GB. MUTATIONS.parquet rose from 548.3 MB to
   1.1 GB. Both need chunked reading.
6. ENA_LOOKUP.parquet no longer exists. WGS_SAMPLES.parquet replaces it.
7. Total record size fell from 3.5 GB to 2.6 GB.

CRyPTIC state explicitly that they do not recommend using v2.1.2 or v3.0.0.

## Schema: the four table groups

The schema divides the tables into four groups.

1. Genetics: WGS_SAMPLES, GENOMES, VARIANTS, MUTATIONS, EFFECTS, PREDICTIONS
2. UKMYC phenotypes only: UKMYC_PLATES, UKMYC_GROWTH, UKMYC_PHENOTYPES, BASHTHEBUG,
   BASHTHEBUG_CLASSIFICATIONS
3. All phenotypes: DST_MEASUREMENTS, DST_SAMPLES. A superset of group 2, including R/S
   results derived using CRyPTIC ECOFFs plus any other phenotypic method (MGIT960, LJ,
   MYCOTB and others). Some samples carry multiple results for the same drug by different
   methods.
4. Reference: SITES, DRUG_CODE, COUNTRIES_LOOKUP, PLATE_LAYOUT

## Schema question 1: what each table holds

WGS_SAMPLES: one row per sample downloaded and processed through EIT Pathogena. Columns
shown: UNIQUEID, study_accession, sample_accession, run_accession, center_name, country,
location, first_public, fastq_ftp, fastq_md5, fastq_bytes. The schema marks this table with
an asterisk noting it also includes further columns for tracking progress through EIT
Pathogena, which are not drawn. The release notes confirm two of those: DATASET and a
processing status field.

GENOMES: per-sample genome metadata. UNIQUEID, SPECIES, N_LINEAGES, LINEAGE, SUBLINEAGE,
MYCOBACTERIAL_READS, TB_READS, TB_COVERAGE, TB_DEPTH, ANTIBIOGRAM, PIPELINE_BUILD.

VARIANTS: per-sample genomic variation in nucleotide space. Index UNIQUEID, GENE, VARIANT.
Then GENOME_POSITION, CODON_IDX, NUCLEOTIDE_INDEX, INDEL_LENGTH, INDEL_NUCLEOTIDES,
VCF_IDX, IS_NULL, IS_MINOR, MINOR_VARIANT, MINOR_READS, COVERAGE.

MUTATIONS: the same variation expressed in gene and protein terms. Index UNIQUEID, GENE,
MUTATION. Then GENE_POSITION, REF, ALT, NUCLEOTIDE_NUMBER, NUCLEOTIDE_INDEX, CODES_PROTEIN,
INDEL_LENGTH, INDEL_NUCLEOTIDES, AMINO_ACID_NUMBER, AMINO_ACID_SEQUENCE,
NUMBER_NUCLEOTIDE_CHANGES, IS_NULL, IS_MINOR, MINOR_MUTATION, MINOR_READS, COVERAGE, FRS.

EFFECTS: per-mutation prediction under a named catalogue. Index UNIQUEID, CATALOGUE_NAME,
CATALOGUE_VERSION, PREDICTION_VALUES, DRUG, GENE, MUTATION. Then PREDICTION and EVIDENCE
(json).

PREDICTIONS: per-sample, per-drug overall call under a named catalogue. Index UNIQUEID,
CATALOGUE_NAME, CATALOGUE_VERSION, CATALOGUE_VALUES, DRUG. Then PREDICTION.

UKMYC_PLATES: one row per plate reading. UNIQUEID, SITEID, SUBJID, LABID, ISOLATENO,
READINGDAY, BELONGS_GPI, PLATEDESIGN, TREE_PATH, IMAGEFILENAME, IMAGE_MD5SUM,
DUPLICATED_IMAGE, IM_IMAGE_DOWNLOADED, IM_IMAGE_FILTERED, IM_WELLS_IDENTIFIED,
IM_POS1GROWTH, IM_POS2GROWTH, IM_POS_AVERAGE, IM_DRUGS_INCONSISTENT_GROWTH,
TRUST_PHENOTYPES.

UKMYC_GROWTH: per-well growth measurement. Index UNIQUEID, READINGDAY, DRUG, DILUTION. Then
PLATEDESIGN, SITEID, WELL_CONC, GROWTH.

UKMYC_PHENOTYPES: the MIC table for UKMYC5/6 plates. Index UNIQUEID, DRUG. Then PLATEDESIGN,
BELONGS_GPI, SITEID, DILUTION, PHENOTYPE_QUALITY, READINGDAY, PRIMARY_DILUTION,
PRIMARY_METHOD, AMYGDA_DILUTION, BASHTHEBUG_DILUTION, BASHTHEBUGPRO_DILUTION,
PHENOTYPE_DESCRIPTION, BASHTHEBUG_NUMBER_CLASSIFICATIONS,
BASHTHEBUGPRO_NUMBER_CLASSIFICATIONS, MIC, LOG2MIC, BINARY_PHENOTYPE.

BASHTHEBUG: aggregated volunteer readings. Index UNIQUEID, READINGDAY, DRUG. Then
IMAGEFILENAME, BB_STATUS, DILUTION_DIFFERENCE, DILUTION_MAX, DILUTION_MEAN, DILUTION_MEDIAN,
DILUTION_MIN, DILUTION_STDEV, LOG2MIC_MEDIAN, MIC_MEDIAN, NUMBER_CANNOT_READ,
NUMBER_CLASSIFICATIONS, NUMBER_FAILED, NUMBER_VALID, NWELLS, PLATEDESIGN, PLATEIMAGE,
SUBJECT_SET.

BASHTHEBUG_CLASSIFICATIONS: the raw individual volunteer classifications, 4.75 million of
them. Index CLASSIFICATION_ID. Then USER_NAME, USER_ID, USER_IP, WORKFLOW_ID, WORKFLOW_NAME,
WORKFLOW_VERSION, CREATED_AT, SUBJECT_IDS, FILENAME, PLATE_IMAGE, PLATE_DESIGN, DRUG, PLATE,
STUDY_ID, READINGDAY, SITEID, BASHTHEBUG_DILUTION.

DST_MEASUREMENTS: all phenotypes from any method. Index UNIQUEID, DRUG. Then SOURCE,
METHOD_1, METHOD_2, METHOD_3, METHOD_CC, METHOD_MIC, PHENOTYPE, QUALITY.

DST_SAMPLES: UNIQUEID, COUNTRY_CODE, NUMBER_DST.

SITES: SITEID, COUNTRY, DESCRIPTION.

COUNTRIES_LOOKUP: COUNTRY_NAME, COUNTRY_CODE_3_LETTER, COUNTRY_CODE_2_LETTER,
COUNTRY_CODE_NUMERIC, LAT, LONG.

DRUG_CODE: DRUG_3_LETTER_CODE, DRUG_NAME.

PLATE_LAYOUT: index PLATEDESIGN, DRUG, DILUTION. Then CONC, ROW, COL, BINARY_PHENOTYPE. This
is where the concentration of each well and the ECOFF-derived R/S assignment live.

## Schema question 2: the join key

UNIQUEID is the sample identifier and the join key across every table in the Genetics, UKMYC
phenotypes and All phenotypes groups. The reference tables join on their own keys instead:
SITES on SITEID, COUNTRIES_LOOKUP on COUNTRY_NAME, DRUG_CODE on DRUG_3_LETTER_CODE,
PLATE_LAYOUT on PLATEDESIGN plus DRUG plus DILUTION.

The single exception inside the phenotype group is BASHTHEBUG_CLASSIFICATIONS, which has no
UNIQUEID. It is keyed on CLASSIFICATION_ID and joins through BASHTHEBUG or through the plate
image filename.

Join integrity has been fragile historically and must be checked empirically rather than
assumed. Two known past failures, both recorded as fixed: in v2.0.1 the UNIQUEID in the
genetics tables carried extra trailing text that prevented joining, and in v2.1.2 Site 07
(PHE) lab identifiers containing "/" and "." were converted to "_" in the phenotype tables
but not in the genetics tables, so those rows would not join.

## Schema question 3: where MICs live

Two tables, for two different purposes.

UKMYC_PHENOTYPES is the UKMYC5/6 MIC table and the substrate for MIC-level work. The drug is
in DRUG. The MIC is in MIC, with LOG2MIC as the numeric form and DILUTION as the well index.

DST_MEASUREMENTS is the superset covering every phenotypic method used, with the drug in
DRUG, an MIC in METHOD_MIC and a binary result in PHENOTYPE. Use it for coverage questions
and for non-UKMYC comparison, not as the primary MIC source.

## Schema question 4: drug codes

The DRUG_CODE table maps DRUG_3_LETTER_CODE to DRUG_NAME. Read from DRUG_CODES.csv.gz: the
codes are the conventional abbreviations, listed under "Drug codes" below.

## Schema question 5: VARIANTS versus MUTATIONS

VARIANTS records nucleotide-level change against the reference genome, positioned by
GENOME_POSITION and traceable to the source VCF row through VCF_IDX. MUTATIONS records the
interpreted, gene-relative consequence: amino acid number and sequence, whether the region
codes for protein, indel description, and the number of nucleotide changes involved. A
resistance catalogue is written against mutations, not raw variants, so MUTATIONS is the
table that joins to EFFECTS and onward to PREDICTIONS.

MUTATIONS also carries FRS, the fraction of reads supporting the call, which VARIANTS does
not. FRS is the column that matters for heteroresistance, which the project's own risk list
names as a known pitfall.

Both tables were generated by gnomonicus using gumpy and piezo. Minor alleles in known
resistance genes are deliberately reported by overriding the MIN_FRS filter in the Clockwork
VCF, so both tables contain a large number of minor alleles flagged through IS_MINOR.

## Schema question 6: WGS_SAMPLES

It holds the ENA study, sample and run accessions, the FTP paths, MD5s and byte sizes of the
FASTQ files, plus centre, country, location and first public date. Undrawn columns track
progress through EIT Pathogena. Two of these are named in the release notes and matter:

DATASET, taking the values CRyPTIC-v1.0 (37,984 samples), CRyPTIC-v2.0 (5,737),
CRyPTIC-v3.0 (8,896) and unknown (25). CRyPTIC-v1.0 is Release One, the set handed to
FIND and Seq&Treat to build the first WHO catalogue. CRyPTIC-v2.0 is samples that gained
WGS, pDST or both after the April 2020 data freeze, including about 1,100 from NICD
enriched for bedaquiline resistance. CRyPTIC-v3.0 is the set identified as a genuine
validation dataset.

A processing status field taking the values complete, cannot assemble, cannot speciate and
not uploaded. Completed counts by dataset: CRyPTIC-v1.0 37,893; CRyPTIC-v2.0 5,523;
CRyPTIC-v3.0 5,827; unknown 21.

## Schema question 7: phenotype quality

UKMYC_PHENOTYPES.PHENOTYPE_QUALITY is the per-MIC quality field. Its three values, defined in
the release notes: HIGH means the measurement was checked by at least two independent
methods that agreed; MEDIUM means only one measurement method was used; LOW means at least
two independent methods were used and none agreed, in which case the MIC reported is the one
read by the laboratory scientist.

DST_MEASUREMENTS.QUALITY carries the same three-level scheme across all phenotype sources.

UKMYC_PLATES.TRUST_PHENOTYPES is a plate-level boolean, sitting alongside the image QC
columns IM_IMAGE_DOWNLOADED, IM_IMAGE_FILTERED, IM_WELLS_IDENTIFIED, IM_POS1GROWTH,
IM_POS2GROWTH, IM_POS_AVERAGE and IM_DRUGS_INCONSISTENT_GROWTH.

The release notes report the HIGH proportion under TMAS as 88.7%, up from 79.2% under
AMyGDA. Computed from the table shipped here it is 78.6% to 80.7% depending on the
denominator, so both figures are recorded as reported rather than confirmed; see the
section on the high-confidence proportion below.

## Schema question 8: off-scale and censored MICs

MIC in UKMYC_PHENOTYPES is typed as categorical, not float, while LOG2MIC alongside it is
float. A categorical MIC field is how the censoring operators are preserved: values carry
prefixes rather than being plain numbers. The release notes confirm this indirectly by
recording a data cleaning pass in DST_MEASUREMENTS that corrected entries such as ">0,12"
into ">=0.12" so they could be joined to PLATE_LAYOUT.

The measurable range per drug is defined by PLATE_LAYOUT, which gives CONC for every
DILUTION of every PLATEDESIGN. The lowest and highest tested dilutions for a drug are where
left and right censoring occur.

PLATE_LAYOUT.BINARY_PHENOTYPE holds the ECOFF-derived R or S call per dilution. The
intermediate category was removed from UKMYC5 and UKMYC6 and the ECOFF applied as defined in
the ECOFF paper, so no intermediate results appear for UKMYC plates in DST_MEASUREMENTS.

Read from UKMYC_PHENOTYPES.parquet. The censoring strings are given under "MIC censoring
convention" below. LOG2MIC is populated on every row that carries a MIC, censored or not,
and carries no marker of censoring, so it cannot distinguish a censored reading from an
interior one. On right-censored rows it is one doubling above the concentration written in
MIC: ">2" carries 2.00 where log2 of 2 is 1. Checked across all 288,904 rows. For interior
and left-censored rows it is log2 of the stated concentration, within the 0.002 that
CRyPTIC's rounded concentration labels introduce. Every censoring decision in this project
therefore parses the MIC string and no step reads LOG2MIC as a measurement.

## Schema question 9: the mutation table's own conventions

Four conventions decide how a row of MUTATIONS.parquet or PLATE_LAYOUT.parquet is read, and
none is stated in the schema document. Each was checked against the 108,773 rows the four
target genes carry.

CODES_PROTEIN describes the mutation rather than the gene. It is True on exactly the rows
whose position is positive: 106,549 SNPs and 1,456 indels inside a coding sequence, against
731 SNPs, 22 indels and 15 gene deletions whose position is negative or absent. For a
protein-coding gene such as Rv0678 the flag is therefore False on every promoter row, so a
parser that reads it as the gene's nature treats those rows as though the gene were rRNA.
This project derives a gene-level flag, GENE_CODES_PROTEIN, true where any mutation in the
gene is placed in its coding sequence, and the grammar parser reads that instead.

An indel's position is the first base it affects, so one written from the promoter can reach
the gene. Twenty-two indel rows across the four genes start at a negative position. Fourteen
stay inside the promoter. Eight reach the coding sequence, all in Rv0678, removing 438 or 483
bases from position -3 or -43. A resolved minor allele can do the same: one sample carries
-3_del_cttgtgag, eight bases from -3, five of them coding. What decides whether such a
deletion shifts the reading frame is the number of bases it takes from the gene, not the sign
of its position.

A large deletion is written twice. Fifteen rows across the four genes carry one deletion as
both del_<fraction> and the explicit sequence removed, thirteen in Rv0678 and two in mmpL5,
with explicit sizes from 282 to 2,829 bases. Every carrier holds exactly one row of each
kind, so counting rows without collapsing them counts one event as two variants, which moves
a sample whose only finding is a gene deletion into a multiple-variant group. This project
marks the explicit row DOUBLE_REPORTED and keeps the fraction row, which carries the deleted
fraction.

PLATE_LAYOUT writes two of its rows per drug and design with operators. The lowest well
reads <=x, and a further row reads >x for the bin above the highest well, repeating that
well's concentration with its own S or R label. Both numbers are tested concentrations, so
the ladder is recovered by stripping the operators and deduplicating, while the R on the bin
row belongs to the bin rather than to the well it repeats.

## Where each published figure is computed

Every figure in docs/RESULTS.md is written by one of these modules into the report named
beside it. The reports are regenerated from the source tables and are not committed, because
outputs/ is ignored.

| Module | Report | What it carries |
| --- | --- | --- |
| code/inspect_mutations.py | mutations_inspection.txt | the mutation table's shape and memory |
| code/check_parsing.py | parsing_check.txt | grammar coverage over six genes |
| code/build_sample_status.py | sample_status_report.txt | per-genome genotype status |
| code/audit_cohort.py | audit_report.txt | join integrity, patient replication, phenotype quality, the discovery and validation halves, rows carrying no MIC, and what the plates can measure |
| code/analyse_groups.py | group_analysis_report.txt | every group against the reference group, the mmpL5 covariate, the gene deletions, stratified estimates, the homogeneity test over the principal lineages, and loss of function against substitution |
| code/model_effects.py | model_report.txt | the lineage effect with site held constant |
| code/cluster_adjust.py | cluster_report.txt | cluster structure, the three treatments of clonality, and the collapsed estimate under several seeds |
| code/mic_model.py | mic_model_report.txt, mic_estimates.csv | the tested ladders, validation against planted parameters, what the excluded rows could do, and the fitted distributions |
| code/heteroresistance.py | heteroresistance_report.txt, heteroresistance_estimates.csv, multi_allele_counts.csv | the uncertain group, resolved minor alleles, and the shifts with site held constant |
| code/build_atlas.py | atlas_report.txt, atlas_evidence.csv | the per-variant layer |
| code/unexplained.py | unexplained_report.txt, unexplained_counts.csv, unexplained_gene_sets.csv | resistance the genotype does not explain |
| code/prediction_metrics.py | prediction_metrics.txt, prediction_metrics.csv, prediction_thresholds.csv | genotype rules as tests, and the threshold sweep |
| code/discovery.py | discovery_report.txt, docs/PRE_REGISTRATION.md | the discovery half and the registered predictions |
| code/validate.py | validation_report.txt | the held-out test |

## Open items raised by the schema and release notes

1. The schema lists AMYGDA_DILUTION in UKMYC_PHENOTYPES and shows no TMAS column, although
   the release notes state TMAS replaced AMyGDA in this version. Resolved against the file:
   a TMAS_DILUTION column exists and the schema does not draw it. PRIMARY_METHOD names the
   laboratory scientist's reading method rather than an automated one, taking VZ and MB.
2. GENOMES carries SPECIES, LINEAGE, SUBLINEAGE and N_LINEAGES. The v2.0.0 notes stated
   lineage was missing because no samples had been run through Mykrobe. The Pathogena
   pipeline used from Release Three onwards does speciate. Resolved against the file: both
   columns are fully populated, with sublineage resolved to four levels, so no separate
   TB-Profiler or Mykrobe run is needed. Counts are under "GENOMES" below.
3. The 425 new National University of Singapore samples have phenotypes with no matching
   genetics, listed as future work.
4. Samples added in groups 6 and 7 have no COUNTRY_CODE in DST_SAMPLES.
5. 2,909 samples from CRyPTIC-v3.0 were never uploaded, most likely because they have three
   FASTQ files in the ENA.

## Licence handling and quarantine

CRyPTIC data is CC BY 4.0 and may be adapted and re-released with attribution.

WHO catalogue content is CC BY-NC-SA 3.0 IGO. It is never merged into released core data or
code, and is pulled at run time inside a separate benchmark module only.

EFFECTS and PREDICTIONS are quarantined under the same rule. Although CRyPTIC distribute
them under CC BY 4.0, their content is the WHO catalogue's classifications applied to
CRyPTIC samples: the CATALOGUE_NAME and CATALOGUE_VERSION columns hold WHO1
(WHO-UCN-GTB-PCI-2021.7) and WHO2 (WHO-UCN-GTB-PCI-2023.5). Re-releasing those columns
inside a core dataset would carry WHO catalogue expression into the output. Use them for
benchmarking, keep them out of the released core tables.

BASHTHEBUG_CLASSIFICATIONS is quarantined for a separate reason. It contains USER_NAME,
USER_ID and USER_IP for the volunteers, which is personal data. It is never redistributed,
in whole or in part, and no derived table retains any of those three columns.

No figure computed from EFFECTS, PREDICTIONS or GENOMES.ANTIBIOGRAM is reproduced anywhere
in this repository. Every count, shift and interval published here derives from mutation
strings in MUTATIONS.parquet, from UKMYC_PHENOTYPES.parquet, from PLATE_LAYOUT.parquet, from
GENOMES.parquet excluding its ANTIBIOGRAM column, and from the reference lookups.

## Interpretive choices already baked into CRyPTIC's predictions

These affect any benchmark computed against PREDICTIONS or EFFECTS.

1. An epistasis rule is applied: a loss of function mutation in mmpL5 overrides any
   resistance-associated mutation in Rv0678, producing an overall prediction of Susceptible.
2. All genetic variation in mmpL5 is reported as Susceptible rather than Unknown, on the
   reasoning that mmpL5 has no variants associated with resistance and is only tracked so
   the epistasis rule can be applied.
3. Minor alleles in resistance genes that are not associated with resistance in WHO2 are
   reported as Susceptible rather than Unknown.
4. CRyPTIC state their interpretation of the WHO2 catalogue is their own and was not
   publicly available at the time the note was written.

## Method papers

CRyPTIC data compendium, PLOS Biology 2022, 20(8):e3001721,
doi 10.1371/journal.pbio.3001721

CRyPTIC ECOFF paper, European Respiratory Journal 2022, 60:2200239,
doi 10.1183/13993003.00239-2022

TMAS, the plate reading model used in v3.4.0: Vo HT, Nguyen S, Tran AT, Nguyen H, Ho Bich H,
Fowler PW, Walker TM, Nguyen TT. Deep learning-based framework for Mycobacterium
tuberculosis bacterial growth detection for antimicrobial susceptibility testing.
Computational and Structural Biotechnology Journal, 2025, doi 10.1016/j.csbj.2025.05.030.
Preprint doi 10.1101/2025.02.14.638231, posted 19 February 2025. Trained on 4,018 CRyPTIC
plate images, benchmarked against AMyGDA and BashTheBug, essential agreement 98.8%.

The compendium and ECOFF papers describe the v1.1.1-era data read by AMyGDA. Where they
differ from v3.4.0, v3.4.0 governs the data and the papers govern the method rationale.

## Related records

- v3.3.0, 24 April 2025, added WGS for 4,832 further samples
- v3.2.0, 8 April 2025, added 9,160 samples and 54,910 pDST measurements
- v3.1.0, 28 February 2025, first Release Three version, WGS_SAMPLES introduced
- v3.0.0, 28 January 2025, https://zenodo.org/records/16041005, not recommended by CRyPTIC
- v2.1.2, 23 July 2024, DOI 10.5281/zenodo.15679886, not recommended by CRyPTIC
- v1.1.1, 25 January 2021, the release used in the primary CRyPTIC publications

Release One tables are mirrored at EMBL-EBI:
https://ftp.ebi.ac.uk/pub/databases/cryptic/release_june2022/reproducibility/data_tables/cryptic-analysis-group/

Note: the CRyPTIC project website lists DOI 10.5281/zenodo.15680920 against both v3.4.0 and
v1.1.1, which cannot both be correct. Version numbers and DOIs here were read off Zenodo,
not off that page.

## Documentation read

DATA_SCHEMA.pdf: read
RELEASE_NOTES.md: read, all versions back to Release One

## Cumulative project spend

USD 0.00 against the USD 50 to 60 ceiling shared across all three trajectories.

# Findings computed from the data (v3.4.0)

Computed from DRUG_CODES.csv.gz, PLATE_LAYOUT.parquet, UKMYC_PHENOTYPES.parquet and
GENOMES.parquet as downloaded from the v3.4.0 record. All figures below are reproducible
from those four files.

## Drug codes (question 4, now answered)

The DRUG_CODE table holds 39 codes. The four of interest:
BDQ bedaquiline, CFZ clofazimine, DLM delamanid, LZD linezolid.
Also present and relevant to the wider new and repurposed group: PAN pretomanid,
DCS D-cycloserine, CYC cycloserine, SZD sutezolid.

## UKMYC_PHENOTYPES: actual shape

288,904 rows, indexed on UNIQUEID and DRUG. 21,685 unique samples across 14 drugs:
AMI, BDQ, CFZ, DLM, EMB, ETH, INH, KAN, LEV, LZD, MXF, PAS, RFB, RIF.
Plate designs: UKMYC6 190,372 rows, UKMYC5 98,532.
Reading day: 14 (268,545), 21 (19,239), 10 (1,120).

Rows, samples and rows carrying a MIC, per drug:

| Drug | Rows | Samples | With a MIC | UKMYC5 | UKMYC6 |
| --- | --- | --- | --- | --- | --- |
| AMI | 21,685 | 21,685 | 21,475 | 7,041 | 14,644 |
| BDQ | 21,681 | 21,681 | 21,437 | 7,037 | 14,644 |
| CFZ | 21,683 | 21,683 | 21,474 | 7,039 | 14,644 |
| DLM | 21,683 | 21,683 | 21,367 | 7,039 | 14,644 |
| EMB | 21,681 | 21,681 | 21,485 | 7,037 | 14,644 |
| ETH | 21,681 | 21,681 | 21,467 | 7,037 | 14,644 |
| INH | 21,682 | 21,682 | 21,316 | 7,038 | 14,644 |
| KAN | 21,684 | 21,684 | 21,461 | 7,040 | 14,644 |
| LEV | 21,681 | 21,681 | 21,519 | 7,037 | 14,644 |
| LZD | 21,681 | 21,681 | 21,520 | 7,037 | 14,644 |
| MXF | 21,681 | 21,681 | 21,516 | 7,037 | 14,644 |
| PAS | 7,038 | 7,038 | 6,882 | 7,038 | 0 |
| RFB | 21,682 | 21,682 | 21,480 | 7,038 | 14,644 |
| RIF | 21,681 | 21,681 | 21,358 | 7,037 | 14,644 |

One row per sample per drug throughout. PAS is carried on UKMYC5 only, which is why its
count is the UKMYC5 total. 3,147 rows carry no MIC, and LOG2MIC is absent on exactly those
rows.

The 21,685 figure supersedes the 21,570 quoted for v2.1.2. Note the count of drugs is 14,
not the 13 quoted in the published papers, because UKMYC5 and UKMYC6 differ in composition.

## The schema PDF is out of date for this table

The actual file contains a TMAS_DILUTION column, which the schema does not draw. The
BASHTHEBUGPRO_DILUTION and BASHTHEBUGPRO_NUMBER_CLASSIFICATIONS columns drawn in the schema
do not exist in the file. Actual columns, in order: PLATEDESIGN, BELONGS_GPI, SITEID,
DILUTION, PHENOTYPE_QUALITY, READINGDAY, PRIMARY_DILUTION, PRIMARY_METHOD, AMYGDA_DILUTION,
BASHTHEBUG_DILUTION, TMAS_DILUTION, PHENOTYPE_DESCRIPTION,
BASHTHEBUG_NUMBER_CLASSIFICATIONS, MIC, LOG2MIC, BINARY_PHENOTYPE.

Where the schema and the data disagree, the data governs.

## The three reading methods

PHENOTYPE_DESCRIPTION records which methods were compared and whether they agreed. Values
and counts: VZ,TM AGREE 195,939; ALL DISAGREE 28,822; VZ ONLY 27,125; BB,TM AGREE 20,662;
VZ,BB AGREE 10,393; BB RUNNING 5,963. VZ is the Vizion reading by the laboratory scientist,
BB is BashTheBug, TM is TMAS.

PRIMARY_METHOD takes two values: VZ 287,292 and MB 1,612. Neither is defined in the schema
or the release notes. The plate-validation paper reads every MIC by three methods, naming
them as the Vizion digital viewing system, a mirrored box and an inverted-light microscope,
and abbreviates the second MB: Rancoita PM, Cugnata F, Cruz AL, Borroni E, Hoosdally SJ,
Walker TM, Grazian C, Davies TJ, Peto TEA, Crook DW, Fowler PW, Cirillo DM, for the CRyPTIC
Consortium. Validating a 14-drug microtiter plate containing bedaquiline and delamanid for
large-scale research susceptibility testing of Mycobacterium tuberculosis. Antimicrobial
Agents and Chemotherapy 2018, 62(9):e00344-18, doi 10.1128/AAC.00344-18. On that reading VZ
is the Vizion and MB the mirrored box, both scientist readings rather than automated ones.
Attributed to that paper, not confirmed by CRyPTIC documentation for this column.

## High-confidence proportion does not reconcile

The v3.4.0 release notes state TMAS raised the proportion of high-confidence MICs from 79.2%
to 88.7%. Computed from this table, HIGH is 226,994 of 288,904 rows, 78.6%. Restricting to
rows carrying a MIC gives 79.1%. Excluding rows still marked BB RUNNING gives 80.2%, and
both restrictions together give 80.7%. None reaches 88.7%.
PLACEHOLDER: the denominator CRyPTIC used for 88.7% is not identified. The TMAS paper cited
under "Method papers" does not carry either figure, reporting 98.8% essential agreement
against ground truth instead, so the release notes remain the only source. Do not cite
88.7% without resolving it. Cite the figure computed from the version used instead.

## MIC censoring convention (question 8, now answered)

MIC is a string. Interior values are plain numbers. Left-censored values carry a "<=" prefix
at the lowest tested concentration. Right-censored values carry a ">" prefix at the highest
tested concentration. PLATE_LAYOUT confirms the same convention: dilution 1 carries a "<=X"
concentration, and the top dilution carries a ">Y" concentration with no ROW or COL, being
the off-plate bin.

The tested range differs by plate design for the same drug, so PLATEDESIGN must be part of
every join and every censoring decision. Never join on DILUTION alone.

| Drug | UKMYC5 range | UKMYC6 range | ECOFF boundary (both designs) |
| --- | --- | --- | --- |
| BDQ | <=0.015 to >2 | <=0.008 to >1 | S at or below 0.25, R at or above 0.5 |
| CFZ | <=0.06 to >4 | <=0.03 to >2 | S at or below 0.25, R at or above 0.5 |
| DLM | <=0.015 to >1 | <=0.008 to >0.5 | S at or below 0.12, R at or above 0.25 |
| LZD | <=0.03 to >2 | <=0.06 to >4 | S at or below 1.0, R at or above 2.0 |

The ECOFF sits at the same concentration on both designs even though the dilution index
differs. PLATE_LAYOUT.BINARY_PHENOTYPE carries the assignment. No intermediate category
appears.

## GENOMES

54,057 rows, one per sample. Every column is fully populated with no nulls.

SPECIES: M. tuberculosis 53,589; M. bovis 191; M. bovis subsp. BCG 98; M. africanum 90;
M. orgyis 73; M. caprae 8; M. microti 4; M. canettii 2; MTB Complex 2.

LINEAGE and SUBLINEAGE are fully populated. Lineage counts: lineage4 25,049; lineage2
16,419; lineage3 6,741; lineage1 4,570; mixed 736; then Bovis, BCG, lineage6, lineage5,
animal clades, Caprae, lineage7. Sublineage is resolved to four levels, for example
lineage4.10 4,599, lineage4.1.2.1 4,113, lineage2.2 4,096.

Consequence: lineage assignment does not require a separate TB-Profiler or Mykrobe run. That
step drops out of the plan.

TB_COVERAGE median 99.30%, TB_DEPTH median 90.0x.
PIPELINE_BUILD takes six values, so the cohort was not processed under a single build.

ANTIBIOGRAM is a 14-character R/S/U string, one character per drug. It is catalogue-derived
and therefore falls under the same quarantine as EFFECTS and PREDICTIONS. No figure in this
repository is computed from it.

The 54,057 row count does not match the 53,897 stated on the record as having both WGS and
pDST. PLACEHOLDER: the 160-sample difference is unexplained and cannot be resolved from the
files on disk. The record's figure counts samples with any pDST result, which lives in
DST_MEASUREMENTS, and that file downloaded as 92 bytes and does not open as parquet.
UKMYC_PHENOTYPES covers only the 96-well plate subset and cannot stand in for it.

## The join: 6,525 MIC samples have no genome

Of the 21,685 samples with UKMYC MICs, 15,160 appear in GENOMES and 6,525 do not. The loss
is not random. It is concentrated by site:

| Site | MIC samples | No genome | Percent lost |
| --- | --- | --- | --- |
| 12 | 425 | 425 | 100.0 |
| 16 | 147 | 147 | 100.0 |
| 02 | 3,709 | 2,622 | 70.7 |
| 03 | 2,285 | 1,019 | 44.6 |
| 17 | 124 | 45 | 36.3 |
| 06 | 2,637 | 637 | 24.2 |
| 10 | 2,537 | 592 | 23.3 |
| 05 | 4,349 | 650 | 14.9 |
| 08 | 2,467 | 237 | 9.6 |
| 04 | 1,581 | 100 | 6.3 |
| 14 | 369 | 23 | 6.2 |
| 01 | 136 | 5 | 3.7 |
| 11 | 457 | 14 | 3.1 |
| 20 | 462 | 9 | 1.9 |

Two sites are lost entirely and one loses seven tenths of its samples. The genome-matched
cohort is therefore geographically biased relative to the full MIC cohort, and any
prevalence figure computed on it inherits that bias. This must be stated in any output.
SITES.csv.gz has been read; the institution and country of every SITEID appear in the site
table below.

## Statistical power for the four target drugs

Restricted to the 15,160 samples that have both a UKMYC MIC and a genome.

| Drug | Rows | Resistant | Percent | HIGH-quality rows | Resistant at HIGH | Left-censored | Right-censored |
| --- | --- | --- | --- | --- | --- | --- | --- |
| BDQ | 15,156 | 171 | 1.13 | 11,652 | 141 | 13.7% | 27 |
| CFZ | 15,158 | 683 | 4.51 | 10,750 | 426 | 40.0% | 25 |
| DLM | 15,158 | 229 | 1.51 | 10,902 | 152 | 60.8% | 100 |
| LZD | 15,156 | 217 | 1.43 | 11,303 | 160 | 1.0% | 71 |

## The bedaquiline and clofazimine cross-resistance axis

Cross-tabulated on the genome-matched cohort:

| | CFZ R | CFZ S |
| --- | --- | --- |
| BDQ R | 88 | 81 |
| BDQ S | 594 | 14,136 |

763 samples are resistant to at least one of the pair, against 171 for bedaquiline alone.

Restricted to samples where both drugs carry a HIGH-quality phenotype: 58 dual-resistant,
28 bedaquiline-only, 310 clofazimine-only.

Resistance by lineage on the genome-matched cohort:

Bedaquiline: lineage3 2.21% (27 of 1,219), lineage2 1.43% (77 of 5,379), mixed 1.35%,
lineage1 0.78%, lineage4 0.77% (58 of 7,525).

Clofazimine: lineage2 6.38% (343 of 5,380), lineage3 5.58% (68 of 1,219), mixed 4.95%,
lineage4 3.15% (237 of 7,525), lineage1 2.85%.

Clofazimine-resistant lineage2 samples by sublineage: lineage2.2.1 113, lineage2.2.3 66,
lineage2.2 63, lineage2.2.7 46, lineage2.2.6 20, then smaller groups.

Lineage2 is the Beijing family, which supplies the regional link to Kazakhstan and Russia.
PLACEHOLDER: mapping from these sublineage labels to the named B0/W148 and Central Asia
Outbreak clones has not been established and must not be assumed.

# Site structure and its consequences

## The genome-match loss, by country

SITES.csv.gz maps SITEID to institution, city and country. The 6,525 MIC samples with no
genome resolve as follows:

| Site | Institution | Country | MIC samples | Matched | Lost | Percent lost |
| --- | --- | --- | --- | --- | --- | --- |
| 02 | Chinese Center for Disease Control and Prevention | China | 3,709 | 1,087 | 2,622 | 70.7 |
| 03 | Institute of Microbiology and Laboratory Medicine, Gauting | Germany | 2,285 | 1,266 | 1,019 | 44.6 |
| 05 | Universidad Peruana Cayetano Heredia | Peru | 4,349 | 3,699 | 650 | 14.9 |
| 06 | San Raffaele Scientific Institute | Italy | 2,637 | 2,000 | 637 | 24.2 |
| 10 | Centre for Tuberculosis, NICD | South Africa | 2,537 | 1,945 | 592 | 23.3 |
| 12 | National University of Singapore | Singapore | 425 | 0 | 425 | 100.0 |
| 08 | Oxford University Clinical Research Unit | Vietnam | 2,467 | 2,230 | 237 | 9.6 |
| 16 | TORCH, Stellenbosch | South Africa | 147 | 0 | 147 | 100.0 |
| 04 | Hinduja Hospital and FMR, Mumbai | India | 1,581 | 1,481 | 100 | 6.3 |
| 17 | African Health Research Institute, Durban | South Africa | 124 | 79 | 45 | 36.3 |
| 14 | Brazil, Sao Paolo | Brazil | 369 | 346 | 23 | 6.2 |
| 11 | Public Health Sweden | Sweden | 457 | 443 | 14 | 3.1 |
| 20 | University of Capetown | South Africa | 462 | 453 | 9 | 1.9 |
| 01 | Centers for Disease Control and Prevention | United States | 136 | 131 | 5 | 3.7 |

Site 12 losing exactly 425 samples confirms the release note stating that the redownload
picked up 425 National University of Singapore samples whose FASTQ files are not yet
processed.

The largest single loss is China at 2,622 samples. Among the Chinese samples that did match,
66.0% are lineage2, the second highest proportion of any site after Vietnam at 68.9%. The
loss therefore falls hardest on the lineage2 population, which is the population the
Beijing-related regional argument depends on.

## Resistance is concentrated in one site

Bedaquiline and clofazimine resistance across the matched cohort, by site:

| Site | Country | n | BDQ resistant | Percent | CFZ resistant | Percent |
| --- | --- | --- | --- | --- | --- | --- |
| 10 | South Africa (NICD) | 1,945 | 61 | 3.14 | 264 | 13.57 |
| 06 | Italy | 2,000 | 24 | 1.20 | 87 | 4.35 |
| 02 | China | 1,087 | 17 | 1.56 | 90 | 8.28 |
| 11 | Sweden | 443 | 16 | 3.61 | 6 | 1.35 |
| 14 | Brazil | 346 | 16 | 4.62 | 41 | 11.85 |
| 04 | India | 1,477 | 11 | 0.74 | 50 | 3.38 |
| 05 | Peru | 3,699 | 9 | 0.24 | 18 | 0.49 |
| 08 | Vietnam | 2,230 | 9 | 0.40 | 90 | 4.04 |
| 20 | South Africa (Capetown) | 453 | 4 | 0.88 | 25 | 5.52 |
| 03 | Germany | 1,266 | 2 | 0.16 | 4 | 0.32 |

NICD alone supplies 61 of 171 bedaquiline-resistant and 264 of 683 clofazimine-resistant
samples. The v3.1.0 release notes state that roughly 1,100 NICD samples were added
specifically because they are enriched for bedaquiline resistance. The concentration is
therefore a sampling decision rather than an epidemiological observation, and no prevalence
figure computed from this dataset can be read as a population estimate.
