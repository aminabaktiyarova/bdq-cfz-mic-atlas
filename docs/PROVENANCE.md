# Trajectory B: data provenance record

Researcher: Amina Baktiyarova, Independent Researcher, ORCID 0009-0007-6265-6493

This file records exactly which data version, licence and documentation this project rests
on, and what the schema says. Every figure published from this project traces back to what
is written here. Facts below were read directly off the Zenodo record, DATA_SCHEMA.pdf and
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
PLACEHOLDER: number of v3.4.0 samples with UKMYC5/UKMYC6 MIC data, per drug, to be computed
from UKMYC_PHENOTYPES.

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

MD5 checksums for every file are shown on the Zenodo record and should be verified after
download.

Naming mismatch to watch: the schema calls the drug lookup table DRUG_CODE, the file on
Zenodo is DRUG_CODES.csv.gz.

## Differences from v2.1.2 that affect this project

1. Plate images in v3.4.0 were read by TMAS, a convolutional neural network, replacing
   AMyGDA. This raised the proportion of MICs held at high confidence (at least two of the
   three measurement methods agreeing) from 79.2% to 88.7%.
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

The DRUG_CODE table maps DRUG_3_LETTER_CODE to DRUG_NAME.
PLACEHOLDER: the actual three-letter codes for bedaquiline, clofazimine, delamanid and
linezolid have not been read from the file. Confirm from DRUG_CODES.csv.gz (385 bytes)
rather than assuming the conventional abbreviations.

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

Under TMAS the HIGH proportion is 88.7%, up from 79.2% under AMyGDA.

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

PLACEHOLDER: the exact set of strings used in MIC for censored values in v3.4.0, and whether
LOG2MIC is populated on censored rows. Read from UKMYC_PHENOTYPES.parquet before any
statistical handling is designed.

## Open items raised by the schema and release notes

1. The schema lists AMYGDA_DILUTION in UKMYC_PHENOTYPES and shows no TMAS column, although
   the release notes state TMAS replaced AMyGDA in this version. PRIMARY_METHOD is the
   likely place a TMAS reading is named. PLACEHOLDER: confirm which column holds the TMAS
   reading and what values PRIMARY_METHOD takes.
2. GENOMES carries SPECIES, LINEAGE, SUBLINEAGE and N_LINEAGES. The v2.0.0 notes stated
   lineage was missing because no samples had been run through Mykrobe. The Pathogena
   pipeline used from Release Three onwards does speciate. PLACEHOLDER: confirm whether
   LINEAGE and SUBLINEAGE are populated in v3.4.0. If they are, lineage assignment does not
   need a separate TB-Profiler or Mykrobe run.
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

PRIMARY_METHOD takes two values: VZ 287,292 and MB 1,612.
PLACEHOLDER: what MB denotes is not documented in the schema or release notes.

## High-confidence proportion does not reconcile

The v3.4.0 release notes state TMAS raised the proportion of high-confidence MICs from 79.2%
to 88.7%. Computed from this table, HIGH is 226,994 of 288,904 rows, 78.6%. Restricting to
rows carrying a MIC gives 79.1%. Excluding rows still marked BB RUNNING gives 80.2%, and
both restrictions together give 80.7%. None reaches 88.7%.
PLACEHOLDER: the denominator CRyPTIC used for 88.7% is not identified. Do not cite that
figure without resolving it. Cite the figure computed from the version used instead.

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
and therefore falls under the same quarantine as EFFECTS and PREDICTIONS.

The 54,057 row count does not match the 53,897 stated on the record as having both WGS and
pDST. PLACEHOLDER: the 160-sample difference is unexplained.

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
PLACEHOLDER: the country and description of each SITEID, from SITES.csv.gz, which has not
been read.

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

# Genotype findings computed from EFFECTS and PREDICTIONS (v3.4.0)

These two tables are catalogue-derived and quarantined from any released output. The figures
below were computed to establish feasibility.

## Structure

EFFECTS: 1,154,127 rows, 53,864 samples, 15 drugs. Predictions: S 917,989; R 127,866;
U 81,122; F 27,150.
PREDICTIONS: 810,615 rows, 54,041 samples. S 621,919; R 120,552; U 58,269; F 9,875.

Only one catalogue is present in v3.4.0: WHO-UCN-GTB-PCI-2023.5, version 2.0. The WHO first
edition results described in the v2.1.0 release notes are not carried in this version.

The EVIDENCE column holds JSON containing WHO catalogue text verbatim, including the fields
FINAL CONFIDENCE GRADING, INITIAL CONFIDENCE GRADING and WHO HGVS, alongside the catalogue's
own solo-sample counts. This is WHO catalogue expression reproduced in full and reinforces
the quarantine on these two tables.

## Genes carried per target drug

BDQ: mmpL5 (101,224 rows, all graded S by the reporting rule), Rv0678 (3,481 rows, 2,832
samples, 733 graded R), pepQ (1,617 rows, 1,584 samples, 17 R), atpE (250 rows, 249 samples,
16 R).

CFZ: mmpL5, Rv0678 (715 R, 1,325 U), pepQ (16 R).

DLM: fgd1 (34,311 samples), fbiC (5,898 samples, 554 R), Rv2983, fbiA, fbiB, ddn (1,201
samples, 93 R).

LZD: rrl (7,704 samples, 74 R), rplC (1,746 samples, 203 R, 80 F).

## Solo mutation analysis, bedaquiline pathway

Solo defined as a sample carrying exactly one mutation graded R or U for the drug across the
catalogue's genes, excluding mmpL5. Restricted to samples with both a genome and a UKMYC MIC.

614 solo samples with paired bedaquiline and clofazimine MICs; 392 with high-quality
phenotypes on both drugs. By gene:

| Gene | Solo samples | BDQ resistant | CFZ resistant | Median log2 BDQ MIC | Median log2 CFZ MIC |
| --- | --- | --- | --- | --- | --- |
| Rv0678 | 383 | 54 | 83 | -3.06 | -3.06 |
| pepQ | 217 | 4 | 17 | -4.06 | -3.06 |
| atpE | 14 | 2 | 0 | -4.56 | -4.06 |

Reference group carrying no R or U graded mutation for bedaquiline: 14,531 samples, median
log2 bedaquiline MIC -5.06, median log2 clofazimine MIC -4.06.

A solo Rv0678 mutation therefore shifts the median bedaquiline MIC by 2 log2 units, a
fourfold rise, and the median clofazimine MIC by 1 log2 unit, a twofold rise. Despite the
smaller clofazimine shift, more solo Rv0678 samples cross the clofazimine ECOFF (83) than the
bedaquiline ECOFF (54), because the clofazimine wild-type population sits closer to its
ECOFF. A binary catalogue cannot express this asymmetry.

atpE is confirmed as descriptive rather than inferential: 14 solo samples in the entire
matched cohort. pepQ is better represented than expected at 217 solo samples, and its
phenotype is consistent with a low-level effect: a 1 log2 shift in bedaquiline MIC with only
4 samples crossing the ECOFF.

## Rv0678 variant sparsity

383 solo Rv0678 samples carry 152 distinct mutations. 112 are seen once and 28 are seen two
to four times. Only nine mutations reach eight or more solo samples:

| Mutation | Solo samples | BDQ resistant | CFZ resistant | Median log2 BDQ | Median log2 CFZ |
| --- | --- | --- | --- | --- | --- |
| 192_ins_g | 52 | 2 | 2 | -6.06 | -4.06 |
| 141_ins_c | 39 | 11 | 17 | -2.00 | -2.00 |
| 138_ins_g | 18 | 4 | 9 | -2.00 | -1.50 |
| E55D | 18 | 0 | 2 | -4.06 | -3.06 |
| M146T | 18 | 0 | 3 | -3.06 | -3.06 |
| L40V | 14 | 0 | 0 | -6.06 | -4.06 |
| R90C | 12 | 1 | 3 | -2.00 | -3.06 |
| N98D | 9 | 0 | 0 | -4.06 | -4.06 |
| N4T | 8 | 0 | 0 | -4.56 | -3.06 |

Per-variant MIC distributions are computable for roughly this handful only. Any resource has
to work at the level of mutation classes and positions, with per-variant detail for the small
well-supported set.

Mutation class alone does not separate the phenotype. Indels: 158 samples, 36 distinct, 27
bedaquiline resistant, median log2 -3.06. Missense and other: 225 samples, 116 distinct, 27
bedaquiline resistant, median log2 -3.06. Two frameshift insertions illustrate the point:
192_ins_g in 52 samples is essentially phenotypically silent, while 141_ins_c in 39 samples
shifts both MICs by around 3 log2 units.

## Lineage structure of the Rv0678 signal

Solo Rv0678 samples by lineage:

| Lineage | Solo samples | BDQ resistant | CFZ resistant |
| --- | --- | --- | --- |
| lineage4 | 179 | 8 | 10 |
| lineage2 | 142 | 30 | 58 |
| lineage3 | 43 | 11 | 10 |
| lineage1 | 19 | 5 | 5 |

The same gene produces resistance at 21% in lineage2 against 4.5% in lineage4. This is
confounded by which mutations occur in which lineage and by site of origin, and is not a
causal claim. It is the observation that motivates a lineage-stratified analysis and it is
the concrete link to the Beijing-dominated regional epidemiology.

## Genotypically unexplained resistance

Counting resistant samples that carry no mutation graded R or U for that drug in any of the
catalogue's genes:

Bedaquiline: 111 of 171 resistant samples, 65%. Restricted to high-quality phenotypes, 85 of
141, 60%.
Clofazimine: 574 of 683 resistant samples, 84%. Restricted to high-quality phenotypes, 358 of
426, 84%.

The unexplained fraction persists at high phenotype quality, so it is not attributable to
measurement error alone. Other candidate explanations are heteroresistance below the
detection limit, mechanisms outside the catalogue's gene set, and genuine phenotypic
variation. This is the quantified gap that a new resource addresses.

## Linezolid solo mutations

Solo samples with a UKMYC linezolid MIC:

| Gene | Solo samples | Distinct mutations | Resistant | Median log2 MIC |
| --- | --- | --- | --- | --- |
| rrl | 998 | 270 | 25 | -1.00 |
| rplC | 193 | 41 | 36 | -1.00 |

rplC C154R alone accounts for 40 solo samples, 29 of them resistant, with a median log2 MIC
of 2.00 against a population median of -1.00: a shift of three doublings, from 0.5 to 4.0.
The next largest signal is rrl g2814t at 3 solo samples, all resistant.

The three most common rrl variants, g1052t (211 solo samples), g2399a (112) and c637g (95),
have median MICs identical to the population median and yield 4, 2 and 5 resistant samples
respectively. Almost all catalogued rrl variation is phenotypically silent in this dataset.

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

## The lineage2 signal is substantially confounded by site

Solo Rv0678 samples by site: NICD South Africa 101 (31 bedaquiline resistant), Peru 86 (1),
Italy 68 (10), Germany 34 (1), China 29 (7), India 29 (1), Vietnam 16 (1), then single
figures elsewhere. NICD alone contributes 31 of the 54 bedaquiline-resistant solo Rv0678
samples.

Stratifying the lineage2 against lineage4 comparison by site:

| Site | Country | lineage2 n | lineage2 BDQ R | lineage4 n | lineage4 BDQ R |
| --- | --- | --- | --- | --- | --- |
| 10 | South Africa | 64 | 22 | 32 | 5 |
| 02 | China | 24 | 6 | 4 | 0 |
| 03 | Germany | 22 | 1 | 10 | 0 |
| 08 | Vietnam | 12 | 1 | 2 | 0 |
| 04 | India | 9 | 0 | 3 | 0 |
| 06 | Italy | 1 | 0 | 33 | 0 |
| 05 | Peru | 1 | 0 | 85 | 1 |

The pooled contrast of 21% against 4.5% is driven substantially by which sites contribute
which lineages. Peru and Italy together supply 118 lineage4 solo samples with one resistant
sample between them, and both are near-zero-lineage2 sites with low resistance throughout.
NICD is the only site with adequate numbers of both, and there the contrast is 34% against
16%, roughly half the pooled gap.

Consequence: site must enter any model of the lineage effect as a covariate. The pooled
lineage figure must not be reported on its own, and the regional argument cannot rest on it.
