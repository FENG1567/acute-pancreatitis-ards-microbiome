# Data dictionary

## Core analysis inputs

### `data/metadata/subject_sample_run_map.tsv`

- `subject_id`: de-identified public participant alias; primary statistical unit.
- `sample_id`, `sample_accession`, `run_accession`: public sample/run identifiers.
- `ARDS_label_original`: public binary outcome used only in the primary cohort.
- `ARDS_definition`, `ARDS_assessment_window`: source-supported endpoint text.
- `collection_time_original`: cohort-level public sampling description.
- `include_primary`: frozen inclusion flag.
- `age`, `sex`, `bmi`: prespecified supportive covariates.
- `mapping_confidence`: strength of the public patient/sample mapping.

### `data/processed/PRJNA893348_genus_counts.tsv.gz`

Rows are public sample accessions; columns are unrarefied bacterial genus counts
after blind technical QC. `Escherichia/Shigella` is retained as the RDP18
compound label. Borreliella is a locked structural-zero feature inserted by the
analysis script and is not replaced by a near neighbour.

### `data/metadata/stage08_patient_sample_ledger.tsv`

`dataset`, `run_accession`, `patient_id`, `specimen`, and `group` identify public
samples across cohorts. For PRJNA428535, every retained patient has one `blood`
and one `neutrophils` record.

### `data/processed/PRJNA428535_genus_counts.tsv.gz`

Rows are run accessions and columns are genus counts for the 124 paired
PRJNA428535 samples.

## Primary derived variables

- `signature_score_sd`: mean of the five positive-direction CLR components,
  each weighted 0.2, standardized within the 65 frozen AP patients.
- `effect_SD`: ARDS minus non-ARDS mean score difference in SD units.
- `one_sided_permutation_p`: `(1 + exceedances)/(10,000 + 1)` under patient-label
  permutation with the prespecified greater-than alternative.
- `bootstrap_CI_low`, `bootstrap_CI_high`: patient-level percentile bootstrap
  interval from 10,000 replicates.
- `OR_firth`: odds ratio per 1-SD score from Firth penalized logistic regression.

## Paired derived variables

- `blood_score`, `neutrophil_score`: compartment-specific standardized relative
  scores; they are not absolute microbial loads.
- `spearman_rho`: within-patient rank association between compartments.
- `neutrophil_minus_blood_mean`: mean paired relative-score difference.
- detection columns: genus count greater than zero in the corresponding sample.

## Host-response variables

- `severity_ordinal`: mild=1, moderately severe=2, severe=3.
- four module columns: equal means of signed, within-gene z scores using the
  frozen module dictionary.
- `BH_FDR`: Benjamini-Hochberg-adjusted module-level P value.

For full headers and row-level values, inspect the TSV files directly. No hidden
recoding or unpublished clinical outcome imputation is used.
