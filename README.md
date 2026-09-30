# Cross-compartment microbiome signals and early ARDS in acute pancreatitis

This repository contains the analysis-ready data, frozen analysis definitions,
statistical code, figure source data, and verification utilities supporting the
manuscript on cross-compartment microbial signals and early ARDS in acute
pancreatitis (AP).

## What can be reproduced from this repository

The default workflow reproduces, without downloading FASTQ files:

1. the locked patient-level PRJNA893348 primary analysis (65 AP patients; 26
   ARDS and 39 non-ARDS);
2. prespecified and reviewer-requested sensitivity analyses, including Firth
   logistic regression;
3. PRJNA428535 paired blood-neutrophil analyses (62 patient pairs);
4. numerical checks against the publication-candidate result tables; and
5. the six publication figures from frozen source-data tables when R and the
   required graphics packages are available.

The confirmatory result is the one-sided 10,000-permutation test in
`results/analysis/primary_test.json` (P = 0.009299070092990702). The value
P = 0.008799120087991202 in `sensitivity_with_ci.tsv` is a separate,
fixed-seed sensitivity run and must not replace the primary P value.

## Quick start

Python 3.12 is recommended.

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python run_all.py
```

If R is not installed, the statistical analyses and validation still run; the
workflow reports that figure regeneration was skipped. To require figure
regeneration:

```bash
python run_all.py --require-figures
```

Run the repository tests separately with:

```bash
python -m unittest discover -s tests -v
```

## Repository layout

```text
data/frozen/              frozen protocol, signature, and software records
data/metadata/            patient/sample/run mappings and cohort flow
data/processed/           analysis-ready genus count matrices
data/host_response/       analysis-ready host-module results and metadata
data/published_results/   frozen publication-candidate numerical tables
data/figure_source_data/  source-data tables mapped to figure panels
scripts/analysis/         public statistical reproduction scripts
scripts/figures/          final multi-panel figure script
scripts/raw_preprocessing/archival raw-to-analysis workflow scripts
results/analysis/         regenerated analysis outputs and validation report
results/figures/          six final PDF figures
docs/                     scope, provenance, accessions, and environment notes
```

## Raw-data boundary

FASTQ files, DADA2 RDS objects, the RDP training database, full host-expression
matrices, and large temporary files are intentionally not included. They must
be retrieved from the public accessions listed in `docs/ACCESSIONS.md` and
processed with the scripts in `scripts/raw_preprocessing/`. Those scripts use
`STAGE0_ROOT` and `STAGE1_ROOT` environment variables rather than machine-
specific paths. The default repository workflow starts from the included
analysis-ready matrices because this is the smallest auditable package that can
recompute the central patient-level statistics.

## Scientific boundaries

- PRJNA893348 is the only public AP cohort in this project with a directly
  mappable patient-level ARDS label. It is the primary outcome cohort and is not
  its own independent replication.
- PRJNA428535 provides paired cellular-compartment context, not ARDS validation.
- GSE194331/PRJNA800337 provides independent AP-severity host-response context,
  not within-patient microbe-host integration.
- No public sequenced negative controls were available for PRJNA893348 or
  PRJNA428535. Microbial reads do not establish viability, migration direction,
  or causality.

## Data and code availability

Code in this repository is released under the MIT License. The public repository
is https://github.com/FENG1567/acute-pancreatitis-ards-microbiome, and the
publication-candidate release is tagged as `v1.0.1`. Public-source and derived
scientific data remain subject to the terms of their originating repositories
and publications; see `DATA_USE_NOTICE.md`. A manuscript-ready Code
Availability statement is provided in `CODE_AVAILABILITY.md`. Archiving the
tagged release in Zenodo remains recommended so that the exact version receives
a DOI.

See `docs/REPRODUCIBILITY_SCOPE.md` for the precise reproduction tiers and
known limitations.
