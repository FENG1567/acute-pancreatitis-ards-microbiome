# Reproducibility scope

## Tier 1: analysis-ready reproduction (fully included)

`python run_all.py` recomputes the primary PRJNA893348 score and inferential
statistics, sensitivity analyses, Firth models, component audits, and the
PRJNA428535 paired analyses. Fixed seeds and 10,000 resamples/permutations are
used exactly as recorded. Newly generated tables are compared numerically with
the frozen publication-candidate tables.

## Tier 2: figure reproduction (fully included; R required)

`scripts/figures/make_figures.R` regenerates six multi-panel figures from
`data/published_results/`. The exact source-data copies used by each figure are
also retained in `data/figure_source_data/`.

## Tier 3: raw-read reconstruction (code and accessions included; large inputs excluded)

Raw preprocessing scripts are provided for auditability, but raw FASTQ files,
large DADA2 intermediates, the RDP18 reference database, and the full host
expression matrix are not bundled. Set `STAGE0_ROOT` and `STAGE1_ROOT` and
reconstruct the expected workspace from the accessions in `ACCESSIONS.md` before
running the archival scripts. Exact raw rebuilding also requires DADA2 1.38.0
and the RDP18 training set with the checksum recorded in
`data/frozen/software_lock.tsv`.

## What is verified automatically

- cohort size: 65 primary AP patients, 26 ARDS and 39 non-ARDS;
- paired support: 62 blood-neutrophil pairs;
- host context: 87 AP transcriptomes;
- fixed five-genus score, pseudocount 0.5, sample-wise CLR, equal weights;
- primary effect, bootstrap interval, AUROC and primary P value;
- Firth estimates and paired-support results within numerical tolerance;
- separation of the primary P value from the secondary random-seed sensitivity;
- required files, checksums, and absence of machine-specific paths or credentials.

## Known limits

The public data do not provide sequenced negative controls for the two
microbiome cohorts, exact patient-level sampling-to-ARDS intervals, a second
independent AP-ARDS cohort, or same-patient microbiome-transcriptome integration.
Accordingly, this repository reproduces an association and triangulated support,
not a causal or prospective prediction claim.
