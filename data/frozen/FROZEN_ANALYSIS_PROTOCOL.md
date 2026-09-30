# Stage 1 Frozen Analysis Protocol

Locked: 2026-09-29T00:52:46+08:00  
Review: author-reviewed

## Confirmatory question

Test whether the Stage 0 outcome-blind, checksum-locked five-genus signature is
higher in the 26 AP patients who subsequently developed ARDS than in the 39 AP
patients who did not. Healthy controls are excluded. The patient is the only
statistical unit.

## Immutable signature

The source lock SHA256 is `0a10ecee753b8c31f00b216d2f457c6795abf5c09abee75714f665f09974b0af` and contains
Borreliella, Escherichia, Staphylococcus, Enterococcus, and Klebsiella, all in
the positive direction. The PRJNA893348 RDP18 table labels Escherichia as the
non-separable compound `Escherichia/Shigella`; this exact compound is the
pre-outcome nomenclature crosswalk. No Borreliella column or classified read is
present; the feature is retained as an all-zero count and no Borrelia or other
near-neighbour substitution is permitted.

For every retained sample, add 0.5 to every bacterial genus count, perform the
sample-wise CLR over the full genus table, take the arithmetic mean of the five
fixed positive-direction components (coefficient 0.2 each), then standardize
that aggregate by the sample standard deviation among all 65 frozen AP
patients. A higher score is prespecified as higher ARDS risk. No feature,
direction, weight, cut-off, transform, or patient can be selected using outcome
results.

## Population and timing

The frozen population is 65 unique AP patients (26 ARDS, 39 nonARDS). All pass
the outcome-blind minimum assigned-depth rule of 1,000 counts. The publication
supports sampling within 24 hours of AP onset and before subsequent ARDS at the
cohort level, but an individual time interval is not public; therefore primary
wording is “early association/identification,” not a patient-level prospective
prediction claim.

## Primary inference

The statistic is mean standardized score in ARDS minus nonARDS. The only
confirmatory P value is a one-sided patient-label permutation test with 10,000
permutations, alternative greater, master seed 20260929, and P<0.05. The
Monte Carlo standard error is reported. The main effect is the mean difference
with a 10,000-replicate stratified patient-bootstrap percentile 95% CI. AUROC
and its bootstrap interval are descriptive and cannot determine success.

## Support and sensitivity

Age, sex, and BMI are the only allowed adjustment variables and enter one
supportive complete-case logistic model without outcome-driven selection.
Critical diagnostics are leave-one-patient-out; log-depth adjustment; exclusion
of the blind bottom 5% depth; fixed pseudocount 1.0; every leave-one-genus-out
score without refitting; and a stringent skin-contaminant analysis omitting
Staphylococcus without redistributing its 0.2 weight. These cannot replace the
primary analysis. Single-genus tests are exploratory and use BH-FDR.

## Contamination and causal boundaries

No re-analyzable negative controls are public, so decontamination cannot be
claimed. Microbial reads do not prove viable organisms, migration direction, or
causality. Different platforms/regions remain separate. Cross-cohort host and
cellular findings are triangular support, never patient-level multi-omics.
