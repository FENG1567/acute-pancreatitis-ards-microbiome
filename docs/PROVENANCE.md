# Provenance and analytic separation

1. PRJNA1031835 and PRJNA771396 were analyzed before the external outcome test.
   PRJNA771396 produced the checksum-locked five-genus signature recorded in
   `data/frozen/stage0_locked_signature.tsv`.
2. The locked score was transferred without outcome-driven refitting to
   PRJNA893348. The primary statistical unit was the patient.
3. PRJNA428535 was analyzed separately as paired blood-neutrophil context. It is
   not an ARDS validation cohort.
4. GSE194331/PRJNA800337 was analyzed separately as an AP-severity host-response
   context. It is not patient-level multi-omics integration.
5. Candidate cohorts without defensible patient-level ARDS mapping were retained
   in eligibility/exclusion tables and were not relabelled.

The public-release copy standardizes machine-specific absolute paths to repository-
relative paths for portability. These packaging changes do not alter samples,
features, transformations, seeds, results, or
scientific interpretations. `protocol_lock_manifest.tsv` retains the sizes and
hashes of the original locked artifacts; its path column identifies the public
repository counterpart.
