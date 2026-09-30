# Reproduction validation record

Validation date: 2026-09-29

## Analysis-ready workflow

Command:

```text
python run_all.py --no-figures
```

Result: PASS.

- primary cohort: 65 AP patients (26 ARDS, 39 non-ARDS);
- primary effect: 0.594743654014341 SD;
- stratified-bootstrap 95% CI: 0.114417893863355 to 1.065313620954965;
- primary one-sided permutation P: 0.009299070092990702;
- age/sex/BMI-adjusted Firth OR: 1.846685973176114;
- leave-Klebsiella-out sensitivity: 0.324354037778469 SD,
  P=0.098390160983902;
- paired blood-neutrophil Spearman rho: 0.101687232435155;
- 19 regenerated result tables and two frozen host-response tables compared;
- five repository tests passed;
- no comparison problem remained.

The Firth optimizer coefficients agreed within absolute tolerance 1e-6 because
last optimizer digits can differ across SciPy builds. The primary effect,
confidence interval, and permutation P value use the stricter 1e-8/1e-9 table
tolerance and explicit value checks.

## Figure workflow

The R figure script was run on R 4.6.1 against a temporary copy of the frozen
publication tables. It completed with `REVISION_FIGURES_COMPLETE` and generated
six non-empty PDF files. The server did not have Arial and substituted Liberation
Sans; this is an environment-specific font warning, not a numerical difference.
PDF hashes are therefore not used as a cross-platform scientific equality test.

The script was made dependency-tolerant: Cairo PDF output is always attempted;
SVG is added when `svglite` is installed, and PNG/TIFF are added when `ragg` is
installed. Missing optional raster packages no longer prevent PDF reproduction.
