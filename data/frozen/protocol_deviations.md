# Protocol deviations

No outcome-informed deviations at lock time.

- Stage 0 QC tracking columns `denoised`, `merged`, and `nonchim` were blank because of an index-name mismatch. Stage 1 reconstructs assigned genus depth from the immutable count table without altering Stage 0. This is a technical-record repair and does not change patients, outcomes, features, or thresholds.

- On 2026-09-29, the first GSE194331 Ensembl REST pass returned transient HTTP errors for five prespecified module genes (NFKB1, NLRP3, OCLN, OLFM4, PADI4). A uniform technical rescue was applied to all five failed symbols after each current GRCh38 Ensembl stable ID was verified through the same public lookup endpoint. Previously successful mappings were cached, all four predeclared modules were recomputed, and no gene was added, removed, or selected using expression direction or significance.

- On 2026-09-29, publication-figure QA identified clipped labels in Figures 2, 4 and 6 and a potentially misleading filled point for Borreliella agreement in Figure 5. Only display labels, margins, text wrapping and device settings were changed. Borreliella is now shown as an open point annotated `both absent`; all numerical inputs, analyses, ordering rules and conclusions are unchanged. A private Stage 1 R library was created for `svglite` so that SVG text remains editable. The dependency was built from CRAN source tarballs downloaded once on the control computer because direct CRAN access from the server failed; no public R library or scientific analysis environment was modified.
