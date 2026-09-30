# R environment

The final figure script was developed with R 4.6.1 and uses base graphics plus:

- `ragg` for optional PNG/TIFF export;
- `svglite` for optional editable SVG export; and
- `systemfonts` through the graphics backends.

Install the packages in a clean R library:

```r
install.packages(c("ragg", "svglite", "systemfonts"))
```

Then run:

```bash
Rscript scripts/figures/make_figures.R /absolute/path/to/repository
```

PDF export requires only a Cairo-capable R installation. When `svglite` or
`ragg` is absent, the corresponding optional formats are skipped without
changing the PDFs. The script reads only `data/published_results/` and writes
figures to `results/figures/`. PDF byte hashes may differ across operating systems because
of device metadata and font embedding; compare page count, dimensions, labels,
source-data tables, and rendered appearance rather than PDF hashes alone.
