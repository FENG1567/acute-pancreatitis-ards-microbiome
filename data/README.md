# Data directory

| Directory | Contents | Access route |
|---|---|---|
| `frozen/` | analysis protocol, signature lock, deviations, variable definitions, software records | included in repository |
| `metadata/` | patient/sample/run mappings, blinded QC, exclusions, population flow | included in repository; derived from public records |
| `processed/` | two compressed genus count matrices used by the public scripts | included in repository; derived from public sequencing accessions |
| `host_response/` | host-module scores, results, gene mapping and coverage | included as analysis-ready derived data |
| `published_results/` | frozen numerical tables used in the manuscript and figures | included in repository |
| `figure_source_data/` | source-data tables mapped to the six main figures | included in repository |

Raw FASTQ files and large intermediate objects are excluded. See
`docs/ACCESSIONS.md` for retrieval routes and `docs/DATA_DICTIONARY.md` for file
and variable definitions. Empty strings or `NA` denote unavailable public data;
they must not be inferred from outcome status.
