#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dada2)
  library(data.table)
})

root <- Sys.getenv("STAGE1_ROOT", unset="")
stage0 <- Sys.getenv("STAGE0_ROOT", unset="")
if (!nzchar(root) || !nzchar(stage0)) stop("Set STAGE1_ROOT and STAGE0_ROOT")
manifest_path <- file.path(root, "01_metadata/PRJNA428535_ena_download_manifest.tsv")
rawdir <- file.path(root, "02_raw/PRJNA428535")
outdir <- file.path(root, "03_processed/PRJNA428535")
filterdir <- file.path(outdir, "filtered")
dir.create(filterdir, recursive = TRUE, showWarnings = FALSE)

meta <- fread(manifest_path)
stopifnot(nrow(meta) == 124L, uniqueN(meta$run_accession) == 124L)
setorder(meta, run_accession)
input <- file.path(rawdir, paste0(meta$run_accession, ".fastq.gz"))
if (!all(file.exists(input))) stop("not all FASTQ files exist")
names(input) <- meta$run_accession
filtered <- file.path(filterdir, paste0(meta$run_accession, ".fastq.gz"))
names(filtered) <- meta$run_accession
threads <- 32L

filtering <- filterAndTrim(
  input, filtered,
  truncLen = 0, minLen = 100, maxLen = 350, maxN = 0,
  maxEE = 2, truncQ = 2, rm.phix = TRUE,
  compress = TRUE, multithread = threads, verbose = TRUE
)
rownames(filtering) <- meta$run_accession

errors <- learnErrors(filtered, nbases = 1e8, multithread = threads,
                      randomize = TRUE, verbose = TRUE)
derep <- derepFastq(filtered, verbose = TRUE)
names(derep) <- meta$run_accession
denoised <- dada(derep, err = errors, multithread = threads, pool = FALSE,
                 BAND_SIZE = 32, HOMOPOLYMER_GAP_PENALTY = -1)
seqtab <- makeSequenceTable(denoised)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = threads, verbose = TRUE)

asv.seq <- colnames(seqtab.nochim)
asv.id <- paste0("ASV", seq_along(asv.seq))
colnames(seqtab.nochim) <- asv.id
taxonomy <- assignTaxonomy(
  asv.seq,
  file.path(stage0, "02_raw/reference/rdp_train_set_18.fa.gz"),
  minBoot = 50, tryRC = TRUE, multithread = threads, verbose = TRUE
)
rownames(taxonomy) <- asv.id
tax.df <- as.data.table(taxonomy, keep.rownames = "ASV")
tax.df[, sequence := asv.seq]
count.df <- as.data.table(seqtab.nochim, keep.rownames = "run_accession")
fwrite(count.df, file.path(outdir, "asv_counts.tsv.gz"), sep = "\t")
fwrite(tax.df, file.path(outdir, "asv_taxonomy.tsv.gz"), sep = "\t")

genus <- taxonomy[, "Genus"]
family <- taxonomy[, "Family"]
genus[is.na(genus) | genus == ""] <- paste0("Unclassified_", family[is.na(genus) | genus == ""])
genus[is.na(genus) | genus == "" | genus == "Unclassified_NA"] <- "Unclassified"
genus.counts <- t(rowsum(t(seqtab.nochim), group = genus, reorder = TRUE))
fwrite(as.data.table(genus.counts, keep.rownames = "run_accession"),
       file.path(outdir, "genus_counts.tsv.gz"), sep = "\t")

denoised.counts <- vapply(denoised, function(x) sum(getUniques(x)), numeric(1))
tracking <- data.table(
  run_accession = meta$run_accession,
  input = as.numeric(filtering[meta$run_accession, "reads.in"]),
  filtered = as.numeric(filtering[meta$run_accession, "reads.out"]),
  denoised = as.numeric(denoised.counts[meta$run_accession]),
  nonchim = as.numeric(rowSums(seqtab.nochim)[meta$run_accession])
)
fwrite(tracking, file.path(root, "04_qc/PRJNA428535_dada2_tracking.tsv"), sep = "\t")
saveRDS(list(sequence_table = seqtab.nochim, taxonomy = tax.df,
             run_metadata = meta, tracking = tracking, errors = errors),
        file.path(outdir, "dada2_result.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
           file.path(root, "logs/PRJNA428535_DADA2_COMPLETE"))
