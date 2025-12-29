#!/usr/bin/env Rscript

# ============================================================
# 08_host_deseq2.R (virus vs mock)
# - Forces condition levels: mock (reference) -> virus
# - Uses explicit contrast: virus vs mock
# - Host-only gene-level DE from Salmon quant.sf + tx2gene
# - No renv
# ============================================================

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(data.table)
})

# --------------------------- Paths ---------------------------
samples_csv  <- "data/metadata/samples.csv"
salmon_dir   <- "results/salmon"
tx2gene_file <- "ref/gencode/tx2gene_gencode_v49.tsv"
out_dir      <- "results/host_gene"

out_tsv <- file.path(out_dir, "deseq2_gene_results.tsv")
out_rds <- file.path(out_dir, "dds_gene.rds")

# --------------------------- Checks --------------------------
stopifnot(file.exists(samples_csv))
stopifnot(dir.exists(salmon_dir))
stopifnot(file.exists(tx2gene_file))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------- Read metadata ------------------------
meta <- readr::read_csv(samples_csv, show_col_types = FALSE)

# Require at least: sample_id, condition
req <- c("sample_id", "condition")
miss <- setdiff(req, names(meta))
if (length(miss) > 0) {
  stop("samples.csv missing columns: ", paste(miss, collapse = ", "),
       "\nAvailable columns: ", paste(names(meta), collapse = ", "))
}

meta <- meta %>%
  mutate(
    sample_id = as.character(sample_id),
    condition = tolower(as.character(condition))
  )

# Enforce expected conditions
expected <- c("mock", "virus")
present <- sort(unique(meta$condition))
if (!all(expected %in% present)) {
  stop("Condition values must include: mock and virus.\n",
       "Found: ", paste(present, collapse = ", "))
}

# Force factor levels: mock (reference) then virus
meta <- meta %>%
  mutate(
    condition = factor(condition, levels = c("mock", "virus"))
  )

# ---------------------- Quant.sf paths -----------------------
files <- file.path(salmon_dir, meta$sample_id, "quant.sf")
names(files) <- meta$sample_id

if (!all(file.exists(files))) {
  missing <- names(files)[!file.exists(files)]
  stop("Missing quant.sf for samples: ", paste(missing, collapse = ", "),
       "\nExpected: results/salmon/<sample_id>/quant.sf")
}

# ----------------------- Read tx2gene ------------------------
tx2gene <- readr::read_tsv(tx2gene_file, show_col_types = FALSE)

# Standardize columns
if (all(c("TXNAME", "GENEID") %in% names(tx2gene))) {
  tx2gene <- dplyr::rename(tx2gene, tx = TXNAME, gene = GENEID)
} else if (all(c("tx", "gene") %in% names(tx2gene))) {
  # ok
} else if (all(c("transcript_id", "gene_id") %in% names(tx2gene))) {
  tx2gene <- dplyr::rename(tx2gene, tx = transcript_id, gene = gene_id)
} else {
  stop("tx2gene columns not recognized.\nAvailable: ",
       paste(names(tx2gene), collapse = ", "))
}

tx2gene <- tx2gene %>%
  transmute(
    tx   = as.character(tx),
    gene = sub("\\..*$", "", as.character(gene))
  ) %>%
  distinct()

# -------------------------- tximport -------------------------
# Viral transcripts will be "missing from tx2gene" and ignored (expected).
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = TRUE)

# -------------------------- DESeq2 ---------------------------
dds <- DESeqDataSetFromTximport(
  txi,
  colData = as.data.frame(meta),
  design  = ~ condition
)

dds <- DESeq(dds)

# Explicit contrast: virus vs mock (log2FC > 0 means higher in virus)
res <- results(dds, contrast = c("condition", "virus", "mock"))

# ------------------------- Save outputs ----------------------
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::select(gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
  arrange(padj)

data.table::fwrite(res_df, out_tsv, sep = "\t")
saveRDS(dds, out_rds)

message("Contrast locked: virus vs mock (log2FC>0 = induced in virus)")
message("Saved: ", out_tsv)
message("Saved: ", out_rds)
message("Done.")
