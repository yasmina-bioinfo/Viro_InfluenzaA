#!/usr/bin/env Rscript

# ============================================================
# Build and save a DESeq2 object (gene-level) for visualization
# from Salmon quant.sf files using tximport + tx2gene mapping.
#
# Robust features:
# - Works with multiple metadata column naming conventions:
#   sample / sample_id / Sample / run / sra_run ...
# - Works with multiple tx2gene column naming conventions:
#   tx/gene, TXNAME/GENEID, transcript_id/gene_id, etc.
# - Validates file paths and quant.sf existence.
#
# Outputs:
#   results/host_gene/dds_gene.rds
#   results/host_gene/vsd_gene.rds
# ============================================================

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(readr)
  library(dplyr)
})

# ----------------------------
# 0) Paths (adapt once per project)
# ----------------------------
meta_file    <- "data/metadata/samples.csv"
salmon_dir   <- "results/salmon"
tx2gene_file <- "ref/gencode/tx2gene_gencode_v49.tsv"

out_dir <- "results/host_gene"
out_dds <- file.path(out_dir, "dds_gene.rds")
out_vsd <- file.path(out_dir, "vsd_gene.rds")

# ----------------------------
# 1) Basic checks
# ----------------------------
stopifnot(file.exists(meta_file))
stopifnot(file.exists(tx2gene_file))
stopifnot(dir.exists(salmon_dir))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("[1/6] Reading metadata: ", meta_file)
meta <- read_csv(meta_file, show_col_types = FALSE)

# ----------------------------
# 2) Robust metadata handling
#    - We need a sample identifier column and a condition column
# ----------------------------
# Candidate names for sample id
sample_candidates <- c("sample", "sample_id", "Sample", "SampleID",
                       "run", "Run", "sra_run", "SRA", "accession", "Accession")

cond_candidates   <- c("condition", "Condition", "group", "Group", "treatment", "Treatment")

pick_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

sample_col <- pick_col(meta, sample_candidates)
cond_col   <- pick_col(meta, cond_candidates)

if (is.na(sample_col)) {
  stop("Could not find a sample ID column in metadata.\nAvailable columns: ",
       paste(names(meta), collapse = ", "),
       "\nExpected one of: ", paste(sample_candidates, collapse = ", "))
}
if (is.na(cond_col)) {
  stop("Could not find a condition/group column in metadata.\nAvailable columns: ",
       paste(names(meta), collapse = ", "),
       "\nExpected one of: ", paste(cond_candidates, collapse = ", "))
}

meta <- meta %>%
  mutate(
    sample = as.character(.data[[sample_col]]),
    condition = as.factor(.data[[cond_col]])
  )

# Ensure sample names are unique
if (any(duplicated(meta$sample))) {
  dup <- unique(meta$sample[duplicated(meta$sample)])
  stop("Duplicate sample identifiers found in metadata: ", paste(dup, collapse = ", "))
}

message("Metadata columns used: sample='", sample_col, "', condition='", cond_col, "'")
message("N samples: ", nrow(meta))

# ----------------------------
# 3) Robust tx2gene handling
# ----------------------------
message("[2/6] Reading tx2gene: ", tx2gene_file)
tx2gene <- read_tsv(tx2gene_file, show_col_types = FALSE)

# Normalize column names for easier matching (case-insensitive)
names_lc <- tolower(names(tx2gene))

# Candidate names for tx and gene columns (case-insensitive)
tx_candidates   <- c("tx", "txname", "transcript", "transcript_id", "target_id", "name")
gene_candidates <- c("gene", "geneid", "gene_id", "ensgid", "gene_name", "geneid_raw")

find_lc <- function(names_lc, candidates_lc) {
  hit <- intersect(candidates_lc, names_lc)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

tx_lc   <- find_lc(names_lc, tx_candidates)
gene_lc <- find_lc(names_lc, gene_candidates)

if (is.na(tx_lc) || is.na(gene_lc)) {
  stop("tx2gene file must contain transcript and gene columns.\nAvailable columns: ",
       paste(names(tx2gene), collapse = ", "),
       "\nTranscript candidates: ", paste(tx_candidates, collapse = ", "),
       "\nGene candidates: ", paste(gene_candidates, collapse = ", "))
}

# Map back to original column names
tx_col   <- names(tx2gene)[which(names_lc == tx_lc)[1]]
gene_col <- names(tx2gene)[which(names_lc == gene_lc)[1]]

tx2gene <- tx2gene %>%
  transmute(
    tx = as.character(.data[[tx_col]]),
    gene = as.character(.data[[gene_col]])
  ) %>%
  filter(!is.na(tx), !is.na(gene))

# Optional: drop Ensembl version suffix in gene IDs (safe for gene_id tables)
tx2gene <- tx2gene %>%
  mutate(gene = sub("\\..*$", "", gene))

message("tx2gene columns used: tx='", tx_col, "', gene='", gene_col, "'")
message("tx2gene rows: ", nrow(tx2gene))

# ----------------------------
# 4) Build quant.sf file list
# ----------------------------
message("[3/6] Building quant.sf file list from: ", salmon_dir)
files <- file.path(salmon_dir, meta$sample, "quant.sf")
names(files) <- meta$sample

missing_q <- names(files)[!file.exists(files)]
if (length(missing_q) > 0) {
  stop("Missing quant.sf for these samples (check folder names under results/salmon/):\n",
       paste(missing_q, collapse = "\n"))
}

# ----------------------------
# 5) tximport + DESeq2 object
# ----------------------------
message("[4/6] Running tximport (Salmon)...")
txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreAfterBar = TRUE
)

message("[5/6] Creating DESeq2 object + VST...")
dds <- DESeqDataSetFromTximport(
  txi,
  colData = as.data.frame(meta),
  design = ~ condition
)

dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = FALSE)

# ----------------------------
# 6) Save objects
# ----------------------------
saveRDS(dds, out_dds)
saveRDS(vsd, out_vsd)

message("[6/6] Saved: ", out_dds)
message("[6/6] Saved: ", out_vsd)
message("Done.")
