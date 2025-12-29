#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
  library(tidyr)
})

# ---------------- Parameters ----------------
samples_csv  <- "data/metadata/samples.csv"
salmon_dir   <- "results/salmon"
tx2gene_file <- "ref/gencode/tx2gene_gencode_v49.tsv"

out_dir <- "results/host_matrix"
out_numreads <- file.path(out_dir, "host_NumReads_matrix.tsv")
out_tpm      <- file.path(out_dir, "host_TPM_matrix.tsv")
# -------------------------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(samples_csv), dir.exists(salmon_dir), file.exists(tx2gene_file))

# Read metadata (your columns)
meta <- read_csv(samples_csv, show_col_types = FALSE) %>%
  mutate(sample = as.character(sample_id))

# Read tx2gene and standardize columns to tx / gene
tx2gene <- read_tsv(tx2gene_file, show_col_types = FALSE)

if (all(c("TXNAME","GENEID") %in% names(tx2gene))) {
  tx2gene <- tx2gene %>% rename(tx = TXNAME, gene = GENEID)
} else if (all(c("tx","gene") %in% names(tx2gene))) {
  # ok
} else if (all(c("transcript_id","gene_id") %in% names(tx2gene))) {
  tx2gene <- tx2gene %>% rename(tx = transcript_id, gene = gene_id)
} else {
  stop("tx2gene columns not recognized. Available: ", paste(names(tx2gene), collapse=", "))
}

tx2gene <- tx2gene %>%
  transmute(tx = as.character(tx),
            gene = as.character(gene)) %>%
  mutate(gene = sub("\\..*$", "", gene)) %>%
  distinct()

# Function to read one quant.sf and keep host-only transcripts (present in tx2gene$tx)
read_one <- function(sample_name) {
  qf <- file.path(salmon_dir, sample_name, "quant.sf")
  stopifnot(file.exists(qf))
  q <- fread(qf) %>%
    select(Name, TPM, NumReads)
  q_host <- q %>% filter(Name %in% tx2gene$tx)
  q_host$sample <- sample_name
  q_host
}

message("[1/3] Reading quant.sf for ", nrow(meta), " samples...")
allq <- bind_rows(lapply(meta$sample, read_one))

message("[2/3] Building matrices...")

numreads_mat <- allq %>%
  select(Name, sample, NumReads) %>%
  pivot_wider(names_from = sample, values_from = NumReads, values_fill = 0) %>%
  arrange(Name)

tpm_mat <- allq %>%
  select(Name, sample, TPM) %>%
  pivot_wider(names_from = sample, values_from = TPM, values_fill = 0) %>%
  arrange(Name)

message("[3/3] Saving...")
fwrite(numreads_mat, out_numreads, sep = "\t")
fwrite(tpm_mat, out_tpm, sep = "\t")

message("Saved: ", out_numreads)
message("Saved: ", out_tpm)
message("Done.")
