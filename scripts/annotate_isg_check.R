#!/usr/bin/env Rscript

# ============================================================
# Annotate DESeq2 gene-level results (ENSG -> gene_name)
# and check canonical ISGs (Interferon-Stimulated Genes)
#
# Inputs:
#   - results/host_gene/deseq2_gene_results.tsv   (must contain gene_id)
#   - ref/gencode.gene.v49.annotation.gtf         (GENCODE v49 GTF)
#
# Outputs:
#   - results/host_gene/deseq2_gene_results_annotated.tsv
#   - results/host_gene/isg_hits.tsv
#   - results/host_gene/top_upregulated.tsv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ----------------------------
# 0) File paths (edit if needed)
# ----------------------------
deseq_file <- "results/host_gene/deseq2_gene_results.tsv"
gtf_file   <- "ref/gencode/gencode.v49.annotation.gtf"

out_annot  <- "results/host_gene/deseq2_gene_results_annotated.tsv"
out_isg    <- "results/host_gene/isg_hits.tsv"
out_topup  <- "results/host_gene/top_upregulated.tsv"

# ----------------------------
# 1) Sanity checks: do files exist?
# ----------------------------
if (!file.exists(deseq_file)) stop("DESeq2 results file not found: ", deseq_file)
if (!file.exists(gtf_file))   stop("GTF file not found: ", gtf_file)

message("[1/5] Reading GTF: ", gtf_file)

# ----------------------------
# 2) Read GTF and extract gene_id / gene_name
#    We keep only 'gene' features to avoid duplicates.
# ----------------------------
gtf <- fread(
  gtf_file,
  sep = "\t",
  header = FALSE,
  data.table = TRUE,
  quote = "",
  fill = TRUE,
  comment.char = "#"
)

# Standard GTF columns (9)
setnames(gtf, c("seqname","source","feature","start","end","score","strand","frame","attribute"))

gtf_gene <- gtf[feature == "gene", .(attribute)]

# Helper to extract an attribute like: gene_id "ENSG..."; gene_name "STAT1";
get_attr <- function(x, key) {
  m <- regmatches(x, regexpr(paste0(key, " \"[^\"]+\""), x))
  ifelse(m == "", NA_character_, gsub(paste0("^", key, " \"|\"$"), "", m))
}

message("[2/5] Parsing gene_id and gene_name from GTF attributes...")

annot <- gtf_gene %>%
  transmute(
    gene_id_raw = get_attr(attribute, "gene_id"),
    gene_name   = get_attr(attribute, "gene_name")
  ) %>%
  filter(!is.na(gene_id_raw)) %>%
  mutate(
    # Remove Ensembl version suffix if present (ENSG... .12 -> ENSG...)
    gene_id = sub("\\..*$", "", gene_id_raw)
  ) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  select(gene_id, gene_name)

message("Annotation rows (genes): ", nrow(annot))

# ----------------------------
# 3) Read DESeq2 results
#    IMPORTANT: your file already contains 'gene_id' column.
# ----------------------------
message("[3/5] Reading DESeq2 results: ", deseq_file)
res <- fread(deseq_file)

if (!("gene_id" %in% names(res))) {
  stop(
    "Expected column 'gene_id' not found in DESeq2 results.\n",
    "Available columns: ", paste(names(res), collapse = ", ")
  )
}

# Ensure gene_id is clean (remove any potential version suffix)
res <- res %>%
  mutate(gene_id = sub("\\..*$", "", gene_id))

# Check required DESeq2 columns exist
need_cols <- c("log2FoldChange", "padj")
missing <- setdiff(need_cols, names(res))
if (length(missing) > 0) {
  stop(
    "Missing required DESeq2 columns: ", paste(missing, collapse = ", "),
    "\nAvailable columns: ", paste(names(res), collapse = ", ")
  )
}

# ----------------------------
# 4) Merge annotation (ENSG -> gene_name)
# ----------------------------
message("[4/5] Merging annotation into DESeq2 results...")

res_annot <- res %>%
  left_join(annot, by = "gene_id") %>%
  relocate(gene_id, gene_name, .before = 1)

fwrite(res_annot, out_annot, sep = "\t")
message("Wrote annotated results: ", out_annot)

# ----------------------------
# 5) ISG check + export top upregulated genes
# ----------------------------
# Canonical ISGs (extend anytime)
isg_list <- c(
  "IFIT1","IFIT2","IFIT3","ISG15","MX1","MX2","OAS1","OAS2","OAS3",
  "IFI6","IFI27","RSAD2","DDX58","IFIH1","STAT1","STAT2","IRF7","IRF9",
  "HERC5","BST2","USP18","CXCL10","CXCL11"
)

# Top upregulated genes (by log2FC, padj not NA)
top_up <- res_annot %>%
  filter(!is.na(padj)) %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange))

# Save a top list (e.g., 200 genes)
fwrite(head(top_up, 200), out_topup, sep = "\t")
message("Wrote top upregulated genes: ", out_topup)

# ISG hits table
isg_hits <- res_annot %>%
  filter(!is.na(gene_name)) %>%
  filter(gene_name %in% isg_list) %>%
  arrange(padj)

fwrite(isg_hits, out_isg, sep = "\t")
message("Wrote ISG hits: ", out_isg)

# Console summary
message("\n===== ISG SUMMARY =====")
message("ISGs found in results: ", nrow(isg_hits), " / ", length(isg_list))

if (nrow(isg_hits) > 0) {
  message("Top ISG hits (up to 10):")
  print(
    isg_hits %>%
      select(gene_id, gene_name, log2FoldChange, padj) %>%
      head(10)
  )
} else {
  message("No ISGs detected via gene_name. Double-check annotation and DESeq2 results.")
}

message("\nDone.")
