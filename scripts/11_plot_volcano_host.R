#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# -------- Parameters (edit only this block) --------
de_tsv      <- "results/host_gene/deseq2_gene_results_annotated.tsv"
lfc_col     <- "log2FoldChange"
padj_col    <- "padj"
label_col   <- "gene_name"     # if missing, will fallback to gene_id
lfc_thresh  <- 1
padj_thresh <- 0.05
label_top_n <- 12              # label top N most significant among Up genes
out_png     <- "results/figures/volcano_host.png"
width       <- 6
height      <- 5
dpi         <- 300
# -----------------------------------------------

stopifnot(file.exists(de_tsv))
dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)

res <- fread(de_tsv)

# Validate required cols
need <- c(lfc_col, padj_col)
miss <- setdiff(need, names(res))
if (length(miss) > 0) {
  stop("Missing required columns: ", paste(miss, collapse=", "),
       "\nAvailable: ", paste(names(res), collapse=", "))
}

# Choose label column safely
if (!(label_col %in% names(res))) {
  if ("gene_id" %in% names(res)) {
    label_col <- "gene_id"
  } else {
    label_col <- NULL
  }
}

df <- res %>%
  mutate(
    padj_safe = ifelse(is.na(.data[[padj_col]]), 1, .data[[padj_col]]),
    neglog10  = -log10(padj_safe),
    sig = case_when(
      .data[[padj_col]] < padj_thresh & .data[[lfc_col]] >=  lfc_thresh ~ "Up",
      .data[[padj_col]] < padj_thresh & .data[[lfc_col]] <= -lfc_thresh ~ "Down",
      TRUE ~ "NS"
    )
  )

p <- ggplot(df, aes(x = .data[[lfc_col]], y = neglog10, color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Volcano plot (DESeq2)",
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)",
    color = "Significance"
  )

# Optional labeling (top N Up genes by padj)
if (!is.null(label_col) && label_top_n > 0) {
  lab <- df %>%
    filter(sig == "Up") %>%
    arrange(.data[[padj_col]]) %>%
    head(label_top_n)

  p <- p + ggrepel::geom_text_repel(
    data = lab,
    aes(label = .data[[label_col]]),
    size = 3,
    max.overlaps = Inf
  )
}

ggsave(out_png, p, width = width, height = height, dpi = dpi)
message("Saved: ", out_png)
