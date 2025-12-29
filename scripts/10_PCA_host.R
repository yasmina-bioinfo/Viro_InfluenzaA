#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

# -------- Parameters (edit only this block) --------
vsd_rds   <- "results/host_gene/vsd_gene.rds"
group_col <- "condition"          # column in colData(dds)
out_png   <- "results/figures/PCA_host.png"
width     <- 6
height    <- 5
dpi       <- 300
# -----------------------------------------------

stopifnot(file.exists(vsd_rds))
dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)

vsd <- readRDS(vsd_rds)

# Validate grouping column
cd <- as.data.frame(colData(vsd))
if (!(group_col %in% colnames(cd))) {
  stop("Grouping column not found in colData: ", group_col,
       "\nAvailable: ", paste(colnames(cd), collapse = ", "))
}

pca_df <- plotPCA(vsd, intgroup = group_col, returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))

p <- ggplot(pca_df, aes(PC1, PC2, color = .data[[group_col]])) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  labs(color = group_col)

ggsave(out_png, p, width = width, height = height, dpi = dpi)
message("Saved: ", out_png)
