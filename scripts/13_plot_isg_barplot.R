# 13_plot_isg_barplot.R
# Purpose: Portfolio-ready ISG barplot (base R, up to 10 genes)

# Load ISG hits table
isg <- read.table(
  "results/host_gene/isg_hits.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Safety checks
required_cols <- c("gene_name", "log2FoldChange", "padj")
missing_cols <- setdiff(required_cols, colnames(isg))
if (length(missing_cols) > 0) {
  stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
}

# Canonical ISG candidates (priority order)
candidate_genes <- c(
  "CXCL10",
  "IFIT2",
  "OAS3",
  "IFIT1",
  "ISG15",
  "IFIT3",
  "IFI27",
  "USP18",
  "IRF7",
  "IFITM3"
)

# Keep only significant UP genes from the candidate list
isg_sel <- isg[
  isg$gene_name %in% candidate_genes &
  isg$log2FoldChange > 1 &   # explicit upregulation threshold
  isg$padj < 0.05,
]

# Order genes by predefined biological order
isg_sel$gene_name <- factor(isg_sel$gene_name, levels = candidate_genes)
isg_sel <- isg_sel[order(isg_sel$gene_name), ]

# Limit to max 10 genes (in case more are present)
isg_sel <- head(isg_sel, 10)

# Stop if nothing to plot (safety)
if (nrow(isg_sel) == 0) {
  stop("No upregulated ISGs found for plotting.")
}

# Plot
png(
  filename = "results/figures/isg_barplot_bulk_influenza.png",
  width = 2200,
  height = 1400,
  res = 300
)

par(mar = c(5, 9, 3, 1))  # larger left margin for gene labels

barplot(
  isg_sel$log2FoldChange,
  names.arg = as.character(isg_sel$gene_name),
  horiz = TRUE,
  las = 1,
  col = "#4C72B0",
  xlab = "log2 fold change (infected vs control)",
  main = "Antiviral interferon-stimulated genes",
  sub = "Canonical ISGs selected with log2FC > 1 and adjusted p-value < 0.05",
  cex.names = 1.2,
  cex.lab = 1.1,
  cex.main = 1.2
)


abline(v = 0, lwd = 2)
dev.off()
