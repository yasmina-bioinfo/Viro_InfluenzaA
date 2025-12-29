suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(dplyr)
  library(pheatmap)
})

# -------- Parameters --------
vsd_rds   <- "results/host_gene/vsd_gene.rds"
de_annot  <- "results/host_gene/deseq2_gene_results_annotated.tsv"
group_col <- "condition"

# ISGs as GENE NAMES (human-readable)
isg_names <- c("IFIT1","IFIT2","IFIT3","ISG15","MX1","MX2","OAS1","OAS2","OAS3","RSAD2","IFI27","STAT1","IRF7")

out_png  <- "results/figures/heatmap_ISGs.png"
width_px <- 1100
height_px <- 900
# ---------------------------

stopifnot(file.exists(vsd_rds))
stopifnot(file.exists(de_annot))
dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)

vsd <- readRDS(vsd_rds)
mat <- assay(vsd)

# Column annotation
cd <- as.data.frame(colData(vsd))
if (!(group_col %in% colnames(cd))) {
  stop("Grouping column not found in colData: ", group_col,
       "\nAvailable: ", paste(colnames(cd), collapse = ", "))
}
ann_col <- cd[, group_col, drop=FALSE]

# Detect whether rownames are gene_id (ENSG...) or gene_name
rn <- rownames(mat)
is_ensg <- any(grepl("^ENSG", rn))

# Load mapping from annotated DE results
map <- fread(de_annot) %>%
  select(gene_id, gene_name) %>%
  distinct() %>%
  filter(!is.na(gene_id), !is.na(gene_name))

if (is_ensg) {
  # rownames are ENSG -> map ISG gene_name to gene_id
  isg_gene_ids <- map %>%
    filter(gene_name %in% isg_names) %>%
    pull(gene_id) %>%
    unique()

  keep <- intersect(isg_gene_ids, rn)
  if (length(keep) == 0) {
    stop("No ISGs found after mapping gene_name -> gene_id.\n",
         "Check that gene_name column matches your ISG list.")
  }

  submat <- mat[keep, , drop=FALSE]

  # Rename rows to gene_name for readability in the heatmap
  id2name <- map %>% distinct(gene_id, gene_name)
  rownames(submat) <- id2name$gene_name[match(rownames(submat), id2name$gene_id)]

} else {
  # rownames are already gene_name
  keep <- intersect(isg_names, rn)
  if (length(keep) == 0) stop("None of the ISGs were found in rownames(vsd assay).")
  submat <- mat[keep, , drop=FALSE]
}

# Save heatmap
png(out_png, width = width_px, height = height_px)
pheatmap(
  submat,
  scale = "row",
  annotation_col = ann_col,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 9,
  main = "Heatmap: ISGs (VST, scaled by gene)"
)
dev.off()

message("Saved: ", out_png)
