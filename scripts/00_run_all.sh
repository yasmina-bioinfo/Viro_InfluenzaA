#!/usr/bin/env bash
set -euo pipefail

echo "[00] Running full pipeline from project root: $(pwd)"

bash scripts/02_qc_fastqc_multiqc.sh
bash scripts/03_build_combined_reference.sh
bash scripts/04_salmon_index.sh
bash scripts/05_salmon_quant.sh
Rscript scripts/06_compute_viral_load.R
Rscript scripts/07_host_deseq2.R
Rscript scripts/08_build_dds_gene.R
Rscript scripts/09_plot_pca.R
Rscript scripts/10_plot_volcano.R
Rscript scripts/11_plot_heatmap_isgs.R

echo "[00] Done. Outputs in results/."
