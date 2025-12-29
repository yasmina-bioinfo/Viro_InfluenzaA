#!/usr/bin/env bash
set -euo pipefail

FASTQ_DIR="data/fastq"
QC_DIR="results/qc"
THREADS=8

mkdir -p "$QC_DIR/fastqc"

echo "[02] Running FastQC..."
fastqc -t "$THREADS" -o "$QC_DIR/fastqc" "$FASTQ_DIR"/*.fastq.gz

echo "[02] Running MultiQC..."
multiqc -o "$QC_DIR" "$QC_DIR/fastqc"

echo "[02] Saved: $QC_DIR/multiqc_report.html"
