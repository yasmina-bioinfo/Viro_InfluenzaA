#!/usr/bin/env bash
set -euo pipefail

# Edit these if needed
FASTQ_DIR="data/fastq"
SAMPLES="data/metadata/samples.csv"

mkdir -p "$FASTQ_DIR"

echo "[01] FASTQ dir: $FASTQ_DIR"
echo "[01] Samples: $SAMPLES"

# Expect a column called sra_run or run or accession; we’ll read sra_run (your file has it)
runs=$(awk -F',' 'NR>1 {print $4}' "$SAMPLES" | tr -d '\r' | grep -v '^$' | sort -u)

echo "[01] Runs found:"
echo "$runs"

# Example strategy:
# 1) Try SRA Toolkit (fasterq-dump)
# 2) If missing/corrupt, fallback to ENA (manual URL)

for r in $runs; do
  echo "[01] Fetching $r"

  # Skip if already present
  if ls "$FASTQ_DIR"/"${r}"*.fastq.gz >/dev/null 2>&1; then
    echo "  -> Already exists, skip."
    continue
  fi

  # SRA Toolkit
  if command -v fasterq-dump >/dev/null 2>&1; then
    echo "  -> Using fasterq-dump"
    fasterq-dump "$r" -O "$FASTQ_DIR" --threads 6
    gzip -f "$FASTQ_DIR"/"${r}"*.fastq
  else
    echo "  -> fasterq-dump not found. Install SRA Toolkit or use ENA method."
  fi
done

echo "[01] Done. If any run failed/corrupt, re-download via ENA for that specific run."
