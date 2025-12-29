#!/usr/bin/env bash
set -euo pipefail

FASTQ_DIR="data/fastq"
SAMPLES="data/metadata/samples.csv"
INDEX_DIR="ref/combined/salmon_index"
OUT_DIR="results/salmon"
THREADS=8
LIBTYPE="A"
MODE="paired"   # "paired" or "single"

mkdir -p "$OUT_DIR"

echo "[05] Quant mode: $MODE"
echo "[05] Reading samples from: $SAMPLES"

while IFS=',' read -r sample_id condition replicate sra_run; do
  # Skip header
  if [[ "$sample_id" == "sample_id" ]]; then
    continue
  fi

  # Clean CRLF (Windows line endings)
  sample_id=$(echo "$sample_id" | tr -d '\r')
  sra_run=$(echo "$sra_run" | tr -d '\r')

  echo "[05] Sample: $sample_id | Run: $sra_run"

  sample_out="$OUT_DIR/$sample_id"
  mkdir -p "$sample_out"
if [ -f "$sample_out/quant.sf" ]; then
  echo "[05] quant.sf exists for $sample_id -> skip"
  continue
fi

if [ "$MODE" = "paired" ]; then
    R1=$(ls "$FASTQ_DIR"/"$sra_run"*"_1".fastq.gz \
        "$FASTQ_DIR"/"$sra_run"*"_R1".fastq.gz 2>/dev/null | head -n 1 || true)

    R2=$(ls "$FASTQ_DIR"/"$sra_run"*"_2".fastq.gz \
        "$FASTQ_DIR"/"$sra_run"*"_R2".fastq.gz 2>/dev/null | head -n 1 || true)


if [ -z "$R1" ] || [ -z "$R2" ]; then
      echo "[ERROR] Missing paired FASTQs for $sra_run"
      echo "        Expected something like: ${sra_run}*_1.fastq.gz and ${sra_run}*_2.fastq.gz in $FASTQ_DIR"
      exit 1
fi

    salmon quant \
      -i "$INDEX_DIR" \
      -l "$LIBTYPE" \
      -1 "$R1" \
      -2 "$R2" \
      -p "$THREADS" \
      --validateMappings \
      -o "$sample_out"

  else
    READ=$(ls "$FASTQ_DIR"/"$sra_run"*".fastq.gz" \
           "$FASTQ_DIR"/"$sra_run"*".fq.gz" 2>/dev/null | head -n 1 || true)


    if [ -z "$READ" ]; then
      echo "[ERROR] Missing single-end FASTQ for $sra_run"
      echo "        Expected something like: ${sra_run}*.fastq.gz in $FASTQ_DIR"
      exit 1
    fi

    salmon quant \
      -i "$INDEX_DIR" \
      -l "$LIBTYPE" \
      -r "$READ" \
      -p "$THREADS" \
      --validateMappings \
      -o "$sample_out"
  fi

done < "$SAMPLES"

echo "[05] Done: quant.sf written to $OUT_DIR/<sample_id>/quant.sf"
