#!/usr/bin/env bash
set -euo pipefail

COMBINED_FA_GZ="ref/combined/combined_transcripts.fa.gz"
INDEX_DIR="ref/combined/salmon_index"
THREADS=8

mkdir -p "$INDEX_DIR"

echo "[04] Building Salmon index..."
salmon index -t "$COMBINED_FA_GZ" -i "$INDEX_DIR" -p "$THREADS"

echo "[04] Done: $INDEX_DIR"
