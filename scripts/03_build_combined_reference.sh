#!/usr/bin/env bash
set -euo pipefail

GTF_GZ="ref/gencode/gencode.v49.annotation.gtf.gz"
TX_FA_GZ="ref/gencode/gencode.v49.transcripts.fa.gz"
VIRAL_FA="ref/iav_wsn33/wsn33_segments.fa"
OUT_DIR="ref/combined"
OUT_FA="$OUT_DIR/combined_transcripts.fa"

mkdir -p "$OUT_DIR"

echo "[03] Decompress host transcripts..."
zcat "$TX_FA_GZ" > "$OUT_DIR/host_transcripts.fa"

echo "[03] Combine host + viral..."
cat "$OUT_DIR/host_transcripts.fa" "$VIRAL_FA" > "$OUT_FA"

gzip -f "$OUT_FA" || true

echo "[03] Combined reference created:"
echo "  - $OUT_FA.gz"
