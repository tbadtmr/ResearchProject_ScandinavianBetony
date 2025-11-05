#!/usr/bin/env bash
# ==============================================================
# Script: 03_clonefilter.sh
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Run Stacks `clone_filter` on all demultiplexed paired FASTQ files.
#
# Inputs:
#   - Demultiplexed reads: 03-analysis/02-demultiplexed/*.1.fq.gz and *.2.fq.gz
# Outputs:
#   - Clone-filtered reads: 03-analysis/03-clonefiltered/
#   - Log:                  05-logs/clone_filter.log
# ==============================================================

set -euo pipefail

INPUT_DIR="03-analysis/02-demultiplexed"
OUTPUT_DIR="03-analysis/03-clonefiltered"
LOG_DIR="05-logs"
LOG_FILE="$LOG_DIR/clone_filter.log"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

echo "Clone filter log - $(date)" > "$LOG_FILE"

# Loop through all forward reads (*.1.fq.gz)
for fwd in "$INPUT_DIR"/*.1.fq.gz; do
  [[ "$fwd" == *.rem.1.fq.gz ]] && continue  # skip remainder files if any

  sample=$(basename "$fwd" .1.fq.gz)
  rev="$INPUT_DIR/${sample}.2.fq.gz"

  if [[ -f "$rev" ]]; then
    echo "[$(date)] Running clone_filter for $sample..." | tee -a "$LOG_FILE"

    clone_filter \
      -1 "$fwd" \
      -2 "$rev" \
      -i gzfastq \
      -o "$OUTPUT_DIR" \
      >> "$LOG_FILE" 2>&1

    echo "[$(date)] Finished $sample." | tee -a "$LOG_FILE"
  else
    echo "[$(date)] Reverse read not found for $sample, skipping..." | tee -a "$LOG_FILE"
  fi
done

echo "Clone filtering complete. Output: $OUTPUT_DIR" | tee -a "$LOG_FILE"
