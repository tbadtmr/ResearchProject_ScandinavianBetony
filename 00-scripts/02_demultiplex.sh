#!/usr/bin/env bash
# ==============================================================
# Script: 02_demultiplex.sh
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Demultiplex RADseq reads for all groups (P2A–P2D) using Stacks
#   `process_radtags` with inline barcodes in R1 and initial QC/adapter trimming.
#
# Inputs:
#   - Raw reads: 02-raw_data/<group>/TC-*.fastq.gz     (groups: P2A,P2B,P2C,P2D)
#   - Barcodes : 01-info_files/barcodes/<group>.txt
# Outputs:
#   - Demultiplexed reads: 03-analysis/02-demultiplexed/
# ==============================================================

set -euo pipefail

RAW_DIR="02-raw_data"
BARCODE_DIR="01-info_files/barcodes"
OUT_DIR="03-analysis/02-demultiplexed"

mkdir -p "$OUT_DIR"

for BATCH in P2A P2B P2C P2D; do
  echo "Demultiplexing ${BATCH} ..."
  IN_DIR="${RAW_DIR}/${BATCH}"
  BARCODES="${BARCODE_DIR}/${BATCH}.txt"

  process_radtags -P \
    -p "$IN_DIR" \
    -o "$OUT_DIR" \
    -b "$BARCODES" \
    -e sbfI -i gzfastq -y gzfastq --inline_null \
    -r -c -q \
    --adapter_1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --adapter_2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --adapter_mm 2 \
    --force-poly-g-check \
    --threads 10

  echo "Finished ${BATCH}"
done

echo "All groups finished. Output: ${OUT_DIR}"
