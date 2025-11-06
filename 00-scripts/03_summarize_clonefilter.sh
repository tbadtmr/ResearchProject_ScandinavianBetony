#!/usr/bin/env bash
# ==============================================================
# Script: 03_summarize_clonefilter.sh
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Summarize read retention after clone filtering.
#   Counts reads before and after filtering and calculates
#   the percentage of reads retained per sample (exact filename).
#
# Inputs:
#   - Unfiltered reads: 03-analysis/02-demultiplexed/*.fq.gz
#   - Filtered reads:   03-analysis/03-clonefiltered/*.fq.gz
# Outputs:
#   - Summary table:    04-results/clone_filter_summary.tsv
# ==============================================================

set -euo pipefail

ORIG_DIR="03-analysis/02-demultiplexed"
FILTERED_DIR="03-analysis/03-clonefiltered"
OUTFILE="04-results/clone_filter_summary.tsv"

mkdir -p "$(dirname "$OUTFILE")"

echo -e "Sample\tOriginal_Reads\tFiltered_Reads\tPct_Retained" > "$OUTFILE"

have_any=false
# loop through all samples
for fq in "$ORIG_DIR"/*.fq.gz; do
  have_any=true
  [[ "$fq" == *.rem.* ]] && continue  # skip remainder files if any

  sample=$(basename "$fq" .fq.gz)  # sample name

  orig_lines=$(zcat -- "$fq" | wc -l)
  orig_reads=$((orig_lines / 4))

  filt_file="$FILTERED_DIR/${sample}.fq.gz"

  if [[ -f "$filt_file" ]]; then
    filt_lines=$(zcat -- "$filt_file" | wc -l)
    filt_reads=$((filt_lines / 4))

    if [[ "$orig_reads" -gt 0 ]]; then
      pct=$(awk "BEGIN { printf \"%.1f\", ($filt_reads / $orig_reads) * 100 }")
    else
      pct="0.0"
    fi

    echo -e "${sample}\t${orig_reads}\t${filt_reads}\t${pct}%" >> "$OUTFILE"
  else
    echo -e "${sample}\t${orig_reads}\tNA\tNA" >> "$OUTFILE"
  fi
done

# if no files were found, exit
if ! $have_any; then
  echo "No input files found in $ORIG_DIR (*.fq.gz). Exiting." >&2
  exit 1
fi

# final statement
echo "Clone filtering summary written to $OUTFILE"

