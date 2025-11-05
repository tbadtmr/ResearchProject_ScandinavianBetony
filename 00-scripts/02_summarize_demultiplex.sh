#!/usr/bin/env bash
# ==============================================================
# Script: 02_summarize_demultiplex.sh
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Parse Stacks process_radtags log files created during demultiplexing,
#   and build a TSV table summarizing per-sample read retention metrics.
#
# Inputs:
#   - Logs: 03-analysis/02-demultiplexed/process_radtags.P2*.log
# Outputs:
#   - Summary TSV: 04-results/demultiplexing_summary.tsv
# ==============================================================

set -euo pipefail

OUTDIR="04-results"
LOGDIR="03-analysis/02-demultiplexed"
mkdir -p "$OUTDIR"
OUTFILE="${OUTDIR}/demultiplexing_summary.tsv"

echo -e "Batch\tSample\tTotal_Reads\tRetained_Reads\tPct_Retained\tProperly_Paired\tPct_Properly_Paired" > "$OUTFILE"

for LOG in ${LOGDIR}/process_radtags.P2*.log; do
  [ -e "$LOG" ] || { echo "No log files found in ${LOGDIR}"; exit 1; }
  BATCH=$(basename "$LOG" | sed 's/.*process_radtags.\(P2[ABCD]\).log/\1/')

  awk -v batch="$BATCH" '
    BEGIN { in_section=0 }
    /^BEGIN per_barcode_raw_read_counts/ { in_section=1; next }
    /^END per_barcode_raw_read_counts/   { in_section=0 }
    in_section == 1 {
      if (NF >= 11) {
        sample=$2
        total=$3
        retained=$7
        pct_retained=$8
        properly_paired=$9
        pct_properly_paired=$10
        gsub("%","",pct_retained)
        gsub("%","",pct_properly_paired)
        print batch "\t" sample "\t" total "\t" retained "\t" pct_retained "\t" properly_paired "\t" pct_properly_paired
      }
    }
  ' "$LOG" >> "$OUTFILE"
done

echo "Demultiplexing summary saved to: $OUTFILE"
