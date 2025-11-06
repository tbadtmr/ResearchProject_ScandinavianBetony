#!/usr/bin/env bash
# ==============================================================
# Script: 05_optimize_M.sh
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Parameter optimization for Stacks 'M' (distance allowed between
#   stacks within individuals). Tests a range of M values while
#   keeping a fixed m value (provided as a command-line argument).
#
#   Each run builds a separate de novo Stacks catalog and runs
#   populations automatically for summary statistics.
#
# Usage:
#   bash 00-scripts/05_optimize_M.sh <m_value>
#   e.g. bash 00-scripts/05_optimize_M.sh 6
#
# Inputs:
#   - Clone-filtered reads: 03-analysis/03-clonefiltered/*.fq.gz
#   - Subset popmap:        01-info_files/subset_samples.txt
# Outputs:
#   - One run per M:        03-analysis/05-parameters/M/M<M>/
#   - Logs:                 05-logs/opt_M_M<M>.log
# ==============================================================

set -euo pipefail

# USER ARGUMENT
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <m_value>"
  echo "Example: $0 6"
  exit 1
fi
m_value="$1"   # fixed m value passed from the command line


# USER SETTINGS
READ_DIR="03-analysis/03-clonefiltered"       # location of input FASTQs
POPMAP="01-info_files/subset_samples.txt"     # subset popmap file

OUT_ROOT="03-analysis/05-parameters/M"        # base directory for outputs
LOG_DIR="05-logs"                             # where to store logs

CORES=8                                       # number of CPU threads
M_VALUES=(2 3 4 5 6 7 8)                      # range of M values to test
RUN_POPULATIONS=true                          # run populations for each M


# Create necessary directories
mkdir -p "$OUT_ROOT" "$LOG_DIR"

# Basic sanity checks
[[ -d "$READ_DIR" ]] || { echo "[ERROR] READ_DIR not found: $READ_DIR"; exit 1; }
[[ -f "$POPMAP"   ]] || { echo "[ERROR] POPMAP not found: $POPMAP"; exit 1; }

# Loop over M values
for M in "${M_VALUES[@]}"; do
  RUN_DIR="${OUT_ROOT}/M${M}"
  LOG="${LOG_DIR}/opt_M_M${M}.log"
  mkdir -p "$RUN_DIR"

  echo "[INFO] Starting denovo_map.pl with m=${m_value}, M=${M} (n = default)" | tee -a "$LOG"
  echo "Output: ${RUN_DIR}" | tee -a "$LOG"

  # Run Stacks assembly for this M
  denovo_map.pl \
    -T "$CORES" \
    --samples "$READ_DIR" \
    --popmap  "$POPMAP" \
    -o "$RUN_DIR" \
    -m "$m_value" \
    -M "$M" \
    &>> "$LOG"

  # Run populations to generate summary stats
  if $RUN_POPULATIONS; then
    echo "[INFO] Running populations for M=${M}..." | tee -a "$LOG"
    populations \
      -P "$RUN_DIR" \
      -M "$POPMAP" \
      -t "$CORES" \
      -r 0.8 \
      --vcf \
      &>> "$LOG"
  fi

  echo "[OK] Finished M=${M}. Results saved in: $RUN_DIR" | tee -a "$LOG"
done

echo "[DONE] All M optimization runs completed."
echo "Check logs in: $LOG_DIR"
echo "Stacks outputs in: $OUT_ROOT"
