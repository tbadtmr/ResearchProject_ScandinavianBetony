#!/usr/bin/env bash
# ==============================================================
# Script: 05_optimize_n.sh
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Optimize Stacks 'n' (distance allowed between stacks when
#   building the catalog across individuals). Runs two settings:
#     - n = M
#     - n = M + 1
#
#   For a given fixed m and M (passed via command line).
#
# Usage:
#   bash 00-scripts/05_optimize_n.sh <m_value> <M_value>
#   e.g. bash 00-scripts/05_optimize_n.sh 6 4
#
# Inputs:
#   - Clone-filtered reads: 03-analysis/03-clonefiltered/*.fq.gz
#   - Subset popmap:        01-info_files/subset_samples.txt
# Outputs:
#   - One run per n:        03-analysis/05-parameters/n/n<n>/
#   - Logs:                 05-logs/opt_n<n>.log
# ==============================================================

set -euo pipefail

# USER ARGUMENTS
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <m_value> <M_value>"
  echo "Example: $0 6 4"
  exit 1
fi
m_value="$1"
M_value="$2"

# USER SETTINGS
READ_DIR="03-analysis/03-clonefiltered"
POPMAP="01-info_files/subset_samples.txt"

OUT_ROOT="03-analysis/05-parameters/n"
LOG_DIR="05-logs"

CORES=8
RUN_POPULATIONS=true


mkdir -p "$OUT_ROOT" "$LOG_DIR"

# sanity checks
[[ -d "$READ_DIR" ]] || { echo "[ERROR] READ_DIR not found: $READ_DIR"; exit 1; }
[[ -f "$POPMAP"   ]] || { echo "[ERROR] POPMAP not found: $POPMAP"; exit 1; }

# The two n values to test
N_VALUES=( "$M_value" $((M_value + 1)) )

for n in "${N_VALUES[@]}"; do
  RUN_DIR="${OUT_ROOT}/n${n}"
  LOG="${LOG_DIR}/opt_n${n}.log"
  mkdir -p "$RUN_DIR"

  echo "[INFO] Running denovo_map.pl with m=${m_value}, M=${M_value}, n=${n}" | tee -a "$LOG"
  echo "Output: ${RUN_DIR}" | tee -a "$LOG"

  denovo_map.pl \
    -T "$CORES" \
    --samples "$READ_DIR" \
    --popmap  "$POPMAP" \
    -o "$RUN_DIR" \
    -m "$m_value" \
    -M "$M_value" \
    -n "$n" \
    &>> "$LOG"

  if $RUN_POPULATIONS; then
    echo "[INFO] Running populations for n=${n}..." | tee -a "$LOG"
    populations \
      -P "$RUN_DIR" \
      -M "$POPMAP" \
      -t "$CORES" \
      -r 0.8 \
      &>> "$LOG"
  fi

  echo "[OK] Finished n=${n}. Results: ${RUN_DIR}" | tee -a "$LOG"
done

echo "[DONE] n optimization complete (tested n=M and n=M+1)."
echo "Logs in: $LOG_DIR"
echo "Outputs in: $OUT_ROOT"

