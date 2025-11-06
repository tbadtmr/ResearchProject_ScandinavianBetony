#!/usr/bin/env bash
# ==============================================================
# Script: 06_rerun_cleaned_samples.sh
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Re-run Stacks `populations` at r=0.7 using cleaned popmaps:
#     - popmap_noTrentino.txt
#     - popmap_noTrentino_noOutgroup.txt
#   Reuses the existing final denovo catalog built with <m> <M> <n>.
#
# Usage:
#   ./00-scripts/06_rerun_cleaned_samples.sh <m> <M> <n>
#   e.g. ./00-scripts/06_rerun_cleaned_samples.sh 6 4 4
#
# Inputs:
#   - Denovo catalog: 03-analysis/06-finalrun/denovo_m<m>_M<M>_n<n>/
#   - Popmaps:
#       01-info_files/popmap_noTrentino.txt
#       01-info_files/popmap_noTrentino_noOutgroup.txt
#
# Outputs:
#   - 03-analysis/06-finalrun/r70_noTrentino/
#   - 03-analysis/06-finalrun/r70_noTrentino_noOutgroup/
#   - Logs in 05-logs/
# ==============================================================

set -euo pipefail

# USER ARGUMENTS
if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <m_value> <M_value> <n_value>"
  echo "Example: $0 6 4 4"
  exit 1
fi
m_value="$1"
M_value="$2"
n_value="$3"

# USER SETTINGS
RESULTS_ROOT="03-analysis/06-finalrun"
LOG_DIR="05-logs"
CORES="${CORES:-8}"
R_VALUE="0.7"

DENOVO_DIR="${RESULTS_ROOT}/denovo_m${m_value}_M${M_value}_n${n_value}"

POPMAP_NOTRENTINO="01-info_files/popmap_noTrentino.txt"
POPMAP_NOTRENTINO_NOOG="01-info_files/popmap_noTrentino_noOutgroup.txt"

declare -A RUNS=(
  ["r70_noTrentino"]="$POPMAP_NOTRENTINO"
  ["r70_noTrentino_noOutgroup"]="$POPMAP_NOTRENTINO_NOOG"
)

# SANITY CHECKS
mkdir -p "$RESULTS_ROOT" "$LOG_DIR"

[[ -d "$DENOVO_DIR" ]] || { echo "[ERROR] Denovo directory not found: $DENOVO_DIR"; exit 1; }
[[ -f "$POPMAP_NOTRENTINO" ]] || { echo "[ERROR] Popmap not found: $POPMAP_NOTRENTINO"; exit 1; }
[[ -f "$POPMAP_NOTRENTINO_NOOG" ]] || { echo "[ERROR] Popmap not found: $POPMAP_NOTRENTINO_NOOG"; exit 1; }

command -v populations >/dev/null 2>&1 || { echo "[ERROR] populations not in PATH"; exit 1; }

# RUN POPULATIONS
for run_name in "${!RUNS[@]}"; do
  popmap="${RUNS[$run_name]}"
  outdir="${RESULTS_ROOT}/${run_name}"
  run_log="${LOG_DIR}/06_rerun_${run_name}.log"

  mkdir -p "$outdir"

  echo "[INFO] Running populations ($run_name) at r=${R_VALUE}"
  {
    echo "=== populations START $(date -u +"%Y-%m-%dT%H:%M:%SZ") ==="
    echo "INPUT:    $DENOVO_DIR"
    echo "POPMAP:   $popmap"
    echo "OUTPUT:   $outdir"
    echo "THREADS:  $CORES"
    echo "r:        $R_VALUE"
  } | tee -a "$run_log"

  populations \
    -P "$DENOVO_DIR" \
    -O "$outdir" \
    -M "$popmap" \
    -t "$CORES" \
    -r "$R_VALUE" \
    --vcf \
    --structure \
    --genepop \
    --plink \
    --phylip \
    &>> "$run_log"

  echo "[OK] populations finished: $run_name -> $outdir" | tee -a "$run_log"
done

echo "[DONE] Re-runs complete. Outputs in: $RESULTS_ROOT; logs in: $LOG_DIR"
