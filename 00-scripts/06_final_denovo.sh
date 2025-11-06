#!/usr/bin/env bash
# ==============================================================
# Script: 06_final_denovo.sh
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Run the final de novo Stacks pipeline on the full dataset
#   with selected m/M/n parameters, then run four populations
#   filtering schemes (different -r and --min-maf).
#
# Usage:
#   bash 00-scripts/06_final_denovo.sh <m_value> <M_value> <n_value>
#   e.g. bash 00-scripts/06_final_denovo.sh 6 4 4
#
# Inputs:
#   - Clone-filtered reads: 03-analysis/03-clonefiltered/*.fq.gz
#   - Popmap (full set):    01-info_files/popmap.txt
#
# Outputs:
#   - denovo_map results:   03-results/06-finalrun/denovo_m<M>_M<M>_n<n>/
#   - populations runs:     03-results/06-finalrun/<run_name>/
#   - Logs:                 05-logs/06_final_denovo_*.log
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

# user SETTINGS
READ_DIR="03-analysis/03-clonefiltered"
POPMAP="01-info_files/popmap.txt"

RESULTS_ROOT="03-results/06-finalrun"
LOG_DIR="05-logs"

CORES=8

# populations runs:
declare -A RUNS=(
  ["r50"]="0.5 0.01"
  ["r60"]="0.6 0.01"
  ["r70"]="0.7 0.01"
  ["r80"]="0.8 0.01"
)

# PRE-RUN CHECKS
mkdir -p "$RESULTS_ROOT" "$LOG_DIR"

[[ -d "$READ_DIR" ]] || { echo "[ERROR] READ_DIR not found: $READ_DIR"; exit 1; }
[[ -f "$POPMAP"   ]] || { echo "[ERROR] POPMAP not found: $POPMAP"; exit 1; }

command -v denovo_map.pl >/dev/null 2>&1 || { echo "[ERROR] denovo_map.pl not in PATH"; exit 1; }
command -v populations   >/dev/null 2>&1 || { echo "[ERROR] populations not in PATH"; exit 1; }

# DENOVO MAP (FULL SET)
DENOVO_DIR="${RESULTS_ROOT}/denovo_m${m_value}_M${M_value}_n${n_value}"
DENOVO_LOG="${LOG_DIR}/06_final_denovo_denovo_m${m_value}_M${M_value}_n${n_value}.log"

mkdir -p "$DENOVO_DIR"

echo "[INFO] Starting denovo_map.pl (m=${m_value}, M=${M_value}, n=${n_value})"
echo "[INFO] Output: $DENOVO_DIR"
{
  echo "=== denovo_map.pl START $(date -u +"%Y-%m-%dT%H:%M:%SZ") ==="
  echo "CORES=${CORES}"
  echo "READ_DIR=${READ_DIR}"
  echo "POPMAP=${POPMAP}"
  echo "OUT=${DENOVO_DIR}"
} | tee -a "$DENOVO_LOG"

# Run denovo_map.pl
denovo_map.pl \
  -T "$CORES" \
  --samples "$READ_DIR" \
  --popmap  "$POPMAP" \
  -o "$DENOVO_DIR" \
  -m "$m_value" \
  -M "$M_value" \
  -n "$n_value" \
  &>> "$DENOVO_LOG"

echo "[OK] denovo_map.pl finished. Results in: $DENOVO_DIR" | tee -a "$DENOVO_LOG"

# POPULATIONS FILTERS
for run in "${!RUNS[@]}"; do
  IFS=' ' read -r r_val maf_val <<< "${RUNS[$run]}"
  OUTDIR="${RESULTS_ROOT}/${run}"
  RUN_LOG="${LOG_DIR}/06_final_denovo_${run}.log"

  mkdir -p "$OUTDIR"

  echo "[INFO] Running populations: ${run}  (r=${r_val}, MAF=${maf_val})"
  {
    echo "=== populations START $(date -u +"%Y-%m-%dT%H:%M:%SZ") ==="
    echo "INP=${DENOVO_DIR}"
    echo "OUT=${OUTDIR}"
    echo "r=${r_val}  min-maf=${maf_val}  CORES=${CORES}"
  } | tee -a "$RUN_LOG"

  populations \
    -P "$DENOVO_DIR" \
    -O "$OUTDIR" \
    -M "$POPMAP" \
    -t "$CORES" \
    -r "$r_val" \
    --min-maf "$maf_val" \
    --vcf \
    --structure \
    --genepop \
    --plink \
    --phylip \
    &>> "$RUN_LOG"

  echo "[OK] populations finished: ${run}. Output: $OUTDIR" | tee -a "$RUN_LOG"
done

echo "[DONE] Final de novo + 4 populations runs complete."
echo "Logs:   $LOG_DIR"
echo "Output: $RESULTS_ROOT"
