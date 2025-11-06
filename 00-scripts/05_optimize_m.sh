#!/usr/bin/env bash
# ==============================================================
# Script: 05_optimize_m.sh
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Parameter optimization for Stacks 'm' (minimum stack depth).
#   Runs denovo_map.pl for m = 3–7 using a representative subset
#   of samples defined in subset_samples.txt.
#
#   Each run builds a separate de novo Stacks catalog while keeping
#   other parameters (M and n) at default values. Populations is run
#   automatically for each m to generate summary statistics (r0.8).
#
# Inputs:
#   - Clone-filtered reads: 03-analysis/03-clonefiltered/*.fq.gz
#   - Subset popmap:        01-info_files/subset_samples.txt
# Outputs:
#   - One run per m:        03-analysis/05-parameters/m/m<m>/
#   - Logs:                 05-logs/opt_m_m<m>.log
# ==============================================================

set -euo pipefail

# USER SETTINGS
READ_DIR="03-analysis/03-clonefiltered"       # location of input FASTQs
POPMAP="01-info_files/subset_samples.txt"     # subset popmap file

OUT_ROOT="03-analysis/05-parameters/m"        # base directory for outputs
LOG_DIR="05-logs"                             # where to store logs

CORES=8                                       # number of CPU threads
M_VALUES=(3 4 5 6 7)                          # range of m values to test
RUN_POPULATIONS=true                          # run populations for each m value

# Create necessary directories
mkdir -p "$OUT_ROOT" "$LOG_DIR"

# Basic sanity checks
[[ -d "$READ_DIR" ]] || { echo "[ERROR] READ_DIR not found: $READ_DIR"; exit 1; }
[[ -f "$POPMAP"   ]] || { echo "[ERROR] POPMAP not found: $POPMAP"; exit 1; }

# Loop over m values
for m in "${M_VALUES[@]}"; do
  # Create unique output and log paths for each test
  RUN_DIR="${OUT_ROOT}/m${m}"
  LOG="${LOG_DIR}/opt_m_m${m}.log"
  mkdir -p "$RUN_DIR"

  echo "[INFO] Starting denovo_map.pl with m=${m} (M and n = defaults)" | tee -a "$LOG"
  echo "Output: ${RUN_DIR}" | tee -a "$LOG"

  # Run Stacks de novo assembly
  denovo_map.pl \
    -T "$CORES" \
    --samples "$READ_DIR" \
    --popmap  "$POPMAP" \
    -o "$RUN_DIR" \
    -m "$m" \
    &>> "$LOG"

  # Run populations for each run
  if $RUN_POPULATIONS; then
    echo "[INFO] Running populations for m=${m}..." | tee -a "$LOG"
    populations \
      -P "$RUN_DIR" \
      -M "$POPMAP" \
      -t "$CORES" \
      -r 0.8 \
      --vcf \
      &>> "$LOG"
  fi

  echo "[OK] Finished m=${m}. Results saved in: $RUN_DIR" | tee -a "$LOG"
done

echo "[DONE] All m optimization runs completed."
echo "Check logs in: $LOG_DIR"
echo "Stacks outputs in: $OUT_ROOT"