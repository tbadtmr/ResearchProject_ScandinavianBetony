#!/usr/bin/env bash
# ==============================================================
# Script: 04_qc_clonefilter.sh
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Run FastQC on all clone-filtered FASTQ files and summarize
#   the results with MultiQC (report goes to 04-results).
#
# Inputs:
#   - Clone-filtered FASTQs: 03-analysis/03-clonefiltered/*.fq.gz
# Outputs:
#   - FastQC reports:        03-analysis/04-qc_clonefiltered/
#   - MultiQC summary:       04-results/multiqc_clonefiltered/multiqc_clonefiltered.html
#   - Log:                   05-logs/qc_clonefilter.log
# ==============================================================

set -euo pipefail

INPUT_DIR="03-analysis/03-clonefiltered"
FASTQC_DIR="03-analysis/04-qc_clonefiltered"
MULTIQC_DIR="04-results/multiqc_clonefiltered"
LOG_DIR="05-logs"
LOG_FILE="${LOG_DIR}/qc_clonefilter.log"

mkdir -p "$FASTQC_DIR" "$MULTIQC_DIR" "$LOG_DIR"

{
  echo "QC (clone-filtered) started at: $(date)"

  # Run FastQC on all clone-filtered FASTQs
  for fq in "${INPUT_DIR}"/*.fq.gz; do
    echo "Running FastQC on: ${fq}"
    fastqc "${fq}" -o "${FASTQC_DIR}"
  done

  # Summarize with MultiQC (report saved in results)
  echo "Running MultiQC..."
  multiqc "${FASTQC_DIR}" -o "${MULTIQC_DIR}" -n "multiqc_clonefiltered"

  echo "QC complete at: $(date)"
  echo "FastQC reports: ${FASTQC_DIR}"
  echo "MultiQC report: ${MULTIQC_DIR}/multiqc_clonefiltered.html"
} > "${LOG_FILE}" 2>&1

echo "QC for clone-filtered reads completed."
echo "Log:        ${LOG_FILE}"
echo "FastQC:     ${FASTQC_DIR}"
echo "MultiQC:    ${MULTIQC_DIR}/multiqc_clonefiltered.html"