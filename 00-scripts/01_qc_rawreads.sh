#!/usr/bin/env bash
# ==============================================================
# Script: 01_qc_rawreads.sh
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Runs FastQC on all raw FASTQ files (per sampling group)
#   and generates a combined MultiQC summary report.
#
# Input:  02-raw_data/P2*/TC-*.fastq.gz
# Output: 03-analysis/01-qc_rawreads/  (FastQC reports)
#         04-results/multiqc_rawreads/ (MultiQC summary)
# ==============================================================

# create output dirs
mkdir -p 03-analysis/01-qc_rawreads
mkdir -p 04-results/multiqc_rawreads

# run FastQC on all FASTQs inside group folders (P2A, P2B, ...)
for fq in 02-raw_data/P2*/TC-*.fastq.gz; do
    fastqc "$fq" -o 03-analysis/01-qc_rawreads/
done

# summarize with MultiQC
multiqc 03-analysis/01-qc_rawreads \
    -o 04-results/multiqc_rawreads \
    -n multiqc_rawreads

echo "QC done."
echo "FastQC reports: 03-analysis/01-qc_rawreads"
echo "MultiQC report: 04-results/multiqc_rawreads/multiqc_rawreads.html"
