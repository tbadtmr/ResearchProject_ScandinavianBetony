#!/usr/bin/env python3
# ==============================================================
# Script: 04_qc_metrics.py
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Extract GC% (FastQC general stats) and mean Phred score for
#   Read 1 from MultiQC outputs produced after clone filtering,
#   and write both metrics to TSV files. Also writes a merged
#   summary table for convenience.
#
# Inputs (from MultiQC):
#   - 04-results/multiqc_clonefiltered/multiqc_data/multiqc_general_stats.txt
#   - 04-results/multiqc_clonefiltered/multiqc_data/fastqc_per_sequence_quality_scores_plot.txt
#
# Outputs:
#   - 04-results/gc_content_read1.tsv          (Sample, GC_Content_Read1)
#   - 04-results/mean_phred_read1.tsv          (Sample, Mean_Phred_Read1)
#   - 04-results/qc_read1_summary.tsv          (merged metrics)
# ==============================================================

import os
import ast
import sys
import pandas as pd

DATA_DIR = "04-results/multiqc_clonefiltered/multiqc_data"
OUT_DIR = "04-results"

GEN_STATS = os.path.join(DATA_DIR, "multiqc_general_stats.txt")
PER_SEQ_QUAL = os.path.join(DATA_DIR, "fastqc_per_sequence_quality_scores_plot.txt")

GC_OUT = os.path.join(OUT_DIR, "gc_content_read1.tsv")
PHRED_OUT = os.path.join(OUT_DIR, "mean_phred_read1.tsv")
MERGED_OUT = os.path.join(OUT_DIR, "qc_read1_summary.tsv")

def ensure_exists(path):
    if not os.path.exists(path):
        sys.stderr.write(f"[ERROR] Missing file: {path}\n")
        sys.exit(1)

def clean_sample(name: str) -> str:
    # MultiQC sample names for clone-filtered R1 look like "<sample>.1.1"
    return name[:-4] if name.endswith(".1.1") else name

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # GC content (from general stats)
    ensure_exists(GEN_STATS)
    df_g = pd.read_csv(GEN_STATS, sep="\t")

    df_g_r1 = df_g[df_g["Sample"].str.endswith(".1.1")].copy()
    df_g_r1["Sample"] = df_g_r1["Sample"].apply(clean_sample)

    # FastQC GC% is usually under 'fastqc-percent_gc' in MultiQC
    if "fastqc-percent_gc" not in df_g_r1.columns:
        sys.stderr.write("[ERROR] Column 'fastqc-percent_gc' not found in general stats.\n")
        sys.exit(1)

    df_gc = df_g_r1[["Sample", "fastqc-percent_gc"]].rename(
        columns={"fastqc-percent_gc": "GC_Content_Read1"}
    )
    df_gc.to_csv(GC_OUT, sep="\t", index=False)

    # Mean Phred (from per-sequence quality plot)
    ensure_exists(PER_SEQ_QUAL)
    df_q = pd.read_csv(PER_SEQ_QUAL, sep="\t")
    df_q_r1 = df_q[df_q["Sample"].str.endswith(".1.1")].copy()

    mean_rows = []
    for _, row in df_q_r1.iterrows():
        sample_raw = row["Sample"]
        sample = clean_sample(sample_raw)

        # Per-sequence quality scores are stored as tuple-like strings "(x, y)"
        # across the remaining columns; ignore NaNs.
        points = []
        for val in row.iloc[1:]:
            if pd.isna(val):
                continue
            try:
                x, y = ast.literal_eval(val)
                points.append((float(x), float(y)))
            except Exception:
                continue

        total_bases = sum(y for _, y in points)
        weighted = sum(x * y for x, y in points)
        mean_phred = round(weighted / total_bases, 2) if total_bases > 0 else 0.0
        mean_rows.append((sample, mean_phred))

    df_phred = pd.DataFrame(mean_rows, columns=["Sample", "Mean_Phred_Read1"])
    df_phred.to_csv(PHRED_OUT, sep="\t", index=False)

    # Merged table
    df_merged = pd.merge(df_gc, df_phred, on="Sample", how="outer")
    df_merged.to_csv(MERGED_OUT, sep="\t", index=False)

    print(f"[OK] GC% (R1):            {GC_OUT}")
    print(f"[OK] Mean Phred (R1):     {PHRED_OUT}")
    print(f"[OK] Merged summary (R1): {MERGED_OUT}")

if __name__ == "__main__":
    main()
