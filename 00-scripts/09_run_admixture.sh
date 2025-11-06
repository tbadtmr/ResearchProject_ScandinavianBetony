#!/usr/bin/env bash
# ==============================================================
# Script: 09_run_admixture.sh
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Run ADMIXTURE (K = 1..KMAX, with REPS replicates each) on the
#   LD-pruned, no-outgroup PLINK files already created for PCA.
#   Collect per-replicate CV error and summarize by K.
#
# Usage:
#   ./00-scripts/09_run_admixture.sh
#
# Requirements: admixture, plink (only for earlier PCA step), awk, grep, sed
# ==============================================================

set -euo pipefail

# USER SETTINGS
RUN="r70_noTrentino_noOutgroup"           # LD-pruned, no outgroup
KMAX="${KMAX:-10}"                        # max K to test
REPS="${REPS:-50}"                        # replicates per K
THREADS="${THREADS:-4}"                   # CPU threads for admixture

# Reuse PLINK binaries from PCA step
PCA_DIR="03-analysis/08-downstream/01-pca"
PREFIX_SRC="${PCA_DIR}/${RUN}_pruned"     # .bed/.bim/.fam already exist

# Output location for ADMIXTURE runs
BASE="03-analysis/08-downstream/02-admixture/${RUN}"
mkdir -p "${BASE}/logs"

# Sanity: ensure PLINK files exist
for ext in bed bim fam; do
  [[ -f "${PREFIX_SRC}.${ext}" ]] || { echo "[ERROR] Missing ${PREFIX_SRC}.${ext}" >&2; exit 1; }
done

# Run from target folder so ADMIXTURE writes outputs there
cd "${BASE}"

# Copy (or symlink) PLINK trio locally for clarity (optional: use ln -s)
cp -n "../../01-pca/${RUN}_pruned."{bed,bim,fam} ./ 2>/dev/null || true
PREFIX="./${RUN}_pruned"

# ADMIXTURE loop
for K in $(seq 1 "${KMAX}"); do
  for R in $(seq 1 "${REPS}"); do
    SEED="${R}"
    echo ">>> K=${K} | replicate=${R} | seed=${SEED}"
    admixture --cv=10 -j"${THREADS}" -s "${SEED}" "${PREFIX}.bed" "${K}" \
      2>&1 | tee "logs/K${K}_rep${R}.log"

    # Rename outputs so replicates don't overwrite each other
    mv "${RUN}_pruned.${K}.Q" "${RUN}_pruned.${K}.rep${R}.Q"
    mv "${RUN}_pruned.${K}.P" "${RUN}_pruned.${K}.rep${R}.P"
  done
done

# Collect CV errors
{
  echo -e "K\treplicate\tcv_error"
  grep -H "CV error (K=" logs/K*.log \
    | sed -E 's#.*/K([0-9]+)_rep([0-9]+)\.log:.*CV error \(K=([0-9]+)\): *([0-9.]+).*#\1\t\2\t\4#'
} > cv_all.tsv

# Summarize CV by K
awk -F'\t' 'NR>1{k=$1; x=$3; n[k]++; s[k]+=x; q[k]+=x*x; if(!(k in mn)||x<mn[k]) mn[k]=x; if(!(k in mx)||x>mx[k]) mx[k]=x}
END{printf "K\tmean\tsd\tmin\tmax\tn\n";
    for(k in n){m=s[k]/n[k]; v=(q[k]/n[k])-(m*m); sd=(v>0?sqrt(v):0);
      printf "%d\t%.6f\t%.6f\t%.6f\t%.6f\t%d\n", k, m, sd, mn[k], mx[k], n[k];}}' \
  cv_all.tsv | sort -n -k1,1 > cv_summary.tsv

echo "[DONE] ADMIXTURE finished for ${RUN}"
echo "  Outputs → ${BASE}"
echo "    - ${RUN}_pruned.bed/.bim/.fam   (copied for convenience)"
echo "    - logs/*.log"
echo "    - cv_all.tsv, cv_summary.tsv"
echo "    - ${RUN}_pruned.<K>.rep<R>.Q / .P"