#!/usr/bin/env bash
# ==============================================================
# Script: 12_fst_bootstrap.sh
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Bootstrap 95% CIs for pairwise FST using per-site .weir.fst
#   outputs produced by 12_diversity_stats.sh.
#
# Inputs:
#   - Per-case dirs:
#       03-analysis/08-downstream/04-diversity/stats_allGroups/
#       03-analysis/08-downstream/04-diversity/stats_scaniaSplit/
#
# Outputs:
#   - stats_<case>/FST_pairwise_with_CI.tsv
#   - FST_pairwise_with_CI.tsv (merged, deduplicated)
#   - Copies of case tables into 04-results/
#
# Usage:
#   B=1000 ./00-scripts/12_fst_bootstrap.sh
# ==============================================================

set -euo pipefail

B=${B:-1000}
BASE="03-analysis/08-downstream/04-diversity"
RDIR="04-results"
mkdir -p "${RDIR}"

CASE_DIRS=("${BASE}/stats_allGroups" "${BASE}/stats_scaniaSplit")

mean_from_log () { awk -F': ' '/Weir and Cockerham mean Fst estimate:/ {print $2}' "$1" | tail -n1; }

declare -a CASE_OUTS=()

for DIR in "${CASE_DIRS[@]}"; do
  if [[ ! -d "${DIR}" ]]; then
    echo "WARNING: ${DIR} not found, skipping."
    continue
  fi

  OUT="${DIR}/FST_pairwise_with_CI.tsv"
  echo -e "pop1\tpop2\tmean_FST\tCI95_low\tCI95_high\tn_sites_boot" > "${OUT}"

  shopt -s nullglob
  for F in "${DIR}"/FST_*_vs_*.weir.fst; do
    base=$(basename "$F" .weir.fst)
    popA=${base#FST_}; popA=${popA%%_vs_*}
    popB=${base##*_vs_}
    LOG="${F%.weir.fst}.log"

    # finite per-site FST values only
    awk 'NR>1 && $3==($3)+0 {print $3}' "$F" > /tmp/fst_vals.$$ || true
    n=$(wc -l < /tmp/fst_vals.$$ || echo 0)

    mean_log=$(mean_from_log "$LOG" || true)
    [ -z "${mean_log:-}" ] && mean_log="NA"

    if [ "$n" -lt 5 ]; then
      echo -e "$popA\t$popB\t$mean_log\tNA\tNA\t$n" >> "${OUT}"
      rm -f /tmp/fst_vals.$$ 2>/dev/null || true
      continue
    fi

    # bootstrap B means
    awk -v B="$B" 'BEGIN{srand()}
      { vals[NR]=$1 }
      END{
        n=NR;
        for(b=1; b<=B; b++){
          s=0;
          for(i=1;i<=n;i++){
            idx=int(rand()*n)+1; s+=vals[idx];
          }
          print s/n;
        }
      }' /tmp/fst_vals.$$ > /tmp/boot_means.$$

    sort -n /tmp/boot_means.$$ -o /tmp/boot_means.$$
    low_idx=$(awk -v B="$B" 'BEGIN{printf "%d\n",(int(0.025*B)+1)}')
    high_idx=$(awk -v B="$B" 'BEGIN{printf "%d\n",(int(0.975*B)+1)}')
    low=$(awk -v i="$low_idx" 'NR==i{printf "%.10f\n",$1}' /tmp/boot_means.$$)
    high=$(awk -v i="$high_idx" 'NR==i{printf "%.10f\n",$1}' /tmp/boot_means.$$)

    mean_out="$mean_log"
    [ "$mean_out" = "NA" ] && mean_out=$(awk '{s+=$1} END{printf "%.10f\n", s/NR}' /tmp/fst_vals.$$)

    echo -e "$popA\t$popB\t$mean_out\t$low\t$high\t$n" >> "${OUT}"
    rm -f /tmp/fst_vals.$$ /tmp/boot_means.$$ 2>/dev/null || true
  done
  shopt -u nullglob

  echo "Wrote ${OUT}"
  CASE_OUTS+=("${OUT}")

  # export a copy to 04-results with the case label in the name
  case_label=$(basename "${DIR}" | sed 's/^stats_//')
  cp -f "${OUT}" "${RDIR}/fst_pairwise_with_CI_${case_label}.tsv"
done

# merge
MERGED="${BASE}/FST_pairwise_with_CI.tsv"
if [ "${#CASE_OUTS[@]}" -gt 0 ]; then
  awk '
    BEGIN{FS=OFS="\t"}
    FNR==1 && NR!=1 { next }
    NR==1 { print; next }
    { key=$1 OFS $2; if(!(key in seen)){seen[key]=1; print} }
  ' "${CASE_OUTS[@]}" > "${MERGED}"
  echo "Wrote ${MERGED}"
  cp -f "${MERGED}" "${RDIR}/fst_pairwise_with_CI.tsv"
else
  echo "No per-case outputs produced; merged table not created."
fi
