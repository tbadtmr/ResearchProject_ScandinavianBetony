#!/usr/bin/env bash
# ==============================================================
# Script: 06_pop_summary.sh
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Summarize key metrics from multiple Stacks `populations` runs.
#   Thresholds:
#     - Low-depth samples: mean DP < 20
#     - High-missing samples: >50% missing genotypes
#
# Usage:
#   ./00-scripts/06_pop_summary.sh [OUT_TSV]
# ==============================================================

set -euo pipefail

RESULTS_ROOT="03-analysis/06-finalrun"
RUNS=( r50 r60 r70 r80 )

OUT_TSV="${1:-}"
LOW_DP_THRESHOLD=20
MISSING_FRAC_THRESHOLD=0.5

# helper functions

# vcf_path: return the path to the VCF file (compressed or plain)
# for a given populations output directory.
vcf_path() {
  local dir="$1"
  if [[ -f "${dir}/populations.snps.vcf.gz" ]]; then
    echo "${dir}/populations.snps.vcf.gz"
  elif [[ -f "${dir}/populations.snps.vcf" ]]; then
    echo "${dir}/populations.snps.vcf"
  else
    echo ""
  fi
}

# count_variants: fast variant count using bcftools index if available.
# Falls back to counting non-header lines if unindexed.
count_variants() {
  local vcf="$1"
  if [[ "$vcf" =~ \.vcf\.gz$ ]]; then
    if bcftools index -n "$vcf" >/dev/null 2>&1; then
      bcftools index -n "$vcf"
      return
    fi
  fi
  bcftools view -H "$vcf" | wc -l
}

# checks
command -v bcftools >/dev/null 2>&1 || { echo "[ERROR] bcftools not found in PATH" >&2; exit 1; }
[[ -d "$RESULTS_ROOT" ]] || { echo "[ERROR] Missing directory: $RESULTS_ROOT" >&2; exit 1; }

# table header
header=$'run\tr\tmaf\tloci\tvariant_sites\tmean_depth\tmean_site_missingness\tsamples_meanDP_lt20\tsamples_gt50pct_missing'
echo "$header"
[[ -n "$OUT_TSV" ]] && { mkdir -p "$(dirname "$OUT_TSV")"; echo "$header" > "$OUT_TSV"; }

# loop over each run
for run in "${RUNS[@]}"; do
  DIR="${RESULTS_ROOT}/${run}"
  VCF="$(vcf_path "$DIR")"

  if [[ -z "$VCF" ]]; then
    echo "[WARN] No VCF found in ${DIR} (skipping)" >&2
    continue
  fi

  # r value: r50 -> 0.50
  r_digits="$(sed -E 's/^r([0-9]+).*/\1/' <<< "$run")"
  r_val="$(awk -v r="$r_digits" 'BEGIN{printf "%.2f", r/100}')"
  maf_val="NA"

  variants="$(count_variants "$VCF")"
  loci="$(bcftools query -f '%CHROM\n' "$VCF" | sort -u | wc -l)"

  mean_dp="$(
    bcftools query -f '[%DP\t]\n' "$VCF" \
    | tr '\t' '\n' \
    | awk '$1!="."{sum+=$1; n++} END{if(n>0) printf "%.2f", sum/n; else print "NA"}'
  )"

  mean_site_miss="$(
    bcftools query -f '[%GT\t]\n' "$VCF" \
    | awk '{
        tot=NF; miss=0;
        for(i=1;i<=NF;i++) if($i=="./." || $i==".") miss++;
        sum+=miss/tot; n++;
      } END{
        if(n>0) printf "%.4f", sum/n; else print "NA"
      }'
  )"

  # per-sample: mean DP and fraction missing; thresholds: mean DP < 20; missing > 0.5
  read -r lowdp gt50miss < <(
    paste \
      <(bcftools query -f '[%DP\t]\n' "$VCF") \
      <(bcftools query -f '[%GT\t]\n' "$VCF") \
    | awk -F'\t' -v TH_DP="$LOW_DP_THRESHOLD" -v TH_MISS="$MISSING_FRAC_THRESHOLD" '
      NR==1{
        nfield=NF/2;
        for(i=1;i<=nfield;i++){sumdp[i]=0; ndp[i]=0; miss[i]=0; tot[i]=0;}
      }
      {
        nfield=NF/2;
        for(i=1;i<=nfield;i++){
          dp=$i; gt=$(i+nfield);
          if(dp!="." && dp!=""){sumdp[i]+=dp; ndp[i]++}
          tot[i]++; if(gt=="./." || gt==".") miss[i]++;
        }
      }
      END{
        low=0; highmiss=0;
        for(i=1;i<=nfield;i++){
          if(ndp[i]>0 && (sumdp[i]/ndp[i])<TH_DP) low++;
          if(tot[i]>0 && (miss[i]/tot[i])>TH_MISS) highmiss++;
        }
        printf "%d %d", low, highmiss;
      }'
  )

  line="${run}\t${r_val}\t${maf_val}\t${loci}\t${variants}\t${mean_dp}\t${mean_site_miss}\t${lowdp}\t${gt50miss}"
  echo -e "$line"
  [[ -n "$OUT_TSV" ]] && echo -e "$line" >> "$OUT_TSV"
done