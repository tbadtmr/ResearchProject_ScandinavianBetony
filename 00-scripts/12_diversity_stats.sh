#!/usr/bin/env bash
# ==============================================================
# Script: 12_diversity_stats.sh
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Compute per-individual heterozygosity (Ho/He/F), per-group summaries,
#   and pairwise FST across groups for two popmap scenarios using VCFtools.
#
# Inputs:
#   - VCF (LD-pruned, no outgroup):
#       03-analysis/07-finalvcfs/<RUN>_pruned.vcf[.gz]
#   - Popmaps (two-column: SAMPLE<TAB>GROUP):
#       01-info_files/popmap_groups.txt
#       01-info_files/popmap_groups_scania.txt
#
# Outputs:
#   - Full: 03-analysis/08-downstream/04-diversity/
#   - Copies of key summaries: 04-results/
#
# Usage:
#   ./00-scripts/12_diversity_stats.sh r70_noTrentino_noOutgroup
# ==============================================================

set -euo pipefail

RUN="${1:-r70_noTrentino_noOutgroup}"
VCF="03-analysis/07-finalvcfs/${RUN}_pruned.vcf"
POPMAP_ALL="01-info_files/popmap_groups.txt"
POPMAP_SPLIT="01-info_files/popmap_groups_scania.txt"

DDIR="03-analysis/08-downstream/04-diversity"
RDIR="04-results"
mkdir -p "${DDIR}" "${RDIR}"

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found"; exit 1; }; }
need vcftools; need awk; need grep; need sed

vcf_flag="--vcf"
[[ "${VCF}" =~ \.vcf\.gz$ ]] && vcf_flag="--gzvcf"

header_line() {
  if [[ "${VCF}" =~ \.vcf\.gz$ ]]; then zgrep -m1 '^#CHROM' "${VCF}"
  else grep -m1 '^#CHROM' "${VCF}"; fi
}

export_case_summaries() {
  local LABEL="$1" OUTDIR="$2"
  # copy concise, well-named tables into 04-results
  cp -f "${OUTDIR}/indiv_HoHe_from_vcftools.tsv"       "${RDIR}/diversity_individual_${LABEL}.tsv"
  cp -f "${OUTDIR}/group_heterozygosity_summary.tsv"   "${RDIR}/diversity_group_summary_${LABEL}.tsv"
  cp -f "${OUTDIR}/group_heterozygosity_with_SE.tsv"   "${RDIR}/diversity_group_SE_${LABEL}.tsv"
  cp -f "${OUTDIR}/FST_pairwise_summary.tsv"           "${RDIR}/fst_pairwise_summary_${LABEL}.tsv"
}

run_case() {
  local LABEL="$1"   # allGroups / scaniaSplit
  local POPMAP="$2"

  echo "Case: ${LABEL}"

  local OUTDIR="${DDIR}/stats_${LABEL}"
  local GROUPDIR="${DDIR}/groups_${LABEL}"
  mkdir -p "${OUTDIR}" "${GROUPDIR}"

  # Clean popmap (CRLF, extra columns)
  sed 's/\r$//' "${POPMAP}" | awk -F'\t' 'NF>=2{print $1"\t"$2}' > "${OUTDIR}/popmap.clean.tsv"

  # Samples present in this VCF
  header_line | tr '\t' '\n' | tail -n +10 > "${OUTDIR}/samples_in_vcf.txt"

  # Per-group sample lists for samples actually present in the VCF
  awk -v GD="${GROUPDIR}" 'BEGIN{FS=OFS="\t"}
    NR==FNR { inv[$1]=1; next }
    ($1 in inv) {
      g=$2; gsub(/[^A-Za-z0-9._-]+/,"_",g);
      print $1 >> (GD "/" g ".txt")
    }' "${OUTDIR}/samples_in_vcf.txt" "${OUTDIR}/popmap.clean.tsv"

  find "${GROUPDIR}" -type f -size 0 -delete || true

  # --het per individual
  vcftools ${vcf_flag} "${VCF}" --het --out "${OUTDIR}/indiv" >/dev/null

  # Map INDV -> GROUP
  awk '{grp[$1]=$2} END{for(k in grp) print k, grp[k]}' "${OUTDIR}/popmap.clean.tsv" > "${OUTDIR}/map.tmp"

  # Compute Ho, He, keep F from vcftools, and join GROUP
  awk 'NR==1{next}{
    INDV=$1; O=$2; E=$3; N=$4; F=$5;
    Ho = 1 - (O/N);
    He = 1 - (E/N);
    Fcalc = (He>0 ? 1 - (Ho/He) : "NA");
    print INDV, Ho, He, F, Fcalc
  }' OFS='\t' "${OUTDIR}/indiv.het" \
  | awk 'BEGIN{OFS="\t"} NR==FNR{g[$1]=$2; next} {print $1,(g[$1]?g[$1]:"NA"),$2,$3,$4,$5}' "${OUTDIR}/map.tmp" - \
  | awk 'BEGIN{OFS="\t"; print "INDV","GROUP","Ho","He","F_vcftools","Fcalc"}1' \
  > "${OUTDIR}/indiv_HoHe_from_vcftools.tsv"

  # Per-group means
  awk 'BEGIN{FS=OFS="\t"} NR==1{next}{
    g=$2; if(g=="NA") next;
    n[g]++; sumHo[g]+=$3; sumHe[g]+=$4; sumF[g]+=$5; sumFc[g]+=$6
  }
  END{
    print "Group","N","mean_Ho","mean_He","mean_F_vcftools","mean_Fcalc";
    for(g in n){
      printf "%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\n",
        g, n[g], sumHo[g]/n[g], sumHe[g]/n[g], sumF[g]/n[g], sumFc[g]/n[g]
    }
  }' "${OUTDIR}/indiv_HoHe_from_vcftools.tsv" \
  | { read -r H; echo "$H"; sort -k1,1; } \
  > "${OUTDIR}/group_heterozygosity_summary.tsv"

  # With SE
  awk 'BEGIN{FS=OFS="\t"} NR==1{next}{
    g=$2; if(g=="NA") next;
    n[g]++; Ho[g]+=$3; He[g]+=$4; Fv[g]+=$5;
    Ho2[g]+=$3*$3; He2[g]+=$4*$4; F2[g]+=$5*$5
  }
  END{
    print "Group","N","mean_Ho","SE_Ho","mean_He","SE_He","mean_F","SE_F";
    for(g in n){
      mHo=Ho[g]/n[g]; mHe=He[g]/n[g]; mF=Fv[g]/n[g];
      seHo=sqrt((Ho2[g]/n[g]-mHo*mHo)/n[g]);
      seHe=sqrt((He2[g]/n[g]-mHe*mHe)/n[g]);
      seF =sqrt((F2[g] /n[g]-mF *mF )/n[g]);
      printf "%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
        g, n[g], mHo, seHo, mHe, seHe, mF, seF
    }
  }' "${OUTDIR}/indiv_HoHe_from_vcftools.tsv" \
  > "${OUTDIR}/group_heterozygosity_with_SE.tsv"

  # Pairwise FST
  mapfile -t groups_arr < <(ls "${GROUPDIR}"/*.txt | sort)
  for ((i=0;i<${#groups_arr[@]}-1;i++)); do
    for ((j=i+1;j<${#groups_arr[@]};j++)); do
      a=${groups_arr[$i]}; b=${groups_arr[$j]}
      an=$(basename "$a" .txt); bn=$(basename "$b" .txt)
      vcftools ${vcf_flag} "${VCF}" \
        --weir-fst-pop "$a" \
        --weir-fst-pop "$b" \
        --out "${OUTDIR}/FST_${an}_vs_${bn}" >/dev/null
    done
  done

  # Summarize weighted mean FST (from .log)
  grep -H "Weir and Cockerham mean Fst" "${OUTDIR}"/FST_*log \
  | sed 's#^.*/FST_##; s/\.log:Weir and Cockerham mean Fst estimate: /\t/; s/_vs_/\t/' \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' \
  | sort > "${OUTDIR}/FST_pairwise_summary.tsv"

  # Export concise copies into 04-results
  export_case_summaries "${LABEL}" "${OUTDIR}"

  echo "Done: ${LABEL}"
}

run_case "allGroups"   "${POPMAP_ALL}"
run_case "scaniaSplit" "${POPMAP_SPLIT}"

echo "Full outputs: ${DDIR}"
echo "Key summaries: ${RDIR}"
