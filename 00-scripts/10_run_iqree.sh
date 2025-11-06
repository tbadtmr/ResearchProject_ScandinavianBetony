#!/usr/bin/env bash
# ==============================================================
# Script: 10_run_iqtree.sh
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Convert a filtered VCF (including the outgroup) to PHYLIP format
#   using vcf2phylip.py, then run IQ-TREE v2 to infer a maximum-
#   likelihood phylogeny with ascertainment-bias correction (+ASC),
#   1,000 ultrafast bootstrap replicates, and SH-aLRT branch supports.
#
# Inputs:
#   - VCF: 03-analysis/07-finalvcfs/<RUN>_filtered.recode.vcf
#   - Outgroup: sample ID (e.g. "Trento_B_hirsuta")
#
# Outputs:
#   - Tree: 03-analysis/08-downstream/03-phylogeny/<RUN>/<RUN>.treefile
#   - Log:  03-analysis/08-downstream/03-phylogeny/<RUN>/<RUN>.log
#   - If invariant sites detected: re-runs on <RUN>.varsites.phy
#
# Usage:
#   OUTGROUP="Trento_B_hirsuta" THREADS=8 ./00-scripts/10_run_iqtree.sh r70_noTrentino_noOutgroup
# ==============================================================

set -euo pipefail

# USER PARAMETERS
RUN="${1:-r70_noTrentino_noOutgroup}"                      # dataset tag
VCF="03-analysis/07-finalvcfs/${RUN}_filtered.recode.vcf"  # filtered VCF WITH outgroup
OUTGROUP="${OUTGROUP:-Trento_B_hirsuta}"                   # outgroup sample ID
THREADS="${THREADS:-8}"                                    # CPU threads
IQTREE_BIN="${IQTREE_BIN:-iqtree2}"                        # IQ-TREE binary name (v2 recommended)

# Paths & setup
ADIR="03-analysis/08-downstream/03-phylogeny/${RUN}"
mkdir -p "${ADIR}"
PHY="${ADIR}/${RUN}.phy"
FINAL_PREFIX="${ADIR}/${RUN}"

# vcf2phylip setup
V2P="00-scripts/vcf2phylip.py"
if [[ ! -f "${V2P}" ]]; then
  echo "Downloading vcf2phylip.py..."
  curl -L -o "${V2P}" https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py
  chmod +x "${V2P}"
fi

# Convert VCF -> PHYLIP
# Keeps only variable sites (typical for RADseq SNP data)
python "${V2P}" -i "${VCF}" -o "${ADIR}" --phylip

# Move output to consistent filename (<RUN>.phy)
mv "${ADIR}/"*.phy "${PHY}"

# RUN IQ-TREE
# Uses model selection (MFP) with ascertainment bias correction (+ASC)
# and estimates support with 1,000 ultrafast bootstrap + SH-aLRT replicates
"${IQTREE_BIN}" \
  -s "${PHY}" \
  -st DNA \
  -m MFP+ASC \
  -B 1000 -alrt 1000 \
  -o "${OUTGROUP}" \
  -nt "${THREADS}" \
  -pre "${FINAL_PREFIX}"

# Handle invariant sites
# If IQ-TREE detects invariant sites with +ASC, it writes a <prefix>.varsites.phy file.
VARSITES="${FINAL_PREFIX}.varsites.phy"
if [[ -f "${VARSITES}" ]]; then
  echo "Detected invariant sites; re-running IQ-TREE on variable-sites alignment..."
  "${IQTREE_BIN}" \
    -s "${VARSITES}" \
    -st DNA \
    -m MFP+ASC \
    -B 1000 -alrt 1000 \
    -o "${OUTGROUP}" \
    -nt "${THREADS}" \
    -pre "${FINAL_PREFIX}"
fi

# Completion message
echo "IQ-TREE run completed."
echo "Tree file: ${FINAL_PREFIX}.treefile"
echo "Support values: UFBoot + SH-aLRT in ${FINAL_PREFIX}.iqtree"
