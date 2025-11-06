#!/usr/bin/env bash
# ==============================================================
# Script: 12_make_splitstree_input.sh
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Convert an LD-pruned VCF (without outgroup) to PHYLIP format
#   for NeighborNet in SplitsTree.
#
# Inputs:
#   - VCF: 03-analysis/07-finalvcfs/<RUN>_pruned.vcf
#
# Outputs:
#   - PHYLIP: 03-analysis/08-downstream/03-phylogeny/<RUN>_splitstree.phy
#
# Usage:
#   ./00-scripts/12_make_splitstree_input.sh r70_noTrentino_noOutgroup
# ==============================================================

set -euo pipefail

RUN="${1:-r70_noTrentino_noOutgroup}"
VCF="03-analysis/07-finalvcfs/${RUN}_pruned.vcf"
ODIR="03-analysis/08-downstream/03-phylogeny/
OFN="${RUN}_splitstree.phy"
mkdir -p "${ODIR}"

# Get vcf2phylip if missing
V2P="00-scripts/vcf2phylip.py"
if [[ ! -f "${V2P}" ]]; then
  echo "Downloading vcf2phylip.py..."
  curl -L -o "${V2P}" https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py
  chmod +x "${V2P}"
fi

# Convert VCF -> PHYLIP (SplitsTree will compute p-distances from the alignment)
python "${V2P}" -i "${VCF}" -o "${ODIR}" --phylip

# Normalize filename
PHY_FOUND=$(ls -1 "${ODIR}"/*.phy | head -n1)
mv -f "${PHY_FOUND}" "${ODIR}/${OFN}"
echo "PHYLIP ready: ${ODIR}/${OFN}"