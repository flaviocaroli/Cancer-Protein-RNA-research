#!/usr/bin/env bash
set -euo pipefail

SCENARIOS=(
  ccle_tmt_to_itraq_labelfree
  ccle_tmt_itraq_to_labelfree
  ccle_tmt_labelfree_to_itraq
  ccle_labelfree_itraq_to_tmt
)
NORMALIZATIONS_CSV="${NORMALIZATIONS:-zscore,rank,binning}"
IFS=',' read -r -a NORMALIZATIONS_ARRAY <<< "${NORMALIZATIONS_CSV}"

for scenario in "${SCENARIOS[@]}"; do
  for norm in "${NORMALIZATIONS_ARRAY[@]}"; do
    echo "Preparing ${scenario} with ${norm}"
    SCENARIO="${scenario}" X_NORM="${norm}" bash scripts/run_prepare.sh
  done
done
