#!/usr/bin/env bash
set -euo pipefail
SCENARIO="${SCENARIO:-ccle_tmt_to_itraq_labelfree}"
Y_NORM="${Y_NORM:-domain_zscore}"
for X_NORM in zscore rank binning; do
  echo "Preparing ${SCENARIO} / ${X_NORM}"
  SCENARIO="${SCENARIO}" X_NORM="${X_NORM}" Y_NORM="${Y_NORM}" bash "$(dirname "$0")/run_prepare.sh"
done
