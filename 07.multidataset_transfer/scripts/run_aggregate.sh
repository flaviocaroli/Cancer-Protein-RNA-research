#!/usr/bin/env bash
set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PROJECT_DIR}"
SCENARIO="${SCENARIO:-ccle_tmt_to_itraq_labelfree}"
X_NORM="${X_NORM:-zscore}"
Y_NORM="${Y_NORM:-domain_zscore}"
MODEL="${MODEL:-enet}"
RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/results/${SCENARIO}/x-${X_NORM}__y-${Y_NORM}/${MODEL}}"
python src/aggregate_results.py --results-dir "${RESULTS_DIR}" --model "${MODEL}"
