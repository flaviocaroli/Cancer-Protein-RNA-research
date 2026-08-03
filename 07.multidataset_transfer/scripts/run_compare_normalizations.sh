#!/usr/bin/env bash
set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PROJECT_DIR}"
SCENARIO="${SCENARIO:-ccle_tmt_to_itraq_labelfree}"
Y_NORM="${Y_NORM:-domain_zscore}"
MODEL="${MODEL:-enet}"
NORMALIZATIONS="${NORMALIZATIONS:-zscore,rank,binning}"
python src/compare_normalizations.py --project-root "${PROJECT_DIR}" --scenario "${SCENARIO}" --y-norm "${Y_NORM}" --model "${MODEL}" --normalizations "${NORMALIZATIONS}"
