#!/usr/bin/env bash
set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "${PROJECT_DIR}"
export PYTHONPATH="${PROJECT_DIR}/src:${PYTHONPATH:-}"
# Build a split via notebook 02 logic in headless form (synthetic fallback data).
python notebooks/_smoke_build_split.py
SPLIT_DIR=data/splits/ccle_tmt_itraq__labelfree MODEL=enet SELECT_K=120 CV_FOLDS=5 \
  MAX_COMBOS=4 PROTEIN_CHUNKS=4 SUBMIT=local bash run_train.sh
SPLIT_DIR=data/splits/ccle_tmt_itraq__labelfree MODEL=enet SELECT_K=120 bash run_predict.sh
echo "Smoke test OK."
