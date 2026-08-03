#!/usr/bin/env bash
set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "${PROJECT_DIR}"
export PYTHONPATH="${PROJECT_DIR}/src:${PYTHONPATH:-}"
SPLIT_DIR="${SPLIT_DIR:?set SPLIT_DIR}"; RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/results}"
MODEL="${MODEL:-enet}"; NORMALIZATION="${NORMALIZATION:-}"; SELECT_K="${SELECT_K:-5000}"
SPLIT_NAME="$(basename "${SPLIT_DIR}")"; NORM_TAG="${NORMALIZATION:-fromcfg}"
OUT_DIR="${RESULTS_DIR}/${SPLIT_NAME}/${MODEL}_${NORM_TAG}_k${SELECT_K}"
python src/predict_transfer.py --split_dir "${SPLIT_DIR}" --models_dir "${OUT_DIR}/models" --out_dir "${OUT_DIR}"
