#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${PROJECT_DIR}"

DATA_DIR="${DATA_DIR:-${PROJECT_DIR}/data}"
RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/results}"
X_PATH="${X_PATH:-}"
Y_PATH="${Y_PATH:-}"
COMBINE_PREDICTIONS="${COMBINE_PREDICTIONS:-0}"

# Choose the final model configs used for the robustness atlas.
PROXY_CONFIG="${PROXY_CONFIG:-proxy:zscore:0}"
ENET_CONFIG="${ENET_CONFIG:-enet:zscore:2000}"
XGB_CONFIG="${XGB_CONFIG:-}"
METRIC="${METRIC:-pearson}"
THRESHOLD="${THRESHOLD:-0.50}"
PASS_RATE="${PASS_RATE:-0.80}"
GAIN_THRESHOLD="${GAIN_THRESHOLD:-0.05}"
N_BOOT="${N_BOOT:-2000}"
REQUIRE_CI_LOW="${REQUIRE_CI_LOW:-0}"

X_ARG=()
Y_ARG=()
if [[ -n "${X_PATH}" ]]; then X_ARG=(--x_path "${X_PATH}"); fi
if [[ -n "${Y_PATH}" ]]; then Y_ARG=(--y_path "${Y_PATH}"); fi

echo "[1/2] Aggregating metrics"
python src/aggregate.py \
  --results_dir "${RESULTS_DIR}" \
  --data_dir "${DATA_DIR}" \
  "${X_ARG[@]}" \
  "${Y_ARG[@]}" \
  --combine_predictions "${COMBINE_PREDICTIONS}"

METRICS_LONG="${RESULTS_DIR}/metrics/metrics_long.parquet"
if [[ ! -f "${METRICS_LONG}" ]]; then
  if [[ -f "${RESULTS_DIR}/metrics/metrics_long.csv.gz" ]]; then
    METRICS_LONG="${RESULTS_DIR}/metrics/metrics_long.csv.gz"
  elif [[ -f "${RESULTS_DIR}/metrics/metrics_long.tsv.gz" ]]; then
    METRICS_LONG="${RESULTS_DIR}/metrics/metrics_long.tsv.gz"
  fi
fi

echo "[2/2] Building robust cascade atlas"
XGB_ARG=()
if [[ -n "${XGB_CONFIG}" ]]; then XGB_ARG=(--xgb_config "${XGB_CONFIG}"); fi
python src/bootstrap_cascade.py \
  --metrics_long "${METRICS_LONG}" \
  --out_dir "${RESULTS_DIR}/cascade" \
  --proxy_config "${PROXY_CONFIG}" \
  --enet_config "${ENET_CONFIG}" \
  "${XGB_ARG[@]}" \
  --metric "${METRIC}" \
  --threshold "${THRESHOLD}" \
  --pass_rate "${PASS_RATE}" \
  --gain_threshold "${GAIN_THRESHOLD}" \
  --n_boot "${N_BOOT}" \
  --require_ci_low "${REQUIRE_CI_LOW}"
