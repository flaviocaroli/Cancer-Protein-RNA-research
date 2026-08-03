#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PROJECT_DIR}"

SCENARIO="${SCENARIO:-ccle_tmt_to_itraq_labelfree}"
TRAIN_DATASETS="${TRAIN_DATASETS:-}"
TEST_DATASETS="${TEST_DATASETS:-}"
OUTPUT_NAME="${OUTPUT_NAME:-}"
X_NORM="${X_NORM:-zscore}"
Y_NORM="${Y_NORM:-domain_zscore}"
N_BINS="${N_BINS:-50}"
EPS="${EPS:-1e-8}"
DOMAIN_BALANCED_NORMALIZER="${DOMAIN_BALANCED_NORMALIZER:-1}"
QUANTILE_FIT_MAX_SAMPLES_PER_DOMAIN="${QUANTILE_FIT_MAX_SAMPLES_PER_DOMAIN:-1000}"
MAKE_FIGURES="${MAKE_FIGURES:-1}"
FORCE="${FORCE:-0}"
PYTHON_EXECUTABLE="${PYTHON_EXECUTABLE:-$(command -v python)}"

ARGS=(
  --project-root "${PROJECT_DIR}"
  --scenario "${SCENARIO}"
  --x-norm "${X_NORM}"
  --y-norm "${Y_NORM}"
  --n-bins "${N_BINS}"
  --eps "${EPS}"
  --domain-balanced-normalizer "${DOMAIN_BALANCED_NORMALIZER}"
  --quantile-fit-max-samples-per-domain "${QUANTILE_FIT_MAX_SAMPLES_PER_DOMAIN}"
  --make-figures "${MAKE_FIGURES}"
)
if [[ -n "${TRAIN_DATASETS}" ]]; then ARGS+=(--train-datasets "${TRAIN_DATASETS}"); fi
if [[ -n "${TEST_DATASETS}" ]]; then ARGS+=(--test-datasets "${TEST_DATASETS}"); fi
if [[ -n "${OUTPUT_NAME}" ]]; then ARGS+=(--output-name "${OUTPUT_NAME}"); fi
if [[ "${FORCE}" == "1" ]]; then ARGS+=(--force); fi

"${PYTHON_EXECUTABLE}" src/prepare_scenario.py "${ARGS[@]}"
