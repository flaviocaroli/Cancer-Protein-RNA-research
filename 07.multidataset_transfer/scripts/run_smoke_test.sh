#!/usr/bin/env bash
set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PROJECT_DIR}"
PYTHON_EXECUTABLE="${PYTHON_EXECUTABLE:-$(command -v python)}"

rm -rf data/Collection_A data/processed data/reference prepared results/smoke figures/preparation/ccle_tmt_to_itraq_labelfree
"${PYTHON_EXECUTABLE}" src/make_synthetic_data.py --project-root "${PROJECT_DIR}" --genes 50 --proteins 4

MAPPING_PATH="${PROJECT_DIR}/data/reference/ensembl_to_symbol.tsv" FORCE=1 bash scripts/run_collection_loader.sh
SCENARIO=ccle_tmt_to_itraq_labelfree X_NORM=zscore Y_NORM=domain_zscore N_BINS=10 MAKE_FIGURES=1 FORCE=1 bash scripts/run_prepare.sh

cat > configs/model_smoke.yaml <<'YAML'
training:
  feature_k: 10
  cv_folds: 2
  cv_metric: domain_mean_pearson
  cv_mode: stratified_domain
  feature_score_mode: domain_rank_mean
  domain_weighting: balanced
  random_seed: 42
  max_param_combos: 1
  enet_grid:
    alpha: [0.1]
    l1_ratio: [0.5]
    max_iter: 2000
    tol: 0.0001
  xgb_grid:
    n_estimators: [10]
    max_depth: [2]
    learning_rate: [0.1]
    min_child_weight: [1]
    subsample: [0.8]
    colsample_bytree: [0.8]
    reg_lambda: [1.0]
    reg_alpha: [0.0]
YAML

for model in enet xgb; do
  rm -rf "results/smoke/${model}"
  SCENARIO=ccle_tmt_to_itraq_labelfree X_NORM=zscore Y_NORM=domain_zscore MODEL="${model}" \
  RESULTS_DIR="${PROJECT_DIR}/results/smoke/${model}" MODEL_CONFIG="${PROJECT_DIR}/configs/model_smoke.yaml" \
  PROTEIN_CHUNKS=2 MAX_HPC_JOBS=20 SUBMIT=local SAVE_PREDICTIONS=1 \
  SAVE_ALL_FEATURE_ROWS=1 bash scripts/run_train.sh
  RESULTS_DIR="${PROJECT_DIR}/results/smoke/${model}" MODEL="${model}" bash scripts/run_aggregate.sh
done

MODEL_FILE=$(find results/smoke/enet/models_per_protein/enet -name '*.joblib' | head -1)
"${PYTHON_EXECUTABLE}" src/predict_saved_model.py \
  --model "${MODEL_FILE}" \
  --rna prepared/ccle_tmt_to_itraq_labelfree/x-zscore__y-domain_zscore/tests/iTRAQ/X_raw.joblib \
  --out results/smoke/reused_model_predictions.parquet

echo "Smoke test completed successfully."
