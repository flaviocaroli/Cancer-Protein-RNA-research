#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${PROJECT_DIR}"

DATA_DIR="${DATA_DIR:-${PROJECT_DIR}/data}"
RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/results}"
X_PATH="${X_PATH:-}"
Y_PATH="${Y_PATH:-}"

MODELS="${MODELS:-proxy,enet}"
NORMALIZATIONS="${NORMALIZATIONS:-zscore}"
SELECT_KS="${SELECT_KS:-2000}"
N_REPEATS="${N_REPEATS:-20}"
TEST_SIZE="${TEST_SIZE:-0.20}"
CV_FOLDS="${CV_FOLDS:-5}"
GRID_PRESET="${GRID_PRESET:-fast}"
MAX_PARAM_COMBOS="${MAX_PARAM_COMBOS:-8}"
KEEP_PREDICTIONS="${KEEP_PREDICTIONS:-0}"
SEED="${SEED:-42}"
N_JOBS="${N_JOBS:-1}"
SKIP_EXISTING="${SKIP_EXISTING:-1}"

# TASK_MODE=packed_by_protein is best for clusters with low job limits: one array task per protein chunk.
# TASK_MODE=grid creates many smaller tasks, mainly useful on clusters without array-size limits.
TASK_MODE="${TASK_MODE:-packed_by_protein}"
PROTEIN_CHUNKS="${PROTEIN_CHUNKS:-20}"
PROTEINS_PER_TASK="${PROTEINS_PER_TASK:-50}"
SPLITS_PER_TASK="${SPLITS_PER_TASK:-5}"

SUBMIT="${SUBMIT:-slurm}"
PARTITION="${PARTITION:-compute}"
TIME="${TIME:-2-00:00:00}"
MEM="${MEM:-20G}"
CPUS_PER_TASK="${CPUS_PER_TASK:-1}"
ARRAY_THROTTLE="${ARRAY_THROTTLE:-20}"

mkdir -p "${RESULTS_DIR}/logs" "${RESULTS_DIR}/split_manifest" "${RESULTS_DIR}/task_tables" "${RESULTS_DIR}/metrics/metrics_shards" "${RESULTS_DIR}/prediction_shards"

X_ARG=()
Y_ARG=()
if [[ -n "${X_PATH}" ]]; then X_ARG=(--x_path "${X_PATH}"); fi
if [[ -n "${Y_PATH}" ]]; then Y_ARG=(--y_path "${Y_PATH}"); fi

echo "Project directory: ${PROJECT_DIR}"
echo "Data directory:    ${DATA_DIR}"
echo "Results directory: ${RESULTS_DIR}"
echo "Submit mode:       ${SUBMIT}"
echo "Task mode:         ${TASK_MODE}"
echo "Protein chunks:    ${PROTEIN_CHUNKS}"
echo "Models:            ${MODELS}"
echo "Normalizations:    ${NORMALIZATIONS}"
echo "Select k:          ${SELECT_KS}"
echo "Repeats:           ${N_REPEATS}"

echo ""
echo "[1/3] Creating repeated 80/20 splits"
python src/make_splits.py \
  --data_dir "${DATA_DIR}" \
  "${X_ARG[@]}" \
  "${Y_ARG[@]}" \
  --out_dir "${RESULTS_DIR}/split_manifest" \
  --n_repeats "${N_REPEATS}" \
  --test_size "${TEST_SIZE}" \
  --seed "${SEED}"

echo ""
echo "[2/3] Creating task table"
python src/make_task_table.py \
  --data_dir "${DATA_DIR}" \
  "${X_ARG[@]}" \
  "${Y_ARG[@]}" \
  --out_path "${RESULTS_DIR}/task_tables/tasks.tsv" \
  --models "${MODELS}" \
  --normalizations "${NORMALIZATIONS}" \
  --select_ks "${SELECT_KS}" \
  --n_repeats "${N_REPEATS}" \
  --task_mode "${TASK_MODE}" \
  --protein_chunks "${PROTEIN_CHUNKS}" \
  --proteins_per_task "${PROTEINS_PER_TASK}" \
  --splits_per_task "${SPLITS_PER_TASK}" \
  --cv_folds "${CV_FOLDS}" \
  --grid_preset "${GRID_PRESET}" \
  --max_param_combos "${MAX_PARAM_COMBOS}" \
  --keep_predictions "${KEEP_PREDICTIONS}"

N_TASKS=$(($(wc -l < "${RESULTS_DIR}/task_tables/tasks.tsv") - 1))
echo ""
echo "[3/3] Task count: ${N_TASKS}"

if [[ "${N_TASKS}" -le 0 ]]; then
  echo "No tasks were created. Exiting." >&2
  exit 1
fi

if [[ "${SUBMIT}" == "slurm" ]]; then
  echo "Submitting one SLURM array with ${N_TASKS} tasks"
  echo "Array throttle: ${ARRAY_THROTTLE} concurrent tasks"

  sbatch \
    --job-name=robust8020 \
    --partition="${PARTITION}" \
    --time="${TIME}" \
    --mem="${MEM}" \
    --cpus-per-task="${CPUS_PER_TASK}" \
    --array=1-"${N_TASKS}"%"${ARRAY_THROTTLE}" \
    --output="${RESULTS_DIR}/logs/robust8020_%A_%a.out" \
    --error="${RESULTS_DIR}/logs/robust8020_%A_%a.err" \
    --wrap="cd ${PROJECT_DIR} && python src/run_task.py \
      --task_table ${RESULTS_DIR}/task_tables/tasks.tsv \
      --task_id \$SLURM_ARRAY_TASK_ID \
      --data_dir ${DATA_DIR} \
      --split_manifest ${RESULTS_DIR}/split_manifest/repeated_8020_splits.tsv \
      --out_dir ${RESULTS_DIR} \
      --seed ${SEED} \
      --n_jobs ${N_JOBS} \
      --skip_existing ${SKIP_EXISTING} \
      ${X_PATH:+--x_path ${X_PATH}} \
      ${Y_PATH:+--y_path ${Y_PATH}}"

elif [[ "${SUBMIT}" == "local" ]]; then
  echo "Running locally. This is intended only for smoke tests or small jobs."
  for TASK_ID in $(seq 1 "${N_TASKS}"); do
    echo "Running task ${TASK_ID}/${N_TASKS}"
    python src/run_task.py \
      --task_table "${RESULTS_DIR}/task_tables/tasks.tsv" \
      --task_id "${TASK_ID}" \
      --data_dir "${DATA_DIR}" \
      --split_manifest "${RESULTS_DIR}/split_manifest/repeated_8020_splits.tsv" \
      --out_dir "${RESULTS_DIR}" \
      --seed "${SEED}" \
      --n_jobs "${N_JOBS}" \
      --skip_existing "${SKIP_EXISTING}" \
      "${X_ARG[@]}" \
      "${Y_ARG[@]}"
  done
else
  echo "Unknown SUBMIT=${SUBMIT}. Use SUBMIT=slurm or SUBMIT=local." >&2
  exit 1
fi

echo ""
echo "run_repeated.sh finished launcher step."
echo "After SLURM jobs finish, run run_after_jobs.sh."
