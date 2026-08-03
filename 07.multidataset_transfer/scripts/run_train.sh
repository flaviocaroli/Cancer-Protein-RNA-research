#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PROJECT_DIR}"

SCENARIO="${SCENARIO:-ccle_tmt_to_itraq_labelfree}"
X_NORM="${X_NORM:-zscore}"
Y_NORM="${Y_NORM:-domain_zscore}"
MODEL="${MODEL:-enet}"
PREPARED_DIR="${PREPARED_DIR:-${PROJECT_DIR}/prepared/${SCENARIO}/x-${X_NORM}__y-${Y_NORM}}"
RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/results/${SCENARIO}/x-${X_NORM}__y-${Y_NORM}/${MODEL}}"
MODEL_CONFIG="${MODEL_CONFIG:-${PROJECT_DIR}/configs/model_defaults.yaml}"
PROTEIN_CHUNKS="${PROTEIN_CHUNKS:-20}"
PROTEIN_SUBSET="${PROTEIN_SUBSET:-}"
MAX_HPC_JOBS="${MAX_HPC_JOBS:-20}"

SUBMIT="${SUBMIT:-slurm}"
PARTITION="${PARTITION:-compute}"
TIME="${TIME:-3-00:00:00}"
MEM="${MEM:-24G}"
CPUS_PER_TASK="${CPUS_PER_TASK:-1}"
N_JOBS="${N_JOBS:-1}"
ARRAY_THROTTLE="${ARRAY_THROTTLE:-20}"
SEED="${SEED:-42}"
SAVE_PREDICTIONS="${SAVE_PREDICTIONS:-1}"
SAVE_ALL_FEATURE_ROWS="${SAVE_ALL_FEATURE_ROWS:-0}"
SKIP_EXISTING="${SKIP_EXISTING:-1}"
PYTHON_EXECUTABLE="${PYTHON_EXECUTABLE:-$(command -v python)}"

if [[ ! -f "${PREPARED_DIR}/manifest.json" ]]; then
  echo "Prepared bundle not found: ${PREPARED_DIR}" >&2
  echo "Run notebook 02 or scripts/run_prepare.sh first." >&2
  exit 1
fi
if [[ "${PROTEIN_CHUNKS}" -gt "${MAX_HPC_JOBS}" ]]; then
  echo "PROTEIN_CHUNKS=${PROTEIN_CHUNKS} exceeds the hard HPC limit MAX_HPC_JOBS=${MAX_HPC_JOBS}." >&2
  exit 1
fi
if [[ "${ARRAY_THROTTLE}" -gt "${MAX_HPC_JOBS}" ]]; then
  ARRAY_THROTTLE="${MAX_HPC_JOBS}"
fi

mkdir -p "${RESULTS_DIR}/logs" "${RESULTS_DIR}/task_tables"
TASK_TABLE="${RESULTS_DIR}/task_tables/tasks.tsv"

MAKE_TASK_ARGS=(
  --prepared-dir "${PREPARED_DIR}"
  --out "${TASK_TABLE}"
  --protein-chunks "${PROTEIN_CHUNKS}"
)
if [[ -n "${PROTEIN_SUBSET}" ]]; then
  MAKE_TASK_ARGS+=(--protein-subset "${PROTEIN_SUBSET}")
fi
"${PYTHON_EXECUTABLE}" src/make_tasks.py "${MAKE_TASK_ARGS[@]}"

N_TASKS=$(( $(wc -l < "${TASK_TABLE}") - 1 ))
if [[ "${N_TASKS}" -lt 1 ]]; then
  echo "No tasks generated." >&2
  exit 1
fi
if [[ "${N_TASKS}" -gt "${MAX_HPC_JOBS}" ]]; then
  echo "Refusing to submit ${N_TASKS} tasks; cluster limit is ${MAX_HPC_JOBS}." >&2
  exit 1
fi

COMMON_ARGS=(
  --prepared-dir "${PREPARED_DIR}"
  --results-dir "${RESULTS_DIR}"
  --task-table "${TASK_TABLE}"
  --model "${MODEL}"
  --model-config "${MODEL_CONFIG}"
  --n-jobs "${N_JOBS}"
  --seed "${SEED}"
  --save-predictions "${SAVE_PREDICTIONS}"
  --save-all-feature-rows "${SAVE_ALL_FEATURE_ROWS}"
  --skip-existing "${SKIP_EXISTING}"
)
if [[ -n "${PROTEIN_SUBSET}" ]]; then
  COMMON_ARGS+=(--protein-subset "${PROTEIN_SUBSET}")
fi

printf -v COMMON_ARGS_ESCAPED ' %q' "${COMMON_ARGS[@]}"

if [[ "${SUBMIT}" == "slurm" ]]; then
  sbatch \
    --job-name="md_${MODEL}" \
    --partition="${PARTITION}" \
    --time="${TIME}" \
    --mem="${MEM}" \
    --cpus-per-task="${CPUS_PER_TASK}" \
    --array="1-${N_TASKS}%${ARRAY_THROTTLE}" \
    --output="${RESULTS_DIR}/logs/${MODEL}_%A_%a.out" \
    --error="${RESULTS_DIR}/logs/${MODEL}_%A_%a.err" \
    --wrap="cd $(printf '%q' "${PROJECT_DIR}") && $(printf '%q' "${PYTHON_EXECUTABLE}") src/train_task.py${COMMON_ARGS_ESCAPED} --task-id \$SLURM_ARRAY_TASK_ID"
elif [[ "${SUBMIT}" == "local" ]]; then
  for TASK_ID in $(seq 1 "${N_TASKS}"); do
    "${PYTHON_EXECUTABLE}" src/train_task.py "${COMMON_ARGS[@]}" --task-id "${TASK_ID}"
  done
else
  echo "SUBMIT must be slurm or local" >&2
  exit 1
fi

cat <<EOF
Submitted/running ${N_TASKS} ${MODEL} jobs (never more than ${MAX_HPC_JOBS}).
Each job loads the prepared bundle once and loops sequentially over its assigned protein chunk.
Prepared bundle: ${PREPARED_DIR}
Results: ${RESULTS_DIR}
After completion run scripts/run_aggregate.sh with the same scenario/norm/model.
EOF
