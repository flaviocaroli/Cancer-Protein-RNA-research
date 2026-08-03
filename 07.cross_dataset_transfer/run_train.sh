#!/usr/bin/env bash
# Train per-protein models for ONE split x model x normalization.
# HPC-FRIENDLY: submits a SINGLE array of at most PROTEIN_CHUNKS (<=20) tasks.
# Each task handles one contiguous protein chunk; train_transfer.py loops over
# the proteins in that chunk. This respects a 20-concurrent-job cap.
set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "${PROJECT_DIR}"
export PYTHONPATH="${PROJECT_DIR}/src:${PYTHONPATH:-}"

SPLIT_DIR="${SPLIT_DIR:?set SPLIT_DIR=data/splits/<split_name>}"
RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/results}"
MODEL="${MODEL:-enet}"                       # enet | xgb
NORMALIZATION="${NORMALIZATION:-}"           # empty -> use split_config model_norm
SELECT_K="${SELECT_K:-5000}"; CV_FOLDS="${CV_FOLDS:-5}"; MAX_COMBOS="${MAX_COMBOS:-8}"; SEED="${SEED:-42}"

PROTEIN_CHUNKS="${PROTEIN_CHUNKS:-20}"       # HARD CAP: keep <= 20
[[ "${PROTEIN_CHUNKS}" -gt 20 ]] && { echo "PROTEIN_CHUNKS>20 exceeds job cap"; exit 1; }
ARRAY_THROTTLE="${ARRAY_THROTTLE:-20}"
PARTITION="${PARTITION:-compute}"; TIME="${TIME:-3-00:00:00}"; MEM="${MEM:-24G}"; CPUS="${CPUS:-1}"
SUBMIT="${SUBMIT:-slurm}"

# total protein count from the split
N_PROT="$(python -c "import sys; sys.path.insert(0,'src'); import io_utils as IO; print(len(IO.protein_list('${SPLIT_DIR}')))")"
CHUNK=$(( (N_PROT + PROTEIN_CHUNKS - 1) / PROTEIN_CHUNKS ))
SPLIT_NAME="$(basename "${SPLIT_DIR}")"
NORM_TAG="${NORMALIZATION:-fromcfg}"
OUT_DIR="${RESULTS_DIR}/${SPLIT_NAME}/${MODEL}_${NORM_TAG}_k${SELECT_K}"
mkdir -p "${OUT_DIR}/logs"
echo "Split=${SPLIT_NAME} model=${MODEL} norm=${NORM_TAG} k=${SELECT_K}"
echo "Proteins=${N_PROT} in ${PROTEIN_CHUNKS} chunks of ${CHUNK} -> ${OUT_DIR}"

NORM_ARG=(); [[ -n "${NORMALIZATION}" ]] && NORM_ARG=(--normalization "${NORMALIZATION}")

if [[ "${SUBMIT}" == "slurm" ]]; then
  sbatch --job-name="tr_${MODEL}" --partition="${PARTITION}" --time="${TIME}" \
    --mem="${MEM}" --cpus-per-task="${CPUS}" \
    --array=0-$(( PROTEIN_CHUNKS - 1 ))%"${ARRAY_THROTTLE}" \
    --output="${OUT_DIR}/logs/%A_%a.out" --error="${OUT_DIR}/logs/%A_%a.err" \
    --wrap="cd ${PROJECT_DIR} && export PYTHONPATH=${PROJECT_DIR}/src:\$PYTHONPATH && \
      i=\$SLURM_ARRAY_TASK_ID; s=\$(( i * ${CHUNK} )); e=\$(( s + ${CHUNK} )); \
      python src/train_transfer.py --split_dir ${SPLIT_DIR} --model ${MODEL} \
      ${NORMALIZATION:+--normalization ${NORMALIZATION}} --select_k ${SELECT_K} \
      --cv_folds ${CV_FOLDS} --max_combos ${MAX_COMBOS} --seed ${SEED} \
      --protein_start \${s} --protein_end \${e} --out_dir ${OUT_DIR}"
else
  # local: explicit loop over the (<=20) chunks, each looping its proteins internally
  for i in $(seq 0 $(( PROTEIN_CHUNKS - 1 ))); do
    s=$(( i * CHUNK )); e=$(( s + CHUNK ))
    python src/train_transfer.py --split_dir "${SPLIT_DIR}" --model "${MODEL}" \
      "${NORM_ARG[@]}" --select_k "${SELECT_K}" --cv_folds "${CV_FOLDS}" \
      --max_combos "${MAX_COMBOS}" --seed "${SEED}" \
      --protein_start "${s}" --protein_end "${e}" --out_dir "${OUT_DIR}"
  done
fi
