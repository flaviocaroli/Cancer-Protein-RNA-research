#!/bin/bash
# submit_chunks.sh
# Runs pipeline.py on proteins listed in ${DATA_DIR}/protein_list.txt using SLURM arrays.
# Saves outputs into ${DATA_DIR}/results (created if missing).

set -euo pipefail

CHUNK_SIZE=500
TOTAL_PROTEINS=5113   # or: TOTAL_PROTEINS=$(wc -l < "${DATA_DIR}/protein_list.txt")
DATA_DIR=~/pipeline_proteins/data
SRC_DIR=~/pipeline_proteins/src

# MUST match pipeline.py choices: fixed_percentile | fixed_p | union
PIPELINE_TYPE="fixed_percentile"

# Optional knobs passed to pipeline.py
PERCENTILE=5
P_VALUE=0.05
N_SPLITS=3

# Default pipeline.py behavior is save-only-best. Set FULL_OUTPUT=1 to add --full_output.
FULL_OUTPUT=0

MODEL_TYPE="elasticnet"

PROTEIN_LIST="${DATA_DIR}/protein_list.txt"
LOG_DIR="${DATA_DIR}/logs"
RESULTS_DIR="${DATA_DIR}/results"

mkdir -p "${LOG_DIR}"
mkdir -p "${RESULTS_DIR}"

if [ ! -f "${PROTEIN_LIST}" ]; then
  echo "Error: Protein list file not found at ${PROTEIN_LIST}"
  exit 1
fi

for START in $(seq 1 "${CHUNK_SIZE}" "${TOTAL_PROTEINS}"); do
  END=$((START + CHUNK_SIZE - 1))
  if [ "${END}" -gt "${TOTAL_PROTEINS}" ]; then
    END="${TOTAL_PROTEINS}"
  fi

  sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=prot_${PIPELINE_TYPE}_%A_%a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --array=${START}-${END}
#SBATCH --output=${LOG_DIR}/pipeline_%A_%a.out
#SBATCH --error=${LOG_DIR}/pipeline_%A_%a.err

set -euo pipefail

module load miniconda3
source activate base
conda activate test_env

cd "${SRC_DIR}"

# Ensure results folder exists on compute node too
mkdir -p "${RESULTS_DIR}"

# Get protein name for this SLURM array task id (1-indexed)
PROTEIN_NAME=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${PROTEIN_LIST}")

if [ -z "\${PROTEIN_NAME}" ]; then
  echo "Error: Empty protein name for SLURM_ARRAY_TASK_ID=\${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Running pipeline for protein: \${PROTEIN_NAME}"
echo "Pipeline: ${PIPELINE_TYPE} | percentile=${PERCENTILE} | p_value=${P_VALUE} | n_splits=${N_SPLITS} | full_output=${FULL_OUTPUT}"
echo "Saving outputs to: ${RESULTS_DIR}"

# Build args
ARGS=(
  --protein "\${PROTEIN_NAME}"
  --pipeline "${PIPELINE_TYPE}"
  --model "${MODEL_TYPE}"
  --percentile "${PERCENTILE}"
  --p_value "${P_VALUE}"
  --n_splits "${N_SPLITS}"
)

if [ "${FULL_OUTPUT}" -eq 1 ]; then
  ARGS+=(--full_output)
fi

# IMPORTANT: output_dir is controlled inside pipeline.py. To save into results/,
# pipeline.py must use output_dir=<.../results> when save_only_best=True.
# If your pipeline.py currently hardcodes output_dir=data_dir, update it to output_dir=results_dir.
python pipeline.py "\${ARGS[@]}"

echo "Done: \${PROTEIN_NAME}"
EOF

done
s