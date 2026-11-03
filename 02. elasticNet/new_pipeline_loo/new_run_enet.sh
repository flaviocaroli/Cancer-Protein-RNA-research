#!/bin/bash
# new_run_enet.sh
# Submits ONE SLURM array job with N_CHUNKS tasks.

set -euo pipefail

# Adjust these to your install location on the cluster
BASE_DIR=${BASE_DIR:-$(pwd)}            # run from new_pipeline_loo/
DATA_DIR="${BASE_DIR}/../data"      # single 02. elasticNet/data folder
SRC_DIR="${BASE_DIR}/src"

PROTEIN_LIST="${DATA_DIR}/protein_list.txt"
NORMALIZATION_TYPE=${NORMALIZATION_TYPE:-zscore}   # zscore | none | ranking | binning
RNA_FILE=${RNA_FILE:-transcriptomics_filtered.joblib}
RESULTS_DIR=${RESULTS_DIR:-${DATA_DIR}/loo/results_${NORMALIZATION_TYPE}}
LOG_DIR="${DATA_DIR}/loo/logs"

mkdir -p "${RESULTS_DIR}" "${LOG_DIR}"

# ---------------------------
# Dynamic chunking controls
# ---------------------------
TARGET_CHUNKS=${TARGET_CHUNKS:-10}
# GROUP_SIZE_OVERRIDE=511

# ---------------------------
# 1. Feature Selection Strategy (PIPELINE)
# ---------------------------
# Choices: fixed_percentile | fixed_p | union | naive_cluster
PIPELINE_TYPE=${PIPELINE_TYPE:-fixed_percentile}

# ---------------------------
# 2. Normalization Method
# ---------------------------
# Choices: zscore | binning | ranking | none
# Use "zscore" for the standard pipeline.
# Use "binning" or "ranking" for Relative Representations (Stahlschmidt et al.).


# ---------------------------
# 3. Hyperparameters
# ---------------------------
PERCENTILE=${PERCENTILE:-10}
P_VALUE=0.05
N_SPLITS=5
N_BINS=${N_BINS:-50}     # Only used if NORMALIZATION_TYPE is "binning"

# Default is best-only. Set FULL_OUTPUT=1 to add --full_output.
FULL_OUTPUT=${FULL_OUTPUT:-1}   # 1 = also write <P>_LOO.parquet (needed for enet_all_results.h5)

if [ ! -f "${PROTEIN_LIST}" ]; then
  echo "Error: protein_list.txt not found at ${PROTEIN_LIST}"
  exit 1
fi

TOTAL_PROTEINS=$(grep -v '^#' "${PROTEIN_LIST}" | sed '/^$/d' | wc -l)
if [ "${TOTAL_PROTEINS}" -eq 0 ]; then
  echo "protein_list.txt is empty."
  exit 0
fi

# Compute dynamic GROUP_SIZE
if [ "${GROUP_SIZE_OVERRIDE:-0}" -gt 0 ]; then
  GROUP_SIZE="${GROUP_SIZE_OVERRIDE}"
else
  GROUP_SIZE=$(( (TOTAL_PROTEINS + TARGET_CHUNKS - 1) / TARGET_CHUNKS ))
  if [ "${GROUP_SIZE}" -lt 1 ]; then GROUP_SIZE=1; fi
fi

N_CHUNKS=$(( (TOTAL_PROTEINS + GROUP_SIZE - 1) / GROUP_SIZE ))
LAST_CHUNK=$(( N_CHUNKS - 1 ))

echo "Submitting array: 0-${LAST_CHUNK}"
echo "Pipeline: ${PIPELINE_TYPE} | Normalization: ${NORMALIZATION_TYPE}"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=np_${NORMALIZATION_TYPE}
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=11:00:00
#SBATCH --mem=32G
#SBATCH --array=0-${LAST_CHUNK}
#SBATCH --output=${LOG_DIR}/newpipe_%A_%a.out
#SBATCH --error=${LOG_DIR}/newpipe_%A_%a.err

set -euo pipefail

module load miniconda3
source activate base
conda activate ${CONDA_ENV:-test_env}

cd "${SRC_DIR}"
mkdir -p "${RESULTS_DIR}"

CHUNK_ID=\${SLURM_ARRAY_TASK_ID}

# Build base arguments
ARGS=(
  --chunk_id "\${CHUNK_ID}"
  --group_size "${GROUP_SIZE}"
  --pipeline "${PIPELINE_TYPE}"
  --normalization "${NORMALIZATION_TYPE}"
  --percentile "${PERCENTILE}"
  --p_value "${P_VALUE}"
  --n_splits "${N_SPLITS}"
  --n_bins "${N_BINS}"
  --data_dir "${DATA_DIR}"
  --rna_file "${RNA_FILE}"
  --results_dir "${RESULTS_DIR}"
  --protein_list "${PROTEIN_LIST}"
)

if [ "${FULL_OUTPUT}" -eq 1 ]; then
  ARGS+=(--full_output)
fi

echo "Running chunk_id=\${CHUNK_ID}..."
python new_pipeline.py "\${ARGS[@]}"

echo "Done chunk_id=\${CHUNK_ID}"
EOF