#!/bin/bash
# new_run_corr.sh

set -euo pipefail

# ✅ Capture current directory to avoid Slurm spool errors
BASE_DIR="$(pwd)"
DATA_DIR="${BASE_DIR}/data"
SRC_DIR="${BASE_DIR}"

PROTEIN_LIST="${DATA_DIR}/protein_list.txt"
RESULTS_DIR="${DATA_DIR}/results"
LOG_DIR="${DATA_DIR}/logs"

# Ensure directories exist before submission
mkdir -p "${RESULTS_DIR}" "${LOG_DIR}"

# ---------------------------
# Strict Chunking Control
# ---------------------------
TARGET_CHUNKS=25

# ---------------------------
# Correlation Knobs
# ---------------------------
P_THRESHOLD=0.05
BOOTSTRAP=1
N_BOOTSTRAPS=100
PER_PROTEIN_FORMAT="parquet"

# Validate input file
if [ ! -f "${PROTEIN_LIST}" ]; then
  echo "Error: protein_list.txt not found at ${PROTEIN_LIST}"
  exit 1
fi

# Count proteins (ignoring comments and empty lines)
TOTAL_PROTEINS=$(grep -v '^#' "${PROTEIN_LIST}" | sed '/^$/d' | wc -l)
if [ "${TOTAL_PROTEINS}" -eq 0 ]; then
  echo "protein_list.txt is empty."
  exit 0
fi

# Calculate GROUP_SIZE to hit exactly TARGET_CHUNKS
# Ceiling math: (Total + Chunks - 1) / Chunks
GROUP_SIZE=$(( (TOTAL_PROTEINS + TARGET_CHUNKS - 1) / TARGET_CHUNKS ))
[ "${GROUP_SIZE}" -lt 1 ] && GROUP_SIZE=1

# Determine the actual number of chunks (last one might be smaller)
N_CHUNKS=$(( (TOTAL_PROTEINS + GROUP_SIZE - 1) / GROUP_SIZE ))
LAST_CHUNK=$(( N_CHUNKS - 1 ))

echo "Total Proteins: ${TOTAL_PROTEINS}"
echo "Group Size: ${GROUP_SIZE}"
echo "Submitting array: 0-${LAST_CHUNK} (Total chunks: ${N_CHUNKS})"

# ---------------------------
# Submit the Job
# ---------------------------
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=corr_chunk
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --array=0-${LAST_CHUNK}
#SBATCH --output=${LOG_DIR}/corr_%A_%a.out
#SBATCH --error=${LOG_DIR}/corr_%A_%a.err
#SBATCH --chdir=${SRC_DIR}

set -euo pipefail

module load miniconda3
source activate base
conda activate test_env

CHUNK_ID=\${SLURM_ARRAY_TASK_ID}

# Construct arguments
ARGS=(
  --chunk_id "\${CHUNK_ID}"
  --group_size "${GROUP_SIZE}"
  --p_threshold "${P_THRESHOLD}"
  --n_bootstraps "${N_BOOTSTRAPS}"
  --data_dir "${DATA_DIR}"
  --results_dir "${RESULTS_DIR}"
  --protein_list "${PROTEIN_LIST}"
  --per_protein_format "${PER_PROTEIN_FORMAT}"
)

if [ "${BOOTSTRAP}" -eq 0 ]; then
  ARGS+=(--no_bootstrap)
fi

echo "Running chunk_id=\${CHUNK_ID} in \$(pwd)"
python corr_pipeline.py "\${ARGS[@]}"
echo "Done chunk_id=\${CHUNK_ID}"
EOF