#!/usr/bin/env bash
set -euo pipefail

# ================= USER SETTINGS (override via env) =================
# Paths (must be visible on all nodes)
DATA_DIR=${DATA_DIR:-$PWD/data}
SRC_DIR=${SRC_DIR:-$PWD/src}
LOG_DIR=${LOG_DIR:-$PWD/logs}

# SLURM resources
PARTITION=${PARTITION:-medium_gpu}
TIME=${TIME:-03:00:00}
CPUS_PER_TASK=${CPUS_PER_TASK:-8}
MEM=${MEM:-32G}
USE_GPU=${USE_GPU:-1}            # 1 = request GPU + pass --force_gpu; 0 = CPU only
GPUS_PER_TASK=${GPUS_PER_TASK:-1}

# Throttle concurrent array tasks (optional)
ARRAY_LIMIT=${ARRAY_LIMIT:-}     # e.g. 16 → "--array=1-N%16"

# Pipeline knobs (forwarded to pipeline.py)
SPLITS=${SPLITS:-3}
SELECT_K=${SELECT_K:-0}        # 0 = no Spearman prefilter (keep all features)
EARLY_STOPPING=${EARLY_STOPPING:-100}
MAX_PARAM_COMBOS=${MAX_PARAM_COMBOS:-12}
MINIMAL_SAVE=${MINIMAL_SAVE:-1}  # 1 = pass --minimal_save
N_SAMPLES=${N_SAMPLES:-}         # leave empty to use all samples
NO_FEATURE_SELECTION=${NO_FEATURE_SELECTION:-0}  # 1 = pass --no_feature_selection
RESUME_SKIP=${RESUME_SKIP:-1}    # 1 = skip proteins with existing outputs
# ====================================================================

mkdir -p "$DATA_DIR" "$LOG_DIR"
# Ensure child processes (Python) see these
export DATA_DIR SRC_DIR LOG_DIR

# --- Ensure required inputs exist ---
if [[ ! -f "$DATA_DIR/proteomics_matched.joblib" || ! -f "$DATA_DIR/transcriptomics_filtered.joblib" ]]; then
  echo "ERROR: Missing data files in $DATA_DIR"
  echo "Expected:"
  echo " - $DATA_DIR/transcriptomics_filtered.joblib"
  echo " - $DATA_DIR/proteomics_matched.joblib"
  exit 1
fi

# --- Build protein list + count (safe to missing env via getenv) ---
echo "Reading proteins from $DATA_DIR/proteomics_matched.joblib ..."
TOTAL_PROTEINS=$(
python - <<'PY'
import os
from joblib import load
data_dir = os.getenv("DATA_DIR", ".")
p = os.path.join(data_dir, "proteomics_matched.joblib")
Y = load(p)
cols = list(getattr(Y, "columns", []))
out = os.path.join(data_dir, "protein_list.txt")
with open(out, "w", encoding="utf-8") as f:
    for c in cols:
        f.write(f"{c}\n")
print(len(cols))
PY
)

if [[ -z "${TOTAL_PROTEINS}" || "${TOTAL_PROTEINS}" -lt 1 ]]; then
  echo "ERROR: No proteins found in proteomics_matched.joblib."
  exit 1
fi

ARRAY_RANGE="1-${TOTAL_PROTEINS}"
[[ -n "${ARRAY_LIMIT}" ]] && ARRAY_RANGE="${ARRAY_RANGE}%${ARRAY_LIMIT}"
echo "Submitting SLURM array: ${ARRAY_RANGE} (TOTAL_PROTEINS=${TOTAL_PROTEINS})"

# --- Submit ONE array job (1 task per protein) ---
sbatch <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=xgb_loo
#SBATCH --partition=${PARTITION}
#SBATCH --time=${TIME}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --mem=${MEM}
#SBATCH --array=${ARRAY_RANGE}
#SBATCH -o ${LOG_DIR}/%x_%A_%a.out
$([ "${USE_GPU}" = "1" ] && echo "#SBATCH --gres=gpu:${GPUS_PER_TASK}")

set -euo pipefail
echo "JobID=\$SLURM_JOB_ID  Task=\$SLURM_ARRAY_TASK_ID  Host=\$(hostname)"
echo "DATA_DIR=${DATA_DIR}"

# ---------- Optional: activate your Python env ----------
# set +u
# module load miniconda3 || true
# source "\$(conda info --base)/etc/profile.d/conda.sh" || true
# conda activate xgb || true
# set -u

# Threading hygiene
export OMP_NUM_THREADS=${CPUS_PER_TASK}
export MKL_NUM_THREADS=${CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${CPUS_PER_TASK}
export NUMEXPR_MAX_THREADS=${CPUS_PER_TASK}

# Node-local scratch (prefer site scratch, then SLURM_TMPDIR, else /tmp)
if [[ -d "/scratch/\$USER" ]]; then
  NODE_SCRATCH="/scratch/\$USER/\${SLURM_JOB_ID:-\$\$}"
  mkdir -p "\$NODE_SCRATCH" || true
else
  NODE_SCRATCH="\${SLURM_TMPDIR:-/tmp}"
fi
WORKDIR="\$NODE_SCRATCH/xgb_\${SLURM_JOB_ID:-\$\$}"
mkdir -p "\$WORKDIR/src" "\$WORKDIR/data"
trap 'rm -rf "\$WORKDIR"' EXIT

# Stage code + data so pipeline.py sees ../data relative to src
rsync -a --delete "${SRC_DIR}/" "\$WORKDIR/src/"
rsync -a "${DATA_DIR}/transcriptomics_filtered.joblib" "\$WORKDIR/data/"
rsync -a "${DATA_DIR}/proteomics_matched.joblib"     "\$WORKDIR/data/"

# Resolve protein by array index
if [[ ! -f "${DATA_DIR}/protein_list.txt" ]]; then
  echo "ERROR: protein_list.txt not found in ${DATA_DIR}"
  exit 1
fi
PROTEIN_NAME=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${DATA_DIR}/protein_list.txt" || true)
if [[ -z "\$PROTEIN_NAME" ]]; then
  echo "No protein for index \$SLURM_ARRAY_TASK_ID; exiting."
  exit 0
fi
echo "Protein: \$PROTEIN_NAME"

# Resume-safe: skip if output already exists on shared FS
OUT_SHARED="${DATA_DIR}/\${PROTEIN_NAME}_loo_xgboost_LOO_results.joblib"
if [[ ${RESUME_SKIP} -eq 1 && -f "\$OUT_SHARED" ]]; then
  echo "Exists, skipping: \$OUT_SHARED"
  exit 0
fi

cd "\$WORKDIR/src"

# Build args for pipeline.py (new simplified CLI)
ARGS=( --protein "\$PROTEIN_NAME"
       --splits ${SPLITS}
       --select_k ${SELECT_K}
       --early_stopping ${EARLY_STOPPING} )

# Optional flags
[[ -n "${MAX_PARAM_COMBOS}" ]]   && ARGS+=( --max_param_combos ${MAX_PARAM_COMBOS} )
[[ -n "${N_SAMPLES}" ]]          && ARGS+=( --n_samples ${N_SAMPLES} )
[[ ${MINIMAL_SAVE} -eq 1 ]]      && ARGS+=( --minimal_save )
[[ ${NO_FEATURE_SELECTION} -eq 1 ]] && ARGS+=( --no_feature_selection )
[[ ${USE_GPU} -eq 1 ]]           && ARGS+=( --force_gpu )

echo "python pipeline.py \${ARGS[*]}"
python pipeline.py "\${ARGS[@]}"

# Copy result back to shared data dir
OUT_LOCAL="\$WORKDIR/data/\${PROTEIN_NAME}_loo_xgboost_LOO_results.joblib"
if [[ -f "\$OUT_LOCAL" ]]; then
  rsync -a "\$OUT_LOCAL" "${DATA_DIR}/"
  echo "Saved -> ${DATA_DIR}/\$(basename "\$OUT_LOCAL")"
else
  echo "WARNING: output not found: \$OUT_LOCAL"
fi

echo "Task \${SLURM_ARRAY_TASK_ID} done."
EOF
