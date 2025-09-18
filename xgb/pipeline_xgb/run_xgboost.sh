#!/bin/bash
set -euo pipefail

# ========= USER SETTINGS =========
CHUNK_SIZE=${CHUNK_SIZE:-500}

DATA_DIR=${DATA_DIR:-$HOME/pipeline_xgb/data}   # must be SHARED & visible to compute nodes
SRC_DIR=${SRC_DIR:-$HOME/pipeline_xgb/src}
LOG_DIR=${LOG_DIR:-$HOME/pipeline_xgb/logs}

PIPELINE_TYPE=${PIPELINE_TYPE:-single}          # aliases -> 'loo'
MODEL_TYPE=${MODEL_TYPE:-xgboost}

PARTITION=${PARTITION:-medium_gpu}              # gpu | medium_gpu | long_gpu | debug_gpu | gpu_plus_8
GPUS_PER_TASK=${GPUS_PER_TASK:-1}
CPUS_PER_TASK=${CPUS_PER_TASK:-20}
MEM=${MEM:-32G}
TIME=${TIME:-3:00:00}

FORCE_GPU=${FORCE_GPU:-1}

# Pipeline knobs (forwarded to pipeline.py)
SPLITS=${SPLITS:-5}
PRESELECT_K=${PRESELECT_K-}                     # leave empty to omit the flag entirely
EARLY_STOPPING=${EARLY_STOPPING:-50}
SELECT_K=${SELECT_K:-0}                         # <=0 disables feature selection entirely
FEATURE_SELECTOR=${FEATURE_SELECTOR:-anova}
MAX_PARAM_COMBOS=${MAX_PARAM_COMBOS:-12}
MINIMAL_SAVE=${MINIMAL_SAVE:-1}
N_SAMPLES=${N_SAMPLES-}                         # empty = all samples
NO_FEATURE_SELECTION=${NO_FEATURE_SELECTION:-1} # 1 = disable FS entirely

# Scheduling / QoS
ARRAY_LIMIT=${ARRAY_LIMIT-}                     # e.g., 2 to limit concurrency; empty = no throttle
RESUME_SKIP=${RESUME_SKIP:-1}                   # 1 = skip tasks whose output already exists

# ========= PREP =========
mkdir -p "$DATA_DIR" "$LOG_DIR"
export DATA_DIR

# Build definitive protein list from joblib (on submit node) and count proteins
if [[ ! -f "$DATA_DIR/proteomics_matched.joblib" ]]; then
  echo "ERROR: $DATA_DIR/proteomics_matched.joblib not found. Set DATA_DIR to your real dataset."
  exit 1
fi

echo "Reading proteins from $DATA_DIR/proteomics_matched.joblib ..."
TOTAL_PROTEINS=$(python - <<'PY'
import os
from joblib import load
p=os.path.join(os.environ["DATA_DIR"],"proteomics_matched.joblib")
Y=load(p)
cols=list(getattr(Y,"columns",[]))
out=os.path.join(os.environ["DATA_DIR"],"protein_list.txt")
with open(out,"w") as f:
    for c in cols:
        f.write(f"{c}\n")
print(len(cols))
PY
)

if [[ -z "${TOTAL_PROTEINS}" || "${TOTAL_PROTEINS}" -lt 1 ]]; then
  echo "ERROR: No proteins found in proteomics_matched.joblib."
  exit 1
fi

# Clamp chunk size to avoid empty tails on tiny datasets
if [[ "${CHUNK_SIZE}" -gt "${TOTAL_PROTEINS}" ]]; then
  CHUNK_SIZE=${TOTAL_PROTEINS}
fi

ARRAY_LIMIT_OPT=""
[[ -n "${ARRAY_LIMIT}" ]] && ARRAY_LIMIT_OPT="%${ARRAY_LIMIT}"

echo "TOTAL_PROTEINS=${TOTAL_PROTEINS}  CHUNK_SIZE=${CHUNK_SIZE}  PARTITION=${PARTITION}  ARRAY_LIMIT=${ARRAY_LIMIT:-none}"

# ========= SUBMIT CHUNKS =========
for START in $(seq 1 "$CHUNK_SIZE" "$TOTAL_PROTEINS"); do
  END=$((START + CHUNK_SIZE - 1))
  [[ $END -gt $TOTAL_PROTEINS ]] && END=$TOTAL_PROTEINS

  sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=xgb_${START}_${END}
#SBATCH --partition=${PARTITION}
#SBATCH --time=${TIME}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --mem=${MEM}
#SBATCH --array=${START}-${END}${ARRAY_LIMIT_OPT}
#SBATCH --gres=gpu:GA100:${GPUS_PER_TASK}
#SBATCH -o ${LOG_DIR}/%x_%A_%a.out

set -euo pipefail
echo "Job \$SLURM_JOB_ID Task \$SLURM_ARRAY_TASK_ID on \$(hostname)"
echo "DATA_DIR=\$DATA_DIR"
echo "CUDA_VISIBLE_DEVICES=\$CUDA_VISIBLE_DEVICES"

# ---- ENV / MODULES (bootstrap conda) ----
set +u
module load miniconda3
source "\$(conda info --base)/etc/profile.d/conda.sh"
# Create env once per user (file lock prevents stampede across array tasks)
LOCKDIR="\$HOME/.env_locks/test_env_install.lock"
mkdir -p "\$HOME/.env_locks"
if ! conda env list | awk '{print \$1}' | grep -qx "test_env"; then
  if mkdir "\$LOCKDIR" 2>/dev/null; then
    trap 'rmdir "\$LOCKDIR"' RETURN
    conda create -y -n test_env python=3.10
    conda activate test_env
    conda install -y -c conda-forge joblib pandas numpy scipy scikit-learn "xgboost>=2.0" "openpyxl>=3.1"
  else
    echo "Waiting for shared env install lock..."
    while [[ -d "\$LOCKDIR" ]]; do sleep 8; done
    conda activate test_env
  fi
else
  conda activate test_env
fi
set -u

# ---- Threading hygiene (avoid oversubscription) ----
export OMP_NUM_THREADS=${CPUS_PER_TASK}
export MKL_NUM_THREADS=${CPUS_PER_TASK}
export NUMEXPR_MAX_THREADS=${CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${CPUS_PER_TASK}

# ---- Node-local scratch with safe fallbacks ----
# Prefer site scratch, else SLURM_TMPDIR, else TMPDIR, else /tmp
if [[ -d "/scratch/\$USER" ]]; then
  NODE_SCRATCH="/scratch/\$USER/\${SLURM_JOB_ID:-\$\$}"
  mkdir -p "\$NODE_SCRATCH" || true
else
  NODE_SCRATCH="\${SLURM_TMPDIR:-\${TMPDIR:-/tmp}}"
fi
WORKDIR="\${NODE_SCRATCH}/xgb_run_\${SLURM_ARRAY_JOB_ID:-\${SLURM_JOB_ID:-\$\$}}"
mkdir -p "\$WORKDIR/src" "\$WORKDIR/data"
echo "Using WORKDIR=\$WORKDIR"
trap 'rm -rf "\$WORKDIR"' EXIT

# ---- Copy only what's needed to local scratch ----
rsync -a --delete "${SRC_DIR}/" "\$WORKDIR/src/"
rsync -a "${DATA_DIR}/transcriptomics_filtered.joblib" "\$WORKDIR/data/"
rsync -a "${DATA_DIR}/proteomics_matched.joblib"     "\$WORKDIR/data/"

cd "\$WORKDIR/src"

# Resolve protein name by array index using the prewritten list (no joblib load here)
if [[ ! -f "\$DATA_DIR/protein_list.txt" ]]; then
  echo "ERROR: protein_list.txt not found in \$DATA_DIR"
  exit 1
fi
PROTEIN_NAME=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "\$DATA_DIR/protein_list.txt" || true)
if [[ -z "\$PROTEIN_NAME" ]]; then
  echo "Task \${SLURM_ARRAY_TASK_ID} exceeds list length; exiting."
  exit 0
fi
echo "Protein: \${PROTEIN_NAME}"

# Optional: skip if output already exists on shared FS (resume-safe)
if [[ ${RESUME_SKIP} -eq 1 ]]; then
  OUT_SHARED="${DATA_DIR}/\${PROTEIN_NAME}_loo_xgboost_LOO_results.joblib"
  if [[ -f "\$OUT_SHARED" ]]; then
    echo "Output exists, skipping: \$OUT_SHARED"
    exit 0
  fi
fi

# ---- Build CLI args ----
ARGS=()
ARGS+=(--protein "\${PROTEIN_NAME}")
ARGS+=(--pipeline ${PIPELINE_TYPE})
ARGS+=(--model ${MODEL_TYPE})
ARGS+=(--splits ${SPLITS})
ARGS+=(--feature_selector ${FEATURE_SELECTOR})
ARGS+=(--select_k ${SELECT_K})
ARGS+=(--early_stopping ${EARLY_STOPPING})
ARGS+=(--max_param_combos ${MAX_PARAM_COMBOS})
[[ ${FORCE_GPU} -eq 1 ]]            && ARGS+=(--force_gpu)
[[ ${MINIMAL_SAVE} -eq 1 ]]         && ARGS+=(--minimal_save)
[[ -n "${N_SAMPLES}" ]]             && ARGS+=(--n_samples ${N_SAMPLES})
[[ ${NO_FEATURE_SELECTION} -eq 1 ]] && ARGS+=(--no_feature_selection)
# Only pass preselect_k if it's numeric AND FS is enabled
if [[ -n "${PRESELECT_K}" && ${NO_FEATURE_SELECTION} -eq 0 && ${SELECT_K} -gt 0 ]]; then
  if [[ "${PRESELECT_K}" =~ ^[0-9]+$ ]]; then
    ARGS+=(--preselect_k ${PRESELECT_K})
  else
    echo "WARN: PRESELECT_K='${PRESELECT_K}' is not numeric; ignoring."
  fi
fi

# ---- Run (pipeline.py reads ../data relative to src) ----
python pipeline.py "\${ARGS[@]}"

# ---- Copy result back to shared DATA_DIR ----
OUT_LOCAL="\$WORKDIR/data/\${PROTEIN_NAME}_loo_xgboost_LOO_results.joblib"
if [[ -f "\$OUT_LOCAL" ]]; then
  rsync -a "\$OUT_LOCAL" "${DATA_DIR}/"
  echo "Saved -> ${DATA_DIR}/\$(basename "\$OUT_LOCAL")"
else
  echo "WARNING: expected output not found: \$OUT_LOCAL"
fi

echo "Task \${SLURM_ARRAY_TASK_ID} finished."
EOF

  echo "Submitted array ${START}-${END}"
done
