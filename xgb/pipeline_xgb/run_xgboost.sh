#!/bin/bash
set -euo pipefail

# ========= USER SETTINGS =========
CHUNK_SIZE=${CHUNK_SIZE:-500}

DATA_DIR=${DATA_DIR:-$HOME/pipeline_xgb/data}
SRC_DIR=${SRC_DIR:-$HOME/pipeline_xgb/src}
LOG_DIR=${LOG_DIR:-$HOME/pipeline_xgb/logs}

PIPELINE_TYPE=${PIPELINE_TYPE:-single}   # aliases: single|sample|protein
MODEL_TYPE=${MODEL_TYPE:-xgboost}

# If 3h isn't enough on your cluster, default to long_gpu. Override at run time if you like.
PARTITION=${PARTITION:-long_gpu}          # gpu | medium_gpu | long_gpu | debug_gpu | gpu_plus_8

# Keep within partition limits (long_gpu: up to 3-00:10:00 on your cluster)
GPUS_PER_TASK=${GPUS_PER_TASK:-1}
CPUS_PER_TASK=${CPUS_PER_TASK:-20}
MEM=${MEM:-32G}
TIME=${TIME:-12:00:00}                    # bumpable; must be <= partition limit

# XGBoost GPU only (1) vs allow CPU fallback (0)
FORCE_GPU=${FORCE_GPU:-1}

# Lighter pipeline knobs (forwarded to pipeline.py)
SPLITS=${SPLITS:-5}                       # fewer inner folds speeds things up
PRESELECT_K=${PRESELECT_K:-5000}          # fast corr prefilter before p-values
EARLY_STOPPING=${EARLY_STOPPING:-50}      # early stopping rounds in inner CV

# Throttle concurrent tasks per array (helps QOSMaxGRESPerUser)
ARRAY_LIMIT=${ARRAY_LIMIT:-}              # empty = no throttle; set 1 or 2 to cap

# If protein_list.txt exists we’ll use it; otherwise submit a safe upper bound
MAX_PROTEINS_GUESS=${MAX_PROTEINS_GUESS:-6000}

# ========= PREP =========
mkdir -p "$DATA_DIR" "$LOG_DIR"
export DATA_DIR

if [[ -f "$DATA_DIR/protein_list.txt" ]]; then
  TOTAL_PROTEINS=$(wc -l < "$DATA_DIR/protein_list.txt" | tr -d ' ')
else
  echo "No protein_list.txt found; using MAX_PROTEINS_GUESS=${MAX_PROTEINS_GUESS}."
  TOTAL_PROTEINS=${MAX_PROTEINS_GUESS}
fi
echo "TOTAL_PROTEINS=${TOTAL_PROTEINS}  CHUNK_SIZE=${CHUNK_SIZE}  PARTITION=${PARTITION}  ARRAY_LIMIT=${ARRAY_LIMIT:-none}"

ARRAY_LIMIT_OPT=""
[[ -n "${ARRAY_LIMIT}" ]] && ARRAY_LIMIT_OPT="%${ARRAY_LIMIT}"

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
# Request a full GA100 GPU (matches your sinfo GRES)
#SBATCH --gres=gpu:GA100:${GPUS_PER_TASK}
#SBATCH -o ${LOG_DIR}/%x_%A_%a.out

set -euo pipefail
echo "Job \$SLURM_JOB_ID Task \$SLURM_ARRAY_TASK_ID on \$(hostname)"
echo "Partition: ${PARTITION}  Time: ${TIME}  CPUs: ${CPUS_PER_TASK}  Mem: ${MEM}  GPUs: ${GPUS_PER_TASK}"
echo "CUDA_VISIBLE_DEVICES=\$CUDA_VISIBLE_DEVICES"

# ---- ENV / MODULES (bootstrap conda + deps) ----
set +u
module load miniconda3
source "\$(conda info --base)/etc/profile.d/conda.sh"

# create env if missing
if ! conda env list | awk '{print \$1}' | grep -qx "test_env"; then
  conda create -y -n test_env python=3.10
fi
conda activate test_env
set -u

# ---- REQUIRED DEPS ----
DEPS=(joblib pandas numpy scipy scikit-learn "xgboost>=2.0" "openpyxl>=3.1")
mkdir -p "\$HOME/.env_locks"
LOCKDIR="\$HOME/.env_locks/test_env_install.lock"

MISSING=\$(python - <<'PY'
import pkgutil
need=["joblib","pandas","numpy","scipy","sklearn","xgboost","openpyxl"]
missing=[m for m in need if pkgutil.find_loader(m) is None]
print(" ".join(missing))
PY
)

if [[ -n "\$MISSING" ]]; then
  if mkdir "\$LOCKDIR" 2>/dev/null; then
    trap 'rmdir "\$LOCKDIR"' EXIT
    echo "Installing missing packages: \$MISSING"
    conda install -y -c conda-forge "\${DEPS[@]}" || {
      echo "Conda install failed or partial; falling back to pip for remaining..."
      python -m pip install -U pip
      python - <<'PY'
import sys, subprocess, pkgutil
deps=["joblib","pandas","numpy","scipy","scikit-learn","xgboost>=2.0","openpyxl>=3.1"]
for d in deps:
    name=d.split("==")[0].split(">=")[0]
    if pkgutil.find_loader(name) is None:
        subprocess.check_call([sys.executable,"-m","pip","install","-U","--no-cache-dir",d])
PY
    }
  else
    echo "Waiting for another task to finish env install..."
    while [[ -d "\$LOCKDIR" ]]; do sleep 10; done
  fi
fi

which python
python -V
python - <<'PY'
import joblib, pandas as pd, numpy as np, sklearn, xgboost as xgb, openpyxl
print("joblib", joblib.__version__,
      "| pandas", pd.__version__,
      "| numpy", np.__version__,
      "| sklearn", sklearn.__version__,
      "| xgboost", xgb.__version__,
      "| openpyxl", openpyxl.__version__)
PY

# ---- RUN ----
cd "${SRC_DIR}"

# Determine protein name by index directly from the joblib (no prebuilt list needed)
PROTEIN_NAME=\$(python - <<'PY'
import os, sys
from joblib import load
data_dir = os.environ["DATA_DIR"]
task_id = int(os.environ["SLURM_ARRAY_TASK_ID"])
Y = load(os.path.join(data_dir, "proteomics_matched.joblib"))
cols = list(Y.columns)
n = len(cols)
if task_id > n:
    print("")  # signal to skip
else:
    print(cols[task_id-1])
PY
)

if [[ -z "\$PROTEIN_NAME" ]]; then
  echo "Task \${SLURM_ARRAY_TASK_ID} exceeds number of proteins; exiting."
  exit 0
fi

echo "Protein: \${PROTEIN_NAME}"

GPU_FLAG=""
[[ ${FORCE_GPU} -eq 1 ]] && GPU_FLAG="--force_gpu"

# Forward the "lighter" flags to pipeline.py
python pipeline.py \
  --protein "\${PROTEIN_NAME}" \
  --pipeline ${PIPELINE_TYPE} \
  --model ${MODEL_TYPE} \
  --splits ${SPLITS} \
  --feature_selector pearson \
  --percentile 5.0 \
  --pvalue 0.05 \
  --preselect_k ${PRESELECT_K} \
  --early_stopping ${EARLY_STOPPING} \
  \${GPU_FLAG}

echo "Task \${SLURM_ARRAY_TASK_ID} finished."
EOF

  echo "Submitted array ${START}-${END}"
done
