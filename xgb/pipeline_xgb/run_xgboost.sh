#!/bin/bash
set -euo pipefail

# ---------- USER SETTINGS ----------
CHUNK_SIZE=${CHUNK_SIZE:-500}
DATA_DIR=${DATA_DIR:-$HOME/pipeline_proteins/data}
SRC_DIR=${SRC_DIR:-$HOME/pipeline_proteins/src}
PIPELINE_TYPE=${PIPELINE_TYPE:-single}   # aliases: single|sample|protein (mapped in pipeline.py)
MODEL_TYPE=${MODEL_TYPE:-xgboost}
PARTITION=${PARTITION:-medium_gpu}       # adjust to your cluster
GPUS_PER_TASK=${GPUS_PER_TASK:-1}
CPUS_PER_TASK=${CPUS_PER_TASK:-20}
MEM=${MEM:-32G}
TIME=${TIME:-03:00:00}
ENV_ACTIVATE_LINE=${ENV_ACTIVATE_LINE:-"module load miniconda3 && source ~/.bashrc && conda activate test_env"}
FORCE_GPU=${FORCE_GPU:-1}                # 1 -> pass --force_gpu to pipeline.py

# ---------- PREP LIST OF PROTEINS ----------
mkdir -p "$DATA_DIR"
if [[ ! -f "$DATA_DIR/protein_list.txt" ]]; then
  echo "protein_list.txt not found — creating it from proteomics_matched.joblib ..."
  python - <<'PY'
from joblib import load
import os
import sys
data_dir = os.environ.get("DATA_DIR")
Y = load(os.path.join(data_dir, "proteomics_matched.joblib"))
with open(os.path.join(data_dir, "protein_list.txt"), "w") as f:
    for p in Y.columns:
        f.write(str(p) + "\n")
print(f"Wrote {len(Y.columns)} proteins to protein_list.txt")
PY
fi

TOTAL_PROTEINS=${TOTAL_PROTEINS:-$(wc -l < "$DATA_DIR/protein_list.txt")}
echo "TOTAL_PROTEINS=$TOTAL_PROTEINS  CHUNK_SIZE=$CHUNK_SIZE"

# ---------- SUBMIT CHUNKS AS ARRAY JOBS ----------
for START in $(seq 1 $CHUNK_SIZE $TOTAL_PROTEINS); do
  END=$((START + CHUNK_SIZE - 1))
  if [[ $END -gt $TOTAL_PROTEINS ]]; then END=$TOTAL_PROTEINS; fi

  sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=xgb_${START}_${END}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --time=${TIME}
#SBATCH --partition=${PARTITION}
#SBATCH --mem=${MEM}
#SBATCH --array=${START}-${END}
#SBATCH --gpus=${GPUS_PER_TASK}
# If your cluster uses gres instead of --gpus, comment the previous line and use:
# #SBATCH --gres=gpu:${GPUS_PER_TASK}

echo "Job \$SLURM_JOB_ID Task \$SLURM_ARRAY_TASK_ID starting on \$(hostname)"
echo "CUDA devices visible: \$CUDA_VISIBLE_DEVICES"

# ---- ENV / MODULES ----
${ENV_ACTIVATE_LINE}
# Optional: load CUDA explicitly if your site requires it, e.g.:
# module load cuda/12.2

export CUDA_VISIBLE_DEVICES=0

# ---- RUN ----
cd "${SRC_DIR}"

PROTEIN_NAME=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${DATA_DIR}/protein_list.txt")
echo "Protein: \${PROTEIN_NAME}"

if [[ ${FORCE_GPU} -eq 1 ]]; then
  GPU_FLAG="--force_gpu"
else
  GPU_FLAG=""
fi

python pipeline.py \
  --protein "\${PROTEIN_NAME}" \
  --pipeline ${PIPELINE_TYPE} \
  --model ${MODEL_TYPE} \
  --splits 5 \
  --feature_selector pearson \
  --percentile 5.0 \
  --pvalue 0.05 \
  \${GPU_FLAG}

echo "Task \$SLURM_ARRAY_TASK_ID finished."
EOF

  echo "Submitted array ${START}-${END}"
done
