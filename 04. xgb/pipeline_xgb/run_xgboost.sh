#!/usr/bin/env bash
set -euo pipefail

# ================= USER SETTINGS (override via env) =================
# Paths (must be visible on all nodes)
DATA_DIR=${DATA_DIR:-$PWD/data}
SRC_DIR=${SRC_DIR:-$PWD/src}
LOG_DIR=${LOG_DIR:-$PWD/logs}

# SLURM resources
# Cluster partitions seen in sinfo:
#   CPU: debug_cpu, medium_cpu, compute, defq
#   GPU: debug_gpuh200, medium_gpuh200, gpuh200, long_gpuh200,
#        debug_gpunew, medium_gpunew, gpunew, long_gpunew
# Default is CPU because this pipeline runs many small/medium XGBoost fits.
PARTITION=${PARTITION:-medium_cpu}
TIME=${TIME:-03:00:00}
CPUS_PER_TASK=${CPUS_PER_TASK:-8}
MEM=${MEM:-32G}
USE_GPU=${USE_GPU:-0}            # 1 = request GPU + pass --force_gpu; 0 = CPU only
GPUS_PER_TASK=${GPUS_PER_TASK:-1}

# Throttle concurrent array tasks.
# You said you can send at most 29 jobs, so default to 29 concurrent protein tasks.
ARRAY_LIMIT=${ARRAY_LIMIT:-29}   # e.g. 29 -> "--array=1-N%29"; empty = no throttle

# Pipeline knobs forwarded to pipeline.py
SPLITS=${SPLITS:-3}
SELECT_K=${SELECT_K:-500}        # 0 = no Spearman prefilter; 500/1000 recommended for speed
EARLY_STOPPING=${EARLY_STOPPING:-30}
MAX_PARAM_COMBOS=${MAX_PARAM_COMBOS:-}  # empty = use whole selected grid preset
MINIMAL_SAVE=${MINIMAL_SAVE:-1}  # 1 = pass --minimal_save
N_SAMPLES=${N_SAMPLES:-}         # leave empty to use all samples
NO_FEATURE_SELECTION=${NO_FEATURE_SELECTION:-0}  # 1 = pass --no_feature_selection
RESUME_SKIP=${RESUME_SKIP:-1}    # 1 = skip proteins with existing outputs

# New pipeline knobs
NORMALIZATION=${NORMALIZATION:-none}       # none | zscore | rank | qbin5 | qbin10
GRID_PRESET=${GRID_PRESET:-fast}           # tiny | fast | current
TUNING_STRATEGY=${TUNING_STRATEGY:-nested} # nested = strict LOO; once = faster but less strict
XGB_THREADS=${XGB_THREADS:-${CPUS_PER_TASK}}

# Output naming
# If OUTPUT_SUFFIX is empty, build a suffix from settings to avoid overwriting comparisons.
OUTPUT_SUFFIX=${OUTPUT_SUFFIX:-}
if [[ -z "${OUTPUT_SUFFIX}" ]]; then
  OUTPUT_SUFFIX="loo_xgboost_${NORMALIZATION}_${GRID_PRESET}_${TUNING_STRATEGY}_k${SELECT_K}"
fi
# ====================================================================

mkdir -p "$DATA_DIR" "$LOG_DIR"
export DATA_DIR SRC_DIR LOG_DIR

# --- Ensure required inputs exist ---
if [[ ! -f "$DATA_DIR/proteomics_matched.joblib" || ! -f "$DATA_DIR/transcriptomics_filtered.joblib" ]]; then
  echo "ERROR: Missing data files in $DATA_DIR"
  echo "Expected:"
  echo " - $DATA_DIR/transcriptomics_filtered.joblib"
  echo " - $DATA_DIR/proteomics_matched.joblib"
  exit 1
fi

# --- Basic validation of options before submitting ---
case "$NORMALIZATION" in
  none|zscore|rank|qbin5|qbin10) ;;
  *) echo "ERROR: NORMALIZATION must be one of: none, zscore, rank, qbin5, qbin10"; exit 1 ;;
esac

case "$GRID_PRESET" in
  tiny|fast|current) ;;
  *) echo "ERROR: GRID_PRESET must be one of: tiny, fast, current"; exit 1 ;;
esac

case "$TUNING_STRATEGY" in
  nested|once) ;;
  *) echo "ERROR: TUNING_STRATEGY must be one of: nested, once"; exit 1 ;;
esac

if [[ "$USE_GPU" != "0" && "$USE_GPU" != "1" ]]; then
  echo "ERROR: USE_GPU must be 0 or 1"
  exit 1
fi

if [[ "$USE_GPU" = "1" ]]; then
  case "$PARTITION" in
    *gpu*|personae*|lcgroup*) ;;
    *)
      echo "WARNING: USE_GPU=1 but PARTITION='$PARTITION' does not look like a GPU partition."
      echo "         Recommended GPU partitions: medium_gpuh200, gpuh200, long_gpuh200, medium_gpunew."
      ;;
  esac
else
  case "$PARTITION" in
    *gpu*|personae*|lcgroup*)
      echo "WARNING: USE_GPU=0 but PARTITION='$PARTITION' looks like a GPU partition."
      echo "         Recommended CPU partitions: medium_cpu, compute, defq, debug_cpu."
      ;;
  esac
fi

# --- Build protein list + count ---
echo "Reading proteins from $DATA_DIR/proteomics_matched.joblib ..."
TOTAL_PROTEINS=$(python - <<'PY'
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
if [[ -n "${ARRAY_LIMIT}" ]]; then
  ARRAY_RANGE="${ARRAY_RANGE}%${ARRAY_LIMIT}"
fi

echo "Submitting SLURM array: ${ARRAY_RANGE} (TOTAL_PROTEINS=${TOTAL_PROTEINS})"
echo "Settings: PARTITION=${PARTITION}, TIME=${TIME}, CPUS_PER_TASK=${CPUS_PER_TASK}, MEM=${MEM}"
echo "Pipeline: USE_GPU=${USE_GPU}, SELECT_K=${SELECT_K}, NORMALIZATION=${NORMALIZATION}, GRID_PRESET=${GRID_PRESET}, TUNING_STRATEGY=${TUNING_STRATEGY}, OUTPUT_SUFFIX=${OUTPUT_SUFFIX}"

# --- Submit ONE array job: 1 task per protein ---
sbatch <<EOF2
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
echo "SRC_DIR=${SRC_DIR}"
echo "LOG_DIR=${LOG_DIR}"
echo "USE_GPU=${USE_GPU}"
echo "SELECT_K=${SELECT_K}"
echo "NORMALIZATION=${NORMALIZATION}"
echo "GRID_PRESET=${GRID_PRESET}"
echo "TUNING_STRATEGY=${TUNING_STRATEGY}"
echo "OUTPUT_SUFFIX=${OUTPUT_SUFFIX}"

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

# Node-local scratch: prefer site scratch, then SLURM_TMPDIR, else /tmp
if [[ -d "/scratch/\$USER" ]]; then
  NODE_SCRATCH="/scratch/\$USER/\${SLURM_JOB_ID:-\$\$}"
  mkdir -p "\$NODE_SCRATCH" || true
else
  NODE_SCRATCH="\${SLURM_TMPDIR:-/tmp}"
fi
WORKDIR="\$NODE_SCRATCH/xgb_\${SLURM_JOB_ID:-\$\$}_\${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "\$WORKDIR/src" "\$WORKDIR/data"
trap 'rm -rf "\$WORKDIR"' EXIT

# Stage code + data.
# pipeline.py supports DATA_DIR, but this preserves the ../data layout too.
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
OUT_SHARED="${DATA_DIR}/\${PROTEIN_NAME}_${OUTPUT_SUFFIX}.joblib"
if [[ ${RESUME_SKIP} -eq 1 && -f "\$OUT_SHARED" ]]; then
  echo "Exists, skipping: \$OUT_SHARED"
  exit 0
fi

cd "\$WORKDIR/src"

# Build args for updated pipeline.py
ARGS=( --protein "\$PROTEIN_NAME"
       --splits ${SPLITS}
       --select_k ${SELECT_K}
       --early_stopping ${EARLY_STOPPING}
       --normalization ${NORMALIZATION}
       --grid_preset ${GRID_PRESET}
       --tuning_strategy ${TUNING_STRATEGY}
       --xgb_threads ${XGB_THREADS}
       --output_suffix ${OUTPUT_SUFFIX} )

# Optional flags
[[ -n "${MAX_PARAM_COMBOS}" ]]      && ARGS+=( --max_param_combos ${MAX_PARAM_COMBOS} )
[[ -n "${N_SAMPLES}" ]]             && ARGS+=( --n_samples ${N_SAMPLES} )
[[ ${MINIMAL_SAVE} -eq 1 ]]         && ARGS+=( --minimal_save )
[[ ${NO_FEATURE_SELECTION} -eq 1 ]] && ARGS+=( --no_feature_selection )
[[ ${USE_GPU} -eq 1 ]]              && ARGS+=( --force_gpu )

echo "python pipeline.py \${ARGS[*]}"
python pipeline.py "\${ARGS[@]}"

# Copy result back to shared data dir
OUT_LOCAL="\$WORKDIR/data/\${PROTEIN_NAME}_${OUTPUT_SUFFIX}.joblib"
if [[ -f "\$OUT_LOCAL" ]]; then
  rsync -a "\$OUT_LOCAL" "${DATA_DIR}/"
  echo "Saved -> ${DATA_DIR}/\$(basename "\$OUT_LOCAL")"
else
  echo "WARNING: output not found: \$OUT_LOCAL"
  echo "Files in work data dir:"
  ls -lh "\$WORKDIR/data" || true
fi

echo "Task \${SLURM_ARRAY_TASK_ID} done."
EOF2