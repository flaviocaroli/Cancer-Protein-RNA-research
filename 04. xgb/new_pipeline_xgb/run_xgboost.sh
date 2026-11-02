#!/usr/bin/env bash
set -euo pipefail

# ================= BASIC PATHS =================
DATA_DIR=${DATA_DIR:-$PWD/../data}
SRC_DIR=${SRC_DIR:-$PWD/src}
RESULTS_DIR=${RESULTS_DIR:-$PWD/../data/results_xgb}
LOG_DIR=${LOG_DIR:-$PWD/logs}

# ================= SLURM SETTINGS =================
PARTITION=${PARTITION:-medium_cpu}
TIME=${TIME:-3:00:00}
CPUS_PER_TASK=${CPUS_PER_TASK:-8}
MEM=${MEM:-32G}
ARRAY_LIMIT=${ARRAY_LIMIT:-29}

USE_GPU=${USE_GPU:-0}
GPUS_PER_TASK=${GPUS_PER_TASK:-1}

# ================= PIPELINE DEFAULTS =================
NORMALIZATIONS=${NORMALIZATIONS:-none,zscore,rank,qbin50}

PROTEINS_PER_TASK=${PROTEINS_PER_TASK:-1000}
PROTEIN_START=${PROTEIN_START:-1}
PROTEIN_END=${PROTEIN_END:-}

SPLITS=${SPLITS:-3}
SELECT_K=${SELECT_K:-2000}
EARLY_STOPPING=${EARLY_STOPPING:-10}
TEST_FRACTION=${TEST_FRACTION:-0.2}
RANDOM_STATE=${RANDOM_STATE:-42}
MAX_PARAM_COMBOS=${MAX_PARAM_COMBOS:-8}
N_SAMPLES=${N_SAMPLES:-}

GRID_PRESET=${GRID_PRESET:-fast}
XGB_THREADS=${XGB_THREADS:-$CPUS_PER_TASK}

H5_COMPRESSION_LEVEL=${H5_COMPRESSION_LEVEL:-9}
RESUME_SKIP=${RESUME_SKIP:-1}
OUTPUT_SUFFIX=${OUTPUT_SUFFIX:-}
RUN_ID=${RUN_ID:-$(date +%Y%m%d_%H%M%S)}
CONDA_ENV=${CONDA_ENV:-xgb_loo}

mkdir -p "$RESULTS_DIR" "$LOG_DIR"

# ================= CHECK INPUTS =================
[[ -f "$SRC_DIR/pipeline.py" ]] || { echo "ERROR: missing $SRC_DIR/pipeline.py"; exit 1; }
[[ -f "$SRC_DIR/functions.py" ]] || { echo "ERROR: missing $SRC_DIR/functions.py"; exit 1; }
[[ -f "$DATA_DIR/transcriptomics_filtered.joblib" ]] || { echo "ERROR: missing $DATA_DIR/transcriptomics_filtered.joblib"; exit 1; }
[[ -f "$DATA_DIR/proteomics_matched.joblib" ]] || { echo "ERROR: missing $DATA_DIR/proteomics_matched.joblib"; exit 1; }

case "$GRID_PRESET" in tiny|fast|current) ;; *) echo "ERROR: bad GRID_PRESET"; exit 1 ;; esac
case "$USE_GPU" in 0|1) ;; *) echo "ERROR: USE_GPU must be 0 or 1"; exit 1 ;; esac

# ================= MAKE TASK FILE =================
TASK_FILE="$RESULTS_DIR/xgb_8020_tasks_${RUN_ID}.tsv"

TOTAL_TASKS=$(DATA_DIR="$DATA_DIR" TASK_FILE="$TASK_FILE" \
NORMALIZATIONS="$NORMALIZATIONS" PROTEINS_PER_TASK="$PROTEINS_PER_TASK" \
PROTEIN_START="$PROTEIN_START" PROTEIN_END="$PROTEIN_END" \
python - <<'PY'
import os
import re
from joblib import load

DATA_DIR = os.environ["DATA_DIR"]
TASK_FILE = os.environ["TASK_FILE"]

X = load(os.path.join(DATA_DIR, "transcriptomics_filtered.joblib"))
Y = load(os.path.join(DATA_DIR, "proteomics_matched.joblib"))

x_samples = set(map(str, X.index))
y_index = set(map(str, Y.index))
y_columns = set(map(str, Y.columns))
index_overlap = len(x_samples & y_index)
column_overlap = len(x_samples & y_columns)

if index_overlap >= column_overlap:
    proteins = list(map(str, Y.columns))
    orientation = "samples_as_rows__proteins_as_columns"
else:
    proteins = list(map(str, Y.index))
    orientation = "proteins_as_rows__samples_as_columns"


def check_norm(n):
    n = n.strip().lower()
    if n in {"none", "zscore", "rank"}:
        return n
    if re.fullmatch(r"qbin[0-9]+", n) and int(n.replace("qbin", "")) >= 2:
        return n
    raise SystemExit(f"ERROR: invalid normalization: {n}")

norms = [check_norm(x) for x in os.environ["NORMALIZATIONS"].split(",") if x.strip()]
if not norms:
    raise SystemExit("ERROR: no normalizations provided")

n_proteins = len(proteins)
start = int(os.environ.get("PROTEIN_START") or 1)
end_env = os.environ.get("PROTEIN_END") or ""
end = int(end_env) if end_env else n_proteins
chunk = int(os.environ.get("PROTEINS_PER_TASK") or 1000)

start = max(1, start)
end = min(n_proteins, end)
if start > end:
    raise SystemExit(f"ERROR: PROTEIN_START={start} > PROTEIN_END={end}")
if chunk < 1:
    raise SystemExit("ERROR: PROTEINS_PER_TASK must be >= 1")

chunks = []
s = start
while s <= end:
    e = min(s + chunk - 1, end)
    chunks.append((s, e))
    s = e + 1

with open(TASK_FILE, "w", encoding="utf-8") as f:
    for norm in norms:
        for s, e in chunks:
            f.write(f"{s}\t{e}\t{norm}\n")

print(f"Detected orientation: {orientation}", flush=True)
print(f"Sample overlap with Y.index: {index_overlap}", flush=True)
print(f"Sample overlap with Y.columns: {column_overlap}", flush=True)
print(f"Proteins detected: {n_proteins}", flush=True)
print(f"First proteins: {proteins[:5]}", flush=True)
print(f"Protein range: {start}-{end}", flush=True)
print(f"Proteins per task: {chunk}", flush=True)
print(f"Chunks: {len(chunks)}", flush=True)
print(f"Normalizations: {norms}", flush=True)
print(f"Task file: {TASK_FILE}", flush=True)
print(len(chunks) * len(norms))
PY
)

TOTAL_TASKS=$(echo "$TOTAL_TASKS" | tail -n 1)
if [[ -z "$TOTAL_TASKS" || "$TOTAL_TASKS" -lt 1 ]]; then
  echo "ERROR: no SLURM tasks created."
  exit 1
fi

ARRAY_RANGE="1-${TOTAL_TASKS}"
if [[ -n "$ARRAY_LIMIT" ]]; then
  ARRAY_RANGE="${ARRAY_RANGE}%${ARRAY_LIMIT}"
fi

echo "Submitting array: $ARRAY_RANGE"
echo "Task file: $TASK_FILE"
echo "Results dir: $RESULTS_DIR"
echo "Partition: $PARTITION"
echo "Time: $TIME"

# ================= SUBMIT ARRAY =================
sbatch <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=xgb_8020
#SBATCH --partition=${PARTITION}
#SBATCH --time=${TIME}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --mem=${MEM}
#SBATCH --array=${ARRAY_RANGE}
#SBATCH -o ${LOG_DIR}/%x_%A_%a.out
$([ "$USE_GPU" = "1" ] && echo "#SBATCH --gres=gpu:${GPUS_PER_TASK}")

set -euo pipefail

echo "Job: \$SLURM_JOB_ID  Task: \$SLURM_ARRAY_TASK_ID  Host: \$(hostname)"

export DATA_DIR="${DATA_DIR}"
export RESULTS_DIR="${RESULTS_DIR}"
export OMP_NUM_THREADS=${CPUS_PER_TASK}
export MKL_NUM_THREADS=${CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${CPUS_PER_TASK}
export NUMEXPR_MAX_THREADS=${CPUS_PER_TASK}

if command -v conda >/dev/null 2>&1; then
  set +u
  source "\$(conda info --base)/etc/profile.d/conda.sh" || true
  conda activate "${CONDA_ENV}" || true
  set -u
fi

echo "Python executable: \$(which python)"
python - <<'PY'
import sys
print("Python version:", sys.version)
PY

TASK_LINE=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${TASK_FILE}" || true)
if [[ -z "\$TASK_LINE" ]]; then
  echo "No task line for array index \$SLURM_ARRAY_TASK_ID"
  exit 0
fi

IFS=\$'\t' read -r CHUNK_START CHUNK_END NORM_CURRENT <<< "\$TASK_LINE"

if [[ -z "${OUTPUT_SUFFIX}" ]]; then
  RUN_OUTPUT_SUFFIX="xgb_8020_\${NORM_CURRENT}_${GRID_PRESET}_k${SELECT_K}"
else
  RUN_OUTPUT_SUFFIX="${OUTPUT_SUFFIX}_\${NORM_CURRENT}"
fi

echo "Chunk: \$CHUNK_START-\$CHUNK_END"
echo "Normalization used for model training: \$NORM_CURRENT"
echo "Output suffix: \$RUN_OUTPUT_SUFFIX"

ARGS=(
  --protein_start "\$CHUNK_START"
  --protein_end "\$CHUNK_END"
  --splits ${SPLITS}
  --select_k ${SELECT_K}
  --early_stopping ${EARLY_STOPPING}
  --test_fraction ${TEST_FRACTION}
  --random_state ${RANDOM_STATE}
  --normalization "\$NORM_CURRENT"
  --grid_preset ${GRID_PRESET}
  --xgb_threads ${XGB_THREADS}
  --max_param_combos ${MAX_PARAM_COMBOS}
  --h5_compression_level ${H5_COMPRESSION_LEVEL}
  --metrics_dir "${RESULTS_DIR}"
  --resume_check_dir "${RESULTS_DIR}"
  --output_suffix "\$RUN_OUTPUT_SUFFIX"
)

[[ -n "${N_SAMPLES}" ]] && ARGS+=(--n_samples ${N_SAMPLES})
[[ "${USE_GPU}" = "1" ]] && ARGS+=(--force_gpu)
[[ "${RESUME_SKIP}" = "1" ]] && ARGS+=(--resume_skip)
[[ "${RESUME_SKIP}" = "0" ]] && ARGS+=(--no_resume_skip)

echo "python -u ${SRC_DIR}/pipeline.py \${ARGS[*]}"
python -u "${SRC_DIR}/pipeline.py" "\${ARGS[@]}"

echo "Done."
EOF
