#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${PROJECT_DIR}"

# This creates a tiny synthetic dataset and runs a local 2-split test.
python - <<'PY'
from pathlib import Path
import numpy as np
import pandas as pd
import joblib
rng = np.random.default_rng(1)
root = Path('data')
root.mkdir(exist_ok=True)
samples = [f'S{i:03d}' for i in range(60)]
genes = [f'G{i:04d}' for i in range(300)] + ['P000', 'P001', 'P002', 'P003', 'P004']
proteins = [f'P{i:03d}' for i in range(10)]
X = pd.DataFrame(rng.normal(size=(len(samples), len(genes))), index=samples, columns=genes)
Y = pd.DataFrame(rng.normal(size=(len(samples), len(proteins))), index=samples, columns=proteins)
# Make a few proteins predictable by same-gene columns.
for p in proteins[:5]:
    if p in X.columns:
        Y[p] = 1.2 * X[p] + rng.normal(scale=0.3, size=len(samples))
joblib.dump(X, root / 'transcriptomics_filtered.joblib')
joblib.dump(Y, root / 'proteomics_matched.joblib')
print('Synthetic data written:', X.shape, Y.shape)
PY

rm -rf results
SUBMIT=local \
MODELS=proxy,enet \
NORMALIZATIONS=zscore \
SELECT_KS=50 \
N_REPEATS=1 \
CV_FOLDS=3 \
GRID_PRESET=tiny \
MAX_PARAM_COMBOS=2 \
PROTEINS_PER_TASK=10 \
SPLITS_PER_TASK=1 \
KEEP_PREDICTIONS=1 \
bash run_repeated.sh

PROXY_CONFIG=proxy:zscore:0 \
ENET_CONFIG=enet:zscore:50 \
N_BOOT=100 \
bash run_after_jobs.sh

echo "Smoke test finished. Check results/cascade/robust_cascade_atlas.tsv"
