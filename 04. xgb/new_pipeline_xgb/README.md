# new_pipeline_xgb

80/20 XGBoost pipeline for protein prediction.

Expected layout:

```text
new_pipeline_xgb/
  data/
    transcriptomics_filtered.joblib
    proteomics_matched.joblib
  src/
    functions.py
    pipeline.py
  run_xgboost.sh
  results_xgb/
  logs/
```

Each saved H5 contains the requested minimal fields only:

- protein
- train_samples
- test_samples
- test_mse
- test_pearson_corr
- test_cosine_similarity
- best_params
- selected_features
- final model trained on all Dataset 1 samples

The H5 intentionally does not store test_size, test_true, test_predicted, rmse, or a standalone normalization field. The output filename includes the normalization used for training.

Small test:

```bash
NORMALIZATIONS=none,qbin50 \
GRID_PRESET=tiny \
MAX_PARAM_COMBOS=2 \
SELECT_K=50 \
EARLY_STOPPING=5 \
PROTEIN_START=1 \
PROTEIN_END=5 \
PROTEINS_PER_TASK=5 \
PARTITION=medium_cpu \
TIME=03:00:00 \
bash run_xgboost.sh
```

Full run defaults:

```bash
bash run_xgboost.sh
```

This defaults to:

```text
NORMALIZATIONS=none,zscore,rank,qbin50
PARTITION=compute
TIME=2-00:00:00
PROTEINS_PER_TASK=1000
SELECT_K=2000
EARLY_STOPPING=10
GRID_PRESET=fast
MAX_PARAM_COMBOS=8
```

Prediction from a saved H5:

```python
from joblib import load
from src.functions import load_minimal_8020_h5, predict_with_artifact

loaded = load_minimal_8020_h5("results_xgb/A2M_xgb_8020_none_fast_k2000.h5")
artifact = loaded["model_artifact"]
X_new = load("data/transcriptomics_filtered.joblib")
pred = predict_with_artifact(artifact, X_new)
print(pred.head())
```
