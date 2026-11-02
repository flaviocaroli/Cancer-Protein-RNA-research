# new_pipeline_xgb – XGBoost, fixed 80/20 split (primary)

Per protein: 80/20 split (295/74), absolute-Spearman top-k features on the training data, and 3-fold CV over at most 8 `fast`-grid combinations with early stopping 10. The final number of trees is max(25, min(n, ⌈1.1·median best iteration⌉)). The model is refit on 295 samples and tested on 74 (Suppl. S1.2).

## Inputs (`data/`)
`transcriptomics_filtered.joblib`, `proteomics_matched.joblib` (01 outputs).

## Runs used in the paper
```bash
# k = 2000, four transformations -> results_xgb/
NORMALIZATIONS=none,zscore,rank,qbin50 SELECT_K=2000 GRID_PRESET=fast MAX_PARAM_COMBOS=8 \
SPLITS=3 EARLY_STOPPING=10 PARTITION=compute TIME=2-00:00:00 bash run_xgboost.sh

# k = 5000, z-score (primary model for external validation) -> res_xgb_zbig/
NORMALIZATIONS=zscore SELECT_K=5000 GRID_PRESET=fast MAX_PARAM_COMBOS=8 SPLITS=3 \
EARLY_STOPPING=10 RESULTS_DIR=res_xgb_zbig PARTITION=compute TIME=2-00:00:00 bash run_xgboost.sh
```
Script defaults: `PARTITION=medium_cpu`, `TIME=3:00:00`, `PROTEINS_PER_TASK=1000`, `RANDOM_STATE=42`, `USE_GPU=0`. Set `USE_GPU=1` for CUDA.

## Outputs
```text
results_xgb/<P>_xgb_8020_{none,zscore,rank,qbin50}_fast_k2000.h5
res_xgb_zbig/<P>_xgb_8020_zscore_fast_k5000.h5
```
Each H5 holds: protein, train/test samples, test MSE / Pearson / cosine, best params, selected features, and the final model artifact (`xgboost.XGBRegressor` plus normalizer).
These folders are inputs to `05. ccle_res_analysis/` and `06. validation/` (6–10).

```python
from joblib import load
from src.functions import load_minimal_8020_h5, predict_with_artifact
art = load_minimal_8020_h5("results_xgb/A2M_xgb_8020_none_fast_k2000.h5")["model_artifact"]
pred = predict_with_artifact(art, load("data/transcriptomics_filtered.joblib"))
```
Small test: `NORMALIZATIONS=none,qbin50 GRID_PRESET=tiny MAX_PARAM_COMBOS=2 SELECT_K=50 EARLY_STOPPING=5 PROTEIN_START=1 PROTEIN_END=5 PROTEINS_PER_TASK=5 bash run_xgboost.sh`
