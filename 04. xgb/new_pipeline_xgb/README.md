# new_pipeline_xgb – XGBoost, fixed 80/20 split (primary)

Per protein: 80/20 split (295/74), absolute-Spearman top-k features, 3-fold CV over at most 8 `fast` combinations with early stopping 10. Final trees = max(25, min(n, ⌈1.1·median best iteration⌉)). Refit on 295 samples, test on 74 (Suppl. S1.2).
Reads and writes **`../data/`**.

```bash
NORMALIZATIONS=none,zscore,rank,qbin50 SELECT_K=2000 GRID_PRESET=fast MAX_PARAM_COMBOS=8 \
SPLITS=3 EARLY_STOPPING=10 PARTITION=compute TIME=2-00:00:00 bash run_xgboost.sh              # -> ../data/results_xgb/

NORMALIZATIONS=zscore SELECT_K=5000 GRID_PRESET=fast MAX_PARAM_COMBOS=8 SPLITS=3 EARLY_STOPPING=10 \
RESULTS_DIR=$PWD/../data/res_xgb_zbig PARTITION=compute TIME=2-00:00:00 bash run_xgboost.sh   # -> ../data/res_xgb_zbig/
```
Defaults: `DATA_DIR=../data`, `RESULTS_DIR=../data/results_xgb`, `PARTITION=medium_cpu`, `TIME=3:00:00`, `USE_GPU=0`.
Outputs: `<P>_xgb_8020_<norm>_fast_k2000.h5` and `<P>_xgb_8020_zscore_fast_k5000.h5`. Copy them to `05. ccle_res_analysis/data/` and `06. validation/data/`.
