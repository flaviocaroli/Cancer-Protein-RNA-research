# new_pipeline_8020 – ElasticNet, fixed 80/20 split (primary)

One model per protein. Each model goes through: 80/20 split (295/74), transformation fitted on the training data, F-regression top-k features, 3-fold CV over at most 8 random `fast`-grid combinations (MSE), refit on the 295 samples, and evaluation on the 74.

## Inputs (`data/`)
`transcriptomics_filtered.joblib`, `proteomics_matched.joblib` – both from `01. pre-processing`.

## Runs used in the paper (settings from Suppl. S1.1)
```bash
# k = 2000, all four transformations  -> results_enet/
NORMALIZATIONS=none,zscore,ranking,binning N_BINS=50 SELECT_K=2000 \
GRID_PRESET=fast MAX_PARAM_COMBOS=8 SPLITS=3 FIT_INTERCEPT=0 bash run_enet.sh

# k = 5000, z-score only (primary model for external validation) -> res_enet_zbig/
NORMALIZATIONS=zscore SELECT_K=5000 GRID_PRESET=fast MAX_PARAM_COMBOS=8 SPLITS=3 \
FIT_INTERCEPT=0 RESULTS_DIR=res_enet_zbig bash run_enet.sh
```
These values are also the script defaults, except `NORMALIZATIONS`, `SELECT_K` and `RESULTS_DIR`.
Cluster settings: `PARTITION` (compute), `TIME` (3-00:00:00), `PROTEINS_PER_TASK` (1000), `CONDA_ENV` (xgb_loo).

## Outputs (one compressed H5 per protein and normalization)
```text
results_enet/<P>_enet_8020_{none,zscore,ranking}_fast_k2000.h5
results_enet/<P>_enet_8020_binning50_fast_k2000.h5     # bin count is in the name
res_enet_zbig/<P>_enet_8020_zscore_fast_k5000.h5
```
Each H5 holds: `protein`, `train_samples`, `test_samples`, `test_mse`, `test_pearson_corr`, `test_cosine_similarity`, `best_params`, `selected_features`, and the final model artifact (an `sklearn` ElasticNet plus its normalizer).
These folders are inputs to `05. ccle_res_analysis/` and `06. validation/` (notebooks 6–10).

## Test and monitor
```bash
NORMALIZATIONS=none,binning GRID_PRESET=tiny MAX_PARAM_COMBOS=2 SELECT_K=50 \
PROTEIN_START=1 PROTEIN_END=5 PROTEINS_PER_TASK=5 PARTITION=medium_cpu TIME=03:00:00 bash run_enet.sh
squeue -u $USER;  find results_enet -name "*.h5" | wc -l
```
`enet_8020_smoke_test.ipynb` runs the same code in-process on a few proteins.
