# new_pipeline_8020 – ElasticNet, fixed 80/20 split (primary)

Per protein: 80/20 split (295/74), transformation fitted on the training data, F-regression top-k features, 3-fold CV over at most 8 random `fast`-grid combinations (MSE), refit on 295 samples, test on 74 (Suppl. S1.1).
Reads and writes **`../data/`** (the single `03. elasticNet/data` folder).

## Runs used in the paper (from this folder)
```bash
NORMALIZATIONS=none,zscore,ranking,binning N_BINS=50 SELECT_K=2000 \
GRID_PRESET=fast MAX_PARAM_COMBOS=8 SPLITS=3 FIT_INTERCEPT=0 bash run_enet.sh                 # -> ../data/results_enet/

NORMALIZATIONS=zscore SELECT_K=5000 GRID_PRESET=fast MAX_PARAM_COMBOS=8 SPLITS=3 FIT_INTERCEPT=0 \
RESULTS_DIR=$PWD/../data/res_enet_zbig bash run_enet.sh                                        # -> ../data/res_enet_zbig/
```
Defaults: `DATA_DIR=../data`, `RESULTS_DIR=../data/results_enet`, logs in `./logs`. Cluster settings: `PARTITION`, `TIME`, `PROTEINS_PER_TASK`, `CONDA_ENV`.

## Outputs
```text
../data/results_enet/<P>_enet_8020_{none,zscore,ranking}_fast_k2000.h5, <P>_enet_8020_binning50_fast_k2000.h5
../data/res_enet_zbig/<P>_enet_8020_zscore_fast_k5000.h5
```
Each H5 holds the split, test metrics, best params, selected features and the final model artifact. These folders are copied to `05. ccle_res_analysis/data/` and `06. validation/data/`.
`enet_8020_smoke_test.ipynb` runs the same code in-process on a few proteins.
