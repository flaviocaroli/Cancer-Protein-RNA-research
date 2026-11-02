# pipeline_xgb_loo – XGBoost leave-one-out (Suppl. S4.2)

## Inputs (`data/`)
`transcriptomics_filtered.joblib`, `proteomics_matched.joblib` (01 outputs). These are the same features as the ElasticNet LOO run.

## Run (from this folder)
```bash
NORMALIZATIONS=none RESULTS_DIR=$PWD/results_xgb_loo bash run_xgboost.sh
```
Output: `results_xgb_loo/<P>_loo_xgboost_<norm>_fast_once_k2000.h5` (+ `models/`).
Next: copy the `*.h5` files into `../model_results/final_data/`.

`xgb_loo_smoke_test.ipynb` tests a few proteins and samples and writes to `data/notebook_trial_outputs/`. `_mock_env/` is a test fixture.
