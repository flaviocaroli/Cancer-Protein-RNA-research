# pipeline_xgb_loo – XGBoost leave-one-out (Suppl. S4.2)

Reads `../data/transcriptomics_filtered.joblib` and `../data/proteomics_matched.joblib` (the same features as the ElasticNet LOO run).
```bash
NORMALIZATIONS=none bash run_xgboost.sh        # -> ../data/results_xgb_loo/<P>_loo_xgboost_<norm>_fast_once_k2000.h5
```
`model_results/01_merge_results_loo.ipynb` reads that folder directly, so there is no copy step.
`xgb_loo_smoke_test.ipynb` runs a quick test, and `_mock_env/` is a test fixture.
