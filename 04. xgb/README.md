# 04. xgb

| Subfolder | Role in paper | Main outputs |
|---|---|---|
| `new_pipeline_xgb/` | **Primary** XGBoost benchmark: 295/74 split, 3-fold CV (Table 1, Figs 2–3) | `results_xgb/`, `res_xgb_zbig/` |
| `pipeline_xgb_loo/` | LOO-CV robustness (Suppl. S4.2, Fig. S10b) | `results_xgb_loo/<P>_loo_xgboost_<norm>_fast_once_k2000.h5` |
| `model_results/` | Merge and analyse the LOO results | `xgb_res/xgb_all_results.h5`, figures |

Both pipelines read `data/transcriptomics_filtered.joblib` and `data/proteomics_matched.joblib` from `01. pre-processing`.
Both default to a folder called `results_xgb/`. Always give the LOO run its own folder (`RESULTS_DIR=results_xgb_loo`).
