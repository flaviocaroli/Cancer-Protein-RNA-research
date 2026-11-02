# 03. elasticNet

| Subfolder | Role in paper | Main outputs |
|---|---|---|
| `new_pipeline_8020/` | **Primary** ElasticNet benchmark: 295/74 split, 3-fold CV (Table 1, Figs 2–3) | `results_enet/`, `res_enet_zbig/` (one H5 per protein) |
| `new_pipeline_loo/` | LOO-CV robustness analysis (Suppl. S4.2, Fig. S10a) | `data/results_<norm>/<P>_LOO.parquet`, `<P>_best.parquet` |
| `model_results/` | LOO analysis + scTranslator CCLE predictions | `res_data/enet_all_results.h5`, `res_data/sc_predictions_*.npz`, figures |

Order: `new_pipeline_8020` → (`model_results/03` for scTranslator) → `new_pipeline_loo` → `model_results/00 → 01 → 02`.
All three use the same RNA matrix, `transcriptomics_filtered.joblib` (22,614 genes).
Install: `pip install -r requirements.txt`.
