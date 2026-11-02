# model_results – ElasticNet LOO analysis and scTranslator on CCLE

## Run order
| # | Notebook | Reads | Writes |
|---|---|---|---|
| 0 | `00_build_enet_all_results.ipynb` | `../new_pipeline_loo/data/results_<NORMALIZATION>/*_LOO.parquet` and `*_best.parquet` (default `zscore`) | `res_data/enet_all_results.h5` (key `data`: protein, LOO sample, true, predicted, params_dict, …), `res_data/enet_best_summary.csv` |
| 1 | `01. res_analysis_loo.ipynb` | `res_data/enet_all_results.h5`, `res_data/correlations.joblib` (01), `res_data/proteomics_matched.joblib` (01), `res_data/protein_clusters_k3.csv` (02 `results/clusters_pg/`), `res_data/xgb_all_results.h5` (optional, 04) | `protein_list.txt` (top 100), `per_protein_metrics_*.csv`, `per_sample_metrics_*.csv`, `exports/*.csv` |
| 2 | `02_enet_scT_comparison.ipynb` | `res_data/enet_all_results.h5`, `res_data/sc_predictions_{2M,160K,10K}.npz`, `res_data/xgb_all_results.h5` | `elasticnet_*_{raw,scaled}.png`, `comparison_plots/`, `outputs/`, `figures_pdf/` (Fig. S4) |
| 3 | `03. compute_sctranslator_predictions.ipynb` (run before 2; GPU recommended) | `data/transcriptomics_filtered.joblib`, `data/proteomics_matched.joblib`, scTranslator code and checkpoints | `res_data/sc_predictions_{2M,160K,10K}.npz` |

`*_best` files alone cannot rebuild `enet_all_results.h5`, because they hold only one LOO split per protein. Run the LOO pipeline with `FULL_OUTPUT=1` (the default).
`res_data/proteins01.csv` and `res_data/pe85.csv` are no longer needed. The XGB set is recomputed from `xgb_all_results.h5`, and `pe85` was never used.

## scTranslator files (notebook 3)
Copy `research_data/00_external/input/scTranslator/` to `<repo>/external/scTranslator/`:
```text
external/scTranslator/scTranslator/code/model/performer_enc_dec.py  + code/model/ID_dic/*.pkl
external/scTranslator/checkpoints/scTranslator_{2M,160K,10K}.pt
```
These are found automatically. You can override them with `SCTRANSLATOR_REPO`, `SCTRANSLATOR_CKPT_DIR`, `CPRR_DATA_DIR` or `CPRR_PRED_OUT`.
Checkpoints are loaded with `weights_only=False`, as required by PyTorch ≥ 2.6. Only use your own trusted files.
Columns with more than 50% NaN are dropped and the rest are filled with 1e-8. Symbols missing from the scTranslator vocabulary are discarded.
The 10K checkpoint is the paper's primary scTranslator reference (Suppl. S3). For a quick test, keep only `10K` in `CONFIG['checkpoints']`.
