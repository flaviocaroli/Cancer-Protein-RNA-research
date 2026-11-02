# 04. xgb / model_results – XGBoost LOO merge and analysis

## Inputs
| Path | Source |
|---|---|
| `final_data/*.h5` | `pipeline_xgb_loo/results_xgb_loo/*.h5` |
| `protein_list.txt` | Same list used for the LOO run (coverage check only) |
| `xgb_res/correlations.joblib`, `xgb_res/proteomics_matched.joblib` | `01. pre-processing` |

## Run order
| # | Notebook | Writes |
|---|---|---|
| 1 | `01_merge_results_loo.ipynb`: run the `final_data → data_h5` cells, then the **last cell** (`merge_h5_folder`) | `data_h5/*.h5`, **`xgb_res/xgb_all_results.h5`** |
| 2 | `02_res_analysis_loo_xgb.ipynb` | `exports/*.csv`, `figures_pdf/` (Fig. S10b) |

The alternative merge cells (`merged.h5`, `xgb_all_results2.h5`) are redundant and can be skipped.
`xgb_res/xgb_all_results.h5` is also an input of `03. elasticNet/model_results` (01 and 02). It is stored there as a copy in `research_data/03_elasticNet/input/model_results/res_data/`, so no manual step is needed.
