# 04. xgb

Code is in the subfolders. **All data lives in one folder: `data/`.**

## `data/` (the only data folder)
| Path | What | Source / producer |
|---|---|---|
| `transcriptomics_filtered.joblib`, `proteomics_matched.joblib`, `correlations.joblib` | CCLE inputs | copy from `01. pre-processing/data/` |
| `protein_list.txt` | Protein list of the LOO run (coverage check) | your list or `02. correlations/data/protein_list.txt` |
| `results_xgb/` | 80/20 k=2000 (none, zscore, rank, qbin50) | `new_pipeline_xgb` |
| `res_xgb_zbig/` | 80/20 z-score k=5000 (**primary model**) | `new_pipeline_xgb` |
| `results_xgb_loo/` | LOO per-protein H5 (+ `models/`) | `pipeline_xgb_loo` |
| `data_h5/`, `xgb_all_results.h5` | Standardised and merged LOO | `model_results/01` |

Order: `new_pipeline_xgb` → `pipeline_xgb_loo` → `model_results/01 → 02`.
