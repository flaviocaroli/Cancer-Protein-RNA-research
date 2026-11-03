# 03. xgb

Code is in the subfolders. **All data lives in one folder: `data/`.**

## `data/` (the only data folder)

> Every entry of `data/` (and every output folder listed below) is a **link** into `research_data/`, created by `link_research_data.py` at the repo root. "Copy from X" in the table means the link points to the dataset stored once under X; nothing is duplicated.

| Path | What | Source / producer |
|---|---|---|
| `transcriptomics_filtered.joblib`, `proteomics_matched.joblib`, `correlations.joblib` | CCLE inputs | copy from `01. pre-processing/data/` |
| `protein_list.txt` | Protein list of the LOO run (coverage check) | your list (same as `02. elasticNet/data/protein_list.txt`) |
| `results_xgb/` | 80/20 k=2000 (none, zscore, rank, qbin50) | `new_pipeline_xgb` |
| `res_xgb_zbig/` | 80/20 z-score k=5000 (**primary model**) | `new_pipeline_xgb` |
| `results_xgb_loo/` | LOO per-protein H5 (+ `models/`) | `pipeline_xgb_loo` |
| `data_h5/`, `xgb_all_results.h5` | Standardised and merged LOO | `model_results/01` |

Order: `new_pipeline_xgb` → `pipeline_xgb_loo` → `model_results/01 → 02`.
