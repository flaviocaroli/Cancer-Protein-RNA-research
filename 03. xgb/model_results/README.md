# 03. xgb / model_results – XGBoost LOO merge and analysis

All inputs and outputs are in **`../data/`**.

| # | Notebook | Reads (`../data/`) | Writes |
|---|---|---|---|
| 1 | `01_merge_results_loo.ipynb`: run the standardisation cells, then the **last cell** | `results_xgb_loo/*.h5`, `protein_list.txt` | `../data/data_h5/*.h5`, **`../data/xgb_all_results.h5`**, coverage `.txt` files |
| 2 | `02_res_analysis_loo_xgb.ipynb` | `xgb_all_results.h5`, `correlations.joblib`, `proteomics_matched.joblib` | `exports/` (Fig. S10b) |

The alternative merge cells (`merged.h5`, `xgb_all_results2.h5`) are optional. `xgb_all_results.h5` is copied to `02. elasticNet/data/`.

## Figures
Every notebook shows its figures inline and saves them to `figures_pdf/<notebook name>/` next to the notebook. The first code cell (*pRNA figure routing*) handles this; keep it.
