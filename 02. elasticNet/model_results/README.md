# model_results – ElasticNet LOO analysis and scTranslator on CCLE

All inputs and outputs are in **`../data/`**.

| # | Notebook | Reads (`../data/`) | Writes |
|---|---|---|---|
| 0 | `00_build_enet_all_results.ipynb` | `loo/results_zscore/*_LOO.parquet`, `*_best.parquet` | `../data/enet_all_results.h5`, `../data/enet_best_summary.csv` |
| 1 | `01. res_analysis_loo.ipynb` | `enet_all_results.h5`, `correlations.joblib`, `proteomics_matched.joblib`, `xgb_all_results.h5` (optional) | `protein_list.txt` (top 100), `per_protein_metrics_*.csv`, `per_sample_metrics_*.csv`, `exports/` |
| 2 | `02_enet_scT_comparison.ipynb` | `enet_all_results.h5`, `sc_predictions_*.npz`, `xgb_all_results.h5` | `comparison_plots/`, `outputs/` (Fig. S4) |
| 3 | `03. compute_sctranslator_predictions.ipynb` (run before 2; GPU recommended) | `transcriptomics_filtered.joblib`, `proteomics_matched.joblib`, scTranslator | `../data/sc_predictions_{2M,160K,10K}.npz` |

`proteins01.csv` and `pe85.csv` are optional legacy files. If present, they are read from `../data/`.

## scTranslator (notebook 3)
Put the scTranslator repo in `<repo>/external/scTranslator/scTranslator/` and the checkpoints in `<repo>/external/scTranslator/checkpoints/` (copy them from `research_data/00_external`). Overrides: `SCTRANSLATOR_REPO`, `SCTRANSLATOR_CKPT_DIR`, `CPRR_DATA_DIR`, `CPRR_PRED_OUT`.

## Figures
Every notebook shows its figures inline and saves them to `figures_pdf/<notebook name>/` next to the notebook. The first code cell (*pRNA figure routing*) handles this; keep it.
