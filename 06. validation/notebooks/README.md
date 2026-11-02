# 06. validation / notebooks

Run from this folder. Paths are relative (`../data`, `../processed_data`, `../val_data`, `../results_*`). Shared helpers are in `src/utils_val.py`.

## Data preparation
| # | Notebook | Reads | Writes |
|---|---|---|---|
| 1 | `1. process_collection_a.ipynb` | `../data/Collection_A/<cohort>/*` | `../processed_data/Collection_A_{TMT,iTRAQ,LabelFree}_{Protein,RNA}.csv`. Samples are renamed `<COHORT>_<i>`, and features with >50% NaN are dropped. |
| 2 | `2. comparison.ipynb` | `../processed_data/` (CCLE joblibs + Collection A CSVs) | `../val_data/{TMT,iTRAQ,LabelFree}_{transcriptomics,proteomics}.joblib` (aligned to CCLE: RNA robust-z, protein asinh + z), `../val_data/{genes,proteins}/overlap_*.txt`, `../val_data/transcriptomics_common.joblib`, `proteomics_common.joblib`, `../processed_data/feature_lists/*`, `../model_results/feature_intersections/*` |

## Exploratory LOO-model transfer (not part of the paper)
| # | Notebook | Reads | Writes |
|---|---|---|---|
| 3 | `3.results_cptac_loo.ipynb` | `../results_val/*_best.(parquet|csv)`, `batch_*`, `../val_data/*` | `../results_val/plots_cptac_unified/*` |
| 4 | `4.results_norms_loo.ipynb` | `../results_{val,none,rk,bin}/*_best.*` | `../results_val/plots_norms_unified/*`, `<results_dir>/plots/` |
| 5 | `5. comparison_models_loo.ipynb` | same four folders | `../results_comparison_all_models/model_comparison_bars.png` |

The `results_*` folders are copies of `03. elasticNet/new_pipeline_loo/data/results_{zscore,none,ranking,binning}/`. Parquet and CSV are both read.

## Paper external validation (Figs 5–7, S6–S9)
Primary analysis, as in the paper's Methods and S5:
- External RNA is completed with CCLE-development (295-sample) gene medians.
- RNA is z-scored with CCLE-development mean and SD for the z-score models.
- Proteins stay on cohort-native scales, without imputation.
- Models are the z-score k=5000 ElasticNet and XGBoost, with no refitting.
- The same-gene reference is the best of raw, CCLE-z and CCLE-train linear regression, chosen per protein, cohort and metric (a descriptive upper envelope).

Clipping, sample-scale alignment and TMT1/TMT2 are sensitivity analyses only.

| # | Notebook | Writes (under `PROJECT_ROOT/results_external_validation_paper_v1/`) |
|---|---|---|
| 6 | `6_normalization.ipynb` | Frozen contract: `validation_manifest.json`, `contract/`, `cohorts/*.joblib`, normalisation figures (Fig. S6) |
| 7 | `7_inference_all_models.ipynb` | Predictions and metrics per cohort (`evaluation/`) |
| 8 | `8_external_plots.ipynb` | Table 2, Fig. 5, S8, S9 (`tables/`, `notebooks/figures_pdf/`) |
| 9 | `9_external_cascade.ipynb` | Fig. 6, S7 (`procedural/`) |
| 10 | `10_xu_contextual_comparison.ipynb` | Fig. 7 (`xu/`); needs `Table0_S1-S8_Combined20220114.xlsx` |

These notebooks need the module **`paper_validation.py`** next to them, which is not yet in the repo. Notebook 8 also needs `results_comparison_all_models_v4/proxy_models_v4/proxy_metrics_by_tissue_protein.csv`, or `PAPER_VALIDATION_PROXY_TISSUE_PROTEIN_PATH`.
Built-in sanity checks: 369 CCLE samples (295/74), 5,113 proteins, 20 cohorts.

## Dependencies
Root `requirements.txt`. `mygene` is optional and only needed for online Ensembl mapping.
