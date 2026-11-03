# 05. validation / notebooks

Run from this folder. Paths are relative (`../data`, `../processed_data`, `../val_data`, `../results_*`). Shared helpers are in `src/utils_val.py`.

## Data preparation
| # | Notebook | Reads | Writes |
|---|---|---|---|
| 1 | `1. process_collection_a.ipynb` | `../data/Collection_A/<cohort>/*` | `../data/processed/Collection_A_{TMT,iTRAQ,LabelFree}_{Protein,RNA}.csv`. Samples are renamed `<COHORT>_<i>`, and features with >50% NaN are dropped. |
| 2 | `2. comparison.ipynb` | `../data/*.joblib` (CCLE) and `../data/processed/` | `../data/val_data/` (aligned platform matrices, overlaps, `*_common.joblib`), `../data/processed/feature_lists/`, `../data/processed/feature_intersections/` |

## Exploratory LOO-model transfer (not part of the paper)
| # | Notebook | Reads | Writes |
|---|---|---|---|
| 3 | `3.results_cptac_loo.ipynb` | `../data/loo/results_zscore/*_best.*`, `batch_*`, `../data/val_data/*` | tables in `../data/loo/results_zscore/plots_cptac_unified/` |
| 4 | `4.results_norms_loo.ipynb` | `../data/loo/results_{zscore,none,ranking,binning}/*_best.*` | tables in `../data/loo/results_zscore/plots_norms_unified/` |
| 5 | `5. comparison_models_loo.ipynb` | same four folders | `figures_pdf/5. comparison_models_loo/…/model_comparison_bars.png` |

`../data/loo/` is a copy of `02. elasticNet/data/loo/`. Parquet and CSV are both read.

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

All five notebooks import **`paper_validation.py`** (in this folder), which holds the shared logic:

| Area | Functions |
|---|---|
| Cohorts | `COHORT_TABLE`, the 20 cohorts with platform and Xu short name |
| Model loading | Reads the saved CCLE models (`data/res_enet_zbig/`, `data/res_xgb_zbig/`) and checks they all share one 80/20 split |
| CCLE reference | Statistics from the 295 training samples: medians, mean/SD, clip quantiles |
| Proxy parameters | Same-gene proxy fits (raw, z-score, linear), training samples only |
| Cohort contracts | Per-cohort matched samples with the >50% missingness filters |
| Inference | `evaluate_all_cohorts_memory_safe` (resumable) |
| Metrics | Pearson, Spearman, MSE, RMSE, normalised RMSE, cosine, calibration slope/intercept |
| Universes and summaries | Matched universe, bootstrap CIs, pooled-platform metrics |

Scenarios are `primary_raw` (paper), `clipped` and `scale_aligned` (sensitivity only).
TMT cohorts whose native protein median is above 5 form the `TMT2` reporting stratum, and the rest form `TMT1`.
Notebook 7 also writes `../data/proxy_models/proxy_metrics_by_tissue_protein.csv`, which notebook 8 reads.
Useful environment variables: `PAPER_VALIDATION_SMOKE_TEST=1` (relaxes the 369/74/5113 checks), `PAPER_VALIDATION_ALLOW_MISSING_MODELS=1`, `PAPER_VALIDATION_FORCE=1` (restart notebook 7), `XU_XLSX_PATH`.
Built-in sanity checks: 369 CCLE samples (295/74), 5,113 proteins, 20 cohorts.

## Dependencies
Root `requirements.txt`. `mygene` is optional and only needed for online Ensembl mapping.

## Figures
Every notebook shows its figures inline and saves them to `figures_pdf/<notebook name>/` next to the notebook. The first code cell (*pRNA figure routing*) redirects every `savefig` there, keeping the original sub-folders. Do not delete that cell.
