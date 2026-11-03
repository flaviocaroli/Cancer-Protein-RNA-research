# 04. ccle_res_analysis – CCLE benchmark, cascade and protein classes

Covers Table 1, Figs 2–4 and Figs S2–S5. All models are compared on the same 74 held-out cell lines. The matched benchmark has 4,100 proteins, and the cascade universe has 5,113.

## `data/` (the only data folder)

> Every entry of `data/` (and every output folder listed below) is a **link** into `research_data/`, created by `link_research_data.py` at the repo root. "Copy from X" in the table means the link points to the dataset stored once under X; nothing is duplicated.

| Path | What | Source |
|---|---|---|
| `transcriptomics_filtered.joblib`, `proteomics_matched.joblib` | CCLE matrices | copy from 01 |
| `results_enet/`, `res_enet_zbig/` | ENet 80/20 H5 | copy from `02. elasticNet/data/` |
| `results_xgb/`, `res_xgb_zbig/` | XGB 80/20 H5 (not LOO) | copy from `03. xgb/data/` |
| `sc_predictions_{2M,160K,10K}.npz` | Optional | copy from `02. elasticNet/data/` |
| `pan_atlas_downloads/` | HPA, UniProt, GO, Reactome (downloaded by `5_hpa_pan_atlas` Part A); `wikipathways_human.gmt` is in the repo | — |

Analysis outputs (tables) are written to their own folders next to the notebooks, as listed below.

## Run order
| # | Notebook | Purpose | Main outputs |
|---|---|---|---|
| 0 | `0_proxies_enet_xgb_plots.ipynb` | Same-gene proxies vs models | `model_performance_test_plots_all_proxies/` |
| 1 | `1. comparison_models.ipynb` | Best model per protein | `comparison_outputs/` (incl. `best_model_per_protein_all_metrics.csv`), `output/` |
| 2 | `2. enet_xgb_plots.ipynb` | Detailed metric plots | `detailed_metric_plots/` |
| 3 | `3_proxy_enet_xgb_all_metrics.ipynb` | **Table 1 / Fig. 2** matched benchmark | `complete_comparison_outputs/tables/`, `figures_pdf/` |
| 4 | `4_ccle_cascade.ipynb` | **Fig. 3** cascade (proxy modes raw, train_zscore, train_linear) | `procedural_unified_tables/<mode>/…/cascade_partition_proteins_long.csv`, `figures_pdf/` |
| 5 | `5_hpa_pan_atlas.ipynb` – Part A (former 6_pan_atlas) | Downloads annotations; pan-atlas analysis | `pan_atlas_downloads/*`, `6_pan_atlas_outputs/` |
| | Part B (former 5_HPA_heatmaps) | HPA class coverage | `4.2_HPA_heatmaps_outputs/` |
| | Part C (former 5_1_hpa_barplots) | **Fig. 4 / S5** HPA enrichment and Pearson-vs-MSE | `4.0_hpa_v4_outputs/` |

`5_hpa_pan_atlas.ipynb` merges the three former notebooks in dependency order. Helpers that were identical are defined once, and duplicate imports are removed.
Parts B and C use the paper's primary proxy mode (`procedural_unified_tables/train_linear/…`) automatically. If that file is missing, they fall back to auto-search.

## Notes
- Same-gene correlations are computed on the 295 training samples inside notebooks 0, 2 and 4 (`same_gene_corr_train`). The precomputed `full_pearson100` / `full_spearman100` matrices are no longer needed.
- The output descriptions in `cascade_partition_outputs_metadata.md` still hold. `procedural_ccle.ipynb` is now `4_ccle_cascade.ipynb`, and `9.procedural_cptac.ipynb` is now `05. validation/notebooks/9_external_cascade.ipynb`.

## Figures
Every notebook shows its figures inline and saves them to `figures_pdf/<notebook name>/` next to the notebook. The first code cell (*pRNA figure routing*) redirects every `savefig` there, keeping the original sub-folders. Do not delete that cell.
