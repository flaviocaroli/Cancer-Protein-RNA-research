# 05. ccle_res_analysis – CCLE benchmark, cascade and protein classes

Covers Table 1, Figs 2–4 and Figs S2–S5. All models are compared on the same 74 held-out cell lines. The matched benchmark has 4,100 proteins, and the cascade universe has 5,113.

## Inputs (paths relative to this folder)
| Path | Source |
|---|---|
| `data/transcriptomics_filtered.joblib`, `data/proteomics_matched.joblib` | 01 |
| `data/full_pearson100.joblib`, `data/full_spearman100.joblib` | 02 `new_corr` (assembled) |
| `results_enet/`, `res_enet_zbig/` | 03 `new_pipeline_8020` |
| `results_xgb/`, `res_xgb_zbig/` | 04 `new_pipeline_xgb` (the 80/20 run, not LOO) |
| `sc_predictions_{2M,160K,10K}.npz` (this folder) | 03 `model_results/res_data/` (optional) |
| `pan_atlas_downloads/` | HPA, UniProt, GO, Reactome (downloaded by notebook 6); `wikipathways_human.gmt` is in the repo |

## Run order
| # | Notebook | Purpose | Main outputs |
|---|---|---|---|
| 0 | `0_proxies_enet_xgb_plots.ipynb` | Same-gene proxies vs models | `model_performance_test_plots_all_proxies/` |
| 1 | `1. comparison_models.ipynb` | Best model per protein | `comparison_outputs/` (incl. `best_model_per_protein_all_metrics.csv`), `output/` |
| 2 | `2. enet_xgb_plots.ipynb` | Detailed metric plots | `detailed_metric_plots/` |
| 3 | `3_proxy_enet_xgb_all_metrics.ipynb` | **Table 1 / Fig. 2** matched benchmark | `complete_comparison_outputs/tables/`, `figures_pdf/` |
| 4 | `4_ccle_cascade.ipynb` | **Fig. 3** cascade (proxy modes raw, train_zscore, train_linear) | `procedural_unified_tables/<mode>/…/cascade_partition_proteins_long.csv`, `figures_pdf/` |
| 6 | `6_pan_atlas.ipynb` | Downloads annotations; pan-atlas analysis | `pan_atlas_downloads/*`, `6_pan_atlas_outputs/` |
| 5 | `5_HPA_heatmaps.ipynb` | HPA class coverage | `4.2_HPA_heatmaps_outputs/` |
| 5_1 | `5_1_hpa_barplots.ipynb` | **Fig. 4 / S5** HPA enrichment and Pearson-vs-MSE | `4.0_hpa_v4_outputs/` |

Run **6 before 5 and 5_1**, because they need `pan_atlas_downloads/hpa_proteinatlas.tsv.zip`.

In notebooks 5 and 5_1, set `CASCADE_LONG_FILES` to the paper's primary mode (train-fitted linear proxy):
```python
CASCADE_LONG_FILES = ['procedural_unified_tables/train_linear/protein_names_by_cascade_partition/cascade_partition_proteins_long.csv']
```
Otherwise the auto-search pools the three proxy modes under one source id.

## Notes
- Proxy eligibility must use training samples only. `4_ccle_cascade` does this (`same_gene_corr_train`). `2. enet_xgb_plots` and `0_proxies` still read the all-sample `full_*100` lookups. Treat their proxy subsets as descriptive (see `src/utils_ccle.py`, LEAKAGE NOTE).
- The output descriptions in `cascade_partition_outputs_metadata.md` still hold. `procedural_ccle.ipynb` is now `4_ccle_cascade.ipynb`, and `9.procedural_cptac.ipynb` is now `06/notebooks/9_external_cascade.ipynb`.
