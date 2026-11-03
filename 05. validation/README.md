# 05. validation – external transfer to 20 proteogenomic cohorts

Data: the 20 cohorts compiled by Xu et al. 2022 (EMBL-EBI BioStudies **S-BSST733**). There are 12 TMT, 3 iTRAQ and 5 label-free cohorts, 2,560 samples in total.

| Track | Notebooks | In paper |
|---|---|---|
| Data preparation | `notebooks/1 → 2` | Collection A matrices, feature overlap (Tables S3–S4) |
| Exploratory LOO-model transfer | `notebooks/3 → 4 → 5` | Not in paper (uses LOO ElasticNet and an asinh/z-aligned target) |
| **Paper external validation** | `notebooks/6 → 7 → 8 → 9 → 10` (module `notebooks/paper_validation.py`) | Figs 5–7, S6–S9, Table 2 |
| scTranslator on external data | `scTranslator-validation-newdata/` | Not in paper (optional) |

## `data/` (the only data folder)

> Every entry of `data/` (and every output folder listed below) is a **link** into `research_data/`, created by `link_research_data.py` at the repo root. "Copy from X" in the table means the link points to the dataset stored once under X; nothing is duplicated.

```text
05. validation/data/
  Collection_A/<cohort>/protein_matrix_*.txt, RNA_matrix_*.txt   # 20 cohorts (S-BSST733), plus brain/rows_metadata.csv
  Table0_S1-S8_Combined20220114.xlsx                            # Xu et al. supplementary (notebook 10)
  transcriptomics_filtered.joblib, proteomics_matched.joblib    # copy from 01
  res_enet_zbig/, res_xgb_zbig/                                 # copy from 03 / 04 data/ (notebooks 6-10)
  loo/results_{zscore,none,ranking,binning}/                    # copy from 02. elasticNet/data/loo/ (notebooks 3-5)
  processed/        <- notebook 1 (Collection_A_<platform>_{RNA,Protein}.csv), notebook 2 (feature_lists/, feature_intersections/)
  val_data/         <- notebook 2 (aligned platform matrices, overlaps, *_common.joblib)
  proxy_models/     <- notebook 7 (proxy_metrics_by_tissue_protein.csv, read by notebook 8)
```
Run outputs: `results_external_validation_paper_v1/` (notebooks 6–10) and `notebooks/figures_pdf/<notebook>/` (all figures).

Cohort folder names must match exactly:

| Platform | Cohorts |
|---|---|
| TMT | colon_2019, renal_2019, Endomet_val, breast2020, pedbrain_2020, lung_2020_smoke, lung_val, lung2021, headneck_2021, pancrea2021, brain_2021, liver_val |
| iTRAQ | breast_2016, ovary_2016, stomach_2019 |
| LabelFree | brain, colon, prostate, liver, lung |

`data/Collection_B/` contains only licence and README files and is not used by any notebook.
