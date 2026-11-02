# 06. validation – external transfer to 20 proteogenomic cohorts

Data: the 20 cohorts compiled by Xu et al. 2022 (EMBL-EBI BioStudies **S-BSST733**). There are 12 TMT, 3 iTRAQ and 5 label-free cohorts, 2,560 samples in total.

| Track | Notebooks | In paper |
|---|---|---|
| Data preparation | `notebooks/1 → 2` | Collection A matrices, feature overlap (Tables S3–S4) |
| Exploratory LOO-model transfer | `notebooks/3 → 4 → 5` | Not in paper (uses LOO ElasticNet and an asinh/z-aligned target) |
| **Paper external validation** | `notebooks/6 → 7 → 8 → 9 → 10` | Figs 5–7, S6–S9, Table 2 |
| scTranslator on external data | `scTranslator-validation-newdata/` | Not in paper (optional) |

## Folder layout (create it from `research_data/06_validation/input/`)
```text
06. validation/
  data/Collection_A/<cohort>/protein_matrix_*.txt, RNA_matrix_*.txt     # 20 cohorts, from S-BSST733
  data/Collection_A/brain/rows_metadata.csv                            # Ensembl -> symbol for brain
  data/Table0_S1-S8_Combined20220114.xlsx                              # Xu et al. supplementary (notebook 10)
  processed_data/transcriptomics_filtered.joblib, proteomics_matched.joblib   # from 01
  res_enet_zbig/, res_xgb_zbig/                                        # from 03 / 04 80/20 runs (notebooks 6-10)
  results_val/ results_none/ results_rk/ results_bin/                  # from 03 LOO (notebooks 3-5)
  notebooks/  scTranslator-validation-newdata/
```
Cohort folder names must match exactly:

| Platform | Cohorts |
|---|---|
| TMT | colon_2019, renal_2019, Endomet_val, breast2020, pedbrain_2020, lung_2020_smoke, lung_val, lung2021, headneck_2021, pancrea2021, brain_2021, liver_val |
| iTRAQ | breast_2016, ovary_2016, stomach_2019 |
| LabelFree | brain, colon, prostate, liver, lung |

`data/Collection_B/` contains only licence and README files and is not used by any notebook.
