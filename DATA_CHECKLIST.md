# pRNA – what goes in each `data/` folder (paper branch)

Every stage has **one** `data/` folder. The tables below list what must be in it before you run that stage's notebooks and scripts, and what the stage writes there itself.

| Source | Meaning |
|---|---|
| **you** | You provide it (download or your own file) |
| **copy from X** | Produced by stage X; copy it here |
| **produced** | Written by this stage; nothing to do |

Before you start, set VS Code **Settings → Jupyter: Notebook File Root** to `${fileDirname}` so every notebook runs from its own folder.

---

## 01. pre-processing/data/ (done)

| File | Source | Used by |
|---|---|---|
| `CCLE_expression_full.csv` | you (figshare 25797014; auto-downloaded) | 01 |
| `Table_S2_Protein_Quant_Normalized.xlsx` | you (Gygi lab; auto-downloaded) | 01 |
| `sample_info.csv` | you (figshare 25494443; auto-downloaded) | 01 |
| `transcriptomics.csv`, `proteomics.csv` | produced | – |
| `transcriptomics_matched.joblib`, `proteomics_matched.joblib` | produced | 01/02_filtering |
| `transcriptomics_filtered.joblib` | produced | **02, 03, 04, 05** |
| `proteomics_matched.joblib` | produced | **02, 03, 04, 05** |
| `correlations.joblib` | produced | **02, 03** |

---

## 02. elasticNet/data/

| File / folder | Source | Used by |
|---|---|---|
| `transcriptomics_filtered.joblib` | copy from 01 | new_pipeline_8020, new_pipeline_loo, model_results/03 |
| `proteomics_matched.joblib` | copy from 01 | new_pipeline_8020, new_pipeline_loo, model_results/01, 03 |
| `correlations.joblib` | copy from 01 | model_results/01 |
| `protein_list.txt` | you (LOO protein list) | new_pipeline_loo |
| `scTranslator/` (contains `code/model/performer_enc_dec.py` and `code/model/ID_dic/`) | you (`git clone https://github.com/TencentAILabHealthcare/scTranslator.git`) | model_results/03 |
| `scTranslator_2M.pt`, `scTranslator_160K.pt`, `scTranslator_10K.pt` | you (Drive links in MANIFEST.xlsx) | model_results/03 |
| `xgb_all_results.h5` | copy from 03 (or let 02_enet_scT_comparison read it from `03. xgb/data/`) | model_results/01 (optional), 02 |
| `proteins01.csv`, `pe85.csv` | you – optional legacy, not needed | model_results/01 |
| `results_enet/` | produced (80/20, k=2000: none, zscore, ranking, binning50) | → copy to 04 |
| `res_enet_zbig/` | produced (80/20, z-score, k=5000, primary model) | → copy to 04, 05 |
| `loo/results_zscore/`, `loo/results_none/`, `loo/results_ranking/`, `loo/results_binning/` | produced by new_pipeline_loo (`*_best`, `*_LOO`, `batch_*`) | model_results/00; → copy to 05 |
| `enet_all_results.h5` | produced by model_results/00, **or your original file** | model_results/01, 02 |
| `enet_best_results.h5` | produced by model_results/00 | reference |
| `sc_predictions_2M.npz`, `sc_predictions_160K.npz`, `sc_predictions_10K.npz` | produced by model_results/03 | model_results/02; → copy to 04 |

**Order:** new_pipeline_8020 → model_results/03 → (new_pipeline_loo) → model_results/00 → 01 → 02

---

## 03. xgb/data/

| File / folder | Source | Used by |
|---|---|---|
| `transcriptomics_filtered.joblib` | copy from 01 | new_pipeline_xgb, pipeline_xgb_loo |
| `proteomics_matched.joblib` | copy from 01 | new_pipeline_xgb, pipeline_xgb_loo, model_results/02 |
| `correlations.joblib` | copy from 01 | model_results/02 |
| `protein_list.txt` | you (same LOO list as 02) | model_results/01 (coverage check) |
| `results_xgb/` | produced (80/20, k=2000: none, zscore, rank, qbin50) | → copy to 04 |
| `res_xgb_zbig/` | produced (80/20, z-score, k=5000, primary model: `<P>_xgb_8020_zscore_fast_k5000.h5`) | → copy to 04, 05 |
| `results_xgb_loo/` | produced by pipeline_xgb_loo (`<P>_loo_xgboost_<norm>_fast_once_k2000.h5`; this is the old `final_data/`) | model_results/01 |
| `data_h5/` | produced by model_results/01 | model_results/01 |
| `xgb_all_results.h5` | produced by model_results/01 (last cell) | model_results/02; → copy to 02 |

**Order:** new_pipeline_xgb → pipeline_xgb_loo → model_results/01 → 02

---

## 04. ccle_res_analysis/data/

| File / folder | Source | Used by |
|---|---|---|
| `transcriptomics_filtered.joblib` | copy from 01 | notebooks 0–4 |
| `proteomics_matched.joblib` | copy from 01 | notebooks 0–4 |
| `results_enet/` | copy from 02 | notebooks 0–4 |
| `res_enet_zbig/` | copy from 02 | notebooks 0–4 |
| `results_xgb/` | copy from 03 (80/20, **not** LOO) | notebooks 0–4 |
| `res_xgb_zbig/` | copy from 03 | notebooks 0–4 |
| `sc_predictions_2M.npz`, `sc_predictions_160K.npz`, `sc_predictions_10K.npz` | copy from 02 (optional: scTranslator columns) | notebooks 3, 4 |
| `pan_atlas_downloads/` | auto-downloaded by 5_hpa_pan_atlas Part A; `wikipathways_human.gmt` is in the repo | 5_hpa_pan_atlas |

**Order:** 0 → 1 → 2 → 3 → 4 → 5_hpa_pan_atlas

---

## 05. validation/data/

| File / folder | Source | Used by |
|---|---|---|
| `Collection_A/<cohort>/protein_matrix_*.txt` + `RNA_matrix_*.txt` (20 cohorts) | you (BioStudies S-BSST733) | notebook 1 |
| `Collection_A/brain/rows_metadata.csv` | you (ships with the brain cohort) | notebook 1 |
| `Table0_S1-S8_Combined20220114.xlsx` | you (Xu et al. 2022 supplementary) | notebook 10 |
| `transcriptomics_filtered.joblib` | copy from 01 | notebooks 2, 6 |
| `proteomics_matched.joblib` | copy from 01 | notebooks 2, 6 |
| `res_enet_zbig/` | copy from 02 (all 5,113 proteins) | notebooks 6–10 |
| `res_xgb_zbig/` | copy from 03 (all 5,113 proteins, same 80/20 split) | notebooks 6–10 |
| `loo/results_zscore/`, `loo/results_none/`, `loo/results_ranking/`, `loo/results_binning/` | copy from 02 (`*_best` files) | notebooks 3, 4, 5 |
| `processed/` | produced by notebooks 1, 2 | notebooks 2, 6 |
| `val_data/` | produced by notebook 2 | notebooks 3–5, scTranslator notebooks |
| `proxy_models/` | produced by notebook 7 | notebook 8 |

Cohort folder names (exact):

| Platform | Folders |
|---|---|
| TMT | colon_2019, renal_2019, Endomet_val, breast2020, pedbrain_2020, lung_2020_smoke, lung_val, lung2021, headneck_2021, pancrea2021, brain_2021, liver_val |
| iTRAQ | breast_2016, ovary_2016, stomach_2019 |
| LabelFree | brain, colon, prostate, liver, lung |

**Order:** 1 → 2 → 3 → 4 → 5 → 6 → 7 → 8 → 9 → 10

**Optional: scTranslator notebooks** (`scTranslator-validation-newdata/`). They need notebook 2 to have run, plus the scTranslator code and checkpoints. Place them in `<repo>/external/scTranslator/scTranslator/` and `<repo>/external/scTranslator/checkpoints/`, or point the environment variables `SCTRANSLATOR_REPO` and `SCTRANSLATOR_CKPT_DIR` at the copies in `02. elasticNet/data/`.

---

## Copy summary (what moves between stages)

| File / folder | From | To |
|---|---|---|
| `transcriptomics_filtered.joblib`, `proteomics_matched.joblib` | 01 | 02, 03, 04, 05 |
| `correlations.joblib` | 01 | 02, 03 |
| `protein_list.txt` | you | 02, 03 |
| `results_enet/` | 02 | 04 |
| `res_enet_zbig/` | 02 | 04, 05 |
| `loo/` | 02 | 05 |
| `sc_predictions_*.npz` | 02 | 04 |
| `results_xgb/` | 03 | 04 |
| `res_xgb_zbig/` | 03 | 04, 05 |
| `xgb_all_results.h5` | 03 | 02 |

## Python environment
Use `numpy==1.26.4`, `tables==3.9.2`, `numexpr==2.8.7`, `h5py==3.11.0`, plus `xgboost`, `pyarrow` and `openpyxl`. For the scTranslator notebooks you also need `torch`, `einops` and `local-attention==1.4.3`.
The first cell of each fixed notebook checks the environment and prints the exact fix if something is broken.
