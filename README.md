# Cancer-Protein-RNA-research

Code for **"Protein Abundance Prediction from Bulk RNA-seq"** (Caroli, Tangherloni, Mazzara, Buffa).

The workflow predicts each protein from bulk RNA-seq with models of increasing complexity and assigns the protein to the simplest one that is good enough:

```text
same-gene proxy -> ElasticNet -> XGBoost (-> scTranslator, reference only) -> unresolved
```

CCLE: 369 cell lines, 22,614 genes, 5,113 proteins. The 80/20 split gives 295 development samples and 74 held-out samples.
External transfer: 20 cohorts (12 TMT, 3 iTRAQ, 5 label-free), 2,560 samples, scored without refitting.

## Repository layout

| Folder | Paper section | What it does |
|---|---|---|
| `01. pre-processing/` | Methods – CCLE preprocessing | Download, match and KDE-filter CCLE. Writes `transcriptomics_filtered.joblib`, `proteomics_matched.joblib`, `correlations.joblib`. |
| `02. elasticNet/` | Methods – ElasticNet; Suppl. S1.1, S4.2 | 80/20 pipeline (primary), LOO pipeline (robustness), LOO analysis, scTranslator CCLE predictions. |
| `03. xgb/` | Methods – XGBoost; Suppl. S1.2, S4.2 | 80/20 pipeline (primary), LOO pipeline, LOO merge and analysis. |
| `04. ccle_res_analysis/` | Results – Table 1, Figs 2–4, S2–S5 | Matched benchmark, same-gene proxies, cascade, HPA and pan-atlas analysis. |
| `05. validation/` | Results – external transfer, Figs 5–7, S6–S9 | 20-cohort external validation, Xu et al. comparison, scTranslator benchmarks. |
| `external/` (git-ignored) | — | scTranslator code and checkpoints. Copy from `research_data/00_external`. |

Every folder has its own README listing its inputs, outputs and run order.
`research_data/` (shared separately, git-ignored) holds all inputs and outputs. Its paths mirror the repo, and `research_data/MANIFEST.xlsx` lists every file.

## Data sources

| Data | Source |
|---|---|
| CCLE RNA-seq, DepMap 20Q2, log2(TPM+1) | https://figshare.com/ndownloader/files/25797014 |
| CCLE sample info | https://figshare.com/ndownloader/files/25494443 |
| CCLE TMT proteomics (Nusinow et al. 2020) | https://gygi.hms.harvard.edu/data/ccle/Table_S2_Protein_Quant_Normalized.xlsx |
| 20 external cohorts ("Collection A", Xu et al. 2022) | EMBL-EBI BioStudies **S-BSST733** |
| Xu et al. cohort statistics | `Table0_S1-S8_Combined20220114.xlsx` (Briefings in Bioinformatics 23(3):bbac091, supplementary) |
| HPA, UniProt, GO, Reactome, WikiPathways | Downloaded automatically by `04. ccle_res_analysis/5_hpa_pan_atlas.ipynb` (Part A) |
| scTranslator code and checkpoints (2M / 160K / 10K) | Liu et al., Nature Biomedical Engineering 2026 |

## Installation

```bash
conda create -n prot_rna python=3.10 -y
conda activate prot_rna
pip install -r requirements.txt
# optional, scTranslator notebooks only:
pip install torch
```

## End-to-end order

The order follows the data dependencies, which is not always the folder number.

| Step | Run | Produces |
|---|---|---|
| 1 | `01. pre-processing` notebooks 01 → 02 (→ 03 optional) | Matched and filtered CCLE matrices, same-gene correlations |
| 2 | `02. elasticNet/new_pipeline_8020` and `03. xgb/new_pipeline_xgb` (k=2000 all normalizations; z-score k=5000) | `results_enet`, `res_enet_zbig`, `results_xgb`, `res_xgb_zbig` |
| 3 | `02. elasticNet/model_results/03` (optional, GPU) | `sc_predictions_{2M,160K,10K}.npz` |
| 4 | `04. ccle_res_analysis` notebooks 0 → 1 → 2 → 3 → 4 → 5_hpa_pan_atlas | Table 1, Figs 2–4, S2–S5 |
| 5 | `03. xgb/pipeline_xgb_loo` → `03. xgb/model_results` 01 → 02 | XGB LOO (Fig. S10b), `xgb_all_results.h5` |
| 6 | `02. elasticNet/new_pipeline_loo` (4 normalizations) → `02. elasticNet/model_results` 00 → 01 → 02 | ENet LOO (Fig. S10a, S4) |
| 7 | `05. validation/notebooks` 1 → 2 → 3 → 4 → 5 | External matrices and exploratory LOO-model transfer |
| 8 | `05. validation/notebooks` 6 → 7 → 8 → 9 → 10 | Paper external validation (Figs 5–7, S6–S9, Table 2) |
| 9 | `05. validation/scTranslator-validation-newdata` (optional) | scTranslator on external platforms |

## Data: one store, linked into the repo
All data lives once in `research_data/`, a separate folder, git-ignored, shared outside git. It is organised by the stage that **produces** each dataset, with your own files in `00_external`.
`link_research_data.py` (repo root) links every stage `data/` entry and output folder to that store. A dataset used by several stages, such as `res_enet_zbig/` used by 02, 04 and 05, is therefore stored once, and runs write straight into the store.

```bash
python link_research_data.py --store ../research_data                 # dry run: what will be linked
python link_research_data.py --store ../research_data --migrate --apply   # move existing repo data into the store and link
python link_research_data.py --store ../research_data --check         # status of every link
```
The full list of datasets and links is in `research_data/MANIFEST.xlsx`. The step-by-step setup is in `research_data/GUIDE.md`.

## Figures
Every notebook displays its figures inline and saves them to `<notebook folder>/figures_pdf/<notebook name>/`. This is handled by the first code cell (*pRNA figure routing*), which redirects every `savefig` call, including calls inside `src/` helpers.

## Reproducibility notes

- Fixed 80/20 split: 295/74. Hyperparameters are chosen by 3-fold CV inside the 295 samples, selecting on mean validation MSE.
- Z-score parameters and feature rankings are fitted on training folds only. Rank and QBin50 transforms are per-sample, so they fit no parameters.
- ElasticNet (paper S1.1): `fit_intercept=False`, `max_iter=10000`, `tol=1e-3`, `selection=random`, at most 8 random grid combinations, `fast` grid.
- XGBoost (paper S1.2): `hist`, `reg:squarederror`, early stopping 10, final n_estimators = max(25, min(n, ⌈1.1·median best iteration⌉)).
- Pipelines skip existing result files, so interrupted runs can be resumed. Use a new results folder whenever you change settings.
- Record the git commit, package versions, seeds and split sample IDs for every run.

## License
MIT – see `LICENSE`.
