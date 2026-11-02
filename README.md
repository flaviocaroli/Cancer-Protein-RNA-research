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
| `01. pre-processing/` | Methods – CCLE preprocessing | Download, match and KDE-filter CCLE. Writes `transcriptomics_filtered.joblib` and `proteomics_matched.joblib`. |
| `02. correlations/` | Supplementary workflow stage 2 | Chunked protein–RNA and protein–protein correlations, plus clustering notebooks. |
| `03. elasticNet/` | Methods – ElasticNet; Suppl. S1.1, S4.2 | 80/20 pipeline (primary), LOO pipeline (robustness), LOO analysis and scTranslator CCLE predictions. |
| `04. xgb/` | Methods – XGBoost; Suppl. S1.2, S4.2 | 80/20 pipeline (primary), LOO pipeline, LOO merge and analysis. |
| `05. ccle_res_analysis/` | Results – Table 1, Figs 2–4, S2–S5 | Matched benchmark, same-gene proxies, cascade, HPA and pan-atlas analysis. |
| `06. validation/` | Results – external transfer, Figs 5–7, S6–S9 | 20-cohort external validation, Xu et al. comparison, scTranslator benchmarks. |
| `external/` (git-ignored) | — | scTranslator code and checkpoints. Copy from `research_data/00_external`. |

Every folder has its own README listing its inputs, outputs and run order.
`research_data/` (shared separately, git-ignored) holds all inputs and outputs. Its paths mirror the repo, and `research_data/MANIFEST.csv` lists every file.

## Data sources

| Data | Source |
|---|---|
| CCLE RNA-seq, DepMap 20Q2, log2(TPM+1) | https://figshare.com/ndownloader/files/25797014 |
| CCLE sample info | https://figshare.com/ndownloader/files/25494443 |
| CCLE TMT proteomics (Nusinow et al. 2020) | https://gygi.hms.harvard.edu/data/ccle/Table_S2_Protein_Quant_Normalized.xlsx |
| 20 external cohorts ("Collection A", Xu et al. 2022) | EMBL-EBI BioStudies **S-BSST733** |
| Xu et al. cohort statistics | `Table0_S1-S8_Combined20220114.xlsx` (Briefings in Bioinformatics 23(3):bbac091, supplementary) |
| HPA, UniProt, GO, Reactome, WikiPathways | Downloaded automatically by `05. ccle_res_analysis/5_hpa_pan_atlas.ipynb` (Part A) |
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
| 1 | `01` notebooks 01 → 02 (→ 03 optional) | Matched and filtered CCLE matrices |
| 2 | `02/new_corr`: Pearson-RNA, Spearman-RNA and Pearson-protein runs, then `02` notebooks 01 → 02 → 03 | `full_pearson100`, `full_spearman100`, `full_prot`, clusters |
| 3 | `03/new_pipeline_8020` and `04/new_pipeline_xgb` (k=2000 all normalizations; z-score k=5000) | `results_enet`, `res_enet_zbig`, `results_xgb`, `res_xgb_zbig` |
| 4 | `03/model_results/03` (optional, GPU) | `sc_predictions_{2M,160K,10K}.npz` |
| 5 | `05` notebooks 0 → 1 → 2 → 3 → 4 → 5_hpa_pan_atlas | Table 1, Figs 2–4, S2–S5 |
| 6 | `03/new_pipeline_loo` (4 normalizations) → `03/model_results` 00 → 01 → 02 | LOO robustness (Fig. S10) |
| 7 | `04/pipeline_xgb_loo` → `04/model_results` 01 → 02 | XGB LOO (Fig. S10b) |
| 8 | `06/notebooks` 1 → 2 → 3 → 4 → 5 | External matrices and exploratory LOO-model transfer |
| 9 | `06/notebooks` 6 → 7 → 8 → 9 → 10 | Paper external validation (Figs 5–7, S6–S9, Table 2) |
| 10 | `06/scTranslator-validation-newdata` (optional) | scTranslator on external platforms |

## Data folders
Each stage folder has exactly **one `data/` folder** holding all of its inputs and shared intermediate files. The code in its subfolders reads `../data/` (or `data/` from the stage root). Files that one stage passes to another (for example `transcriptomics_filtered.joblib`) are copied or symlinked between stage `data/` folders. `research_data/` mirrors this layout.

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
