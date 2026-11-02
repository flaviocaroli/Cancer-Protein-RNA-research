# 02. correlations

Protein-level correlation structure (Supplementary workflow stage 2). Code is in `new_corr/` and the notebooks. **All data lives in one folder: `data/`.**

## `data/` (the only data folder)
| Path | What | Source |
|---|---|---|
| `transcriptomics_filtered.joblib`, `proteomics_matched.joblib` | CCLE matrices | copy from `01. pre-processing/data/` |
| `protein_list.txt` | 5,113 proteins | in repo |
| `corr_results/<method>_<target>/per_protein/*.parquet`, `batch_*` | Chunk outputs | `new_corr` |
| `corr_results/logs/` | SLURM logs | `new_corr/new_run_corr.sh` |
| `full_pearson100.joblib` | Protein × gene Pearson | `new_corr` assembly |
| `full_spearman100.joblib` | Protein × gene Spearman (used by 05) | `new_corr` assembly |
| `full_prot.joblib` | Protein × protein Pearson | `new_corr` assembly |

## Step 1 – compute the matrices (from `new_corr/`)
```bash
bash new_run_corr.sh                     # -> data/full_pearson100.joblib
METHOD=spearman bash new_run_corr.sh     # -> data/full_spearman100.joblib
TARGET=protein bash new_run_corr.sh      # -> data/full_prot.joblib
```
Each submission runs a chunk array plus an automatic assembly job. The assembled matrix is written straight into `data/`, so there is no copy step.

## Step 2 – notebooks (run in order, from this folder)
| # | Notebook | Reads (`data/`) | Writes |
|---|---|---|---|
| 1 | `01. pearson_anlaysis.ipynb` | `full_pearson100.joblib`, `transcriptomics_filtered.joblib` | `pearson proteins/protein_list_0{5,6,7,8}.txt`, `results/clusters_pg/*` |
| 2 | `02. prot_clustering.ipynb` | `full_prot.joblib`, `proteomics_matched.joblib` | `results/clusters_pp/*`, `results/clusters_pe/*` |
| 3 | `03. cluster_anlaysis.ipynb` | `full_pearson100.joblib`, `transcriptomics_filtered.joblib`, `results/clusters_*` | `results/overlap_analysis/*` |

`results/clusters_pg/protein_clusters_k3.csv` is copied to `03. elasticNet/data/`.

## Figures
Every notebook shows its figures inline and saves them to `figures_pdf/<notebook name>/` next to the notebook. The first code cell (*pRNA figure routing*) handles this; keep it.
