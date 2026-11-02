# 02. correlations

Protein-level correlation structure (paper, Supplementary workflow stage 2). The `new_corr/` pipeline computes the full matrices, and the three notebooks analyse them.

## Step 1 – compute matrices (`new_corr/`, see its README)
| Run | Command (from `new_corr/`) | Final file |
|---|---|---|
| Protein × RNA, Pearson, 100 bootstraps | `bash new_run_corr.sh` | `data/results/full_pearson100.joblib` |
| Protein × RNA, Spearman, 100 bootstraps | `METHOD=spearman bash new_run_corr.sh` | `data/results_spearman/full_spearman100.joblib` |
| Protein × protein, Pearson | `TARGET=protein bash new_run_corr.sh` | `data/results_prot/full_prot.joblib` |

Each submission runs a SLURM array of chunks plus an automatic assembly job.

## Step 2 – copy into `data/`
```text
data/full_pearson100.joblib        <- new_corr/data/results/full_pearson100.joblib
data/full_prot.joblib              <- new_corr/data/results_prot/full_prot.joblib
data/transcriptomics_filtered.joblib, data/proteomics_matched.joblib   <- 01 outputs
```
`full_spearman100.joblib` is not used here. It goes to `05. ccle_res_analysis/data/`.

## Step 3 – notebooks (run in order)
| # | Notebook | Reads | Writes |
|---|---|---|---|
| 1 | `01. pearson_anlaysis.ipynb` | `data/full_pearson100.joblib`, `data/transcriptomics_filtered.joblib` | `pearson proteins/protein_list_0{5,6,7,8}.txt`, `results/clustermaps/map_pg_std.png`, `results/clusters_pg/*` (protein × gene clusters, communities, network tables) |
| 2 | `02. prot_clustering.ipynb` | `data/full_prot.joblib`, `data/proteomics_matched.joblib` | `results/clustermaps/map_pp_raw.png`, `map_pe_raw.png`, `results/clusters_pp/*`, `results/clusters_pe/*` |
| 3 | `03. cluster_anlaysis.ipynb` | `data/full_pearson100.joblib`, `data/transcriptomics_filtered.joblib`, `results/clusters_{pg,pp,pe}/` | `results/overlap_analysis/*` |

`results/clusters_pg/protein_clusters_k3.csv` is also used by `03. elasticNet/model_results/01`.
Helpers: `src/corr_utils.py`.
