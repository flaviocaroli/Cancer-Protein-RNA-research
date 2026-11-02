# 03. elasticNet

Code is in the subfolders. **All data lives in one folder: `data/`.**

## `data/` (the only data folder)
| Path | What | Source / producer |
|---|---|---|
| `transcriptomics_filtered.joblib`, `proteomics_matched.joblib`, `correlations.joblib` | CCLE inputs | copy from `01. pre-processing/data/` |
| `protein_list.txt` | Protein list for the LOO run | your local list |
| `protein_clusters_k3.csv` | Clusters | copy from `02. correlations/results/clusters_pg/` |
| `xgb_all_results.h5` | Merged XGB LOO | copy from `04. xgb/data/` |
| `results_enet/` | 80/20 k=2000 (none, zscore, ranking, binning50) | `new_pipeline_8020` |
| `res_enet_zbig/` | 80/20 z-score k=5000 (**primary model**) | `new_pipeline_8020` |
| `loo/results_{zscore,none,ranking,binning}/` | LOO `<P>_best.parquet`, `<P>_LOO.parquet`, `batch_*` | `new_pipeline_loo` |
| `loo/logs/` | LOO SLURM logs | `new_pipeline_loo` |
| `enet_all_results.h5`, `enet_best_summary.csv` | Merged LOO | `model_results/00` |
| `sc_predictions_{2M,160K,10K}.npz` | scTranslator on CCLE | `model_results/03` |

## Order
`new_pipeline_8020` → (`model_results/03`) → `new_pipeline_loo` → `model_results/00 → 01 → 02`.
Install with `pip install -r requirements.txt`. Each subfolder README explains how to run it.
