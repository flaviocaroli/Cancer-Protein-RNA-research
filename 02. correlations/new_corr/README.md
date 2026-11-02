# new_corr – chunked correlation pipeline

Correlates each protein with every RNA gene (`--target rna`) or with every protein (`--target protein`), using Pearson or Spearman, then assembles one matrix.
Code only: inputs and outputs are in **`../data/`** (the single data folder of `02. correlations`).

## Inputs (`../data/`)
`transcriptomics_filtered.joblib`, `proteomics_matched.joblib` (from 01) and `protein_list.txt` (in repo).

## Significance rule
| Mode | Kept value | Otherwise |
|---|---|---|
| Bootstrap (default, 100 draws over the protein's non-missing samples) | Bootstrap median, if the 95% percentile CI excludes 0 and there are at least max(10, n/10) usable draws | 0 |
| `--no_bootstrap` | r, if the parametric p-value is below `--p_threshold` | 0 |

Parametric results equal `scipy.stats.pearsonr` / `spearmanr`.

## Outputs (`../data/`)
| Run | Chunk files | Assembled matrix |
|---|---|---|
| `--method pearson --target rna` (default) | `corr_results/pearson_rna/` | `full_pearson100.joblib` (proteins × genes) |
| `--method spearman --target rna` | `corr_results/spearman_rna/` | `full_spearman100.joblib` |
| `--method pearson --target protein` | `corr_results/pearson_protein/` | `full_prot.joblib` (square, diagonal 1) |

## Run
```bash
bash new_run_corr.sh                     # SLURM: chunks + automatic assembly
METHOD=spearman bash new_run_corr.sh
TARGET=protein bash new_run_corr.sh
# locally
python corr_pipeline.py --chunk_id 0 --group_size 500
python corr_pipeline.py --assemble
```
Optional environment variables: `TARGET_CHUNKS` (25), `N_BOOTSTRAPS` (100), `BOOTSTRAP` (1), `P_THRESHOLD` (0.05), `CONDA_ENV`, `RESULTS_DIR`.
`smoke_test.ipynb` runs a synthetic test.
