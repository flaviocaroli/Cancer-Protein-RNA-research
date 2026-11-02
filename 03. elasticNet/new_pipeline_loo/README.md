# new_pipeline_loo – ElasticNet leave-one-out (Suppl. S4.2)

Reads and writes **`../data/`** (the single `03. elasticNet/data` folder).
Inputs: `transcriptomics_filtered.joblib` (default `--rna_file`), `proteomics_matched.joblib` and `protein_list.txt` (your list).

## Run (from this folder), one submission per normalization
```bash
for N in zscore none ranking binning; do NORMALIZATION_TYPE=$N bash new_run_enet.sh; done
```
| Variable | Default |
|---|---|
| `NORMALIZATION_TYPE` | zscore |
| `RESULTS_DIR` | `../data/loo/results_<norm>` |
| `FULL_OUTPUT` | 1 (also writes `<P>_LOO.parquet`) |
| `RNA_FILE` | transcriptomics_filtered.joblib |
| `PIPELINE_TYPE`, `PERCENTILE`, `N_BINS`, `TARGET_CHUNKS`, `CONDA_ENV` | fixed_percentile, 10, 50, 10, test_env |

Local: `python src/new_pipeline.py --normalization zscore --full_output --group_size 50`.

## Outputs (`../data/loo/results_<norm>/`)
| File | Content | Used by |
|---|---|---|
| `<P>_best.parquet` | Best LOO split: coefficients, params, z-score parameters | `06. validation` notebooks 3–5 (copy the four folders to `06. validation/data/loo/`) |
| `<P>_LOO.parquet` | One row per LOO sample | `model_results/00_build_enet_all_results` |
| `batch_<s>_<e>_fixed_percentile_<norm>.parquet` | Chunk summary | `06` notebook 3 |

`new_smoke.ipynb` tests 2 proteins and writes to `../data/loo/results_smoke/`.
