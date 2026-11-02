# new_pipeline_loo – ElasticNet leave-one-out (robustness, Suppl. S4.2)

For each protein and each LOO split, the pipeline uses inner 5-fold CV and ANOVA percentile feature selection on the training samples, then predicts the held-out sample.

## Inputs (`data/`)
| File | Source |
|---|---|
| `transcriptomics_filtered.joblib` | 01 output. This is the default `--rna_file`, the same features as the 80/20 and XGB pipelines. |
| `proteomics_matched.joblib` | 01 output |
| `protein_list.txt` | **Your local list** (in `research_data/03_elasticNet/input/new_pipeline_loo/data/protein_list.txt`) |

Optional sensitivity run on genes shared with all external sets: put `transcriptomics_common.joblib` (from `06/notebooks/2`) in `data/` and pass `RNA_FILE=transcriptomics_common.joblib`.

## Run (from this folder) – one submission per normalization
```bash
for N in zscore none ranking binning; do NORMALIZATION_TYPE=$N bash new_run_enet.sh; done
```
| Variable | Default | Meaning |
|---|---|---|
| `NORMALIZATION_TYPE` | zscore | zscore, none, ranking or binning |
| `RESULTS_DIR` | `data/results_<norm>` | One folder per normalization, so runs no longer overwrite each other |
| `FULL_OUTPUT` | 1 | Also write per-sample `<P>_LOO.parquet`. Needed for `model_results/00`. |
| `RNA_FILE` | transcriptomics_filtered.joblib | RNA matrix in `data/` |
| `PIPELINE_TYPE`, `PERCENTILE`, `N_BINS`, `TARGET_CHUNKS`, `CONDA_ENV`, `BASE_DIR` | fixed_percentile, 10, 50, 10, test_env, current folder | |

Local: `python src/new_pipeline.py --normalization zscore --full_output --group_size 50 --results_dir data/results_zscore`

## Outputs (`data/results_<norm>/`)
| File | Content | Used by |
|---|---|---|
| `<P>_best.parquet` | Best LOO split: coefficients (`feat_coeffs_dict`), params, `mse_all_others`, and z-score parameters for the z-score run | `06/notebooks` 3, 4, 5 |
| `<P>_LOO.parquet` | One row per LOO sample: `protein`, `LOO sample`, `true`, `predicted`, `params_dict`, `features_coef_dict` | `model_results/00_build_enet_all_results` |
| `batch_<s>_<e>_fixed_percentile_<norm>.parquet` | Chunk summary (dict columns stored as JSON) | `06/notebooks/3` |

Where the folders go next:

| Folder | Copy to | Also used by |
|---|---|---|
| `results_zscore` | `06. validation/results_val/` | `model_results/00` |
| `results_none` | `06. validation/results_none/` | |
| `results_ranking` | `06. validation/results_rk/` | |
| `results_binning` | `06. validation/results_bin/` | |

`new_smoke.ipynb` runs 2 proteins in-process and through the CLI, writing to `data/results_smoke/`.
