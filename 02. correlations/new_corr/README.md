# new_corr – chunked correlation pipeline

Correlates each protein with every RNA gene (`--target rna`) or with every other protein (`--target protein`), using Pearson or Spearman. It then assembles the per-protein results into one matrix.

## Inputs (`data/`)
| File | Source |
|---|---|
| `transcriptomics_filtered.joblib` | 01 output (samples × genes) |
| `proteomics_matched.joblib` | 01 output (samples × proteins) |
| `protein_list.txt` | in repo (5,113 proteins, one per line) |

## Significance rule
| Mode | Kept value | Otherwise |
|---|---|---|
| Bootstrap (default, 100 draws over the protein's non-missing samples) | Bootstrap median, if the 95% percentile CI excludes 0 and there are at least max(10, n/10) usable draws | 0 |
| `--no_bootstrap` | r, if the parametric p-value is below `--p_threshold` (default 0.05) | 0 |

The computation is vectorised over features. Parametric results equal `scipy.stats.pearsonr` / `spearmanr`. Bootstrap medians follow the same rule as the original per-gene loop, but the random streams differ, so they are not bit-identical.

## Outputs
| Run | Folder (default) | Per chunk | Assembled |
|---|---|---|---|
| `--method pearson --target rna` | `data/results/` | `per_protein/<P>_corr.parquet`, `batch_<s>_<e>_corr_summary.parquet` | `full_pearson100.joblib` (proteins × genes) |
| `--method spearman --target rna` | `data/results_spearman/` | same | `full_spearman100.joblib` |
| `--method pearson --target protein` | `data/results_prot/` | same | `full_prot.joblib` (square proteins × proteins, diagonal 1) |

Per-protein files are 1 row × features, and non-significant cells are 0. Without a parquet engine the files fall back to CSV.

## Run on SLURM (from this folder)
```bash
bash new_run_corr.sh                       # Pearson vs RNA  -> full_pearson100
METHOD=spearman bash new_run_corr.sh       # Spearman vs RNA -> full_spearman100
TARGET=protein bash new_run_corr.sh        # protein-protein -> full_prot
```
Optional environment variables: `TARGET_CHUNKS` (25), `N_BOOTSTRAPS` (100), `BOOTSTRAP` (1), `P_THRESHOLD` (0.05), `CONDA_ENV` (test_env), `RESULTS_DIR`.
Each call submits the chunk array plus an assembly job that runs after the array (`--dependency=afterok`).

## Run locally
```bash
python corr_pipeline.py --chunk_id 0 --group_size 500             # one chunk
python corr_pipeline.py --method spearman --chunk_id 0 --group_size 500
python corr_pipeline.py --assemble                                # after all chunks
python corr_pipeline.py --method spearman --assemble
python corr_pipeline.py --target protein --assemble
```
`--assemble` warns about proteins in the list that have no per-protein file. Re-run those chunks before you use the matrix.

`smoke_test.ipynb` runs a synthetic test without writing to the real results folders.
