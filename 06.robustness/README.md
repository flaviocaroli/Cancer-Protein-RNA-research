# 06.robustness_repeated_8020

Repeated 80/20 robustness pipeline for the CCLE RNA -> protein benchmark.

This version is designed for clusters with low SLURM job limits. By default it uses
`TASK_MODE=packed_by_protein`, which creates one SLURM array task per protein chunk.
With `PROTEIN_CHUNKS=20`, only 20 jobs are submitted once.

## Expected data

Put the preprocessed matrices in `data/`:

```text
transcriptomics_filtered.joblib   # samples x genes
proteomics_matched.joblib         # samples x proteins
```

or pass explicit paths with `X_PATH=/path/to/X.joblib` and `Y_PATH=/path/to/Y.joblib`.

For the main CCLE benchmark you should see approximately:

```text
X = 369 samples x ~22,000 genes
Y = 369 samples x ~5,113 proteins
```

If you see `4843 samples x 369 genes`, your input matrices are wrong or transposed before this pipeline.

## Smoke test

```bash
bash run_smoke_test.sh
```

## First real 20-split ENet run on HPC

```bash
rm -rf results

DATA_DIR=$PWD/data \
SUBMIT=slurm \
PARTITION=compute \
TIME=3-00:00:00 \
MEM=24G \
CPUS_PER_TASK=1 \
N_JOBS=1 \
TASK_MODE=packed_by_protein \
PROTEIN_CHUNKS=20 \
ARRAY_THROTTLE=20 \
MODELS=proxy,enet \
NORMALIZATIONS=zscore \
SELECT_KS=2000 \
N_REPEATS=20 \
CV_FOLDS=5 \
GRID_PRESET=fast \
MAX_PARAM_COMBOS=8 \
KEEP_PREDICTIONS=0 \
bash run_repeated.sh
```

This submits exactly 20 array tasks by default. Each task handles one protein chunk and loops over all splits/configs for that chunk.

After all jobs finish:

```bash
DATA_DIR=$PWD/data \
PROXY_CONFIG=proxy:zscore:0 \
ENET_CONFIG=enet:zscore:2000 \
METRIC=pearson \
THRESHOLD=0.50 \
PASS_RATE=0.80 \
GAIN_THRESHOLD=0.05 \
N_BOOT=2000 \
bash run_after_jobs.sh
```

Main outputs:

```text
results/metrics/metrics_long.parquet      # or metrics_long.csv.gz if pyarrow is missing
results/cascade/robust_cascade_atlas.tsv
results/cascade/robust_cascade_counts.tsv
results/cascade/paired_gain_tests.tsv
results/cascade/bootstrap_metric_ci.tsv
```

## Final 50-split run

```bash
rm -rf results_final_enet

DATA_DIR=$PWD/data \
RESULTS_DIR=$PWD/results_final_enet \
SUBMIT=slurm \
PARTITION=compute \
TIME=3-00:00:00 \
MEM=24G \
CPUS_PER_TASK=1 \
N_JOBS=1 \
TASK_MODE=packed_by_protein \
PROTEIN_CHUNKS=20 \
ARRAY_THROTTLE=20 \
MODELS=proxy,enet \
NORMALIZATIONS=zscore \
SELECT_KS=2000 \
N_REPEATS=50 \
CV_FOLDS=5 \
GRID_PRESET=fast \
MAX_PARAM_COMBOS=8 \
KEEP_PREDICTIONS=0 \
bash run_repeated.sh
```

Then:

```bash
DATA_DIR=$PWD/data \
RESULTS_DIR=$PWD/results_final_enet \
PROXY_CONFIG=proxy:zscore:0 \
ENET_CONFIG=enet:zscore:2000 \
METRIC=pearson \
THRESHOLD=0.50 \
PASS_RATE=0.80 \
GAIN_THRESHOLD=0.05 \
N_BOOT=2000 \
bash run_after_jobs.sh
```

## Notes

- `KEEP_PREDICTIONS=0` keeps storage light. The robust atlas is built from metric shards.
- Install `pyarrow` if possible. Without it, the code writes `.csv.gz` fallbacks.
- Use `TASK_MODE=grid` only on clusters without strict job-count limits.
