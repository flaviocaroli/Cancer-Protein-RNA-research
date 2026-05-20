# new_pipeline_enet

Fast 80/20 ElasticNet pipeline for protein prediction from transcriptomics.

Expected layout:

```text
new_pipeline_enet/
  run_enet.sh
  enet_8020_smoke_test.ipynb
  src/
    functions.py
    pipeline.py
  data/
    transcriptomics_filtered.joblib
    proteomics_matched.joblib
  results_enet/
  logs/
```

The pipeline auto-detects the proteomics orientation. For your current data it expects:

```text
transcriptomics_filtered.joblib: rows = samples, columns = genes
proteomics_matched.joblib:      rows = samples, columns = proteins
```

Each protein-normalization run uses an 80/20 internal train/test split, tunes ElasticNet on the train split using inner CV, evaluates on the test split, then trains a final model on all Dataset 1 samples and saves it to H5.

Saved H5 contains:

```text
protein
train_samples
test_samples
test_mse
test_pearson_corr
test_cosine_similarity
best_params
selected_features
final all-samples ElasticNet model artifact
```

It intentionally does not store test true/predicted values or RMSE.

## Small cluster test

```bash
NORMALIZATIONS=none,binning \
GRID_PRESET=tiny \
MAX_PARAM_COMBOS=2 \
SELECT_K=50 \
PROTEIN_START=1 \
PROTEIN_END=5 \
PROTEINS_PER_TASK=5 \
PARTITION=medium_cpu \
TIME=03:00:00 \
bash run_enet.sh
```

## Full run

```bash
NORMALIZATIONS=none,zscore,ranking,binning \
GRID_PRESET=fast \
MAX_PARAM_COMBOS=8 \
SELECT_K=2000 \
N_BINS=50 \
PROTEINS_PER_TASK=1000 \
PARTITION=compute \
TIME=2-00:00:00 \
bash run_enet.sh
```

## Monitor

```bash
squeue -u $USER
sacct -j JOBID --format=JobID,JobName%20,State,Elapsed,Start,End,ExitCode
find results_enet -name "*.h5" | wc -l
```
