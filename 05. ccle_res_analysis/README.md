# CCLE Protein Prediction Runs

This folder contains model outputs from four 80/20 train-test pipelines used to predict protein abundance from transcriptomics.

## Data used

All runs use the same input files:

```text
data/transcriptomics_filtered.joblib
data/proteomics_matched.joblib
```

Expected orientation:

```text
transcriptomics_filtered.joblib:
  rows    = samples
  columns = genes / RNA features

proteomics_matched.joblib:
  rows    = samples
  columns = proteins
```

Each pipeline performs:

```text
1. Random 80/20 train-test split of samples
2. Feature selection on the training set
3. Hyperparameter tuning using inner CV on the 80% training set
4. Evaluation on the 20% held-out test set
5. Saving one compressed H5 file per protein per run
6. Saving the final model artifact for future prediction
```

The main metrics saved per protein are:

```text
test_mse
test_pearson_corr
test_cosine_similarity
selected_features
best_params
train_samples
test_samples
model artifact
```

---

# Run 1: XGBoost baseline, K=2000

## Output folder

```text
results_xgb/
```

## Model

```text
XGBoost regressor
```

## Main settings

```text
Pipeline type:       80/20 train-test
Normalizations:      none, zscore, rank, qbin50
Feature selection:   top 2000 genes
SELECT_K:            2000
GRID_PRESET:         fast
MAX_PARAM_COMBOS:    8
SPLITS:              3
EARLY_STOPPING:      10
TUNING_STRATEGY:     once
Output mode:         minimal H5
```

## Typical command

```bash
NORMALIZATIONS=none,zscore,rank,qbin50 \
GRID_PRESET=fast \
MAX_PARAM_COMBOS=8 \
SELECT_K=2000 \
EARLY_STOPPING=10 \
SPLITS=3 \
TUNING_STRATEGY=once \
PROTEINS_PER_TASK=1000 \
PARTITION=compute \
TIME=2-00:00:00 \
bash run_xgboost.sh
```

## Expected output files

```text
results_xgb/<PROTEIN>_xgb_8020_none_fast_k2000.h5
results_xgb/<PROTEIN>_xgb_8020_zscore_fast_k2000.h5
results_xgb/<PROTEIN>_xgb_8020_rank_fast_k2000.h5
results_xgb/<PROTEIN>_xgb_8020_qbin50_fast_k2000.h5
```

---

# Run 2: XGBoost z-score big, K=5000, 16 params

## Output folder

```text
res_xgb_zbig/
```

## Model

```text
XGBoost regressor
```

## Main settings

```text
Pipeline type:       80/20 train-test
Normalization:       zscore only
Feature selection:   top 5000 genes
SELECT_K:            5000
GRID_PRESET:         fast
MAX_PARAM_COMBOS:    16
SPLITS:              5
EARLY_STOPPING:      20
TUNING_STRATEGY:     once
Output mode:         minimal H5
```

## Typical command

```bash
NORMALIZATIONS=zscore \
GRID_PRESET=fast \
MAX_PARAM_COMBOS=16 \
SELECT_K=5000 \
EARLY_STOPPING=20 \
SPLITS=5 \
TUNING_STRATEGY=once \
PROTEINS_PER_TASK=500 \
PARTITION=compute \
TIME=2-00:00:00 \
RESULTS_DIR=res_xgb_zbig \
bash run_xgboost.sh
```

## Expected output files

```text
res_xgb_zbig/<PROTEIN>_xgb_8020_zscore_fast_k5000.h5
```

---

# Run 3: ElasticNet baseline, K=2000

## Output folder

```text
results_enet/
```

## Model

```text
ElasticNet regression
```

## Main settings

```text
Pipeline type:       80/20 train-test
Normalizations:      none, zscore, ranking, binning
Feature selection:   top 2000 genes
SELECT_K:            2000
GRID_PRESET:         fast
MAX_PARAM_COMBOS:    16
N_BINS:              10 or 50 depending on the run
FIT_INTERCEPT:       1
Output mode:         minimal H5
```

## Typical command

```bash
CONDA_ENV=xgb_loo \
NORMALIZATIONS=none,zscore,ranking,binning \
N_BINS=10 \
GRID_PRESET=fast \
MAX_PARAM_COMBOS=16 \
SELECT_K=2000 \
FIT_INTERCEPT=1 \
PROTEINS_PER_TASK=1000 \
PARTITION=compute \
TIME=2-00:00:00 \
bash run_enet.sh
```

## Expected output files

```text
results_enet/<PROTEIN>_enet_8020_none_fast_k2000.h5
results_enet/<PROTEIN>_enet_8020_zscore_fast_k2000.h5
results_enet/<PROTEIN>_enet_8020_ranking_fast_k2000.h5
results_enet/<PROTEIN>_enet_8020_binning_fast_k2000.h5
```

---

# Run 4: ElasticNet z-score big, K=5000, 16 params

## Output folder

```text
res_enet_zbig/
```

## Model

```text
ElasticNet regression
```

## Main settings

```text
Pipeline type:       80/20 train-test
Normalization:       zscore only
Feature selection:   top 5000 genes
SELECT_K:            5000
GRID_PRESET:         fast
MAX_PARAM_COMBOS:    16
FIT_INTERCEPT:       1
Output mode:         minimal H5
```

## Typical command

```bash
CONDA_ENV=xgb_loo \
NORMALIZATIONS=zscore \
GRID_PRESET=fast \
MAX_PARAM_COMBOS=16 \
SELECT_K=5000 \
FIT_INTERCEPT=1 \
PROTEINS_PER_TASK=1000 \
PARTITION=compute \
TIME=2-00:00:00 \
RESULTS_DIR=res_enet_zbig \
bash run_enet.sh
```

## Expected output files

```text
res_enet_zbig/<PROTEIN>_enet_8020_zscore_fast_k5000.h5
```

---

# Comparison notebook

The final comparison notebook expects these folders:

```text
results_xgb/
res_xgb_zbig/
results_enet/
res_enet_zbig/
```

It also expects same-gene correlation inputs:

```text
data/full_pearson100.joblib
data/full_spearman100.joblib
```

The comparison notebook compares:

```text
1. XGB K=2000
2. XGB z-score K=5000
3. ENet K=2000
4. ENet z-score K=5000
5. Same-name gene Pearson baseline
6. Same-name gene Spearman baseline
```

Metrics compared:

```text
test_pearson_corr
test_cosine_similarity
test_mse
```

Thresholds used for well-predicted proteins:

```text
Pearson: good if >= 0.70
Cosine:  good if >= 0.70
MSE:     good if <= selected MSE percentile threshold
```

---

# Notes

## XGBoost

The XGBoost models are saved as fitted `xgboost.XGBRegressor` objects inside each H5 model artifact.

The artifact includes:

```text
model
selected_features
normalizer
best_params
protein
```

## ElasticNet

The ElasticNet models are saved as fitted `sklearn.linear_model.ElasticNet` objects.

The feature weights are available as:

```python
artifact["model"].coef_
```

The matching gene order is:

```python
artifact["selected_features"]
```

## Binning note

For ENet binning, the number of bins is controlled by:

```bash
N_BINS=10
```

or:

```bash
N_BINS=50
```

For XGBoost qbin, the bin count is encoded directly in the normalization name:

```bash
NORMALIZATIONS=qbin50
```

## Re-running failed jobs

Both pipelines support resume behavior. Existing H5 files are skipped when rerunning the same command with the same output suffix.

This means failed or timed-out runs can usually be resumed by running the same command again, possibly with a smaller `PROTEINS_PER_TASK`.
