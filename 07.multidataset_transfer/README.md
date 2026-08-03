# 07. Multi-dataset platform transfer

This folder implements strict cross-platform RNA-to-protein transfer experiments using CCLE and the Collection A tumor cohorts.

It is deliberately split into three stages:

1. **Notebook 01 / Collection A loader**: merge raw Collection A cohorts into TMT, iTRAQ, and LabelFree matrices.
2. **Notebook 02 / scenario preparation**: choose train/test datasets, match features, learn training-only normalization parameters, save prepared matrices and diagnostic figures.
3. **HPC training**: train one ElasticNet or XGBoost model per protein using five-fold inner CV, then predict the untouched held-out platform(s).

The HPC launcher submits at most **20 jobs**. Each job loads the prepared bundle once and loops sequentially over one protein chunk.

## Implemented transfer scenarios

| Scenario | Training datasets | Held-out test datasets |
|---|---|---|
| `ccle_tmt_to_itraq_labelfree` | CCLE + TMT | iTRAQ and LabelFree, evaluated separately |
| `ccle_tmt_itraq_to_labelfree` | CCLE + TMT + iTRAQ | LabelFree |
| `ccle_tmt_labelfree_to_itraq` | CCLE + TMT + LabelFree | iTRAQ |
| `ccle_labelfree_itraq_to_tmt` | CCLE + LabelFree + iTRAQ | TMT |

There is intentionally no TMT+iTRAQ or TMT+LabelFree two-dataset training scenario. The only two-dataset training scenario is CCLE+TMT, as requested.

## Folder structure

```text
07.multidataset_transfer/
├── configs/
│   ├── collection_a_groups.yaml
│   ├── datasets.yaml
│   ├── scenarios.yaml
│   └── model_defaults.yaml
├── data/
│   ├── Collection_A/                  # raw Collection A cohort folders
│   ├── processed/
│   │   ├── CCLE/                      # user-supplied CCLE matrices
│   │   └── platforms/                 # notebook 01 outputs
│   └── reference/
│       └── ensembl_to_symbol.tsv      # recommended local mapping
├── notebooks/
│   ├── 01_load_collection_a.ipynb
│   └── 02_match_split_normalize.ipynb
├── prepared/                          # notebook 02 outputs
├── results/                           # per-protein models and metrics
├── figures/                           # PDF/SVG/PNG diagnostics
├── scripts/
│   ├── run_collection_loader.sh
│   ├── run_prepare.sh
│   ├── run_prepare_all_scenarios.sh
│   ├── run_train.sh
│   ├── run_aggregate.sh
│   ├── run_compare_normalizations.sh
│   └── run_smoke_test.sh
└── src/
    ├── utils.py
    ├── build_collection_platforms.py
    ├── prepare_scenario.py
    ├── make_tasks.py
    ├── train_task.py
    ├── aggregate_results.py
    ├── compare_normalizations.py
    └── predict_saved_model.py
```

## 1. Environment

```bash
cd 07.multidataset_transfer
conda env create -f environment.yml
conda activate multidataset_transfer
```

A virtual environment is also fine:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## 2. Input data

### Collection A

Copy the raw collection exactly as used by `05. validation`:

```text
data/Collection_A/<cohort>/protein_matrix_*.txt
data/Collection_A/<cohort>/RNA_matrix_*.txt
```

The platform assignments are in `configs/collection_a_groups.yaml`.

### CCLE

Place the matched CCLE matrices here:

```text
data/processed/CCLE/transcriptomics.joblib
data/processed/CCLE/proteomics.joblib
data/processed/CCLE/metadata.tsv        # optional but recommended
```

Matrices may be samples × features or features × samples; the dataset loader chooses the orientation with the largest RNA/protein sample-ID overlap. After loading, all downstream matrices are samples × features.

### Ensembl-to-symbol mapping

For reproducibility, provide:

```text
data/reference/ensembl_to_symbol.tsv
```

Accepted columns include `ensembl_gene_id` plus `gene_symbol`. If unavailable, notebook 01 can query MyGene once and save a cache, but a fixed local mapping is preferable.

## 3. Notebook 01: Collection A loader only

Open:

```bash
jupyter lab notebooks/01_load_collection_a.ipynb
```

The notebook mirrors only the Collection A loader from `05. validation`:

- groups cohorts into TMT, iTRAQ, and LabelFree;
- aligns RNA/protein samples within each cohort;
- gives samples deterministic cohort-prefixed IDs;
- maps RNA and protein Ensembl IDs to gene symbols;
- averages duplicate symbols;
- filters features and samples at the configured non-missing threshold;
- saves platform joblibs and metadata.

Equivalent command-line run:

```bash
MAPPING_PATH=$PWD/data/reference/ensembl_to_symbol.tsv \
FORCE=1 \
bash scripts/run_collection_loader.sh
```

Outputs:

```text
data/processed/platforms/TMT/transcriptomics.joblib
data/processed/platforms/TMT/proteomics.joblib
data/processed/platforms/TMT/metadata.tsv

... and the corresponding iTRAQ and LabelFree folders.
```

## 4. Notebook 02: split, match, normalize, save

Open:

```bash
jupyter lab notebooks/02_match_split_normalize.ipynb
```

The global configuration cell controls:

```python
SCENARIO = "ccle_tmt_to_itraq_labelfree"
TRAIN_DATASETS = None
TEST_DATASETS = None
OUTPUT_NAME = None
X_NORMALIZATION = "zscore"   # zscore, rank, binning, none
Y_NORMALIZATION = "domain_zscore"
N_BINS = 50
EPS = 1e-8
DOMAIN_BALANCED_NORMALIZER = True
FORCE_REBUILD = False
```

To use one of the four predefined splits, change only `SCENARIO`. For a custom split, provide `TRAIN_DATASETS`, `TEST_DATASETS`, and a unique `OUTPUT_NAME`.

### Leakage control

For every scenario:

- the held-out test dataset does not contribute to fitted normalization parameters;
- z-score means/SDs and bin edges are learned on training datasets only;
- test RNA is transformed using those saved parameters unchanged;
- rank normalization is sample-wise and parameter-free;
- the matched feature names may use the test feature list, but no test abundance values are used to fit a transformation or model;
- five-fold model CV refits normalization and feature selection inside each training fold;
- the full-training normalized matrices produced by notebook 02 are used only for the final fit after hyperparameter selection.

Equivalent command:

```bash
SCENARIO=ccle_tmt_to_itraq_labelfree \
X_NORM=zscore \
Y_NORM=domain_zscore \
N_BINS=50 \
MAKE_FIGURES=1 \
FORCE=1 \
bash scripts/run_prepare.sh
```

Prepare all four scenarios for z-score, rank, and binning:

```bash
NORMALIZATIONS=zscore,rank,binning \
bash scripts/run_prepare_all_scenarios.sh
```

## 5. Normalization definitions

### RNA z-score

A single gene-wise mean and SD is learned from the pooled training datasets using domain-balanced sample weights. Near-zero SDs are replaced by one. Missing values are imputed with the training mean before scaling. The same parameters are applied to the held-out platform.

### RNA rank

Each sample is replaced by within-sample percentile ranks across the matched gene universe. This does not fit test-cohort parameters and is robust to platform-wide scale shifts, but it removes absolute expression amplitude.

### RNA binning

For every gene, quantile edges are learned from a domain-balanced training view. Values are converted to normalized bin indices in `[0,1]`. The held-out platform uses the training edges unchanged.

### No RNA normalization

Missing values are imputed with training medians; otherwise values remain on the original scale. This is a control and is unlikely to be reliable across heterogeneous RNA pipelines.

### Protein target normalization

`domain_zscore` is the recommended target definition. Each protein is standardized separately within every training domain before models are fitted. This prevents TMT, iTRAQ, LabelFree, and CCLE numerical scales from being treated as directly interchangeable.

Consequences:

- predictions are dimensionless relative-protein-abundance scores;
- held-out Pearson and Spearman against raw protein abundance are the primary strict-transfer metrics;
- held-out MSE/RMSE/R²/calibration are calculated against a clearly labeled evaluation-only standardized test target;
- evaluation-only test scaling is never used to fit the model or generate predictions and must not be described as deployment calibration.

## 6. Model training

### Feature selection

Feature selection is repeated inside every CV fold and once on the complete training set for the final model.

- ElasticNet: top 5,000 genes by univariate F-regression/ANOVA score.
- XGBoost: top 5,000 genes by absolute Spearman correlation.
- Default `domain_rank_mean`: calculate feature scores within each training domain, convert them to percentile ranks, and average the ranks so the largest domain does not dominate.

### Five-fold CV

Both model families use five-fold inner CV on the combined training samples. The default splitter is stratified by training dataset/platform. The default selection score is mean Pearson correlation across training domains.

### ElasticNet HPC run

```bash
SCENARIO=ccle_tmt_to_itraq_labelfree \
X_NORM=zscore \
Y_NORM=domain_zscore \
MODEL=enet \
SUBMIT=slurm \
PARTITION=compute \
TIME=3-00:00:00 \
MEM=24G \
CPUS_PER_TASK=1 \
N_JOBS=1 \
PROTEIN_CHUNKS=20 \
MAX_HPC_JOBS=20 \
ARRAY_THROTTLE=20 \
SAVE_PREDICTIONS=1 \
SAVE_ALL_FEATURE_ROWS=0 \
bash scripts/run_train.sh
```

This submits at most 20 array tasks. If the protein universe contains 5,000 proteins, each task loops through approximately 250 proteins sequentially.

### XGBoost HPC run

```bash
SCENARIO=ccle_tmt_to_itraq_labelfree \
X_NORM=zscore \
Y_NORM=domain_zscore \
MODEL=xgb \
SUBMIT=slurm \
PARTITION=compute \
TIME=3-00:00:00 \
MEM=32G \
CPUS_PER_TASK=4 \
N_JOBS=4 \
PROTEIN_CHUNKS=20 \
MAX_HPC_JOBS=20 \
ARRAY_THROTTLE=5 \
SAVE_PREDICTIONS=1 \
SAVE_ALL_FEATURE_ROWS=0 \
bash scripts/run_train.sh
```

`ARRAY_THROTTLE=5` does not create extra jobs; it limits concurrent XGB tasks while all 20 remain in one array submission.

### Protein subset

To train only selected proteins:

```bash
PROTEIN_SUBSET=/path/to/proteins.tsv \
PROTEIN_CHUNKS=20 \
... \
bash scripts/run_train.sh
```

The subset file needs a `protein` column or a first column containing protein names.

## 7. Per-protein model artifacts

Every protein gets a reusable joblib:

```text
results/<scenario>/x-<norm>__y-<target>/<model>/models_per_protein/<model>/<protein>.joblib
```

The artifact contains:

- fitted estimator;
- selected genes and feature-selection scores;
- all 5,000 ElasticNet coefficients or XGB importances aligned to selected genes;
- intercept and nonzero count for ElasticNet;
- best hyperparameters and complete CV table;
- training/test dataset names and preparation manifest;
- paths to the shared RNA normalizer and target transformer;
- held-out metrics and, when enabled, held-out predictions.

To keep aggregate files manageable, feature shards contain nonzero ElasticNet coefficients/nonzero XGB importances by default. Set `SAVE_ALL_FEATURE_ROWS=1` to also write all 5,000 rows per protein to shard tables. The per-protein joblib always stores all coefficients/importances.

## 8. Aggregate a completed run

```bash
SCENARIO=ccle_tmt_to_itraq_labelfree \
X_NORM=zscore \
Y_NORM=domain_zscore \
MODEL=enet \
bash scripts/run_aggregate.sh
```

Important outputs:

```text
metrics_all.parquet
predictions_all.parquet
selected_features_all.parquet
summary_by_model_test.tsv
figures/transfer_pearson_by_test_dataset.pdf
```

CSV-gzip fallbacks are used when Parquet is unavailable.

## 9. Compare z-score, rank, and binning

After preparing, training, and aggregating all three normalizations:

```bash
SCENARIO=ccle_tmt_to_itraq_labelfree \
Y_NORM=domain_zscore \
MODEL=enet \
NORMALIZATIONS=zscore,rank,binning \
bash scripts/run_compare_normalizations.sh
```

Outputs are saved under:

```text
results/normalization_comparisons/<scenario>/<model>/
```

## 10. Reuse a saved protein model

```bash
python src/predict_saved_model.py \
  --model results/.../models_per_protein/enet/TP53.joblib \
  --rna /path/to/new_samples_by_genes.joblib \
  --out TP53_predictions.parquet
```

The normalizer path is read from the model artifact. Supply `--normalizer` only when the prepared folder has moved.

## 11. Recommended execution order

For each scenario:

1. Run notebook 01 once for Collection A.
2. Run notebook 02 separately for `zscore`, `rank`, and `binning`.
3. Run ElasticNet for each normalization.
4. Aggregate and compare normalizations.
5. Lock one RNA normalization before interpreting XGBoost.
6. Run XGBoost with the locked normalization.
7. Compare model transfer per held-out dataset, per protein, and per source cohort where metadata allows.

## 12. Smoke test

The smoke test generates synthetic CCLE and Collection A files, builds the platform matrices, prepares a scenario, fits per-protein models, aggregates results, and reuses one saved model:

```bash
bash scripts/run_smoke_test.sh
```
