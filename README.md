# Cancer-Protein-RNA-research

Interpretable prediction and analysis of protein abundance from RNA expression in cancer.

This repository studies how much of the proteome can be inferred from transcriptomics using a staged, protein-level framework:

```text
same-gene RNA proxy -> ElasticNet -> XGBoost -> unresolved proteins
```

The goal is not only to compare model performance, but to understand **which proteins are predictable by simple RNA-protein coupling**, which require **small explainable multigene models**, which benefit from **larger nonlinear models**, and which remain difficult to predict from RNA alone.

The project is organized as a reproducible analysis pipeline:

```text
pre-processing -> correlations -> model training -> CCLE comparison -> external validation
```

---

## Table of Contents

- [Project Overview](#project-overview)
- [Repository Structure](#repository-structure)
- [Datasets](#datasets)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Main Pipelines](#main-pipelines)
- [CCLE Comparison](#ccle-comparison)
- [External Validation](#external-validation)
- [Main Outputs](#main-outputs)
- [Reproducibility Notes](#reproducibility-notes)
- [Troubleshooting](#troubleshooting)
- [License](#license)

---

## Project Overview

This repository focuses on transcriptome-to-proteome prediction in cancer cell lines and external tumor proteomics datasets.

The main analyses include:

- preprocessing and matching CCLE transcriptomics/proteomics;
- computing RNA-protein correlation structure;
- training per-protein ElasticNet models;
- training per-protein XGBoost models;
- comparing models against same-gene RNA proxy baselines;
- building protein-level cascade partitions;
- validating trained models on external CPTAC-derived datasets;
- comparing tissue-level validation results against Xu-style atlas benchmarks;
- evaluating optional scTranslator predictions where available.

The central scientific idea is that protein prediction should be **protein-specific, metric-dependent, and interpretable**. Small models are used first. Larger models are only considered useful when they add predictive value beyond simpler baselines.

---

## Repository Structure

```text
.
├── 01. pre-processing/
│   ├── src/
│   ├── 01. pre-processing.ipynb
│   ├── 02. filtering_matched.ipynb
│   ├── 03.normalizations_plots.ipynb
│   └── README.md
│
├── 02. correlations/
│   ├── data/
│   ├── new_corr/
│   │   ├── corr_functions.py
│   │   ├── corr_pipeline.py
│   │   ├── new_run_corr.sh
│   │   └── smoke_test.ipynb
│   ├── src/
│   ├── 01. pearson_anlaysis.ipynb
│   ├── 02. prot_clustering.ipynb
│   ├── 03. cluster_anlaysis.ipynb
│   └── README.md
│
├── 03. elasticNet/
│   ├── model_results/
│   ├── new_pipeline_8020/
│   │   ├── src/
│   │   ├── run_enet.sh
│   │   └── enet_8020_smoke_test.ipynb
│   ├── new_pipeline_loo/
│   └── README.md
│
├── 04. xgb/
│   ├── model_results/
│   ├── new_pipeline_xgb/
│   │   ├── src/
│   │   ├── run_xgboost.sh
│   │   └── README.md
│   ├── pipeline_xgb_loo/
│   └── README.md
│
├── 04.5 ccle comparison/
│   ├── data/
│   ├── src/
│   ├── pan_atlas_downloads/
│   ├── 1. comparison_models.ipynb
│   ├── 2. final_comparison_ccle.ipynb
│   ├── 3.procedural_ccle.ipynb
│   ├── 3.1.proxy_z.ipynb
│   ├── 3.2. proxy_linear.ipynb
│   ├── 4. pan_atlas.ipynb
│   ├── 4.1 atlas_group.ipynb
│   ├── cascade_partition_outputs_metadata.md
│   └── README.md
│
├── 05. validation/
│   ├── data/
│   │   ├── Collection_A/ brain/
│   │   └── Collection_B/
│   │       └── CPTAC/
│   │           ├── brca_cptac_2020/
│   │           ├── gbm_cptac_2021/
│   │           ├── luad_cptac_2020/
│   │           ├── lusc_cptac_2021/
│   │           ├── paad_cptac_2021/
│   │           └── ucec_cptac_2020/
│   ├── notebooks/
│   │   ├── src/
│   │   ├── 7_test_8020.ipynb
│   │   ├── 8.prot_model_val.ipynb
│   │   ├── 9.procedural_cptac.ipynb
│   │   ├── 10.xu_comparison.ipynb
│   │   └── README.md
│   └── scTrasnlator-validation-newdata/
│
├── .gitignore
├── LICENSE
├── requirements.txt
└── README.md
```

---

## Datasets

### CCLE data

The main development benchmark uses public CCLE transcriptomics, proteomics, and metadata.

- Transcriptomics: <https://figshare.com/ndownloader/files/25797014>
- Proteomics: <https://gygi.hms.harvard.edu/data/ccle/Table_S2_Protein_Quant_Normalized.xlsx>
- Metadata: <https://figshare.com/ndownloader/files/25797011>

The preprocessing notebooks align CCLE sample identifiers, clean features, handle missingness, and save matched matrices for downstream analysis.

Expected processed files:

```text
transcriptomics_filtered.joblib
proteomics_matched.joblib
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

### External validation data

External validation is organized under:

```text
05. validation/data/
```

This includes brain/atlas-style data and CPTAC-derived cancer cohorts such as BRCA, GBM, LUAD, LUSC, PAAD, and UCEC.

---

## Installation

### 1. Clone the repository

```bash
git clone --branch draft https://github.com/flaviocaroli/Cancer-Protein-RNA-research.git
cd Cancer-Protein-RNA-research
```

### 2. Create a Python environment

```bash
conda create -n prot_rna python=3.10 -y
conda activate prot_rna
```

### 3. Install dependencies

The repository contains a `requirements.txt` file with pinned packages.

```bash
pip install -r requirements.txt
```

If pip has trouble because the requirements are stored on one line, create a cleaned requirements file:

```bash
python - <<'PY'
from pathlib import Path
txt = Path("requirements.txt").read_text().split()
Path("requirements.clean.txt").write_text("\n".join(txt) + "\n")
PY

pip install -r requirements.clean.txt
```

Main dependencies include:

```text
pandas
numpy
scipy
scikit-learn
matplotlib
seaborn
joblib
networkx
jupyterlab
ipykernel
tables
h5py
numexpr
```

For validation notebooks, install additional packages listed in the validation folder if needed.

---

## Quick Start

The safest way to reproduce the project is to run the pipeline in order.

### 1. Preprocess CCLE data

Open and run:

```text
01. pre-processing/01. pre-processing.ipynb
01. pre-processing/02. filtering_matched.ipynb
```

Optional normalization plots:

```text
01. pre-processing/03.normalizations_plots.ipynb
```

The expected output is a pair of matched matrices:

```text
transcriptomics_filtered.joblib
proteomics_matched.joblib
```

### 2. Compute RNA-protein correlations

For small tests, use:

```text
02. correlations/new_corr/smoke_test.ipynb
```

For local execution:

```bash
cd "02. correlations/new_corr"

python corr_pipeline.py \
  --data_dir /path/to/data \
  --results_dir /path/to/results/correlations \
  --chunk_id 0 \
  --group_size 500 \
  --p_threshold 0.05 \
  --n_bootstraps 100
```

For faster parametric correlations without bootstrap:

```bash
python corr_pipeline.py \
  --data_dir /path/to/data \
  --results_dir /path/to/results/correlations \
  --chunk_id 0 \
  --group_size 500 \
  --no_bootstrap \
  --p_threshold 0.05
```

For cluster execution, edit paths inside:

```text
02. correlations/new_corr/new_run_corr.sh
```

then run:

```bash
bash new_run_corr.sh
```

### 3. Train ElasticNet models

Small test:

```bash
cd "03. elasticNet/new_pipeline_8020"

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

Full-style run:

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

### 4. Train XGBoost models

Small test:

```bash
cd "04. xgb/new_pipeline_xgb"

NORMALIZATIONS=none,qbin50 \
GRID_PRESET=tiny \
MAX_PARAM_COMBOS=2 \
SELECT_K=50 \
EARLY_STOPPING=5 \
PROTEIN_START=1 \
PROTEIN_END=5 \
PROTEINS_PER_TASK=5 \
PARTITION=medium_cpu \
TIME=03:00:00 \
bash run_xgboost.sh
```

Full-style run:

```bash
bash run_xgboost.sh
```

Default XGBoost settings include:

```text
NORMALIZATIONS=none,zscore,rank,qbin50
PARTITION=compute
TIME=2-00:00:00
PROTEINS_PER_TASK=1000
SELECT_K=2000
EARLY_STOPPING=10
GRID_PRESET=fast
MAX_PARAM_COMBOS=8
```

### 5. Compare models on CCLE

Open:

```text
04.5 ccle comparison/3.procedural_ccle.ipynb
```

This notebook compares:

```text
1. same-gene Pearson proxy
2. same-gene Spearman proxy
3. ElasticNet
4. XGBoost
5. optional scTranslator predictions
6. unresolved proteins
```

The main cascade is:

```text
Proxy -> ElasticNet -> XGBoost -> scTranslator optional -> still not good
```

### 6. Validate externally on CPTAC-style datasets

Open notebooks in:

```text
05. validation/notebooks/
```

Important notebooks include:

```text
7_test_8020.ipynb
8.prot_model_val.ipynb
9.procedural_cptac.ipynb
10.xu_comparison.ipynb
```

These notebooks evaluate saved CCLE-trained models on external datasets, build CPTAC cascade partitions, and compare tissue-level outputs to Xu-style atlas benchmarks.

---

## Main Pipelines

### Correlation pipeline

Folder:

```text
02. correlations/new_corr/
```

Purpose:

- compute each protein against all RNA features;
- support chunked execution;
- support SLURM array jobs;
- save per-protein wide correlation outputs;
- save compact batch summaries.

Outputs:

```text
per_protein/
batch_*_corr_summary.(parquet|csv)
```

If parquet support is unavailable, CSV fallback is used.

---

### ElasticNet 80/20 pipeline

Folder:

```text
03. elasticNet/new_pipeline_8020/
```

Purpose:

- train one model per protein;
- perform feature selection on the training split;
- tune ElasticNet hyperparameters using inner cross-validation;
- evaluate on a held-out 20% test split;
- save a final all-samples model artifact.

Typical output folder:

```text
results_enet/
```

Saved H5 fields include:

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

ElasticNet models are saved as fitted `sklearn.linear_model.ElasticNet` objects. Coefficients can be inspected through:

```python
artifact["model"].coef_
artifact["selected_features"]
```

---

### XGBoost 80/20 pipeline

Folder:

```text
04. xgb/new_pipeline_xgb/
```

Purpose:

- train one XGBoost model per protein;
- use top-k feature selection;
- support several normalizations;
- tune hyperparameters using cross-validation;
- optionally use GPU support;
- save minimal H5 model artifacts.

Typical output folder:

```text
results_xgb/
```

Saved H5 fields include:

```text
protein
train_samples
test_samples
test_mse
test_pearson_corr
test_cosine_similarity
best_params
selected_features
final all-samples XGBoost model artifact
```

Example prediction from a saved H5:

```python
from joblib import load
from src.functions import load_minimal_8020_h5, predict_with_artifact

loaded = load_minimal_8020_h5("results_xgb/A2M_xgb_8020_none_fast_k2000.h5")
artifact = loaded["model_artifact"]

X_new = load("data/transcriptomics_filtered.joblib")
pred = predict_with_artifact(artifact, X_new)

print(pred.head())
```

---

## CCLE Comparison

Folder:

```text
04.5 ccle comparison/
```

This folder combines model outputs and baselines to produce protein-level comparisons.

Main comparison notebooks:

```text
1. comparison_models.ipynb
2. final_comparison_ccle.ipynb
3.procedural_ccle.ipynb
3.1.proxy_z.ipynb
3.2. proxy_linear.ipynb
```

Main model runs compared:

```text
XGBoost baseline, K=2000
XGBoost z-score big, K=5000
ElasticNet baseline, K=2000
ElasticNet z-score big, K=5000
same-gene Pearson proxy
same-gene Spearman proxy
```

Main metrics:

```text
test_pearson_corr
test_cosine_similarity
test_mse
```

The procedural CCLE notebook builds interpretable partitions such as:

```text
Proxy Pearson
Proxy Spearman
ENet
XGB
Still not good
Still not good after scTranslator
```

---

## External Validation

Folder:

```text
05. validation/
```

The validation notebooks apply saved model artifacts to external transcriptomics/proteomics datasets.

Main tasks:

- process external datasets;
- evaluate CCLE-trained models on LabelFree, TMT, and iTRAQ-style data;
- compute protein-level and tissue-level metrics;
- compare validation performance across platforms;
- build CPTAC cascade partitions;
- compare tissue-level results to Xu-style atlas metrics.

Important notebooks:

```text
05. validation/notebooks/7_test_8020.ipynb
05. validation/notebooks/8.prot_model_val.ipynb
05. validation/notebooks/9.procedural_cptac.ipynb
05. validation/notebooks/10.xu_comparison.ipynb
```

---

## Main Outputs

### Model artifacts

ElasticNet:

```text
03. elasticNet/new_pipeline_8020/results_enet/
03. elasticNet/model_results/
```

XGBoost:

```text
04. xgb/new_pipeline_xgb/results_xgb/
04. xgb/model_results/
```

### CCLE cascade outputs

The procedural notebooks can export protein-level cascade partitions.

Metadata for these outputs is described in:

```text
04.5 ccle comparison/cascade_partition_outputs_metadata.md
```

Typical exported files:

```text
cascade_partition_proteins_long.csv
cascade_partition_proteins_grouped.csv
cascade_partition_proteins_files_manifest.csv
protein_lists_one_file_per_partition/
```

Use cases:

- find proteins rescued by ElasticNet at a chosen Pearson threshold;
- identify proteins still unresolved after all models;
- compare cascade stages across CCLE and CPTAC validation datasets;
- create protein lists for downstream biological annotation.

### Validation outputs

Validation notebooks can generate:

```text
eval_predictions_long
eval_metrics_by_protein.csv
eval_metrics_by_tissue.csv
model_inventory.csv
feature_coverage.csv
```

Exact filenames may depend on notebook configuration and output directory settings.

---

## Reproducibility Notes

- Most notebooks and pipelines use relative paths; check the first configuration cells before running.
- Large raw data, model outputs, and H5 artifacts may not be committed to the repository.
- Feature selection is performed on the training data before test evaluation.
- The 80/20 pipelines save one compressed H5 file per protein per run.
- Minimal H5 files intentionally do not always store test predictions, RMSE, or full intermediate objects.
- Existing H5 output files are skipped when rerunning the same command, allowing failed jobs to be resumed.
- For large runs, reduce `PROTEINS_PER_TASK` if jobs fail because of memory or wall-time limits.
- For parquet outputs, install `pyarrow` or `fastparquet`; otherwise, many scripts fall back to CSV.

---

## Troubleshooting

### `pip install -r requirements.txt` fails

The requirements file may be stored as a single line. Create a cleaned version:

```bash
python - <<'PY'
from pathlib import Path
txt = Path("requirements.txt").read_text().split()
Path("requirements.clean.txt").write_text("\n".join(txt) + "\n")
PY

pip install -r requirements.clean.txt
```

### Missing `transcriptomics_filtered.joblib` or `proteomics_matched.joblib`

Run the preprocessing notebooks first:

```text
01. pre-processing/01. pre-processing.ipynb
01. pre-processing/02. filtering_matched.ipynb
```

Then copy or symlink the processed files into the `data/` folder expected by the pipeline.

### Correlation jobs are too slow

Use chunked execution:

```text
02. correlations/new_corr/corr_pipeline.py
```

Reduce:

```text
group_size
n_bootstraps
```

or use:

```bash
--no_bootstrap
```

### Cluster jobs fail or time out

Try smaller chunks:

```text
PROTEINS_PER_TASK=100
```

or use a shorter protein range for testing:

```text
PROTEIN_START=1
PROTEIN_END=5
```

### H5 outputs already exist

The ENet and XGB pipelines skip existing H5 files with the same output suffix. To rerun from scratch, delete or move the old output files.

### Validation paths do not resolve

Validation notebooks expect the original relative directory structure. Run them from:

```text
05. validation/notebooks/
```

and check paths such as:

```text
../val_data
../processed_data
../results_comparison_all_models
```

---

## Recommended Workflow for New Users

For a first pass, run only small tests:

```text
1. Preprocessing notebooks
2. Correlation smoke test
3. ElasticNet 5-protein test
4. XGBoost 5-protein test
5. CCLE procedural comparison notebook
```

After the small tests work, scale to full protein runs on a cluster.

---

## Research Interpretation

This repository supports a protein-centered interpretation of RNA-to-protein prediction.

The key modeling question is:

```text
What is the simplest model that predicts each protein well?
```

This allows proteins to be grouped into interpretable categories:

```text
same-gene RNA predictable
ElasticNet rescued
XGBoost rescued
scTranslator rescued
still unresolved
```

This is different from a standard leaderboard. The goal is not only to identify the highest-performing model, but to map the complexity needed to predict each protein.

---

## License

This repository is released under the MIT License. See:

```text
LICENSE
```
