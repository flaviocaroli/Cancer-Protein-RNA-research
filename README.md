```markdown
# Cancer-Proteins-RNA-Research

Prediction and statistical analysis of protein abundances in cancer cell lines using correlations, linear models (Elastic Net), gradient boosting (XGBoost), and simple baselines — organized as a reproducible pipeline: **pre-processing → correlations → modeling (EN/XGB) → baselines → results**, with optional HPC/SLURM execution.

---

## Table of Contents

- [Repository Structure](#repository-structure)
- [Datasets](#datasets)
- [Installation](#installation)
- [Quick Start (Local)](#quick-start-local)
- [HPC / SLURM Workflow](#hpc--slurm-workflow)
- [Outputs](#outputs)
- [Reproducibility Notes](#reproducibility-notes)
- [Troubleshooting](#troubleshooting)
- [License](#license)

---

## Repository Structure

```

.
├─ 01. pre-processing/        # data alignment/cleaning utilities & notebooks
├─ 02. correlations/          # correlation notebooks/scripts + intermediate artifacts
├─ 03. elasticNet/            # Elastic Net training, configs, and (optionally) SLURM scripts
├─ 04. xgb/
│  ├─ clustering/             # clustering & overlap/consensus analysis outputs
│  ├─ model_results/          # XGBoost model artifacts, metrics, summaries
│  └─ pipeline_xgb/           # XGBoost pipeline code and (optionally) SLURM scripts
├─ 05. baselines/             # simpler baselines / 
├─ LICENSE
├─ requirements.txt
└─ README.md


## Datasets

We work with public CCLE resources:

- **Transcriptomics:** https://figshare.com/ndownloader/files/25797014  
- **Proteomics:** https://gygi.hms.harvard.edu/data/ccle/Table_S2_Protein_Quant_Normalized.xlsx  
- **Metadata:** https://figshare.com/ndownloader/files/25797011

We align CCLE sample IDs, remove excessive missing values, and persist filtered/matched matrices for downstream modeling.

---

## Installation

We create a clean environment and install dependencies from `requirements.txt`.

```bash
# Conda (recommended)
conda create -n prot_pred python=3.10 -y
conda activate prot_pred
pip install -r requirements.txt
````

**Pinned requirements (for reproducibility):**

```
pandas==2.3.2
numpy==1.24.4
seaborn==0.12.2
matplotlib==3.7.2
scikit-learn==1.7.2
scipy==1.11.1
joblib==1.2.0
networkx==3.1
jupyterlab>=4.0
ipykernel>=6.25
tables>=3.8
h5py>=3.8
numexpr>=2.8.4
ipython>=8.15.0
```

---

## Quick Start (Local)

### 1) Clone the repository

```bash
git clone https://github.com/yourusername/Cancer-Proteins-Prediction-Research.git
cd Cancer-Proteins-Prediction-Research
```

### 2) Create a data directory for intermediates

We keep large artifacts out of version control:

```bash
mkdir -p "02. correlations/data"
```

### 3) Run the correlations notebook

This step downloads/ingests datasets, filters/matches samples, and writes intermediates.

```bash
jupyter notebook "02. correlations/correlations.ipynb"
```

**Expected outputs (written to `02. correlations/data/`):**

* `transcriptomics_filtered.joblib`
* `proteomics_matched.joblib`

### 4) Stage filtered files for modeling pipelines

We copy intermediates into the modeling folders:

```bash
# Elastic Net
mkdir -p "03. elasticNet/data"
cp "02. correlations/data/transcriptomics_filtered.joblib" "03. elasticNet/data/"
cp "02. correlations/data/proteomics_matched.joblib"     "03. elasticNet/data/"

# XGBoost
mkdir -p "04. xgb/pipeline_xgb/data"
cp "02. correlations/data/transcriptomics_filtered.joblib" "04. xgb/pipeline_xgb/data/"
cp "02. correlations/data/proteomics_matched.joblib"       "04. xgb/pipeline_xgb/data/"
```

### 5) (Optional) Generate protein lists / split work

If a helper like `generate_protein_list.py` exists, we run it to create a batched protein list:

```bash
# adjust path if different in the repo
python "04. xgb/pipeline_xgb/src/generate_protein_list.py"
```

This writes a text file under the corresponding `data/` directory to parallelize jobs.

### 6) Train models locally

* **Elastic Net**: notebooks/scripts in `03. elasticNet/`
* **XGBoost**: notebooks/scripts in `04. xgb/pipeline_xgb/`

Artifacts and metrics are written under `03. elasticNet/...` and `04. xgb/model_results/`. Clustering/overlap outputs (if enabled) land in `04. xgb/clustering/`.

---

## HPC / SLURM Workflow

### 1) Upload modeling folders and prepared data

We copy only what’s needed:

```bash
# from repo root (local)
scp -r "03. elasticNet" "04. xgb" your_username@hpc_address:/path/to/work/
```

### 2) Set up the environment on the cluster

```bash
ssh your_username@hpc_address
module load miniconda3
conda create --name prot_env python=3.8 -y   # or python=3.10 if available site-wide
conda activate prot_env
pip install -r /path/to/work/requirements.txt
```

### 3) Submit jobs (if SLURM scripts are present)

```bash
# Elastic Net
cd /path/to/work/"03. elasticNet"
sbatch run_enet.sh

# XGBoost
cd /path/to/work/"04. xgb"/pipeline_xgb
sbatch run_xgboost.sh
```

### 4) Monitor & retrieve results

```bash
squeue -u your_username

# after completion, on local:
scp -r your_username@hpc_address:/path/to/work/"04. xgb"/model_results ./local_results_xgb/
# and similarly for Elastic Net outputs if saved to a results folder
```

---

## Outputs

* **Intermediates:** `02. correlations/data/`
* **Elastic Net results:** under `03. elasticNet/` (metrics, predictions, saved models where applicable)
* **XGBoost results:** `04. xgb/model_results/` (per-protein metrics, predictions, summaries)
* **Clustering & overlaps:** `04. xgb/clustering/` (contingency tables, overlap summaries, and—if enabled—plots)
* **Baselines:** `05. baselines/` for simplified checks and sanity tests

CSV summaries (overlaps/consensus and per-protein/per-sample metrics) are exported alongside figures for easy inspection whenever those components are enabled.

---

## Reproducibility Notes

* We set random seeds in modeling scripts/notebooks where applicable.
* We perform feature selection **inside** each cross-validation fold to avoid leakage.
* Paths are relative to the repository; large raw files remain out of version control.
* For HDF workflows, we ensure `tables` and `numexpr` are installed (included in `requirements.txt`).

---

## Troubleshooting

* **OOM / slow correlations**
  We reduce bootstrap iterations, switch correlation type (Pearson/Spearman/Kendall), or chunk computations.
* **Cluster labels shift between runs**
  We fix random seeds and persist the label map; PCA sanity checks help validate clustering stability.
* **Environment mismatches (local vs HPC)**
  We keep Python and package versions aligned; we recreate the env using the same `requirements.txt`.
* **Missing intermediates**
  We re-run the correlations notebook (Quick Start step 3) and confirm the two `.joblib` files exist before modeling.

---

## License

See [`LICENSE`](LICENSE) for details.

```
```
