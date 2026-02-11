# pipeline_proteins

## Files
This folder contains a **driver notebook** plus supporting scripts. Run the notebook first; the CLI mirrors the same logic.

### 01. pipeline_test.ipynb
- **What it does:** Orchestrates the four pipelines in `src/pipeline_functions.py` on a **controllable list of proteins**:
  - `fixed_percentile`, `variable_p`, `fixed_p`, and `fixed_features`.
- **Typical outputs:** Per‑protein Joblib files named
  `data/{protein}_{pipeline}_{model}_LOO_results.joblib`, and a combined CSV/Joblib (`ALL_{model}_LOO_results.*`).
- **Use cases:** Quick experimentation on selected proteins, side‑by‑side comparison of selection strategies and model hyper‑grids.

### 02. results_to_h5.ipynb
- **What it does:** (Optional) Loads the per‑protein Joblib results and consolidates them to an **HDF5** file (and/or CSV) for downstream analysis/plotting.
- **Typical outputs:** `data/pipeline_results.h5` with tables per pipeline, plus helper summaries.
- **Use cases:** Archival, fast reloads, and external visualization tools that prefer HDF5/CSV.

### 03. src/pipeline.py (CLI)
- **What it does:** Command‑line runner for a **single protein** with **leave‑one‑out** across all samples. Mirrors the notebook’s logic.
- **How to run (examples):**
  ```bash
  # ElasticNet + fixed percentile feature selection
  python src/pipeline.py --protein PTPN11 --pipeline fixed_percentile --model elasticnet

  # XGBoost + variable p-values (CPU)
  python src/pipeline.py --protein EGFR --pipeline variable_p --model xgboost
  ```
- **Notes:** GPU is optional for XGBoost. If you don’t have one, ensure `tree_method` uses `'hist'` (the notebook already toggles this via `XGB_USE_GPU=False`).

---

### Quick start
1. **Install dependencies** (from repo root):
   ```bash
   pip install -r requirements.txt
   ```

2. **Place required data** in `data/` (same level as `src/`):
   - `transcriptomics_filtered.joblib` → shape **(samples × features)**; index = sample IDs.
   - `proteomics_matched.joblib` → shape **(samples × proteins)**; index must match transcriptomics index one‑to‑one.
   - *(Optional, only for `fixed_features`)* precomputed features per protein:
     `{PROTEIN}_feature_selection_LOO.joblib` containing, for each LOO sample, the list of selected features.

3. **Run the notebook** `pipeline_test.ipynb`:
   - Edit the **Controls** cell: set `PROTEINS`, choose `MODEL_TYPE` (`'elasticnet'` or `'xgboost'`), and tweak CV/selection params.
   - Execute each **chapter** cell to run the corresponding pipeline.
   - Outputs are written to `data/` automatically; a final cell combines all results.
   - If using XGBoost without a GPU, set `XGB_USE_GPU=False`.

4. **(Alternative) Run via CLI** for one protein:
   ```bash
   python src/pipeline.py --protein <NAME> --pipeline <fixed_percentile|variable_p|fixed_p|fixed_features> --model <elasticnet|xgboost>
   ```

### Expected layout
```
pipeline_proteins/
├─ src/
│  ├─ pipeline_functions.py
│  ├─ pipeline.py
│  ├─ generate_protein_list.py
│  └─ __pycache__/...
├─ data/
│  ├─ transcriptomics_filtered.joblib
│  ├─ proteomics_matched.joblib
│  └─ {PROTEIN}_feature_selection_LOO.joblib        # optional, per-protein, for fixed_features
├─ pipeline_test.ipynb
├─ results_to_h5.ipynb
├─ requirements.txt
├─ run_enet.sh
└─ README.md
```

### Notes on models & parameters
- **ElasticNet**: Grid defaults in both notebook and CLI (`alpha`, `l1_ratio`). Adjust as needed for data scale.
- **XGBoost**: Defaults include `max_depth`, `learning_rate`, `n_estimators`, `subsample`, `colsample_bytree`.  
  - GPU path uses `tree_method='gpu_hist'` and `gpu_id=0`.
  - CPU path uses `tree_method='hist'`.
- **Reproducibility**: Seeds are fixed for NumPy and KFold; models are initialized with deterministic seeds when applicable.

### Outputs
- Per‑protein results: `data/{protein}_{pipeline}_{model}_LOO_results.joblib`
- Combined results (when using the notebook): `data/ALL_{model}_LOO_results.joblib` and `.csv`
- (Optional) HDF5 archive from `results_to_h5.ipynb`.

