# pipeline_proteins_v2

## Files
This folder contains a **script-driven pipeline** (with SLURM support) plus supporting modules in `src/`. The workflow is designed for **chunked execution over many proteins** with flexible feature selection and normalization strategies.

---

### 01. src/new_pipeline.py (main CLI)
- **What it does:** Core entry point to run the pipeline on a **group of proteins** (chunked or manual selection).  
- Supports:
  - Feature selection pipelines:  
    `fixed_percentile`, `fixed_p`, `union`, `naive_cluster`
  - Normalization methods:  
    `zscore`, `binning`, `ranking`, `none`
- **Key feature:** Handles **chunking logic** for scalable runs.
- **Typical outputs:**
  - Per-protein best results:  
    `data/results/{PROTEIN}_best.parquet`
  - Batch results:  
    `data/results/batch_{start}_{end}_{pipeline}_{normalization}.parquet`

---

### 02. src/new_functions.py
- Implements core pipelines (ElasticNet + LOO + CV).
- Includes: `fixed_percentile`, `fixed_p`, `union`, `fixed_features`
- New: `naive_cluster` pipeline

---

### 03. src/norm_functions.py
- Optimized pipelines for:
  - `binning`
  - `ranking`
- Transformation computed once per run (faster).

---

### 04. new_run_enet.sh
- SLURM array launcher for parallel execution.
- Handles chunking, logging, and environment setup.

---

### 05. new_smoke.ipynb
- Optional notebook for quick testing/debugging.

---

## Quick start

### Install
```bash
pip install -r requirements.txt
```

### Run locally
```bash
python src/new_pipeline.py --group_size 50 --pipeline fixed_percentile --normalization zscore
```

### Run with chunking
```bash
python src/new_pipeline.py --chunk_id 0 --target_chunks 10 --pipeline fixed_percentile --normalization zscore
```

### Run on SLURM
```bash
bash new_run_enet.sh
```

---

## Expected layout
```
pipeline_proteins_v2/
├─ src/
│  ├─ new_pipeline.py
│  ├─ new_functions.py
│  ├─ norm_functions.py
├─ data/
│  ├─ transcriptomics_common.joblib
│  ├─ proteomics_matched.joblib
│  ├─ protein_list.txt
│  ├─ results/
│  └─ logs/
├─ new_run_enet.sh
├─ new_smoke.ipynb
├─ requirements.txt
└─ README.md
```

---

## Pipelines
- fixed_percentile
- fixed_p
- union
- naive_cluster

---

## Normalization
- zscore
- binning
- ranking
- none

---

## Outputs
- `data/results/{PROTEIN}_best.parquet`
- `data/results/batch_*.parquet`
- logs in `data/logs/`

---

## Notes
- Model: ElasticNet
- Outer CV: Leave-One-Out
- Inner CV: KFold
- Seeds fixed for reproducibility
