# feature_selection

## Files
This folder contains three files for single-protein **LOO** feature selection and batch execution.  

### selection_functions.py
- **What it does:** Implements two feature selectors for a *single protein* across leave-one-out splits:  
  - `feature_selection(...)`: ANOVA (F-test), Pearson (r), and Spearman (ρ) with p-value percentile thresholds; returns per-LOO dicts of selected genes (with weights/scores).  
  - `naive_cluster_feature_selection(...)`: Pearson-only; selects top-N genes based on the protein’s **cluster** from `protein_clusters_k3.csv` using a `max_features` mapping; returns per-LOO lists.
- **Outputs:** In-memory Python dicts/lists keyed by LOO fold; intended to be serialized by the CLI script.

### selection_script.py
- **What it does:** CLI wrapper to run feature selection for **one protein** across **all samples (LOO)**.  
  Loads `data/transcriptomics_kde_filtered.joblib` (X) and `data/proteomics_matched.joblib` (Y); saves `data/{PROTEIN}_feature_selection_LOO.joblib`.  
  Uses `--protein` (required) and `--percentile` (default 5). Currently calls the **naive** selector.
- **Outputs:** `data/{PROTEIN}_feature_selection_LOO.joblib`, console logs.

### run_feat.sh
- **What it does:** SLURM batch launcher to process many proteins in chunks.  
  Reads `data/protein_list.txt`, submits an array job `[START..END]` with `CHUNK_SIZE`, and calls `selection_script.py --protein "$PROTEIN_NAME" --percentile $PERCENTILE`.  
  Adjust paths, log dirs, conda env, and resources as needed.
- **Outputs:** SLURM logs (`slurm-*.out`), one `*_feature_selection_LOO.joblib` per processed protein.

### feature_selection_test.ipynb
- **What it does:** Jupyter notebook to test and demonstrate the feature selection functions on a small subset of data.
- **Outputs:** In-notebook results and visualizations with possibility to save outputs.
---

### Quick start
1. Install dependencies (from repo root):
   ```bash
   pip install -r requirements.txt
   ```
2. Ensure required inputs exist in `data/`:
   - `transcriptomics_kde_filtered.joblib` (features, X)
   - `proteomics_matched.joblib` (targets, Y)
   - (for naive selector) `protein_clusters_k3.csv`
3. Run a single protein (example):
   ```bash
   python selection_script.py --protein TP53 --percentile 5
   ```
4. Batch on SLURM (edit paths/env in `run_feat.sh`, then):
   ```bash
   sbatch run_feat.sh
   ```
5. Also it is possible to run the selection functions directly from the jupyter notebook.


### Expected layout
```
.
├─ selection_functions.py
├─ selection_script.py
├─ run_feat.sh
└─ README.md
```
