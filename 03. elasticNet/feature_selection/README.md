# Feature Selection (LOO) — Quick Guide

## What’s here
- **`selection_functions.py`**: Implements two LOO feature selectors for a *single protein*:
  - `feature_selection(...)`: ANOVA (F-test), Pearson (r), and Spearman (ρ) with p-value percentile thresholds; returns per-LOO dicts of selected genes with weights. :contentReference[oaicite:0]{index=0}
  - `naive_cluster_feature_selection(...)`: Pearson-only; selects the top-N genes based on the protein’s **cluster** from `protein_clusters_k3.csv` using the `max_features` mapping. Returns per-LOO lists of selected genes. :contentReference[oaicite:1]{index=1}

- **`selection_script.py`**: CLI wrapper to run the pipeline for **one protein** across **all samples (LOO)**.  
  Loads `data/transcriptomics_kde_filtered.joblib` (X) and `data/proteomics_matched.joblib` (Y), then saves results to `data/{PROTEIN}_feature_selection_LOO.joblib`. Uses `--protein` (required) and `--percentile` (default 5). Currently calls the **naive** selector.

- **`run_feat.sh`**: SLURM batch launcher to process many proteins in chunks.  
  Reads `data/protein_list.txt`, submits an array job `[START..END]` with `CHUNK_SIZE`, and runs `selection_script.py --protein "$PROTEIN_NAME" --percentile $PERCENTILE`. Adjusts paths, log dirs, conda env, and resources.
