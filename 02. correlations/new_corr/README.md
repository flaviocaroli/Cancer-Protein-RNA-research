# new_corr (chunked correlation runner)

This folder contains a chunked / SLURM-array-friendly rewrite of your original correlation workflow.

Instead of running **one protein per SLURM task**, we run **one chunk (group) of proteins per task** (same pattern as your `new_pipeline.py` + `new_run_enet.sh` style).

---

## Contents

- `corr_functions.py`
  - `compute_correlation_for_protein(...)`: protein vs all genes correlations (Pearson) with either:
    - **bootstrap CI** significance (default), or
    - **parametric p-value** threshold (`bootstrap=False`)
  - `run_correlation_group(...)`: loop over a list of proteins, optionally saving per-protein files, and returning a compact batch summary.

- `corr_pipeline.py`
  - CLI entrypoint to run **one chunk**:
    - chooses proteins from `protein_list.txt` based on `--chunk_id` + `--group_size` (or `--start_idx`)
    - loads your joblibs
    - runs `run_correlation_group`
    - saves:
      - per-protein wide outputs to `<results_dir>/per_protein/`
      - one batch summary file per chunk to `<results_dir>/`

- `new_run_corr.sh`
  - Submits **ONE SLURM array**: one task per chunk.
  - Chunking is controlled by:
    - `CHUNK_SIZE_OVERRIDE` (fixed proteins per task), OR
    - `TARGET_CHUNKS` (auto-compute group size)

- `smoke_test.ipynb`
  - Small synthetic test:
    - builds fake transcriptomics/proteomics
    - runs both single-protein and group functions
    - writes outputs to `./smoke_out/`

- `requirements.txt`
  - minimal Python deps

---

## Expected data layout

By default `corr_pipeline.py` assumes:

```
<data_dir>/
  transcriptomics_filtered.joblib
  proteomics_matched.joblib
  protein_list.txt
```

Outputs go to:

```
<results_dir>/
  per_protein/
    <protein>_corr.(parquet|csv)
  batch_<start>_<end>_corr_summary.(parquet|csv)
```

Notes:
- Per-protein files are **wide**: 1 row, genes as columns (values are correlation or 0).
- The batch summary is **compact**: it stores only significant (non-zero) correlations in `corr_json` (a JSON dict).

---

## Local usage (no SLURM)

Run the first chunk of 500 proteins:

```bash
python corr_pipeline.py \
  --data_dir /path/to/data \
  --results_dir /path/to/data/correlations \
  --chunk_id 0 \
  --group_size 500 \
  --p_threshold 0.05 \
  --n_bootstraps 100
```

Disable bootstrap and use parametric p-values:

```bash
python corr_pipeline.py \
  --data_dir /path/to/data \
  --results_dir /path/to/data/correlations \
  --chunk_id 0 \
  --group_size 500 \
  --no_bootstrap \
  --p_threshold 0.05
```

---

## SLURM usage

1) Copy this `new_corr/` folder into your project, and set paths in `new_run_corr.sh`:

- `BASE_DIR`
- `DATA_DIR`
- `SRC_DIR`

2) Submit:

```bash
bash new_run_corr.sh
```

---

## Install

```bash
pip install -r requirements.txt
```

If you want parquet outputs, install a parquet engine:

```bash
pip install pyarrow
# or
pip install fastparquet
```

(If parquet isn’t available, the code falls back to CSV automatically in a few places.)
