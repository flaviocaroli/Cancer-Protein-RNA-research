# pre-process

Short analysis and EDA for preparing data before modeling/analysis.

## Files

- **01. pre-processing.ipynb** — Loads raw data and performs initial cleaning/transforms (type fixes, renaming/standardizing columns, missing-value handling, basic deduping/quality checks). Saves a tidy/curated dataset for downstream steps.
- **02. filtering_matched.ipynb** — Applies your matching & filtering rules on the curated data (e.g., keeping only records that meet match criteria, removing ambiguous/unwanted rows). Produces a final filtered dataset ready for exploration or modeling.

## Quick start

1. Install dependencies (from the repo root):  
   ```bash
   pip install -r requirements.txt
   ```
2. Open the notebooks **in order** (01 ➜ 02) and run top‑to‑bottom.
3. Adjust any input/output paths in the first cells of each notebook as needed.

## Expected layout

```
pre-process/
├─ 01. pre-processing.ipynb
├─ 02. filtering_matched.ipynb
└─ README.md
```

> Notes: The notebooks assume you have a `data/` directory available and sufficient disk space for intermediate files, if not, in the first notebook 01. pre-processing.ipynb it is ensured that the `data/` directory is created if it doesn't exist.
