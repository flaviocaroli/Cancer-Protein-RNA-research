# pre-processing


## Files
This folder contains two notebooks. Run them **in order** (01 → 02).  
Filenames are kept as-is.

### 01. pre-processing.ipynb
- **What it does:** Loads raw data and performs initial cleaning/transforms (type fixes, renaming/standardizing columns, missing-value handling, deduping/quality checks).
- **Typical outputs:** A tidy/curated dataset (e.g., CSV/Parquet/Joblib) and optional summary logs/reports on schema and missingness.
- **Use cases:** Create a consistent, analysis-ready table for downstream steps; ensure reproducibility of early cleaning decisions.

### 02. filtering_matched.ipynb
- **What it does:** Applies matching & filtering rules on the curated data (e.g., enforce match criteria, remove ambiguous/unwanted rows).
- **Typical outputs:** Final filtered dataset ready for exploration/modeling, plus optional audit tables of kept/dropped rows.
- **Use cases:** Produce a high-confidence subset for EDA/modeling; make selection criteria explicit and repeatable.
---

### Quick start
1. Install dependencies (from repo root):
   ```bash
   pip install -r requirements.txt
   ```
2. Open each notebook and set input/output paths in the first configuration cell (defaults usually point to a `data/` and `figures/` folder).
3. Run all cells top-to-bottom.

### Expected layout
```
pre-process/
├─ 01. pre-processing.ipynb
├─ 02. filtering_matched.ipynb
└─ README.md
```
