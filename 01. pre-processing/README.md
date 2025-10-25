# pre-processing

## Files
This folder contains two notebooks. Run them **in order** (01 → 02).  

### 01. pre-processing.ipynb
- **What it does:** Loads raw data and performs initial cleaning and transforms (type fixes, renaming/standardizing columns, missing-value handling, deduping/quality checks).
- **Outputs:** Curated datasets and optional summary reports on schema and missingness and distributions.

### 02. filtering_matched.ipynb
- **What it does:** Applies matching & filtering rules on the curated data (e.g., enforce match criteria, remove ambiguous/unwanted rows).
- **Outputs:** Final filtered dataset ready for exploration/modeling, plus optional audit tables of kept/dropped rows.
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
