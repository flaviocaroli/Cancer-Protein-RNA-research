# Cleaned CPTAC validation notebooks

This folder contains the refactored notebooks and shared helper module requested.

## Main changes

- `3.results_cptac.ipynb` replaces the three split CPTAC-method notebooks:
  - `3. results_tmt.ipynb`
  - `3. results_labelfree.ipynb`
  - `3. results_itraq.ipynb`
- `4.results_norms.ipynb` replaces the four split normalization notebooks:
  - `4. results_z.ipynb`
  - `4. results_none.ipynb`
  - `4. results_rk.ipynb`
  - `4. results_bin.ipynb`
- Long reusable helper functions were moved into `src/utils_val.py`.
- Original notebooks 1, 2, 5, 6, 7, 8, and 9 were updated to import helpers from `src/utils_val.py` and have outputs stripped.

## Expected run location

These notebooks keep the original relative paths such as `../val_data`, `../processed_data`, and `../results_comparison_all_models`. Run them from the same notebooks directory level as before, with `src/` adjacent to the notebooks.

## Dependencies

Install the package list in `requirements_validation.txt`. `pyarrow` is used for parquet files when available; the notebooks keep CSV fallbacks where the original code had them.
