# scTranslator on external platforms (optional, not in the paper)

Zero-shot scTranslator 2M / 160K / 10K predictions on the TMT, iTRAQ and label-free matrices.

## Setup
1. Copy `research_data/00_external/input/scTranslator/` to `<repo>/external/scTranslator/`:
   - `external/scTranslator/scTranslator/` – the scTranslator repo (old path `C:\Users\avsd8\OneDrive\Desktop\thesis\scTranslator stuff\scTranslator`)
   - `external/scTranslator/checkpoints/scTranslator_{2M,160K,10K}.pt` (old folder `…\scTranslator stuff\`)
2. Run `05. validation/notebooks/2. comparison.ipynb` first. The benchmarks read `05. validation/data/val_data/` directly.
3. Optional overrides: `SCTRANSLATOR_REPO`, `SCTRANSLATOR_CKPT_DIR`, `VAL_DATA_DIR`. `pip install torch`; a GPU is recommended.

## Notebooks and outputs
| Notebook | Inputs (`data/val_data/`) | Output folder |
|---|---|---|
| `benchmark_scTranslator_tmt.ipynb` | `TMT_*.joblib` | `benchmark_results_tmt/` |
| `benchmark_scTranslator_itraq.ipynb` | `iTRAQ_*.joblib` | `benchmark_results_itraq/` |
| `benchmark_scTranslator_label.ipynb` | `LabelFree_*.joblib` | `benchmark_results_labelfree/` |
| `overlap/benchmark_scTranslator_overlap_{tmt,itraq,labelfree}.ipynb` | same, restricted to `genes/overlap_ALL_4.txt` and `proteins/overlap_ALL_4.txt` | `overlap/benchmark_results_OVERLAP_{tmt,itraq,labelfree}/` |

Each output folder holds `predictions_<ckpt>.npz/.npy`, `results_<ckpt>.npz`, `summary.txt` and plots. The overlap folders also have tissue-level CSVs.

## Figures
Every notebook shows its figures inline and saves them to `figures_pdf/<notebook name>/` next to the notebook. The first code cell (*pRNA figure routing*) redirects every `savefig` there, keeping the original sub-folders. Do not delete that cell.
