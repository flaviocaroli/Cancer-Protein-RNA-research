# 01. pre-processing

Builds the matched CCLE matrices used by every later stage (paper: *CCLE preprocessing and feature filtering*, Fig. S1).

## Run order
| # | Notebook | Reads (`data/`) | Writes |
|---|---|---|---|
| 1 | `01_pre-processing.ipynb` | `CCLE_expression_full.csv`, `Table_S2_Protein_Quant_Normalized.xlsx` (or `.csv`), `sample_info.csv`. All three are downloaded automatically if missing. | `data/transcriptomics.csv`, `data/proteomics.csv` (harmonised), `data/correlations.joblib` (same-gene Pearson/Spearman/Kendall), `data/transcriptomics_matched.joblib`, `data/proteomics_matched.joblib`, `figures_pdf/*.pdf` |
| 2 | `02_filtering_matched.ipynb` | `data/transcriptomics_matched.joblib`, `data/proteomics_matched.joblib` | `data/transcriptomics_filtered.joblib` (KDE mean/SD filter → 22,614 genes), `figures_pdf/*.pdf` |
| 3 (optional) | `03.normalizations_plots.ipynb` | `data/transcriptomics_filtered.joblib` | `outputs/zscored.*`, `outputs/binned_B20.*`, `outputs/ranked.*` (illustration only) |

## Key outputs (used downstream)
| File | Shape | Used by |
|---|---|---|
| `transcriptomics_filtered.joblib` | 369 samples × 22,614 genes | 02, 03, 04, 05, 06 |
| `proteomics_matched.joblib` | 369 samples × 5,113 proteins | 02, 03, 04, 05, 06 |
| `correlations.joblib` | proteins × {pearson, spearman, kendall} | `03/model_results/01`, `04/model_results/02` |

## How to run
1. Copy `research_data/01_pre-processing/input/data/*` into `data/`, or let the notebook download the files.
2. Run the notebooks in the order above, from inside this folder.
3. Copy the outputs into `research_data/01_pre-processing/output/` and into the `data/` folders of the later stages (see their READMEs).

Helper code: `src/standardised_pipeline_utils.py`.

## Figures
Every notebook shows its figures inline and saves them to `figures_pdf/<notebook name>/` next to the notebook. The first code cell (*pRNA figure routing*) redirects every `savefig` there, keeping the original sub-folders. Do not delete that cell.
