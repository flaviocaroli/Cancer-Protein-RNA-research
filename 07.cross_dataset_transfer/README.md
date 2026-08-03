# 07. Cross-Dataset Transfer (revised)

Train per-protein RNA→protein models on a **union of datasets** and predict on a **held-out dataset**,
with a **literature-grade, leakage-safe harmonization** fit on the training union, a selectable model
representation (**z-score / binning / ranking**), 5-fold CV tuning, ANOVA/Spearman feature selection
(k=5000), and **self-contained per-protein models that save ElasticNet gene coefficients**.

Datasets: **CCLE, TMT, iTRAQ, LabelFree**. The four splits:

| split name | train | test |
|---|---|---|
| `ccle_tmt__itraq_labelfree` | CCLE + TMT | iTRAQ **and** LabelFree |
| `ccle_tmt_itraq__labelfree` | CCLE + TMT + iTRAQ | LabelFree |
| `ccle_tmt_labelfree__itraq` | CCLE + TMT + LabelFree | iTRAQ |
| `ccle_labelfree_itraq__tmt` | CCLE + LabelFree + iTRAQ | TMT |

## Why this harmonization (not the 05.validation recipe)

For training one ML model on mixed platforms, the validated approach in the cross-platform literature
(Thompson et al. 2020, *Genome Biology*; PRONE 2025 for proteomics NMs) is a **reference-anchored
transform fit on the training data and applied frozen to the held-out set** — quantile→normal mapping
(default), robust-z, or rank. Batch-correction tools (ComBat/limma) are deliberately avoided: they can't
be frozen and projected onto a held-out platform. Here the reference is the **train union**, so no
held-out labels or distributions inform the transform. The model representation (z-score/binning/
ranking) is a separate, sweepable axis applied to the harmonized genes.

## Pipeline

**Notebook 00 — `notebooks/00_load_collection_a.ipynb`** (loader; imitates `05.validation` only here).
Groups Collection A cohorts into TMT/iTRAQ/LabelFree, concatenates, maps RNA Ensembl→symbol, filters
>50% missing, transposes to samples×features, saves per-platform joblibs + CCLE into `data/loaded/`.

**Notebook 02 — `notebooks/02_match_and_normalize.ipynb`** (split-controllable). Globals:
`TRAIN_DATASETS`, `TEST_DATASETS`, `HARMONIZATION_METHOD` (`quantile_normal`|`robust_z`|`standard_z`),
`MODEL_NORM` (`zscore`|`binning`|`ranking`|`none`). Builds the common feature space, **fits
harmonization on the train union**, applies it to all datasets, applies the chosen model
representation, prints dimensions, renders **Nature-style before/after** figures (PDF+PNG @300dpi), and
saves the split to `data/splits/<split_name>/`: `X_train, Y_train, X_test__<D>, Y_test__<D>`,
`common_genes/proteins`, `harmonization_params.joblib`, `split_config.json`.

**HPC training/prediction** (`src/` + `run_train.sh`/`run_predict.sh`). Reads a split dir; per protein:
select k=5000 genes (ANOVA-F for ENet, |Spearman| for XGB) → fit model normalizer on train, apply to
test → 5-fold CV → refit on train union → **save a self-contained artifact**.

## 20-job HPC cap

`run_train.sh` submits **one SLURM array of at most `PROTEIN_CHUNKS` (≤20) tasks**. Each task gets a
contiguous protein chunk and `train_transfer.py` **loops over the proteins in that chunk**. Default
`PROTEIN_CHUNKS=20`; the script refuses values >20. Run one (split, model, normalization) config at a
time to stay under the cap.

## Per-protein artifact (what gets saved)

`results/<split>/<model>_<norm>_k<k>/models/<PROTEIN>.joblib`:
- `features_coef_dict` (**gene → coefficient**) + `intercept` + `n_nonzero_coef`  ← ENet
- `booster` (fitted XGBRegressor)  ← XGB
- `selected_features`, `normalizer` (train-fitted), `best_params`, `cv_pearson_mean`, provenance.

Matches the `03. elasticNet` `features_coef_dict` format; `src/models.predict_with_artifact` re-predicts
on any dataset from the artifact alone.

## Run

```bash
pip install -r requirements.txt
# 0. load Collection A + CCLE  -> data/loaded/
jupyter nbconvert --to notebook --execute notebooks/00_load_collection_a.ipynb
# 1. choose a split + harmonization + model_norm in the globals, run:
jupyter nbconvert --to notebook --execute notebooks/02_match_and_normalize.ipynb   # -> data/splits/<split>/
# 2. train (<=20 jobs) then predict:
SPLIT_DIR=data/splits/ccle_tmt_itraq__labelfree MODEL=enet SELECT_K=5000 CV_FOLDS=5 PROTEIN_CHUNKS=20 bash run_train.sh
SPLIT_DIR=data/splits/ccle_tmt_itraq__labelfree MODEL=enet bash run_predict.sh
```

Sweep: `MODEL` in {enet, xgb}; `MODEL_NORM` in {zscore, binning, ranking} (re-run notebook 02 per norm);
all four splits. Local dry run: add `SUBMIT=local` to `run_train.sh`.

## Smoke test (no real data)

```bash
bash run_smoke_test.sh   # synthetic split -> chunked train (4 chunks) -> predict
```

## Notes to set on real data
- `RNA_TRANSFORM`: `None` if CCLE RNA is log2(TPM+1), else `"log2p1"`.
- `HARMONIZATION_METHOD`: `quantile_normal` is the recommended default; compare against `robust_z`.
- Feature intersection shrinks as you add platforms; notebook 02 prints the usable gene/protein counts
  per split — report them.
