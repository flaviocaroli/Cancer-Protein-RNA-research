# model_results

Analysis notebooks for the ElasticNet leave-one-out (LOO) results and the
ElasticNet-vs-scTranslator comparison. Notebooks are numbered in the order you
would normally run them; all of them read from / write to a sibling `res_data/`
folder (not committed — it holds the large result artifacts).

## Files

### 01. res_analysis_loo.ipynb
- **What it does:** Loads the consolidated ElasticNet LOO results and produces the
  per-protein / per-sample performance analysis and plots.
- **Inputs:** `res_data/enet_all_results.h5`, `res_data/correlations.joblib`.
- **Outputs:** In-notebook tables and figures (optionally saved).

### 02. enet_scT_comparison.ipynb
- **What it does:** Runs the full ElasticNet analysis twice (**RAW** and **SCALED**,
  per-sample min–max scaling scTranslator-style), then compares ElasticNet against
  the three scTranslator checkpoints (2M, 160K, 10K).
- **Inputs:** `res_data/enet_all_results.h5` and the prediction bundles
  `res_data/sc_predictions_{2M,160K,10K}.npz` (produced by notebook 03).
- **Outputs:** Metric figures (`elasticnet_metrics_*.png`, `comparison_plots/…`).

### 03. compute_sctranslator_predictions.ipynb
- **What it does:** Computes the scTranslator zero-shot protein-abundance predictions
  that notebook 02 consumes. For each checkpoint (2M, 160K, 10K) it maps genes/proteins
  to the scTranslator vocabulary, applies the paper's per-cell min–max normalization,
  runs inference (proteins are chunked to the decoder's 1000 limit) and saves one bundle.
- **Inputs:** matched `transcriptomics_filtered.joblib` + `proteomics_matched.joblib`,
  plus the scTranslator repo, checkpoints, and ID dictionaries (see *Prerequisites*).
- **Outputs:** `res_data/sc_predictions_{2M,160K,10K}.npz`, each with keys
  `predictions` (cells × proteins), `protein_names`, `protein_ids`, `cell_names`.

---

## Prerequisites for notebook 03

The scTranslator model code, checkpoints and ID dictionaries are **not** in this
repo — download them from the scTranslator release and point the notebook at them
via the config cell (or environment variables):

| What | Default env var | Contents |
|------|-----------------|----------|
| scTranslator source repo | `SCTRANSLATOR_REPO` | provides `code/model/performer_enc_dec.py` and `code/model/ID_dic/` |
| Checkpoint folder | `SCTRANSLATOR_CKPT_DIR` | `scTranslator_2M.pt`, `scTranslator_160K.pt`, `scTranslator_10K.pt` |
| Input data folder | `CPRR_DATA_DIR` | `transcriptomics_filtered.joblib`, `proteomics_matched.joblib` |
| Prediction output | `CPRR_PRED_OUT` | defaults to `res_data/` |

`transcriptomics_filtered.joblib` / `proteomics_matched.joblib` follow the same
convention as the ElasticNet pipeline: **rows = samples, columns = gene/protein symbols**,
with matched sample indices.

### Quick start
1. Install dependencies (from repo root): `pip install -r requirements.txt`
   (needs `torch`, plus `joblib numpy pandas scipy scikit-learn tqdm`).
2. Set the paths in the **Paths & Configuration** cell (or export the env vars above).
3. Run top to bottom. Predictions are written to `res_data/sc_predictions_<checkpoint>.npz`.

### Smoke test
Inference cost scales with `cells × proteins`; the full three-checkpoint run is heavy
(and very slow on CPU). For a quick check, keep only the `10K` entry in
`CONFIG['checkpoints']` — it produces one bundle and a per-cell / per-protein PCC readout.

## Notes / nuances
- **`torch.load` + PyTorch ≥ 2.6:** the checkpoints are full pickled model objects, so
  the notebook loads them with `weights_only=False`. Under the newer default
  (`weights_only=True`) the load fails with an *"Unsupported global scPerformerEncDec"*
  error. The load is only safe because these are our own trusted checkpoints.
- **NaN handling:** columns with **> 50 %** NaN across samples are dropped; the rest are
  filled with `1e-8` (0 is reserved as scTranslator's prediction token).
- **Unmapped symbols** (not in the scTranslator vocabulary) are discarded, so the mapped
  protein count can be smaller than the raw column count — the saved `protein_names`
  reflect exactly what was predicted, and notebook 02 realigns them to the ground truth.
- A GPU is strongly recommended; the notebook auto-selects CUDA when available.
