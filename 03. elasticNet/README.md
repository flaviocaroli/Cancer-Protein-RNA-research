# 03. elasticNet 

## Files
This folder contains three folders for   **LOO** feature selection and regression execution.  

### /new_pipeline_loo
- **What it does:**: pipeline adapted to run chunks of proteins Elasticnet training and LOO testing. It includes nested feature selection (anova) and choice of normalization (z-score, none, ranking, binning) and 5 fold Cross validation for hyperparameter tuning.
- **Outputs:** `results/{PROTEIN}_LOO.csv`, predictions for each sample for each protein. Or best model for each protein (features coefficients, and hyperparameters).

### /new_pipeline_8020
- **What it does:** CLI wrapper for the 80/20 ElasticNet pipeline: runs CV-based hyperparameter tuning and trains either one protein or a protein chunk on the cluster.
- **Outputs:** compressed H5 files in results_enet/, named like {PROTEIN}_enet_8020_{normalization}_{grid}_k{select_k}.h5, containing split metrics, predictions, selected features, and final all-samples model artifact.

### /model_results
- **What it does:** Contains a notebook to analyze the results of the loo regression models. Also it contains notebooks with comparisons with scTranslator on CCLE.
- **Outputs:** In-notebook results and visualizations with possibility to save outputs.
---

### Quick start
1. Install dependencies (from repo root):
   ```bash
   pip install -r requirements.txt
   ```
2. Run each step as described in their respective READMEs, each folder contains a test to run feature selection or model on a subset of proteins/samples/genes of choice.

### Layout
```
.
new_pipeline_loo/
├─ src/
├─ new_run_enet.sh
└─ new_smoke.ipynb
new_pipeline_8020/
├─ src/
├─ enet_8020_smoke_test.ipynb
├─ run_enet.sh
└─ README.md
model_results/
├─ 01.res_analysis_loo.ipynb
├─ 02.enet_scT_comparison.ipynb
```
