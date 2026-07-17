# 04. xgb

## Files
This folder contains two folders for **XGBoost** model execution and one folder for result analysis.

### /new_pipeline_xgb
- **What it does:** CLI wrapper for the 80/20 XGBoost protein prediction pipeline. It runs one protein or protein chunks on the cluster, with Spearman top-k feature selection, optional normalization (`none`, `zscore`, `rank`, `qbin50`), CV-based hyperparameter tuning, and optional GPU support.
- **Outputs:** compressed H5 files in `results_xgb/`, named like `{PROTEIN}_xgb_8020_{normalization}_{grid}_k{select_k}.h5`, containing train/test samples, test MSE, Pearson correlation, cosine similarity, best parameters, selected features, and the final all-samples XGBoost model artifact.

### /pipeline_xgb_loo
- **What it does:** XGBoost pipeline/testing folder for LOO or subset-based execution. It contains smoke-test notebooks to validate small protein/sample runs and minimal H5 result loading before full cluster execution.
- **Outputs:** test pipeline outputs and minimal H5 result files, plus in-notebook validation of predictions, metrics, selected features, and saved model artifacts.

### /model_results
- **What it does:** Contains notebooks to merge and analyze XGBoost LOO regression outputs.
- **Outputs:** merged result tables, in-notebook result summaries, and visualizations.

---

### Quick start
1. Install dependencies from repo root:
   ```bash
   pip install -r requirements.txt