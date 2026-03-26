# 03. elasticNet 

## Files
This folder contains three folders for   **LOO** feature selection and regression execution.  

### /new_pipeline
- **What it does:**: pipeline adapted to run chunks of proteins Elasticnet training and LOO testing. It includes nested feature selection (anova) and choice of normalization (z-score, none, ranking, binning) and 5 fold Cross validation for hyperparameter tuning.
- **Outputs:** `results/{PROTEIN}_LOO.csv`, predictions for each sample for each protein. Or best model for each protein (features coefficients, and hyperparameters).

### /single pipeline
- **What it does:** CLI wrapper to run feature CV and regression training for **one protein** at time (in cluster) across **all samples (LOO)**, provides also .sh file to run everything on a cluster.  
- **Outputs:** `data/{PROTEIN}_LOO.joblib`, predictions.

### /model_results
- **What it does:** Contains a notebook to analyze the results of the regression model. Also it contains notebooks with comparisons with
- **Outputs:** In-notebook results and visualizations with possibility to save outputs.
---

### Quick start
1. Install dependencies (from repo root):
   ```bash
   pip install -r requirements.txt
   ```
2. Run each step as described in their respective READMEs, each folder contains a test to run feature selection or model on a subset of proteins/samples/genes of choice.

### Expected layout
```
.
new_pipeline/
├─ src/
├─ new_run_enet.sh
└─ new_smoke.ipynb
pipeline_proteins/
├─ src/
├─ enet_pipeline_test.py
├─ run_enet.sh
└─ README.md
├─ requirements.txt
model_results/
├─ res_analysis.ipynb
├─ scTrasnslator-eNet-comparison_results/
└─ scTrasnlator-validation-newdata/

```
