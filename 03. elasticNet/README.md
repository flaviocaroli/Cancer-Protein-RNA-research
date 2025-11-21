# 03. elasticNet 

## Files
This folder contains three folders for   **LOO** feature selection and regression execution.  

### /feature_selection
- **What it does:** Implements two feature selectors for a *single protein* across leave-one-out splits:  
- **Outputs:** In-memory Python dicts/lists keyed by LOO fold; intended to be serialized by the CLI script.

### /pipeline_proteins
- **What it does:** CLI wrapper to run feature CV and regression training for **one protein** across **all samples (LOO)**, provides also .sh file to run everything on a cluster.  
- **Outputs:** `data/{PROTEIN}_LOO.joblib`, predictions.

### /model_results
- **What it does:** Contains a notebook to analyze the results of the regression model.
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
feature_selection/
├─ selection_functions.py
├─ selection_script.py
├─ run_feat.sh
└─ README.md
pipeline_proteins/
├─ src/
├─ enet_pipeline_test.py
├─ run_enet.sh
└─ README.md
├─ requirements.txt
model_results/
├─ res_analysis.ipynb
```
