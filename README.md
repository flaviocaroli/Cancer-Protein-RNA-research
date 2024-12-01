```markdown
# Cancer-Proteins-Prediction-Research

Prediction and statistical analysis of protein abundances in cancer cell lines using correlations, linear regressions, boosting methods, and deep learning.

## Datasets Used

- **Transcriptomics Dataset:** https://figshare.com/ndownloader/files/25797014  
- **Proteomics Dataset:** https://gygi.hms.harvard.edu/data/ccle/Table_S2_Protein_Quant_Normalized.xlsx 

## Usage Instructions

### Step 1: Clone the Repository
Clone the repository to your local machine:
```bash
git clone https://github.com/yourusername/Cancer-Proteins-Prediction-Research.git
cd Cancer-Proteins-Prediction-Research
```

### Step 2: Create a Data Directory
Since some datasets are larger than GitHub's file size limit, create an empty data folder in the `notebooks` directory to store these files:
```bash
mkdir notebooks/data
```

### Step 3: Run `correlations.ipynb`
This notebook downloads the datasets, performs filtering, and generates intermediate processed files. Open and run the notebook:
```bash
jupyter notebook notebooks/correlations.ipynb
```

**Output Files:**
- `transcriptomics_filtered.joblib`
- `proteomics_matched.joblib`

These files will be saved in the `notebooks/data` folder.

### Step 4: Copy Filtered Files
Copy the processed datasets to the appropriate directory for the pipeline:
```bash
cp notebooks/data/transcriptomics_filtered.joblib pipeline_proteins/data/
cp notebooks/data/proteomics_matched.joblib pipeline_proteins/data/
```

### Step 5: Split the Proteins into Subsets
Run the `split_prot.py` script to generate 10 subsets of proteins for parallel processing:
```bash
cd pipeline_proteins/src
python split_prot.py
```

**Output:**  
This script creates 10 text files (`protein_list_1.txt`, `protein_list_2.txt`, ..., `protein_list_10.txt`) in the `pipeline_proteins/data` folder.

### Step 6: Upload Files to Your HPC
Transfer the entire `pipeline_proteins` folder to your HPC:
```bash
scp -r ../pipeline_proteins your_username@hpc_address:/path/to/your/hpc/directory/
```

### Step 7: Set Up the Environment on the HPC
Log in to your HPC and set up a Python environment:
```bash
ssh your_username@hpc_address
module load miniconda3
conda create --name prot_env python=3.8 -y
conda activate prot_env
pip install -r pipeline_proteins/requirements.txt
```

### Step 8: Run SLURM Scripts for the Models
Submit the SLURM scripts to execute the models:

- **Run Elastic Net Regression:**
  ```bash
  cd pipeline_proteins
  sbatch run_enet.sh
  ```

- **Run XGBoost Regression:**
  ```bash
  sbatch run_xgboost.sh
  ```

### Step 9: Monitor Job Status
To check the status of your jobs:
```bash
squeue -u your_username
```

### Step 10: Retrieve Results
Once the jobs are complete, download the results back to your local machine:
```bash
scp -r your_username@hpc_address:/path/to/results_folder ./local_results_folder/
```
```
