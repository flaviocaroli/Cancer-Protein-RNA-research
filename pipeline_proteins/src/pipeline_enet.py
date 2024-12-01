import os
import argparse
import pandas as pd
import numpy as np
from joblib import load
from sklearn.linear_model import ElasticNet
from models_utils import full_pipeline2  # Import your custom pipeline function

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Run pipeline on a subset of proteins.')
parser.add_argument('--subset', type=int, required=True, help='Protein subset number (1-10)')
args = parser.parse_args()

# Validate subset number
subset_number = args.subset
if not 1 <= subset_number <= 10:
    raise ValueError("Subset number must be between 1 and 10.")

# Paths
current_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(current_dir, "..", "data"))
protein_lists_dir = os.path.join(data_dir, "protein_lists")
transcriptomics_file = os.path.join(data_dir, "transcriptomics_filtered.joblib")
proteomics_file = os.path.join(data_dir, "proteomics_matched.joblib")
protein_list_file = os.path.join(protein_lists_dir, f"protein_list_{subset_number}.txt")

# Load data
print("Loading data...")
try:
    transcriptomics = load(transcriptomics_file)
    proteomics = load(proteomics_file)
except Exception as e:
    print(f"Error loading data: {e}")
    sys.exit(1)

# Load protein list for this subset
try:
    with open(protein_list_file, 'r') as f:
        protein_list = [line.strip() for line in f]
except FileNotFoundError:
    print(f"Protein list file not found: {protein_list_file}")
    sys.exit(1)

# Ensure protein list matches proteomics data columns
protein_list = [protein.strip() for protein in protein_list]  # Remove extra spaces
missing_proteins = [p for p in protein_list if p not in proteomics.columns]
if missing_proteins:
    print(f"Warning: The following proteins are missing in the proteomics dataset: {missing_proteins}")
    protein_list = [p for p in protein_list if p in proteomics.columns]

# Subset the proteomics data
try:
    proteomics_subset = proteomics.loc[:, protein_list]
    if proteomics_subset.empty:
        raise ValueError("Proteomics subset is empty. Check protein list or dataset.")
except Exception as e:
    print(f"Error subsetting proteomics data: {e}")
    sys.exit(1)

# Define the ElasticNet model and parameter grid
elasticnet_model = ElasticNet()
param_grid_elasticnet = {
    'alpha': [0.1, 0.3],            # Regularization strength
    'l1_ratio': [0.1, 0.5, 0.9],    # Mix of L1 and L2 regularization
    'p_value': [0.1, 0.05, 0.01]    # P-values for feature selection
}

# Run the pipeline
print(f"Running ElasticNet pipeline on protein subset {subset_number}...")
try:
    elasticnet_scores = full_pipeline2(
        model=elasticnet_model,
        param_grid=param_grid_elasticnet,
        X=transcriptomics,
        y=proteomics_subset,
        feature_selector='anova',
        n_splits=5,
        n_samples=None,
        n_proteins=None
    )
except Exception as e:
    print(f"Error during pipeline execution: {e}")
    sys.exit(1)

# Save results
output_file = os.path.join(data_dir, f"elasticnet_scores_part{subset_number}.csv")
try:
    print(f"Saving results to {output_file}...")
    elasticnet_scores.to_csv(output_file, index=True)
    print("Pipeline finished successfully.")
except Exception as e:
    print(f"Error saving results: {e}")
    sys.exit(1)
