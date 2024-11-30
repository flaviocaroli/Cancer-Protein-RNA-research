
import os
import sys
import argparse
import pandas as pd
import numpy as np
from joblib import load
from sklearn.linear_model import ElasticNet
from models_utils import full_pipeline2  # Import your pipeline function

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Run pipeline on a subset of proteins.')
parser.add_argument('--subset', type=int, required=True, help='Protein subset number (1-10)')
args = parser.parse_args()

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
transcriptomics = load(transcriptomics_file)
proteomics = load(proteomics_file)

# Load protein list for this subset
with open(protein_list_file, 'r') as f:
    protein_list = [line.strip() for line in f]

# Subset the proteomics data
proteomics_subset = proteomics.loc[protein_list]

# Define the ElasticNet model and parameter grid
elasticnet_model = ElasticNet()
param_grid_elasticnet = {
    'alpha': [0.1, 0.3],            # Regularization strength
    'l1_ratio': [0.1, 0.5, 0.9],    # Mix of L1 and L2 regularization
    'p_value': [0.1, 0.05, 0.01]    # P-values for feature selection
}

# Run the pipeline
print(f"Running ElasticNet pipeline on protein subset {subset_number}...")
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

# Save results
output_file = os.path.join(data_dir, f"elasticnet_scores_part{subset_number}.csv")
print(f"Saving results to {output_file}...")
elasticnet_scores.to_csv(output_file, index=True)
print("Pipeline finished successfully.")
