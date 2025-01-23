import os
import sys
import argparse
import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNet
from joblib import load
from single_pipeline import single_protein_pipeline

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Run pipeline on a single protein (LOO for all samples).')
parser.add_argument('--protein', type=str, required=True, help='Protein name')
args = parser.parse_args()

protein_name = args.protein

# Paths
current_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(current_dir, "..", "data"))
transcriptomics_file = os.path.join(data_dir, "transcriptomics_filtered.joblib")
proteomics_file = os.path.join(data_dir, "proteomics_matched.joblib")

# Load data
print("Loading data...")
try:
    transcriptomics = load(transcriptomics_file)  # samples x features
    proteomics = load(proteomics_file)            # samples x proteins
except Exception as e:
    print(f"Error loading data: {e}")
    sys.exit(1)

# Check if protein exists
if protein_name not in proteomics.columns:
    print(f"Protein {protein_name} not found in proteomics data.")
    sys.exit(1)

# Define the ElasticNet model and parameter grid
elasticnet_model = ElasticNet()
param_grid_elasticnet = {
    'alpha': [0.1, 0.3],          # Regularization strengths
    'l1_ratio': [0.1, 0.5, 0.9],  # Mix of L1 and L2
    'p_value': [0.1, 0.05, 0.01]  # P-values for feature selection
}

print(f"Running ElasticNet pipeline for protein {protein_name} over all samples (LOO)...")
try:
    # Run LOO pipeline for the single protein
    results_df = single_protein_pipeline(
        model=elasticnet_model,
        param_grid=param_grid_elasticnet,
        X=transcriptomics,
        y=proteomics,
        protein_name=protein_name,
        feature_selector='anova',
        n_splits=5
    )
except Exception as e:
    print(f"Error during pipeline execution: {e}")
    sys.exit(1)

# Save results
output_file = os.path.join(data_dir, f"{protein_name}_LOO_results.csv")
try:
    print(f"Saving results to {output_file}...")
    results_df.to_csv(output_file, index=False)
    print("Pipeline finished successfully.")
except Exception as e:
    print(f"Error saving results: {e}")
    sys.exit(1)
