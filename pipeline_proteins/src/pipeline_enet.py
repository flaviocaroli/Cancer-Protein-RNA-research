import os
import sys
import argparse
import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNet
from joblib import load
from single_pipeline import single_protein_pipeline, sample_selection_pipeline, protein_selection_pipeline

# Utility function to select the pipeline
def select_pipeline(pipeline_name):
    """
    Selects the appropriate pipeline based on the input name.
    """
    pipelines = {
        'single': single_protein_pipeline,
        'sample': sample_selection_pipeline,
        'protein': protein_selection_pipeline
    }
    if pipeline_name not in pipelines:
        raise ValueError(f"Invalid pipeline name: {pipeline_name}. Choose from 'single', 'sample', or 'protein'.")
    return pipelines[pipeline_name]

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Run pipeline on a single protein (LOO for all samples).')
parser.add_argument('--protein', type=str, required=True, help='Protein name')
parser.add_argument('--pipeline', type=str, required=True, choices=['single', 'sample', 'protein'],
                    help="Pipeline type: 'single', 'sample', or 'protein'")
args = parser.parse_args()

protein_name = args.protein
pipeline_name = args.pipeline

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
param_grid_elasticnet2 = {
    'alpha': [0.1, 0.5, 1, 10],   # Regularization strengths
    'l1_ratio': [0.1, 0.5, 0.9],  # Mix of L1 and L2
    # p-value is fixed at 0.1 for sample and protein pipelines
}

n_splits = 5
n_splits2 = 10

# Select the pipeline
pipeline = select_pipeline(pipeline_name)

# Choose the appropriate parameter grid based on the pipeline
if pipeline_name in ['sample', 'protein']:
    param_grid = param_grid_elasticnet2
else:
    param_grid = param_grid_elasticnet

print(f"Running {pipeline_name} pipeline for protein {protein_name} over all samples (LOO)...")
try:
    # Run the selected pipeline for the single protein
    results_df = pipeline(
        model=elasticnet_model,
        param_grid=param_grid,
        X=transcriptomics,
        y=proteomics,
        protein_name=protein_name,
        feature_selector='anova',
        n_splits=n_splits
    )
except Exception as e:
    print(f"Error during pipeline execution: {e}")
    sys.exit(1)

# Save results
output_file = os.path.join(data_dir, f"{protein_name}_{pipeline_name}_LOO_results.csv")
try:
    print(f"Saving results to {output_file}...")
    results_df.to_csv(output_file, index=False)
    print("Pipeline finished successfully.")
except Exception as e:
    print(f"Error saving results: {e}")
    sys.exit(1)