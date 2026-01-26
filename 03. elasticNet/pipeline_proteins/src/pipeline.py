import os
import sys
import argparse
import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNet
from joblib import load, dump
from pipeline_functions import cv_pipeline_fixed_percentile, cv_pipeline_fixed_p, cv_pipeline_fixed_features

# Utility function to select the pipeline
def select_pipeline(pipeline_name, model_type):
    """
    Selects the appropriate pipeline based on the input name and model type.
    """
    pipelines = {
        'fixed_percentile': cv_pipeline_fixed_percentile,  # Pipeline with fixed percentile
        'fixed_p': cv_pipeline_fixed_p,                    # Pipeline with fixed p-value
        'fixed_features': cv_pipeline_fixed_features       # Pipeline with fixed features
    }

    if pipeline_name not in pipelines:
        raise ValueError(f"Invalid pipeline name: {pipeline_name}. Choose from 'fixed_percentile', 'fixed_p', or 'fixed_features'.")
    
    if model_type not in ['elasticnet', 'xgboost']:
        raise ValueError(f"Invalid model type: {model_type}. Choose from 'elasticnet' or 'xgboost'.")
    
    return pipelines[pipeline_name], model_type

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Run pipeline on a single protein (LOO for all samples).')
parser.add_argument('--protein', type=str, required=True, help='Protein name')
parser.add_argument('--pipeline', type=str, required=True, choices=['fixed_percentile', 'variable_p', 'fixed_p'],
                    help="Pipeline type: 'fixed_percentile', 'variable_p' (variable p-values), or 'fixed_p' (fixed p-value)")
parser.add_argument('--model', type=str, required=True, choices=['elasticnet', 'xgboost'],
                    help="Model type: 'elasticnet' or 'xgboost'")
args = parser.parse_args()

protein_name = args.protein
pipeline_name = args.pipeline
model_type = args.model

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
    transcriptomics.index.name = "samples"
    proteomics.index.name = "samples"
    proteomics.columns.name = "proteins"
    transcriptomics.columns.name = "genes"
except Exception as e:
    print(f"Error loading data: {e}")
    sys.exit(1)

# Check if protein exists
if protein_name not in proteomics.columns:
    print(f"Protein {protein_name} not found in proteomics data.")
    sys.exit(1)

# Define models and parameter grids
if model_type == 'elasticnet':
    model = ElasticNet()
    param_grid = {
        'alpha': [0.1, 0.5, 1, 10],   # Regularization strengths
        'l1_ratio': [0.1, 0.5, 0.9],  # Mix of L1 and L2
    }

n_splits = 5

# Select the pipeline
pipeline, model_type = select_pipeline(pipeline_name, model_type)

print(f"Running {pipeline_name} pipeline for protein {protein_name} over all samples (LOO)...")
try:
    # Run the selected pipeline for the single protein
    if pipeline_name == 'fixed_percentile':
        results_df = pipeline(
            model=model,
            param_grid=param_grid,
            X=transcriptomics,
            y=proteomics,
            protein_name=protein_name,
            feature_selector='pearson',
            n_splits=n_splits,
            percentile=5  # Fixed percentile
        )
    elif pipeline_name == 'fixed_p':
        results_df = pipeline(
            model=model,
            param_grid=param_grid,
            X=transcriptomics,
            y=proteomics,
            protein_name=protein_name,
            feature_selector='pearson',
            n_splits=n_splits,
            p_value=0.05  # Fixed p-value
        )
    elif pipeline_name == 'fixed_features':
        results_df = pipeline(
            model=model,
            param_grid=param_grid,
            X=transcriptomics,
            y=proteomics,
            protein_name=protein_name,
            features=load(f"data/{protein_name}_feature_selection_LOO.joblib"),
            n_splits=n_splits
        )
except Exception as e:
    print(f"Error during pipeline execution: {e}")
    sys.exit(1)

# Save results using Joblib
output_file = os.path.join(data_dir, f"{protein_name}_{pipeline_name}_{model_type}_LOO_results.joblib")
try:
    print(f"Saving results to {output_file}...")
    dump(results_df, output_file)
    print("Pipeline finished successfully.")
except Exception as e:
    print(f"Error saving results: {e}")
    sys.exit(1)