import os
import sys
import argparse
import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNet
from joblib import load, dump
import xgboost as xgb
from pipeline_functions import cv_pipeline_variable_p, cv_pipeline_fixed_p

# Utility function to select the pipeline
def select_pipeline(pipeline_name, model_type):
    """
    Selects the appropriate pipeline based on the input name and model type.
    """
    pipelines = {
        'variable_p': cv_pipeline_variable_p,  # Pipeline with variable p-values
        'fixed_p': cv_pipeline_fixed_p         # Pipeline with fixed p-value
    }
    if pipeline_name not in pipelines:
        raise ValueError(f"Invalid pipeline name: {pipeline_name}. Choose from 'variable_p' or 'fixed_p'.")
    
    if model_type not in ['elasticnet', 'xgboost']:
        raise ValueError(f"Invalid model type: {model_type}. Choose from 'elasticnet' or 'xgboost'.")
    
    return pipelines[pipeline_name], model_type

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Run pipeline on a single protein (LOO for all samples).')
parser.add_argument('--protein', type=str, required=True, help='Protein name')
parser.add_argument('--pipeline', type=str, required=True, choices=['variable_p', 'fixed_p'],
                    help="Pipeline type: 'variable_p' (variable p-values) or 'fixed_p' (fixed p-value)")
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
elif model_type == 'xgboost':
    model = xgb.XGBRegressor()
    param_grid = {
        'max_depth': [3, 6, 9],
        'learning_rate': [0.01, 0.1, 0.3],
        'n_estimators': [100, 200],
        'subsample': [0.8, 1.0],
        'colsample_bytree': [0.8, 1.0],
        'tree_method': ['gpu_hist'],
        'gpu_id': [0]
    }

n_splits = 5

# Select the pipeline
pipeline, model_type = select_pipeline(pipeline_name, model_type)

print(f"Running {pipeline_name} pipeline for protein {protein_name} over all samples (LOO)...")
try:
    # Run the selected pipeline for the single protein
    if pipeline_name == 'variable_p':
        results_df = pipeline(
            model=model,
            param_grid=param_grid,
            X=transcriptomics,
            y=proteomics,
            protein_name=protein_name,
            feature_selector='pearson',
            n_splits=n_splits
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