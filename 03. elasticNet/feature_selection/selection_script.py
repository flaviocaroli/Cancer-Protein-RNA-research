import os
import sys
import argparse
from joblib import load, dump
from selection_functions import feature_selection, naive_cluster_feature_selection

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description='Run feature selection pipeline on a single protein (LOO for all samples).'
)
parser.add_argument('--protein', type=str, required=True, help='Protein name')
parser.add_argument('--percentile', type=float, default=5, 
                    help='Percentile threshold (0-100) for feature selection (default: 5)')
args = parser.parse_args()

protein_name = args.protein.rstrip()
percentile = args.percentile

# Define paths
current_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(current_dir, "data")) ## add '..' if there is a src folder
transcriptomics_file = os.path.join(data_dir, "transcriptomics_kde_filtered.joblib")
proteomics_file = os.path.join(data_dir, "proteomics_matched.joblib")

# Load data
print("Loading data...")
try:
    transcriptomics = load(transcriptomics_file)  # DataFrame: samples x features
    proteomics = load(proteomics_file)              # DataFrame: samples x proteins
    transcriptomics.index.name = "samples"
    proteomics.index.name = "samples"
    proteomics.columns.name = "proteins"
    transcriptomics.columns.name = "genes"
except Exception as e:
    print(f"Error loading data: {e}")
    sys.exit(1)

# Check if the requested protein exists in the proteomics data
if protein_name not in proteomics.columns:
    print(f"Protein {protein_name} not found in proteomics data.")
    sys.exit(1)

print(f"Running feature selection pipeline for protein {protein_name} over all samples (LOO)...")
try:
    results_df = naive_cluster_feature_selection(
        X=transcriptomics,
        y=proteomics,
        protein_name=protein_name
    )
except Exception as e:
    print(f"Error during pipeline execution: {e}")
    sys.exit(1)

# Save results using Joblib
output_file = os.path.join(data_dir, f"{protein_name}_feature_selection_LOO.joblib")
try:
    print(f"Saving results to {output_file}...")
    dump(results_df, output_file)
    print("Pipeline finished successfully.")
except Exception as e:
    print(f"Error saving results: {e}")
    sys.exit(1)
