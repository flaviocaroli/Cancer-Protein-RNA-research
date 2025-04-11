#!/usr/bin/env python
import os
import sys
import argparse
import pandas as pd
from joblib import load
from corr_functions import compute_correlation_for_protein

def main():
    parser = argparse.ArgumentParser(
        description='Compute correlation row for a single protein and save to an individual file.'
    )
    parser.add_argument('--protein', type=str, required=True, help='Protein name')
    args = parser.parse_args()
    protein_name = args.protein

    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(current_dir, "data") ## add '..' if there is a src folder
    transcriptomics_file = os.path.join(data_dir, "proteomics_matched.joblib")
    proteomics_file = os.path.join(data_dir, "proteomics_matched.joblib")
    output_dir = os.path.join(data_dir, "correlations")

    os.makedirs(output_dir, exist_ok=True)

    try:
        print("Loading data...")
        transcriptomics = load(transcriptomics_file)
        proteomics = load(proteomics_file)
        transcriptomics.index.name = "samples"
        proteomics.index.name = "samples"
        proteomics.columns.name = "proteins"
        transcriptomics.columns.name = "genes"
    except Exception as e:
        print(f"Error loading data: {e}")
        sys.exit(1)

    print(f"Computing correlation row for protein: {protein_name}")
    try:
        correlation_row = compute_correlation_for_protein(
            transcriptomics, proteomics, protein_name, p_threshold=0.05
        )
    except Exception as e:
        print(f"Error during computation for {protein_name}: {e}")
        sys.exit(1)

    # Save the individual correlation row to its own file.
    output_file = os.path.join(output_dir, f"{protein_name}_corr.csv")
    try:
        # Convert Series to a one-row DataFrame and save to CSV.
        correlation_row.to_frame().T.to_csv(output_file)
        print(f"Saved correlation row for protein {protein_name} to {output_file}")
    except Exception as e:
        print(f"Error saving result for protein {protein_name}: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
