import os
from joblib import load
import numpy as np

# Paths
current_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(current_dir, "..", "data"))
proteomics_file = os.path.join(data_dir, "proteomics_matched.joblib")
protein_lists_dir = os.path.join(data_dir, "protein_lists")

# Load proteomics data
print("Loading proteomics data...")
proteomics = load(proteomics_file)

# Get the list of proteins
proteins = proteomics.T.index.tolist()
total_proteins = len(proteins)
print(f"Total number of proteins: {total_proteins}")

# Create protein_lists directory if it doesn't exist
if not os.path.exists(protein_lists_dir):
    os.makedirs(protein_lists_dir)

# Split proteins into 10 subsets (this can be changed according to the number of jobs available)
num_splits = 10
protein_subsets = np.array_split(proteins, num_splits)

# Save each subset to a separate text file
for i, subset in enumerate(protein_subsets, 1):
    list_filename = os.path.join(protein_lists_dir, f"protein_list_{i}.txt")
    with open(list_filename, 'w') as f:
        for protein in subset:
            f.write(f"{protein}\n")
    print(f"Saved protein subset {i} to {list_filename}")
