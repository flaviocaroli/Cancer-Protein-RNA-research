import os
import sys
from joblib import load

current_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(current_dir, "..", "data"))
proteomics_file = os.path.join(data_dir, "proteomics_matched.joblib")

try:
    proteomics = load(proteomics_file)
except Exception as e:
    print(f"Error loading proteomics data: {e}")
    sys.exit(1)

# Extract protein names from the columns of the proteomics DataFrame
protein_list = proteomics.columns.tolist()

protein_list_file = os.path.join(data_dir, "protein_list.txt")

try:
    with open(protein_list_file, "w") as f:
        for protein in protein_list[0]: # REMOVE [0] TO CONSIDER ALL PROTEINS
            f.write(f"{protein}\n")
    print(f"Protein list written to {protein_list_file}, total proteins: {len(protein_list)}")
except Exception as e:
    print(f"Error writing protein list: {e}")
    sys.exit(1)
