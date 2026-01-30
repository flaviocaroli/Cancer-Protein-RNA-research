#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np
from joblib import load

from selection_functions import naive_cluster_feature_selection_indices


def main():
    p = argparse.ArgumentParser(
        description="Run naive cluster-based LOO feature selection for a single protein, "
                    "saving a compact NPZ of indices."
    )
    p.add_argument("--protein", required=True, type=str, help="Protein name (column in proteomics).")
    p.add_argument("--percentile", type=float, default=10.0,
                   help="(Ignored for naive selection; kept for CLI compatibility.)")
    p.add_argument("--n_samples", type=int, default=None,
                   help="If set, only the first N samples are used.")
    p.add_argument("--data_dir", type=str, default=None,
                   help="Override data directory (default: $DATA_DIR or ./data)")
    args = p.parse_args()

    protein_name = args.protein.strip()

    # --- Resolve data dir ---
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = args.data_dir or os.getenv("DATA_DIR", os.path.join(current_dir, "data"))
    if not os.path.isdir(data_dir):
        print(f"ERROR: data_dir not found: {data_dir}")
        sys.exit(1)

    X_path = os.path.join(data_dir, "transcriptomics_filtered.joblib")
    Y_path = os.path.join(data_dir, "proteomics_matched.joblib")
    clusters_csv = os.path.join(data_dir, "protein_clusters_k3.csv")

    print("Loading data...")
    try:
        X = load(X_path)  # samples x genes
        Y = load(Y_path)  # samples x proteins
        X.index.name, X.columns.name = "samples", "genes"
        Y.index.name, Y.columns.name = "samples", "proteins"
    except Exception as e:
        print(f"Error loading data ({X_path} / {Y_path}): {e}")
        sys.exit(1)

    if protein_name not in Y.columns:
        print(f"Protein {protein_name} not found in proteomics data.")
        sys.exit(1)

    # Persist the gene list once (order matters; matches X.columns)
    genes_txt = os.path.join(data_dir, "gene_list.txt")
    if not os.path.exists(genes_txt):
        try:
            with open(genes_txt, "w", encoding="utf-8") as f:
                for g in X.columns:
                    f.write(f"{g}\n")
            print(f"Wrote gene list: {genes_txt}")
        except Exception as e:
            print(f"WARNING: Failed to write gene_list.txt ({e})")

    print(f"Running naive feature selection (cluster-aware) for '{protein_name}'...")
    try:
        res = naive_cluster_feature_selection_indices(
            X=X,
            y=Y,
            protein_name=protein_name,
            clusters_csv_path=clusters_csv,
            n_samples=args.n_samples,
        )
    except Exception as e:
        print(f"Error during feature selection: {e}")
        sys.exit(1)

    # Save compact NPZ
    out_npz = os.path.join(data_dir, f"{protein_name}_feature_selection_LOO_indices.npz")
    try:
        np.savez_compressed(
            out_npz,
            samples=res["samples"],
            indptr=res["indptr"],
            indices=res["indices"],
            n_genes=res["n_genes"],
        )
        print(f"Saved compact indices to: {out_npz}")
    except Exception as e:
        print(f"Error saving NPZ: {e}")
        sys.exit(1)

    # Optional: tiny counts CSV for quick QA (no heavy objects)
    try:
        counts = np.diff(res["indptr"].astype(np.int64))
        counts_csv = os.path.join(data_dir, f"{protein_name}_feature_selection_counts.csv")
        with open(counts_csv, "w", encoding="utf-8") as f:
            f.write("sample,selected_count\n")
            for s, c in zip(res["samples"], counts):
                f.write(f"{s},{int(c)}\n")
        print(f"Wrote counts CSV: {counts_csv}")
    except Exception as e:
        print(f"WARNING: Failed to write counts CSV ({e})")

    print("Done.")


if __name__ == "__main__":
    main()
