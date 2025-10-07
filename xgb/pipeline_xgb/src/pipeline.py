#!/usr/bin/env python3
import os
import sys
import argparse
from joblib import load, dump

from pipeline_functions import cv_pipeline_loo

# ----------------------
# Arg parsing
# ----------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Run XGBoost LOO pipeline (outer LOO over samples) for a single protein."
    )
    # Required
    p.add_argument("--protein", required=True, type=str,
                   help="Protein name (column in proteomics).")

    # CV / speed controls
    p.add_argument("--splits", type=int, default=3,
                   help="Inner CV folds (default: 3).")
    p.add_argument("--select_k", type=int, default=500,
                   help="Keep top-K features via Spearman prefilter per LOO split (0 = keep all).")
    p.add_argument("--early_stopping", type=int, default=50,
                   help="Early stopping rounds in inner CV (0 disables).")
    p.add_argument("--n_samples", type=int, default=None,
                   help="If set, run only on the first N samples (outer LOO).")
    p.add_argument("--max_param_combos", type=int, default=None,
                   help="Randomly cap the number of hyperparam combos explored per inner CV.")

    # Convenience flags
    p.add_argument("--minimal_save", action="store_true",
                   help="Save only [protein, LOO sample, true, predicted].")
    p.add_argument("--force_gpu", action="store_true",
                   help="Use GPU (gpu_hist/gpu_predictor) if available.")
    p.add_argument("--no_feature_selection", action="store_true",
                   help="Disable Spearman prefilter (equivalent to --select_k 0).")
    return p.parse_args()

# ----------------------
# Small helpers
# ----------------------
def data_paths():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.abspath(os.path.join(current_dir, "..", "data"))
    X_path = os.path.join(data_dir, "transcriptomics_filtered.joblib")
    Y_path = os.path.join(data_dir, "proteomics_matched.joblib")
    return data_dir, X_path, Y_path

def build_param_grid(use_gpu: bool):
    common = {
        "max_depth": [3, 6, 8],
        "min_child_weight": [1, 3],          # stronger when larger
        "gamma": [0.0, 0.1],                  # minimum loss reduction to split
        "learning_rate": [0.05, 0.1],
        "n_estimators": [600, 1200, 2000],    # early stopping will clip this

        # sampling
        "subsample": [0.7, 0.9],
        "colsample_bytree": [0.6, 0.8, 1.0],

        # regularization
        "reg_alpha": [0.0, 0.1, 0.5],
        "reg_lambda": [0.5, 1.0, 2.0],

        # histogram params
        "max_bin": [256],
    }
    if use_gpu:
        return {**common,
                "tree_method": ["gpu_hist"],
                "predictor": ["gpu_predictor"]}
    else:
        return {**common,
                "tree_method": ["hist"]}

# ----------------------
# Main
# ----------------------
def main():
    args = parse_args()
    pipeline_name = "loo"

    # Handle convenience flag
    if args.no_feature_selection:
        args.select_k = 0

    # Load data
    data_dir, X_path, Y_path = data_paths()
    print("Loading data...")
    try:
        X = load(X_path)
        Y = load(Y_path)
        X.index.name, X.columns.name = "samples", "genes"
        Y.index.name, Y.columns.name = "samples", "proteins"
    except Exception as e:
        print(f"Error loading data ({X_path} / {Y_path}): {e}")
        sys.exit(1)

    protein = args.protein
    if protein not in Y.columns:
        print(f"Protein '{protein}' not found in proteomics.")
        sys.exit(1)

    use_gpu = bool(args.force_gpu)
    param_grid = build_param_grid(use_gpu=use_gpu)

    print(f"Running {pipeline_name} for protein '{protein}' (outer LOO across samples)...")
    try:
        results_df = cv_pipeline_loo(
            model=None,                         # kept for API compatibility; unused
            param_grid=param_grid,
            X=X,
            y=Y,
            protein_name=protein,
            n_splits=args.splits,
            n_samples=args.n_samples,
            select_k=args.select_k,            # Spearman prefilter size per outer split (0 = all)
            early_stopping_rounds=args.early_stopping,
            minimal_save=args.minimal_save,
            max_param_combos=args.max_param_combos,
            use_gpu=use_gpu,                   # works with the updated pipeline_functions.py
        )
    except Exception as e:
        print(f"Error during pipeline execution: {e}")
        sys.exit(1)

    # Save
    out_path = os.path.join(data_dir, f"{protein}_{pipeline_name}_xgboost_LOO_results.joblib")
    try:
        print(f"Saving: {out_path}")
        try:
            dump(results_df, out_path, compress=("zstd", 3))
        except Exception:
            dump(results_df, out_path, compress=3)
        print("Done.")
    except Exception as e:
        print(f"Error saving results: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
