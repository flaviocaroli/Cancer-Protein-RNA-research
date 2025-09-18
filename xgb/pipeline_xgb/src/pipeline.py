#!/usr/bin/env python3
import os
import sys
import argparse
from joblib import load, dump
import xgboost as xgb

from pipeline_functions import cv_pipeline_loo

# Back-compat aliases
ALIAS_MAP = {"single": "loo", "sample": "loo", "protein": "loo"}

def parse_args():
    p = argparse.ArgumentParser(
        description="Run XGBoost LOO pipeline (outer LOO over samples) for a single protein."
    )
    p.add_argument("--protein", required=True, type=str, help="Protein name (column in proteomics).")
    p.add_argument("--pipeline", required=True, type=str,
                   choices=["loo", "single", "sample", "protein"],
                   help="Pipeline type (only 'loo' supported; aliases map to 'loo').")
    p.add_argument("--model", required=True, type=str, choices=["xgboost"],
                   help="Model type (only xgboost supported).")

    # Speed / selection knobs
    p.add_argument("--splits", type=int, default=3, help="Inner CV folds (default: 3).")
    p.add_argument("--feature_selector", type=str, default="spearman",
                   choices=["pearson", "spearman", "anova"],
                   help="Univariate selector to rank features (if selection enabled).")
    p.add_argument("--select_k", type=int, default=500,
                   help="Keep top-K features inside each inner fold; <=0 disables feature selection.")
    p.add_argument("--preselect_k", type=int, default=5000,
                   help="Fast prefilter by |Pearson corr| before p-values (None to disable).")
    p.add_argument("--early_stopping", type=int, default=50,
                   help="Early stopping rounds in inner CV (0 disables).")
    p.add_argument("--n_samples", type=int, default=None,
                   help="If set, run only on the first N samples (outer LOO).")
    p.add_argument("--max_param_combos", type=int, default=None,
                   help="Randomly cap the number of hyperparam combos explored per inner CV.")

    # Saving knobs
    p.add_argument("--minimal_save", action="store_true",
                   help="Save only [protein, LOO sample, true, predicted].")
    p.add_argument("--force_gpu", action="store_true", help="Use GPU (gpu_hist).")
    p.add_argument("--no_feature_selection", action="store_true",
                   help="Disable preselect + p-value feature selection (use all features).")

    # Backward-compat NO-OP flags (accepted but ignored)
    p.add_argument("--percentile", type=float, default=None, help=argparse.SUPPRESS)
    p.add_argument("--pvalue", type=float, default=None, help=argparse.SUPPRESS)
    p.add_argument("--features_path", type=str, default=None, help=argparse.SUPPRESS)
    return p.parse_args()

def main():
    args = parse_args()
    pipeline_name = ALIAS_MAP.get(args.pipeline, args.pipeline)  # -> 'loo'

    # Optionally disable FS entirely
    if args.no_feature_selection:
        args.select_k = 0
        args.preselect_k = None

    # Paths relative to src/
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.abspath(os.path.join(current_dir, "..", "data"))
    X_path = os.path.join(data_dir, "transcriptomics_filtered.joblib")
    Y_path = os.path.join(data_dir, "proteomics_matched.joblib")

    print("Loading data...")
    try:
        X = load(X_path)
        Y = load(Y_path)
        X.index.name = "samples"; X.columns.name = "genes"
        Y.index.name = "samples"; Y.columns.name = "proteins"
    except Exception as e:
        print(f"Error loading data: {e}")
        sys.exit(1)

    protein = args.protein
    if protein not in Y.columns:
        print(f"Protein '{protein}' not in proteomics.")
        sys.exit(1)

    # Base model (unused directly but kept for signature)
    model = xgb.XGBRegressor(
        objective="reg:squarederror",
        n_jobs=-1,
        random_state=42,
        verbosity=0,
    )

    # Compact grids; rely on early stopping
    if args.force_gpu:
        param_grid = {
            "max_depth": [3, 6],
            "learning_rate": [0.1],
            "n_estimators": [1200],
            "subsample": [0.7, 0.9],
            "colsample_bytree": [0.5, 0.8],
            "reg_alpha": [0.0, 0.1],
            "reg_lambda": [1.0],
            "max_bin": [128],
            "tree_method": ["gpu_hist"],
            "predictor": ["gpu_predictor"],
        }
    else:
        param_grid = {
            "max_depth": [3, 6],
            "learning_rate": [0.1],
            "n_estimators": [1200],
            "subsample": [0.7, 0.9],
            "colsample_bytree": [0.5, 0.8],
            "reg_alpha": [0.0, 0.1],
            "reg_lambda": [1.0],
            "max_bin": [128],
            "tree_method": ["hist"],  # CPU fast histogram
        }

    print(f"Running {pipeline_name} for protein '{protein}' (LOO across samples)...")
    try:
        results_df = cv_pipeline_loo(
            model=model,
            param_grid=param_grid,
            X=X,
            y=Y,
            protein_name=protein,
            feature_selector=args.feature_selector,
            n_splits=args.splits,
            n_samples=args.n_samples,
            select_k=args.select_k,
            preselect_k=args.preselect_k,
            early_stopping_rounds=args.early_stopping,
            minimal_save=args.minimal_save,
            max_param_combos=args.max_param_combos,
        )
    except Exception as e:
        print(f"Error during pipeline execution: {e}")
        sys.exit(1)

    # Save compactly (try zstd -> fallback)
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
