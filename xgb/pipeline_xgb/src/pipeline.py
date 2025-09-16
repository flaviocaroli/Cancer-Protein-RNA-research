#!/usr/bin/env python3
import os
import sys
import argparse
from joblib import load, dump
import xgboost as xgb

# --- import pipeline functions from the same folder ---
# (we run from SRC_DIR, so a simple relative import works)
from pipeline_functions import (
    cv_pipeline_fixed_percentile,
    cv_pipeline_variable_p,
    cv_pipeline_fixed_p,
    cv_pipeline_fixed_features,
)

ALIAS_MAP = {
    "single": "fixed_percentile",
    "sample": "variable_p",
    "protein": "fixed_p",
}

def parse_args():
    p = argparse.ArgumentParser(
        description="Run XGBoost pipeline (LOO over samples) for a single protein."
    )
    p.add_argument("--protein", required=True, type=str, help="Protein name (column in proteomics).")
    p.add_argument(
        "--pipeline",
        required=True,
        type=str,
        choices=["fixed_percentile", "variable_p", "fixed_p", "fixed_features", "single", "sample", "protein"],
        help="Pipeline type."
    )
    p.add_argument("--model", required=True, type=str, choices=["xgboost"], help="Model type (only xgboost supported).")
    p.add_argument("--splits", type=int, default=5, help="Inner CV splits.")
    p.add_argument("--feature_selector", type=str, default="pearson", choices=["pearson", "spearman", "anova"])
    p.add_argument("--percentile", type=float, default=5.0, help="For fixed_percentile.")
    p.add_argument("--pvalue", type=float, default=0.05, help="For fixed_p.")
    p.add_argument("--features_path", type=str, default=None, help="For fixed_features (joblib list of features).")
    p.add_argument("--force_gpu", action="store_true", help="Force GPU (restrict param grid to gpu_hist only).")
    return p.parse_args()

def main():
    args = parse_args()

    pipeline_name = ALIAS_MAP.get(args.pipeline, args.pipeline)

    # Paths relative to src/
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.abspath(os.path.join(current_dir, "..", "data"))
    X_path = os.path.join(data_dir, "transcriptomics_filtered.joblib")
    Y_path = os.path.join(data_dir, "proteomics_matched.joblib")

    print("Loading data...")
    try:
        X = load(X_path)   # samples x genes
        Y = load(Y_path)   # samples x proteins
        X.index.name = "samples"; X.columns.name = "genes"
        Y.index.name = "samples"; Y.columns.name = "proteins"
    except Exception as e:
        print(f"Error loading data: {e}")
        sys.exit(1)

    protein = args.protein
    if protein not in Y.columns:
        print(f"Protein '{protein}' not in proteomics.")
        sys.exit(1)

    model = xgb.XGBRegressor(
        objective="reg:squarederror",
        n_jobs=-1,
        random_state=42,
        verbosity=0,
    )

    if args.force_gpu:
        param_grid = {
            "max_depth": [3, 6, 9],
            "learning_rate": [0.01, 0.1, 0.3],
            "n_estimators": [100, 200],
            "subsample": [0.8, 1.0],
            "colsample_bytree": [0.8, 1.0],
            "tree_method": ["gpu_hist"],
        }
    else:
        param_grid = {
            "max_depth": [3, 6, 9],
            "learning_rate": [0.01, 0.1, 0.3],
            "n_estimators": [100, 200],
            "subsample": [0.8, 1.0],
            "colsample_bytree": [0.8, 1.0],
            "tree_method": ["gpu_hist", "hist"],
        }

    print(f"Running {pipeline_name} for protein '{protein}' (LOO across samples)...")
    try:
        if pipeline_name == "fixed_percentile":
            results_df = cv_pipeline_fixed_percentile(
                model=model, param_grid=param_grid, X=X, y=Y, protein_name=protein,
                feature_selector=args.feature_selector, n_splits=args.splits,
                percentile=float(args.percentile),
            )
        elif pipeline_name == "variable_p":
            results_df = cv_pipeline_variable_p(
                model=model, param_grid=param_grid, X=X, y=Y, protein_name=protein,
                feature_selector=args.feature_selector, n_splits=args.splits,
                p_values=[0.1, 0.05, 0.01],
            )
        elif pipeline_name == "fixed_p":
            results_df = cv_pipeline_fixed_p(
                model=model, param_grid=param_grid, X=X, y=Y, protein_name=protein,
                feature_selector=args.feature_selector, n_splits=args.splits,
                p_value=float(args.pvalue),
            )
        elif pipeline_name == "fixed_features":
            features_path = args.features_path or os.path.join(data_dir, f"{protein}_feature_selection_LOO.joblib")
            if not os.path.exists(features_path):
                print(f"Fixed-features file not found: {features_path}")
                sys.exit(1)
            results_df = cv_pipeline_fixed_features(
                model=model, param_grid=param_grid, X=X, y=Y, protein_name=protein,
                features=load(features_path), n_splits=args.splits,
            )
        else:
            raise ValueError(f"Unknown pipeline: {pipeline_name}")
    except Exception as e:
        print(f"Error during pipeline execution: {e}")
        sys.exit(1)

    out_path = os.path.join(data_dir, f"{protein}_{pipeline_name}_xgboost_LOO_results.joblib")
    try:
        print(f"Saving: {out_path}")
        dump(results_df, out_path)
        print("Done.")
    except Exception as e:
        print(f"Error saving results: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
