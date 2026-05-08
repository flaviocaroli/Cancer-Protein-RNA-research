#!/usr/bin/env python3
import os
import sys
import argparse
from joblib import load, dump

from pipeline_functions import cv_pipeline_loo


def parse_args():
    p = argparse.ArgumentParser(
        description="Run XGBoost LOO pipeline for a single protein."
    )

    # Required
    p.add_argument("--protein", required=True, type=str,
                   help="Protein name (column in proteomics).")

    # CV / speed controls
    p.add_argument("--splits", type=int, default=3,
                   help="Inner CV folds. Default: 3.")
    p.add_argument("--select_k", type=int, default=500,
                   help="Keep top-K features via Spearman prefilter per LOO split. 0 keeps all features.")
    p.add_argument("--early_stopping", type=int, default=30,
                   help="Early stopping rounds in inner CV. 0 disables early stopping.")
    p.add_argument("--n_samples", type=int, default=None,
                   help="If set, run only on the first N samples.")
    p.add_argument("--max_param_combos", type=int, default=None,
                   help="Randomly cap the number of hyperparameter combos explored.")
    p.add_argument("--tuning_strategy", choices=["nested", "once"], default="nested",
                   help=(
                       "nested = strict LOO: tune hyperparameters inside every LOO split; "
                       "once = faster: tune once per protein, then run LOO with fixed params. "
                       "Default: nested."
                   ))
    p.add_argument("--grid_preset", choices=["tiny", "fast", "current"], default="fast",
                   help="Hyperparameter grid size. tiny=4 combos, fast=16 combos, current=original 1944 combos.")
    p.add_argument(
        "--normalization",
        choices=["none", "zscore", "rank", "qbin5", "qbin10"],
        default="none",
        help=(
            "Transcriptomics normalization. "
            "none = raw data; "
            "zscore = fit mean/std inside each LOO train split; "
            "rank = rank-transform once outside LOO; "
            "qbin5/qbin10 = quantile-bin once outside LOO."
        ),
    )
    p.add_argument("--xgb_threads", type=int, default=None,
                   help="Threads passed to XGBoost. Default: SLURM_CPUS_PER_TASK or all visible CPUs.")

    # Convenience flags
    p.add_argument("--minimal_save", action="store_true",
                   help="Save only [protein, LOO sample, true, predicted].")
    p.add_argument("--force_gpu", action="store_true",
                   help="Use GPU with tree_method='hist', device='cuda'.")
    p.add_argument("--no_feature_selection", action="store_true",
                   help="Disable Spearman prefilter. Equivalent to --select_k 0.")
    p.add_argument("--output_suffix", default="loo_xgboost_LOO_results",
                   help="Output suffix before .joblib. Use different suffixes when comparing normalizations.")
    return p.parse_args()


def data_paths():
    """Use DATA_DIR if exported; otherwise fall back to ../data relative to this file."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    default_data_dir = os.path.abspath(os.path.join(current_dir, "..", "data"))
    data_dir = os.path.abspath(os.getenv("DATA_DIR", default_data_dir))
    X_path = os.path.join(data_dir, "transcriptomics_filtered.joblib")
    Y_path = os.path.join(data_dir, "proteomics_matched.joblib")
    return data_dir, X_path, Y_path


def build_param_grid(use_gpu: bool, preset: str = "fast"):
    """
    Backend params are fixed per run; they are not searched.

    current: original large grid, useful only for final confirmation.
    fast:    compact grid focused on the parameters that matter most.
    tiny:    very small pilot grid.
    """
    if preset == "current":
        common = {
            "max_depth": [3, 6, 8],
            "min_child_weight": [1, 3],
            "gamma": [0.0, 0.1],
            "learning_rate": [0.05, 0.1],
            "n_estimators": [600, 1200, 2000],
            "subsample": [0.7, 0.9],
            "colsample_bytree": [0.6, 0.8, 1.0],
            "reg_alpha": [0.0, 0.1, 0.5],
            "reg_lambda": [0.5, 1.0, 2.0],
            "max_bin": [256],
        }
    elif preset == "tiny":
        common = {
            "max_depth": [3, 4],
            "min_child_weight": [1, 3],
            "learning_rate": [0.05],
            "n_estimators": [600],
            "subsample": [0.9],
            "colsample_bytree": [0.8],
            "gamma": [0.0],
            "reg_alpha": [0.1],
            "reg_lambda": [1.0],
            "max_bin": [256],
        }
    else:  # fast
        common = {
            "max_depth": [3, 4],
            "min_child_weight": [1, 3],
            "learning_rate": [0.05, 0.1],
            "n_estimators": [600],
            "subsample": [0.9],
            "colsample_bytree": [0.6, 0.8],
            "gamma": [0.0],
            "reg_alpha": [0.1],
            "reg_lambda": [1.0],
            "max_bin": [256],
        }

    backend = {"tree_method": ["hist"], "device": ["cuda" if use_gpu else "cpu"]}
    return {**common, **backend}


def main():
    args = parse_args()

    if args.no_feature_selection:
        args.select_k = 0

    data_dir, X_path, Y_path = data_paths()
    print("Loading data...")
    print(f"DATA_DIR: {data_dir}")
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
    param_grid = build_param_grid(use_gpu=use_gpu, preset=args.grid_preset)

    print(
        "Running LOO XGBoost for protein "
        f"'{protein}' | grid={args.grid_preset} | tuning={args.tuning_strategy} | "
        f"normalization={args.normalization} | select_k={args.select_k} | gpu={use_gpu}"
    )

    try:
        results_df = cv_pipeline_loo(
            model=None,
            param_grid=param_grid,
            X=X,
            y=Y,
            protein_name=protein,
            n_splits=args.splits,
            n_samples=args.n_samples,
            select_k=args.select_k,
            early_stopping_rounds=args.early_stopping,
            minimal_save=args.minimal_save,
            max_param_combos=args.max_param_combos,
            use_gpu=use_gpu,
            normalization=args.normalization,
            tuning_strategy=args.tuning_strategy,
            xgb_threads=args.xgb_threads,
            grid_preset=args.grid_preset,
        )
    except Exception as e:
        print(f"Error during pipeline execution: {e}")
        sys.exit(1)

    out_path = os.path.join(data_dir, f"{protein}_{args.output_suffix}.joblib")
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
