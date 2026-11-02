#!/usr/bin/env python3
import argparse
import gc
import os
import re
import sys
from typing import List, Optional, Tuple

import pandas as pd
from joblib import load

from functions import (
    save_minimal_8020_h5,
    train_test_pipeline_for_protein,
    validate_normalization,
)


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Run an 80/20 XGBoost protein pipeline. The saved H5 contains only the requested "
            "minimal metrics/metadata plus a final model trained on all Dataset 1 samples."
        )
    )

    p.add_argument("--protein", default=None, type=str)
    p.add_argument("--protein_start", type=int, default=None, help="1-based inclusive start index for protein chunk mode.")
    p.add_argument("--protein_end", type=int, default=None, help="1-based inclusive end index for protein chunk mode.")

    p.add_argument("--splits", type=int, default=3)
    p.add_argument("--select_k", type=int, default=2000)
    p.add_argument("--early_stopping", type=int, default=10)
    p.add_argument("--test_fraction", type=float, default=0.2)
    p.add_argument("--random_state", type=int, default=42)
    p.add_argument("--n_samples", type=int, default=None)
    p.add_argument("--max_param_combos", type=int, default=8)
    p.add_argument("--grid_preset", choices=["tiny", "fast", "current"], default="fast")
    p.add_argument("--normalization", default="none", help="none, zscore, rank, or qbinN such as qbin50.")
    p.add_argument("--xgb_threads", type=int, default=None)
    p.add_argument("--force_gpu", action="store_true")
    p.add_argument("--no_feature_selection", action="store_true")

    p.add_argument("--metrics_dir", default=None, help="Output directory. Default: RESULTS_DIR or ./results_xgb.")
    p.add_argument("--resume_check_dir", default=None, help="Directory checked for existing H5 outputs. Default: metrics_dir.")
    p.add_argument("--resume_skip", action="store_true", default=True)
    p.add_argument("--no_resume_skip", action="store_false", dest="resume_skip")
    p.add_argument("--h5_compression_level", type=int, default=9)
    p.add_argument("--output_suffix", default="xgb_8020_fast_k2000")

    return p.parse_args()


def data_paths():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.abspath(os.path.join(current_dir, ".."))
    default_data_dir = os.path.join(project_dir, "..", "data")
    default_results_dir = os.path.join(project_dir, "..", "data", "results_xgb")

    data_dir = os.path.abspath(os.getenv("DATA_DIR", default_data_dir))
    results_dir = os.path.abspath(os.getenv("RESULTS_DIR", default_results_dir))
    x_path = os.path.join(data_dir, "transcriptomics_filtered.joblib")
    y_path = os.path.join(data_dir, "proteomics_matched.joblib")
    return project_dir, data_dir, results_dir, x_path, y_path


def safe_path_component(value: str) -> str:
    out = str(value).replace(os.sep, "_").replace("/", "_").replace("\\", "_")
    return re.sub(r"[^A-Za-z0-9._+=@-]+", "_", out)


def build_param_grid(use_gpu: bool, preset: str = "fast"):
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
    else:
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


def detect_protein_axis(X: pd.DataFrame, Y: pd.DataFrame) -> Tuple[List[str], str, int, int]:
    x_samples = pd.Index(map(str, X.index))
    y_index = pd.Index(map(str, Y.index))
    y_columns = pd.Index(map(str, Y.columns))
    index_overlap = len(x_samples.intersection(y_index))
    column_overlap = len(x_samples.intersection(y_columns))

    if index_overlap >= column_overlap:
        proteins = list(map(str, Y.columns))
        orientation = "samples_as_rows__proteins_as_columns"
    else:
        proteins = list(map(str, Y.index))
        orientation = "proteins_as_rows__samples_as_columns"
    return proteins, orientation, index_overlap, column_overlap


def select_proteins(
    all_proteins: List[str],
    *,
    protein: Optional[str],
    protein_start: Optional[int],
    protein_end: Optional[int],
) -> List[str]:
    if protein:
        return [protein]

    n = len(all_proteins)
    start = 1 if protein_start is None else int(protein_start)
    end = n if protein_end is None else int(protein_end)

    if start < 1:
        raise ValueError("--protein_start must be >= 1.")
    if end < start:
        raise ValueError("--protein_end must be >= --protein_start.")
    if start > n:
        raise ValueError(f"--protein_start={start} is beyond number of proteins={n}.")

    end = min(end, n)
    return all_proteins[start - 1:end]


def expected_output_path(protein: str, *, output_suffix: str, metrics_dir: str) -> str:
    return os.path.join(metrics_dir, f"{safe_path_component(protein)}_{output_suffix}.h5")


def run_one_protein(
    *,
    protein: str,
    X: pd.DataFrame,
    Y: pd.DataFrame,
    param_grid: dict,
    args,
    metrics_dir: str,
    use_gpu: bool,
):
    out_path = expected_output_path(protein, output_suffix=args.output_suffix, metrics_dir=metrics_dir)
    check_dir = os.path.abspath(args.resume_check_dir or metrics_dir)
    check_path = expected_output_path(protein, output_suffix=args.output_suffix, metrics_dir=check_dir)

    if args.resume_skip and os.path.exists(check_path):
        print(f"Skipping existing output for protein: {protein}")
        return "skipped"

    print("=" * 100)
    print(
        f"Running 80/20 XGBoost for protein '{protein}' | "
        f"grid={args.grid_preset} | normalization={args.normalization} | "
        f"select_k={args.select_k} | gpu={use_gpu}"
    )

    try:
        result = train_test_pipeline_for_protein(
            X=X,
            y=Y,
            protein_name=protein,
            param_grid=param_grid,
            normalization=args.normalization,
            test_fraction=args.test_fraction,
            random_state=args.random_state,
            n_splits=args.splits,
            select_k=args.select_k,
            early_stopping_rounds=args.early_stopping,
            max_param_combos=args.max_param_combos,
            n_samples=args.n_samples,
            use_gpu=use_gpu,
            xgb_threads=args.xgb_threads,
        )

        print(f"Saving H5: {out_path}")
        save_minimal_8020_h5(result, out_path, compression_level=args.h5_compression_level)

        summary = result["summary"]
        print("Summary metrics:")
        print(f"  protein: {summary['protein']}")
        print(f"  test_mse: {summary['test_mse']}")
        print(f"  test_pearson_corr: {summary['test_pearson_corr']}")
        print(f"  test_cosine_similarity: {summary['test_cosine_similarity']}")
        print(f"  train_samples: {len(summary['train_samples'])}")
        print(f"  test_samples: {len(summary['test_samples'])}")
        print(f"  selected_features: {len(summary['selected_features'])}")
        print(f"Saved output: {out_path}")
        return "done"

    except Exception as e:
        print(f"ERROR for protein {protein}: {type(e).__name__}: {e}", file=sys.stderr)
        return "failed"
    finally:
        gc.collect()


def main():
    args = parse_args()
    args.normalization = validate_normalization(args.normalization)

    if args.no_feature_selection:
        args.select_k = 0

    if not (0 < float(args.test_fraction) < 1):
        print("ERROR: --test_fraction must be between 0 and 1.", file=sys.stderr)
        sys.exit(1)

    project_dir, data_dir, default_results_dir, x_path, y_path = data_paths()
    metrics_dir = os.path.abspath(args.metrics_dir or default_results_dir)
    os.makedirs(metrics_dir, exist_ok=True)

    print("Loading data...")
    print(f"PROJECT_DIR: {project_dir}")
    print(f"DATA_DIR: {data_dir}")
    print(f"RESULTS_DIR: {metrics_dir}")
    print(f"X path: {x_path}")
    print(f"Y path: {y_path}")

    try:
        X = load(x_path)
        Y = load(y_path)
    except Exception as e:
        print(f"Error loading data ({x_path} / {y_path}): {e}", file=sys.stderr)
        sys.exit(1)

    X.index.name = "samples"
    X.columns.name = "genes"

    all_proteins, orientation, index_overlap, column_overlap = detect_protein_axis(X, Y)
    print(f"Transcriptomics shape: {getattr(X, 'shape', None)}")
    print(f"Proteomics shape: {getattr(Y, 'shape', None)}")
    print(f"Detected proteomics orientation: {orientation}")
    print(f"Sample overlap using Y rows/index: {index_overlap}")
    print(f"Sample overlap using Y columns: {column_overlap}")
    print(f"Number of detected proteins: {len(all_proteins)}")
    print(f"First detected proteins: {all_proteins[:5]}")

    try:
        proteins_to_run = select_proteins(
            all_proteins,
            protein=args.protein,
            protein_start=args.protein_start,
            protein_end=args.protein_end,
        )
    except Exception as e:
        print(f"Error selecting proteins: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Selected protein count: {len(proteins_to_run)}")
    print(f"Selected protein preview: {proteins_to_run[:5]}")

    use_gpu = bool(args.force_gpu)
    param_grid = build_param_grid(use_gpu=use_gpu, preset=args.grid_preset)
    statuses = {"done": 0, "skipped": 0, "failed": 0}

    for idx, protein in enumerate(proteins_to_run, start=1):
        print("=" * 100)
        print(f"Protein {idx}/{len(proteins_to_run)} in this pipeline call: {protein}")
        status = run_one_protein(
            protein=protein,
            X=X,
            Y=Y,
            param_grid=param_grid,
            args=args,
            metrics_dir=metrics_dir,
            use_gpu=use_gpu,
        )
        statuses[status] = statuses.get(status, 0) + 1

    print("=" * 100)
    print("Chunk finished.")
    print(f"Done: {statuses.get('done', 0)}")
    print(f"Skipped: {statuses.get('skipped', 0)}")
    print(f"Failed: {statuses.get('failed', 0)}")

    if statuses.get("failed", 0) > 0:
        sys.exit(2)


if __name__ == "__main__":
    main()
