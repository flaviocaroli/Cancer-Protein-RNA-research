#!/usr/bin/env python3
import argparse
import gc
import os
import re
import sys
from typing import List, Optional, Tuple

import pandas as pd
from joblib import load, dump

from pipeline_functions import cv_pipeline_loo, save_minimal_protein_h5


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Run XGBoost LOO pipeline for one protein or a chunk of proteins. "
            "Default output mode is one compressed H5 per protein."
        )
    )

    p.add_argument("--protein", default=None, type=str)
    p.add_argument("--protein_start", type=int, default=None)
    p.add_argument("--protein_end", type=int, default=None)

    p.add_argument("--splits", type=int, default=3)
    p.add_argument("--select_k", type=int, default=2000)
    p.add_argument("--early_stopping", type=int, default=10)
    p.add_argument("--n_samples", type=int, default=None)
    p.add_argument("--max_param_combos", type=int, default=None)

    p.add_argument("--tuning_strategy", choices=["nested", "once"], default="once")
    p.add_argument("--grid_preset", choices=["tiny", "fast", "current"], default="fast")

    p.add_argument(
        "--normalization",
        default="none",
        help="none, zscore, rank, or qbinN such as qbin50.",
    )

    p.add_argument("--xgb_threads", type=int, default=None)

    p.add_argument("--minimal_save", action="store_true", default=True)
    p.add_argument("--save_models", action="store_true", default=True)
    p.add_argument("--no_save_models", action="store_false", dest="save_models")

    p.add_argument("--save_final_model", action="store_true", default=False)
    p.add_argument("--no_save_final_model", action="store_false", dest="save_final_model")

    p.add_argument("--model_dir", default=None)
    p.add_argument("--metrics_dir", default=None)
    p.add_argument("--resume_check_dir", default=None)

    p.add_argument("--resume_skip", action="store_true", default=True)
    p.add_argument("--no_resume_skip", action="store_false", dest="resume_skip")

    p.add_argument("--compression", choices=["xz", "gzip", "zlib", "zstd", "none"], default="xz")

    p.add_argument("--output_mode", choices=["full_joblib", "minimal_h5"], default="minimal_h5")
    p.add_argument("--h5_compression_level", type=int, default=9)
    p.add_argument("--h5_model", choices=["best_loo", "final"], default="best_loo")

    p.add_argument("--force_gpu", action="store_true")
    p.add_argument("--no_feature_selection", action="store_true")

    p.add_argument("--output_suffix", default="loo_xgboost_none_fast_once_k2000")

    return p.parse_args()


def validate_normalization(normalization: str) -> str:
    normalization = (normalization or "none").lower().strip()

    if normalization in {"none", "zscore", "rank"}:
        return normalization

    if re.fullmatch(r"qbin[0-9]+", normalization):
        n_bins = int(normalization.replace("qbin", ""))
        if n_bins < 2:
            raise ValueError("qbin normalization requires at least 2 bins.")
        return normalization

    raise ValueError(
        f"Invalid normalization: {normalization}. "
        "Allowed: none, zscore, rank, or qbinN such as qbin50."
    )


def data_paths():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.abspath(os.path.join(current_dir, ".."))

    default_data_dir = os.path.join(project_dir, "..", "data")
    default_results_dir = os.path.join(project_dir, "..", "data", "results_xgb_loo")

    data_dir = os.path.abspath(os.getenv("DATA_DIR", default_data_dir))
    results_dir = os.path.abspath(os.getenv("RESULTS_DIR", default_results_dir))

    X_path = os.path.join(data_dir, "transcriptomics_filtered.joblib")
    Y_path = os.path.join(data_dir, "proteomics_matched.joblib")

    return project_dir, data_dir, results_dir, X_path, Y_path


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


def _compression_arg(kind: str):
    kind = (kind or "xz").lower()
    if kind == "none":
        return 0
    if kind == "xz":
        return ("xz", 9)
    if kind == "gzip":
        return ("gzip", 9)
    if kind == "zlib":
        return 9
    if kind == "zstd":
        return ("zstd", 19)
    raise ValueError(f"Unknown compression: {kind}")


def dump_compressed(obj, path: str, compression: str = "xz"):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    comp = _compression_arg(compression)

    try:
        dump(obj, path, compress=comp)
    except Exception:
        try:
            dump(obj, path, compress=("gzip", 9))
        except Exception:
            dump(obj, path, compress=9)


def safe_path_component(value: str) -> str:
    out = str(value).replace(os.sep, "_").replace("/", "_").replace("\\", "_")
    out = re.sub(r"[^A-Za-z0-9._+=@-]+", "_", out)
    return out


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


def expected_output_path(
    protein: str,
    *,
    output_mode: str,
    output_suffix: str,
    metrics_dir: str,
    model_dir: str,
    kind: str = "primary",
) -> str:
    safe_protein = safe_path_component(protein)

    if output_mode == "minimal_h5":
        return os.path.join(metrics_dir, f"{safe_protein}_{output_suffix}.h5")

    if kind == "predictions":
        return os.path.join(metrics_dir, f"{safe_protein}_{output_suffix}_predictions.joblib")
    if kind == "metrics":
        return os.path.join(metrics_dir, f"{safe_protein}_{output_suffix}_metrics.joblib")
    if kind == "best_model":
        return os.path.join(model_dir, f"{safe_protein}_{output_suffix}_best_loo_model.joblib")
    if kind == "final_model":
        return os.path.join(model_dir, f"{safe_protein}_{output_suffix}_final_all_samples_model.joblib")

    return os.path.join(metrics_dir, f"{safe_protein}_{output_suffix}_metrics.joblib")


def output_exists_for_resume(
    protein: str,
    *,
    output_mode: str,
    output_suffix: str,
    metrics_dir: str,
    model_dir: str,
    resume_check_dir: Optional[str],
) -> bool:
    check_metrics_dir = os.path.abspath(resume_check_dir or metrics_dir)
    check_model_dir = os.path.join(check_metrics_dir, "models")

    path = expected_output_path(
        protein,
        output_mode=output_mode,
        output_suffix=output_suffix,
        metrics_dir=check_metrics_dir,
        model_dir=check_model_dir,
        kind="primary",
    )

    return os.path.exists(path)


def save_result(
    result,
    *,
    protein: str,
    args,
    metrics_dir: str,
    model_dir: str,
):
    predictions_df = result["predictions"]
    metrics = result["metrics"]
    best_model_artifact = result.get("best_loo_model_artifact")
    final_model_artifact = result.get("final_model_artifact")

    if args.output_mode == "minimal_h5":
        h5_path = expected_output_path(
            protein,
            output_mode="minimal_h5",
            output_suffix=args.output_suffix,
            metrics_dir=metrics_dir,
            model_dir=model_dir,
        )
        print(f"Saving minimal compressed H5: {h5_path}")
        save_minimal_protein_h5(
            result,
            h5_path,
            protein_name=protein,
            model_choice=args.h5_model,
            compression_level=args.h5_compression_level,
        )
        return h5_path

    pred_path = expected_output_path(
        protein,
        output_mode="full_joblib",
        output_suffix=args.output_suffix,
        metrics_dir=metrics_dir,
        model_dir=model_dir,
        kind="predictions",
    )
    metrics_path = expected_output_path(
        protein,
        output_mode="full_joblib",
        output_suffix=args.output_suffix,
        metrics_dir=metrics_dir,
        model_dir=model_dir,
        kind="metrics",
    )

    print(f"Saving predictions: {pred_path}")
    dump_compressed(predictions_df, pred_path, args.compression)

    print(f"Saving metrics: {metrics_path}")
    dump_compressed(metrics, metrics_path, args.compression)

    if args.save_models and best_model_artifact is not None:
        best_path = expected_output_path(
            protein,
            output_mode="full_joblib",
            output_suffix=args.output_suffix,
            metrics_dir=metrics_dir,
            model_dir=model_dir,
            kind="best_model",
        )
        print(f"Saving best LOO model artifact: {best_path}")
        dump_compressed(best_model_artifact, best_path, args.compression)

    if args.save_models and args.save_final_model and final_model_artifact is not None:
        final_path = expected_output_path(
            protein,
            output_mode="full_joblib",
            output_suffix=args.output_suffix,
            metrics_dir=metrics_dir,
            model_dir=model_dir,
            kind="final_model",
        )
        print(f"Saving final all-samples model artifact: {final_path}")
        dump_compressed(final_model_artifact, final_path, args.compression)

    return metrics_path


def print_summary(metrics: dict):
    print("Summary metrics:")
    for k, v in metrics.items():
        if k not in {"best_loo_model_summary", "final_model_summary"}:
            print(f"  {k}: {v}")


def run_one_protein(
    *,
    protein: str,
    X: pd.DataFrame,
    Y: pd.DataFrame,
    param_grid: dict,
    args,
    metrics_dir: str,
    model_dir: str,
    use_gpu: bool,
):
    if args.resume_skip and output_exists_for_resume(
        protein,
        output_mode=args.output_mode,
        output_suffix=args.output_suffix,
        metrics_dir=metrics_dir,
        model_dir=model_dir,
        resume_check_dir=args.resume_check_dir,
    ):
        print(f"Skipping existing output for protein: {protein}")
        return "skipped"

    print("=" * 100)
    print(
        f"Running LOO XGBoost for protein '{protein}' | "
        f"grid={args.grid_preset} | tuning={args.tuning_strategy} | "
        f"normalization={args.normalization} | select_k={args.select_k} | gpu={use_gpu}"
    )

    try:
        result = cv_pipeline_loo(
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
            save_models=args.save_models,
            save_final_model=args.save_final_model,
        )

        saved_path = save_result(
            result,
            protein=protein,
            args=args,
            metrics_dir=metrics_dir,
            model_dir=model_dir,
        )

        print_summary(result["metrics"])
        print(f"Saved output: {saved_path}")
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

    if args.output_mode == "minimal_h5":
        args.minimal_save = True
        args.save_models = True

        if args.h5_model == "best_loo":
            args.save_final_model = False
        elif args.h5_model == "final":
            args.save_final_model = True

    project_dir, data_dir, default_results_dir, X_path, Y_path = data_paths()

    metrics_dir = os.path.abspath(args.metrics_dir or default_results_dir)
    model_dir = os.path.abspath(args.model_dir or os.path.join(metrics_dir, "models"))

    os.makedirs(metrics_dir, exist_ok=True)
    os.makedirs(model_dir, exist_ok=True)

    print("Loading data...")
    print(f"PROJECT_DIR: {project_dir}")
    print(f"DATA_DIR: {data_dir}")
    print(f"METRICS_DIR: {metrics_dir}")
    print(f"MODEL_DIR: {model_dir}")
    print(f"X path: {X_path}")
    print(f"Y path: {Y_path}")

    try:
        X = load(X_path)
        Y = load(Y_path)
    except Exception as e:
        print(f"Error loading data ({X_path} / {Y_path}): {e}", file=sys.stderr)
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
    print(f"Selected protein range preview: {proteins_to_run[:5]}")

    if not proteins_to_run:
        print("No proteins selected. Exiting.")
        return

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
            model_dir=model_dir,
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