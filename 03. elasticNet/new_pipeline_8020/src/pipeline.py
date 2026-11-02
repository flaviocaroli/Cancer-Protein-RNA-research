#!/usr/bin/env python3
"""Chunk-aware 80/20 ElasticNet pipeline.

Reads:
  data/transcriptomics_filtered.joblib
  data/proteomics_matched.joblib

Writes:
  results_enet/<PROTEIN>_enet_8020_<normalization>_<grid>_k<select_k>.h5
"""

import argparse
import gc
import os
import sys
from typing import List, Optional, Tuple

import pandas as pd
from joblib import load

from functions import (
    detect_protein_axis,
    run_enet_8020_for_protein,
    safe_path_component,
    validate_normalization,
)


def parse_args():
    p = argparse.ArgumentParser(description="Run 80/20 ElasticNet for one protein or a protein chunk.")

    p.add_argument("--protein", default=None, type=str)
    p.add_argument("--protein_start", type=int, default=None, help="1-based inclusive protein start index.")
    p.add_argument("--protein_end", type=int, default=None, help="1-based inclusive protein end index.")

    p.add_argument("--normalization", default="zscore", help="none, zscore, ranking, binning. Aliases: rank, qbin50.")
    p.add_argument("--n_bins", type=int, default=50)
    p.add_argument("--select_k", type=int, default=2000)
    p.add_argument("--test_size", type=float, default=0.2)
    p.add_argument("--splits", type=int, default=3)
    p.add_argument("--max_param_combos", type=int, default=8)
    p.add_argument("--grid_preset", choices=["tiny", "fast", "current"], default="fast")

    p.add_argument("--enet_max_iter", type=int, default=10000)
    p.add_argument("--enet_tol", type=float, default=1e-3)
    p.add_argument("--enet_selection", choices=["cyclic", "random"], default="random")
    p.add_argument("--fit_intercept", action="store_true")
    p.add_argument("--n_samples", type=int, default=None)

    p.add_argument("--metrics_dir", default=None)
    p.add_argument("--resume_check_dir", default=None)
    p.add_argument("--resume_skip", action="store_true", default=True)
    p.add_argument("--no_resume_skip", action="store_false", dest="resume_skip")
    p.add_argument("--h5_compression_level", type=int, default=9)
    p.add_argument("--output_suffix", default=None)

    return p.parse_args()


def data_paths():
    here = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.abspath(os.path.join(here, ".."))
    data_dir = os.path.abspath(os.getenv("DATA_DIR", os.path.join(project_dir, "..", "data")))
    results_dir = os.path.abspath(os.getenv("RESULTS_DIR", os.path.join(project_dir, "..", "data", "results_enet")))
    x_path = os.path.join(data_dir, "transcriptomics_filtered.joblib")
    y_path = os.path.join(data_dir, "proteomics_matched.joblib")
    return project_dir, data_dir, results_dir, x_path, y_path


def build_param_grid(preset: str = "fast"):
    if preset == "tiny":
        return {
            "alpha": [0.1, 1.0],
            "l1_ratio": [0.1, 0.9],
        }
    if preset == "current":
        return {
            "alpha": [0.0001, 0.001, 0.01, 0.1, 0.5, 1.0, 5.0, 10.0],
            "l1_ratio": [0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95],
        }
    return {
        "alpha": [0.001, 0.01, 0.1, 1.0],
        "l1_ratio": [0.1, 0.5, 0.9],
    }


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
        raise ValueError("--protein_start must be >= 1")
    if end < start:
        raise ValueError("--protein_end must be >= --protein_start")
    if start > n:
        raise ValueError(f"--protein_start={start} is beyond number of proteins={n}")
    end = min(end, n)
    return all_proteins[start - 1:end]


def output_path_for(protein: str, metrics_dir: str, output_suffix: str) -> str:
    return os.path.join(metrics_dir, f"{safe_path_component(protein)}_{output_suffix}.h5")


def exists_for_resume(protein: str, metrics_dir: str, resume_check_dir: Optional[str], output_suffix: str) -> bool:
    check_dir = os.path.abspath(resume_check_dir or metrics_dir)
    return os.path.exists(output_path_for(protein, check_dir, output_suffix))


def main():
    args = parse_args()
    args.normalization = validate_normalization(args.normalization)

    project_dir, data_dir, default_results_dir, x_path, y_path = data_paths()
    metrics_dir = os.path.abspath(args.metrics_dir or default_results_dir)
    os.makedirs(metrics_dir, exist_ok=True)

    if args.output_suffix is None:
        norm_label = args.normalization
        if args.normalization == "binning":
            norm_label = f"binning{int(args.n_bins)}"
        args.output_suffix = f"enet_8020_{norm_label}_{args.grid_preset}_k{int(args.select_k)}"

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
        print(f"ERROR loading data: {type(e).__name__}: {e}", file=sys.stderr)
        sys.exit(1)

    X.index.name = "samples"
    X.columns.name = "genes"

    all_proteins, orientation, index_overlap, column_overlap = detect_protein_axis(X, Y)
    print(f"Transcriptomics shape: {getattr(X, 'shape', None)}")
    print(f"Proteomics shape: {getattr(Y, 'shape', None)}")
    print(f"Detected orientation: {orientation}")
    print(f"Sample overlap with Y.index: {index_overlap}")
    print(f"Sample overlap with Y.columns: {column_overlap}")
    print(f"Proteins detected: {len(all_proteins)}")
    print(f"First proteins: {all_proteins[:5]}")

    try:
        proteins = select_proteins(
            all_proteins,
            protein=args.protein,
            protein_start=args.protein_start,
            protein_end=args.protein_end,
        )
    except Exception as e:
        print(f"ERROR selecting proteins: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Selected proteins: {len(proteins)}")
    print(f"Preview: {proteins[:5]}")
    print(f"Normalization: {args.normalization}")
    print(f"Output suffix: {args.output_suffix}")

    param_grid = build_param_grid(args.grid_preset)
    statuses = {"done": 0, "skipped": 0, "failed": 0}

    for i, protein in enumerate(proteins, start=1):
        print("=" * 100)
        print(f"Protein {i}/{len(proteins)}: {protein}")
        out_path = output_path_for(protein, metrics_dir, args.output_suffix)

        if args.resume_skip and exists_for_resume(protein, metrics_dir, args.resume_check_dir, args.output_suffix):
            print(f"Skipping existing: {out_path}")
            statuses["skipped"] += 1
            continue

        try:
            result = run_enet_8020_for_protein(
                X,
                Y,
                protein,
                param_grid=param_grid,
                normalization=args.normalization,
                n_bins=args.n_bins,
                select_k=args.select_k,
                test_size=args.test_size,
                n_splits=args.splits,
                max_param_combos=args.max_param_combos,
                enet_max_iter=args.enet_max_iter,
                enet_tol=args.enet_tol,
                enet_selection=args.enet_selection,
                fit_intercept=args.fit_intercept,
                n_samples=args.n_samples,
                h5_path=out_path,
                h5_compression_level=args.h5_compression_level,
            )
            m = result["metrics"]
            print("Saved:", out_path)
            print(f"  test_mse: {m.get('test_mse')}")
            print(f"  test_pearson_corr: {m.get('test_pearson_corr')}")
            print(f"  test_cosine_similarity: {m.get('test_cosine_similarity')}")
            print(f"  selected_features: {len(result.get('selected_features', []))}")
            print(f"  best_params: {result.get('best_params')}")
            statuses["done"] += 1
        except Exception as e:
            print(f"ERROR for protein {protein}: {type(e).__name__}: {e}", file=sys.stderr)
            statuses["failed"] += 1
        finally:
            gc.collect()

    print("=" * 100)
    print("Chunk finished.")
    print(f"Done: {statuses['done']}")
    print(f"Skipped: {statuses['skipped']}")
    print(f"Failed: {statuses['failed']}")
    if statuses["failed"] > 0:
        sys.exit(2)


if __name__ == "__main__":
    main()
