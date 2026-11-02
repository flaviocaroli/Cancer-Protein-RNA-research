# new_pipeline.py
# Runs a *group* of proteins per call.
# Writes per-protein best files and a concatenated batch parquet.
# Compatible with Python 3.8+.

import argparse
import os
import sys
from pathlib import Path

import json

import numpy as np
import pandas as pd
from joblib import load

# Import both standard and normalization runners
from new_functions import best_only_row_for_protein, run_protein_group
from norm_functions import _best_only_row_from_tracking, _safe_to_parquet_or_csv, run_protein_group_norm


def read_protein_list(protein_list_path):
    """Return proteins, ignoring blank lines and comments (#...)."""
    with open(protein_list_path, "r") as f:
        return [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]


def _jsonify_dict_columns(df):
    """Dict/list columns -> JSON strings so batch files stay compact and CSV/parquet round-trip."""
    df = df.copy()
    for col in df.columns:
        if df[col].map(lambda v: isinstance(v, (dict, list, np.ndarray))).any():
            df[col] = df[col].map(lambda v: json.dumps(v.tolist() if isinstance(v, np.ndarray) else v)
                                  if isinstance(v, (dict, list, np.ndarray)) else v)
    return df


def _write_full_output_files(out, X, Y, args, results_dir):
    """--full_output: write <P>_LOO.parquet (one row per LOO sample) AND the usual <P>_best.parquet,
    so the same run feeds both 03/model_results (enet_all_results.h5) and 06 validation."""
    best_rows = []
    for prot, raw in out.groupby("protein", sort=False):
        raw = raw.dropna(subset=["LOO sample"]).reset_index(drop=True)
        if raw.empty:
            continue
        loo = _jsonify_dict_columns(raw)
        loo["normalization"] = args.normalization
        _safe_to_parquet_or_csv(loo, Path(results_dir) / f"{prot}_LOO.parquet")
        cols = ["LOO sample", "true", "predicted", "params_dict", "features_coef_dict"]
        if args.normalization in ("binning", "ranking"):
            err = raw["true"].astype(float) - raw["predicted"].astype(float)
            i = int(err.abs().values.argmin())
            best = _best_only_row_from_tracking(
                protein_name=prot, sample_ids=raw["LOO sample"].astype(str).to_numpy(), best_i=i,
                best_abs_error=float(abs(err.iloc[i])), total_sse=float((err ** 2).sum()),
                best_sq_error=float(err.iloc[i] ** 2), params_dict=dict(raw.loc[i, "params_dict"]),
                coef_dict=dict(raw.loc[i, "features_coef_dict"]), fit_intercept=bool(args.fit_intercept),
                output_dir=str(results_dir), normalization="RelativeRepresentation")
        else:
            best = best_only_row_for_protein(
                raw[cols], X, Y[prot], prot, zscore=(args.normalization == "zscore"), zscore_eps=1e-8,
                fit_intercept=bool(args.fit_intercept), output_dir=str(results_dir))
        best_rows.append(best)
    print(f"full_output: wrote {len(best_rows)} *_LOO and *_best files")
    return out


def main():
    parser = argparse.ArgumentParser()

    # ---------------------------
    # 1. Chunking / Grouping
    # ---------------------------
    parser.add_argument("--chunk_id", type=int, default=None, help="0-based chunk id.")
    parser.add_argument("--start_idx", type=int, default=None, help="0-based start index.")
    parser.add_argument("--group_size", type=int, default=None, help="Proteins per run.")
    parser.add_argument("--target_chunks", type=int, default=10, help="Target number of chunks if group_size not set.")

    # ---------------------------
    # 2. Pipeline Strategy (Feature Selection logic)
    # ---------------------------
    parser.add_argument(
        "--pipeline",
        type=str,
        default="fixed_percentile",
        choices=["fixed_percentile", "fixed_p", "union", "naive_cluster"],
        help="The feature selection strategy to use.",
    )

    # ---------------------------
    # 3. Normalization Method (Data Transformation)
    # ---------------------------
    parser.add_argument(
        "--normalization",
        type=str,
        default="zscore",
        choices=["zscore", "binning", "ranking", "none"],
        help="The data normalization method. 'zscore' is standard. 'binning'/'ranking' use Relative Representations.",
    )

    # ---------------------------
    # 4. Hyperparameters & Settings
    # ---------------------------
    parser.add_argument("--percentile", type=int, default=5)
    parser.add_argument("--p_value", type=float, default=0.05)
    parser.add_argument("--n_splits", type=int, default=5)
    parser.add_argument("--enet_max_iter", type=int, default=5000)
    parser.add_argument("--enet_tol", type=float, default=1e-3)
    parser.add_argument("--enet_selection", type=str, default="random")
    parser.add_argument("--fit_intercept", action="store_true", help="Fit ElasticNet intercept.")
    
    # Binning specific
    parser.add_argument("--n_bins", type=int, default=50, help="Bins for --normalization binning.")

    # Naive Cluster specific
    parser.add_argument("--clusters_csv", type=str, default=None, help="Required for --pipeline naive_cluster.")
    parser.add_argument("--naive_n_samples", type=int, default=None, help="Debug limit for naive_cluster.")

    # Output / Paths
    parser.add_argument("--full_output", action="store_true", help="Save full LOO outputs.")
    parser.add_argument("--data_dir", type=str, default=None)
    parser.add_argument("--rna_file", type=str, default="transcriptomics_filtered.joblib",
                        help="RNA matrix inside data_dir. Default transcriptomics_filtered.joblib (same features as "
                             "the 80/20 and XGB pipelines). Use transcriptomics_common.joblib only to reproduce the "
                             "common-gene sensitivity run.")
    parser.add_argument("--results_dir", type=str, default=None)
    parser.add_argument("--protein_list", type=str, default=None)
    parser.add_argument("--no_batch_file", action="store_true")
    parser.add_argument("--consume_list", action="store_true")

    args = parser.parse_args()

    # ---------------------------
    # Setup Directories
    # ---------------------------
    if args.data_dir is None:
        here = Path(__file__).resolve().parent
        data_dir = (here.parent.parent / "data").resolve()   # single 03. elasticNet/data folder
    else:
        data_dir = Path(args.data_dir).resolve()

    results_dir = Path(args.results_dir).resolve() if args.results_dir else (data_dir / "loo" / f"results_{args.normalization}").resolve()
    results_dir.mkdir(parents=True, exist_ok=True)

    protein_list_path = args.protein_list or str((data_dir / "protein_list.txt").resolve())
    transcript_path = str((data_dir / args.rna_file).resolve())
    protein_path = str((data_dir / "proteomics_matched.joblib").resolve())

    # Path Validations
    if not os.path.exists(transcript_path):
        sys.exit(f"Error: missing {transcript_path}")
    if not os.path.exists(protein_path):
        sys.exit(f"Error: missing {protein_path}")
    if not os.path.exists(protein_list_path):
        sys.exit(f"Error: missing {protein_list_path}")

    # Pipeline Specific Validations
    if args.pipeline == "naive_cluster":
        if args.clusters_csv is None:
            sys.exit("Error: --pipeline naive_cluster requires --clusters_csv.")
        if not os.path.exists(args.clusters_csv):
            sys.exit(f"Error: missing {args.clusters_csv}")

    # Binning/Ranking Validity Check
    # norm_functions.py currently only supports fixed_percentile logic
    if args.normalization in ["binning", "ranking"]:
        if args.pipeline != "fixed_percentile":
            sys.exit(f"Error: --normalization {args.normalization} currently only supports --pipeline fixed_percentile.")

    print("Loading data...")
    X = load(transcript_path)
    Y = load(protein_path)

    all_proteins = read_protein_list(protein_list_path)
    total = len(all_proteins)
    if total == 0:
        print("Protein list is empty.")
        return

    # ---------------------------
    # Grouping Logic
    # ---------------------------
    if args.group_size is not None:
        group_size = int(args.group_size)
    else:
        target_chunks = int(args.target_chunks) if args.target_chunks else 10
        group_size = max((total + target_chunks - 1) // target_chunks, 1)

    if args.start_idx is not None:
        start_idx = int(args.start_idx)
    elif args.chunk_id is not None:
        start_idx = int(args.chunk_id) * group_size
    else:
        start_idx = 0

    end_idx = min(start_idx + group_size, total)
    proteins = all_proteins[start_idx:end_idx]

    if not proteins:
        print("No proteins selected.")
        return

    print(f"Proteins: {start_idx}..{end_idx - 1} ({len(proteins)}/{total})")
    print(f"Pipeline: {args.pipeline} | Normalization: {args.normalization}")

    # ---------------------------
    # Execution Dispatch
    # ---------------------------
    param_grid = {"alpha": [0.1, 0.5, 1.0], "l1_ratio": [0.1, 0.5, 0.9]}

    if args.normalization in ["binning", "ranking"]:
        # Use norm_functions
        # Note: We pass args.normalization as the 'pipeline' arg to run_protein_group_norm 
        # because that function expects 'binning' or 'ranking' as the string key.
        out = run_protein_group_norm(
            proteins=proteins,
            pipeline=args.normalization,  # 'binning' or 'ranking'
            X=X,
            Y=Y,
            param_grid=param_grid,
            feature_selector="anova",
            percentile=args.percentile,
            n_splits=args.n_splits,
            n_bins=args.n_bins,
            enet_max_iter=args.enet_max_iter,
            enet_tol=args.enet_tol,
            enet_selection=args.enet_selection,
            fit_intercept=bool(args.fit_intercept),
            save_only_best=(not args.full_output),
            output_dir=str(results_dir),
        )
    else:
        # Use new_functions (Standard / Z-score)
        # Map normalization arg to boolean
        use_zscore = (args.normalization == "zscore")
        
        out = run_protein_group(
            proteins=proteins,
            pipeline=args.pipeline,
            X=X,
            Y=Y,
            param_grid=param_grid,
            percentile=args.percentile,
            p_value=args.p_value,
            n_splits=args.n_splits,
            n_splits_stage1=args.n_splits,
            n_splits_stage2=args.n_splits,
            enet_max_iter=args.enet_max_iter,
            enet_tol=args.enet_tol,
            enet_selection=args.enet_selection,
            fit_intercept=bool(args.fit_intercept),
            zscore=use_zscore,
            zscore_eps=1e-8,
            save_only_best=(not args.full_output),
            output_dir=str(results_dir),
            clusters_csv_path=args.clusters_csv,
            naive_n_samples=args.naive_n_samples,
        )

    # ---------------------------
    # Saving Results
    # ---------------------------
    if args.full_output and not out.empty and "LOO sample" in out.columns:
        out = _write_full_output_files(out, X, Y, args, results_dir)

    out = _jsonify_dict_columns(out)
    if not args.no_batch_file and not out.empty:
        # Include normalization in filename for clarity
        batch_name = f"batch_{start_idx}_{end_idx-1}_{args.pipeline}_{args.normalization}.parquet"
        batch_path = results_dir / batch_name
        try:
            out.to_parquet(batch_path, compression="snappy", index=False)
        except Exception:
            out.to_csv(batch_path.with_suffix(".csv"), index=False)
        print(f"Saved batch: {batch_path}")

    if args.consume_list:
        kept_prots = all_proteins[:start_idx] + all_proteins[end_idx:]
        tmp = protein_list_path + ".tmp"
        with open(tmp, "w") as f:
            for p in kept_prots:
                f.write(p + "\n")
        os.replace(tmp, protein_list_path)
        print(f"Consumed {len(proteins)} proteins.")

    print("Done.")


if __name__ == "__main__":
    main()