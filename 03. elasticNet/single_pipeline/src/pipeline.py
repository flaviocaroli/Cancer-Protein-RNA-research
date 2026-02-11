# pipeline.py
import os
import sys
import argparse
import pandas as pd
from joblib import load, dump

from pipeline_functions import (
    cv_pipeline_fixed_percentile,
    cv_pipeline_fixed_p,
    cv_union,
)


def select_pipeline(pipeline_name):
    pipelines = {
        "fixed_percentile": cv_pipeline_fixed_percentile,
        "fixed_p": cv_pipeline_fixed_p,
        "union": cv_union,
    }
    if pipeline_name not in pipelines:
        raise ValueError("Invalid pipeline: {}".format(pipeline_name))
    return pipelines[pipeline_name]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein", type=str, required=True)
    parser.add_argument(
        "--pipeline",
        type=str,
        default="fixed_percentile",
        choices=["fixed_percentile", "fixed_p", "union"],
        help="Which pipeline to run. Default: fixed_percentile",
    )
    parser.add_argument("--model", type=str, default="elasticnet")

    parser.add_argument(
        "--full_output",
        action="store_true",
        help="If set, outputs full LOO results + consensus model + loo_stats. Otherwise saves only {protein}_best.parquet.",
    )

    # Optional knobs
    parser.add_argument("--percentile", type=int, default=5)
    parser.add_argument("--p_value", type=float, default=0.05)
    parser.add_argument("--n_splits", type=int, default=3)

    args = parser.parse_args()

    # sanitize protein name (important for SLURM + sed + CRLF issues)
    args.protein = args.protein.strip().replace("\r", "")

    save_only_best = (not args.full_output)  # default True

    # Paths
    data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))
    results_dir = os.path.join(data_dir, "results")
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    transcript_path = os.path.join(data_dir, "transcriptomics_filtered.joblib")
    protein_path = os.path.join(data_dir, "proteomics_matched.joblib")

    print("Loading Data...")
    print("data_dir:", data_dir)
    print("results_dir:", results_dir)
    print("transcript_path:", transcript_path, "exists:", os.path.exists(transcript_path))
    print("protein_path:", protein_path, "exists:", os.path.exists(protein_path))

    X = load(transcript_path)
    Y = load(protein_path)

    print("Running protein:", repr(args.protein))
    print("Proteomics shape:", getattr(Y, "shape", None))

    if args.protein not in Y.columns:
        # print a few close matches to debug typos/whitespace
        print("Error: {} not found in Y.columns.".format(repr(args.protein)))
        print("First 10 Y columns:", list(Y.columns[:10]))
        sys.exit(1)

    param_grid = {"alpha": [0.1, 0.5, 1.0], "l1_ratio": [0.1, 0.5, 0.9]}
    pipeline_func = select_pipeline(args.pipeline)

    print("Starting LOO for {} using pipeline={} ...".format(args.protein, args.pipeline))
    print("save_only_best =", save_only_best)
    print("Outputs will be saved to:", results_dir)

    # ---------------------------------------------------------
    # Run pipeline
    # ---------------------------------------------------------
    if args.pipeline == "fixed_percentile":
        raw_results = pipeline_func(
            None,
            param_grid,
            X,
            Y,
            args.protein,
            feature_selector="anova",
            n_splits=args.n_splits,
            percentile=args.percentile,
            save_only_best=save_only_best,
            output_dir=results_dir,
        )

    elif args.pipeline == "fixed_p":
        raw_results = pipeline_func(
            None,
            param_grid,
            X,
            Y,
            args.protein,
            feature_selector="anova",
            n_splits=args.n_splits,
            p_value=args.p_value,
            save_only_best=save_only_best,
            output_dir=results_dir,
        )

    else:  # "union"
        raw_results = pipeline_func(
            None,
            param_grid,
            X,
            Y,
            args.protein,
            feature_selector="anova",
            percentile=args.percentile,
            n_splits_stage1=args.n_splits,
            n_splits_stage2=args.n_splits,
            save_only_best=save_only_best,
            output_dir=results_dir,
        )

    # ---------------------------------------------------------
    # Output logic
    # ---------------------------------------------------------
    if save_only_best:
        # pipeline_functions already wrote {results_dir}/{protein}_best.parquet (or CSV fallback)
        out_best = os.path.join(results_dir, "{}_best.parquet".format(args.protein))
        out_best_csv = os.path.join(results_dir, "{}_best.csv".format(args.protein))

        if os.path.exists(out_best):
            print("Done. Saved:", out_best)
        elif os.path.exists(out_best_csv):
            print("Done. Parquet unavailable; saved:", out_best_csv)
        else:
            print("Done, but best output file was not found. Check results_dir permissions/path.")
        return

    # Full output mode: consensus + loo_stats
    coefs_df = pd.DataFrame(raw_results["features_coef_dict"].tolist()).fillna(0.0)
    mean_weights = coefs_df.mean(axis=0)
    stability = (coefs_df != 0).mean(axis=0)

    consensus_df = pd.DataFrame(
        {"gene": mean_weights.index, "weight": mean_weights.values, "stability": stability.values}
    ).query("weight != 0")

    output_path = os.path.join(results_dir, "{}_consensus_model.parquet".format(args.protein))
    consensus_df.to_parquet(output_path, compression="snappy")

    perf_path = os.path.join(results_dir, "{}_loo_stats.joblib".format(args.protein))
    dump(raw_results[["LOO sample", "true", "predicted"]], perf_path)

    print("Process complete. Consensus saved to {}".format(output_path))
    print("LOO stats saved to {}".format(perf_path))


if __name__ == "__main__":
    main()
