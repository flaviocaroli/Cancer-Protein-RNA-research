# new_pipeline.py
# Runs a *group* of proteins per call (e.g., 511 proteins).
# Writes per-protein best files and a concatenated batch parquet into data/results/.
# Compatible with Python 3.8+.

import argparse
import os
import sys
from pathlib import Path

from joblib import load

from new_functions import run_protein_group


def read_protein_list(protein_list_path):
    """Return proteins, ignoring blank lines and comments (#...)."""
    with open(protein_list_path, "r") as f:
        return [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]


def main():
    parser = argparse.ArgumentParser()

    # chunking
    parser.add_argument("--chunk_id", type=int, default=None, help="0-based chunk id. start_idx = chunk_id*group_size")
    parser.add_argument("--start_idx", type=int, default=None, help="0-based start index into protein_list.txt")

    # NEW: group_size can be computed from target_chunks
    parser.add_argument(
        "--group_size",
        type=int,
        default=None,
        help="How many proteins to run in this call (overrides --target_chunks).",
    )
    parser.add_argument(
        "--target_chunks",
        type=int,
        default=10,
        help="If --group_size is not set, group_size = ceil(total_proteins/target_chunks).",
    )

    # pipeline
    parser.add_argument(
        "--pipeline",
        type=str,
        default="fixed_percentile",
        choices=["fixed_percentile", "fixed_p", "union"],
    )
    parser.add_argument(
        "--full_output",
        action="store_true",
        help="If set, save full LOO outputs (large). Default: best-only.",
    )

    # knobs
    parser.add_argument("--percentile", type=int, default=5)
    parser.add_argument("--p_value", type=float, default=0.05)
    parser.add_argument("--n_splits", type=int, default=5)
    parser.add_argument("--enet_max_iter", type=int, default=5000)
    parser.add_argument("--enet_tol", type=float, default=1e-3)
    parser.add_argument("--enet_selection", type=str, default="random")
    parser.add_argument("--fit_intercept", action="store_true", help="Fit ElasticNet intercept. Default False.")
    parser.add_argument("--no_zscore", action="store_true", help="Disable z-scoring. Default is z-score (no leakage).")

    # paths
    parser.add_argument("--data_dir", type=str, default=None, help="Path to data folder (joblibs + protein_list.txt).")
    parser.add_argument("--results_dir", type=str, default=None, help="Where to write outputs. Default: <data_dir>/results")
    parser.add_argument("--protein_list", type=str, default=None, help="Path to protein_list.txt")

    # output controls
    parser.add_argument("--no_batch_file", action="store_true", help="Do not save concatenated batch parquet.")
    parser.add_argument(
        "--consume_list",
        action="store_true",
        help="DANGEROUS with parallel jobs: delete processed lines from protein_list.txt (sequential use only).",
    )

    args = parser.parse_args()

    # Resolve dirs (assumes src/ is next to ../data/)
    if args.data_dir is None:
        here = Path(__file__).resolve().parent
        data_dir = (here.parent / "data").resolve()
    else:
        data_dir = Path(args.data_dir).resolve()

    results_dir = Path(args.results_dir).resolve() if args.results_dir else (data_dir / "results").resolve()
    results_dir.mkdir(parents=True, exist_ok=True)

    protein_list_path = args.protein_list or str((data_dir / "protein_list.txt").resolve())

    transcript_path = str((data_dir / "transcriptomics_common.joblib").resolve())
    protein_path = str((data_dir / "proteomics_matched.joblib").resolve())

    if not os.path.exists(transcript_path):
        print("Error: missing transcriptomics joblib at {}".format(transcript_path))
        sys.exit(1)
    if not os.path.exists(protein_path):
        print("Error: missing proteomics joblib at {}".format(protein_path))
        sys.exit(1)
    if not os.path.exists(protein_list_path):
        print("Error: missing protein list at {}".format(protein_list_path))
        sys.exit(1)

    print("Loading data...")
    X = load(transcript_path)
    Y = load(protein_path)

    # Load protein list and compute total
    all_proteins = read_protein_list(protein_list_path)
    total = len(all_proteins)
    if total == 0:
        print("protein_list.txt is empty (after removing comments/blanks). Nothing to do.")
        return

    # Decide group_size
    if args.group_size is not None:
        group_size = int(args.group_size)
    else:
        target_chunks = int(args.target_chunks) if args.target_chunks else 10
        group_size = (total + target_chunks - 1) // target_chunks  # ceil(total/target_chunks)
        group_size = max(group_size, 1)

    # Determine start index
    if args.start_idx is not None:
        start_idx = int(args.start_idx)
    elif args.chunk_id is not None:
        start_idx = int(args.chunk_id) * group_size
    else:
        start_idx = 0

    end_idx = min(start_idx + group_size, total)
    proteins = all_proteins[start_idx:end_idx]

    if len(proteins) == 0:
        print("No proteins selected (start_idx={}, group_size={}, total={}). Nothing to do.".format(start_idx, group_size, total))
        return

    print("Selected proteins: {}..{} (count={}, total_list={})".format(start_idx, end_idx - 1, len(proteins), total))
    print("Pipeline={}, save_only_best={}, results_dir={}".format(args.pipeline, (not args.full_output), str(results_dir)))
    print("group_size={} (arg_group_size={}, target_chunks={})".format(group_size, args.group_size, args.target_chunks))

    # grid
    param_grid = {"alpha": [0.1, 0.5, 1.0], "l1_ratio": [0.1, 0.5, 0.9]}

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
        zscore=(not args.no_zscore),
        zscore_eps=1e-8,
        save_only_best=(not args.full_output),
        output_dir=str(results_dir),
    )

    # Save batch file
    if not args.no_batch_file:
        batch_name = "batch_{}_{}_{}.parquet".format(start_idx, end_idx - 1, args.pipeline)
        batch_path = results_dir / batch_name
        try:
            out.to_parquet(batch_path, compression="snappy", index=False)
        except Exception:
            out.to_csv(batch_path.with_suffix(".csv"), index=False)
        print("Saved batch results to:", str(batch_path))

    if args.consume_list:
        # Since we parsed proteins ignoring comments/blanks, consume_list should operate on that same filtered list.
        kept_prots = all_proteins[:start_idx] + all_proteins[end_idx:]

        tmp = protein_list_path + ".tmp"
        with open(tmp, "w") as f:
            for p in kept_prots:
                f.write(p.strip() + "\n")
        os.replace(tmp, protein_list_path)
        print("Consumed {} proteins from protein_list.txt (deleted indices {}..{})".format(len(proteins), start_idx, end_idx - 1))

    print("Done.")


if __name__ == "__main__":
    main()
