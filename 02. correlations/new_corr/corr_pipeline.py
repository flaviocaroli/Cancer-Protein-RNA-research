# corr_pipeline.py  (UPDATED: assumes data/ lives INSIDE new_corr/)
#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path

from joblib import load
from corr_functions import run_correlation_group


def read_protein_list(protein_list_path: str):
    """Return proteins, ignoring blank lines and comments (#...)."""
    with open(protein_list_path, "r") as f:
        return [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]


def main():
    parser = argparse.ArgumentParser(description="Run correlation for a chunk (group) of proteins.")

    # chunking
    parser.add_argument("--chunk_id", type=int, default=None, help="0-based chunk id. start_idx = chunk_id*group_size")
    parser.add_argument("--start_idx", type=int, default=None, help="0-based start index into protein_list.txt")

    parser.add_argument("--group_size", type=int, default=None, help="How many proteins to run in this call.")
    parser.add_argument(
        "--target_chunks",
        type=int,
        default=10,
        help="If --group_size is not set, group_size = ceil(total_proteins/target_chunks).",
    )

    # correlation knobs
    parser.add_argument("--p_threshold", type=float, default=0.05)
    parser.add_argument("--no_bootstrap", action="store_true", help="Disable bootstrap; use parametric p-value.")
    parser.add_argument("--n_bootstraps", type=int, default=100)
    parser.add_argument("--random_state", type=int, default=42)

    # paths
    parser.add_argument(
        "--data_dir",
        type=str,
        default=None,
        help="Path to data folder (joblibs + protein_list.txt). Default: <this_script_dir>/data",
    )
    parser.add_argument(
        "--results_dir",
        type=str,
        default=None,
        help="Where to write outputs. Default: <data_dir>/results",
    )
    parser.add_argument("--protein_list", type=str, default=None, help="Path to protein_list.txt")

    # I/O controls
    parser.add_argument("--per_protein_format", choices=["parquet", "csv"], default="parquet")
    parser.add_argument("--no_batch_file", action="store_true", help="Do not save batch summary parquet/csv.")

    args = parser.parse_args()

    # ✅ NEW DEFAULT: data_dir is inside new_corr/
    here = Path(__file__).resolve().parent
    data_dir = Path(args.data_dir).resolve() if args.data_dir else (here / "data").resolve()

    # ✅ NEW DEFAULT: results inside data/results
    results_dir = Path(args.results_dir).resolve() if args.results_dir else (data_dir / "results").resolve()
    results_dir.mkdir(parents=True, exist_ok=True)

    protein_list_path = args.protein_list or str((data_dir / "protein_list.txt").resolve())
    transcript_path = str((data_dir / "transcriptomics_filtered.joblib").resolve())
    protein_path = str((data_dir / "proteomics_matched.joblib").resolve())

    if not os.path.exists(transcript_path):
        print(f"Error: missing transcriptomics joblib at {transcript_path}")
        sys.exit(1)
    if not os.path.exists(protein_path):
        print(f"Error: missing proteomics joblib at {protein_path}")
        sys.exit(1)
    if not os.path.exists(protein_list_path):
        print(f"Error: missing protein list at {protein_list_path}")
        sys.exit(1)

    print("Loading data...")
    transcriptomics = load(transcript_path)
    proteomics = load(protein_path)

    transcriptomics.index.name = "samples"
    proteomics.index.name = "samples"
    transcriptomics.columns.name = "genes"
    proteomics.columns.name = "proteins"

    all_proteins = read_protein_list(protein_list_path)
    total = len(all_proteins)
    if total == 0:
        print("protein_list.txt is empty (after removing comments/blanks). Nothing to do.")
        return

    # Decide group_size
    if args.group_size is not None:
        group_size = max(int(args.group_size), 1)
    else:
        target_chunks = int(args.target_chunks) if args.target_chunks else 10
        group_size = (total + target_chunks - 1) // target_chunks
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
        print(f"No proteins selected (start_idx={start_idx}, group_size={group_size}, total={total}). Nothing to do.")
        return

    print(f"Selected proteins: {start_idx}..{end_idx - 1} (count={len(proteins)}, total_list={total})")
    print(f"Bootstrap={'False' if args.no_bootstrap else 'True'} | p_threshold={args.p_threshold} | n_bootstraps={args.n_bootstraps}")
    print(f"✅ Writing outputs under: {results_dir}")

    batch_df = run_correlation_group(
        proteins=proteins,
        transcriptomics=transcriptomics,
        proteomics=proteomics,
        p_threshold=float(args.p_threshold),
        bootstrap=(not args.no_bootstrap),
        n_bootstraps=int(args.n_bootstraps),
        random_state=int(args.random_state) if args.random_state is not None else None,
        output_dir=str(results_dir),          # => new_corr/data/results
        per_protein_subdir="per_protein",     # => new_corr/data/results/per_protein
        per_protein_format=args.per_protein_format,
    )

    if not args.no_batch_file:
        batch_name = f"batch_{start_idx}_{end_idx - 1}_corr_summary"
        batch_path = results_dir / f"{batch_name}.parquet"
        try:
            batch_df.to_parquet(batch_path, compression="snappy", index=False)
            print("Saved batch summary to:", str(batch_path))
        except Exception:
            csv_path = results_dir / f"{batch_name}.csv"
            batch_df.to_csv(csv_path, index=False)
            print("Saved batch summary to:", str(csv_path))

    print("Done.")


if __name__ == "__main__":
    main()
