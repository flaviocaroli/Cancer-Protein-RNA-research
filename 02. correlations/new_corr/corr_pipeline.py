#!/usr/bin/env python3
"""Run one chunk of the correlation pipeline (or assemble the full matrix).

Defaults reproduce the original run (Pearson, protein vs RNA, results in data/results/).
  --method   pearson | spearman
  --target   rna     (features = transcriptomics_filtered genes)  -> full_pearson100 / full_spearman100
             protein (features = proteomics_matched proteins)      -> full_prot
  --assemble          skip computation; stack per_protein/* into <results_dir>/<name>.joblib
"""
import argparse
import os
import sys
from pathlib import Path

from joblib import dump, load
from corr_functions import assemble_full_matrix, run_correlation_group

DEFAULT_SUBDIR = {("pearson", "rna"): "results",
                  ("spearman", "rna"): "results_spearman",
                  ("pearson", "protein"): "results_prot",
                  ("spearman", "protein"): "results_prot_spearman"}


def full_name(method, target, n_boot, bootstrap):
    if target == "protein":
        return "full_prot" if method == "pearson" else "full_prot_spearman"
    return f"full_{method}{n_boot if bootstrap else 'param'}"


def read_protein_list(protein_list_path: str):
    with open(protein_list_path, "r") as f:
        return [l.strip() for l in f if l.strip() and not l.strip().startswith("#")]


def main():
    p = argparse.ArgumentParser(description="Chunked protein correlation pipeline.")
    p.add_argument("--chunk_id", type=int, default=None)
    p.add_argument("--start_idx", type=int, default=None)
    p.add_argument("--group_size", type=int, default=None)
    p.add_argument("--target_chunks", type=int, default=10)
    p.add_argument("--method", choices=["pearson", "spearman"], default="pearson")
    p.add_argument("--target", choices=["rna", "protein"], default="rna")
    p.add_argument("--p_threshold", type=float, default=0.05)
    p.add_argument("--no_bootstrap", action="store_true")
    p.add_argument("--n_bootstraps", type=int, default=100)
    p.add_argument("--random_state", type=int, default=42)
    p.add_argument("--data_dir", type=str, default=None, help="Default: <script_dir>/data")
    p.add_argument("--results_dir", type=str, default=None, help="Default: <data_dir>/<results|results_spearman|results_prot>")
    p.add_argument("--protein_list", type=str, default=None)
    p.add_argument("--per_protein_format", choices=["parquet", "csv"], default="parquet")
    p.add_argument("--no_batch_file", action="store_true")
    p.add_argument("--assemble", action="store_true", help="Assemble per_protein/* into one joblib and exit.")
    p.add_argument("--assemble_out", type=str, default=None, help="Output .joblib (default <results_dir>/<full_name>.joblib)")
    args = p.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = Path(args.data_dir).resolve() if args.data_dir else (here / "data").resolve()
    results_dir = (Path(args.results_dir).resolve() if args.results_dir
                   else (data_dir / DEFAULT_SUBDIR[(args.method, args.target)]).resolve())
    results_dir.mkdir(parents=True, exist_ok=True)
    protein_list_path = args.protein_list or str(data_dir / "protein_list.txt")
    if not os.path.exists(protein_list_path):
        sys.exit(f"Error: missing protein list at {protein_list_path}")
    all_proteins = read_protein_list(protein_list_path)

    if args.assemble:
        name = full_name(args.method, args.target, args.n_bootstraps, not args.no_bootstrap)
        out = Path(args.assemble_out) if args.assemble_out else results_dir / f"{name}.joblib"
        full = assemble_full_matrix(results_dir / "per_protein", protein_order=all_proteins,
                                    square=(args.target == "protein"))
        missing = [x for x in all_proteins if x not in full.index]
        dump(full, out, compress=3)
        print(f"Assembled {full.shape[0]} proteins x {full.shape[1]} features -> {out}")
        if missing:
            print(f"WARNING: {len(missing)} proteins from the list have no per-protein file "
                  f"(first: {missing[:5]}). Re-run their chunks before using the matrix.")
        return

    transcript_path = data_dir / "transcriptomics_filtered.joblib"
    protein_path = data_dir / "proteomics_matched.joblib"
    for path in (transcript_path, protein_path):
        if not path.exists():
            sys.exit(f"Error: missing {path}")
    print("Loading data...")
    proteomics = load(protein_path)
    features = proteomics if args.target == "protein" else load(transcript_path)
    features.index.name = "samples"
    proteomics.index.name = "samples"

    total = len(all_proteins)
    if total == 0:
        print("protein_list.txt is empty. Nothing to do.")
        return
    group_size = max(int(args.group_size), 1) if args.group_size else max((total + args.target_chunks - 1) // args.target_chunks, 1)
    start_idx = args.start_idx if args.start_idx is not None else (args.chunk_id * group_size if args.chunk_id is not None else 0)
    end_idx = min(start_idx + group_size, total)
    proteins = all_proteins[start_idx:end_idx]
    if not proteins:
        print(f"No proteins selected (start_idx={start_idx}, total={total}).")
        return
    print(f"Proteins {start_idx}..{end_idx - 1} | method={args.method} target={args.target} "
          f"bootstrap={not args.no_bootstrap} n_boot={args.n_bootstraps} -> {results_dir}")

    batch_df = run_correlation_group(
        proteins, features, proteomics, p_threshold=args.p_threshold,
        bootstrap=not args.no_bootstrap, n_bootstraps=args.n_bootstraps,
        random_state=args.random_state, output_dir=str(results_dir),
        per_protein_subdir="per_protein", per_protein_format=args.per_protein_format,
        method=args.method)

    if not args.no_batch_file:
        stem = results_dir / f"batch_{start_idx}_{end_idx - 1}_corr_summary"
        try:
            batch_df.to_parquet(stem.with_suffix(".parquet"), compression="snappy", index=False)
        except Exception:
            batch_df.to_csv(stem.with_suffix(".csv"), index=False)
    print("Done.")


if __name__ == "__main__":
    main()
