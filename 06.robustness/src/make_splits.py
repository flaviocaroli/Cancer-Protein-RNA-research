#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from utils import load_aligned_data, make_repeated_splits, resolve_input_paths, load_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Create repeated 80/20 CCLE split manifest.")
    parser.add_argument("--data_dir", required=True)
    parser.add_argument("--x_path", default=None)
    parser.add_argument("--y_path", default=None)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--n_repeats", type=int, default=50)
    parser.add_argument("--test_size", type=float, default=0.20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--sample_info_path", default=None)
    parser.add_argument("--stratify_column", default=None, help="Optional sample_info column to stratify random splits.")
    args = parser.parse_args()

    X, Y = load_aligned_data(args.data_dir, x_path=args.x_path, y_path=args.y_path)
    stratify = None
    if args.stratify_column:
        _, _, sample_info_resolved = resolve_input_paths(args.data_dir, args.x_path, args.y_path, args.sample_info_path)
        if sample_info_resolved is None:
            raise FileNotFoundError("--stratify_column was supplied but no sample_info file was found.")
        sample_info = load_table(sample_info_resolved)
        sample_info.index = sample_info.index.astype(str)
        if args.stratify_column not in sample_info.columns:
            raise ValueError(f"Column {args.stratify_column} not in sample_info. Available: {list(sample_info.columns)[:20]}")
        stratify = sample_info.reindex(X.index)[args.stratify_column].astype(str).fillna("missing").values
        counts = pd.Series(stratify).value_counts()
        if (counts < 2).any():
            bad = counts[counts < 2].index.tolist()
            raise ValueError(f"Cannot stratify: classes with <2 samples: {bad[:10]}")

    manifest = make_repeated_splits(
        sample_ids=X.index.tolist(),
        n_repeats=args.n_repeats,
        test_size=args.test_size,
        seed=args.seed,
        stratify=stratify,
    )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "repeated_8020_splits.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    summary = manifest.groupby(["split_id", "set"]).size().unstack(fill_value=0).reset_index()
    summary.to_csv(out_dir / "split_summary.tsv", sep="\t", index=False)
    print(f"[OK] wrote {manifest_path}")
    print(summary.head())


if __name__ == "__main__":
    main()
