#!/usr/bin/env python
from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd

from utils import chunk_ranges, load_aligned_data


def parse_csv(value: str):
    return [v.strip() for v in str(value).split(",") if v.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a SLURM-friendly task table for repeated 80/20 runs.")
    parser.add_argument("--out_path", required=True)
    parser.add_argument("--data_dir", default=None, help="Optional; if supplied, protein count/order is read from aligned proteomics.")
    parser.add_argument("--x_path", default=None)
    parser.add_argument("--y_path", default=None)
    parser.add_argument("--models", default="proxy,enet")
    parser.add_argument("--normalizations", default="zscore")
    parser.add_argument("--select_ks", default="2000")
    parser.add_argument("--n_repeats", type=int, default=20)
    parser.add_argument("--proteins_total", type=int, default=None)
    parser.add_argument("--protein_list", default=None, help="Optional newline-separated protein list; overrides proteins_total.")
    parser.add_argument("--proteins_per_task", type=int, default=50)
    parser.add_argument("--protein_chunks", type=int, default=20, help="Used by --task_mode packed_by_protein. Creates this many protein chunks/tasks.")
    parser.add_argument("--splits_per_task", type=int, default=5)
    parser.add_argument("--cv_folds", type=int, default=5)
    parser.add_argument("--grid_preset", default="fast", choices=["tiny", "fast", "wide"])
    parser.add_argument("--max_param_combos", type=int, default=8)
    parser.add_argument("--keep_predictions", type=int, default=1)
    parser.add_argument(
        "--task_mode",
        default="packed_by_protein",
        choices=["packed_by_protein", "grid"],
        help=(
            "packed_by_protein creates one SLURM task per protein chunk; each task loops over all requested models, "
            "normalizations, k values and splits. This is best for clusters with a low job limit. "
            "grid creates one task per model/normalization/k/split/protein chunk."
        ),
    )
    args = parser.parse_args()

    if args.protein_list:
        proteins = [line.strip() for line in Path(args.protein_list).read_text().splitlines() if line.strip()]
        proteins_total = len(proteins)
    elif args.data_dir:
        _, Y = load_aligned_data(args.data_dir, x_path=args.x_path, y_path=args.y_path)
        proteins_total = Y.shape[1]
    elif args.proteins_total is not None:
        proteins_total = args.proteins_total
    else:
        raise ValueError("Supply --data_dir, --protein_list, or --proteins_total.")

    models = parse_csv(args.models)
    normalizations = parse_csv(args.normalizations)
    select_ks = [int(x) for x in parse_csv(args.select_ks)]

    rows = []
    task_id = 1

    if args.task_mode == "packed_by_protein":
        if args.protein_chunks < 1:
            raise ValueError("--protein_chunks must be >= 1")
        proteins_per_chunk = int(math.ceil(proteins_total / args.protein_chunks))
        protein_chunks = chunk_ranges(proteins_total, proteins_per_chunk)
        for protein_start, protein_end in protein_chunks:
            rows.append(
                {
                    "task_id": task_id,
                    "task_mode": "packed_by_protein",
                    "model": "packed",
                    "normalization": "packed",
                    "select_k": -1,
                    "split_start": 0,
                    "split_end": args.n_repeats,
                    "protein_start": protein_start,
                    "protein_end": protein_end,
                    "cv_folds": args.cv_folds,
                    "grid_preset": args.grid_preset,
                    "max_param_combos": args.max_param_combos,
                    "keep_predictions": args.keep_predictions,
                    "packed_models": ",".join(models),
                    "packed_normalizations": ",".join(normalizations),
                    "packed_select_ks": ",".join(str(k) for k in select_ks),
                }
            )
            task_id += 1
    else:
        split_chunks = chunk_ranges(args.n_repeats, args.splits_per_task)
        protein_chunks = chunk_ranges(proteins_total, args.proteins_per_task)
        for model in models:
            model = model.lower()
            if model == "proxy":
                k_values = [0]
            else:
                k_values = select_ks
            for normalization in normalizations:
                for k in k_values:
                    for split_start, split_end in split_chunks:
                        for protein_start, protein_end in protein_chunks:
                            rows.append(
                                {
                                    "task_id": task_id,
                                    "task_mode": "grid",
                                    "model": model,
                                    "normalization": normalization,
                                    "select_k": int(k),
                                    "split_start": split_start,
                                    "split_end": split_end,
                                    "protein_start": protein_start,
                                    "protein_end": protein_end,
                                    "cv_folds": args.cv_folds,
                                    "grid_preset": args.grid_preset,
                                    "max_param_combos": args.max_param_combos,
                                    "keep_predictions": args.keep_predictions,
                                    "packed_models": "",
                                    "packed_normalizations": "",
                                    "packed_select_ks": "",
                                }
                            )
                            task_id += 1

    out = pd.DataFrame(rows)
    out_path = Path(args.out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    print(f"[OK] wrote {len(out)} tasks to {out_path}")
    print(out.head())


if __name__ == "__main__":
    main()
