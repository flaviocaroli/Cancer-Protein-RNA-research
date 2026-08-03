#!/usr/bin/env python
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from utils import list_metric_shards, load_aligned_data, read_parquet_or_csv, write_parquet_or_csv


CONFIG_RE = re.compile(r"model=(?P<model>.+?)_norm=(?P<norm>.+?)_k=(?P<k>\d+)$")


def aggregate_metrics(metrics_dir: Path, out_path: Path) -> pd.DataFrame:
    files = list_metric_shards(metrics_dir)
    if not files:
        raise FileNotFoundError(f"No metric shards found in {metrics_dir}")
    frames = []
    for f in files:
        try:
            frames.append(read_parquet_or_csv(f))
        except Exception as exc:
            print(f"[WARN] could not read {f}: {exc}")
    if not frames:
        raise RuntimeError("No readable metric shards.")
    df = pd.concat(frames, axis=0, ignore_index=True)
    df = df.sort_values(["model", "normalization", "select_k", "split_id", "protein"])
    actual_path = write_parquet_or_csv(df, out_path)
    csv_path = out_path.with_suffix(".tsv.gz")
    df.to_csv(csv_path, sep="\t", index=False)
    print(f"[OK] wrote metrics_long: {actual_path} rows={len(df)}")
    print(f"[OK] wrote metrics_long TSV backup: {csv_path}")
    return df


def combine_prediction_config(config_dir: Path, out_dir: Path, sample_order, protein_order) -> None:
    match = CONFIG_RE.search(config_dir.name)
    if not match:
        print(f"[WARN] skipping unrecognized prediction config dir: {config_dir}")
        return
    model = match.group("model")
    norm = match.group("norm")
    k = int(match.group("k"))
    pred_files = sorted(config_dir.glob("*_pred.parquet")) + sorted(config_dir.glob("*_pred.csv.gz")) + sorted(config_dir.glob("*_pred.csv"))
    true_files = sorted(config_dir.glob("*_true.parquet")) + sorted(config_dir.glob("*_true.csv.gz")) + sorted(config_dir.glob("*_true.csv"))
    if not pred_files:
        print(f"[WARN] no prediction files for {config_dir}")
        return

    sample_to_i = {s: i for i, s in enumerate(sample_order)}
    protein_to_j = {p: j for j, p in enumerate(protein_order)}
    split_ids = set()
    for f in pred_files:
        df = read_parquet_or_csv(f)
        split_ids.update(df["split_id"].astype(int).unique().tolist())
    n_splits = max(split_ids) + 1 if split_ids else 0
    pred = np.full((n_splits, len(sample_order), len(protein_order)), np.nan, dtype=np.float32)
    ytrue = np.full_like(pred, np.nan)

    for pred_file in pred_files:
        true_file = Path(str(pred_file).replace("_pred", "_true"))
        if not true_file.exists():
            # try other suffixes if parquet fallback happened
            candidate = config_dir / pred_file.name.replace("_pred", "_true")
            true_file = candidate
        if not true_file.exists():
            print(f"[WARN] missing true file for {pred_file}; skipping")
            continue
        p_df = read_parquet_or_csv(pred_file)
        t_df = read_parquet_or_csv(true_file)
        value_cols = [c for c in p_df.columns if c not in {"split_id", "sample_id"}]
        for col in value_cols:
            if col not in protein_to_j:
                continue
            j = protein_to_j[col]
            p_vals = p_df[col].astype(float).values
            t_vals = t_df[col].astype(float).values if col in t_df.columns else np.full(len(p_df), np.nan)
            splits = p_df["split_id"].astype(int).values
            sample_ids = p_df["sample_id"].astype(str).values
            sample_indices = np.array([sample_to_i.get(s, -1) for s in sample_ids], dtype=int)
            keep = sample_indices >= 0
            pred[splits[keep], sample_indices[keep], j] = p_vals[keep]
            ytrue[splits[keep], sample_indices[keep], j] = t_vals[keep]

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{config_dir.name}.npz"
    np.savez_compressed(
        out_path,
        pred=pred,
        y_true=ytrue,
        sample_ids=np.array(sample_order, dtype=object),
        protein_ids=np.array(protein_order, dtype=object),
        model=model,
        normalization=norm,
        select_k=k,
    )
    print(f"[OK] wrote compressed prediction cube: {out_path} shape={pred.shape}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Aggregate metric shards and optionally combine prediction shards into compressed NPZ cubes.")
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--data_dir", default=None, help="Required for prediction cube sample/protein order.")
    parser.add_argument("--x_path", default=None)
    parser.add_argument("--y_path", default=None)
    parser.add_argument("--combine_predictions", type=int, default=0)
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    metrics_dir = results_dir / "metrics" / "metrics_shards"
    out_path = results_dir / "metrics" / "metrics_long.parquet"
    aggregate_metrics(metrics_dir, out_path)

    if args.combine_predictions:
        if not args.data_dir:
            raise ValueError("--data_dir is required when --combine_predictions 1")
        X, Y = load_aligned_data(args.data_dir, x_path=args.x_path, y_path=args.y_path)
        sample_order = X.index.astype(str).tolist()
        protein_order = Y.columns.astype(str).tolist()
        pred_root = results_dir / "prediction_shards"
        cube_dir = results_dir / "prediction_cubes"
        for config_dir in sorted([p for p in pred_root.iterdir() if p.is_dir()]):
            combine_prediction_config(config_dir, cube_dir, sample_order, protein_order)


if __name__ == "__main__":
    main()
