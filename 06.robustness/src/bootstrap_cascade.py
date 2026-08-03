#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

from utils import bootstrap_ci


def parse_config(text: str) -> Tuple[str, str, int]:
    parts = text.split(":")
    if len(parts) != 3:
        raise ValueError(f"Config must look like model:normalization:k, got {text}")
    return parts[0].lower(), parts[1].lower(), int(parts[2])


def select_config(df: pd.DataFrame, config: Tuple[str, str, int]) -> pd.DataFrame:
    model, norm, k = config
    sub = df[
        (df["model"].astype(str).str.lower() == model)
        & (df["normalization"].astype(str).str.lower() == norm)
        & (df["select_k"].astype(int) == int(k))
    ].copy()
    if sub.empty:
        raise ValueError(f"No rows found for config {model}:{norm}:{k}")
    return sub


def summarize_model(sub: pd.DataFrame, metric: str, threshold: float, n_boot: int, seed: int, model_label: str) -> pd.DataFrame:
    rows = []
    for idx, (protein, g) in enumerate(sub.groupby("protein", sort=False)):
        values = pd.to_numeric(g[metric], errors="coerce").values.astype(float)
        values = values[np.isfinite(values)]
        center, low, high = bootstrap_ci(values, n_boot=n_boot, seed=seed + idx, statistic="median")
        rows.append(
            {
                "protein": protein,
                "model_label": model_label,
                "n_splits_valid": int(len(values)),
                "median_metric": center,
                "ci_low": low,
                "ci_high": high,
                "q25": float(np.nanpercentile(values, 25)) if len(values) else np.nan,
                "q75": float(np.nanpercentile(values, 75)) if len(values) else np.nan,
                "pass_rate": float(np.nanmean(values >= threshold)) if len(values) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def paired_gain(a: pd.DataFrame, b: pd.DataFrame, metric: str, n_boot: int, seed: int, label: str) -> pd.DataFrame:
    # b - a on matched split/protein rows
    keep_cols = ["split_id", "protein", metric]
    ma = a[keep_cols].rename(columns={metric: "metric_a"})
    mb = b[keep_cols].rename(columns={metric: "metric_b"})
    merged = ma.merge(mb, on=["split_id", "protein"], how="inner")
    merged["delta"] = pd.to_numeric(merged["metric_b"], errors="coerce") - pd.to_numeric(merged["metric_a"], errors="coerce")
    rows = []
    for idx, (protein, g) in enumerate(merged.groupby("protein", sort=False)):
        values = g["delta"].values.astype(float)
        values = values[np.isfinite(values)]
        center, low, high = bootstrap_ci(values, n_boot=n_boot, seed=seed + idx, statistic="median")
        rows.append(
            {
                "protein": protein,
                "comparison": label,
                "n_splits_valid": int(len(values)),
                "median_delta": center,
                "ci_low": low,
                "ci_high": high,
                "fraction_delta_gt_0": float(np.nanmean(values > 0)) if len(values) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build robust cascade atlas from repeated split metrics.")
    parser.add_argument("--metrics_long", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--proxy_config", default="proxy:zscore:0")
    parser.add_argument("--enet_config", default="enet:zscore:2000")
    parser.add_argument("--xgb_config", default=None, help="Optional, e.g. xgb:zscore:2000")
    parser.add_argument("--metric", default="pearson")
    parser.add_argument("--threshold", type=float, default=0.50)
    parser.add_argument("--pass_rate", type=float, default=0.80)
    parser.add_argument("--gain_threshold", type=float, default=0.05)
    parser.add_argument("--n_boot", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--require_ci_low", type=int, default=0, help="If 1, model median CI low must exceed threshold.")
    args = parser.parse_args()

    metrics_path = Path(args.metrics_long)
    if metrics_path.suffix == ".parquet":
        df = pd.read_parquet(metrics_path)
    else:
        sep = "\t" if metrics_path.name.endswith(".tsv") or metrics_path.name.endswith(".tsv.gz") else ","
        df = pd.read_csv(metrics_path, sep=sep)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    proxy_cfg = parse_config(args.proxy_config)
    enet_cfg = parse_config(args.enet_config)
    xgb_cfg = parse_config(args.xgb_config) if args.xgb_config else None

    proxy = select_config(df, proxy_cfg)
    enet = select_config(df, enet_cfg)
    xgb = select_config(df, xgb_cfg) if xgb_cfg else None

    summaries = [summarize_model(proxy, args.metric, args.threshold, args.n_boot, args.seed, "proxy")]
    summaries.append(summarize_model(enet, args.metric, args.threshold, args.n_boot, args.seed + 100000, "enet"))
    if xgb is not None:
        summaries.append(summarize_model(xgb, args.metric, args.threshold, args.n_boot, args.seed + 200000, "xgb"))
    summary_long = pd.concat(summaries, axis=0, ignore_index=True)
    summary_long.to_csv(out_dir / "bootstrap_metric_ci.tsv", sep="\t", index=False)

    gain_frames = [paired_gain(proxy, enet, args.metric, args.n_boot, args.seed + 300000, "enet_minus_proxy")]
    if xgb is not None:
        gain_frames.append(paired_gain(enet, xgb, args.metric, args.n_boot, args.seed + 400000, "xgb_minus_enet"))
    gains = pd.concat(gain_frames, axis=0, ignore_index=True)
    gains.to_csv(out_dir / "paired_gain_tests.tsv", sep="\t", index=False)

    # Wide protein-level table
    wide = None
    for label in ["proxy", "enet", "xgb"]:
        sub = summary_long[summary_long["model_label"] == label].copy()
        if sub.empty:
            continue
        sub = sub.rename(
            columns={
                "n_splits_valid": f"{label}_n_splits_valid",
                "median_metric": f"{label}_median_{args.metric}",
                "ci_low": f"{label}_ci_low",
                "ci_high": f"{label}_ci_high",
                "q25": f"{label}_q25",
                "q75": f"{label}_q75",
                "pass_rate": f"{label}_pass_rate",
            }
        )
        sub = sub.drop(columns=["model_label"])
        wide = sub if wide is None else wide.merge(sub, on="protein", how="outer")
    assert wide is not None

    for comp in ["enet_minus_proxy", "xgb_minus_enet"]:
        sub = gains[gains["comparison"] == comp].copy()
        if sub.empty:
            continue
        prefix = comp.replace("_minus_", "_vs_")
        sub = sub.rename(
            columns={
                "n_splits_valid": f"{prefix}_n_splits_valid",
                "median_delta": f"{prefix}_median_delta",
                "ci_low": f"{prefix}_ci_low",
                "ci_high": f"{prefix}_ci_high",
                "fraction_delta_gt_0": f"{prefix}_fraction_delta_gt_0",
            }
        ).drop(columns=["comparison"])
        wide = wide.merge(sub, on="protein", how="left")

    def model_good(row, label: str) -> bool:
        pr = row.get(f"{label}_pass_rate", np.nan)
        low = row.get(f"{label}_ci_low", np.nan)
        if not np.isfinite(pr) or pr < args.pass_rate:
            return False
        if args.require_ci_low and (not np.isfinite(low) or low < args.threshold):
            return False
        return True

    classes = []
    for _, row in wide.iterrows():
        proxy_good = model_good(row, "proxy")
        enet_good = model_good(row, "enet")
        xgb_good = model_good(row, "xgb") if xgb is not None else False
        enet_gain = row.get("enet_vs_proxy_ci_low", np.nan)
        xgb_gain = row.get("xgb_vs_enet_ci_low", np.nan)

        if proxy_good:
            cls = "robust_proxy_sufficient"
        elif enet_good and np.isfinite(enet_gain) and enet_gain > args.gain_threshold:
            cls = "robust_linear_improvable"
        elif xgb_good and np.isfinite(xgb_gain) and xgb_gain > args.gain_threshold:
            cls = "robust_nonlinear_improvable"
        else:
            pass_rates = [row.get("proxy_pass_rate", np.nan), row.get("enet_pass_rate", np.nan)]
            if xgb is not None:
                pass_rates.append(row.get("xgb_pass_rate", np.nan))
            finite = [p for p in pass_rates if np.isfinite(p)]
            if finite and max(finite) > 0.0:
                cls = "unstable_or_borderline"
            else:
                cls = "unresolved"
        classes.append(cls)
    wide["robust_class"] = classes
    wide["threshold"] = args.threshold
    wide["pass_rate_rule"] = args.pass_rate
    wide["gain_threshold_rule"] = args.gain_threshold

    out_path = out_dir / "robust_cascade_atlas.tsv"
    wide.to_csv(out_path, sep="\t", index=False)
    counts = wide["robust_class"].value_counts().rename_axis("robust_class").reset_index(name="n_proteins")
    counts.to_csv(out_dir / "robust_cascade_counts.tsv", sep="\t", index=False)
    print(f"[OK] wrote atlas: {out_path}")
    print(counts)


if __name__ == "__main__":
    main()
