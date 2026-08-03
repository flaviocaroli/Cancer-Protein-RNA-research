#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import sys

import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from utils import read_table, write_table, setup_paper_style, save_figure


def collect(pattern: str, shard_dir: Path) -> pd.DataFrame:
    files = sorted(shard_dir.glob(pattern))
    frames = []
    for path in files:
        try:
            frames.append(read_table(path))
        except Exception as exc:
            print(f"[WARN] Could not read {path}: {exc}")
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def make_summary_plots(metrics: pd.DataFrame, out_dir: Path) -> None:
    if metrics.empty:
        return
    import matplotlib.pyplot as plt
    import seaborn as sns
    setup_paper_style()
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    sns.boxplot(data=metrics, x="test_dataset", y="pearson", hue="model", showfliers=False, ax=ax)
    sns.stripplot(data=metrics.sample(min(len(metrics), 4000), random_state=42), x="test_dataset", y="pearson", hue="model", dodge=True, alpha=0.12, size=1.5, legend=False, ax=ax)
    ax.axhline(0.5, color="0.35", linewidth=1, linestyle="--")
    ax.set_xlabel("")
    ax.set_ylabel("Protein-wise Pearson correlation")
    ax.set_title("Held-out platform transfer performance")
    ax.legend(frameon=False, title="Model")
    save_figure(fig, fig_dir / "transfer_pearson_by_test_dataset")
    plt.close(fig)

    summary = metrics.groupby(["model", "test_dataset"], as_index=False).agg(
        proteins=("protein", "nunique"),
        median_pearson=("pearson", "median"),
        mean_pearson=("pearson", "mean"),
        fraction_ge_0_3=("pearson", lambda x: (x >= 0.3).mean()),
        fraction_ge_0_5=("pearson", lambda x: (x >= 0.5).mean()),
        fraction_ge_0_7=("pearson", lambda x: (x >= 0.7).mean()),
        median_rmse_eval_scaled=("rmse_eval_scaled", "median"),
        median_calibration_slope=("calibration_slope", "median"),
    )
    summary.to_csv(out_dir / "summary_by_model_test.tsv", sep="\t", index=False)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--results-dir", required=True)
    p.add_argument("--model", choices=["enet", "xgb", "all"], default="all")
    args = p.parse_args()
    results_dir = Path(args.results_dir)
    shard_dir = results_dir / "shards"
    models = ["enet", "xgb"] if args.model == "all" else [args.model]

    metrics_frames, pred_frames, feat_frames = [], [], []
    for model in models:
        metrics_frames.append(collect(f"metrics_{model}_task_*.parquet", shard_dir))
        metrics_frames.append(collect(f"metrics_{model}_task_*.csv.gz", shard_dir))
        pred_frames.append(collect(f"predictions_{model}_task_*.parquet", shard_dir))
        pred_frames.append(collect(f"predictions_{model}_task_*.csv.gz", shard_dir))
        feat_frames.append(collect(f"features_{model}_task_*.parquet", shard_dir))
        feat_frames.append(collect(f"features_{model}_task_*.csv.gz", shard_dir))

    metrics = pd.concat([x for x in metrics_frames if not x.empty], ignore_index=True) if any(not x.empty for x in metrics_frames) else pd.DataFrame()
    preds = pd.concat([x for x in pred_frames if not x.empty], ignore_index=True) if any(not x.empty for x in pred_frames) else pd.DataFrame()
    features = pd.concat([x for x in feat_frames if not x.empty], ignore_index=True) if any(not x.empty for x in feat_frames) else pd.DataFrame()

    if not metrics.empty:
        metrics = metrics.drop_duplicates(["scenario", "model", "protein", "test_dataset"], keep="last")
        write_table(metrics, results_dir / "metrics_all.parquet")
    if not preds.empty:
        preds = preds.drop_duplicates(["scenario", "model", "test_dataset", "sample_id", "protein"], keep="last")
        write_table(preds, results_dir / "predictions_all.parquet")
    if not features.empty:
        key = [c for c in ["scenario", "protein", "gene", "coefficient", "importance"] if c in features.columns]
        if key:
            features = features.drop_duplicates(key, keep="last")
        write_table(features, results_dir / "selected_features_all.parquet")

    make_summary_plots(metrics, results_dir)
    print(f"[DONE] Aggregate outputs in {results_dir}")

if __name__ == "__main__":
    main()
