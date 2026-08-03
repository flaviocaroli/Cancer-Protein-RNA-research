#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from utils import read_table, write_table, setup_paper_style, save_figure


def find_metrics(path: Path) -> pd.DataFrame:
    for candidate in [path / "metrics_all.parquet", path / "metrics_all.csv.gz"]:
        if candidate.exists():
            return read_table(candidate)
    raise FileNotFoundError(f"No aggregated metrics in {path}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--project-root", default=str(Path(__file__).resolve().parents[1]))
    p.add_argument("--scenario", required=True)
    p.add_argument("--y-norm", default="domain_zscore")
    p.add_argument("--model", choices=["enet", "xgb"], required=True)
    p.add_argument("--normalizations", default="zscore,rank,binning")
    args = p.parse_args()

    root = Path(args.project_root)
    norms = [x.strip() for x in args.normalizations.split(",") if x.strip()]
    frames = []
    for norm in norms:
        path = root / "results" / args.scenario / f"x-{norm}__y-{args.y_norm}" / args.model
        try:
            df = find_metrics(path)
        except FileNotFoundError as exc:
            print(f"[WARN] {exc}")
            continue
        df = df.copy()
        df["normalization"] = norm
        frames.append(df)
    if not frames:
        raise SystemExit("No normalization results found")
    all_df = pd.concat(frames, ignore_index=True)
    out_dir = root / "results" / "normalization_comparisons" / args.scenario / args.model
    out_dir.mkdir(parents=True, exist_ok=True)
    write_table(all_df, out_dir / "protein_metrics_all_normalizations.parquet")

    rows = []
    for (norm, test), sub in all_df.groupby(["normalization", "test_dataset"]):
        rows.append({
            "normalization": norm,
            "test_dataset": test,
            "n_proteins": sub["protein"].nunique(),
            "median_pearson": sub["pearson"].median(),
            "mean_pearson": sub["pearson"].mean(),
            "fraction_pearson_ge_0_3": (sub["pearson"] >= 0.3).mean(),
            "fraction_pearson_ge_0_5": (sub["pearson"] >= 0.5).mean(),
            "fraction_pearson_ge_0_7": (sub["pearson"] >= 0.7).mean(),
            "median_rmse_eval_scaled": sub["rmse_eval_scaled"].median(),
            "median_calibration_slope": sub["calibration_slope"].median(),
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(out_dir / "normalization_summary.tsv", sep="\t", index=False)

    import matplotlib.pyplot as plt
    import seaborn as sns
    setup_paper_style()

    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    sns.boxplot(data=all_df, x="test_dataset", y="pearson", hue="normalization", showfliers=False, ax=ax)
    ax.axhline(0.5, color="0.35", linestyle="--", linewidth=1)
    ax.set_xlabel("")
    ax.set_ylabel("Protein-wise Pearson correlation")
    ax.set_title(f"{args.scenario}: RNA normalization comparison")
    ax.legend(frameon=False, title="RNA normalization")
    save_figure(fig, out_dir / "normalization_pearson_distributions")
    plt.close(fig)

    frac = summary.melt(
        id_vars=["normalization", "test_dataset"],
        value_vars=["fraction_pearson_ge_0_3", "fraction_pearson_ge_0_5", "fraction_pearson_ge_0_7"],
        var_name="threshold", value_name="fraction",
    )
    frac["threshold"] = frac["threshold"].str.extract(r"([0-9]_[0-9])")[0].str.replace("_", ".", regex=False)
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    sns.barplot(data=frac, x="test_dataset", y="fraction", hue="normalization", errorbar=None, ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("Fraction of proteins above Pearson threshold")
    ax.set_title("Transferable protein fraction by normalization")
    ax.legend(frameon=False, title="RNA normalization")
    save_figure(fig, out_dir / "normalization_fraction_predictable")
    plt.close(fig)

    print(f"[DONE] {out_dir}")

if __name__ == "__main__":
    main()
