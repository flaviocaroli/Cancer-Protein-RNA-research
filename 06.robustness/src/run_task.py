#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from utils import (
    apply_preprocessor,
    build_gene_lookup,
    compute_metrics,
    config_name,
    ensure_dir,
    fit_preprocessor,
    get_split_samples,
    linear_proxy_predict,
    load_aligned_data,
    same_gene_column,
    tune_and_predict_enet,
    tune_and_predict_xgb,
    valid_xy_for_protein,
    write_parquet_or_csv,
)


def parse_csv(value: object):
    if pd.isna(value):
        return []
    return [v.strip() for v in str(value).split(",") if v.strip()]


def parse_int_csv(value: object):
    return [int(v) for v in parse_csv(value)]


def run_proxy_for_protein(X, y, train_samples, test_samples, normalization, gene_lookup, seed):
    X_train_raw, y_train, kept_train = valid_xy_for_protein(X, y, train_samples)
    X_test_raw, y_test, kept_test = valid_xy_for_protein(X, y, test_samples)
    info = {"n_train": len(y_train), "n_test": len(y_test), "test_samples": kept_test}
    gene_col = same_gene_column(y.name, gene_lookup)
    if gene_col is None:
        info.update({"status": "missing_same_gene", "gene_col": None})
        return y_test, np.full(len(y_test), np.nan), info
    if len(y_train) < 3 or len(y_test) < 3 or np.std(y_train) == 0:
        info.update({"status": "too_few_samples_or_constant_y", "gene_col": gene_col})
        return y_test, np.full(len(y_test), np.nan), info
    prep = fit_preprocessor(X_train_raw[[gene_col]], normalization)
    X_train = apply_preprocessor(X_train_raw[[gene_col]], prep)
    X_test = apply_preprocessor(X_test_raw[[gene_col]], prep)
    try:
        y_pred = linear_proxy_predict(X_train, X_test, y_train, gene_col)
        info.update({"status": "ok", "gene_col": gene_col})
    except Exception as exc:
        y_pred = np.full(len(y_test), np.nan)
        info.update({"status": f"proxy_failed: {exc}", "gene_col": gene_col})
    return y_test, y_pred, info


def run_one_config(
    *,
    task_id: int,
    X: pd.DataFrame,
    Y: pd.DataFrame,
    split_manifest: pd.DataFrame,
    gene_lookup: dict,
    model: str,
    normalization: str,
    select_k: int,
    split_start: int,
    split_end: int,
    protein_start: int,
    protein_end: int,
    cv_folds: int,
    grid_preset: str,
    max_param_combos: int,
    keep_predictions: int,
    out_dir: Path,
    seed_base: int,
    n_jobs: int,
    skip_existing: int,
) -> None:
    model = str(model).lower()
    normalization = str(normalization).lower()
    select_k = int(select_k)

    metric_dir = ensure_dir(out_dir / "metrics" / "metrics_shards")
    pred_dir = ensure_dir(out_dir / "prediction_shards" / config_name(model, normalization, select_k))

    cfg = config_name(model, normalization, select_k)
    metric_path = metric_dir / f"task_{task_id:06d}_{cfg}.parquet"
    pred_path = pred_dir / f"task_{task_id:06d}_pred.parquet"
    true_path = pred_dir / f"task_{task_id:06d}_true.parquet"

    if skip_existing and metric_path.exists() and (not keep_predictions or (pred_path.exists() and true_path.exists())):
        print(f"[SKIP] existing task {task_id} config {cfg}: {metric_path}")
        return

    proteins = Y.columns.tolist()[protein_start:protein_end]
    metrics_rows = []
    pred_frames = []
    true_frames = []

    for split_id in range(int(split_start), int(split_end)):
        train_samples, test_samples = get_split_samples(split_manifest, split_id)
        pred_df = pd.DataFrame({"split_id": int(split_id), "sample_id": test_samples}).set_index("sample_id")
        true_df = pd.DataFrame({"split_id": int(split_id), "sample_id": test_samples}).set_index("sample_id")

        for local_i, protein in enumerate(proteins):
            y = Y[protein]
            seed = seed_base + int(split_id) * 100000 + int(protein_start) + local_i
            if model == "proxy":
                y_test, y_pred, info = run_proxy_for_protein(
                    X, y, train_samples, test_samples, normalization, gene_lookup, seed
                )
                best_params_json = "{}"
                best_cv_mse = np.nan
                n_selected = 1 if info.get("gene_col") is not None else 0
                n_nonzero = np.nan
            elif model == "enet":
                y_test, y_pred, info = tune_and_predict_enet(
                    X=X,
                    y=y,
                    train_samples=train_samples,
                    test_samples=test_samples,
                    normalization=normalization,
                    select_k=select_k,
                    cv_folds=cv_folds,
                    grid_preset=grid_preset,
                    max_param_combos=max_param_combos,
                    seed=seed,
                )
                best_params_json = info.get("best_params_json", "{}")
                best_cv_mse = info.get("best_cv_mse", np.nan)
                n_selected = info.get("n_selected_features", np.nan)
                n_nonzero = info.get("n_nonzero_features", np.nan)
            elif model == "xgb":
                y_test, y_pred, info = tune_and_predict_xgb(
                    X=X,
                    y=y,
                    train_samples=train_samples,
                    test_samples=test_samples,
                    normalization=normalization,
                    select_k=select_k,
                    cv_folds=cv_folds,
                    grid_preset=grid_preset,
                    max_param_combos=max_param_combos,
                    seed=seed,
                    n_jobs=n_jobs,
                )
                best_params_json = info.get("best_params_json", "{}")
                best_cv_mse = info.get("best_cv_mse", np.nan)
                n_selected = info.get("n_selected_features", np.nan)
                n_nonzero = info.get("n_nonzero_features", np.nan)
            else:
                raise ValueError(f"Unknown model: {model}")

            met = compute_metrics(y_test, y_pred)
            met.update(
                {
                    "split_id": int(split_id),
                    "model": model,
                    "normalization": normalization,
                    "select_k": int(select_k),
                    "protein": protein,
                    "n_train": int(info.get("n_train", np.nan)) if not pd.isna(info.get("n_train", np.nan)) else np.nan,
                    "status": str(info.get("status", "unknown")),
                    "best_params_json": best_params_json,
                    "best_cv_mse": best_cv_mse,
                    "n_selected_features": n_selected,
                    "n_nonzero_features": n_nonzero,
                    "protein_start": protein_start,
                    "protein_end": protein_end,
                    "task_id": task_id,
                }
            )
            if model == "proxy":
                met["same_gene_col"] = info.get("gene_col", None)
            metrics_rows.append(met)

            if keep_predictions:
                kept_test = info.get("test_samples", [])
                if len(kept_test) == len(y_test):
                    true_df.loc[kept_test, protein] = y_test.astype(np.float32)
                    pred_df.loc[kept_test, protein] = np.asarray(y_pred, dtype=np.float32)

        if keep_predictions:
            pred_frames.append(pred_df.reset_index())
            true_frames.append(true_df.reset_index())

    metrics_df = pd.DataFrame(metrics_rows)
    metric_written = write_parquet_or_csv(metrics_df, metric_path)
    print(f"[OK] wrote metrics: {metric_written} rows={len(metrics_df)}")

    if keep_predictions:
        pred_all = pd.concat(pred_frames, axis=0, ignore_index=True) if pred_frames else pd.DataFrame()
        true_all = pd.concat(true_frames, axis=0, ignore_index=True) if true_frames else pd.DataFrame()
        pred_written = write_parquet_or_csv(pred_all, pred_path)
        true_written = write_parquet_or_csv(true_all, true_path)
        print(f"[OK] wrote predictions: {pred_written} shape={pred_all.shape}")
        print(f"[OK] wrote y_true: {true_written} shape={true_all.shape}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run one repeated-8020 robustness task.")
    parser.add_argument("--task_table", required=True)
    parser.add_argument("--task_id", type=int, required=True)
    parser.add_argument("--data_dir", required=True)
    parser.add_argument("--x_path", default=None)
    parser.add_argument("--y_path", default=None)
    parser.add_argument("--split_manifest", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n_jobs", type=int, default=1)
    parser.add_argument("--skip_existing", type=int, default=1)
    args = parser.parse_args()

    task_table = pd.read_csv(args.task_table, sep="\t")
    row = task_table.loc[task_table["task_id"].astype(int) == int(args.task_id)]
    if row.empty:
        raise ValueError(f"task_id {args.task_id} not found in {args.task_table}")
    task = row.iloc[0].to_dict()

    task_mode = str(task.get("task_mode", "grid"))
    split_start = int(task["split_start"])
    split_end = int(task["split_end"])
    protein_start = int(task["protein_start"])
    protein_end = int(task["protein_end"])
    cv_folds = int(task["cv_folds"])
    grid_preset = str(task["grid_preset"])
    max_param_combos = int(task["max_param_combos"])
    keep_predictions = int(task["keep_predictions"])
    out_dir = Path(args.out_dir)

    X, Y = load_aligned_data(args.data_dir, x_path=args.x_path, y_path=args.y_path)
    split_manifest = pd.read_csv(args.split_manifest, sep="\t")
    split_manifest["sample_id"] = split_manifest["sample_id"].astype(str)
    gene_lookup = build_gene_lookup(X.columns)

    if task_mode == "packed_by_protein" or str(task["model"]).lower() == "packed":
        models = [m.lower() for m in parse_csv(task.get("packed_models", "proxy,enet"))]
        normalizations = [n.lower() for n in parse_csv(task.get("packed_normalizations", "zscore"))]
        select_ks = parse_int_csv(task.get("packed_select_ks", "2000"))
        configs = []
        for model in models:
            k_values = [0] if model == "proxy" else select_ks
            for normalization in normalizations:
                for k in k_values:
                    configs.append((model, normalization, k))
        print(f"[INFO] packed task {args.task_id}: {len(configs)} configs, splits {split_start}:{split_end}, proteins {protein_start}:{protein_end}")
    else:
        configs = [(str(task["model"]).lower(), str(task["normalization"]).lower(), int(task["select_k"]))]

    for model, normalization, select_k in configs:
        run_one_config(
            task_id=int(args.task_id),
            X=X,
            Y=Y,
            split_manifest=split_manifest,
            gene_lookup=gene_lookup,
            model=model,
            normalization=normalization,
            select_k=select_k,
            split_start=split_start,
            split_end=split_end,
            protein_start=protein_start,
            protein_end=protein_end,
            cv_folds=cv_folds,
            grid_preset=grid_preset,
            max_param_combos=max_param_combos,
            keep_predictions=keep_predictions,
            out_dir=out_dir,
            seed_base=args.seed,
            n_jobs=args.n_jobs,
            skip_existing=args.skip_existing,
        )


if __name__ == "__main__":
    main()
