#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import traceback
from pathlib import Path
import sys

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNet

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from utils import (
    compute_metrics,
    domain_balanced_sample_weights,
    dump_json,
    fit_x_normalizer,
    fit_y_transformer,
    load_yaml,
    make_cv_splits,
    parameter_grid,
    safe_name,
    score_predictions,
    select_feature_indices,
    transform_x,
    transform_y,
    write_table,
)


def build_model(model_type: str, params: dict, n_jobs: int, seed: int):
    if model_type == "enet":
        return ElasticNet(
            alpha=float(params["alpha"]),
            l1_ratio=float(params["l1_ratio"]),
            max_iter=int(params.get("max_iter", 20000)),
            tol=float(params.get("tol", 1e-4)),
            random_state=seed,
            selection="cyclic",
        )
    if model_type == "xgb":
        try:
            from xgboost import XGBRegressor
        except ImportError as exc:
            raise ImportError("xgboost is required for MODEL=xgb") from exc
        return XGBRegressor(
            objective="reg:squarederror",
            random_state=seed,
            n_jobs=n_jobs,
            tree_method="hist",
            verbosity=0,
            **params,
        )
    raise ValueError(model_type)


def fit_with_weights(model, X, y, weights):
    try:
        model.fit(X, y, sample_weight=weights)
    except TypeError:
        model.fit(X, y)
    return model


def cv_select_params(
    model_type: str,
    X_raw: pd.DataFrame,
    y_raw: pd.Series,
    meta: pd.DataFrame,
    prep_manifest: dict,
    train_cfg: dict,
    param_combos: list[dict],
    n_jobs: int,
    protein_seed: int,
) -> tuple[dict, pd.DataFrame]:
    domains = meta["dataset"].astype(str).to_numpy()
    cohorts = meta["cohort"].astype(str).to_numpy() if "cohort" in meta else None
    splits = make_cv_splits(domains, cohorts, train_cfg["cv_folds"], train_cfg["cv_mode"], protein_seed)
    records = []
    best_score = -np.inf
    best_params = None

    for combo_id, params in enumerate(param_combos):
        fold_scores = []
        for fold_id, (tr, va) in enumerate(splits):
            Xtr_raw, Xva_raw = X_raw.iloc[tr], X_raw.iloc[va]
            ytr_raw, yva_raw = y_raw.iloc[tr], y_raw.iloc[va]
            dtr = meta.iloc[tr]["dataset"]
            dva = meta.iloc[va]["dataset"]

            x_bundle = fit_x_normalizer(
                Xtr_raw,
                dtr,
                method=prep_manifest["x_normalization"],
                n_bins=int(prep_manifest["n_bins"]),
                eps=float(prep_manifest["eps"]),
                domain_balanced=bool(prep_manifest["domain_balanced_normalizer"]),
                quantile_fit_max_samples_per_domain=1000,
                seed=protein_seed + fold_id,
            )
            Xtr = transform_x(Xtr_raw, x_bundle)
            Xva = transform_x(Xva_raw, x_bundle)

            ytr_df = ytr_raw.to_frame(name=y_raw.name)
            yva_df = yva_raw.to_frame(name=y_raw.name)
            y_bundle = fit_y_transformer(ytr_df, dtr, prep_manifest["y_normalization"], eps=float(prep_manifest["eps"]))
            ytr = transform_y(ytr_df, dtr, y_bundle).iloc[:, 0].to_numpy(dtype=float)
            yva = transform_y(yva_df, dva, y_bundle).iloc[:, 0].to_numpy(dtype=float)

            selected_idx, _ = select_feature_indices(
                Xtr.to_numpy(dtype=np.float32),
                ytr,
                dtr.astype(str).to_numpy(),
                model_type=model_type,
                k=int(train_cfg["feature_k"]),
                score_mode=train_cfg["feature_score_mode"],
            )
            weights = domain_balanced_sample_weights(dtr) if train_cfg["domain_weighting"] == "balanced" else np.ones(len(tr))
            model = build_model(model_type, params, n_jobs=n_jobs, seed=protein_seed + combo_id * 101 + fold_id)
            fit_with_weights(model, Xtr.to_numpy()[:, selected_idx], ytr, weights)
            pred = model.predict(Xva.to_numpy()[:, selected_idx])
            fold_score = score_predictions(yva, pred, dva.astype(str).to_numpy(), train_cfg["cv_metric"])
            fold_scores.append(fold_score)

        mean_score = float(np.mean(fold_scores))
        records.append({
            "combo_id": combo_id,
            "params_json": json.dumps(params, sort_keys=True),
            "mean_cv_score": mean_score,
            "sd_cv_score": float(np.std(fold_scores, ddof=0)),
            "fold_scores_json": json.dumps([float(x) for x in fold_scores]),
        })
        if mean_score > best_score:
            best_score = mean_score
            best_params = params

    if best_params is None:
        raise RuntimeError("No valid hyperparameter combination")
    return best_params, pd.DataFrame(records)




def load_prepared_data(prepared_dir: Path) -> tuple[dict, dict]:
    manifest = json.loads((prepared_dir / "manifest.json").read_text())
    data = {
        "X_raw_all": joblib.load(prepared_dir / "train_X_raw.joblib"),
        "X_norm_all": joblib.load(prepared_dir / "train_X_normalized.joblib"),
        "Y_raw_all": joblib.load(prepared_dir / "train_Y_raw.joblib"),
        "Y_norm_all": joblib.load(prepared_dir / "train_Y_normalized.joblib"),
    }
    data["meta_all"] = pd.read_csv(prepared_dir / "train_metadata.tsv", sep="\t", index_col=0).reindex(data["X_raw_all"].index)
    tests = {}
    for test_name in manifest["test_datasets"]:
        test_dir = prepared_dir / "tests" / test_name
        tests[test_name] = {
            "X_norm": joblib.load(test_dir / "X_normalized.joblib"),
            "Y_raw": joblib.load(test_dir / "Y_raw.joblib"),
            "Y_eval": joblib.load(test_dir / "Y_evaluation_scaled.joblib"),
        }
    data["tests"] = tests
    return manifest, data

def train_one_protein(
    protein: str,
    model_type: str,
    prepared_dir: Path,
    manifest: dict,
    data: dict,
    model_cfg: dict,
    n_jobs: int,
    seed: int,
    save_predictions: bool,
    save_all_feature_rows: bool,
) -> tuple[list[dict], list[dict], list[dict], dict]:
    X_raw_all: pd.DataFrame = data["X_raw_all"]
    X_norm_all: pd.DataFrame = data["X_norm_all"]
    Y_raw_all: pd.DataFrame = data["Y_raw_all"]
    Y_norm_all: pd.DataFrame = data["Y_norm_all"]
    meta_all: pd.DataFrame = data["meta_all"]

    valid = Y_raw_all[protein].notna() & Y_norm_all[protein].notna()
    X_raw = X_raw_all.loc[valid]
    X_norm = X_norm_all.loc[valid]
    y_raw = Y_raw_all.loc[valid, protein]
    y_norm = Y_norm_all.loc[valid, protein].to_numpy(dtype=float)
    meta = meta_all.loc[valid]
    domains = meta["dataset"].astype(str).to_numpy()

    if len(y_raw) < max(25, int(model_cfg["training"]["cv_folds"]) * 4) or np.nanstd(y_norm) == 0:
        raise ValueError(f"Insufficient/non-variable training observations for {protein}: n={len(y_raw)}")

    train_cfg = dict(model_cfg["training"])
    grid = dict(train_cfg[f"{model_type}_grid"])
    combos = parameter_grid(grid, int(train_cfg["max_param_combos"]), seed)
    best_params, cv_df = cv_select_params(
        model_type=model_type,
        X_raw=X_raw,
        y_raw=y_raw,
        meta=meta,
        prep_manifest=manifest,
        train_cfg=train_cfg,
        param_combos=combos,
        n_jobs=n_jobs,
        protein_seed=seed,
    )

    selected_idx, selected_scores = select_feature_indices(
        X_norm.to_numpy(dtype=np.float32),
        y_norm,
        domains,
        model_type=model_type,
        k=int(train_cfg["feature_k"]),
        score_mode=train_cfg["feature_score_mode"],
    )
    selected_features = X_norm.columns[selected_idx].astype(str).tolist()
    weights = domain_balanced_sample_weights(meta["dataset"]) if train_cfg["domain_weighting"] == "balanced" else np.ones(len(meta))
    model = build_model(model_type, best_params, n_jobs=n_jobs, seed=seed)
    fit_with_weights(model, X_norm.to_numpy()[:, selected_idx], y_norm, weights)

    metric_rows: list[dict] = []
    prediction_rows: list[dict] = []
    coefficient_rows: list[dict] = []

    for test_name in manifest["test_datasets"]:
        X_test: pd.DataFrame = data["tests"][test_name]["X_norm"]
        Y_test_raw: pd.DataFrame = data["tests"][test_name]["Y_raw"]
        Y_test_eval: pd.DataFrame = data["tests"][test_name]["Y_eval"]
        pred = model.predict(X_test[selected_features].to_numpy(dtype=np.float32)).astype(np.float32)
        y_raw_test = Y_test_raw[protein].to_numpy(dtype=float)
        y_eval_test = Y_test_eval[protein].to_numpy(dtype=float)
        metrics = compute_metrics(y_raw_test, y_eval_test, pred)
        metric_rows.append({
            "scenario": manifest["scenario"],
            "model": model_type,
            "x_normalization": manifest["x_normalization"],
            "y_normalization": manifest["y_normalization"],
            "protein": protein,
            "test_dataset": test_name,
            "n_train": len(y_raw),
            "n_selected_features": len(selected_features),
            "best_params_json": json.dumps(best_params, sort_keys=True),
            "best_cv_score": float(cv_df["mean_cv_score"].max()),
            **metrics,
        })
        if save_predictions:
            prediction_rows.extend({
                "scenario": manifest["scenario"],
                "model": model_type,
                "x_normalization": manifest["x_normalization"],
                "y_normalization": manifest["y_normalization"],
                "test_dataset": test_name,
                "sample_id": sample_id,
                "protein": protein,
                "y_true_raw": float(yr) if np.isfinite(yr) else np.nan,
                "y_true_evaluation_scaled": float(ye) if np.isfinite(ye) else np.nan,
                "y_pred": float(yp),
            } for sample_id, yr, ye, yp in zip(X_test.index, y_raw_test, y_eval_test, pred))

    if model_type == "enet":
        coeffs = np.asarray(model.coef_, dtype=np.float32)
        for gene, score, coef in zip(selected_features, selected_scores, coeffs):
            if save_all_feature_rows or coef != 0:
                coefficient_rows.append({
                    "scenario": manifest["scenario"],
                    "protein": protein,
                    "gene": gene,
                    "coefficient": float(coef),
                    "abs_coefficient": abs(float(coef)),
                    "is_nonzero": bool(coef != 0),
                    "feature_score": float(score),
                })
        model_details = {
            "intercept": float(model.intercept_),
            "coefficients": coeffs,
            "n_nonzero_coefficients": int(np.count_nonzero(coeffs)),
        }
    else:
        importances = np.asarray(model.feature_importances_, dtype=np.float32)
        for gene, score, importance in zip(selected_features, selected_scores, importances):
            if save_all_feature_rows or importance > 0:
                coefficient_rows.append({
                    "scenario": manifest["scenario"],
                    "protein": protein,
                    "gene": gene,
                    "importance": float(importance),
                    "is_nonzero": bool(importance > 0),
                    "feature_score": float(score),
                })
        model_details = {
            "feature_importances": importances,
            "n_nonzero_importances": int(np.count_nonzero(importances)),
        }

    model_artifact = {
        "scenario": manifest["scenario"],
        "model_type": model_type,
        "protein": protein,
        "model": model,
        "selected_features": selected_features,
        "selected_feature_scores": selected_scores,
        "best_params": best_params,
        "cv_results": cv_df,
        "x_normalization": manifest["x_normalization"],
        "y_normalization": manifest["y_normalization"],
        "shared_normalization_bundle": str((prepared_dir / "normalization/x_normalizer.joblib").resolve()),
        "shared_target_transformer": str((prepared_dir / "normalization/y_train_transformer.joblib").resolve()),
        "prepared_manifest": manifest,
        "training_datasets": manifest["train_datasets"],
        "test_datasets": manifest["test_datasets"],
        "feature_k": int(train_cfg["feature_k"]),
        "feature_score_mode": train_cfg["feature_score_mode"],
        "cv_metric": train_cfg["cv_metric"],
        "target_interpretation": manifest["important_target_note"],
        "evaluation_metrics": metric_rows,
        "heldout_predictions": prediction_rows if save_predictions else [],
        "feature_table": coefficient_rows,
        **model_details,
    }
    return metric_rows, prediction_rows, coefficient_rows, model_artifact


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--prepared-dir", required=True)
    p.add_argument("--results-dir", required=True)
    p.add_argument("--task-table", required=True)
    p.add_argument("--task-id", type=int, required=True)
    p.add_argument("--model", choices=["enet", "xgb"], required=True)
    p.add_argument("--model-config", required=True)
    p.add_argument("--protein-subset", default="")
    p.add_argument("--n-jobs", type=int, default=1)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--save-predictions", type=int, choices=[0, 1], default=1)
    p.add_argument("--save-all-feature-rows", type=int, choices=[0, 1], default=0,
                   help="Write all 5,000 feature rows to shard tables. All coefficients/importances are always stored in each model artifact.")
    p.add_argument("--skip-existing", type=int, choices=[0, 1], default=1)
    args = p.parse_args()

    prepared_dir = Path(args.prepared_dir).resolve()
    results_dir = Path(args.results_dir).resolve()
    task_table = pd.read_csv(args.task_table, sep="\t")
    row = task_table.loc[task_table.task_id.eq(args.task_id)].iloc[0]
    model_cfg = load_yaml(args.model_config)
    manifest, prepared_data = load_prepared_data(prepared_dir)

    proteins = pd.read_csv(prepared_dir / "protein_universe.tsv", sep="\t")["protein"].astype(str).tolist()
    if args.protein_subset:
        sub_path = Path(args.protein_subset)
        sub = pd.read_csv(sub_path, sep=None, engine="python")
        col = "protein" if "protein" in sub.columns else sub.columns[0]
        wanted = set(sub[col].astype(str))
        proteins = [x for x in proteins if x in wanted]
    proteins = proteins[int(row.protein_start):int(row.protein_end)]

    shard_dir = results_dir / "shards"
    model_dir = results_dir / "models_per_protein" / args.model
    log_dir = results_dir / "protein_status" / args.model
    for d in [shard_dir, model_dir, log_dir]:
        d.mkdir(parents=True, exist_ok=True)

    metrics_all, preds_all, coefs_all, status_all = [], [], [], []
    for offset, protein in enumerate(proteins):
        protein_seed = args.seed + args.task_id * 100000 + offset
        model_path = model_dir / f"{safe_name(protein)}.joblib"
        if args.skip_existing and model_path.exists():
            try:
                existing = joblib.load(model_path)
                metrics_all.extend(existing.get("evaluation_metrics", []))
                if args.save_predictions:
                    preds_all.extend(existing.get("heldout_predictions", []))
                coefs_all.extend(existing.get("feature_table", []))
                status_all.append({"protein": protein, "status": "skipped_existing_recovered", "model_path": str(model_path)})
                continue
            except Exception as exc:
                print(f"[WARN] Could not recover existing artifact for {protein}: {exc}; retraining.", flush=True)
        try:
            metric_rows, prediction_rows, coefficient_rows, artifact = train_one_protein(
                protein=protein,
                model_type=args.model,
                prepared_dir=prepared_dir,
                manifest=manifest,
                data=prepared_data,
                model_cfg=model_cfg,
                n_jobs=args.n_jobs,
                seed=protein_seed,
                save_predictions=bool(args.save_predictions),
                save_all_feature_rows=bool(args.save_all_feature_rows),
            )
            joblib.dump(artifact, model_path, compress=3)
            metrics_all.extend(metric_rows)
            preds_all.extend(prediction_rows)
            coefs_all.extend(coefficient_rows)
            status_all.append({"protein": protein, "status": "completed", "model_path": str(model_path)})
            print(f"[OK] {protein}", flush=True)
        except Exception as exc:
            status_all.append({
                "protein": protein,
                "status": "failed",
                "error_type": type(exc).__name__,
                "error": str(exc),
                "traceback": traceback.format_exc(limit=8),
            })
            print(f"[FAIL] {protein}: {type(exc).__name__}: {exc}", flush=True)

    tag = f"{args.model}_task_{args.task_id:03d}"
    if metrics_all:
        write_table(pd.DataFrame(metrics_all), shard_dir / f"metrics_{tag}.parquet")
    if args.save_predictions and preds_all:
        write_table(pd.DataFrame(preds_all), shard_dir / f"predictions_{tag}.parquet")
    if coefs_all:
        write_table(pd.DataFrame(coefs_all), shard_dir / f"features_{tag}.parquet")
    pd.DataFrame(status_all).to_csv(log_dir / f"status_task_{args.task_id:03d}.tsv", sep="\t", index=False)
    dump_json({"task_id": args.task_id, "model": args.model, "n_requested": len(proteins), "n_status": len(status_all)}, log_dir / f"task_{args.task_id:03d}.json")

if __name__ == "__main__":
    main()
