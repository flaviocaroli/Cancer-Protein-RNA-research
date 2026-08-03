"""Shared utilities for repeated 80/20 robustness benchmarking.

This module intentionally contains most of the implementation to keep the src/
folder small. The pipeline is evaluation-only: it saves metrics and optional
lightweight prediction shards, but does not save one model object per protein
per split.
"""
from __future__ import annotations

import itertools
import json
import math
import os
import pickle
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import joblib
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.base import clone
from sklearn.feature_selection import f_regression
from sklearn.linear_model import ElasticNet, LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold, train_test_split


# -----------------------------
# I/O and data alignment
# -----------------------------


def load_table(path: str | Path) -> pd.DataFrame:
    """Load a tabular object as a pandas DataFrame.

    Supported: joblib, pkl/pickle, parquet, csv, tsv/txt. For delimited files,
    the first column is treated as the index.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    suffix = path.suffix.lower()
    if suffix in {".joblib", ".jl"}:
        obj = joblib.load(path)
    elif suffix in {".pkl", ".pickle"}:
        with open(path, "rb") as handle:
            obj = pickle.load(handle)
    elif suffix == ".parquet":
        obj = pd.read_parquet(path)
    elif suffix in {".tsv", ".txt"}:
        obj = pd.read_csv(path, sep="\t", index_col=0)
    elif suffix == ".csv":
        obj = pd.read_csv(path, index_col=0)
    else:
        raise ValueError(f"Unsupported file suffix for {path}")

    if not isinstance(obj, pd.DataFrame):
        obj = pd.DataFrame(obj)
    obj.index = obj.index.astype(str)
    obj.columns = obj.columns.astype(str)
    return obj


def find_first_existing(candidates: Sequence[str | Path]) -> Optional[Path]:
    for candidate in candidates:
        p = Path(candidate)
        if p.exists():
            return p
    return None


def resolve_input_paths(
    data_dir: str | Path,
    x_path: Optional[str] = None,
    y_path: Optional[str] = None,
    sample_info_path: Optional[str] = None,
) -> Tuple[Path, Path, Optional[Path]]:
    """Resolve transcriptomics/proteomics paths with sensible defaults."""
    data_dir = Path(data_dir)
    x_candidates = [
        data_dir / "transcriptomics_filtered.joblib",
        data_dir / "transcriptomics.joblib",
        data_dir / "rna.joblib",
        data_dir / "X.joblib",
        data_dir / "transcriptomics_filtered.parquet",
        data_dir / "transcriptomics.parquet",
        data_dir / "rna.parquet",
        data_dir / "X.parquet",
        data_dir / "transcriptomics_filtered.tsv",
        data_dir / "transcriptomics.tsv",
        data_dir / "rna.tsv",
        data_dir / "X.tsv",
    ]
    y_candidates = [
        data_dir / "proteomics_matched.joblib",
        data_dir / "proteomics.joblib",
        data_dir / "protein.joblib",
        data_dir / "Y.joblib",
        data_dir / "proteomics_matched.parquet",
        data_dir / "proteomics.parquet",
        data_dir / "protein.parquet",
        data_dir / "Y.parquet",
        data_dir / "proteomics_matched.tsv",
        data_dir / "proteomics.tsv",
        data_dir / "protein.tsv",
        data_dir / "Y.tsv",
    ]
    sample_candidates = [
        data_dir / "sample_info.joblib",
        data_dir / "sample_info.parquet",
        data_dir / "sample_info.tsv",
        data_dir / "metadata.joblib",
        data_dir / "metadata.parquet",
        data_dir / "metadata.tsv",
    ]
    x = Path(x_path) if x_path else find_first_existing(x_candidates)
    y = Path(y_path) if y_path else find_first_existing(y_candidates)
    s = Path(sample_info_path) if sample_info_path else find_first_existing(sample_candidates)
    if x is None:
        raise FileNotFoundError(f"Could not find transcriptomics file in {data_dir}")
    if y is None:
        raise FileNotFoundError(f"Could not find proteomics file in {data_dir}")
    return x, y, s


def _deduplicate_axis(df: pd.DataFrame, axis: int, name: str) -> pd.DataFrame:
    if axis == 0:
        duplicated = df.index.duplicated(keep="first")
        if duplicated.any():
            print(f"[WARN] Dropping {duplicated.sum()} duplicated {name} rows, keeping first.")
            df = df.loc[~duplicated]
    else:
        duplicated = df.columns.duplicated(keep="first")
        if duplicated.any():
            print(f"[WARN] Dropping {duplicated.sum()} duplicated {name} columns, keeping first.")
            df = df.loc[:, ~duplicated]
    return df


def _numeric_df(df: pd.DataFrame, label: str) -> pd.DataFrame:
    out = df.apply(pd.to_numeric, errors="coerce")
    all_na_cols = out.columns[out.isna().all()].tolist()
    if all_na_cols:
        print(f"[WARN] {label}: dropping {len(all_na_cols)} all-NA columns.")
        out = out.drop(columns=all_na_cols)
    return out


def orient_and_align(X_raw: pd.DataFrame, Y_raw: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Orient transcriptomics and proteomics so rows are shared samples.

    The function tries four orientations and chooses the pair with the largest
    sample overlap between row labels.
    """
    X_options = [("X", X_raw), ("X.T", X_raw.T)]
    Y_options = [("Y", Y_raw), ("Y.T", Y_raw.T)]
    best = None
    for x_name, X in X_options:
        for y_name, Y in Y_options:
            overlap = len(pd.Index(X.index.astype(str)).intersection(pd.Index(Y.index.astype(str))))
            candidate = (overlap, x_name, y_name, X, Y)
            if best is None or candidate[0] > best[0]:
                best = candidate
    assert best is not None
    overlap, x_name, y_name, X, Y = best
    if overlap == 0:
        raise ValueError("No overlapping sample IDs found between transcriptomics and proteomics in any orientation.")
    if x_name == "X.T":
        print("[INFO] Transposing transcriptomics so rows are samples.")
    if y_name == "Y.T":
        print("[INFO] Transposing proteomics so rows are samples.")
    X = X.copy()
    Y = Y.copy()
    X.index = X.index.astype(str)
    X.columns = X.columns.astype(str)
    Y.index = Y.index.astype(str)
    Y.columns = Y.columns.astype(str)
    X = _deduplicate_axis(X, 0, "transcriptomics sample")
    X = _deduplicate_axis(X, 1, "gene")
    Y = _deduplicate_axis(Y, 0, "proteomics sample")
    Y = _deduplicate_axis(Y, 1, "protein")
    shared_samples = X.index.intersection(Y.index)
    X = X.loc[shared_samples]
    Y = Y.loc[shared_samples]
    X = _numeric_df(X, "transcriptomics")
    Y = _numeric_df(Y, "proteomics")
    print(f"[INFO] Aligned data: X={X.shape[0]} samples x {X.shape[1]} genes; Y={Y.shape[0]} samples x {Y.shape[1]} proteins")
    return X, Y


def load_aligned_data(
    data_dir: str | Path,
    x_path: Optional[str] = None,
    y_path: Optional[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    x_resolved, y_resolved, _ = resolve_input_paths(data_dir, x_path=x_path, y_path=y_path)
    print(f"[INFO] Loading transcriptomics: {x_resolved}")
    print(f"[INFO] Loading proteomics: {y_resolved}")
    X_raw = load_table(x_resolved)
    Y_raw = load_table(y_resolved)
    return orient_and_align(X_raw, Y_raw)


# -----------------------------
# Splits and task tables
# -----------------------------


def make_repeated_splits(
    sample_ids: Sequence[str],
    n_repeats: int,
    test_size: float,
    seed: int,
    stratify: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    rows = []
    sample_ids = np.array([str(s) for s in sample_ids])
    strat = None if stratify is None else np.array(stratify)
    for split_id in range(n_repeats):
        rs = seed + split_id
        train_ids, test_ids = train_test_split(
            sample_ids,
            test_size=test_size,
            random_state=rs,
            shuffle=True,
            stratify=strat,
        )
        rows.extend({"split_id": split_id, "sample_id": s, "set": "train"} for s in train_ids)
        rows.extend({"split_id": split_id, "sample_id": s, "set": "test"} for s in test_ids)
    return pd.DataFrame(rows)


def get_split_samples(split_manifest: pd.DataFrame, split_id: int) -> Tuple[List[str], List[str]]:
    sub = split_manifest.loc[split_manifest["split_id"].astype(int) == int(split_id)]
    train = sub.loc[sub["set"] == "train", "sample_id"].astype(str).tolist()
    test = sub.loc[sub["set"] == "test", "sample_id"].astype(str).tolist()
    return train, test


def chunk_ranges(total: int, chunk_size: int) -> List[Tuple[int, int]]:
    return [(start, min(start + chunk_size, total)) for start in range(0, total, chunk_size)]


# -----------------------------
# Normalization and feature selection
# -----------------------------


@dataclass
class Preprocessor:
    method: str
    impute_values: pd.Series
    mean: Optional[pd.Series] = None
    std: Optional[pd.Series] = None


def fit_preprocessor(X_train: pd.DataFrame, method: str) -> Preprocessor:
    method = method.lower()
    if method not in {"none", "zscore"}:
        raise ValueError(f"Unknown normalization method: {method}")
    impute_values = X_train.mean(axis=0, skipna=True).replace([np.inf, -np.inf], np.nan).fillna(0.0)
    if method == "zscore":
        filled = X_train.fillna(impute_values)
        mean = filled.mean(axis=0)
        std = filled.std(axis=0, ddof=0).replace(0.0, 1.0).fillna(1.0)
        return Preprocessor(method=method, impute_values=impute_values, mean=mean, std=std)
    return Preprocessor(method=method, impute_values=impute_values)


def apply_preprocessor(X: pd.DataFrame, prep: Preprocessor) -> pd.DataFrame:
    X2 = X.fillna(prep.impute_values).replace([np.inf, -np.inf], np.nan).fillna(0.0)
    if prep.method == "zscore":
        X2 = (X2 - prep.mean) / prep.std
        X2 = X2.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return X2


def select_f_regression_features(X_train: pd.DataFrame, y_train: np.ndarray, k: int) -> List[str]:
    k = min(int(k), X_train.shape[1])
    if k <= 0 or k >= X_train.shape[1]:
        return X_train.columns.tolist()
    try:
        scores, _ = f_regression(X_train.values, y_train, center=True, force_finite=True)
    except TypeError:
        scores, _ = f_regression(X_train.values, y_train, center=True)
        scores = np.nan_to_num(scores, nan=-np.inf, posinf=np.inf, neginf=-np.inf)
    order = np.argsort(scores)[::-1][:k]
    return X_train.columns[order].tolist()


def select_spearman_features(X_train: pd.DataFrame, y_train: np.ndarray, k: int) -> List[str]:
    k = min(int(k), X_train.shape[1])
    if k <= 0 or k >= X_train.shape[1]:
        return X_train.columns.tolist()
    X_rank = X_train.rank(axis=0, method="average").values.astype(float)
    y_rank = stats.rankdata(y_train).astype(float)
    Xc = X_rank - X_rank.mean(axis=0, keepdims=True)
    yc = y_rank - y_rank.mean()
    denom = np.sqrt(np.sum(Xc * Xc, axis=0) * np.sum(yc * yc))
    denom[denom == 0] = np.nan
    corr = np.abs(np.sum(Xc * yc[:, None], axis=0) / denom)
    corr = np.nan_to_num(corr, nan=-np.inf, posinf=np.inf, neginf=-np.inf)
    order = np.argsort(corr)[::-1][:k]
    return X_train.columns[order].tolist()


# -----------------------------
# Metrics
# -----------------------------


def safe_pearson(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    if len(y_true) < 3 or np.nanstd(y_true) == 0 or np.nanstd(y_pred) == 0:
        return np.nan
    return float(stats.pearsonr(y_true, y_pred)[0])


def safe_spearman(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    if len(y_true) < 3 or np.nanstd(y_true) == 0 or np.nanstd(y_pred) == 0:
        return np.nan
    return float(stats.spearmanr(y_true, y_pred).correlation)


def safe_cosine(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    denom = np.linalg.norm(y_true) * np.linalg.norm(y_pred)
    if denom == 0:
        return np.nan
    return float(np.dot(y_true, y_pred) / denom)


def compute_metrics(y_true: Sequence[float], y_pred: Sequence[float]) -> Dict[str, float]:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    y_true = y_true[mask]
    y_pred = y_pred[mask]
    out: Dict[str, float] = {"n_test": int(len(y_true))}
    if len(y_true) == 0:
        for key in ["pearson", "spearman", "mse", "rmse", "nrmse", "r2", "cosine", "calibration_slope", "calibration_intercept", "pred_sd", "obs_sd", "sd_ratio"]:
            out[key] = np.nan
        return out
    mse = mean_squared_error(y_true, y_pred)
    rmse = math.sqrt(mse)
    obs_sd = float(np.std(y_true, ddof=1)) if len(y_true) > 1 else np.nan
    pred_sd = float(np.std(y_pred, ddof=1)) if len(y_pred) > 1 else np.nan
    out["pearson"] = safe_pearson(y_true, y_pred)
    out["spearman"] = safe_spearman(y_true, y_pred)
    out["mse"] = float(mse)
    out["rmse"] = float(rmse)
    out["nrmse"] = float(rmse / obs_sd) if np.isfinite(obs_sd) and obs_sd > 0 else np.nan
    try:
        out["r2"] = float(r2_score(y_true, y_pred))
    except Exception:
        out["r2"] = np.nan
    out["cosine"] = safe_cosine(y_true, y_pred)
    if np.isfinite(pred_sd) and pred_sd > 0 and len(y_true) >= 3:
        slope, intercept = np.polyfit(y_pred, y_true, 1)
        out["calibration_slope"] = float(slope)
        out["calibration_intercept"] = float(intercept)
    else:
        out["calibration_slope"] = np.nan
        out["calibration_intercept"] = np.nan
    out["pred_sd"] = pred_sd
    out["obs_sd"] = obs_sd
    out["sd_ratio"] = float(pred_sd / obs_sd) if np.isfinite(obs_sd) and obs_sd > 0 else np.nan
    return out


# -----------------------------
# Same-gene proxy helpers
# -----------------------------


def clean_gene_like(x: str) -> str:
    x = str(x)
    for sep in ["|", ";", ",", " ", "_"]:
        if sep in x:
            x = x.split(sep)[0]
    # remove common version suffixes but keep gene symbols with hyphens
    if "." in x:
        left, right = x.rsplit(".", 1)
        if right.isdigit():
            x = left
    return x.upper()


def build_gene_lookup(gene_columns: Sequence[str]) -> Dict[str, str]:
    lookup = {}
    for col in gene_columns:
        key = clean_gene_like(col)
        if key not in lookup:
            lookup[key] = str(col)
    return lookup


def same_gene_column(protein: str, gene_lookup: Dict[str, str]) -> Optional[str]:
    key = clean_gene_like(protein)
    if key in gene_lookup:
        return gene_lookup[key]
    return None


def linear_proxy_predict(X_train: pd.DataFrame, X_test: pd.DataFrame, y_train: np.ndarray, gene_col: str) -> np.ndarray:
    xtr = X_train[[gene_col]].values.astype(float)
    xte = X_test[[gene_col]].values.astype(float)
    model = LinearRegression()
    model.fit(xtr, y_train)
    return model.predict(xte)


# -----------------------------
# Model grids and fitting
# -----------------------------


def enet_param_grid(preset: str = "fast") -> List[Dict[str, float]]:
    if preset == "tiny":
        alphas = [0.01, 0.1]
        l1_ratios = [0.5]
    elif preset == "wide":
        alphas = [1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0]
        l1_ratios = [0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95]
    else:
        alphas = [1e-3, 1e-2, 1e-1, 1.0]
        l1_ratios = [0.1, 0.5, 0.9]
    return [{"alpha": a, "l1_ratio": l} for a, l in itertools.product(alphas, l1_ratios)]


def xgb_param_grid(preset: str = "fast") -> List[Dict[str, float]]:
    if preset == "tiny":
        raw = {
            "n_estimators": [100],
            "max_depth": [2],
            "learning_rate": [0.1],
            "subsample": [0.8],
            "colsample_bytree": [0.5],
            "min_child_weight": [1],
            "reg_lambda": [1.0],
        }
    elif preset == "wide":
        raw = {
            "n_estimators": [100, 300, 600],
            "max_depth": [2, 3, 4],
            "learning_rate": [0.03, 0.05, 0.1],
            "subsample": [0.6, 0.8, 1.0],
            "colsample_bytree": [0.3, 0.5, 0.8],
            "min_child_weight": [1, 5, 10],
            "reg_lambda": [1.0, 5.0, 10.0],
        }
    else:
        raw = {
            "n_estimators": [100, 300, 600],
            "max_depth": [2, 3],
            "learning_rate": [0.03, 0.1],
            "subsample": [0.8],
            "colsample_bytree": [0.5, 0.8],
            "min_child_weight": [1, 5],
            "reg_lambda": [1.0, 5.0],
        }
    keys = list(raw.keys())
    return [dict(zip(keys, values)) for values in itertools.product(*[raw[k] for k in keys])]


def subset_param_grid(grid: List[Dict[str, float]], max_param_combos: int, seed: int) -> List[Dict[str, float]]:
    if max_param_combos <= 0 or len(grid) <= max_param_combos:
        return grid
    rng = random.Random(seed)
    return rng.sample(grid, max_param_combos)


def valid_xy_for_protein(X: pd.DataFrame, y: pd.Series, samples: Sequence[str]) -> Tuple[pd.DataFrame, np.ndarray, List[str]]:
    samples = [s for s in samples if s in X.index and s in y.index]
    y_sub = y.loc[samples]
    keep = y_sub.replace([np.inf, -np.inf], np.nan).notna()
    kept_samples = y_sub.index[keep].astype(str).tolist()
    return X.loc[kept_samples], y_sub.loc[kept_samples].values.astype(float), kept_samples


def tune_and_predict_enet(
    X: pd.DataFrame,
    y: pd.Series,
    train_samples: Sequence[str],
    test_samples: Sequence[str],
    normalization: str,
    select_k: int,
    cv_folds: int,
    grid_preset: str,
    max_param_combos: int,
    seed: int,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, object]]:
    X_train_raw, y_train, kept_train = valid_xy_for_protein(X, y, train_samples)
    X_test_raw, y_test, kept_test = valid_xy_for_protein(X, y, test_samples)
    info: Dict[str, object] = {"n_train": len(y_train), "n_test": len(y_test), "status": "ok"}
    if len(y_train) < max(cv_folds, 10) or len(y_test) < 3 or np.std(y_train) == 0:
        info["status"] = "too_few_samples_or_constant_y"
        info["test_samples"] = kept_test
        return y_test, np.full(len(y_test), np.nan), info

    n_folds = min(cv_folds, len(y_train))
    kfold = KFold(n_splits=n_folds, shuffle=True, random_state=seed)
    grid = subset_param_grid(enet_param_grid(grid_preset), max_param_combos, seed)
    fold_scores = {json.dumps(p, sort_keys=True): [] for p in grid}

    for tr_idx, val_idx in kfold.split(X_train_raw):
        X_tr_raw = X_train_raw.iloc[tr_idx]
        X_val_raw = X_train_raw.iloc[val_idx]
        y_tr = y_train[tr_idx]
        y_val = y_train[val_idx]
        prep = fit_preprocessor(X_tr_raw, normalization)
        X_tr = apply_preprocessor(X_tr_raw, prep)
        X_val = apply_preprocessor(X_val_raw, prep)
        selected = select_f_regression_features(X_tr, y_tr, select_k)
        X_tr_sel = X_tr[selected].values
        X_val_sel = X_val[selected].values
        for params in grid:
            model = ElasticNet(
                alpha=float(params["alpha"]),
                l1_ratio=float(params["l1_ratio"]),
                max_iter=10000,
                random_state=seed,
                selection="cyclic",
            )
            try:
                model.fit(X_tr_sel, y_tr)
                pred = model.predict(X_val_sel)
                fold_scores[json.dumps(params, sort_keys=True)].append(mean_squared_error(y_val, pred))
            except Exception:
                fold_scores[json.dumps(params, sort_keys=True)].append(np.inf)

    mean_scores = {k: float(np.nanmean(v)) for k, v in fold_scores.items()}
    best_key = min(mean_scores, key=mean_scores.get)
    best_params = json.loads(best_key)

    prep = fit_preprocessor(X_train_raw, normalization)
    X_train = apply_preprocessor(X_train_raw, prep)
    X_test = apply_preprocessor(X_test_raw, prep)
    selected = select_f_regression_features(X_train, y_train, select_k)
    model = ElasticNet(
        alpha=float(best_params["alpha"]),
        l1_ratio=float(best_params["l1_ratio"]),
        max_iter=10000,
        random_state=seed,
        selection="cyclic",
    )
    model.fit(X_train[selected].values, y_train)
    y_pred = model.predict(X_test[selected].values)
    info.update(
        {
            "best_params_json": json.dumps(best_params, sort_keys=True),
            "best_cv_mse": mean_scores[best_key],
            "n_selected_features": len(selected),
            "n_nonzero_features": int(np.sum(np.abs(model.coef_) > 1e-12)),
            "test_samples": kept_test,
        }
    )
    return y_test, y_pred, info


def tune_and_predict_xgb(
    X: pd.DataFrame,
    y: pd.Series,
    train_samples: Sequence[str],
    test_samples: Sequence[str],
    normalization: str,
    select_k: int,
    cv_folds: int,
    grid_preset: str,
    max_param_combos: int,
    seed: int,
    n_jobs: int = 1,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, object]]:
    try:
        from xgboost import XGBRegressor
    except Exception as exc:
        X_test_raw, y_test, kept_test = valid_xy_for_protein(X, y, test_samples)
        return y_test, np.full(len(y_test), np.nan), {"status": f"xgboost_import_failed: {exc}", "n_train": 0, "n_test": len(y_test), "test_samples": kept_test}

    X_train_raw, y_train, kept_train = valid_xy_for_protein(X, y, train_samples)
    X_test_raw, y_test, kept_test = valid_xy_for_protein(X, y, test_samples)
    info: Dict[str, object] = {"n_train": len(y_train), "n_test": len(y_test), "status": "ok"}
    if len(y_train) < max(cv_folds, 10) or len(y_test) < 3 or np.std(y_train) == 0:
        info["status"] = "too_few_samples_or_constant_y"
        info["test_samples"] = kept_test
        return y_test, np.full(len(y_test), np.nan), info

    n_folds = min(cv_folds, len(y_train))
    kfold = KFold(n_splits=n_folds, shuffle=True, random_state=seed)
    grid = subset_param_grid(xgb_param_grid(grid_preset), max_param_combos, seed)
    fold_scores = {json.dumps(p, sort_keys=True): [] for p in grid}

    for tr_idx, val_idx in kfold.split(X_train_raw):
        X_tr_raw = X_train_raw.iloc[tr_idx]
        X_val_raw = X_train_raw.iloc[val_idx]
        y_tr = y_train[tr_idx]
        y_val = y_train[val_idx]
        prep = fit_preprocessor(X_tr_raw, normalization)
        X_tr = apply_preprocessor(X_tr_raw, prep)
        X_val = apply_preprocessor(X_val_raw, prep)
        selected = select_spearman_features(X_tr, y_tr, select_k)
        X_tr_sel = X_tr[selected].values
        X_val_sel = X_val[selected].values
        for params in grid:
            model = XGBRegressor(
                objective="reg:squarederror",
                eval_metric="rmse",
                random_state=seed,
                n_jobs=n_jobs,
                verbosity=0,
                **params,
            )
            try:
                model.fit(X_tr_sel, y_tr)
                pred = model.predict(X_val_sel)
                fold_scores[json.dumps(params, sort_keys=True)].append(mean_squared_error(y_val, pred))
            except Exception:
                fold_scores[json.dumps(params, sort_keys=True)].append(np.inf)

    mean_scores = {k: float(np.nanmean(v)) for k, v in fold_scores.items()}
    best_key = min(mean_scores, key=mean_scores.get)
    best_params = json.loads(best_key)

    prep = fit_preprocessor(X_train_raw, normalization)
    X_train = apply_preprocessor(X_train_raw, prep)
    X_test = apply_preprocessor(X_test_raw, prep)
    selected = select_spearman_features(X_train, y_train, select_k)
    model = XGBRegressor(
        objective="reg:squarederror",
        eval_metric="rmse",
        random_state=seed,
        n_jobs=n_jobs,
        verbosity=0,
        **best_params,
    )
    model.fit(X_train[selected].values, y_train)
    y_pred = model.predict(X_test[selected].values)
    info.update(
        {
            "best_params_json": json.dumps(best_params, sort_keys=True),
            "best_cv_mse": mean_scores[best_key],
            "n_selected_features": len(selected),
            "n_nonzero_features": np.nan,
            "test_samples": kept_test,
        }
    )
    return y_test, y_pred, info


# -----------------------------
# Saving and aggregation helpers
# -----------------------------


def config_name(model: str, normalization: str, select_k: int) -> str:
    return f"model={model}_norm={normalization}_k={int(select_k)}"


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def write_parquet_or_csv(df: pd.DataFrame, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".parquet":
        try:
            df.to_parquet(path, index=False)
            return path
        except Exception as exc:
            fallback = path.with_suffix(".csv.gz")
            print(f"[WARN] Could not write parquet ({exc}); writing {fallback}")
            df.to_csv(fallback, index=False)
            return fallback
    df.to_csv(path, index=False)
    return path


def read_parquet_or_csv(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path)


def list_metric_shards(metrics_dir: str | Path) -> List[Path]:
    metrics_dir = Path(metrics_dir)
    return sorted(list(metrics_dir.glob("*.parquet")) + list(metrics_dir.glob("*.csv.gz")) + list(metrics_dir.glob("*.csv")))


def bootstrap_ci(values: np.ndarray, n_boot: int, seed: int, statistic: str = "median") -> Tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    if statistic == "mean":
        center = float(np.mean(values))
        boot = [np.mean(rng.choice(values, size=len(values), replace=True)) for _ in range(n_boot)]
    else:
        center = float(np.median(values))
        boot = [np.median(rng.choice(values, size=len(values), replace=True)) for _ in range(n_boot)]
    low, high = np.nanpercentile(boot, [2.5, 97.5])
    return center, float(low), float(high)
