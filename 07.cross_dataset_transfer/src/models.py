"""
models.py -- per-protein training with inner 5-fold CV, and a self-contained
artifact per protein that carries everything needed to predict on another dataset.

Key design points requested:
  * Train on PROTEIN accuracy: one model per protein, target = that protein.
  * ENet: save gene coefficients (feat_coeffs_dict) + intercept, matching the
    `03. elasticNet` artifact so `_predict_with_coef_dict`-style code still works.
  * XGB: save the fitted booster.
  * Both: save selected_features (order), the fitted feature-normalizer, best CV
    hyper-parameters, and provenance (train/test datasets, k, method).
  * Feature selection + normalization are refit INSIDE each CV fold on the fold's
    training rows only, so hyper-parameter scores are leakage-safe. The final
    model is refit on the full train union with selection/normalization fit on
    all train rows.

Artifact schema (joblib dict), per protein:
{
  "protein": str,
  "model_family": "enet"|"xgb",
  "scenario": str,
  "train_datasets": [..], "test_datasets": [..],
  "select_k": int, "selector": "anova_f"|"spearman",
  "normalization": "zscore"|"binning"|"ranking"|"none",
  "normalizer": {...},                 # from normalization.fit_feature_normalizer
  "selected_features": [gene,...],     # order matches normalizer columns
  # ENet only:
  "features_coef_dict": {gene: coef},  # <-- the coefficients you asked to keep
  "intercept": float,
  "best_params": {"alpha":.., "l1_ratio":..},
  # XGB only:
  "booster": <xgboost.XGBRegressor>,
  "best_params": {..},
  "cv_pearson_mean": float, "cv_pearson_std": float, "n_train": int,
}
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import KFold

from feature_selection import select_features
from normalization import fit_feature_normalizer, apply_feature_normalizer

ENET_GRID = {
    "alpha": [0.001, 0.01, 0.1, 1.0],
    "l1_ratio": [0.1, 0.5, 0.9],
}
XGB_GRID = {
    "n_estimators": [300, 600],
    "max_depth": [2, 3],
    "learning_rate": [0.03, 0.1],
    "subsample": [0.8],
    "colsample_bytree": [0.6],
    "min_child_weight": [5],
    "reg_lambda": [1.0, 5.0],
}


def _iter_param_grid(grid: Dict, max_combos: Optional[int] = None, seed: int = 0):
    import itertools
    keys = list(grid.keys())
    combos = [dict(zip(keys, vals)) for vals in itertools.product(*[grid[k] for k in keys])]
    if max_combos is not None and len(combos) > max_combos:
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(combos), size=max_combos, replace=False)
        combos = [combos[i] for i in idx]
    return combos


def _safe_pearson(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.size < 3 or np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return np.nan
    return float(pearsonr(a, b)[0])


def _cv_score_config(
    X_tr: pd.DataFrame, y_tr: np.ndarray, model_family: str, params: Dict,
    normalization: str, select_k: int, cv_folds: int, n_bins: int, seed: int,
) -> float:
    """Mean CV Pearson for one hyper-parameter config, with in-fold selection+norm."""
    kf = KFold(n_splits=cv_folds, shuffle=True, random_state=seed)
    scores = []
    for tr_idx, va_idx in kf.split(X_tr):
        Xf, Xv = X_tr.iloc[tr_idx], X_tr.iloc[va_idx]
        yf, yv = y_tr[tr_idx], y_tr[va_idx]
        sel = select_features(Xf, yf, model_family, k=select_k)
        if not sel:
            scores.append(np.nan); continue
        norm = fit_feature_normalizer(Xf[sel], method=normalization, n_bins=n_bins)
        Af = apply_feature_normalizer(Xf[sel], norm)
        Av = apply_feature_normalizer(Xv[sel], norm)
        if model_family == "enet":
            est = ElasticNet(alpha=params["alpha"], l1_ratio=params["l1_ratio"],
                             max_iter=5000, random_state=seed)
        else:
            from xgboost import XGBRegressor
            est = XGBRegressor(**params, tree_method="hist", n_jobs=1,
                               random_state=seed, verbosity=0)
        est.fit(Af, yf)
        scores.append(_safe_pearson(yv, est.predict(Av)))
    return float(np.nanmean(scores)) if np.any(np.isfinite(scores)) else np.nan


def train_protein_model(
    X_train: pd.DataFrame,
    y_train: np.ndarray,
    protein: str,
    model_family: str,
    normalization: str = "zscore",
    select_k: int = 5000,
    cv_folds: int = 5,
    n_bins: int = 50,
    max_combos: Optional[int] = 8,
    seed: int = 42,
    scenario: str = "",
    train_datasets: Optional[List[str]] = None,
    test_datasets: Optional[List[str]] = None,
) -> Optional[Dict]:
    """
    Full per-protein routine: CV-tune -> refit on full train union -> build a
    self-contained artifact (with ENet coefficients) ready to predict elsewhere.
    """
    y_train = np.asarray(y_train, float)
    ok = np.isfinite(y_train)
    if ok.sum() < max(10, cv_folds + 1):
        return None
    Xt = X_train.loc[ok] if X_train.index.equals(pd.RangeIndex(len(X_train))) is False else X_train.iloc[np.where(ok)[0]]
    Xt = X_train.iloc[np.where(ok)[0]]
    yt = y_train[ok]

    grid = ENET_GRID if model_family == "enet" else XGB_GRID
    combos = _iter_param_grid(grid, max_combos=max_combos, seed=seed)

    best_params, best_score = None, -np.inf
    for p in combos:
        s = _cv_score_config(Xt, yt, model_family, p, normalization,
                             select_k, cv_folds, n_bins, seed)
        if np.isfinite(s) and s > best_score:
            best_score, best_params = s, p
    if best_params is None:
        return None

    # Final refit on full train union: selection + normalization on ALL train rows.
    sel = select_features(Xt, yt, model_family, k=select_k)
    if not sel:
        return None
    norm = fit_feature_normalizer(Xt[sel], method=normalization, n_bins=n_bins)
    A = apply_feature_normalizer(Xt[sel], norm)

    artifact = {
        "protein": protein,
        "model_family": model_family,
        "scenario": scenario,
        "train_datasets": train_datasets or [],
        "test_datasets": test_datasets or [],
        "select_k": int(select_k),
        "selector": "anova_f" if model_family == "enet" else "spearman",
        "normalization": normalization,
        "normalizer": norm,
        "selected_features": list(map(str, sel)),
        "best_params": best_params,
        "cv_pearson_mean": float(best_score),
        "n_train": int(len(yt)),
    }

    if model_family == "enet":
        est = ElasticNet(alpha=best_params["alpha"], l1_ratio=best_params["l1_ratio"],
                         max_iter=10000, random_state=seed)
        est.fit(A, yt)
        artifact["features_coef_dict"] = {str(g): float(c) for g, c in zip(sel, est.coef_)}
        artifact["intercept"] = float(est.intercept_)
        artifact["n_nonzero_coef"] = int(np.sum(np.asarray(est.coef_) != 0))
    else:
        from xgboost import XGBRegressor
        est = XGBRegressor(**best_params, tree_method="hist", n_jobs=1,
                           random_state=seed, verbosity=0)
        est.fit(A, yt)
        artifact["booster"] = est

    return artifact


def predict_with_artifact(artifact: Dict, X_new: pd.DataFrame) -> np.ndarray:
    """
    Apply a saved per-protein artifact to a new dataset (the held-out test set).
    Uses the stored selected_features + fitted normalizer, so the test data is
    transformed with train-learned parameters only.
    """
    sel = artifact["selected_features"]
    norm = artifact["normalizer"]
    A = apply_feature_normalizer(X_new.reindex(columns=sel), norm)
    if artifact["model_family"] == "enet":
        w = np.array([artifact["features_coef_dict"][g] for g in sel], float)
        return A.astype(float) @ w + artifact["intercept"]
    return artifact["booster"].predict(A)
