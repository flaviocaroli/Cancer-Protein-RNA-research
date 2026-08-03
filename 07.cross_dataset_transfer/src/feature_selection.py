"""
feature_selection.py -- per-protein univariate feature selection.

ENet : ANOVA / F-regression, keep top-k by F score          (default k=5000)
XGB  : absolute Spearman correlation, keep top-k             (default k=5000)

Both are fit on the TRAINING union only (labels = that protein's abundance),
then the identical gene list is applied to the held-out test dataset. Selection
happens once per protein per scenario on the full train union; inside the 5-fold
CV the same train-only selection discipline is preserved by refitting on each CV
train fold (see models.py). Genes with zero/constant variance are dropped first.
"""
from __future__ import annotations

from typing import List

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.feature_selection import f_regression


def _finite_variable_genes(X: pd.DataFrame) -> List[str]:
    A = X.to_numpy(np.float64)
    ok = np.isfinite(A).all(axis=0) & (np.nanstd(A, axis=0) > 1e-12)
    return list(X.columns[ok])


def select_anova_f(X_train: pd.DataFrame, y_train: np.ndarray, k: int = 5000) -> List[str]:
    """Top-k genes by ANOVA F-regression score (ENet selector)."""
    genes = _finite_variable_genes(X_train)
    if not genes:
        return []
    Xg = X_train[genes].to_numpy(np.float64)
    y = np.asarray(y_train, np.float64)
    with np.errstate(all="ignore"):
        F, _ = f_regression(Xg, y, center=True)
    F = np.nan_to_num(F, nan=0.0, posinf=0.0, neginf=0.0)
    k = min(k, len(genes))
    top = np.argsort(F)[::-1][:k]
    return [genes[i] for i in top]


def select_spearman(X_train: pd.DataFrame, y_train: np.ndarray, k: int = 5000) -> List[str]:
    """Top-k genes by absolute Spearman correlation with the target (XGB selector)."""
    genes = _finite_variable_genes(X_train)
    if not genes:
        return []
    Xg = X_train[genes].to_numpy(np.float64)
    y = np.asarray(y_train, np.float64)
    yr = rankdata(y)
    yr = (yr - yr.mean())
    yden = np.sqrt((yr ** 2).sum()) + 1e-12
    # rank each gene column, then Pearson on ranks == Spearman
    Xr = np.apply_along_axis(rankdata, 0, Xg)
    Xr = Xr - Xr.mean(axis=0, keepdims=True)
    xden = np.sqrt((Xr ** 2).sum(axis=0)) + 1e-12
    rho = (Xr * yr[:, None]).sum(axis=0) / (xden * yden)
    rho = np.nan_to_num(np.abs(rho), nan=0.0)
    k = min(k, len(genes))
    top = np.argsort(rho)[::-1][:k]
    return [genes[i] for i in top]


def select_features(X_train, y_train, model_family: str, k: int = 5000) -> List[str]:
    if model_family == "enet":
        return select_anova_f(X_train, y_train, k)
    if model_family == "xgb":
        return select_spearman(X_train, y_train, k)
    raise ValueError(model_family)
