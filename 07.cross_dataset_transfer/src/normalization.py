"""
normalization.py -- model-level feature normalization for the selected genes.

Stage B, inner normalizer. Operates on the (already harmonized) selected-feature
matrix that goes into ENet/XGB. Three methods, each fit on the TRAIN union and
applied unchanged to the held-out test dataset:

    zscore   : per-gene (train mean, train std). Stored as mean/std -> matches
               the `zscore` branch of utils_val.apply_artifact_normalizer.
    binning  : per-gene quantile bins; edges learned on train.
    ranking  : per-gene rank mapped to [0,1] via train sorted values.

The stored normalizer dicts use the SAME keys your 05.validation
apply_artifact_normalizer already understands (method/columns/mean/std,
edges, sorted_values), so trained artifacts stay compatible with existing code.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

EPS = 1e-8


def fit_feature_normalizer(
    X_train: pd.DataFrame, method: str = "zscore", n_bins: int = 50
) -> Dict:
    cols = list(map(str, X_train.columns))
    A = X_train.to_numpy(dtype=np.float64, copy=True)
    method = method.lower()

    if method in {"none", "raw"}:
        return {"method": "none", "columns": cols}

    if method == "zscore":
        mean = np.nanmean(A, axis=0)
        std = np.nanstd(A, axis=0)
        std[~np.isfinite(std) | (std < EPS)] = 1.0
        return {"method": "zscore", "columns": cols,
                "mean": mean.astype(np.float32), "std": std.astype(np.float32)}

    if method in {"binning", "qbin"}:
        qs = np.linspace(0.0, 1.0, int(n_bins) + 1)[1:-1]
        edges = {}
        for j, c in enumerate(cols):
            col = A[:, j]
            col = col[np.isfinite(col)]
            edges[c] = (np.quantile(col, qs).astype(np.float32)
                        if col.size else np.array([], np.float32))
        return {"method": "qbin", "columns": cols, "n_bins": int(n_bins), "edges": edges}

    if method in {"ranking", "rank"}:
        sorted_values = {}
        for j, c in enumerate(cols):
            col = A[:, j]
            col = np.sort(col[np.isfinite(col)])
            sorted_values[c] = col.astype(np.float32)
        return {"method": "rank", "columns": cols, "sorted_values": sorted_values}

    raise ValueError(f"Unknown normalization method: {method}")


def apply_feature_normalizer(X: pd.DataFrame, normalizer: Dict) -> np.ndarray:
    """Apply a fitted feature normalizer, returning a float32 matrix in column order."""
    method = normalizer.get("method", "none")
    cols = list(map(str, normalizer["columns"]))
    Xr = X.copy()
    Xr.columns = Xr.columns.map(str)
    Xr = Xr.reindex(columns=cols)
    A = Xr.to_numpy(dtype=np.float32, copy=True)

    if method == "none":
        A[~np.isfinite(A)] = 0.0
        return A

    if method == "zscore":
        mean = np.asarray(normalizer["mean"], np.float32)
        std = np.asarray(normalizer["std"], np.float32)
        std[~np.isfinite(std) | (std == 0)] = 1.0
        for j in range(A.shape[1]):
            bad = ~np.isfinite(A[:, j])
            if bad.any():
                A[bad, j] = mean[j]  # -> maps to 0 after centering
        return ((A - mean) / std).astype(np.float32)

    if method == "qbin":
        n_bins = int(normalizer.get("n_bins", 50))
        denom = max(1, n_bins - 1)
        out = np.zeros_like(A, dtype=np.float32)
        for j, c in enumerate(cols):
            edges = np.asarray(normalizer["edges"].get(c, []), np.float32)
            vals = A[:, j]
            m = np.isfinite(vals)
            if edges.size:
                out[m, j] = np.searchsorted(edges, vals[m], side="right") / float(denom)
        return out

    if method == "rank":
        out = np.zeros_like(A, dtype=np.float32)
        for j, c in enumerate(cols):
            sv = np.asarray(normalizer["sorted_values"].get(c, []), np.float32)
            vals = A[:, j]
            m = np.isfinite(vals)
            if sv.size:
                out[m, j] = np.searchsorted(sv, vals[m], side="right") / float(sv.size)
        return out

    raise ValueError(method)
