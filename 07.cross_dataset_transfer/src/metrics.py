"""metrics.py -- per-protein metrics, matching 05.validation naming."""
from __future__ import annotations

import numpy as np
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error


def _cosine(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    na = np.linalg.norm(a); nb = np.linalg.norm(b)
    if na < 1e-12 or nb < 1e-12:
        return np.nan
    return float(np.dot(a, b) / (na * nb))


def protein_metrics(y_true, y_pred) -> dict:
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    m = np.isfinite(y_true) & np.isfinite(y_pred)
    y_true, y_pred = y_true[m], y_pred[m]
    out = {
        "n_test_used": int(y_true.size),
        "test_mse": np.nan,
        "test_pearson_corr": np.nan,
        "test_cosine_similarity": np.nan,
    }
    if y_true.size >= 2:
        out["test_mse"] = float(mean_squared_error(y_true, y_pred))
        out["test_cosine_similarity"] = _cosine(y_true, y_pred)
        if np.std(y_true) > 1e-12 and np.std(y_pred) > 1e-12 and y_true.size >= 3:
            out["test_pearson_corr"] = float(pearsonr(y_true, y_pred)[0])
    return out
