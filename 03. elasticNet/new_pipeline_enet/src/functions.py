#!/usr/bin/env python3
"""ElasticNet 80/20 train-test pipeline utilities.

Designed for project layout:
  new_pipeline_enet/
    src/functions.py
    src/pipeline.py
    data/transcriptomics_filtered.joblib
    data/proteomics_matched.joblib
    results_enet/

This file intentionally implements a faster 80/20 workflow rather than the old
LOO ElasticNet workflow. It saves compact H5 files containing internal test
metrics and a final all-samples model artifact for external prediction.
"""

import itertools
import json
import lzma
import os
import pickle
import re
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from sklearn.feature_selection import f_regression
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import KFold, train_test_split

RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)


# -----------------------------------------------------------------------------
# General helpers
# -----------------------------------------------------------------------------

def safe_path_component(value: str) -> str:
    out = str(value).replace(os.sep, "_").replace("/", "_").replace("\\", "_")
    out = re.sub(r"[^A-Za-z0-9._+=@-]+", "_", out)
    return out


def all_param_combinations(param_grid: Dict[str, Sequence[Any]]) -> List[Dict[str, Any]]:
    if not param_grid:
        return [dict()]
    keys, vals = zip(*param_grid.items())
    return [dict(zip(keys, combo)) for combo in itertools.product(*vals)]


def sample_param_combos(
    param_grid: Dict[str, Sequence[Any]],
    max_param_combos: Optional[int] = None,
) -> List[Dict[str, Any]]:
    combos = list(all_param_combinations(param_grid))
    if max_param_combos is not None and len(combos) > int(max_param_combos):
        rng = np.random.default_rng(RANDOM_STATE)
        idx = rng.choice(len(combos), size=int(max_param_combos), replace=False)
        combos = [combos[i] for i in idx]
    return combos


def validate_normalization(normalization: str) -> str:
    n = (normalization or "none").lower().strip()
    aliases = {
        "rank": "ranking",
        "ranks": "ranking",
        "qbin": "binning",
        "qbin50": "binning",
        "bin": "binning",
        "bins": "binning",
    }
    n = aliases.get(n, n)
    if n not in {"none", "zscore", "ranking", "binning"}:
        raise ValueError("normalization must be one of: none, zscore, ranking, binning")
    return n


def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, Any]:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    out = {
        "test_mse": np.nan,
        "test_pearson_corr": np.nan,
        "test_cosine_similarity": np.nan,
        "n_test_samples": int(len(y_true)),
        "n_valid_test_samples": int(mask.sum()),
    }
    if mask.sum() == 0:
        return out
    yt = y_true[mask]
    yp = y_pred[mask]
    out["test_mse"] = float(np.mean((yt - yp) ** 2))
    if len(yt) >= 2 and np.nanstd(yt) > 0 and np.nanstd(yp) > 0:
        out["test_pearson_corr"] = float(np.corrcoef(yt, yp)[0, 1])
    denom = float(np.linalg.norm(yt) * np.linalg.norm(yp))
    if denom > 0:
        out["test_cosine_similarity"] = float(np.dot(yt, yp) / denom)
    return out


def json_safe(obj: Any) -> Any:
    if obj is None or isinstance(obj, (str, bool, int)):
        return obj
    if isinstance(obj, float):
        return None if not np.isfinite(obj) else obj
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        val = float(obj)
        return None if not np.isfinite(val) else val
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, pd.Timestamp):
        return obj.isoformat()
    if isinstance(obj, dict):
        return {str(k): json_safe(v) for k, v in obj.items() if k != "model"}
    if isinstance(obj, (list, tuple, set)):
        return [json_safe(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return json_safe(obj.tolist())
    try:
        if pd.isna(obj):
            return None
    except Exception:
        pass
    return str(obj)


def fixed_utf8_array(values: Sequence[Any]) -> np.ndarray:
    encoded = [str(v).encode("utf-8") for v in values]
    width = max([len(v) for v in encoded] + [1])
    arr = np.empty(len(encoded), dtype=f"S{width}")
    arr[:] = encoded
    return arr


def write_dataset_compressed(group, name: str, data, *, compression_level: int = 9):
    compression_level = max(0, min(9, int(compression_level)))
    if compression_level == 0:
        return group.create_dataset(name, data=data)
    return group.create_dataset(
        name,
        data=data,
        compression="gzip",
        compression_opts=compression_level,
        shuffle=True,
        fletcher32=True,
    )


def pickle_xz_bytes(obj: Any) -> np.ndarray:
    payload = pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL)
    compressed = lzma.compress(payload, preset=9 | lzma.PRESET_EXTREME)
    return np.frombuffer(compressed, dtype=np.uint8)


def unpickle_xz_bytes(payload: np.ndarray) -> Any:
    return pickle.loads(lzma.decompress(np.asarray(payload, dtype=np.uint8).tobytes()))


# -----------------------------------------------------------------------------
# Data orientation and alignment
# -----------------------------------------------------------------------------

def detect_protein_axis(X: pd.DataFrame, Y: pd.DataFrame) -> Tuple[List[str], str, int, int]:
    x_samples = pd.Index(map(str, X.index))
    y_index = pd.Index(map(str, Y.index))
    y_columns = pd.Index(map(str, Y.columns))

    index_overlap = len(x_samples.intersection(y_index))
    column_overlap = len(x_samples.intersection(y_columns))

    if index_overlap >= column_overlap:
        proteins = list(map(str, Y.columns))
        orientation = "samples_as_rows__proteins_as_columns"
    else:
        proteins = list(map(str, Y.index))
        orientation = "proteins_as_rows__samples_as_columns"

    return proteins, orientation, index_overlap, column_overlap


def extract_protein_series(Y: pd.DataFrame, protein_name: str) -> Tuple[pd.Series, str]:
    if protein_name in Y.columns:
        s = Y[protein_name].copy()
        s.name = protein_name
        return s, "samples_as_rows__proteins_as_columns"
    if protein_name in Y.index:
        s = Y.loc[protein_name].copy()
        s.name = protein_name
        return s, "proteins_as_rows__samples_as_columns"
    raise ValueError(f"Protein {protein_name!r} not found in proteomics rows or columns.")


def align_X_y_for_protein(
    X: pd.DataFrame,
    Y: pd.DataFrame,
    protein_name: str,
    *,
    n_samples: Optional[int] = None,
) -> Tuple[pd.DataFrame, pd.Series, str]:
    y_protein, orientation = extract_protein_series(Y, protein_name)
    common = X.index.intersection(y_protein.index)

    if len(common) == 0:
        X2 = X.copy()
        y2 = y_protein.copy()
        X2.index = X2.index.map(str)
        y2.index = y2.index.map(str)
        common = X2.index.intersection(y2.index)
        if len(common) == 0:
            raise ValueError("No overlapping sample IDs between transcriptomics and proteomics.")
        X = X2
        y_protein = y2

    if n_samples is not None:
        common = common[: int(n_samples)]

    X_aligned = X.loc[common]
    y_aligned = y_protein.loc[common]
    valid = y_aligned.notna()
    X_aligned = X_aligned.loc[valid]
    y_aligned = y_aligned.loc[valid]
    return X_aligned, y_aligned, orientation


# -----------------------------------------------------------------------------
# Normalization
# -----------------------------------------------------------------------------

def fit_zscore_normalizer(X_train: pd.DataFrame) -> Dict[str, Any]:
    arr = X_train.to_numpy(dtype=np.float32, copy=False)
    mean = np.nanmean(arr, axis=0).astype(np.float32)
    std = np.nanstd(arr, axis=0).astype(np.float32)
    mean[~np.isfinite(mean)] = 0.0
    std[~np.isfinite(std) | (std == 0)] = 1.0
    return {"method": "zscore", "columns": list(map(str, X_train.columns)), "mean": mean, "std": std}


def apply_zscore(X: pd.DataFrame, normalizer: Dict[str, Any]) -> pd.DataFrame:
    cols = list(map(str, normalizer["columns"]))
    X_aligned = X.reindex(columns=cols)
    arr = X_aligned.to_numpy(dtype=np.float32, copy=False)
    out = ((arr - normalizer["mean"]) / normalizer["std"]).astype(np.float32, copy=False)
    return pd.DataFrame(out, index=X_aligned.index, columns=cols)


def apply_ranking_rows(X: pd.DataFrame) -> pd.DataFrame:
    arr = X.to_numpy(dtype=np.float32, copy=False)
    ranked = np.argsort(np.argsort(arr, axis=1), axis=1).astype(np.float32, copy=False)
    denom = max(1, arr.shape[1] - 1)
    ranked = ranked / float(denom)
    return pd.DataFrame(ranked, index=X.index, columns=X.columns)


def apply_binning_rows(X: pd.DataFrame, n_bins: int = 50) -> pd.DataFrame:
    if int(n_bins) < 2:
        raise ValueError("n_bins must be >= 2")
    x = X.to_numpy(dtype=np.float64, copy=True)
    x[x == 0] = np.nan
    qs = np.linspace(0.0, 1.0, int(n_bins), dtype=np.float64)
    q = np.nanquantile(x, qs, axis=1).T
    out = np.zeros_like(x, dtype=np.float32)
    for i in range(x.shape[0]):
        t = q[i]
        if np.isnan(t).all():
            continue
        row = x[i]
        m = ~np.isnan(row)
        if not m.any():
            continue
        idx = np.searchsorted(t, row[m], side="left")
        out[i, m] = (float(n_bins) - idx.astype(np.float32)) / float(n_bins)
    return pd.DataFrame(out, index=X.index, columns=X.columns)


def fit_normalizer(X_train: pd.DataFrame, normalization: str, *, n_bins: int = 50) -> Dict[str, Any]:
    normalization = validate_normalization(normalization)
    if normalization == "zscore":
        return fit_zscore_normalizer(X_train)
    if normalization == "ranking":
        return {"method": "ranking", "columns": list(map(str, X_train.columns))}
    if normalization == "binning":
        return {"method": "binning", "columns": list(map(str, X_train.columns)), "n_bins": int(n_bins)}
    return {"method": "none", "columns": list(map(str, X_train.columns))}


def apply_normalizer(X: pd.DataFrame, normalizer: Dict[str, Any]) -> pd.DataFrame:
    method = normalizer.get("method", "none")
    cols = list(map(str, normalizer.get("columns", X.columns)))
    X_aligned = X.reindex(columns=cols)
    if method == "none":
        return X_aligned.astype(np.float32, copy=False)
    if method == "zscore":
        return apply_zscore(X_aligned, normalizer)
    if method == "ranking":
        return apply_ranking_rows(X_aligned)
    if method == "binning":
        return apply_binning_rows(X_aligned, n_bins=int(normalizer.get("n_bins", 50)))
    raise ValueError(f"Unknown normalizer method: {method}")


# -----------------------------------------------------------------------------
# Feature selection and ElasticNet fitting
# -----------------------------------------------------------------------------

def select_top_k_features(X_train: pd.DataFrame, y_train: pd.Series, select_k: int) -> List[str]:
    if select_k is None or int(select_k) <= 0 or X_train.shape[1] <= int(select_k):
        return list(map(str, X_train.columns))

    y_arr = y_train.to_numpy(dtype=float, copy=False)
    X_arr = X_train.to_numpy(dtype=float, copy=False)
    scores, _ = f_regression(X_arr, y_arr)
    scores = np.asarray(scores, dtype=float)
    np.nan_to_num(scores, copy=False, nan=-np.inf, posinf=-np.inf, neginf=-np.inf)

    k = min(int(select_k), X_train.shape[1])
    if not np.isfinite(scores).any():
        return list(map(str, X_train.columns[:k]))

    idx = np.argpartition(scores, X_train.shape[1] - k)[X_train.shape[1] - k:]
    idx = idx[np.argsort(scores[idx])[::-1]]
    return list(map(str, X_train.columns[idx]))


def fit_enet(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    params: Dict[str, Any],
    *,
    enet_max_iter: int = 10000,
    enet_tol: float = 1e-3,
    enet_selection: str = "random",
    fit_intercept: bool = False,
) -> ElasticNet:
    est = ElasticNet(
        max_iter=int(enet_max_iter),
        tol=float(enet_tol),
        selection=str(enet_selection),
        random_state=RANDOM_STATE,
        fit_intercept=bool(fit_intercept),
        **params,
    )
    est.fit(X_train.to_numpy(dtype=np.float32, copy=False), y_train.to_numpy(dtype=float, copy=False))
    return est


def tune_enet_params(
    X_train_raw: pd.DataFrame,
    y_train: pd.Series,
    param_grid: Dict[str, Sequence[Any]],
    *,
    normalization: str,
    n_bins: int = 50,
    select_k: int = 2000,
    n_splits: int = 3,
    max_param_combos: Optional[int] = 8,
    enet_max_iter: int = 10000,
    enet_tol: float = 1e-3,
    enet_selection: str = "random",
    fit_intercept: bool = False,
) -> Dict[str, Any]:
    combos = sample_param_combos(param_grid, max_param_combos=max_param_combos)
    if not combos:
        return {}
    n_splits = min(int(n_splits), len(y_train))
    if n_splits < 2:
        return dict(combos[0])

    kf = KFold(n_splits=n_splits, shuffle=True, random_state=RANDOM_STATE)
    avg_mse = np.zeros(len(combos), dtype=float)

    for tr_idx, va_idx in kf.split(X_train_raw):
        X_fit_raw = X_train_raw.iloc[tr_idx]
        X_val_raw = X_train_raw.iloc[va_idx]
        y_fit = y_train.iloc[tr_idx]
        y_val = y_train.iloc[va_idx]

        norm = fit_normalizer(X_fit_raw, normalization, n_bins=n_bins)
        X_fit_norm = apply_normalizer(X_fit_raw, norm)
        X_val_norm = apply_normalizer(X_val_raw, norm)

        selected = select_top_k_features(X_fit_norm, y_fit, select_k)
        X_fit_sel = X_fit_norm[selected]
        X_val_sel = X_val_norm[selected]

        for i, params in enumerate(combos):
            try:
                est = fit_enet(
                    X_fit_sel,
                    y_fit,
                    params,
                    enet_max_iter=enet_max_iter,
                    enet_tol=enet_tol,
                    enet_selection=enet_selection,
                    fit_intercept=fit_intercept,
                )
                pred = est.predict(X_val_sel.to_numpy(dtype=np.float32, copy=False))
                avg_mse[i] += mean_squared_error(y_val.to_numpy(dtype=float, copy=False), pred) / n_splits
            except Exception as e:
                print(f"WARNING: ElasticNet CV failed for params={params}: {type(e).__name__}: {e}")
                avg_mse[i] += 1e6 / n_splits

    best_idx = int(np.argmin(avg_mse))
    best = dict(combos[best_idx])
    best["_cv_mse"] = float(avg_mse[best_idx])
    return best


# -----------------------------------------------------------------------------
# Core 80/20 pipeline
# -----------------------------------------------------------------------------

def run_enet_8020_for_protein(
    X: pd.DataFrame,
    Y: pd.DataFrame,
    protein_name: str,
    *,
    param_grid: Dict[str, Sequence[Any]],
    normalization: str = "zscore",
    n_bins: int = 50,
    select_k: int = 2000,
    test_size: float = 0.2,
    n_splits: int = 3,
    max_param_combos: Optional[int] = 8,
    enet_max_iter: int = 10000,
    enet_tol: float = 1e-3,
    enet_selection: str = "random",
    fit_intercept: bool = False,
    n_samples: Optional[int] = None,
    h5_path: Optional[str] = None,
    h5_compression_level: int = 9,
) -> Dict[str, Any]:
    normalization = validate_normalization(normalization)
    X_work, y_work, orientation = align_X_y_for_protein(X, Y, protein_name, n_samples=n_samples)
    if len(y_work) < max(10, int(n_splits) + 2):
        raise ValueError(f"Not enough samples for {protein_name}: n={len(y_work)}")

    train_idx, test_idx = train_test_split(
        np.arange(len(y_work)),
        test_size=float(test_size),
        random_state=RANDOM_STATE,
        shuffle=True,
    )
    train_samples = list(map(str, X_work.index[train_idx]))
    test_samples = list(map(str, X_work.index[test_idx]))

    X_train_raw = X_work.iloc[train_idx]
    X_test_raw = X_work.iloc[test_idx]
    y_train = y_work.iloc[train_idx]
    y_test = y_work.iloc[test_idx]

    best_params = tune_enet_params(
        X_train_raw,
        y_train,
        param_grid,
        normalization=normalization,
        n_bins=n_bins,
        select_k=select_k,
        n_splits=n_splits,
        max_param_combos=max_param_combos,
        enet_max_iter=enet_max_iter,
        enet_tol=enet_tol,
        enet_selection=enet_selection,
        fit_intercept=fit_intercept,
    )
    clean_best_params = {k: v for k, v in best_params.items() if not str(k).startswith("_")}

    # Internal 20% test model: trained only on train split.
    train_norm = fit_normalizer(X_train_raw, normalization, n_bins=n_bins)
    X_train_norm = apply_normalizer(X_train_raw, train_norm)
    X_test_norm = apply_normalizer(X_test_raw, train_norm)
    selected_train = select_top_k_features(X_train_norm, y_train, select_k)
    test_model = fit_enet(
        X_train_norm[selected_train],
        y_train,
        clean_best_params,
        enet_max_iter=enet_max_iter,
        enet_tol=enet_tol,
        enet_selection=enet_selection,
        fit_intercept=fit_intercept,
    )
    test_pred = test_model.predict(X_test_norm[selected_train].to_numpy(dtype=np.float32, copy=False))
    metrics = compute_metrics(y_test.to_numpy(dtype=float, copy=False), test_pred)

    # Final model: trained on all Dataset 1 samples for external prediction.
    final_norm = fit_normalizer(X_work, normalization, n_bins=n_bins)
    X_all_norm = apply_normalizer(X_work, final_norm)
    selected_final = select_top_k_features(X_all_norm, y_work, select_k)
    final_model = fit_enet(
        X_all_norm[selected_final],
        y_work,
        clean_best_params,
        enet_max_iter=enet_max_iter,
        enet_tol=enet_tol,
        enet_selection=enet_selection,
        fit_intercept=fit_intercept,
    )

    coef = final_model.coef_.astype(float, copy=False)
    nonzero = int(np.sum(np.abs(coef) > 0))

    artifact = {
        "artifact_kind": "final_all_samples_enet_model",
        "protein": str(protein_name),
        "model": final_model,
        "normalization": normalization,
        "normalizer": final_norm,
        "selected_features": list(map(str, selected_final)),
        "params": clean_best_params,
        "params_with_cv_metadata": dict(best_params),
        "fit_intercept": bool(fit_intercept),
        "intercept": float(final_model.intercept_) if fit_intercept else 0.0,
        "coef_nonzero_count": nonzero,
        "n_training_samples_final_model": int(len(y_work)),
        "training_samples_final_model": list(map(str, X_work.index)),
    }

    result = {
        "protein": str(protein_name),
        "proteomics_orientation": orientation,
        "train_samples": train_samples,
        "test_samples": test_samples,
        "test_true": y_test.to_numpy(dtype=float, copy=False),
        "test_predicted": np.asarray(test_pred, dtype=float),
        "metrics": metrics,
        "best_params": dict(best_params),
        "selected_features": list(map(str, selected_final)),
        "model_artifact": artifact,
        "normalization": normalization,
        "n_bins": int(n_bins),
        "select_k": int(select_k),
        "test_size": float(test_size),
    }

    if h5_path is not None:
        save_enet_8020_h5(result, h5_path, compression_level=h5_compression_level)

    return result


def predict_with_artifact(artifact: Dict[str, Any], X_new: pd.DataFrame) -> pd.Series:
    protein = artifact.get("protein", "prediction")
    features = list(map(str, artifact["selected_features"]))
    normalizer = artifact.get("normalizer", {"method": "none", "columns": features})
    X_norm = apply_normalizer(X_new, normalizer)
    X_sel = X_norm.reindex(columns=features)
    pred = artifact["model"].predict(X_sel.to_numpy(dtype=np.float32, copy=False))
    return pd.Series(pred, index=X_new.index, name=str(protein))


# -----------------------------------------------------------------------------
# H5 IO
# -----------------------------------------------------------------------------

def save_enet_8020_h5(result: Dict[str, Any], path: str, *, compression_level: int = 9) -> str:
    try:
        import h5py
    except ImportError as e:
        raise ImportError("Saving H5 output requires h5py to be installed.") from e

    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    metrics = dict(result["metrics"])
    artifact = dict(result["model_artifact"])
    protein = str(result["protein"])

    # Minimal metrics JSON intentionally excludes normalization, test_size, true/predicted and RMSE.
    metrics_json_obj = {
        "protein": protein,
        "train_samples": list(map(str, result["train_samples"])),
        "test_samples": list(map(str, result["test_samples"])),
        "test_mse": metrics.get("test_mse", np.nan),
        "test_pearson_corr": metrics.get("test_pearson_corr", np.nan),
        "test_cosine_similarity": metrics.get("test_cosine_similarity", np.nan),
        "best_params": result.get("best_params", {}),
        "selected_features": list(map(str, result.get("selected_features", []))),
    }
    metrics_json = json.dumps(json_safe(metrics_json_obj), ensure_ascii=False, sort_keys=True).encode("utf-8")
    params_json = json.dumps(json_safe(result.get("best_params", {})), ensure_ascii=False, sort_keys=True).encode("utf-8")

    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "minimal_protein_8020_elasticnet_h5"
        h5.attrs["format_version"] = "1.0"
        h5.attrs["protein"] = protein
        h5.attrs["model_label"] = "final_all_samples_enet_model"

        g = h5.create_group("metrics")
        write_dataset_compressed(g, "json", np.frombuffer(metrics_json, dtype=np.uint8), compression_level=compression_level)
        g.attrs["test_mse"] = np.nan if metrics.get("test_mse") is None else float(metrics.get("test_mse", np.nan))
        g.attrs["test_pearson_corr"] = np.nan if metrics.get("test_pearson_corr") is None else float(metrics.get("test_pearson_corr", np.nan))
        g.attrs["test_cosine_similarity"] = np.nan if metrics.get("test_cosine_similarity") is None else float(metrics.get("test_cosine_similarity", np.nan))

        g = h5.create_group("split")
        write_dataset_compressed(g, "train_samples", fixed_utf8_array(result["train_samples"]), compression_level=compression_level)
        write_dataset_compressed(g, "test_samples", fixed_utf8_array(result["test_samples"]), compression_level=compression_level)

        g = h5.create_group("model")
        g.attrs["artifact_kind"] = str(artifact.get("artifact_kind"))
        g.attrs["selected_feature_count"] = int(len(artifact.get("selected_features", [])))
        g.attrs["coef_nonzero_count"] = int(artifact.get("coef_nonzero_count", 0))
        write_dataset_compressed(g, "artifact_pickle_xz", pickle_xz_bytes(artifact), compression_level=compression_level)
        write_dataset_compressed(g, "selected_features", fixed_utf8_array(artifact.get("selected_features", [])), compression_level=compression_level)
        write_dataset_compressed(g, "params_json", np.frombuffer(params_json, dtype=np.uint8), compression_level=compression_level)

    return path


def load_enet_8020_h5(path: str) -> Dict[str, Any]:
    try:
        import h5py
    except ImportError as e:
        raise ImportError("Loading H5 output requires h5py to be installed.") from e

    with h5py.File(path, "r") as h5:
        metrics = json.loads(bytes(h5["metrics/json"][:]).decode("utf-8"))
        train_samples = [x.decode("utf-8") for x in h5["split/train_samples"][:]]
        test_samples = [x.decode("utf-8") for x in h5["split/test_samples"][:]]
        artifact = unpickle_xz_bytes(h5["model/artifact_pickle_xz"][:])
    return {"metrics": metrics, "train_samples": train_samples, "test_samples": test_samples, "model_artifact": artifact}
