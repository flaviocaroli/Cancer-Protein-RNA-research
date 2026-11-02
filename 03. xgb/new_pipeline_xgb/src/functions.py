import itertools
import json
import lzma
import os
import pickle
import re
import warnings
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import KFold, train_test_split
import xgboost as xgb

np.random.seed(42)


# -----------------------------
# Normalization helpers
# -----------------------------
def is_qbin_normalization(normalization: str) -> bool:
    normalization = (normalization or "").lower().strip()
    return re.fullmatch(r"qbin[0-9]+", normalization) is not None


def qbin_count(normalization: str) -> int:
    if not is_qbin_normalization(normalization):
        raise ValueError(f"Not a qbin normalization: {normalization}")
    n_bins = int(normalization.lower().replace("qbin", ""))
    if n_bins < 2:
        raise ValueError("qbin normalization requires at least 2 bins.")
    return n_bins


def validate_normalization(normalization: str) -> str:
    normalization = (normalization or "none").lower().strip()
    if normalization in {"none", "zscore", "rank"}:
        return normalization
    if is_qbin_normalization(normalization):
        qbin_count(normalization)
        return normalization
    raise ValueError(
        f"Invalid normalization: {normalization}. Allowed: none, zscore, rank, or qbinN such as qbin50."
    )


def _fit_none_normalizer(columns: Sequence[str]) -> Dict[str, Any]:
    return {"method": "none", "columns": list(map(str, columns))}


def _fit_zscore_normalizer(X_train: pd.DataFrame) -> Dict[str, Any]:
    arr = X_train.to_numpy(dtype=np.float32, copy=False)
    mean = np.nanmean(arr, axis=0).astype(np.float32)
    std = np.nanstd(arr, axis=0).astype(np.float32)
    mean[~np.isfinite(mean)] = 0.0
    std[~np.isfinite(std) | (std == 0)] = 1.0
    return {
        "method": "zscore",
        "columns": list(map(str, X_train.columns)),
        "mean": mean,
        "std": std,
    }


def _fit_rank_normalizer(X_train: pd.DataFrame) -> Dict[str, Any]:
    sorted_values = {}
    for col in X_train.columns:
        arr = X_train[col].to_numpy(dtype=np.float32, copy=False)
        arr = arr[np.isfinite(arr)]
        sorted_values[str(col)] = np.sort(arr.astype(np.float32, copy=False))
    return {
        "method": "rank",
        "columns": list(map(str, X_train.columns)),
        "sorted_values": sorted_values,
    }


def _fit_qbin_normalizer(X_train: pd.DataFrame, n_bins: int) -> Dict[str, Any]:
    if n_bins < 2:
        raise ValueError("qbin normalization requires at least 2 bins.")
    qs = np.linspace(0.0, 1.0, n_bins + 1)[1:-1]
    edges = {}
    for col in X_train.columns:
        arr = X_train[col].to_numpy(dtype=np.float32, copy=False)
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            edges[str(col)] = np.array([], dtype=np.float32)
        else:
            edges[str(col)] = np.unique(np.nanquantile(arr, qs).astype(np.float32))
    return {
        "method": f"qbin{int(n_bins)}",
        "n_bins": int(n_bins),
        "columns": list(map(str, X_train.columns)),
        "edges": edges,
    }


def fit_normalizer(X_train: pd.DataFrame, normalization: str) -> Dict[str, Any]:
    normalization = validate_normalization(normalization)
    if normalization == "none":
        return _fit_none_normalizer(X_train.columns)
    if normalization == "zscore":
        return _fit_zscore_normalizer(X_train)
    if normalization == "rank":
        return _fit_rank_normalizer(X_train)
    if is_qbin_normalization(normalization):
        return _fit_qbin_normalizer(X_train, qbin_count(normalization))
    raise ValueError(f"Unknown normalization: {normalization}")


def apply_normalizer_to_frame(X: pd.DataFrame, normalizer: Dict[str, Any]) -> pd.DataFrame:
    method = normalizer.get("method", "none")
    cols = list(map(str, normalizer.get("columns", X.columns)))
    X_aligned = X.reindex(columns=cols)

    if method == "none":
        return X_aligned.astype(np.float32, copy=False)

    if method == "zscore":
        arr = X_aligned.to_numpy(dtype=np.float32, copy=False)
        out = ((arr - normalizer["mean"]) / normalizer["std"]).astype(np.float32, copy=False)
        return pd.DataFrame(out, index=X_aligned.index, columns=cols)

    if method == "rank":
        out = np.full((X_aligned.shape[0], len(cols)), np.nan, dtype=np.float32)
        for j, col in enumerate(cols):
            sorted_col = normalizer["sorted_values"].get(col, np.array([], dtype=np.float32))
            if sorted_col.size == 0:
                continue
            vals = X_aligned.iloc[:, j].to_numpy(dtype=np.float32, copy=False)
            mask = np.isfinite(vals)
            out[mask, j] = np.searchsorted(sorted_col, vals[mask], side="right") / sorted_col.size
        return pd.DataFrame(out, index=X_aligned.index, columns=cols)

    if is_qbin_normalization(method):
        n_bins = int(normalizer["n_bins"])
        denom = max(1, n_bins - 1)
        out = np.full((X_aligned.shape[0], len(cols)), np.nan, dtype=np.float32)
        for j, col in enumerate(cols):
            edges = normalizer["edges"].get(col, np.array([], dtype=np.float32))
            vals = X_aligned.iloc[:, j].to_numpy(dtype=np.float32, copy=False)
            mask = np.isfinite(vals)
            if edges.size == 0:
                out[mask, j] = 0.0
            else:
                out[mask, j] = np.searchsorted(edges, vals[mask], side="right") / denom
        return pd.DataFrame(out, index=X_aligned.index, columns=cols)

    raise ValueError(f"Unknown normalizer method: {method}")


# -----------------------------
# XGBoost helpers
# -----------------------------
def all_param_combinations(param_grid: dict) -> Sequence[dict]:
    if not param_grid:
        return [dict()]
    keys, vals = zip(*param_grid.items())
    return [dict(zip(keys, combo)) for combo in itertools.product(*vals)]


def sample_param_combos(param_grid: dict, max_param_combos: Optional[int]) -> Sequence[dict]:
    combos = list(all_param_combinations(param_grid))
    if max_param_combos is not None and len(combos) > int(max_param_combos):
        rng = np.random.default_rng(42)
        idx = rng.choice(len(combos), size=int(max_param_combos), replace=False)
        combos = [combos[i] for i in idx]
    return combos


def default_xgb_threads(xgb_threads: Optional[int]) -> int:
    if xgb_threads is not None and int(xgb_threads) > 0:
        return int(xgb_threads)
    slurm_threads = os.getenv("SLURM_CPUS_PER_TASK")
    if slurm_threads and slurm_threads.isdigit():
        return max(1, int(slurm_threads))
    return max(1, os.cpu_count() or 1)


def clean_params(params: Optional[dict]) -> dict:
    return {k: v for k, v in (params or {}).items() if not str(k).startswith("_")}


def clone_xgb_with_params(
    params: Optional[dict],
    *,
    seed: int = 42,
    use_gpu: bool = False,
    xgb_threads: Optional[int] = None,
) -> xgb.XGBRegressor:
    base = dict(
        objective="reg:squarederror",
        n_jobs=default_xgb_threads(xgb_threads),
        random_state=seed,
        verbosity=0,
        tree_method="hist",
        device="cuda" if use_gpu else "cpu",
    )
    base.update(clean_params(params))
    return xgb.XGBRegressor(**base)


def fit_xgb(
    params: Optional[dict],
    Xtr: np.ndarray,
    ytr: np.ndarray,
    *,
    Xva: Optional[np.ndarray] = None,
    yva: Optional[np.ndarray] = None,
    early_stopping_rounds: int = 0,
    use_gpu: bool = False,
    xgb_threads: Optional[int] = None,
) -> xgb.XGBRegressor:
    use_es = Xva is not None and yva is not None and int(early_stopping_rounds or 0) > 0
    if not use_es:
        est = clone_xgb_with_params(params, use_gpu=use_gpu, xgb_threads=xgb_threads)
        est.fit(Xtr, ytr)
        return est

    est = clone_xgb_with_params(params, use_gpu=use_gpu, xgb_threads=xgb_threads)
    try:
        est.fit(
            Xtr,
            ytr,
            eval_set=[(Xva, yva)],
            early_stopping_rounds=int(early_stopping_rounds),
            verbose=False,
        )
        return est
    except TypeError:
        params2 = dict(clean_params(params))
        params2["early_stopping_rounds"] = int(early_stopping_rounds)
        est = clone_xgb_with_params(params2, use_gpu=use_gpu, xgb_threads=xgb_threads)
        est.fit(Xtr, ytr, eval_set=[(Xva, yva)], verbose=False)
        return est


def best_iteration_count(est: xgb.XGBRegressor, fallback: int) -> int:
    best = getattr(est, "best_iteration", None)
    if best is None:
        try:
            best = est.get_booster().best_iteration
        except Exception:
            best = None
    return int(fallback) if best is None else max(1, int(best) + 1)


# -----------------------------
# Feature selection and metrics
# -----------------------------
def spearman_abs_corr_topk(X: pd.DataFrame, y: pd.Series, k: Optional[int]) -> pd.Index:
    if k is None or int(k) <= 0 or X.shape[1] <= int(k):
        return X.columns

    k = int(k)
    y_num = pd.to_numeric(y, errors="coerce")
    valid_rows = y_num.notna()
    X_use = X.loc[valid_rows]
    y_use = y_num.loc[valid_rows]

    if len(y_use) < 2 or y_use.nunique(dropna=True) < 2:
        return X.columns[: min(k, X.shape[1])]

    valid_cols = []
    for col in X_use.columns:
        s = pd.to_numeric(X_use[col], errors="coerce")
        if s.notna().sum() >= 2 and s.nunique(dropna=True) >= 2:
            valid_cols.append(col)

    if not valid_cols:
        return X.columns[: min(k, X.shape[1])]

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        scores = (
            X_use[valid_cols]
            .rank(axis=0, method="average")
            .corrwith(y_use.rank(method="average"))
            .abs()
            .replace([np.inf, -np.inf], np.nan)
            .fillna(0.0)
        )

    kept = list(scores.nlargest(k).index)
    if len(kept) < k:
        kept.extend([c for c in X.columns if c not in kept][: k - len(kept)])
    return pd.Index(kept)


def prediction_metrics(y_true: Sequence[float], y_pred: Sequence[float]) -> Dict[str, Any]:
    yt0 = np.asarray(y_true, dtype=float)
    yp0 = np.asarray(y_pred, dtype=float)
    mask = np.isfinite(yt0) & np.isfinite(yp0)
    out = {
        "n_test_predictions": int(len(yt0)),
        "n_valid_test_predictions": int(mask.sum()),
        "test_mse": np.nan,
        "test_pearson_corr": np.nan,
        "test_cosine_similarity": np.nan,
    }
    if mask.sum() == 0:
        return out
    yt = yt0[mask]
    yp = yp0[mask]
    err = yt - yp
    out["test_mse"] = float(np.mean(err ** 2))
    if len(yt) >= 2 and np.nanstd(yt) > 0 and np.nanstd(yp) > 0:
        out["test_pearson_corr"] = float(np.corrcoef(yt, yp)[0, 1])
    denom = float(np.linalg.norm(yt) * np.linalg.norm(yp))
    if denom > 0:
        out["test_cosine_similarity"] = float(np.dot(yt, yp) / denom)
    return out


# -----------------------------
# Training pipeline
# -----------------------------
def fit_and_score_cv(
    param_grid: dict,
    X_train_sel: pd.DataFrame,
    y_train: pd.Series,
    n_splits: int,
    *,
    early_stopping_rounds: int = 10,
    max_param_combos: Optional[int] = None,
    use_gpu: bool = False,
    xgb_threads: Optional[int] = None,
    use_best_iteration: bool = True,
) -> dict:
    params_list = list(sample_param_combos(param_grid, max_param_combos))
    if not params_list:
        return {}

    n_splits = min(int(n_splits), len(y_train))
    if n_splits < 2:
        return dict(params_list[0])

    avg_mse = np.zeros(len(params_list), dtype=float)
    best_rounds: List[List[int]] = [[] for _ in params_list]
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

    X_values = X_train_sel.to_numpy(dtype=np.float32, copy=False)
    y_values = y_train.to_numpy(dtype=np.float32, copy=False)

    for tr_idx, va_idx in kf.split(X_values):
        Xtr = X_values[tr_idx]
        Xva = X_values[va_idx]
        ytr = y_values[tr_idx]
        yva = y_values[va_idx]
        for j, params in enumerate(params_list):
            try:
                est = fit_xgb(
                    params,
                    Xtr,
                    ytr,
                    Xva=Xva,
                    yva=yva,
                    early_stopping_rounds=early_stopping_rounds,
                    use_gpu=use_gpu,
                    xgb_threads=xgb_threads,
                )
                pred = est.predict(Xva)
                avg_mse[j] += mean_squared_error(yva, pred) / n_splits
                if use_best_iteration and int(early_stopping_rounds or 0) > 0:
                    best_rounds[j].append(best_iteration_count(est, int(params.get("n_estimators", 100))))
            except Exception as e:
                print(f"WARNING: CV fit failed for params={params}: {type(e).__name__}: {e}")
                avg_mse[j] += 1e6 / n_splits

    best_idx = int(np.argmin(avg_mse))
    best_params = dict(params_list[best_idx])
    if use_best_iteration and best_rounds[best_idx]:
        original = int(best_params.get("n_estimators", 100))
        median_best = int(np.median(best_rounds[best_idx]))
        best_params["n_estimators"] = max(25, min(original, int(np.ceil(median_best * 1.10))))
        best_params["_cv_median_best_iteration"] = median_best
    best_params["_cv_mse"] = float(avg_mse[best_idx])
    return best_params


def _extract_gain_importance(est: xgb.XGBRegressor, columns: Sequence[str]) -> Optional[Dict[str, float]]:
    try:
        gain = est.get_booster().get_score(importance_type="gain")
        fmap = {f"f{i}": str(col) for i, col in enumerate(columns)}
        return {fmap.get(k, k): float(v) for k, v in gain.items()}
    except Exception:
        return None


def extract_protein_series(y: pd.DataFrame, protein_name: str) -> Tuple[pd.Series, str]:
    if protein_name in y.columns:
        s = y[protein_name].copy()
        s.name = protein_name
        return s, "samples_as_rows__proteins_as_columns"
    if protein_name in y.index:
        s = y.loc[protein_name].copy()
        s.name = protein_name
        return s, "proteins_as_rows__samples_as_columns"
    raise ValueError(f"Protein {protein_name!r} not found in proteomics rows or columns.")


def align_X_and_y_for_protein(
    X: pd.DataFrame,
    y_protein: pd.Series,
    *,
    n_samples: Optional[int] = None,
) -> Tuple[pd.DataFrame, pd.Series]:
    common = X.index.intersection(y_protein.index)
    if len(common) == 0:
        X2 = X.copy()
        y2 = y_protein.copy()
        X2.index = X2.index.map(str)
        y2.index = y2.index.map(str)
        common = X2.index.intersection(y2.index)
        if len(common) == 0:
            raise ValueError("No overlapping sample IDs between transcriptomics rows and proteomics samples.")
        X = X2
        y_protein = y2
    if n_samples is not None:
        common = common[: int(n_samples)]
    return X.loc[common], y_protein.loc[common]


def train_test_pipeline_for_protein(
    *,
    X: pd.DataFrame,
    y: pd.DataFrame,
    protein_name: str,
    param_grid: dict,
    normalization: str = "none",
    test_fraction: float = 0.2,
    random_state: int = 42,
    n_splits: int = 3,
    select_k: int = 2000,
    early_stopping_rounds: int = 10,
    max_param_combos: Optional[int] = 8,
    n_samples: Optional[int] = None,
    use_gpu: bool = False,
    xgb_threads: Optional[int] = None,
) -> Dict[str, Any]:
    normalization = validate_normalization(normalization)
    y_protein, proteomics_orientation = extract_protein_series(y, protein_name)
    X_raw, y_aligned = align_X_and_y_for_protein(X, y_protein, n_samples=n_samples)

    valid = y_aligned.notna()
    X_raw = X_raw.loc[valid]
    y_aligned = y_aligned.loc[valid]

    if len(y_aligned) < max(10, int(n_splits) + 3):
        raise ValueError(f"Not enough non-missing samples for protein {protein_name}: n={len(y_aligned)}")

    train_idx, test_idx = train_test_split(
        np.arange(len(y_aligned)),
        test_size=float(test_fraction),
        random_state=int(random_state),
        shuffle=True,
    )

    X_train_raw = X_raw.iloc[train_idx]
    X_test_raw = X_raw.iloc[test_idx]
    y_train = y_aligned.iloc[train_idx]
    y_test = y_aligned.iloc[test_idx]

    eval_normalizer = fit_normalizer(X_train_raw, normalization)
    X_train_norm = apply_normalizer_to_frame(X_train_raw, eval_normalizer)
    X_test_norm = apply_normalizer_to_frame(X_test_raw, eval_normalizer)

    kept_eval = spearman_abs_corr_topk(X_train_norm, y_train, int(select_k))
    X_train_sel = X_train_norm[kept_eval]
    X_test_sel = X_test_norm[kept_eval]

    best_params = fit_and_score_cv(
        param_grid,
        X_train_sel,
        y_train,
        int(n_splits),
        early_stopping_rounds=int(early_stopping_rounds),
        max_param_combos=max_param_combos,
        use_gpu=use_gpu,
        xgb_threads=xgb_threads,
        use_best_iteration=True,
    )

    eval_model = fit_xgb(
        clean_params(best_params),
        X_train_sel.to_numpy(dtype=np.float32, copy=False),
        y_train.to_numpy(dtype=np.float32, copy=False),
        use_gpu=use_gpu,
        xgb_threads=xgb_threads,
    )
    y_pred = eval_model.predict(X_test_sel.to_numpy(dtype=np.float32, copy=False))
    metrics = prediction_metrics(y_test.to_numpy(dtype=float), y_pred)

    # Final model trained on all Dataset 1 samples using parameters selected from the train split.
    final_normalizer = fit_normalizer(X_raw, normalization)
    X_all_norm = apply_normalizer_to_frame(X_raw, final_normalizer)
    kept_final = spearman_abs_corr_topk(X_all_norm, y_aligned, int(select_k))
    X_all_sel = X_all_norm[kept_final]

    final_model = fit_xgb(
        clean_params(best_params),
        X_all_sel.to_numpy(dtype=np.float32, copy=False),
        y_aligned.to_numpy(dtype=np.float32, copy=False),
        use_gpu=use_gpu,
        xgb_threads=xgb_threads,
    )

    artifact = {
        "artifact_kind": "final_all_samples_model_after_80_20_evaluation",
        "protein": protein_name,
        "model": final_model,
        "params": clean_params(best_params),
        "params_with_cv_metadata": dict(best_params or {}),
        "selected_features": list(map(str, kept_final)),
        "normalizer": final_normalizer,
        "feature_selector": "spearman_prefilter",
        "n_training_samples": int(len(y_aligned)),
        "training_samples": list(map(str, X_raw.index)),
        "xgboost_version": xgb.__version__,
        "use_gpu": bool(use_gpu),
        "xgb_threads": default_xgb_threads(xgb_threads),
        "feature_importance_gain": _extract_gain_importance(final_model, kept_final),
    }

    summary = {
        "protein": str(protein_name),
        "train_samples": list(map(str, X_train_raw.index)),
        "test_samples": list(map(str, X_test_raw.index)),
        "test_mse": metrics["test_mse"],
        "test_pearson_corr": metrics["test_pearson_corr"],
        "test_cosine_similarity": metrics["test_cosine_similarity"],
        "n_test_predictions": metrics["n_test_predictions"],
        "n_valid_test_predictions": metrics["n_valid_test_predictions"],
        "best_params": clean_params(best_params),
        "selected_features": list(map(str, kept_final)),
        "proteomics_orientation": proteomics_orientation,
        "model_note": "Final model is trained on all Dataset 1 samples after 80/20 test metrics were computed.",
    }

    return {
        "summary": summary,
        "model_artifact": artifact,
    }


# -----------------------------
# Minimal H5 save/load
# -----------------------------
def _json_safe(obj: Any) -> Any:
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
        return {str(k): _json_safe(v) for k, v in obj.items() if k != "model"}
    if isinstance(obj, (list, tuple, set)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _json_safe(obj.tolist())
    try:
        if pd.isna(obj):
            return None
    except Exception:
        pass
    return str(obj)


def _fixed_utf8_array(values: Sequence[Any]) -> np.ndarray:
    encoded = [str(v).encode("utf-8") for v in values]
    width = max([len(v) for v in encoded] + [1])
    arr = np.empty(len(encoded), dtype=f"S{width}")
    arr[:] = encoded
    return arr


def _write_dataset_compressed(group, name: str, data, *, compression_level: int = 9):
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


def _pickle_xz_bytes(obj: Any) -> np.ndarray:
    payload = pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL)
    compressed = lzma.compress(payload, preset=9 | lzma.PRESET_EXTREME)
    return np.frombuffer(compressed, dtype=np.uint8)


def _unpickle_xz_bytes(payload: np.ndarray) -> Any:
    return pickle.loads(lzma.decompress(np.asarray(payload, dtype=np.uint8).tobytes()))


def save_minimal_8020_h5(
    result: Dict[str, Any],
    path: str,
    *,
    compression_level: int = 9,
) -> str:
    try:
        import h5py
    except ImportError as e:
        raise ImportError("H5 output requires h5py to be installed.") from e

    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    summary = result["summary"]
    model_artifact = result["model_artifact"]

    # Deliberately omitted from the H5 summary as requested:
    # normalization, test_size, test_true, test_predicted, rmse.
    summary_for_json = {
        "protein": summary["protein"],
        "train_samples": summary["train_samples"],
        "test_samples": summary["test_samples"],
        "test_mse": summary["test_mse"],
        "test_pearson_corr": summary["test_pearson_corr"],
        "test_cosine_similarity": summary["test_cosine_similarity"],
        "best_params": summary["best_params"],
        "selected_features": summary["selected_features"],
        "model_note": summary["model_note"],
    }
    summary_json = json.dumps(_json_safe(summary_for_json), ensure_ascii=False, sort_keys=True).encode("utf-8")
    params_json = json.dumps(_json_safe(summary["best_params"]), ensure_ascii=False, sort_keys=True).encode("utf-8")

    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "minimal_8020_xgboost_h5"
        h5.attrs["format_version"] = "1.0"
        h5.attrs["protein"] = str(summary["protein"])
        h5.attrs["model_kind"] = "final_all_samples_model"
        h5.attrs["compression"] = f"hdf5_gzip_{int(compression_level)}_shuffle_fletcher32_plus_model_pickle_xz_extreme"

        g = h5.create_group("metrics")
        g.attrs["test_mse"] = np.nan if summary["test_mse"] is None else float(summary["test_mse"])
        g.attrs["test_pearson_corr"] = np.nan if summary["test_pearson_corr"] is None else float(summary["test_pearson_corr"])
        g.attrs["test_cosine_similarity"] = np.nan if summary["test_cosine_similarity"] is None else float(summary["test_cosine_similarity"])
        _write_dataset_compressed(g, "summary_json", np.frombuffer(summary_json, dtype=np.uint8), compression_level=compression_level)

        g = h5.create_group("samples")
        _write_dataset_compressed(g, "train_samples", _fixed_utf8_array(summary["train_samples"]), compression_level=compression_level)
        _write_dataset_compressed(g, "test_samples", _fixed_utf8_array(summary["test_samples"]), compression_level=compression_level)

        g = h5.create_group("model")
        g.attrs["artifact_kind"] = str(model_artifact.get("artifact_kind"))
        g.attrs["selected_feature_count"] = int(len(summary["selected_features"]))
        _write_dataset_compressed(g, "artifact_pickle_xz", _pickle_xz_bytes(model_artifact), compression_level=compression_level)
        _write_dataset_compressed(g, "selected_features", _fixed_utf8_array(summary["selected_features"]), compression_level=compression_level)
        _write_dataset_compressed(g, "best_params_json", np.frombuffer(params_json, dtype=np.uint8), compression_level=compression_level)

    return path


def load_minimal_8020_h5(path: str) -> Dict[str, Any]:
    try:
        import h5py
    except ImportError as e:
        raise ImportError("Loading H5 output requires h5py to be installed.") from e

    with h5py.File(path, "r") as h5:
        summary = json.loads(bytes(h5["metrics/summary_json"][:]).decode("utf-8"))
        model_artifact = _unpickle_xz_bytes(h5["model/artifact_pickle_xz"][:])
    return {"summary": summary, "model_artifact": model_artifact}


def predict_with_artifact(artifact: Dict[str, Any], X_new: pd.DataFrame) -> pd.Series:
    protein = artifact.get("protein", "prediction")
    features = list(map(str, artifact["selected_features"]))
    X_sel = X_new.reindex(columns=features)
    X_trans = apply_normalizer_to_frame(X_sel, artifact["normalizer"])
    pred = artifact["model"].predict(X_trans.to_numpy(dtype=np.float32, copy=False))
    return pd.Series(pred, index=X_new.index, name=str(protein))
