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
from sklearn.model_selection import KFold
import xgboost as xgb

np.random.seed(42)


def _is_qbin_normalization(normalization: str) -> bool:
    normalization = (normalization or "").lower().strip()
    return re.fullmatch(r"qbin[0-9]+", normalization) is not None


def _qbin_count(normalization: str) -> int:
    if not _is_qbin_normalization(normalization):
        raise ValueError(f"Not a qbin normalization: {normalization}")
    n_bins = int(normalization.lower().replace("qbin", ""))
    if n_bins < 2:
        raise ValueError("qbin normalization requires at least 2 bins.")
    return n_bins


def _all_param_combinations(param_grid: dict) -> Sequence[dict]:
    if not param_grid:
        return [dict()]
    keys, vals = zip(*param_grid.items())
    return [dict(zip(keys, combo)) for combo in itertools.product(*vals)]


def _sample_param_combos(param_grid: dict, max_param_combos: Optional[int]) -> Sequence[dict]:
    combos = list(_all_param_combinations(param_grid))
    if max_param_combos is not None and len(combos) > int(max_param_combos):
        rng = np.random.default_rng(42)
        idx = rng.choice(len(combos), size=int(max_param_combos), replace=False)
        combos = [combos[i] for i in idx]
    return combos


def _default_xgb_threads(xgb_threads: Optional[int]) -> int:
    if xgb_threads is not None and int(xgb_threads) > 0:
        return int(xgb_threads)

    slurm_threads = os.getenv("SLURM_CPUS_PER_TASK")
    if slurm_threads and slurm_threads.isdigit():
        return max(1, int(slurm_threads))

    return max(1, os.cpu_count() or 1)


def _clean_params(params: Optional[dict]) -> dict:
    return {k: v for k, v in (params or {}).items() if not str(k).startswith("_")}


def _clone_xgb_with_params(
    params: Optional[dict],
    *,
    seed: int = 42,
    use_gpu: bool = False,
    xgb_threads: Optional[int] = None,
) -> xgb.XGBRegressor:
    base = dict(
        objective="reg:squarederror",
        n_jobs=_default_xgb_threads(xgb_threads),
        random_state=seed,
        verbosity=0,
        tree_method="hist",
        device="cuda" if use_gpu else "cpu",
    )
    base.update(_clean_params(params))
    return xgb.XGBRegressor(**base)


def _fit_xgb(
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
        est = _clone_xgb_with_params(params, use_gpu=use_gpu, xgb_threads=xgb_threads)
        est.fit(Xtr, ytr)
        return est

    est = _clone_xgb_with_params(params, use_gpu=use_gpu, xgb_threads=xgb_threads)

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
        params2 = dict(_clean_params(params))
        params2["early_stopping_rounds"] = int(early_stopping_rounds)
        est = _clone_xgb_with_params(params2, use_gpu=use_gpu, xgb_threads=xgb_threads)
        est.fit(Xtr, ytr, eval_set=[(Xva, yva)], verbose=False)
        return est


def _best_iteration_count(est: xgb.XGBRegressor, fallback: int) -> int:
    best = getattr(est, "best_iteration", None)

    if best is None:
        try:
            best = est.get_booster().best_iteration
        except Exception:
            best = None

    return int(fallback) if best is None else max(1, int(best) + 1)


def _to_f32_series(s: pd.Series) -> np.ndarray:
    return s.to_numpy(dtype=np.float32, copy=False)


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


def _fit_rank_normalizer(X: pd.DataFrame) -> Dict[str, Any]:
    sorted_values = {}

    for col in X.columns:
        arr = X[col].to_numpy(dtype=np.float32, copy=False)
        arr = arr[np.isfinite(arr)]
        sorted_values[str(col)] = np.sort(arr.astype(np.float32, copy=False))

    return {
        "method": "rank",
        "columns": list(map(str, X.columns)),
        "sorted_values": sorted_values,
    }


def _fit_qbin_normalizer(X: pd.DataFrame, n_bins: int) -> Dict[str, Any]:
    if n_bins < 2:
        raise ValueError("qbin normalization requires at least 2 bins.")

    qs = np.linspace(0.0, 1.0, n_bins + 1)[1:-1]
    edges = {}

    for col in X.columns:
        arr = X[col].to_numpy(dtype=np.float32, copy=False)
        arr = arr[np.isfinite(arr)]

        if arr.size == 0:
            edges[str(col)] = np.array([], dtype=np.float32)
        else:
            edges[str(col)] = np.unique(np.nanquantile(arr, qs).astype(np.float32))

    return {
        "method": f"qbin{int(n_bins)}",
        "n_bins": int(n_bins),
        "columns": list(map(str, X.columns)),
        "edges": edges,
    }


def _fit_global_normalizer(X: pd.DataFrame, normalization: str) -> Dict[str, Any]:
    normalization = (normalization or "none").lower().strip()

    if normalization in {"none", "zscore"}:
        return _fit_none_normalizer(X.columns)

    if normalization == "rank":
        return _fit_rank_normalizer(X)

    if _is_qbin_normalization(normalization):
        return _fit_qbin_normalizer(X, _qbin_count(normalization))

    raise ValueError(f"Unknown normalization: {normalization}")


def _slice_normalizer(normalizer: Dict[str, Any], columns: Sequence[str]) -> Dict[str, Any]:
    cols = list(map(str, columns))
    method = normalizer.get("method", "none")

    if method == "none":
        return _fit_none_normalizer(cols)

    if method == "zscore":
        old_cols = list(map(str, normalizer["columns"]))
        idx = [old_cols.index(c) for c in cols]
        return {
            "method": "zscore",
            "columns": cols,
            "mean": np.asarray(normalizer["mean"], dtype=np.float32)[idx],
            "std": np.asarray(normalizer["std"], dtype=np.float32)[idx],
        }

    if method == "rank":
        return {
            "method": "rank",
            "columns": cols,
            "sorted_values": {
                c: normalizer["sorted_values"].get(c, np.array([], dtype=np.float32))
                for c in cols
            },
        }

    if _is_qbin_normalization(method):
        return {
            "method": method,
            "n_bins": int(normalizer["n_bins"]),
            "columns": cols,
            "edges": {
                c: normalizer["edges"].get(c, np.array([], dtype=np.float32))
                for c in cols
            },
        }

    raise ValueError(f"Unknown normalizer method: {method}")


def _apply_normalizer_to_frame(X: pd.DataFrame, normalizer: Dict[str, Any]) -> pd.DataFrame:
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

    if _is_qbin_normalization(method):
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


def _make_train_eval_arrays(
    X_train: pd.DataFrame,
    X_eval: pd.DataFrame,
    *,
    normalization: str,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    normalization = (normalization or "none").lower().strip()

    if normalization == "zscore":
        norm = _fit_zscore_normalizer(X_train)
        Xtr = _apply_normalizer_to_frame(X_train, norm)
        Xev = _apply_normalizer_to_frame(X_eval, norm)
        return (
            Xtr.to_numpy(dtype=np.float32, copy=False),
            Xev.to_numpy(dtype=np.float32, copy=False),
            norm,
        )

    norm = _fit_none_normalizer(X_train.columns)

    return (
        X_train.to_numpy(dtype=np.float32, copy=False),
        X_eval.to_numpy(dtype=np.float32, copy=False),
        norm,
    )


def _spearman_abs_corr_topk(X: pd.DataFrame, y: pd.Series, k: Optional[int]) -> pd.Index:
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


def _compute_prediction_metrics(predictions_df: pd.DataFrame) -> Dict[str, Any]:
    y_true = predictions_df["true"].to_numpy(dtype=float)
    y_pred = predictions_df["predicted"].to_numpy(dtype=float)
    mask = np.isfinite(y_true) & np.isfinite(y_pred)

    metrics: Dict[str, Any] = {
        "n_loo_predictions": int(len(y_true)),
        "n_valid_predictions": int(mask.sum()),
        "mse": np.nan,
        "rmse": np.nan,
        "pearson_corr": np.nan,
        "vector_cosine_similarity_all_loo": np.nan,
        "mean_cosine_similarity_across_loo_iters": np.nan,
        "mean_absolute_error": np.nan,
    }

    if mask.sum() == 0:
        return metrics

    yt = y_true[mask]
    yp = y_pred[mask]
    err = yt - yp

    metrics["mse"] = float(np.mean(err ** 2))
    metrics["rmse"] = float(np.sqrt(metrics["mse"]))
    metrics["mean_absolute_error"] = float(np.mean(np.abs(err)))

    if len(yt) >= 2 and np.nanstd(yt) > 0 and np.nanstd(yp) > 0:
        metrics["pearson_corr"] = float(np.corrcoef(yt, yp)[0, 1])

    denom = float(np.linalg.norm(yt) * np.linalg.norm(yp))
    if denom > 0:
        metrics["vector_cosine_similarity_all_loo"] = float(np.dot(yt, yp) / denom)

    scalar_cos = []

    for a, b in zip(yt, yp):
        d = abs(a) * abs(b)
        if d > 0:
            scalar_cos.append((a * b) / d)

    if scalar_cos:
        metrics["mean_cosine_similarity_across_loo_iters"] = float(np.mean(scalar_cos))

    return metrics


def _fit_and_score_cv(
    param_grid: dict,
    X_train: pd.DataFrame,
    y_train: pd.Series,
    n_splits: int,
    *,
    select_k: int = 0,
    normalization: str = "none",
    early_stopping_rounds: int = 30,
    max_param_combos: Optional[int] = None,
    use_gpu: bool = False,
    xgb_threads: Optional[int] = None,
    use_best_iteration: bool = True,
) -> dict:
    params_list = list(_sample_param_combos(param_grid, max_param_combos))

    if not params_list:
        return {}

    n_splits = min(int(n_splits), len(y_train))

    if n_splits < 2:
        return dict(params_list[0])

    avg_mse = np.zeros(len(params_list), dtype=float)
    best_rounds: List[List[int]] = [[] for _ in params_list]

    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

    for tr_idx, va_idx in kf.split(X_train):
        Xtr_df, Xva_df = X_train.iloc[tr_idx], X_train.iloc[va_idx]
        ytr_s, yva_s = y_train.iloc[tr_idx], y_train.iloc[va_idx]

        if select_k and int(select_k) > 0:
            kept = _spearman_abs_corr_topk(Xtr_df, ytr_s, int(select_k))
            Xtr_df = Xtr_df[kept]
            Xva_df = Xva_df[kept]

        Xtr, Xva, _ = _make_train_eval_arrays(Xtr_df, Xva_df, normalization=normalization)
        ytr, yva = _to_f32_series(ytr_s), _to_f32_series(yva_s)

        for j, params in enumerate(params_list):
            try:
                est = _fit_xgb(
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
                    best_rounds[j].append(
                        _best_iteration_count(est, int(params.get("n_estimators", 100)))
                    )

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


def _extract_gain_importance(
    est: xgb.XGBRegressor,
    columns: Sequence[str],
) -> Optional[Dict[str, float]]:
    try:
        gain = est.get_booster().get_score(importance_type="gain")
        fmap = {f"f{i}": str(col) for i, col in enumerate(columns)}
        return {fmap.get(k, k): float(v) for k, v in gain.items()}
    except Exception:
        return None


def _fit_model_and_make_artifact(
    params: dict,
    X_train_sel: pd.DataFrame,
    y_train: pd.Series,
    *,
    protein_name: str,
    normalization: str,
    global_normalizer: Dict[str, Any],
    use_gpu: bool = False,
    xgb_threads: Optional[int] = None,
    artifact_kind: str = "loo_model",
    extra_metadata: Optional[Dict[str, Any]] = None,
) -> Tuple[xgb.XGBRegressor, Optional[Dict[str, float]], Dict[str, Any]]:
    normalization = (normalization or "none").lower().strip()

    if normalization == "zscore":
        artifact_normalizer = _fit_zscore_normalizer(X_train_sel)
        Xtr_df = _apply_normalizer_to_frame(X_train_sel, artifact_normalizer)

    elif normalization == "rank" or _is_qbin_normalization(normalization):
        artifact_normalizer = _slice_normalizer(global_normalizer, X_train_sel.columns)
        Xtr_df = X_train_sel.astype(np.float32, copy=False)

    else:
        artifact_normalizer = _fit_none_normalizer(X_train_sel.columns)
        Xtr_df = X_train_sel.astype(np.float32, copy=False)

    est = _fit_xgb(
        _clean_params(params),
        Xtr_df.to_numpy(dtype=np.float32, copy=False),
        _to_f32_series(y_train),
        use_gpu=use_gpu,
        xgb_threads=xgb_threads,
    )

    feat_imp = _extract_gain_importance(est, X_train_sel.columns)

    artifact = {
        "artifact_kind": artifact_kind,
        "protein": protein_name,
        "model": est,
        "params": _clean_params(params),
        "params_with_cv_metadata": dict(params or {}),
        "selected_features": list(map(str, X_train_sel.columns)),
        "normalization": normalization,
        "normalizer": artifact_normalizer,
        "feature_selector": "spearman_prefilter",
        "n_training_samples": int(len(y_train)),
        "training_samples": list(map(str, X_train_sel.index)),
        "xgboost_version": xgb.__version__,
        "use_gpu": bool(use_gpu),
        "xgb_threads": _default_xgb_threads(xgb_threads),
        "feature_importance_gain": feat_imp,
    }

    if extra_metadata:
        artifact.update(extra_metadata)

    return est, feat_imp, artifact


def _fit_final_and_predict(
    params: dict,
    X_train_sel: pd.DataFrame,
    y_train: pd.Series,
    X_test_sel: pd.DataFrame,
    *,
    protein_name: str,
    normalization: str,
    global_normalizer: Dict[str, Any],
    use_gpu: bool = False,
    xgb_threads: Optional[int] = None,
    extra_metadata: Optional[Dict[str, Any]] = None,
):
    normalization = (normalization or "none").lower().strip()

    try:
        est, feat_imp, artifact = _fit_model_and_make_artifact(
            params,
            X_train_sel,
            y_train,
            protein_name=protein_name,
            normalization=normalization,
            global_normalizer=global_normalizer,
            use_gpu=use_gpu,
            xgb_threads=xgb_threads,
            artifact_kind="best_loo_candidate_model",
            extra_metadata=extra_metadata,
        )

        if normalization == "zscore":
            Xte_df = _apply_normalizer_to_frame(X_test_sel, artifact["normalizer"])
        else:
            Xte_df = X_test_sel.astype(np.float32, copy=False)

        y_pred = est.predict(Xte_df.to_numpy(dtype=np.float32, copy=False))

        return y_pred, feat_imp, artifact

    except Exception as e:
        print(f"WARNING: final fit failed: {type(e).__name__}: {e}")
        return np.array([np.nan], dtype=np.float32), None, None


def _train_final_all_samples_artifact(
    *,
    param_grid: dict,
    best_params_from_once: Optional[dict],
    X_work: pd.DataFrame,
    y_protein: pd.Series,
    protein_name: str,
    n_splits: int,
    select_k: int,
    normalization: str,
    global_normalizer: Dict[str, Any],
    early_stopping_rounds: int,
    max_param_combos: Optional[int],
    use_gpu: bool,
    xgb_threads: Optional[int],
    grid_preset: str,
    tuning_strategy: str,
    proteomics_orientation: str,
) -> Optional[Dict[str, Any]]:
    try:
        if best_params_from_once is not None:
            final_params = dict(best_params_from_once)
        else:
            print("Tuning final all-samples model parameters...")
            final_params = _fit_and_score_cv(
                param_grid,
                X_work,
                y_protein,
                int(n_splits),
                select_k=int(select_k),
                normalization=normalization,
                early_stopping_rounds=int(early_stopping_rounds),
                max_param_combos=max_param_combos,
                use_gpu=use_gpu,
                xgb_threads=xgb_threads,
            )

        kept = _spearman_abs_corr_topk(X_work, y_protein, int(select_k))

        _, _, artifact = _fit_model_and_make_artifact(
            final_params,
            X_work[kept],
            y_protein,
            protein_name=protein_name,
            normalization=normalization,
            global_normalizer=global_normalizer,
            use_gpu=use_gpu,
            xgb_threads=xgb_threads,
            artifact_kind="final_all_samples_model",
            extra_metadata={
                "select_k": int(select_k),
                "proteomics_orientation": proteomics_orientation,
                "grid_preset": grid_preset,
                "tuning_strategy_used_for_loo": tuning_strategy,
                "final_model_note": "Trained on all available samples for this protein after LOO evaluation.",
            },
        )

        return artifact

    except Exception as e:
        print(f"WARNING: final all-samples model training failed: {type(e).__name__}: {e}")
        return None


def predict_with_artifact(artifact: Dict[str, Any], X_new: pd.DataFrame) -> pd.Series:
    protein = artifact.get("protein", "prediction")
    features = list(map(str, artifact["selected_features"]))

    X_sel = X_new.reindex(columns=features)
    X_trans = _apply_normalizer_to_frame(
        X_sel,
        artifact.get("normalizer", _fit_none_normalizer(features)),
    )

    pred = artifact["model"].predict(X_trans.to_numpy(dtype=np.float32, copy=False))

    return pd.Series(pred, index=X_new.index, name=str(protein))


def _extract_protein_series(y: pd.DataFrame, protein_name: str) -> Tuple[pd.Series, str]:
    if protein_name in y.index:
        s = y.loc[protein_name].copy()
        s.name = protein_name
        return s, "proteins_as_rows"

    if protein_name in y.columns:
        s = y[protein_name].copy()
        s.name = protein_name
        return s, "proteins_as_columns_legacy"

    raise ValueError(
        f"Protein {protein_name!r} not found in proteomics rows or columns. "
        "Expected proteomics_matched orientation: rows=proteins, columns=samples, "
        "or rows=samples, columns=proteins."
    )


def _align_X_and_y_for_protein(
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
            raise ValueError(
                "No overlapping sample IDs between transcriptomics rows and proteomics samples. "
                "Check that transcriptomics_filtered rows are samples and proteomics_matched has matching samples."
            )

        X = X2
        y_protein = y2

    if n_samples is not None:
        common = common[: int(n_samples)]

    return X.loc[common], y_protein.loc[common]


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


def save_minimal_protein_h5(
    result: Dict[str, Any],
    path: str,
    *,
    protein_name: str,
    model_choice: str = "best_loo",
    compression_level: int = 9,
) -> str:
    try:
        import h5py
    except ImportError as e:
        raise ImportError("minimal_h5 output requires h5py to be installed.") from e

    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    predictions_df = result["predictions"]
    metrics = dict(result["metrics"])

    model_choice = (model_choice or "best_loo").lower()

    if model_choice not in {"best_loo", "final"}:
        raise ValueError("model_choice must be 'best_loo' or 'final'.")

    if model_choice == "final":
        model_artifact = result.get("final_model_artifact")
        model_label = "final_model_artifact"
    else:
        model_artifact = result.get("best_loo_model_artifact")
        model_label = "best_loo_model_artifact"

    if model_artifact is None:
        raise ValueError("No model artifact is available. minimal_h5 output requires save_models=True.")

    needed = ["LOO sample", "true", "predicted", "squared_error", "absolute_error"]
    missing = [c for c in needed if c not in predictions_df.columns]

    if missing:
        raise ValueError(f"Predictions dataframe is missing required columns: {missing}")

    compact_artifact = dict(model_artifact)
    compact_artifact.pop("final_model_summary", None)

    metrics_json_obj = dict(metrics)
    metrics_json_obj.update(
        {
            "protein": protein_name,
            "h5_model_label": model_label,
            "h5_contents": {
                "metrics": "protein-level LOO metrics including pearson_corr, mse and cosine similarity",
                "loo": "LOO sample IDs, true values, predictions and errors",
                "model": "model artifact stored as pickle+xz bytes",
            },
        }
    )

    metrics_json = json.dumps(
        _json_safe(metrics_json_obj),
        ensure_ascii=False,
        sort_keys=True,
    ).encode("utf-8")

    params_json = json.dumps(
        _json_safe(compact_artifact.get("params", {})),
        ensure_ascii=False,
        sort_keys=True,
    ).encode("utf-8")

    summary_json = json.dumps(
        _json_safe(compact_artifact.get("best_loo_model_summary", metrics.get("best_loo_model_summary"))),
        ensure_ascii=False,
        sort_keys=True,
    ).encode("utf-8")

    pred = predictions_df[needed].copy()

    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "minimal_protein_loo_xgboost_h5"
        h5.attrs["format_version"] = "1.1"
        h5.attrs["protein"] = str(protein_name)
        h5.attrs["model_label"] = model_label
        h5.attrs["compression"] = (
            f"hdf5_gzip_{int(compression_level)}_shuffle_fletcher32_plus_model_pickle_xz_extreme"
        )

        g = h5.create_group("metrics")
        _write_dataset_compressed(
            g,
            "json",
            np.frombuffer(metrics_json, dtype=np.uint8),
            compression_level=compression_level,
        )
        g.attrs["pearson_corr"] = np.nan if metrics.get("pearson_corr") is None else float(metrics.get("pearson_corr", np.nan))
        g.attrs["mse"] = np.nan if metrics.get("mse") is None else float(metrics.get("mse", np.nan))
        g.attrs["rmse"] = np.nan if metrics.get("rmse") is None else float(metrics.get("rmse", np.nan))
        g.attrs["cosine_similarity"] = np.nan if metrics.get("vector_cosine_similarity_all_loo") is None else float(metrics.get("vector_cosine_similarity_all_loo", np.nan))
        g.attrs["n_valid_predictions"] = int(metrics.get("n_valid_predictions", 0) or 0)

        g = h5.create_group("loo")
        _write_dataset_compressed(
            g,
            "samples",
            _fixed_utf8_array(pred["LOO sample"].astype(str)),
            compression_level=compression_level,
        )
        _write_dataset_compressed(
            g,
            "true",
            pred["true"].to_numpy(dtype=np.float64),
            compression_level=compression_level,
        )
        _write_dataset_compressed(
            g,
            "predicted",
            pred["predicted"].to_numpy(dtype=np.float64),
            compression_level=compression_level,
        )
        _write_dataset_compressed(
            g,
            "squared_error",
            pred["squared_error"].to_numpy(dtype=np.float64),
            compression_level=compression_level,
        )
        _write_dataset_compressed(
            g,
            "absolute_error",
            pred["absolute_error"].to_numpy(dtype=np.float64),
            compression_level=compression_level,
        )

        g = h5.create_group("model")
        g.attrs["artifact_label"] = model_label
        g.attrs["artifact_kind"] = str(compact_artifact.get("artifact_kind"))
        g.attrs["selected_feature_count"] = int(len(compact_artifact.get("selected_features", [])))

        _write_dataset_compressed(
            g,
            "artifact_pickle_xz",
            _pickle_xz_bytes(compact_artifact),
            compression_level=compression_level,
        )
        _write_dataset_compressed(
            g,
            "selected_features",
            _fixed_utf8_array(compact_artifact.get("selected_features", [])),
            compression_level=compression_level,
        )
        _write_dataset_compressed(
            g,
            "params_json",
            np.frombuffer(params_json, dtype=np.uint8),
            compression_level=compression_level,
        )
        _write_dataset_compressed(
            g,
            "best_loo_model_summary_json",
            np.frombuffer(summary_json, dtype=np.uint8),
            compression_level=compression_level,
        )

    return path


def load_minimal_protein_h5(path: str) -> Dict[str, Any]:
    try:
        import h5py
    except ImportError as e:
        raise ImportError("Loading minimal H5 output requires h5py to be installed.") from e

    with h5py.File(path, "r") as h5:
        metrics = json.loads(bytes(h5["metrics/json"][:]).decode("utf-8"))

        samples = [x.decode("utf-8") for x in h5["loo/samples"][:]]

        predictions = pd.DataFrame(
            {
                "LOO sample": samples,
                "true": h5["loo/true"][:],
                "predicted": h5["loo/predicted"][:],
                "squared_error": h5["loo/squared_error"][:],
                "absolute_error": h5["loo/absolute_error"][:],
            }
        )
        predictions.insert(0, "protein", metrics.get("protein", h5.attrs.get("protein", "")))

        model_artifact = _unpickle_xz_bytes(h5["model/artifact_pickle_xz"][:])

    return {
        "metrics": metrics,
        "predictions": predictions,
        "model_artifact": model_artifact,
    }


def _make_summary_row(
    predictions_df: pd.DataFrame,
    metrics: Dict[str, Any],
    protein_name: str,
    model_artifact: Optional[Dict[str, Any]],
    label: str,
) -> pd.DataFrame:
    pred = predictions_df[["LOO sample", "true", "predicted", "squared_error", "absolute_error"]].copy()

    return pd.DataFrame(
        [
            {
                "protein": protein_name,
                "pearson_corr": metrics.get("pearson_corr", np.nan),
                "mse": metrics.get("mse", np.nan),
                "rmse": metrics.get("rmse", np.nan),
                "cosine_similarity": metrics.get("vector_cosine_similarity_all_loo", np.nan),
                "n_loo_predictions": metrics.get("n_loo_predictions", len(predictions_df)),
                "n_valid_predictions": metrics.get("n_valid_predictions", np.nan),
                "loo_samples": pred["LOO sample"].astype(str).tolist(),
                "loo_true": pred["true"].astype(float).tolist(),
                "loo_predicted": pred["predicted"].astype(float).tolist(),
                "model_artifact_label": label,
                "model_artifact": model_artifact,
                "model_artifact_kind": None if model_artifact is None else model_artifact.get("artifact_kind"),
                "model_selected_feature_count": None if model_artifact is None else len(model_artifact.get("selected_features", [])),
                "model_params": None if model_artifact is None else model_artifact.get("params"),
            }
        ]
    )


def cv_pipeline_loo(
    model,
    param_grid: dict,
    X: pd.DataFrame,
    y: pd.DataFrame,
    protein_name: str,
    *,
    feature_selector: str = "spearman",
    n_splits: int = 3,
    n_samples: Optional[int] = None,
    select_k: int = 500,
    preselect_k: Optional[int] = 5000,
    early_stopping_rounds: int = 30,
    minimal_save: bool = True,
    max_param_combos: Optional[int] = None,
    use_gpu: bool = False,
    normalization: str = "none",
    tuning_strategy: str = "nested",
    xgb_threads: Optional[int] = None,
    grid_preset: str = "fast",
    save_models: bool = True,
    save_final_model: bool = True,
    minimal_summary_save: bool = False,
    summary_model: str = "final",
) -> Dict[str, Any]:
    y_protein, proteomics_orientation = _extract_protein_series(y, protein_name)

    normalization = (normalization or "none").lower().strip()
    tuning_strategy = (tuning_strategy or "nested").lower().strip()
    summary_model = (summary_model or "final").lower().strip()

    if normalization not in {"none", "zscore", "rank"} and not _is_qbin_normalization(normalization):
        raise ValueError(
            f"Invalid normalization: {normalization}. "
            "Allowed: none, zscore, rank, or qbinN such as qbin50."
        )

    if tuning_strategy not in {"nested", "once"}:
        raise ValueError("tuning_strategy must be 'nested' or 'once'.")

    if summary_model not in {"final", "best_loo"}:
        raise ValueError("summary_model must be 'final' or 'best_loo'.")

    X_raw, y_protein = _align_X_and_y_for_protein(X, y_protein, n_samples=n_samples)

    valid = y_protein.notna()
    X_raw = X_raw.loc[valid]
    y_protein = y_protein.loc[valid]

    sample_idx = X_raw.index

    if len(sample_idx) < max(3, int(n_splits) + 1):
        raise ValueError(
            f"Not enough non-missing samples for protein {protein_name}: "
            f"n={len(sample_idx)}, n_splits={n_splits}"
        )

    global_normalizer = _fit_global_normalizer(X_raw, normalization)

    if normalization == "rank" or _is_qbin_normalization(normalization):
        print(f"Applying global {normalization} transform once outside LOO...")
        X_work = _apply_normalizer_to_frame(X_raw, global_normalizer)
    else:
        X_work = X_raw

    global_best_params = None

    if tuning_strategy == "once":
        print("Tuning hyperparameters once for this protein...")

        global_best_params = _fit_and_score_cv(
            param_grid,
            X_work,
            y_protein,
            int(n_splits),
            select_k=int(select_k),
            normalization=normalization,
            early_stopping_rounds=int(early_stopping_rounds),
            max_param_combos=max_param_combos,
            use_gpu=use_gpu,
            xgb_threads=xgb_threads,
        )

        print(f"Best once-tuned params: {global_best_params}")

    rows = []
    best_loo_model_artifact = None
    best_loo_error = np.inf
    best_loo_summary = None

    for i, test_id in enumerate(sample_idx, start=1):
        print(f"LOO {i}/{len(sample_idx)}: held-out sample={test_id}")

        X_train = X_work.drop(index=[test_id])
        y_train = y_protein.drop(index=[test_id])
        X_test = X_work.loc[[test_id]]
        y_test = y_protein.loc[[test_id]]

        kept = _spearman_abs_corr_topk(X_train, y_train, int(select_k))

        Xtr_sel = X_train[kept]
        Xte_sel = X_test[kept]

        if tuning_strategy == "nested":
            best_params = _fit_and_score_cv(
                param_grid,
                Xtr_sel,
                y_train,
                int(n_splits),
                select_k=0,
                normalization=normalization,
                early_stopping_rounds=int(early_stopping_rounds),
                max_param_combos=max_param_combos,
                use_gpu=use_gpu,
                xgb_threads=xgb_threads,
            )
        else:
            best_params = dict(global_best_params or {})

        y_pred, feat_imp, candidate_artifact = _fit_final_and_predict(
            best_params,
            Xtr_sel,
            y_train,
            Xte_sel,
            protein_name=protein_name,
            normalization=normalization,
            global_normalizer=global_normalizer,
            use_gpu=use_gpu,
            xgb_threads=xgb_threads,
            extra_metadata={
                "held_out_sample": str(test_id),
                "proteomics_orientation": proteomics_orientation,
                "select_k": int(select_k),
                "grid_preset": grid_preset,
                "tuning_strategy": tuning_strategy,
            },
        )

        pred_val = float(y_pred[0]) if isinstance(y_pred, np.ndarray) else float(y_pred)
        true_val = float(y_test.values[0])

        squared_error = float((true_val - pred_val) ** 2) if np.isfinite(pred_val) else np.inf
        absolute_error = float(abs(true_val - pred_val)) if np.isfinite(pred_val) else np.inf

        if save_models and candidate_artifact is not None and squared_error < best_loo_error:
            best_loo_error = squared_error
            best_loo_model_artifact = candidate_artifact
            best_loo_summary = {
                "held_out_sample": str(test_id),
                "true": true_val,
                "predicted": pred_val,
                "squared_error": squared_error,
                "absolute_error": absolute_error,
                "params": dict(best_params or {}),
                "selected_feature_count": int(len(kept)),
            }
            best_loo_model_artifact["best_loo_model_summary"] = best_loo_summary

        row = {
            "protein": protein_name,
            "LOO sample": test_id,
            "true": true_val,
            "predicted": pred_val,
            "squared_error": squared_error if np.isfinite(squared_error) else np.nan,
            "absolute_error": absolute_error if np.isfinite(absolute_error) else np.nan,
        }

        if not minimal_save:
            row.update(
                {
                    "params_dict": {
                        **best_params,
                        "feature_selector": "spearman_prefilter",
                        "select_k": int(select_k),
                        "normalization": normalization,
                        "tuning_strategy": tuning_strategy,
                        "grid_preset": grid_preset,
                        "early_stopping_rounds": int(early_stopping_rounds),
                        "inner_cv_folds": int(n_splits),
                        "xgb_threads": _default_xgb_threads(xgb_threads),
                        "use_gpu": bool(use_gpu),
                    },
                    "features_coef_dict": feat_imp,
                    "features_selected": list(map(str, kept)),
                }
            )

        rows.append(row)

    predictions_df = pd.DataFrame(rows)
    metrics = _compute_prediction_metrics(predictions_df)

    metrics.update(
        {
            "protein": protein_name,
            "proteomics_orientation": proteomics_orientation,
            "normalization": normalization,
            "select_k": int(select_k),
            "grid_preset": grid_preset,
            "tuning_strategy": tuning_strategy,
            "n_splits_inner_cv": int(n_splits),
            "early_stopping_rounds": int(early_stopping_rounds),
            "max_param_combos": None if max_param_combos is None else int(max_param_combos),
            "use_gpu": bool(use_gpu),
            "xgb_threads": _default_xgb_threads(xgb_threads),
            "minimal_summary_save": bool(minimal_summary_save),
            "summary_model": summary_model,
            "best_loo_model_summary": best_loo_summary,
        }
    )

    final_model_artifact = None

    if save_models and save_final_model:
        final_model_artifact = _train_final_all_samples_artifact(
            param_grid=param_grid,
            best_params_from_once=global_best_params,
            X_work=X_work,
            y_protein=y_protein,
            protein_name=protein_name,
            n_splits=int(n_splits),
            select_k=int(select_k),
            normalization=normalization,
            global_normalizer=global_normalizer,
            early_stopping_rounds=int(early_stopping_rounds),
            max_param_combos=max_param_combos,
            use_gpu=use_gpu,
            xgb_threads=xgb_threads,
            grid_preset=grid_preset,
            tuning_strategy=tuning_strategy,
            proteomics_orientation=proteomics_orientation,
        )

        if final_model_artifact is not None:
            metrics["final_model_summary"] = {
                "artifact_kind": final_model_artifact.get("artifact_kind"),
                "n_training_samples": final_model_artifact.get("n_training_samples"),
                "selected_feature_count": len(final_model_artifact.get("selected_features", [])),
                "params": final_model_artifact.get("params"),
            }

    if summary_model == "best_loo":
        summary_artifact = best_loo_model_artifact
        summary_label = "best_loo_model_artifact"
    else:
        summary_artifact = final_model_artifact
        summary_label = "final_model_artifact"

    summary_row = _make_summary_row(
        predictions_df,
        metrics,
        protein_name,
        summary_artifact,
        summary_label,
    )

    return {
        "predictions": predictions_df,
        "metrics": metrics,
        "summary_row": summary_row,
        "best_loo_model_artifact": best_loo_model_artifact,
        "final_model_artifact": final_model_artifact,
    }