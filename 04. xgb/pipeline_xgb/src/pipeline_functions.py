import itertools
import os
from typing import Optional, Sequence, Tuple, Dict, Any, List

import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import KFold
import xgboost as xgb

np.random.seed(42)


# ----------------------
# Hyperparameter utilities
# ----------------------
def _all_param_combinations(param_grid: dict) -> Sequence[dict]:
    """Cartesian product of a parameter grid. Empty grid -> [{}]."""
    if not param_grid:
        return [dict()]
    keys, vals = zip(*param_grid.items())
    return [dict(zip(keys, combo)) for combo in itertools.product(*vals)]


def _sample_param_combos(param_grid: dict, max_param_combos: Optional[int]) -> Sequence[dict]:
    """Optionally sample a reproducible subset of parameter combinations."""
    combos = list(_all_param_combinations(param_grid))
    if max_param_combos is not None and len(combos) > int(max_param_combos):
        rng = np.random.default_rng(42)
        picks = rng.choice(len(combos), size=int(max_param_combos), replace=False)
        combos = [combos[i] for i in picks]
    return combos


def _default_xgb_threads(xgb_threads: Optional[int]) -> int:
    """Pick a thread count that respects SLURM when available."""
    if xgb_threads is not None and int(xgb_threads) > 0:
        return int(xgb_threads)
    slurm_threads = os.getenv("SLURM_CPUS_PER_TASK")
    if slurm_threads and slurm_threads.isdigit():
        return max(1, int(slurm_threads))
    return max(1, os.cpu_count() or 1)


def _clean_params(params: Optional[dict]) -> dict:
    """Remove metadata keys before passing parameters to XGBoost."""
    return {k: v for k, v in (params or {}).items() if not str(k).startswith("_")}


def _clone_xgb_with_params(
    params: Optional[dict],
    *,
    seed: int = 42,
    use_gpu: bool = False,
    xgb_threads: Optional[int] = None,
) -> xgb.XGBRegressor:
    """Create XGBRegressor with explicit threads and a fixed CPU/GPU backend."""
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
    """
    Fit XGBoost while tolerating sklearn-API differences across XGBoost versions.

    Some versions accept early_stopping_rounds in fit(); newer versions prefer it
    in the constructor. This helper tries fit() first, then constructor fallback.
    """
    use_es = Xva is not None and yva is not None and early_stopping_rounds and early_stopping_rounds > 0

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
    """Return number of trees selected by early stopping, or fallback."""
    best = getattr(est, "best_iteration", None)
    if best is None:
        try:
            best = est.get_booster().best_iteration
        except Exception:
            best = None
    if best is None:
        return int(fallback)
    return max(1, int(best) + 1)  # best_iteration is zero-based.


# ----------------------
# Normalization utilities
# ----------------------
def _to_f32_series(s: pd.Series) -> np.ndarray:
    return s.to_numpy(dtype=np.float32, copy=False)


def _global_normalize_X(X: pd.DataFrame, normalization: str) -> pd.DataFrame:
    """
    Global unsupervised transformations done once before the LOO loop.

    none:
        Return X unchanged.
    zscore:
        Return X unchanged here. Z-score is fitted inside each LOO/train fold.
    rank:
        Rank-transform each gene across samples once outside LOO.
    qbin5/qbin10:
        Quantile-bin each gene once outside LOO and scale bins to [0, 1].
    """
    normalization = (normalization or "none").lower()

    if normalization in {"none", "zscore"}:
        return X

    if normalization == "rank":
        return X.rank(axis=0, method="average", pct=True).astype(np.float32)

    if normalization.startswith("qbin"):
        n_bins = int(normalization.replace("qbin", ""))
        if n_bins < 2:
            raise ValueError("qbin normalization requires at least 2 bins.")

        Xb = pd.DataFrame(index=X.index, columns=X.columns, dtype=np.float32)
        for col in X.columns:
            s = X[col]
            ranks = s.rank(method="average", pct=True)
            b = np.floor(ranks.to_numpy(dtype=np.float64) * n_bins)
            b = np.clip(b, 0, n_bins - 1)
            b = b / max(1, n_bins - 1)
            b[pd.isna(s).to_numpy()] = np.nan
            Xb[col] = b.astype(np.float32)
        return Xb

    raise ValueError(f"Unknown normalization: {normalization}")


def _zscore_fit_transform_train_test(
    X_train: pd.DataFrame,
    X_eval: pd.DataFrame,
) -> Tuple[np.ndarray, np.ndarray]:
    """Fit z-score on X_train only, then apply to X_train and X_eval."""
    Xtr = X_train.to_numpy(dtype=np.float32, copy=False)
    Xev = X_eval.to_numpy(dtype=np.float32, copy=False)

    mean = np.nanmean(Xtr, axis=0).astype(np.float32)
    std = np.nanstd(Xtr, axis=0).astype(np.float32)
    mean[~np.isfinite(mean)] = 0.0
    std[~np.isfinite(std) | (std == 0)] = 1.0

    Xtr_z = ((Xtr - mean) / std).astype(np.float32, copy=False)
    Xev_z = ((Xev - mean) / std).astype(np.float32, copy=False)
    return Xtr_z, Xev_z


def _make_train_eval_arrays(
    X_train: pd.DataFrame,
    X_eval: pd.DataFrame,
    *,
    normalization: str = "none",
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert train/eval data to float32 arrays.

    zscore is local: fit on X_train only.
    none/rank/qbin are already represented in the dataframe.
    """
    normalization = (normalization or "none").lower()
    if normalization == "zscore":
        return _zscore_fit_transform_train_test(X_train, X_eval)
    return (
        X_train.to_numpy(dtype=np.float32, copy=False),
        X_eval.to_numpy(dtype=np.float32, copy=False),
    )


# ----------------------
# Spearman prefilter
# ----------------------
def _spearman_abs_corr_topk(X: pd.DataFrame, y: pd.Series, k: Optional[int]) -> pd.Index:
    """Keep top-k columns by absolute Spearman correlation with y."""
    if k is None or int(k) <= 0 or X.shape[1] <= int(k):
        return X.columns
    Xr = X.rank(axis=0, method="average")
    yr = y.rank(method="average")
    s = Xr.corrwith(yr).abs().replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return s.nlargest(int(k)).index


# ----------------------
# Inner CV tuning
# ----------------------
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
    """
    KFold CV for hyperparameter tuning.

    If select_k > 0, Spearman filtering is performed inside each inner fold using
    only that fold's training samples. If X_train is already prefiltered by the
    outer split, pass select_k=0.
    """
    params_list = list(_sample_param_combos(param_grid, max_param_combos))
    if len(params_list) == 0:
        return {}

    avg_mse = np.zeros(len(params_list), dtype=float)
    best_rounds: List[List[int]] = [[] for _ in params_list]
    kf = KFold(n_splits=int(n_splits), shuffle=True, random_state=42)

    for tr_idx, va_idx in kf.split(X_train):
        Xtr_df, Xva_df = X_train.iloc[tr_idx], X_train.iloc[va_idx]
        ytr_s, yva_s = y_train.iloc[tr_idx], y_train.iloc[va_idx]

        if select_k and int(select_k) > 0:
            kept = _spearman_abs_corr_topk(Xtr_df, ytr_s, int(select_k))
            Xtr_df = Xtr_df[kept]
            Xva_df = Xva_df[kept]

        Xtr, Xva = _make_train_eval_arrays(Xtr_df, Xva_df, normalization=normalization)
        ytr, yva = _to_f32_series(ytr_s), _to_f32_series(yva_s)

        for k_idx, params in enumerate(params_list):
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
                avg_mse[k_idx] += mean_squared_error(yva, pred) / int(n_splits)
                if use_best_iteration and early_stopping_rounds and early_stopping_rounds > 0:
                    fallback = int(params.get("n_estimators", 100))
                    best_rounds[k_idx].append(_best_iteration_count(est, fallback))
            except Exception as e:
                print(f"WARNING: CV fit failed for params={params}: {type(e).__name__}: {e}")
                avg_mse[k_idx] += 1e6 / int(n_splits)

    best_idx = int(np.argmin(avg_mse))
    best_params = dict(params_list[best_idx])

    if use_best_iteration and best_rounds[best_idx]:
        original = int(best_params.get("n_estimators", 100))
        median_best = int(np.median(best_rounds[best_idx]))
        adjusted = max(25, min(original, int(np.ceil(median_best * 1.10))))
        best_params["n_estimators"] = adjusted
        best_params["_cv_median_best_iteration"] = median_best

    best_params["_cv_mse"] = float(avg_mse[best_idx])
    return best_params


# ----------------------
# Final fit and predict
# ----------------------
def _fit_final_and_predict(
    params: dict,
    X_train_sel: pd.DataFrame,
    y_train: pd.Series,
    X_test_sel: pd.DataFrame,
    *,
    normalization: str = "none",
    use_gpu: bool = False,
    xgb_threads: Optional[int] = None,
):
    Xtr, Xte = _make_train_eval_arrays(X_train_sel, X_test_sel, normalization=normalization)
    ytr = _to_f32_series(y_train)

    try:
        est = _fit_xgb(_clean_params(params), Xtr, ytr, use_gpu=use_gpu, xgb_threads=xgb_threads)
        y_pred = est.predict(Xte)

        feat_imp = None
        try:
            gain = est.get_booster().get_score(importance_type="gain")
            fmap = {f"f{i}": col for i, col in enumerate(X_train_sel.columns)}
            feat_imp = {fmap.get(k, k): float(v) for k, v in gain.items()}
        except Exception:
            pass

        return y_pred, feat_imp
    except Exception as e:
        print(f"WARNING: final fit failed: {type(e).__name__}: {e}")
        return np.array([np.nan], dtype=np.float32), None


# ----------------------
# Unified LOO pipeline
# ----------------------
def cv_pipeline_loo(
    model,  # kept for signature compatibility; not used
    param_grid: dict,
    X: pd.DataFrame,
    y: pd.DataFrame,
    protein_name: str,
    *,
    feature_selector: str = "spearman",  # kept for compatibility
    n_splits: int = 3,
    n_samples: Optional[int] = None,
    select_k: int = 500,
    preselect_k: Optional[int] = 5000,  # kept for compatibility; ignored
    early_stopping_rounds: int = 30,
    minimal_save: bool = True,
    max_param_combos: Optional[int] = None,
    use_gpu: bool = False,
    normalization: str = "none",
    tuning_strategy: str = "nested",
    xgb_threads: Optional[int] = None,
    grid_preset: str = "fast",
) -> pd.DataFrame:
    """
    Outer evaluation is always Leave-One-Out CV over samples.

    tuning_strategy='nested':
        Strict original-style LOO evaluation. Each held-out sample gets its own
        hyperparameter tuning on the outer training set.

    tuning_strategy='once':
        Faster mode. Hyperparameters are tuned once per protein and then reused
        for all LOO fits. The predictions are still LOO predictions, but the
        hyperparameter choice has seen all samples, so this is less strict for
        unbiased benchmarking.

    normalization:
        none   = raw transcriptomics.
        zscore = fitted inside each train/eval split, including each LOO split.
        rank   = rank-transformed once outside LOO.
        qbin5/qbin10 = binned once outside LOO.
    """
    if protein_name not in y.columns:
        raise ValueError(f"Protein {protein_name} not found in y columns.")

    normalization = (normalization or "none").lower()
    if normalization not in {"none", "zscore", "rank", "qbin5", "qbin10"}:
        raise ValueError(f"Invalid normalization: {normalization}")

    tuning_strategy = (tuning_strategy or "nested").lower()
    if tuning_strategy not in {"nested", "once"}:
        raise ValueError("tuning_strategy must be 'nested' or 'once'.")

    y_protein = y[protein_name]
    sample_idx = X.index if n_samples is None else X.index[:int(n_samples)]
    X = X.loc[sample_idx]
    y_protein = y_protein.loc[sample_idx]

    # Drop samples with missing target for this protein.
    valid = y_protein.notna()
    X = X.loc[valid]
    y_protein = y_protein.loc[valid]
    sample_idx = X.index

    if len(sample_idx) < max(3, int(n_splits) + 1):
        raise ValueError(
            f"Not enough non-missing samples for protein {protein_name}: "
            f"n={len(sample_idx)}, n_splits={n_splits}"
        )

    # rank/qbin are intentionally computed once outside LOO.
    # zscore is intentionally not done here; it is fitted in each train/eval split.
    X = _global_normalize_X(X, normalization)

    global_best_params = None
    if tuning_strategy == "once":
        print("Tuning hyperparameters once for this protein...")
        global_best_params = _fit_and_score_cv(
            param_grid=param_grid,
            X_train=X,
            y_train=y_protein,
            n_splits=int(n_splits),
            select_k=int(select_k),
            normalization=normalization,
            early_stopping_rounds=int(early_stopping_rounds),
            max_param_combos=max_param_combos,
            use_gpu=use_gpu,
            xgb_threads=xgb_threads,
            use_best_iteration=True,
        )
        print(f"Best once-tuned params: {global_best_params}")

    rows = []
    for i, test_id in enumerate(sample_idx, start=1):
        print(f"LOO {i}/{len(sample_idx)}: held-out sample={test_id}")

        X_train = X.drop(index=[test_id])
        y_train = y_protein.drop(index=[test_id])
        X_test = X.loc[[test_id]]
        y_test = y_protein.loc[[test_id]]

        # Spearman prefilter on the outer training set only.
        kept_features = _spearman_abs_corr_topk(X_train, y_train, int(select_k))
        Xtr_sel = X_train[kept_features]
        Xte_sel = X_test[kept_features]

        if tuning_strategy == "nested":
            # Original-style: the outer training set is prefiltered once, then
            # inner CV tunes hyperparameters only on this fixed feature set.
            best_params = _fit_and_score_cv(
                param_grid=param_grid,
                X_train=Xtr_sel,
                y_train=y_train,
                n_splits=int(n_splits),
                select_k=0,
                normalization=normalization,
                early_stopping_rounds=int(early_stopping_rounds),
                max_param_combos=max_param_combos,
                use_gpu=use_gpu,
                xgb_threads=xgb_threads,
                use_best_iteration=True,
            )
        else:
            best_params = dict(global_best_params or {})

        y_pred, feat_imp = _fit_final_and_predict(
            best_params,
            Xtr_sel,
            y_train,
            Xte_sel,
            normalization=normalization,
            use_gpu=use_gpu,
            xgb_threads=xgb_threads,
        )
        pred_val = float(y_pred[0]) if isinstance(y_pred, np.ndarray) else float(y_pred)

        if minimal_save:
            rows.append({
                "protein": protein_name,
                "LOO sample": test_id,
                "true": float(y_test.values[0]),
                "predicted": pred_val,
            })
        else:
            rows.append({
                "protein": protein_name,
                "LOO sample": test_id,
                "true": float(y_test.values[0]),
                "predicted": pred_val,
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
                "features_selected": list(map(str, kept_features)),
            })

    return pd.DataFrame(rows)
