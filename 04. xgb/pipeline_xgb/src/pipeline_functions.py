import itertools
from typing import Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import KFold
import xgboost as xgb

np.random.seed(42)

# ----------------------
# Small utilities
# ----------------------
def _all_param_combinations(param_grid: dict) -> Sequence[dict]:
    """Cartesian product of a param grid (empty -> [{}])."""
    if not param_grid:
        return [dict()]
    keys, vals = zip(*param_grid.items())
    return [dict(zip(keys, combo)) for combo in itertools.product(*vals)]

def _sample_param_combos(param_grid: dict, max_param_combos: Optional[int]) -> Sequence[dict]:
    """Optionally sample a subset of param combos for speed."""
    combos = _all_param_combinations(param_grid)
    if max_param_combos is not None and len(combos) > max_param_combos:
        rng = np.random.default_rng(42)
        picks = rng.choice(len(combos), size=max_param_combos, replace=False)
        combos = [combos[i] for i in picks]
    return combos

def _clone_xgb_with_params(params: dict, seed: int = 42, use_gpu: bool = False) -> xgb.XGBRegressor:
    base = dict(objective="reg:squarederror", n_jobs=-1, random_state=seed, verbosity=0)
    # Enable GPU unless user explicitly overrides in param_grid
    if use_gpu and "tree_method" not in params and "device" not in params:
        base.update(dict(tree_method="gpu_hist", predictor="gpu_predictor"))
    base.update(params or {})
    return xgb.XGBRegressor(**base)

def _to_f32(df: pd.DataFrame) -> np.ndarray:
    return df.to_numpy(dtype=np.float32, copy=False)

def _to_f32_series(s: pd.Series) -> np.ndarray:
    return s.to_numpy(dtype=np.float32, copy=False)

# ----------------------
# Spearman prefilter (run ONCE per outer split)
# ----------------------
def _spearman_abs_corr_topk(X: pd.DataFrame, y: pd.Series, k: Optional[int]) -> pd.Index:
    """
    Keep top-k columns by absolute Spearman correlation with y.
    Vectorized via rank transform: Spearman = Pearson on ranks.
    """
    if k is None or k <= 0 or X.shape[1] <= k:
        return X.columns
    Xr = X.rank(axis=0, method="average")
    yr = y.rank(method="average")
    s = Xr.corrwith(yr).abs().fillna(0.0)
    return s.nlargest(int(k)).index

# ----------------------
# Inner CV: tune hyperparams ONLY (no feature selection here)
# ----------------------
def _fit_and_score_cv(
    param_grid: dict,
    X_train: pd.DataFrame,
    y_train: pd.Series,
    n_splits: int,
    early_stopping_rounds: int = 50,
    max_param_combos: Optional[int] = None,
    use_gpu: bool = False,
) -> dict:
    """
    KFold CV for hyperparameter tuning.
    Assumes columns of X_train are already fixed (prefiltered).
    """
    params_list = _sample_param_combos(param_grid, max_param_combos)
    if len(params_list) == 0:
        return {}

    avg_mse = np.zeros(len(params_list), dtype=float)
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

    for tr_idx, va_idx in kf.split(X_train):
        Xtr_df, Xva_df = X_train.iloc[tr_idx], X_train.iloc[va_idx]
        ytr_s, yva_s = y_train.iloc[tr_idx], y_train.iloc[va_idx]

        Xtr, Xva = _to_f32(Xtr_df), _to_f32(Xva_df)
        ytr, yva = _to_f32_series(ytr_s), _to_f32_series(yva_s)

        for k_idx, params in enumerate(params_list):
            est = _clone_xgb_with_params(params, seed=42, use_gpu=use_gpu)
            try:
                if early_stopping_rounds and early_stopping_rounds > 0:
                    est.fit(Xtr, ytr, eval_set=[(Xva, yva)],
                            early_stopping_rounds=int(early_stopping_rounds),
                            verbose=False)
                else:
                    est.fit(Xtr, ytr)
                pred = est.predict(Xva)
                avg_mse[k_idx] += mean_squared_error(yva, pred) / n_splits
            except Exception:
                # Penalize failures heavily to avoid selecting broken configs
                avg_mse[k_idx] += 1e6 / n_splits

    best_idx = int(np.argmin(avg_mse))
    return params_list[best_idx]

# ----------------------
# Final fit and predict on held-out sample
# ----------------------
def _fit_final_and_predict(
    params: dict,
    X_train_sel: pd.DataFrame,
    y_train: pd.Series,
    X_test_sel: pd.DataFrame,
    use_gpu: bool = False,
):
    est = _clone_xgb_with_params(params, seed=42, use_gpu=use_gpu)
    try:
        est.fit(_to_f32(X_train_sel), _to_f32_series(y_train))
        y_pred = est.predict(_to_f32(X_test_sel))

        # Importances via booster (gain)
        feat_imp = None
        try:
            gain = est.get_booster().get_score(importance_type="gain")  # keys like "f0"
            fmap = {f"f{i}": col for i, col in enumerate(X_train_sel.columns)}
            feat_imp = {fmap.get(k, k): float(v) for k, v in gain.items()}
        except Exception:
            pass

        return y_pred, feat_imp
    except Exception:
        return np.array([np.nan]), None

# ----------------------
# Single, unified LOO pipeline (outer LOO only)
# ----------------------
def cv_pipeline_loo(
    model,  # kept for signature compatibility; not used
    param_grid: dict,
    X: pd.DataFrame,
    y: pd.DataFrame,
    protein_name: str,
    *,
    feature_selector: str = "spearman",   # ignored; always Spearman prefilter
    n_splits: int = 3,
    n_samples: Optional[int] = None,
    select_k: int = 500,                  # number of features to keep via Spearman prefilter
    preselect_k: Optional[int] = 5000,    # ignored; kept for compatibility
    early_stopping_rounds: int = 50,
    minimal_save: bool = False,
    max_param_combos: Optional[int] = None,
    use_gpu: bool = False,                # set True to use GPU
) -> pd.DataFrame:
    """
    Outer loop: Leave-One-Out CV across samples.
    Inside each outer split:
      1) Spearman prefilter ONCE on the outer training set (top-k by |rho|).
      2) Inner KFold CV tunes hyperparams using the fixed feature set.
      3) Train on full outer training set with best params; predict held-out sample.

    Notes:
      - No feature selection happens inside inner CV.
      - 'feature_selector' and 'preselect_k' are ignored (kept for API compatibility).
    """
    if protein_name not in y.columns:
        raise ValueError(f"Protein {protein_name} not found in y columns.")

    y_protein = y[protein_name]
    sample_idx = X.index if n_samples is None else X.index[:n_samples]
    X, y_protein = X.loc[sample_idx], y_protein.loc[sample_idx]

    rows = []
    for test_id in sample_idx:
        # Outer LOO split
        X_train = X.drop(index=[test_id])
        y_train = y_protein.drop(index=[test_id])
        X_test = X.loc[[test_id]]
        y_test = y_protein.loc[[test_id]]

        # (1) Spearman prefilter ONCE on outer training set
        kept_features = _spearman_abs_corr_topk(X_train, y_train, select_k)
        Xtr_sel = X_train[kept_features]
        Xte_sel = X_test[kept_features]

        # (2) Inner CV: choose params ONLY (no FS inside folds)
        best_params = _fit_and_score_cv(
            param_grid=param_grid,
            X_train=Xtr_sel,
            y_train=y_train,
            n_splits=int(n_splits),
            early_stopping_rounds=int(early_stopping_rounds),
            max_param_combos=max_param_combos,
            use_gpu=use_gpu,
        )

        # (3) Final fit on full outer train with chosen params; predict held-out
        y_pred, feat_imp = _fit_final_and_predict(
            best_params,
            Xtr_sel,
            y_train,
            Xte_sel,
            use_gpu=use_gpu,
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
                    "preselect_k": None,  # ignored
                    "early_stopping_rounds": int(early_stopping_rounds),
                    "inner_cv_folds": int(n_splits),
                },
                "features_coef_dict": feat_imp,                      # may be None
                "features_selected": list(map(str, kept_features)),   # explicit list
            })

    return pd.DataFrame(rows)
