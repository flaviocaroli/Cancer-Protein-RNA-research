import os
import numpy as np
import pandas as pd
import itertools
from typing import Sequence, Tuple, Optional

from sklearn.model_selection import KFold
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr, pearsonr
import xgboost as xgb

np.random.seed(42)

# ----------------------
# Utilities
# ----------------------
def _all_param_combinations(param_grid: dict) -> Sequence[dict]:
    """Cartesian product of param grid (empty -> [{}])."""
    if not param_grid:
        return [dict()]
    keys, values = zip(*param_grid.items())
    return [dict(zip(keys, v)) for v in itertools.product(*values)]

def _sample_param_combos(param_grid: dict, max_param_combos: Optional[int]) -> Sequence[dict]:
    """Optionally sample a subset of param combos for speed."""
    combos = _all_param_combinations(param_grid)
    if max_param_combos is not None and len(combos) > max_param_combos:
        rng = np.random.default_rng(42)
        idx = rng.choice(len(combos), size=max_param_combos, replace=False)
        combos = [combos[i] for i in idx]
    return combos

def _clone_xgb_with_params(params: dict, seed: int = 42) -> xgb.XGBRegressor:
    # Make XGBoost honor Slurm CPU allocation:
    #   - run_xgboost.sh exports XGB_NJOBS=$CPUS_PER_TASK (default 1)
    n_jobs_env = int(os.environ.get("XGB_NJOBS", "1"))
    base = dict(objective="reg:squarederror", n_jobs=n_jobs_env, random_state=seed, verbosity=0)
    base.update(params or {})
    return xgb.XGBRegressor(**base)

def _compute_pvalues(X: pd.DataFrame, y: pd.Series, method: str) -> np.ndarray:
    """Return per-feature p-values for association with y."""
    if method == "anova":
        sel = SelectKBest(score_func=f_regression, k="all")
        sel.fit(X, y)
        pvals = sel.pvalues_
    elif method == "spearman":
        pvals = np.array([spearmanr(X.iloc[:, i], y)[1] for i in range(X.shape[1])])
    elif method == "pearson":
        pvals = np.array([pearsonr(X.iloc[:, i], y)[1] for i in range(X.shape[1])])
    else:
        raise ValueError(f"Unknown feature selection method: {method}")
    # Replace non-finite with 1.0 so they rank last
    return np.where(np.isfinite(pvals), pvals, 1.0)

def _select_topk_by_pvalue(X: pd.DataFrame, y: pd.Series, method: str, k: int) -> pd.Index:
    """Keep top-k features with smallest p-values (guarantee at least 1)."""
    k_eff = max(1, min(int(k), X.shape[1]))
    pvals = _compute_pvalues(X, y, method)
    top_idx = np.argsort(pvals)[:k_eff]
    return X.columns[top_idx]

def _fast_prefilter_corr(X: pd.DataFrame, y: pd.Series, k: Optional[int]) -> pd.Index:
    """Keep top-k columns by absolute Pearson corr (vectorized & fast)."""
    if k is None or X.shape[1] <= k:
        return X.columns
    s = X.corrwith(y).fillna(0.0).abs()
    return s.nlargest(k).index

def _to_f32(df: pd.DataFrame) -> np.ndarray:
    return df.to_numpy(dtype=np.float32, copy=False)

def _to_f32_series(s: pd.Series) -> np.ndarray:
    return s.to_numpy(dtype=np.float32, copy=False)

# ----------------------
# Inner CV: pick params & features (no leakage)
# ----------------------
def _fit_and_score_cv(
    param_grid: dict,
    X_train: pd.DataFrame,
    y_train: pd.Series,
    feature_selector: str,
    select_k: int,
    n_splits: int,
    preselect_k: Optional[int] = None,
    early_stopping_rounds: int = 50,
    max_param_combos: Optional[int] = None,
) -> Tuple[dict, pd.Index]:
    """
    Inner KFold CV to choose best params. Feature selection happens *inside each fold*:
      1) optional fast prefilter by |corr| (top preselect_k)
      2) select top-k by p-value
    After picking best params, recompute features on the *full* outer training set.
    If select_k <= 0, disables all feature selection (keeps all columns).
    """
    params_list = _sample_param_combos(param_grid, max_param_combos)
    avg_mse = np.zeros(len(params_list), dtype=float)
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

    for fit_idx, val_idx in kf.split(X_train):
        X_fit_full, X_val_full = X_train.iloc[fit_idx], X_train.iloc[val_idx]
        y_fit, y_val = y_train.iloc[fit_idx], y_train.iloc[val_idx]

        if select_k <= 0:
            X_fit, X_val = X_fit_full, X_val_full
            sel_features = X_fit.columns
        else:
            keep = _fast_prefilter_corr(X_fit_full, y_fit, preselect_k)
            X_fit, X_val = X_fit_full[keep], X_val_full[keep]
            sel_features = _select_topk_by_pvalue(X_fit, y_fit, feature_selector, select_k)

        X_fit_sel, X_val_sel = X_fit[sel_features], X_val[sel_features]

        X_fit_np = _to_f32(X_fit_sel)
        X_val_np = _to_f32(X_val_sel)
        y_fit_np = _to_f32_series(y_fit)
        y_val_np = _to_f32_series(y_val)

        for k_idx, params in enumerate(params_list):
            est = _clone_xgb_with_params(params, seed=42)
            try:
                if early_stopping_rounds and early_stopping_rounds > 0:
                    est.fit(
                        X_fit_np, y_fit_np,
                        eval_set=[(X_val_np, y_val_np)],
                        early_stopping_rounds=int(early_stopping_rounds),
                        verbose=False,
                    )
                else:
                    est.fit(X_fit_np, y_fit_np)
                pred = est.predict(X_val_np)
                mse = mean_squared_error(y_val_np, pred)
                avg_mse[k_idx] += mse / n_splits
            except Exception:
                avg_mse[k_idx] += 1e6 / n_splits

    best_idx = int(np.argmin(avg_mse))
    best_params = params_list[best_idx] if len(params_list) > 0 else {}

    if select_k <= 0:
        best_features = X_train.columns
    else:
        keep_full = _fast_prefilter_corr(X_train, y_train, preselect_k)
        X_train_pref = X_train[keep_full]
        best_features = _select_topk_by_pvalue(X_train_pref, y_train, feature_selector, select_k)

    return best_params, best_features

# ----------------------
# Final fit and predict on held-out sample
# ----------------------
def _fit_final_and_predict(params: dict, X_train_sel: pd.DataFrame, y_train: pd.Series, X_test_sel: pd.DataFrame):
    est = _clone_xgb_with_params(params, seed=42)
    try:
        Xtr = _to_f32(X_train_sel); ytr = _to_f32_series(y_train)
        Xte = _to_f32(X_test_sel)
        est.fit(Xtr, ytr)
        y_pred = est.predict(Xte)

        imp = None
        try:
            booster = est.get_booster()
            gain = booster.get_score(importance_type="gain")
            fmap = {f"f{i}": col for i, col in enumerate(X_train_sel.columns)}
            imp = {fmap.get(k, k): float(v) for k, v in gain.items()}
        except Exception:
            pass

        return y_pred, imp
    except Exception:
        return np.array([np.nan]), None

# ----------------------
# Single, unified LOO pipeline (outer LOO only)
# ----------------------
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
    early_stopping_rounds: int = 50,
    minimal_save: bool = False,
    max_param_combos: Optional[int] = None,
) -> pd.DataFrame:
    """
    Outer loop: Leave-One-Out CV across samples.
    Inner loop: KFold for hyperparameter tuning + in-fold feature selection.
    Returns a tidy DataFrame with one row per held-out sample.
    """
    if protein_name not in y.columns:
        raise ValueError(f"Protein {protein_name} not found in y columns.")

    y_protein = y[protein_name]
    sample_list = X.index if n_samples is None else X.index[:n_samples]
    X = X.loc[sample_list]
    y_protein = y_protein.loc[sample_list]

    results = []

    for test_sample in sample_list:
        X_train = X.drop(index=[test_sample])
        y_train = y_protein.drop(index=[test_sample])
        X_test = X.loc[[test_sample]]
        y_test = y_protein.loc[[test_sample]]

        best_params, best_features = _fit_and_score_cv(
            param_grid=param_grid,
            X_train=X_train,
            y_train=y_train,
            feature_selector=feature_selector,
            select_k=int(select_k),
            n_splits=int(n_splits),
            preselect_k=preselect_k,
            early_stopping_rounds=int(early_stopping_rounds),
            max_param_combos=max_param_combos,
        )

        Xtr = X_train[best_features]
        Xte = X_test[best_features]
        y_pred, feat_imp = _fit_final_and_predict(best_params, Xtr, y_train, Xte)

        pred_val = float(y_pred[0]) if isinstance(y_pred, np.ndarray) else float(y_pred)

        if minimal_save:
            results.append({
                "protein": protein_name,
                "LOO sample": test_sample,
                "true": float(y_test.values[0]),
                "predicted": pred_val,
            })
        else:
            results.append({
                "protein": protein_name,
                "LOO sample": test_sample,
                "true": float(y_test.values[0]),
                "predicted": pred_val,
                "params_dict": {
                    **best_params,
                    "feature_selector": feature_selector,
                    "select_k": int(select_k),
                    "preselect_k": None if preselect_k is None else int(preselect_k),
                    "early_stopping_rounds": int(early_stopping_rounds),
                    "inner_cv_folds": int(n_splits),
                },
                "features_coef_dict": feat_imp,
                "features_selected": list(map(str, best_features))
            })

    return pd.DataFrame(results)
