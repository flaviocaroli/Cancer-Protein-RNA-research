import numpy as np
import pandas as pd
import itertools
from typing import Iterable, Sequence, Tuple, Optional

from sklearn.model_selection import KFold
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr, pearsonr
import xgboost as xgb

np.random.seed(42)

def _all_param_combinations(param_grid: dict) -> Sequence[dict]:
    keys, values = zip(*param_grid.items())
    return [dict(zip(keys, v)) for v in itertools.product(*values)]

def _clone_xgb_with_params(params: dict, seed: int = 42) -> xgb.XGBRegressor:
    base = dict(objective="reg:squarederror", n_jobs=-1, random_state=seed, verbosity=0)
    base.update(params or {})
    return xgb.XGBRegressor(**base)

def _compute_pvalues(X: pd.DataFrame, y: pd.Series, method: str) -> np.ndarray:
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
    return np.where(np.isfinite(pvals), pvals, 1.0)

def feature_selection_percentile(X: pd.DataFrame, y: pd.Series, percentile: float = 10.0, method: str = "pearson") -> pd.Index:
    pvals = _compute_pvalues(X, y, method)
    thresh = np.percentile(pvals, percentile)
    mask = pvals <= thresh
    if not np.any(mask):
        mask = pvals <= np.sort(pvals)[0]
    return X.columns[mask]

def _select_by_pvalue(X: pd.DataFrame, y: pd.Series, p_value: float, method: str) -> pd.Index:
    pvals = _compute_pvalues(X, y, method)
    mask = pvals <= p_value
    if not np.any(mask):
        k = max(1, int(0.01 * X.shape[1]))
        top_idx = np.argsort(pvals)[:k]
        return X.columns[top_idx]
    return X.columns[mask]

def _fast_prefilter_corr(X: pd.DataFrame, y: pd.Series, k: Optional[int]) -> pd.Index:
    """Keep top-k columns by absolute Pearson corr (vectorized & fast)."""
    if k is None or X.shape[1] <= k:
        return X.columns
    s = X.corrwith(y).fillna(0.0).abs()
    return s.nlargest(k).index

def _fit_and_score_cv(
    model: xgb.XGBRegressor,
    param_grid: dict,
    X_train: pd.DataFrame,
    y_train: pd.Series,
    feature_selector: str,
    select_mode: str,
    select_value: float,
    n_splits: int,
    preselect_k: Optional[int] = None,
    early_stopping_rounds: int = 50,
) -> Tuple[dict, pd.Index]:
    params_list = _all_param_combinations(param_grid)
    avg_mse = np.zeros(len(params_list), dtype=float)
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

    # Preselect once per full X_train for consistency across folds
    keep_full = _fast_prefilter_corr(X_train, y_train, preselect_k)
    X_train = X_train[keep_full]

    for fit_idx, val_idx in kf.split(X_train):
        X_fit, X_val = X_train.iloc[fit_idx], X_train.iloc[val_idx]
        y_fit, y_val = y_train.iloc[fit_idx], y_train.iloc[val_idx]

        if select_mode == "percentile":
            sel_features = feature_selection_percentile(X_fit, y_fit, float(select_value), feature_selector)
        elif select_mode == "pvalue":
            sel_features = _select_by_pvalue(X_fit, y_fit, float(select_value), feature_selector)
        else:
            raise ValueError("select_mode must be 'percentile' or 'pvalue'.")

        if len(sel_features) == 0:
            continue

        X_fit_sel = X_fit[sel_features]
        X_val_sel = X_val[sel_features]

        for k_idx, params in enumerate(params_list):
            est = _clone_xgb_with_params(params, seed=42)
            try:
                if early_stopping_rounds and early_stopping_rounds > 0:
                    est.fit(
                        X_fit_sel, y_fit,
                        eval_set=[(X_val_sel, y_val)],
                        early_stopping_rounds=int(early_stopping_rounds),
                        verbose=False,
                    )
                else:
                    est.fit(X_fit_sel, y_fit)
                pred = est.predict(X_val_sel)
                mse = mean_squared_error(y_val, pred)
                avg_mse[k_idx] += mse / n_splits
            except Exception:
                avg_mse[k_idx] += 1e6 / n_splits

    best_idx = int(np.argmin(avg_mse))
    best_params = params_list[best_idx]

    # Recompute features on full training data (with the same preselection universe)
    if select_mode == "percentile":
        best_features = feature_selection_percentile(X_train, y_train, float(select_value), feature_selector)
    else:
        best_features = _select_by_pvalue(X_train, y_train, float(select_value), feature_selector)

    return best_params, best_features

def _fit_final_and_predict(params: dict, X_train_sel: pd.DataFrame, y_train: pd.Series, X_test_sel: pd.DataFrame):
    est = _clone_xgb_with_params(params, seed=42)
    try:
        est.fit(X_train_sel, y_train)
        y_pred = est.predict(X_test_sel)
        imp = None
        if hasattr(est, "feature_importances_") and len(X_train_sel.columns) == len(est.feature_importances_):
            imp = dict(zip(X_train_sel.columns, est.feature_importances_.tolist()))
        return y_pred, imp
    except Exception:
        return np.array([np.nan]), None

def cv_pipeline_fixed_percentile(
    model, param_grid, X, y, protein_name,
    feature_selector="pearson", n_splits=3, n_samples=None, percentile=10.0,
    preselect_k: Optional[int] = 5000, early_stopping_rounds: int = 50,
):
    if protein_name not in y.columns:
        raise ValueError(f"Protein {protein_name} not found in y columns.")
    y_protein = y[protein_name]
    sample_list = X.index if n_samples is None else X.index[:n_samples]
    X = X.loc[sample_list]; y_protein = y_protein.loc[sample_list]
    results = []

    for test_sample in sample_list:
        X_train = X.drop(index=[test_sample]); y_train = y_protein.drop(index=[test_sample])
        X_test = X.loc[[test_sample]]; y_test = y_protein.loc[[test_sample]]

        best_params, best_features = _fit_and_score_cv(
            model, param_grid, X_train, y_train, feature_selector,
            "percentile", float(percentile), n_splits,
            preselect_k=preselect_k, early_stopping_rounds=early_stopping_rounds
        )
        Xtr = X_train[best_features]; Xte = X_test[best_features]
        y_pred, feat_imp = _fit_final_and_predict(best_params, Xtr, y_train, Xte)

        results.append({
            "protein": protein_name,
            "LOO sample": test_sample,
            "true": y_test.values[0],
            "predicted": y_pred[0] if isinstance(y_pred, np.ndarray) else y_pred,
            "params_dict": {**best_params, "feature_selector": feature_selector, "percentile": float(percentile)},
            "features_coef_dict": feat_imp,
            "percentile": float(percentile),
        })
    return pd.DataFrame(results)

def cv_pipeline_variable_p(
    model, param_grid, X, y, protein_name,
    feature_selector="pearson", n_splits=3, n_samples=None, p_values=None,
    preselect_k: Optional[int] = 5000, early_stopping_rounds: int = 50,
):
    if protein_name not in y.columns:
        raise ValueError(f"Protein {protein_name} not found in y columns.")
    if p_values is None:
        p_values = [0.1, 0.05, 0.01]

    y_protein = y[protein_name]
    sample_list = X.index if n_samples is None else X.index[:n_samples]
    X = X.loc[sample_list]; y_protein = y_protein.loc[sample_list]
    results = []

    for test_sample in sample_list:
        X_train = X.drop(index=[test_sample]); y_train = y_protein.drop(index=[test_sample])
        X_test = X.loc[[test_sample]]; y_test = y_protein.loc[[test_sample]]

        for pv in p_values:
            best_params, best_features = _fit_and_score_cv(
                model, param_grid, X_train, y_train, feature_selector,
                "pvalue", float(pv), n_splits,
                preselect_k=preselect_k, early_stopping_rounds=early_stopping_rounds
            )
            Xtr = X_train[best_features]; Xte = X_test[best_features]
            y_pred, feat_imp = _fit_final_and_predict(best_params, Xtr, y_train, Xte)
            results.append({
                "protein": protein_name,
                "LOO sample": test_sample,
                "true": y_test.values[0],
                "predicted": y_pred[0] if isinstance(y_pred, np.ndarray) else y_pred,
                "params_dict": {**best_params, "feature_selector": feature_selector, "p_value": float(pv)},
                "features_coef_dict": feat_imp,
                "p_value": float(pv),
            })
    return pd.DataFrame(results)

def cv_pipeline_fixed_p(
    model, param_grid, X, y, protein_name,
    feature_selector="pearson", n_splits=3, n_samples=None, p_value=0.05,
    preselect_k: Optional[int] = 5000, early_stopping_rounds: int = 50,
):
    if protein_name not in y.columns:
        raise ValueError(f"Protein {protein_name} not found in y columns.")
    y_protein = y[protein_name]
    sample_list = X.index if n_samples is None else X.index[:n_samples]
    X = X.loc[sample_list]; y_protein = y_protein.loc[sample_list]
    results = []

    for test_sample in sample_list:
        X_train = X.drop(index=[test_sample]); y_train = y_protein.drop(index=[test_sample])
        X_test = X.loc[[test_sample]]; y_test = y_protein.loc[[test_sample]]

        best_params, best_features = _fit_and_score_cv(
            model, param_grid, X_train, y_train, feature_selector,
            "pvalue", float(p_value), n_splits,
            preselect_k=preselect_k, early_stopping_rounds=early_stopping_rounds
        )
        Xtr = X_train[best_features]; Xte = X_test[best_features]
        y_pred, feat_imp = _fit_final_and_predict(best_params, Xtr, y_train, Xte)
        results.append({
            "protein": protein_name,
            "LOO sample": test_sample,
            "true": y_test.values[0],
            "predicted": y_pred[0] if isinstance(y_pred, np.ndarray) else y_pred,
            "params_dict": {**best_params, "feature_selector": feature_selector, "p_value": float(p_value)},
            "features_coef_dict": feat_imp,
            "p_value": float(p_value),
        })
    return pd.DataFrame(results)

def cv_pipeline_fixed_features(
    model, param_grid, X, y, protein_name, features: Iterable[str], n_splits=3, n_samples=None
):
    if protein_name not in y.columns:
        raise ValueError(f"Protein {protein_name} not found in y columns.")
    features = pd.Index(features)
    missing = features.difference(X.columns)
    if len(missing) > 0:
        raise ValueError(f"{len(missing)} features not found in X. Example: {missing[:5].tolist()}")

    y_protein = y[protein_name]
    sample_list = X.index if n_samples is None else X.index[:n_samples]
    X = X.loc[sample_list]; y_protein = y_protein.loc[sample_list]
    results = []
    params_list = _all_param_combinations(param_grid)

    for test_sample in sample_list:
        X_train = X.drop(index=[test_sample]); y_train = y_protein.drop(index=[test_sample])
        X_test = X.loc[[test_sample]]; y_test = y_protein.loc[[test_sample]]

        avg_mse = np.zeros(len(params_list), dtype=float)
        kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
        for fit_idx, val_idx in kf.split(X_train):
            X_fit, X_val = X_train.iloc[fit_idx][features], X_train.iloc[val_idx][features]
            y_fit, y_val = y_train.iloc[fit_idx], y_train.iloc[val_idx]
            for k_idx, params in enumerate(params_list):
                est = _clone_xgb_with_params(params, seed=42)
                try:
                    est.fit(X_fit, y_fit)
                    pred = est.predict(X_val)
                    avg_mse[k_idx] += mean_squared_error(y_val, pred) / n_splits
                except Exception:
                    avg_mse[k_idx] += 1e6 / n_splits

        best_idx = int(np.argmin(avg_mse))
        best_params = params_list[best_idx]
        y_pred, feat_imp = _fit_final_and_predict(best_params, X_train[features], y_train, X_test[features])

        results.append({
            "protein": protein_name,
            "LOO sample": test_sample,
            "true": y_test.values[0],
            "predicted": y_pred[0] if isinstance(y_pred, np.ndarray) else y_pred,
            "params_dict": {**best_params, "feature_selector": "fixed_features"},
            "features_coef_dict": feat_imp,
        })
    return pd.DataFrame(results)
