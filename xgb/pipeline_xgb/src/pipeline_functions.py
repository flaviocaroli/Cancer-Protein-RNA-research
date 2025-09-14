import numpy as np
import pandas as pd
import itertools
from typing import Iterable, Optional, Sequence, Tuple

from sklearn.model_selection import KFold
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr, pearsonr
import xgboost as xgb

np.random.seed(42)


def _all_param_combinations(param_grid: dict) -> Sequence[dict]:
    """
    Expand a scikit-learn–style parameter grid into a list of dictionaries.

    Parameters
    ----------
    param_grid : dict
        Dictionary mapping hyperparameter names to iterables of candidate values,
        e.g. {"max_depth": [3, 6], "learning_rate": [0.1, 0.3]}.

    Returns
    -------
    list of dict
        Each dict is one unique combination across the Cartesian product of
        all provided hyperparameter values. The order is deterministic.
    """
    keys, values = zip(*param_grid.items())
    return [dict(zip(keys, v)) for v in itertools.product(*values)]


def _clone_xgb_with_params(params: dict, seed: int = 42) -> xgb.XGBRegressor:
    """
    Instantiate a fresh XGBRegressor with sane defaults plus provided params.

    Parameters
    ----------
    params : dict
        Hyperparameters to pass into `xgboost.XGBRegressor`. These override the
        built-in base defaults defined here.
    seed : int, default=42
        Random seed for reproducibility (mapped to `random_state`).

    Returns
    -------
    xgboost.XGBRegressor
        New, unfitted estimator with combined parameters:
        - objective="reg:squarederror"
        - n_jobs=-1
        - random_state=seed
        - verbosity=0
        plus any keys in `params`.
    """
    base = dict(objective="reg:squarederror", n_jobs=-1, random_state=seed, verbosity=0)
    base.update(params or {})
    return xgb.XGBRegressor(**base)


def _compute_pvalues(X: pd.DataFrame, y: pd.Series, method: str) -> np.ndarray:
    """
    Compute per-feature p-values for association with the target.

    Parameters
    ----------
    X : DataFrame, shape (n_samples, n_features)
        Feature matrix.
    y : Series or array-like, shape (n_samples,)
        Target values.
    method : {"anova", "spearman", "pearson"}
        Statistical test used to score features:
        - "anova": f_regression (linear relationship; same as scikit-learn)
        - "spearman": Spearman rank correlation
        - "pearson": Pearson correlation

    Returns
    -------
    np.ndarray, shape (n_features,)
        P-values aligned to `X.columns`. Any non-finite values are replaced by
        1.0 (worst).
    """
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


def feature_selection_percentile(
    X: pd.DataFrame,
    y: pd.Series,
    percentile: float = 10.0,
    method: str = "pearson",
) -> pd.Index:
    """
    Select features by taking the best (lowest p-value) `percentile` percent.

    Parameters
    ----------
    X : DataFrame, shape (n_samples, n_features)
        Feature matrix.
    y : Series, shape (n_samples,)
        Target values.
    percentile : float, default=10.0
        Percentage of features to keep based on lowest p-values (0–100).
    method : {"pearson", "spearman", "anova"}, default="pearson"
        Statistical test used to score features.

    Returns
    -------
    pd.Index
        Names of selected features. If no feature passes due to ties/NaNs,
        fall back to keeping the single best feature.
    """
    pvals = _compute_pvalues(X, y, method)
    thresh = np.percentile(pvals, percentile)
    mask = pvals <= thresh
    if not np.any(mask):
        mask = pvals <= np.sort(pvals)[0]
    return X.columns[mask]


def _select_by_pvalue(X: pd.DataFrame, y: pd.Series, p_value: float, method: str) -> pd.Index:
    """
    Select features whose p-values are <= a fixed threshold.

    Parameters
    ----------
    X : DataFrame, shape (n_samples, n_features)
        Feature matrix.
    y : Series, shape (n_samples,)
        Target values.
    p_value : float
        P-value threshold (e.g., 0.05).
    method : {"pearson", "spearman", "anova"}
        Statistical test used to score features.

    Returns
    -------
    pd.Index
        Names of selected features. If none pass, fall back to the top 1% by
        lowest p-values (at least one feature).
    """
    pvals = _compute_pvalues(X, y, method)
    mask = pvals <= p_value
    if not np.any(mask):
        k = max(1, int(0.01 * X.shape[1]))
        top_idx = np.argsort(pvals)[:k]
        return X.columns[top_idx]
    return X.columns[mask]


def _fit_and_score_cv(
    model: xgb.XGBRegressor,
    param_grid: dict,
    X_train: pd.DataFrame,
    y_train: pd.Series,
    feature_selector: str,
    select_mode: str,
    select_value: float,
    n_splits: int,
) -> Tuple[dict, pd.Index]:
    """
    Inner CV loop to choose best hyperparameters and derive final feature set.

    In each fold:
      1) Perform feature selection on the fold's training split using either a
         percentile cutoff or a fixed p-value threshold.
      2) Fit an XGBoost model for every parameter combination on selected
         features and evaluate MSE on the validation split.
    After CV:
      • Pick the best parameter set by lowest average MSE across folds.
      • Recompute selected features on the *full* X_train using the chosen
        thresholding mode and value. Return these as the final feature set.

    Parameters
    ----------
    model : xgboost.XGBRegressor
        Prototype estimator (only used for its class/type).
    param_grid : dict
        Hyperparameter grid to search. See `_all_param_combinations`.
    X_train : DataFrame, shape (n_samples_train, n_features)
        Training data (all LOO training samples).
    y_train : Series, shape (n_samples_train,)
        Training targets.
    feature_selector : {"pearson", "spearman", "anova"}
        Statistical method for scoring features.
    select_mode : {"percentile", "pvalue"}
        Selection strategy.
    select_value : float
        If `percentile`, interpreted as a percentage (0–100);
        if `pvalue`, interpreted as a threshold (e.g., 0.05).
    n_splits : int
        Number of KFold splits for inner cross-validation.

    Returns
    -------
    (best_params, best_features) : (dict, pd.Index)
        The best-performing hyperparameter dict and the final feature set
        recomputed on full X_train with the chosen selection rule.

    Notes
    -----
    • Any training failure for a param combo is penalized with a large MSE to
      ensure it is not selected.
    """
    params_list = _all_param_combinations(param_grid)
    avg_mse = np.zeros(len(params_list), dtype=float)
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

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
                est.fit(X_fit_sel, y_fit)
                pred = est.predict(X_val_sel)
                mse = mean_squared_error(y_val, pred)
                avg_mse[k_idx] += mse / n_splits
            except Exception:
                avg_mse[k_idx] += 1e6 / n_splits

    best_idx = int(np.argmin(avg_mse))
    best_params = params_list[best_idx]

    # Recompute features on full training data
    if select_mode == "percentile":
        best_features = feature_selection_percentile(X_train, y_train, float(select_value), feature_selector)
    else:
        best_features = _select_by_pvalue(X_train, y_train, float(select_value), feature_selector)

    return best_params, best_features


def _fit_final_and_predict(
    params: dict,
    X_train_sel: pd.DataFrame,
    y_train: pd.Series,
    X_test_sel: pd.DataFrame
):
    """
    Fit final XGBoost on selected features and predict the held-out sample.

    Parameters
    ----------
    params : dict
        Best hyperparameters to instantiate the model.
    X_train_sel : DataFrame, shape (n_train, n_features_selected)
        Training matrix restricted to selected features.
    y_train : Series, shape (n_train,)
        Training targets.
    X_test_sel : DataFrame, shape (1, n_features_selected)
        Single test sample (LOO) restricted to selected features.

    Returns
    -------
    (y_pred, feature_importances) : (np.ndarray, dict or None)
        - y_pred : array of shape (1,) with the predicted value for the test sample.
        - feature_importances : dict mapping feature -> importance if available,
          else None.

    Notes
    -----
    Uses `_clone_xgb_with_params` to ensure a fresh estimator.
    """
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
    model,
    param_grid,
    X,
    y,
    protein_name,
    feature_selector="pearson",
    n_splits=5,
    n_samples=None,
    percentile=10.0
):
    """
    Leave-one-out (LOO) pipeline using a fixed percentile for feature selection.

    For each sample:
      • Hold out the sample.
      • Inner CV selects the best hyperparameters while performing feature
        selection at a fixed percentile (lowest p-values).
      • Refit on full training set (selected features) and predict the held-out
        sample.

    Parameters
    ----------
    model : xgboost.XGBRegressor
        Prototype model (its class is used to instantiate per-try estimators).
    param_grid : dict
        Hyperparameter grid to search.
    X : DataFrame, shape (n_samples, n_features)
        Transcriptomics (or other) feature matrix.
    y : DataFrame, shape (n_samples, n_proteins)
        Proteomics (or other targets). Must contain `protein_name` column.
    protein_name : str
        The protein to model (target column in `y`).
    feature_selector : {"pearson", "spearman", "anova"}, default="pearson"
        Method to compute per-feature p-values.
    n_splits : int, default=5
        Number of inner CV folds.
    n_samples : int or None, default=None
        If set, restrict to the first `n_samples` rows for quick runs.
    percentile : float, default=10.0
        Percent of features to keep (lowest p-values).

    Returns
    -------
    DataFrame
        One row per LOO test sample with:
        ["protein", "LOO sample", "true", "predicted", "params_dict",
         "features_coef_dict", "percentile"].
        `features_coef_dict` holds feature importances (if available).
    """
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
            model, param_grid, X_train, y_train, feature_selector, "percentile", float(percentile), n_splits
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
    model,
    param_grid,
    X,
    y,
    protein_name,
    feature_selector="pearson",
    n_splits=5,
    n_samples=None,
    p_values=None
):
    """
    Leave-one-out (LOO) pipeline sweeping multiple p-value thresholds.

    For each sample and each p-value in `p_values`:
      • Hold out the sample.
      • Inner CV selects hyperparameters with feature selection done using a
        fixed p-value threshold.
      • Refit on full training set with the chosen params and selected features,
        predict held-out sample, and record results.

    Parameters
    ----------
    model : xgboost.XGBRegressor
        Prototype model (its class is used to instantiate per-try estimators).
    param_grid : dict
        Hyperparameter grid to search.
    X : DataFrame, shape (n_samples, n_features)
        Feature matrix.
    y : DataFrame, shape (n_samples, n_targets)
        Target matrix containing `protein_name`.
    protein_name : str
        Column name in `y` to model.
    feature_selector : {"pearson", "spearman", "anova"}, default="pearson"
        Method to compute per-feature p-values.
    n_splits : int, default=5
        Number of inner CV folds.
    n_samples : int or None, default=None
        If set, restrict to the first `n_samples` rows for quick runs.
    p_values : iterable of float or None, default=None
        List of p-value thresholds to evaluate. Defaults to [0.1, 0.05, 0.01].

    Returns
    -------
    DataFrame
        One row per (LOO sample, p-value) with the same columns as the fixed
        percentile pipeline, plus the `p_value` used.
    """
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
                model, param_grid, X_train, y_train, feature_selector, "pvalue", float(pv), n_splits
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
    model,
    param_grid,
    X,
    y,
    protein_name,
    feature_selector="pearson",
    n_splits=5,
    n_samples=None,
    p_value=0.05
):
    """
    Leave-one-out (LOO) pipeline using a single fixed p-value threshold.

    This is a convenience specialization of the variable-p pipeline.

    Parameters
    ----------
    model : xgboost.XGBRegressor
        Prototype model.
    param_grid : dict
        Hyperparameter grid to search.
    X : DataFrame, shape (n_samples, n_features)
        Feature matrix.
    y : DataFrame, shape (n_samples, n_targets)
        Target matrix containing `protein_name`.
    protein_name : str
        Column name in `y` to model.
    feature_selector : {"pearson", "spearman", "anova"}, default="pearson"
        Method to compute per-feature p-values.
    n_splits : int, default=5
        Number of inner CV folds.
    n_samples : int or None, default=None
        If set, restrict to the first `n_samples` rows for quick runs.
    p_value : float, default=0.05
        P-value threshold used for feature selection.

    Returns
    -------
    DataFrame
        One row per LOO test sample with selected hyperparameters, prediction,
        and metadata (including the fixed `p_value` used).
    """
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
            model, param_grid, X_train, y_train, feature_selector, "pvalue", float(p_value), n_splits
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
    model,
    param_grid,
    X,
    y,
    protein_name,
    features: Iterable[str],
    n_splits=5,
    n_samples=None
):
    """
    Leave-one-out (LOO) pipeline using a pre-defined fixed set of features.

    Feature selection is *not* performed inside CV. Instead, the provided
    `features` list is used throughout. Inner CV searches only over model
    hyperparameters.

    Parameters
    ----------
    model : xgboost.XGBRegressor
        Prototype model.
    param_grid : dict
        Hyperparameter grid to search.
    X : DataFrame, shape (n_samples, n_features)
        Feature matrix.
    y : DataFrame, shape (n_samples, n_targets)
        Target matrix containing `protein_name`.
    protein_name : str
        Column name in `y` to model.
    features : iterable of str
        Fixed feature names to use (must be present in `X.columns`).
    n_splits : int, default=5
        Number of inner CV folds.
    n_samples : int or None, default=None
        If set, restrict to the first `n_samples` rows for quick runs.

    Returns
    -------
    DataFrame
        One row per LOO test sample with chosen hyperparameters, prediction,
        and feature importances (if available). No feature-selection metadata
        is attached other than the string "fixed_features" in `params_dict`.

    Raises
    ------
    ValueError
        If any requested `features` are missing from `X.columns`.
    """
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
