import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.feature_selection import SelectPercentile, f_regression, SelectKBest
from sklearn.metrics import mean_squared_error
from sklearn.linear_model import ElasticNet
from scipy.stats import spearmanr, pearsonr
import itertools

# Fix seed globally for NumPy
np.random.seed(42)


# -----------------------------
# Utilities
# -----------------------------

def _safe_corr_pvalue(x: np.ndarray, y: np.ndarray, method: str) -> float:
    """Return a p-value for correlation; fall back to 1.0 on failure/NaN."""
    try:
        if method == 'spearman':
            _, p = spearmanr(x, y)
        else:
            _, p = pearsonr(x, y)
        if np.isnan(p):
            return 1.0
        return float(p)
    except Exception:
        return 1.0


def _pvalues_series(X_feat_by_sample: pd.DataFrame, y: pd.Series, method: str) -> pd.Series:
    """Compute per-feature p-values (index: feature names).
    X_feat_by_sample is expected to be features x samples (rows=features, cols=samples).
    y is aligned to the columns (samples).
    """
    pv = [
        _safe_corr_pvalue(X_feat_by_sample.loc[f].values, y.values, method)
        for f in X_feat_by_sample.index
    ]
    return pd.Series(pv, index=X_feat_by_sample.index, dtype=float)


def _select_features_percentile(X_feat_by_sample: pd.DataFrame,
                                y: pd.Series,
                                feature_selector: str,
                                percentile: int) -> list:
    """Return a list of selected features given a percentile and selector."""
    if feature_selector == 'anova':
        # anova requires shape (samples x features)
        selector = SelectPercentile(score_func=f_regression, percentile=int(percentile))
        selector.fit(X_feat_by_sample.T, y.values)
        mask = selector.get_support()
        return list(X_feat_by_sample.index[mask])

    # Correlation-based (pearson/spearman): select the lowest p-values by percentile
    method = 'spearman' if feature_selector == 'spearman' else 'pearson'
    pv = _pvalues_series(X_feat_by_sample, y, method)
    pv = pv.replace([np.inf, -np.inf], 1.0).fillna(1.0)
    # threshold at the given percentile of p-values (lower is better)
    thresh = np.percentile(pv.values, percentile)
    selected = pv.index[pv.values <= thresh].tolist()

    # Guarantee at least one feature (pick the best available)
    if len(selected) == 0 and len(pv) > 0:
        selected = [pv.idxmin()]
    return selected


def _select_features_pvalue(X_feat_by_sample: pd.DataFrame,
                            y: pd.Series,
                            feature_selector: str,
                            p_value: float) -> list:
    """Return a list of selected features given a p-value cutoff and selector."""
    if feature_selector == 'anova':
        selector = SelectKBest(score_func=f_regression, k='all')
        selector.fit(X_feat_by_sample.T, y.values)
        pv = pd.Series(selector.pvalues_, index=X_feat_by_sample.index, dtype=float)
    else:
        method = 'spearman' if feature_selector == 'spearman' else 'pearson'
        pv = _pvalues_series(X_feat_by_sample, y, method)

    pv = pv.replace([np.inf, -np.inf], 1.0).fillna(1.0)
    selected = pv.index[pv.values <= p_value].tolist()

    if len(selected) == 0 and len(pv) > 0:
        selected = [pv.idxmin()]
    return selected


# -----------------------------
# Feature selection (saved per LOO sample)
# -----------------------------

def feature_selection_percentile(X: pd.DataFrame,
                                 y: pd.Series,
                                 percentile: int = 10) -> pd.DataFrame:
    """
    For each LOO sample, select features using the remaining samples only (Pearson by default).
    Returns a DataFrame indexed by sample with one column 'features' (list of feature names).

    Parameters
    ----------
    X : DataFrame
        Samples x Features (rows are samples, columns are features).
    y : Series
        Target aligned to X.index.
    percentile : int
        Percent of features to keep (lowest p-values / highest F-scores).
    """
    sample_list = X.index
    records = []

    for sample in sample_list:
        X_train = X.drop(index=[sample])
        y_train = y.drop(index=[sample])

        # features x samples
        X_train_fx = X_train.T
        selected = _select_features_percentile(X_train_fx, y_train, 'pearson', percentile)
        records.append({'LOO sample': sample, 'features': selected})

    return pd.DataFrame.from_records(records).set_index('LOO sample')


# -----------------------------
# Cross-validation pipelines (ElasticNet only)
# -----------------------------

def cv_pipeline_fixed_percentile(model,
                                 param_grid: dict,
                                 X: pd.DataFrame,
                                 y: pd.DataFrame,
                                 protein_name: str,
                                 feature_selector: str = 'pearson',
                                 n_splits: int = 5,
                                 n_samples=None,
                                 percentile: int = 10) -> pd.DataFrame:
    """
    LOO CV for a single protein with a fixed percentile for feature selection.
    Model: ElasticNet only. `model` is ignored unless it's an ElasticNet instance.
    """
    if protein_name not in y.columns:
        raise ValueError(f"Protein {protein_name} not found in y columns.")

    if (model is not None) and not isinstance(model, ElasticNet):
        raise ValueError("Only ElasticNet is supported. Pass an ElasticNet instance or None.")

    y_protein = y[protein_name]
    sample_list = X.index if n_samples is None else X.index[:n_samples]
    X = X.loc[sample_list]
    y_protein = y_protein.loc[sample_list]

    results = []

    keys, values = zip(*param_grid.items()) if len(param_grid) else ([], [])
    param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)] if len(keys) else [{}]

    for sample in sample_list:
        X_train = X.drop(index=[sample])
        X_test = X.loc[[sample]]
        y_train = y_protein.drop(index=[sample])
        y_test = y_protein.loc[[sample]]

        columns_train = X_train.index
        avg_mse = np.zeros(len(param_combinations), dtype=float)

        # Inner CV for hyperparameters
        kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

        for fit_idx, val_idx in kf.split(columns_train):
            fit_samples = columns_train[fit_idx]
            val_samples = columns_train[val_idx]

            X_fit = X_train.loc[fit_samples].T  # features x samples
            y_fit = y_train.loc[fit_samples]
            X_val = X_train.loc[val_samples].T
            y_val = y_train.loc[val_samples]

            # Select features on the fold's training data
            selected_features = _select_features_percentile(X_fit, y_fit, feature_selector, percentile)
            if len(selected_features) == 0:
                continue

            X_fit_selected = X_fit.loc[selected_features]
            X_val_selected = X_val.loc[selected_features]

            for k_idx, params in enumerate(param_combinations):
                estimator = ElasticNet(max_iter=10000, random_state=42, **params)
                try:
                    estimator.fit(X_fit_selected.T, y_fit.values)
                    y_pred_val = estimator.predict(X_val_selected.T)
                    mse = mean_squared_error(y_val.values, y_pred_val)
                    avg_mse[k_idx] += mse / n_splits
                except Exception:
                    avg_mse[k_idx] += 1e6 / n_splits

        # Best params
        min_mse_idx = int(np.argmin(avg_mse)) if np.isfinite(avg_mse).any() else 0
        best_params = {**param_combinations[min_mse_idx],
                       'percentile': percentile,
                       'feature_selector': feature_selector}

        # Re-select features on the full training set (post inner-CV)
        X_train_fx = X_train.T  # features x samples
        selected_features = _select_features_percentile(X_train_fx, y_train, feature_selector, percentile)

        # Final train and test
        features_coef = None
        y_pred = np.nan

        if len(selected_features) > 0:
            X_train_selected = X_train_fx.loc[selected_features]
            X_test_selected = X_test.T.loc[selected_features]

            estimator = ElasticNet(max_iter=10000, random_state=42, **param_combinations[min_mse_idx])
            try:
                estimator.fit(X_train_selected.T, y_train.values)
                features_coef = dict(zip(selected_features, estimator.coef_.tolist()))
                y_pred_arr = estimator.predict(X_test_selected.T)
                y_pred = float(y_pred_arr[0])
            except Exception:
                y_pred = np.nan

        results.append({
            'protein': protein_name,
            'LOO sample': sample,
            'true': float(y_test.values[0]),
            'predicted': y_pred,
            'params_dict': best_params,
            'features_coef_dict': features_coef,
            'percentile': percentile
        })

    return pd.DataFrame(results)


def cv_pipeline_fixed_p(model,
                        param_grid: dict,
                        X: pd.DataFrame,
                        y: pd.DataFrame,
                        protein_name: str,
                        feature_selector: str = 'pearson',
                        n_splits: int = 5,
                        n_samples=None,
                        p_value: float = 0.05) -> pd.DataFrame:
    """
    LOO CV for a single protein with a fixed p-value cutoff for feature selection.
    Model: ElasticNet only. `model` is ignored unless it's an ElasticNet instance.
    """
    if protein_name not in y.columns:
        raise ValueError(f"Protein {protein_name} not found in y columns.")

    if (model is not None) and not isinstance(model, ElasticNet):
        raise ValueError("Only ElasticNet is supported. Pass an ElasticNet instance or None.")

    y_protein = y[protein_name]
    sample_list = X.index if n_samples is None else X.index[:n_samples]
    X = X.loc[sample_list]
    y_protein = y_protein.loc[sample_list]

    results = []

    keys, values = zip(*param_grid.items()) if len(param_grid) else ([], [])
    param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)] if len(keys) else [{}]

    for sample in sample_list:
        X_train = X.drop(index=[sample])
        X_test = X.loc[[sample]]
        y_train = y_protein.drop(index=[sample])
        y_test = y_protein.loc[[sample]]

        columns_train = X_train.index
        avg_mse = np.zeros(len(param_combinations), dtype=float)

        kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

        for fit_idx, val_idx in kf.split(columns_train):
            fit_samples = columns_train[fit_idx]
            val_samples = columns_train[val_idx]

            X_fit = X_train.loc[fit_samples].T  # features x samples
            y_fit = y_train.loc[fit_samples]
            X_val = X_train.loc[val_samples].T
            y_val = y_train.loc[val_samples]

            selected_features = _select_features_pvalue(X_fit, y_fit, feature_selector, p_value)
            if len(selected_features) == 0:
                continue

            X_fit_selected = X_fit.loc[selected_features]
            X_val_selected = X_val.loc[selected_features]

            for k_idx, params in enumerate(param_combinations):
                estimator = ElasticNet(max_iter=10000, random_state=42, **params)
                try:
                    estimator.fit(X_fit_selected.T, y_fit.values)
                    y_pred_val = estimator.predict(X_val_selected.T)
                    mse = mean_squared_error(y_val.values, y_pred_val)
                    avg_mse[k_idx] += mse / n_splits
                except Exception:
                    avg_mse[k_idx] += 1e6 / n_splits

        min_mse_idx = int(np.argmin(avg_mse)) if np.isfinite(avg_mse).any() else 0
        best_params = {**param_combinations[min_mse_idx],
                       'p_value': p_value,
                       'feature_selector': feature_selector}

        # Re-select features on the full training set
        X_train_fx = X_train.T
        selected_features = _select_features_pvalue(X_train_fx, y_train, feature_selector, p_value)

        features_coef = None
        y_pred = np.nan

        if len(selected_features) > 0:
            X_train_selected = X_train_fx.loc[selected_features]
            X_test_selected = X_test.T.loc[selected_features]

            estimator = ElasticNet(max_iter=10000, random_state=42, **param_combinations[min_mse_idx])
            try:
                estimator.fit(X_train_selected.T, y_train.values)
                features_coef = dict(zip(selected_features, estimator.coef_.tolist()))
                y_pred_arr = estimator.predict(X_test_selected.T)
                y_pred = float(y_pred_arr[0])
            except Exception:
                y_pred = np.nan

        results.append({
            'protein': protein_name,
            'LOO sample': sample,
            'true': float(y_test.values[0]),
            'predicted': y_pred,
            'params_dict': best_params,
            'features_coef_dict': features_coef,
            'p_value': p_value
        })

    return pd.DataFrame(results)


def cv_pipeline_fixed_features(model,
                               param_grid: dict,
                               X: pd.DataFrame,
                               y: pd.DataFrame,
                               protein_name: str,
                               features: pd.DataFrame,
                               n_splits: int = 5,
                               n_samples=None) -> pd.DataFrame:
    """
    LOO CV for a single protein with a fixed set of features (provided per sample).
    `features` must be a DataFrame indexed by sample with a column 'features' (list of feature names).
    Model: ElasticNet only. `model` is ignored unless it's an ElasticNet instance.
    """
    if protein_name not in y.columns:
        raise ValueError(f"Protein {protein_name} not found in y columns.")

    if (model is not None) and not isinstance(model, ElasticNet):
        raise ValueError("Only ElasticNet is supported. Pass an ElasticNet instance or None.")

    if 'features' not in features.columns:
        raise ValueError("`features` DataFrame must contain a 'features' column with lists of feature names per sample.")

    y_protein = y[protein_name]
    sample_list = X.index if n_samples is None else X.index[:n_samples]
    X = X.loc[sample_list]
    y_protein = y_protein.loc[sample_list]

    results = []

    keys, values = zip(*param_grid.items()) if len(param_grid) else ([], [])
    param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)] if len(keys) else [{}]

    for sample in sample_list:
        X_train = X.drop(index=[sample])
        X_test = X.loc[[sample]]
        y_train = y_protein.drop(index=[sample])
        y_test = y_protein.loc[[sample]]

        columns_train = X_train.index
        avg_mse = np.zeros(len(param_combinations), dtype=float)

        # features for this LOO sample
        selected_features = features.loc[sample, 'features']
        if not isinstance(selected_features, (list, tuple, np.ndarray)):
            raise ValueError("Each entry in `features['features']` must be a list/tuple/array of feature names.")

        # Inner CV for hyperparameters
        kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

        for fit_idx, val_idx in kf.split(columns_train):
            fit_samples = columns_train[fit_idx]
            val_samples = columns_train[val_idx]

            X_fit = X_train.loc[fit_samples].T  # features x samples
            y_fit = y_train.loc[fit_samples]
            X_val = X_train.loc[val_samples].T
            y_val = y_train.loc[val_samples]

            if len(selected_features) == 0:
                continue

            # Align to available features
            sel = [f for f in selected_features if f in X_fit.index]
            if len(sel) == 0:
                continue

            X_fit_selected = X_fit.loc[sel]
            X_val_selected = X_val.loc[sel]

            for k_idx, params in enumerate(param_combinations):
                estimator = ElasticNet(max_iter=10000, random_state=42, **params)
                try:
                    estimator.fit(X_fit_selected.T, y_fit.values)
                    y_pred_val = estimator.predict(X_val_selected.T)
                    mse = mean_squared_error(y_val.values, y_pred_val)
                    avg_mse[k_idx] += mse / n_splits
                except Exception:
                    avg_mse[k_idx] += 1e6 / n_splits

        min_mse_idx = int(np.argmin(avg_mse)) if np.isfinite(avg_mse).any() else 0
        best_params = {**param_combinations[min_mse_idx]}

        # Final train/test
        y_pred = np.nan
        features_coef = None

        if len(selected_features) > 0:
            X_train_fx = X_train.T
            sel = [f for f in selected_features if f in X_train_fx.index]
            if len(sel) > 0:
                X_train_selected = X_train_fx.loc[sel]
                X_test_selected = X_test.T.loc[sel]
                estimator = ElasticNet(max_iter=10000, random_state=42, **param_combinations[min_mse_idx])
                try:
                    estimator.fit(X_train_selected.T, y_train.values)
                    features_coef = dict(zip(sel, estimator.coef_.tolist()))
                    y_pred_arr = estimator.predict(X_test_selected.T)
                    y_pred = float(y_pred_arr[0])
                except Exception:
                    y_pred = np.nan

        results.append({
            'protein': protein_name,
            'LOO sample': sample,
            'true': float(y_test.values[0]),
            'predicted': y_pred,
            'params_dict': best_params,
            'features_coef_dict': features_coef
        })

    return pd.DataFrame(results)
