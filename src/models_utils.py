import pandas as pd
import numpy as np
import itertools
from sklearn.model_selection import KFold
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr, pearsonr

def full_adaptive_pipeline(model,
                           param_grid,
                           X,
                           y,
                           feature_selector='anova',
                           p_values=[0.1, 0.05, 0.01],
                           n_splits=5,
                           n_proteins=None,
                           n_samples=None):
    """
    Perform adaptive pipeline for regression model tuning with different feature selection methods.

    Parameters:
    - model: sklearn-like estimator
    - param_grid: dict
    - X: pandas DataFrame
    - y: pandas DataFrame
    - feature_selector: str, optional (default='anova')
        The feature selection method to use. Options are 'anova', 'spearman', or 'pearson'.
    - p_values: list of float, optional (default=[0.1, 0.05, 0.01])
        List of p-values to use for feature selection during hyperparameter tuning.
    - n_splits: int, optional (default=5)
        Number of splits for K-fold cross-validation.
    - n_proteins: int, optional (default=None)
        Number of proteins to use for training. If None, all proteins are used.
    - n_samples: int, optional (default=None)
        Number of samples to use for training. If None, all samples are used.

    Returns:
    - df: pandas DataFrame
    """
    sample_list = y.index
    protein_list = y.columns

    if n_samples is not None:
        sample_list = sample_list[:n_samples]
    if n_proteins is not None:
        protein_list = protein_list[:n_proteins]

    df = pd.DataFrame(columns=['protein', 'LOO sample', 'true', 'predicted', 'params'])

    keys, values = zip(*param_grid.items())
    param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)]

    for protein in protein_list:
        yp = y[protein]

        for sample in sample_list:
            # Leave the sample out
            X_train = X.drop(sample)
            X_test = X.loc[[sample]]
            y_train = yp.drop(sample)
            y_test = yp.loc[[sample]]

            kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
            avg_mse = np.zeros(len(param_combinations) * len(p_values))

            for fold_idx, (fit_idx, val_idx) in enumerate(kf.split(X_train)):
                X_fit = X_train.iloc[fit_idx]
                y_fit = y_train.iloc[fit_idx]
                X_val = X_train.iloc[val_idx]
                y_val = y_train.iloc[val_idx]

                for k, params in enumerate(param_combinations):
                    for p_idx, p_value in enumerate(p_values):
                        if feature_selector == 'anova':
                            selector = SelectKBest(score_func=f_regression, k='all')
                            selector.fit(X_fit, y_fit)
                            mask = selector.pvalues_ <= p_value
                        elif feature_selector == 'spearman':
                            scores = np.array([spearmanr(X_fit.iloc[:, i], y_fit)[0] for i in range(X_fit.shape[1])])
                            pvalues = np.array([spearmanr(X_fit.iloc[:, i], y_fit)[1] for i in range(X_fit.shape[1])])
                            mask = pvalues <= p_value
                        elif feature_selector == 'pearson':
                            scores = np.array([pearsonr(X_fit.iloc[:, i], y_fit)[0] for i in range(X_fit.shape[1])])
                            pvalues = np.array([pearsonr(X_fit.iloc[:, i], y_fit)[1] for i in range(X_fit.shape[1])])
                            mask = pvalues <= p_value
                        else:
                            raise ValueError("feature_selector must be 'anova', 'spearman', or 'pearson'")

                        # If no features are selected, continue to next iteration
                        if not np.any(mask):
                            continue

                        X_fit_selected = X_fit.loc[:, mask]
                        X_val_selected = X_val.loc[:, mask]

                        estimator = model.__class__()  # Create a new instance of the model
                        estimator.set_params(**params)

                        try:
                            estimator.fit(X_fit_selected, y_fit)
                            y_pred = estimator.predict(X_val_selected)
                            mse = mean_squared_error(y_val, y_pred)
                            avg_mse[p_idx * len(param_combinations) + k] += mse / n_splits
                        except Exception:
                            # In case of any error during fitting/predicting, set mse to a large number
                            avg_mse[p_idx * len(param_combinations) + k] += 1e6 / n_splits

            # Select the best parameters and p_value
            best_index = np.argmin(avg_mse)
            best_p_value_idx = best_index // len(param_combinations)
            best_param_idx = best_index % len(param_combinations)
            best_params = param_combinations[best_param_idx]
            best_p_value = p_values[best_p_value_idx]

            # Perform feature selection on full X_train
            if feature_selector == 'anova':
                selector = SelectKBest(score_func=f_regression, k='all')
                selector.fit(X_train, y_train)
                mask = selector.pvalues_ <= best_p_value
            elif feature_selector == 'spearman':
                scores = np.array([spearmanr(X_train.iloc[:, i], y_train)[0] for i in range(X_train.shape[1])])
                pvalues = np.array([spearmanr(X_train.iloc[:, i], y_train)[1] for i in range(X_train.shape[1])])
                mask = pvalues <= best_p_value
            elif feature_selector == 'pearson':
                scores = np.array([pearsonr(X_train.iloc[:, i], y_train)[0] for i in range(X_train.shape[1])])
                pvalues = np.array([pearsonr(X_train.iloc[:, i], y_train)[1] for i in range(X_train.shape[1])])
                mask = pvalues <= best_p_value

            # If no features are selected, skip to next iteration
            if not np.any(mask):
                continue

            X_train_selected = X_train.loc[:, mask]
            X_test_selected = X_test.loc[:, mask]

            estimator = model.__class__()  # Create a new instance of the model
            estimator.set_params(**best_params)

            try:
                estimator.fit(X_train_selected, y_train)
                y_pred = estimator.predict(X_test_selected)
            except Exception:
                y_pred = np.nan  # If fitting fails, set prediction to NaN

            df.loc[len(df)] = {
                'protein': protein,
                'LOO sample': sample,
                'true': y_test.values[0],
                'predicted': y_pred[0] if isinstance(y_pred, np.ndarray) else y_pred,
                'params': {**best_params, 'p_value': best_p_value, 'feature_selector': feature_selector}
            }

    return df
