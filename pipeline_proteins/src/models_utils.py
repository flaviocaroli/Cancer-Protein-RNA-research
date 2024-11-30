import pandas as pd
import numpy as np
import itertools
from sklearn.model_selection import KFold
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr, pearsonr


def full_pipeline2(model,
                            param_grid,
                            X,
                            y,
                            feature_selector='anova',
                            n_splits=5,
                            n_proteins=None,
                            n_samples=None):
    """
    Perform adaptive pipeline for regression model tuning with different feature selection methods.

    Parameters:
    - model: sklearn-like estimator
    - param_grid: dict (should include 'p_value' as one of the hyperparameters)
    - X: pandas DataFrame
    - y: pandas DataFrame
    - feature_selector: str, optional (default='anova')
        The feature selection method to use. Options are 'anova', 'spearman', or 'pearson'.
    - n_splits: int, optional (default=5)
        Number of splits for K-fold cross-validation.
    - n_proteins: int, optional (default=None)
        Number of proteins to use for training. If None, all proteins are used.
    - n_samples: int, optional (default=None)
        Number of samples to use for training. If None, all samples are used.

    Returns:
    - df: pandas DataFrame
    """
    # Transpose X and y to use labels instead of indices for better performance
    X = X.T
    y = y.T

    sample_list = y.columns
    protein_list = y.index

    if n_samples is not None:
        sample_list = sample_list[:n_samples]
    if n_proteins is not None:
        protein_list = protein_list[:n_proteins]

    df = pd.DataFrame(columns=['protein', 'LOO sample', 'true', 'predicted', 'params'])

    # Generate all combinations of hyperparameters, including p_value
    keys, values = zip(*param_grid.items())
    param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)]

    for protein in protein_list:
        yp = y.loc[protein]

        for sample in sample_list:
            # Leave the sample out
            X_train = X.drop(columns=[sample])
            X_test = X[[sample]]
            y_train = yp.drop(index=[sample])
            y_test = yp[[sample]]

            kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
            avg_mse = np.zeros(len(param_combinations))

            for fold_idx, (fit_idx, val_idx) in enumerate(kf.split(X_train.columns)):
                fit_columns = X_train.columns[fit_idx]
                val_columns = X_train.columns[val_idx]

                X_fit = X_train[fit_columns]
                y_fit = y_train[fit_columns]
                X_val = X_train[val_columns]
                y_val = y_train[val_columns]

                # Compute p-values once per fold
                if feature_selector == 'anova':
                    selector = SelectKBest(score_func=f_regression, k='all')
                    selector.fit(X_fit.T, y_fit)
                    pvalues = selector.pvalues_
                elif feature_selector == 'spearman':
                    pvalues = np.array([spearmanr(X_fit.loc[i], y_fit)[1] for i in X_fit.index])
                elif feature_selector == 'pearson':
                    pvalues = np.array([pearsonr(X_fit.loc[i], y_fit)[1] for i in X_fit.index])
                else:
                    raise ValueError("feature_selector must be 'anova', 'spearman', or 'pearson'")

                for k, params in enumerate(param_combinations):
                    # Make a copy of params to avoid modifying the original
                    params_copy = params.copy()
                    # Extract p_value and remove it from params_copy
                    p_value = params_copy.pop('p_value')

                    mask = pvalues <= p_value

                    # If no features are selected, penalize and continue
                    if not np.any(mask):
                        avg_mse[k] += 1e6 / n_splits  # Penalize for no features selected
                        continue

                    X_fit_selected = X_fit.loc[mask]
                    X_val_selected = X_val.loc[mask]

                    estimator = model.__class__(max_iter=10)  # Increase max_iter if needed
                    estimator.set_params(**params_copy)

                    try:
                        estimator.fit(X_fit_selected.T, y_fit)
                        y_pred = estimator.predict(X_val_selected.T)
                        mse = mean_squared_error(y_val, y_pred)
                        avg_mse[k] += mse / n_splits
                    except Exception as e:
                        # In case of any error during fitting/predicting, set mse to a large number
                        avg_mse[k] += 1e6 / n_splits
                        # Optionally, print or log the exception e

            # Select the best parameters based on avg_mse
            best_index = np.argmin(avg_mse)
            # Make a copy of best_params to avoid modifying the original
            best_params = param_combinations[best_index].copy()
            best_p_value = best_params.pop('p_value')

            # Perform feature selection on full X_train using best_p_value
            if feature_selector == 'anova':
                selector = SelectKBest(score_func=f_regression, k='all')
                selector.fit(X_train.T, y_train)
                pvalues = selector.pvalues_
                mask = pvalues <= best_p_value
            elif feature_selector == 'spearman':
                pvalues = np.array([spearmanr(X_train.loc[i], y_train)[1] for i in X_train.index])
                mask = pvalues <= best_p_value
            elif feature_selector == 'pearson':
                pvalues = np.array([pearsonr(X_train.loc[i], y_train)[1] for i in X_train.index])
                mask = pvalues <= best_p_value

            # If no features are selected, skip to next iteration
            if not np.any(mask):
                continue

            X_train_selected = X_train.loc[mask]
            X_test_selected = X_test.loc[mask]

            
            estimator = model.__class__(max_iter=100)  
            estimator.set_params(**best_params)

            try:
                estimator.fit(X_train_selected.T, y_train)
                y_pred = estimator.predict(X_test_selected.T)
            except Exception as e:
                y_pred = np.nan  # If fitting fails, set prediction to NaN
                # Optionally, print or log the exception e

            # Add p_value back to params for recording
            best_params['p_value'] = best_p_value
            best_params['feature_selector'] = feature_selector

            df.loc[len(df)] = {
                'protein': protein,
                'LOO sample': sample,
                'true': y_test.values[0],
                'predicted': y_pred[0] if isinstance(y_pred, np.ndarray) else y_pred,
                'params': best_params
            }

    return df
