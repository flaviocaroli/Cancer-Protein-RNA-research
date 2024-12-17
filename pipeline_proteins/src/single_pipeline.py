import pandas as pd
import numpy as np
import itertools
from sklearn.model_selection import KFold
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr, pearsonr

def single_protein_pipeline(model,
                            param_grid,
                            X,
                            y,
                            protein_name,
                            feature_selector='anova',
                            n_splits=5,
                            n_samples=None):
    """
    Perform a leave-one-out (LOO) scenario for a single protein.
    For each sample in the dataset:
      - That sample is held out as test.
      - Parameter selection is done via K-fold cross-validation on the remaining samples.
      - After finding the best parameters, retrain on all remaining samples and predict for the held-out sample.
    
    Parameters:
    - model: sklearn-like estimator (e.g., ElasticNet)
    - param_grid: dict of parameters (including 'p_value') to try
    - X: pandas DataFrame of transcriptomics data (samples x features)
    - y: pandas DataFrame of proteomics data (samples x proteins)
    - protein_name: str, the name of the protein to analyze
    - feature_selector: 'anova', 'spearman', or 'pearson'
    - n_splits: int, for K-fold in parameter selection
    - n_samples: int, optional, to limit number of samples (for debugging)

    Returns:
    - df: DataFrame with columns ['protein', 'LOO sample', 'true', 'predicted', 'params']
    """

    if protein_name not in y.columns:
        raise ValueError(f"Protein {protein_name} not found in y columns.")

    y_protein = y[protein_name]

    sample_list = X.index
    if n_samples is not None:
        sample_list = sample_list[:n_samples]
        X = X.loc[sample_list]
        y_protein = y_protein.loc[sample_list]

    # Prepare results dataframe
    df = pd.DataFrame(columns=['protein', 'LOO sample', 'true', 'predicted', 'params'])

    # Generate all parameter combinations
    keys, values = zip(*param_grid.items())
    param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)]

    # Loop over each sample as LOO
    for sample in sample_list:
        # Split into training and test
        X_train = X.drop(index=[sample])
        X_test = X.loc[[sample]]
        y_train = y_protein.drop(index=[sample])
        y_test = y_protein.loc[[sample]]

        # Perform parameter selection using K-fold on X_train, y_train
        columns_train = X_train.index
        avg_mse = np.zeros(len(param_combinations))

        kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

        for fit_idx, val_idx in kf.split(columns_train):
            fit_samples = columns_train[fit_idx]
            val_samples = columns_train[val_idx]

            X_fit = X_train.loc[fit_samples].T  # features x fit_samples
            y_fit = y_train.loc[fit_samples]
            X_val = X_train.loc[val_samples].T   # features x val_samples
            y_val = y_train.loc[val_samples]

            # Compute p-values once per fold
            if feature_selector == 'anova':
                selector = SelectKBest(score_func=f_regression, k='all')
                selector.fit(X_fit.T, y_fit)
                pvalues = selector.pvalues_
                feature_index = X_fit.index
            elif feature_selector == 'spearman':
                feature_index = X_fit.index
                pvalues = np.array([spearmanr(X_fit.loc[i], y_fit)[1] for i in feature_index])
            elif feature_selector == 'pearson':
                feature_index = X_fit.index
                pvalues = np.array([pearsonr(X_fit.loc[i], y_fit)[1] for i in feature_index])
            else:
                raise ValueError("feature_selector must be 'anova', 'spearman', or 'pearson'")

            for k_idx, params in enumerate(param_combinations):
                params_copy = params.copy()
                p_value = params_copy.pop('p_value')

                mask = pvalues <= p_value

                if not np.any(mask):
                    # No features selected
                    avg_mse[k_idx] += 1e6 / n_splits
                    continue

                X_fit_selected = X_fit.loc[mask]
                X_val_selected = X_val.loc[mask]

                estimator = model.__class__(max_iter=100)
                estimator.set_params(**params_copy)

                try:
                    estimator.fit(X_fit_selected.T, y_fit)
                    y_pred_val = estimator.predict(X_val_selected.T)
                    mse = mean_squared_error(y_val, y_pred_val)
                    avg_mse[k_idx] += mse / n_splits
                except:
                    avg_mse[k_idx] += 1e6 / n_splits

        # Select the best parameters based on avg_mse
        best_index = np.argmin(avg_mse)
        best_params = param_combinations[best_index].copy()
        best_p_value = best_params.pop('p_value')

        # Final training on full X_train with selected p-value
        X_train_t = X_train.T
        if feature_selector == 'anova':
            selector = SelectKBest(score_func=f_regression, k='all')
            selector.fit(X_train_t.T, y_train)
            pvalues = selector.pvalues_
            mask = pvalues <= best_p_value
        elif feature_selector == 'spearman':
            pvalues = np.array([spearmanr(X_train_t.loc[i], y_train)[1] for i in X_train_t.index])
            mask = pvalues <= best_p_value
        elif feature_selector == 'pearson':
            pvalues = np.array([pearsonr(X_train_t.loc[i], y_train)[1] for i in X_train_t.index])
            mask = pvalues <= best_p_value

        if not np.any(mask):
            # No features selected for final training
            y_pred = np.nan
        else:
            X_train_selected = X_train_t.loc[mask]
            X_test_t = X_test.T.loc[mask]

            estimator = model.__class__(max_iter=100)
            estimator.set_params(**best_params)

            try:
                estimator.fit(X_train_selected.T, y_train)
                y_pred = estimator.predict(X_test_t.T)
            except:
                y_pred = np.nan

        # Add p_value and feature_selector back to params for recording
        best_params['p_value'] = best_p_value
        best_params['feature_selector'] = feature_selector

        df.loc[len(df)] = {
            'protein': protein_name,
            'LOO sample': sample,
            'true': y_test.values[0],
            'predicted': y_pred[0] if isinstance(y_pred, np.ndarray) else y_pred,
            'params': best_params
        }

    return df
