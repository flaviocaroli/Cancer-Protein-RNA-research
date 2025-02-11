import pandas as pd
import numpy as np
import itertools
from sklearn.model_selection import KFold
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr, pearsonr

def cv_pipeline_variable_p(model, param_grid, X, y, protein_name, feature_selector='anova', n_splits=5, n_samples=None, p_values=None):
    """
    Perform leave-one-out (LOO) cross-validation for a single protein with variable p-values for feature selection.
    For each sample:
      - Hold out the sample as the test set.
      - Perform K-fold CV on the remaining samples to select the best model parameters and features.
      - Retrain on the remaining samples and predict the held-out sample.
    
    Parameters:
    - model: sklearn-like estimator (e.g., ElasticNet)
    - param_grid: dict of parameters to try
    - X: pandas DataFrame of transcriptomics data (samples x features)
    - y: pandas DataFrame of proteomics data (samples x proteins)
    - protein_name: str, the name of the protein to analyze
    - feature_selector: 'anova', 'spearman', or 'pearson' for feature selection
    - n_splits: int, number of folds for K-fold CV
    - n_samples: int, optional, to limit the number of samples (for debugging)
    - p_values: list of floats, p-values for feature selection (default: [0.1, 0.05, 0.01])

    Returns:
    - df: DataFrame with results for each LOO sample, including true and predicted values, parameters, and selected features.
    """
    if protein_name not in y.columns:
        raise ValueError(f"Protein {protein_name} not found in y columns.")
    
    if p_values is None:
        p_values = [0.1, 0.05, 0.01]
    
    y_protein = y[protein_name]
    sample_list = X.index if n_samples is None else X.index[:n_samples]
    X = X.loc[sample_list]
    y_protein = y_protein.loc[sample_list]
    
    df = pd.DataFrame(columns=['protein', 'LOO sample', 'true', 'predicted', 
                               'params_dict', 'features_coef_dict', 'p_value'])
    
    keys, values = zip(*param_grid.items())
    param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)]
    
    for sample in sample_list:
        X_train = X.drop(index=[sample])
        X_test = X.loc[[sample]]
        y_train = y_protein.drop(index=[sample])
        y_test = y_protein.loc[[sample]]
        
        for p_value in p_values:
            kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
            avg_mse = np.zeros(len(param_combinations))
            best_features = None
            
            for fit_idx, val_idx in kf.split(X_train):
                X_fit, X_val = X_train.iloc[fit_idx], X_train.iloc[val_idx]
                y_fit, y_val = y_train.iloc[fit_idx], y_train.iloc[val_idx]
                
                # Feature selection based on the chosen method
                if feature_selector == 'anova':
                    selector = SelectKBest(score_func=f_regression, k='all')
                    selector.fit(X_fit, y_fit)
                    pvalues = selector.pvalues_
                elif feature_selector == 'spearman':
                    pvalues = np.array([spearmanr(X_fit.iloc[:, i], y_fit)[1] for i in range(X_fit.shape[1])])
                elif feature_selector == 'pearson':
                    pvalues = np.array([pearsonr(X_fit.iloc[:, i], y_fit)[1] for i in range(X_fit.shape[1])])
                
                mask = pvalues <= p_value
                if not np.any(mask):
                    continue
                
                selected_features = X_train.columns[mask]
                
                # Evaluate each parameter combination
                for k_idx, params in enumerate(param_combinations):
                    estimator = model.__class__(max_iter=100)
                    estimator.set_params(**params)
                    try:
                        estimator.fit(X_fit[selected_features], y_fit)
                        y_pred_val = estimator.predict(X_val[selected_features])
                        mse = mean_squared_error(y_val, y_pred_val)
                        avg_mse[k_idx] += mse / n_splits
                    except:
                        avg_mse[k_idx] += 1e6 / n_splits
            
            # Select the best parameters based on the lowest MSE
            min_mse_idx = np.argmin(avg_mse)
            best_params = {**param_combinations[min_mse_idx], 'p_value': p_value, 'feature_selector': feature_selector}
            
            if best_features is None:
                best_features = selected_features
            
            X_train_selected = X_train[best_features]
            estimator = model.__class__(max_iter=100)
            estimator.set_params(**param_combinations[min_mse_idx])
            
            features_coef = None
            y_pred = np.nan
            
            try:
                estimator.fit(X_train_selected, y_train)
                features_coef = dict(zip(best_features, estimator.coef_))
                X_test_selected = X_test[best_features]
                y_pred = estimator.predict(X_test_selected)
            except:
                y_pred = np.nan
            
            # Store results in the DataFrame
            df.loc[len(df)] = {
                'protein': protein_name,
                'LOO sample': sample,
                'true': y_test.values[0],
                'predicted': y_pred[0] if isinstance(y_pred, np.ndarray) else y_pred,
                'params_dict': best_params,
                'features_coef_dict': features_coef,
                'p_value': p_value
            }
    
    return df


def cv_pipeline_fixed_p(model, param_grid, X, y, protein_name, feature_selector='anova', n_splits=5, n_samples=None, p_value=0.05):
    """
    Perform leave-one-out (LOO) cross-validation for a single protein with a fixed p-value for feature selection.
    For each sample:
      - Hold out the sample as the test set.
      - Perform K-fold CV on the remaining samples to select the best model parameters.
      - Retrain on the remaining samples and predict the held-out sample.
    
    Parameters:
    - model: sklearn-like estimator (e.g., ElasticNet)
    - param_grid: dict of parameters to try
    - X: pandas DataFrame of transcriptomics data (samples x features)
    - y: pandas DataFrame of proteomics data (samples x proteins)
    - protein_name: str, the name of the protein to analyze
    - feature_selector: 'anova', 'spearman', or 'pearson' for feature selection
    - n_splits: int, number of folds for K-fold CV
    - n_samples: int, optional, to limit the number of samples (for debugging)
    - p_value: float, fixed p-value for feature selection (default: 0.05)

    Returns:
    - df: DataFrame with results for each LOO sample, including true and predicted values, parameters, and selected features.
    """
    if protein_name not in y.columns:
        raise ValueError(f"Protein {protein_name} not found in y columns.")
    
    y_protein = y[protein_name]
    
    sample_list = X.index if n_samples is None else X.index[:n_samples]
    X = X.loc[sample_list]
    y_protein = y_protein.loc[sample_list]

    df = pd.DataFrame(columns=['protein', 'LOO sample', 'true', 'predicted', 
                             'params_dict', 'features_coef_dict', 'p_value'])

    keys, values = zip(*param_grid.items())
    param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)]

    for sample in sample_list:
        X_train = X.drop(index=[sample])
        X_test = X.loc[[sample]]
        y_train = y_protein.drop(index=[sample])
        y_test = y_protein.loc[[sample]]
        
        columns_train = X_train.index
        avg_mse = np.zeros(len(param_combinations))
        
        kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
        
        for fit_idx, val_idx in kf.split(columns_train):
            fit_samples = columns_train[fit_idx]
            val_samples = columns_train[val_idx]
            
            X_fit = X_train.loc[fit_samples].T
            y_fit = y_train.loc[fit_samples]
            X_val = X_train.loc[val_samples].T
            y_val = y_train.loc[val_samples]
            
            # Feature selection based on the chosen method
            if feature_selector == 'anova':
                selector = SelectKBest(score_func=f_regression, k='all')
                selector.fit(X_fit.T, y_fit)
                pvalues = selector.pvalues_
            elif feature_selector == 'spearman':
                pvalues = np.array([spearmanr(X_fit.loc[i], y_fit)[1] for i in X_fit.index])
            elif feature_selector == 'pearson':
                pvalues = np.array([pearsonr(X_fit.loc[i], y_fit)[1] for i in X_fit.index])
            
            mask = pvalues <= p_value
            
            if not np.any(mask):
                continue
                
            X_fit_selected = X_fit.loc[mask]
            X_val_selected = X_val.loc[mask]
            
            # Evaluate each parameter combination
            for k_idx, params in enumerate(param_combinations):
                estimator = model.__class__(max_iter=100)
                estimator.set_params(**params)
                
                try:
                    estimator.fit(X_fit_selected.T, y_fit)
                    y_pred_val = estimator.predict(X_val_selected.T)
                    mse = mean_squared_error(y_val, y_pred_val)
                    avg_mse[k_idx] += mse / n_splits
                except:
                    avg_mse[k_idx] += 1e6 / n_splits
        
        # Select the best parameters based on the lowest MSE
        min_mse_idx = np.argmin(avg_mse)
        best_params = {**param_combinations[min_mse_idx], 
                     'p_value': p_value,
                     'feature_selector': feature_selector}
        
        selected_features = X_train.T.index[pvalues <= p_value]
        
        features_coef = None
        y_pred = np.nan
        
        if np.any(selected_features):
            X_train_selected = X_train.T.loc[selected_features]
            estimator = model.__class__(max_iter=100)
            estimator.set_params(**param_combinations[min_mse_idx])
            try:
                estimator.fit(X_train_selected.T, y_train)
                features_coef = dict(zip(selected_features, estimator.coef_))
                
                X_test_selected = X_test.T.loc[selected_features]
                y_pred = estimator.predict(X_test_selected.T)
            except:
                y_pred = np.nan

        # Store results in the DataFrame
        df.loc[len(df)] = {
            'protein': protein_name,
            'LOO sample': sample,
            'true': y_test.values[0],
            'predicted': y_pred[0] if isinstance(y_pred, np.ndarray) else y_pred,
            'params_dict': best_params,
            'features_coef_dict': features_coef,
            'p_value': p_value
        }

    return df