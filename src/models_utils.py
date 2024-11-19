
from sklearn.model_selection import KFold, train_test_split
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.metrics import mean_squared_error
import pandas as pd
import numpy as np
import itertools

## adaptive loop for the model

def full_adaptive_pipeline(model,
                           param_grid,
                           X,
                           y,
                           anova_features,
                           n_splits=5,
                           n_proteins=None,
                           n_samples=None):
    
    # Get samples and proteins
    sample_list = y.index
    protein_list = y.columns

    if n_samples != None:
        sample_list = sample_list[:n_samples]
    if n_proteins != None:
        protein_list = protein_list[:n_proteins]
        

    # Initialize the results table
    df = pd.DataFrame(columns=['protein', 'LOO sample', 'true', 'predicted', 'params'])

    # Generate all parameter combinations
    keys, values = zip(*param_grid.items())
    param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)]
    del keys, values

    # Nested loops for proteins and samples
    for protein in protein_list:
        # Target variable for the current protein
        yp = y[protein]

        for sample in sample_list:
            # Leave-One-Out (LOO) split
            X_train, X_test, y_train, y_test = train_test_split(X, yp, test_size=1)

            # K-Fold cross-validation
            kf = KFold(n_splits=n_splits, shuffle=True)
            folds = list(kf.split(X_train.index))

            # Initialize average MSE list
            avg_mse = np.zeros(len(param_combinations))

            for fit_idx, val_idx in folds:
                # Construct training and validation sets
                X_fit = X_train.iloc[fit_idx]
                y_fit = y_train.iloc[fit_idx]
                X_val = X_train.iloc[val_idx]
                y_val = y_train.iloc[val_idx]

                # Feature selection using ANOVA F-test
                selector = SelectKBest(score_func=f_regression, k=min(anova_features, X_fit.shape[1]))
                X_fit = selector.fit_transform(X_fit, y_fit)
                X_val = selector.transform(X_val)

                for k, params in enumerate(param_combinations):
                    # Set model parameters
                    estimator = model
                    estimator.set_params(**params)

                    # Fit the model and predict
                    estimator.fit(X_fit, y_fit)
                    y_pred = estimator.predict(X_val)

                    # Calculate MSE and update the average
                    mse = mean_squared_error(y_val, y_pred)
                    avg_mse[k] += mse / n_splits

            # Find the best parameter combination (minimum MSE)
            best_index = np.argmin(avg_mse)
            best_params = param_combinations[best_index]

            # Train the final model using the best parameters
            estimator = model
            estimator.set_params(**best_params)

            # Feature selection for the final model
            X_train_selected = selector.fit_transform(X_train, y_train)
            X_test_selected = selector.transform(X_test)

            estimator.fit(X_train_selected, y_train)
            y_pred = estimator.predict(X_test_selected)

            # Add the result to the DataFrame
            df.loc[len(df)] = [
                protein,
                sample,
                y_test.values[0],
                y_pred[0],
                best_params
            ]

    return df

# OLD ONE - Samples outside  - don't use

def boosting_pipeline(model, param_grid, X, y, anova_features):

    # get samples and proteins
    sample_list = y.index
    protein_list = y.columns

    # initialize the results table
    df = pd.DataFrame(columns=['protein',
                               'LOO sample',
                               'true',
                               'predicted',
                               'params'])

    for sample in sample_list:
        print(sample)

        # leave one out split
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=1)

        # get folds
        kf = KFold(n_splits=5, shuffle=True)
        folds = list(kf.split(X_train.index))

        # get parameter combinations
        keys, values = zip(*param_grid.items())
        param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)]

        for protein in protein_list:
            
            # remove all proteins but one from y
            yp_train = y_train[protein]
            yp_test = y_test[protein]

            # select the best features by F-value
            selector = SelectKBest(score_func=f_regression, k=anova_features)
            selector.fit_transform(X_train, yp_train)
            
            # remove all other features from X
            selected_indices = np.where(selector.get_support())[0]
            selected_features = X_train.columns[selected_indices]
            Xp_train = X_train[selected_features]
            Xp_test = X_test[selected_features]

            # initialize best_params variable
            best_params = param_combinations[0]
            best_params_mse = np.inf

            for params in param_combinations:

                avg_mse = 0

                for fold in folds:

                    # fit and validation indices
                    fit_idx, val_idx = fold

                    # construct sets
                    X_fit = Xp_train.iloc[fit_idx]
                    y_fit = yp_train.iloc[fit_idx]
                    X_val = Xp_train.iloc[val_idx]
                    y_val = yp_train.iloc[val_idx]

                    # build model for this fold and get predictions
                    estimator = model
                    estimator.set_params(**params)

                    estimator.fit(X_fit, y_fit)
                    y_pred = estimator.predict(X_val)

                    # get mse on this fold
                    mse = mean_squared_error(y_val, y_pred)
                    avg_mse += mse
                
                # get mse for parameter configuration
                avg_mse /= len(folds)

                # update best config if these parameters have the lowest mean mse
                if avg_mse < best_params_mse:
                    best_params_mse = avg_mse
                    best_params = params
                
            # train the final model with the best parameters
            estimator = model
            estimator.set_params(**best_params)

            estimator.fit(Xp_train, yp_train)
            y_pred = model.predict(Xp_test)

            # add entry to the results table
            df.loc[len(df)] = [
                protein,
                sample,
                yp_test.values[0],
                y_pred[0],
                best_params
            ]
        
    return df