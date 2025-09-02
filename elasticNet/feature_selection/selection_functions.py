import numpy as np
import pandas as pd
from sklearn.feature_selection import f_regression
from scipy.stats import pearsonr, spearmanr

clusters = pd.read_csv("protein_clusters_k3.csv").set_index("Protein")['Cluster']

max_features = {
    1:8600,
    2:10400,
    3:6700
}

def feature_selection(X, y, protein_name, percentile=5, n_samples=None):
    """
    Perform leave-one-out (LOO) feature selection for a single protein using a fixed percentile threshold.
    
    For each left-out sample, the training set is used to compute three types of feature selection:
      - 'anova': Uses f_regression to obtain F scores and p-values.
      - 'pearson': Computes Pearson correlation coefficients and p-values for each feature.
      - 'spearman': Computes Spearman correlation coefficients and p-values for each feature.
      
    For each method, the threshold is defined as the given percentile of the p-values computed over the training set.
    Only features with p-values below or equal to that threshold are selected, and their weight is defined as:
      - anova: the F score (higher means more predictive)
      - pearson: the Pearson correlation coefficient
      - spearman: the Spearman correlation coefficient

    Parameters:
        X (pd.DataFrame): DataFrame of predictors (samples x features).
        y (pd.DataFrame): DataFrame of responses (samples x proteins).
        protein_name (str): The protein column name in y to analyze.
        percentile (float): The percentile threshold (0–100) to select features based on p-values.
        n_samples (int, optional): If provided, only the first n_samples in X and y are used.

    Returns:
        pd.DataFrame: A DataFrame whose index corresponds to the left-out sample and three columns:
                      'anova', 'pearson', 'spearman'. Each cell contains a dict mapping selected
                      feature names to their weight.
    """
    if protein_name not in y.columns:
        raise ValueError(f"Protein {protein_name} not found in y columns.")
    
    y_protein = y[protein_name]
    sample_list = X.index if n_samples is None else X.index[:n_samples]
    X = X.loc[sample_list]
    y_protein = y_protein.loc[sample_list]
    
    results = []

    # Loop over each sample to leave out one at a time
    for sample in sample_list:
        # Define training set (all samples except the one left out)
        X_train = X.drop(index=[sample])
        y_train = y_protein.drop(index=[sample])
        
        features = X_train.columns
        
        # --- ANOVA Feature Selection ---
        # f_regression returns F scores and p-values.
        F_scores, pvals_anova = f_regression(X_train, y_train)
        threshold_anova = np.percentile(pvals_anova, percentile)
        selected_anova = {feat: F_scores[i] 
                          for i, feat in enumerate(features) 
                          if pvals_anova[i] <= threshold_anova}
        
        # --- Pearson Feature Selection ---
        pearson_results = [pearsonr(X_train[col].values, y_train.values) for col in features]
        r_pearson = [res[0] for res in pearson_results]
        pvals_pearson = [res[1] for res in pearson_results]
        threshold_pearson = np.percentile(pvals_pearson, percentile)
        selected_pearson = {feat: r_pearson[i]
                            for i, feat in enumerate(features)
                            if pvals_pearson[i] <= threshold_pearson}
        
        # --- Spearman Feature Selection ---
        spearman_results = [spearmanr(X_train[col].values, y_train.values) for col in features]
        r_spearman = [res[0] for res in spearman_results]
        pvals_spearman = [res[1] for res in spearman_results]
        threshold_spearman = np.percentile(pvals_spearman, percentile)
        selected_spearman = {feat: r_spearman[i]
                             for i, feat in enumerate(features)
                             if pvals_spearman[i] <= threshold_spearman}
        
        results.append({
            'LOO_sample': sample,
            'anova': selected_anova,
            'pearson': selected_pearson,
            'spearman': selected_spearman
        })
    
    # Create DataFrame with left-out sample as index and columns for each selection method
    results_df = pd.DataFrame(results).set_index('LOO_sample')
    return results_df


def naive_cluster_feature_selection(X, y, protein_name):
    """
    Perform leave-one-out (LOO) feature selection for a single protein using a fixed percentile threshold.
    
    For each left-out sample, the training set is used to compute Pearson's correlation
    coefficient between the protein profile and each gene expression profile.

    The number of selected features varies from protein to protein based on their cluster.
    The features with the highest correletion coefficient are kept and saved.

    Parameters:
        X (pd.DataFrame): DataFrame of predictors (samples x features).
        y (pd.DataFrame): DataFrame of responses (samples x proteins).
        protein_name (str): The protein column name in y to analyze.

    Returns:
        pd.DataFrame: A DataFrame whose index corresponds to the left-out sample and one column:
                      'pearson'. Each cell contains a list of the selected genes, per sample.
    """
    if protein_name not in y.columns:
        raise ValueError(f"Protein {protein_name} not found in y columns.")
    
    y_protein = y[protein_name]
    sample_list = X.index
    X = X.loc[sample_list]
    y_protein = y_protein.loc[sample_list]
    
    results = []

    # Loop over each sample to leave out one at a time
    for sample in sample_list:

        # Define training set (all samples except the one left out)
        X_train = X.drop(index=[sample])
        y_train = y_protein.drop(index=[sample])
        
        features = X_train.columns
        
        # --- Pearson Feature Selection ---
        results_pearson = {gene : pearsonr(X_train[gene].values, y_train.values)[0] for gene in features}
        threshold_pearson = max_features[clusters[protein_name]]
        selected_pearson = sorted(results_pearson, key=results_pearson.get, reverse=True)[:threshold_pearson]
        
        results.append({
            'sample': sample,
            'pearson': selected_pearson
        })
    
    # Create DataFrame with left-out sample as index and columns for each selection method
    results_df = pd.DataFrame(results).set_index('sample')
    return results_df