import numpy as np
import pandas as pd
import os  # Added os to handle paths dynamically
from sklearn.feature_selection import f_regression
from scipy.stats import pearsonr, spearmanr
from typing import Dict, Optional, Sequence

# Per-cluster caps for number of selected features
MAX_FEATURES: Dict[int, int] = {
    1: 8600,
    2: 10400,
    3: 6700,
}

# --- FIX: Load clusters using an absolute path relative to this script ---
# This prevents FileNotFoundError when importing this module from other folders.
current_dir = os.path.dirname(os.path.abspath(__file__))
clusters_path = os.path.join(current_dir, "data", "protein_clusters_k3.csv")

if os.path.exists(clusters_path):
    clusters = pd.read_csv(clusters_path).set_index("Protein")['Cluster']
else:
    # Fallback or empty if file is missing (prevents crash on import, allows debugging)
    print(f"Warning: {clusters_path} not found.")
    clusters = pd.Series(dtype=int)


def feature_selection_anova_only(X, y, protein_name, percentile=5, n_samples=None):
    """
    Perform leave-one-out (LOO) feature selection for a single protein using ONLY ANOVA.
    
    For each left-out sample, the training set is used to compute:
      - 'anova': Uses f_regression to obtain F scores and p-values.
      
    The threshold is defined as the given percentile of the p-values computed over the training set.
    Only features with p-values below or equal to that threshold are selected.

    Returns:
        pd.DataFrame: A DataFrame indexed by the left-out sample with a single column 'anova'.
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
        
        # --- ANOVA Feature Selection ONLY ---
        # f_regression returns F scores and p-values.
        F_scores, pvals_anova = f_regression(X_train, y_train)
        
        # Calculate threshold based on percentile
        threshold_anova = np.percentile(pvals_anova, percentile)
        
        # Select features where p-value is <= threshold
        selected_anova = {feat: F_scores[i] 
                          for i, feat in enumerate(features) 
                          if pvals_anova[i] <= threshold_anova}
        
        # Append only the ANOVA results
        results.append({
            'LOO_sample': sample,
            'anova': selected_anova
        })
    
    # Create DataFrame with left-out sample as index
    results_df = pd.DataFrame(results).set_index('LOO_sample')
    return results_df


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


def _load_clusters_map(clusters_csv_path: str) -> pd.Series:
    """
    Load a CSV with columns: Protein, Cluster
    Returns a Series indexed by Protein with integer cluster ids.
    """
    if not os.path.isfile(clusters_csv_path):
        raise FileNotFoundError(
            f"protein_clusters_k3.csv not found at: {clusters_csv_path}"
        )
    df = pd.read_csv(clusters_csv_path)
    if "Protein" not in df.columns or "Cluster" not in df.columns:
        raise ValueError(
            f"Expected columns ['Protein','Cluster'] in {clusters_csv_path}"
        )
    return df.set_index("Protein")["Cluster"]


def naive_cluster_feature_selection_indices(
    X: pd.DataFrame,
    y: pd.DataFrame,
    protein_name: str,
    *,
    clusters_csv_path: str,
    n_samples: Optional[int] = None,
) -> dict:
    """
    Highly optimized LOO feature selection using vectorized ANOVA.
    
    OPTIMIZATION NOTES:
    - Replaced slow 'pearsonr' loop with sklearn's 'f_regression'.
    - f_regression is vectorized (calculates all genes at once).
    - It ranks features identically to Pearson correlation |r|.
    """
    if protein_name not in y.columns:
        raise ValueError(f"Protein {protein_name} not found in y columns.")

    # 1. Resolve cluster and N features
    clusters_map = _load_clusters_map(clusters_csv_path)
    if protein_name not in clusters_map.index:
        raise ValueError(f"Protein {protein_name} not present in clusters file.")
    
    cluster_id = int(clusters_map.loc[protein_name])
    if cluster_id not in MAX_FEATURES:
        raise ValueError(f"Cluster id {cluster_id} not in MAX_FEATURES.")
        
    keep_top = int(MAX_FEATURES[cluster_id])

    # 2. Prepare Data (Slicing)
    sample_idx = X.index if n_samples is None else X.index[:n_samples]
    X_subset = X.loc[sample_idx]
    y_subset = y.loc[sample_idx, protein_name]

    # Convert to numpy arrays ONCE to avoid DataFrame overhead inside loop
    X_values = X_subset.values  # shape: (n_samples, n_genes)
    y_values = y_subset.values  # shape: (n_samples,)
    
    n_genes = X_values.shape[1]
    samples: Sequence[str] = list(sample_idx)

    # 3. Output buffers
    indptr = [0]
    indices = []

    # 4. Optimized LOO Loop
    # We loop over indices (i) to mask rows using boolean indexing or slicing
    n_rows = X_values.shape[0]
    
    for i in range(n_rows):
        # Create a boolean mask to select everything EXCEPT the current sample i
        # This is faster than dropping rows with pandas
        mask = np.ones(n_rows, dtype=bool)
        mask[i] = False
        
        X_train = X_values[mask]
        y_train = y_values[mask]

        # --- SPEED OPTIMIZATION HERE ---
        # Instead of looping 20,000 times with pearsonr, we run ONE f_regression call.
        # f_regression returns (F-scores, p-values). Higher F-score == Stronger Correlation.
        f_scores, _ = f_regression(X_train, y_train)
        
        # Handle potential NaNs (e.g. constant features) by setting them to -1
        np.nan_to_num(f_scores, copy=False, nan=-1.0)

        # Select top-N indices efficiently
        k = min(keep_top, n_genes)
        if k <= 0:
            chosen_idx = np.empty(0, dtype=np.int64)
        else:
            # np.argpartition puts the top K elements at the end of the array (unsorted)
            # This is O(n), much faster than full sort O(n log n)
            part_indices = np.argpartition(f_scores, n_genes - k)[n_genes - k:]
            
            # Optional: Sort these top K to have them in descending order of importance
            # (If order doesn't matter for your downstream task, you can skip this sort)
            top_scores = f_scores[part_indices]
            sorted_order = np.argsort(top_scores)[::-1]
            chosen_idx = part_indices[sorted_order]

        # Append to efficient list
        indices.extend(chosen_idx)
        indptr.append(len(indices))

    return {
        "samples": np.asarray(samples),
        "indptr": np.asarray(indptr, dtype=np.uint32),
        "indices": np.asarray(indices, dtype=np.uint32),
        "n_genes": n_genes,
    }