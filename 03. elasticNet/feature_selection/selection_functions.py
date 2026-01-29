import numpy as np
import pandas as pd
from sklearn.feature_selection import f_regression
from scipy.stats import pearsonr, spearmanr
from typing import Dict, Optional, Sequence

# Per-cluster caps for number of selected features
MAX_FEATURES: Dict[int, int] = {
    1: 8600,
    2: 10400,
    3: 6700,
}


clusters = pd.read_csv("protein_clusters_k3.csv").set_index("Protein")['Cluster']


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
    Naive, cluster-aware LOO feature selection that returns a **compact** index layout.

    For each left-out sample:
      - Compute Pearson r between y_protein(train) and each gene X(train, j)
      - Keep top-N genes, where N depends on the protein's cluster:
            N = MAX_FEATURES[ clusters[protein_name] ]
      - Store only the global **gene indices** (not strings) for memory efficiency.

    Returns a dict suitable for np.savez_compressed:
        {
          "samples": np.ndarray[str]     # sample ids in order
          "indptr":  np.ndarray[uint32]  # CSR pointers (len = n_samples + 1)
          "indices": np.ndarray[uint32]  # concatenated selected gene indices
          "n_genes": int                 # total gene count (len(X.columns))
        }
    """
    if protein_name not in y.columns:
        raise ValueError(f"Protein {protein_name} not found in y columns.")

    # Resolve cluster and cap
    clusters = _load_clusters_map(clusters_csv_path)
    if protein_name not in clusters.index:
        raise ValueError(
            f"Protein {protein_name} not present in clusters file: {clusters_csv_path}"
        )
    cluster_id = int(clusters.loc[protein_name])
    if cluster_id not in MAX_FEATURES:
        raise ValueError(
            f"Cluster id {cluster_id} for {protein_name} not in MAX_FEATURES mapping."
        )
    keep_top = int(MAX_FEATURES[cluster_id])

    # Restrict to a subset of samples if requested (keeps order)
    sample_idx = X.index if n_samples is None else X.index[:n_samples]
    X = X.loc[sample_idx]
    y_protein = y.loc[sample_idx, protein_name]

    n_genes = X.shape[1]
    samples: Sequence[str] = list(sample_idx)
    genes_list: Sequence[str] = list(X.columns)
    gene_to_idx = {g: i for i, g in enumerate(genes_list)}

    # Output buffers (CSR format)
    indptr = [0]
    indices: list[int] = []

    # LOO over samples
    for s in samples:
        X_train = X.drop(index=[s])
        y_train = y_protein.drop(index=[s])

        # Compute Pearson r for each gene; collect as (idx, r)
        # NOTE: leave as a simple loop for clarity; vectorization is possible but not required here.
        r_vals = np.empty(n_genes, dtype=np.float32)
        xvals = X_train.values  # (n_train, n_genes)
        yvals = y_train.values  # (n_train,)

        # Avoid repeated .values in loop
        for j, col in enumerate(genes_list):
            # pearsonr returns (r, p); handle NaNs by setting r = -inf so they sink
            try:
                r, _ = pearsonr(xvals[:, j], yvals)
                if np.isnan(r):
                    r = -np.inf
            except Exception:
                r = -np.inf
            r_vals[j] = r

        # Select top-N by r (descending). Use argpartition for speed.
        # Bound keep_top to n_genes to be safe.
        k = min(keep_top, n_genes)
        if k <= 0:
            chosen_idx = np.empty(0, dtype=np.int64)
        else:
            # largest k -> partition at n-k, then sort that slice descending
            part = np.argpartition(r_vals, n_genes - k)[n_genes - k:]
            # sort those by r descending to get stable order (optional)
            chosen_idx = part[np.argsort(r_vals[part])[::-1]]

        # Record global gene indices
        indices.extend(int(j) for j in chosen_idx)
        indptr.append(len(indices))

    return {
        "samples": np.asarray(samples),
        "indptr": np.asarray(indptr, dtype=np.uint32),
        "indices": np.asarray(indices, dtype=np.uint32),
        "n_genes": n_genes,
    }
