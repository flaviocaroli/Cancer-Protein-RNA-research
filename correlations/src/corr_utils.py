

### FUNCTIONS USED IN "corr_analysis.ipynb" ###


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os


from scipy.spatial.distance import pdist, squareform
import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import dendrogram

from sklearn.metrics import silhouette_score

# Colormap (used in heatmaps)
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap



def silhouette_method(pearson_clean, linkage_matrix, k_range=(2, 15)):
    """Computes the best number of clusters using the Silhouette Score."""
    best_k = k_range[0]
    best_score = -1

    for k in range(k_range[0], k_range[1]):
        labels = sch.fcluster(linkage_matrix, k, criterion='maxclust')
        score = silhouette_score(pearson_clean, labels, metric='euclidean')

        if score > best_score:
            best_score = score
            best_k = k

    return best_k


def elbow_method(pearson_clean, linkage_matrix, k_range=(2, 15)):
    """Computes the Elbow Method to find optimal clusters."""
    wcss = []
    for k in range(k_range[0], k_range[1]):
        labels = sch.fcluster(linkage_matrix, k, criterion='maxclust')
        cluster_var = sum(np.var(pearson_clean[labels == i], axis=0).sum() for i in np.unique(labels))
        wcss.append(cluster_var)

    # Find the "elbow" (biggest drop in variance)
    diffs = np.diff(wcss)  # First derivative
    elbow_k = k_range[0] + np.argmin(diffs)  # Best k is where WCSS drops the most

    # Plot the elbow curve
    plt.figure(figsize=(8, 4))
    plt.plot(range(k_range[0], k_range[1]), wcss, marker='o', label="WCSS")
    plt.axvline(x=elbow_k, color='r', linestyle='--', label=f"Elbow at k={elbow_k}")
    plt.xlabel("Number of Clusters")
    plt.ylabel("Within-cluster Variance (WCSS)")
    plt.title("Elbow Method for Optimal Clusters")
    plt.legend()
    plt.show()

    return elbow_k


# Calculate the CDF for consensus values (useful for choosing k)
def plot_consensus_cdf(consensus_matrices, figsize=(10, 6)):
    """
    Plot the empirical CDF of consensus values for each k.
    
    Parameters:
    consensus_matrices: List or array of consensus matrices for each k
    figsize: Size of the figure (width, height)
    """
    plt.figure(figsize=figsize)
    
    # For each consensus matrix
    for k_idx, consensus in enumerate(consensus_matrices):
        k = k_idx + 2  # k starts at 2
        
        # Flatten the upper triangular part (excluding diagonal)
        mask = np.triu(np.ones(consensus.shape), k=1).astype(bool)
        values = consensus[mask]
        
        # Sort values for CDF
        values_sorted = np.sort(values)
        
        # Calculate empirical CDF
        cdf = np.arange(1, len(values_sorted) + 1) / len(values_sorted)
        
        # Plot CDF
        plt.plot(values_sorted, cdf, label=f'k={k}')
    
    plt.xlabel('Consensus Value')
    plt.ylabel('CDF')
    plt.title('Empirical CDF of Consensus Values', fontsize=16)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    return plt

# Function to compute Dunn Index
def dunn_index(data, labels):
    distances = squareform(pdist(data, metric="euclidean"))
    unique_clusters = np.unique(labels)
    
    intercluster_distances = np.min([
        np.min(distances[labels == i][:, labels == j]) 
        for i in unique_clusters for j in unique_clusters if i != j
    ])
    
    intracluster_distances = np.max([
        np.max(distances[labels == i][:, labels == i]) 
        for i in unique_clusters
    ])
    
    return intercluster_distances / intracluster_distances


def visualize_all_k(consensus_matrices, k_range=range(2, 7)):
    n = len(k_range)
    cols = 3  # You can adjust this based on your preference
    rows = (n + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))

    for ax, k in zip(axes.flat, k_range):
        matrix = consensus_matrices[k]
        im = ax.imshow(matrix, cmap='viridis', vmin=0, vmax=1)
        ax.set_title(f"Consensus Matrix (k={k})")
        ax.axis("off")

    # Turn off any unused subplots
    for i in range(n, rows * cols):
        axes.flat[i].axis("off")

    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.6)
    plt.suptitle("Consensus Matrices for k = 2 to 6", fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()


# Example usage - visualize each consensus matrix
def visualize_all_consensus_matrices(consensus_matrices):
    """
    Visualize all consensus matrices.
    
    Parameters:
    consensus_matrices: Array of consensus matrices (k-2, n_proteins, n_proteins)
    """
    for k_idx, consensus in enumerate(consensus_matrices):
        k = k_idx + 2  # k starts at 2
        plot_consensus_matrix(consensus, k)
        plt.savefig(f'consensus_matrix_k{k}.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # Plot CDFs
    plot_consensus_cdf(consensus_matrices)
    plt.savefig('consensus_cdf.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Calculate and plot stability metrics
    metrics = calculate_stability_metrics(consensus_matrices)
    
    # Plot PAC values (lower is better)
    plt.figure(figsize=(10, 6))
    plt.plot(metrics['k'], metrics['PAC'], 'o-', linewidth=2)
    plt.xlabel('Number of Clusters (k)')
    plt.ylabel('PAC Value')
    plt.title('Proportion of Ambiguous Clustering (lower is better)')
    plt.grid(True, alpha=0.3)
    plt.xticks(metrics['k'])
    plt.savefig('pac_values.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Print metrics table
    print("Cluster Stability Metrics:")
    print(metrics)
    
    return metrics


# If memory is an issue, you can visualize one at a time:
def visualize_one_k(consensus_matrices, k):
    """Visualize consensus matrix for a specific k"""
    consensus = consensus_matrices[k-2]  # k starts at 2, so k=2 is at index 0
    plot_consensus_matrix(consensus, k)
    plt.show()

# Example for k=3:
# visualize_one_k(consensus_matrices, 3)

def plot_dendrogram(linkage_matrix, threshold=None):
    """Plots the dendrogram and suggests the number of clusters based on visual gap."""
    plt.figure(figsize=(10, 5))
    dendrogram(linkage_matrix)
    if threshold:
        plt.axhline(y=threshold, color='r', linestyle='--')
    plt.xlabel("Proteins")
    plt.ylabel("Cluster Distance")
    plt.title("Dendrogram for Hierarchical Clustering")
    plt.show()



# Function to calculate cluster stability metrics
def calculate_stability_metrics(consensus_matrices):
    """
    Calculate stability metrics for different values of k.
    
    Parameters:
    consensus_matrices: List or array of consensus matrices for each k
    
    Returns:
    DataFrame with stability metrics for each k
    """
    results = []
    
    for k_idx, consensus in enumerate(consensus_matrices):
        k = k_idx + 2  # k starts at 2
        
        # Calculate PAC (Proportion of Ambiguous Clustering)
        # PAC is defined as the proportion of sample pairs with consensus values in [0.1, 0.9]
        mask = np.triu(np.ones(consensus.shape), k=1).astype(bool)
        values = consensus[mask]
        pac = np.sum((values >= 0.1) & (values <= 0.9)) / len(values)
        
        # Calculate mean consensus value
        mean_consensus = np.mean(values)
        
        # Calculate median consensus value
        median_consensus = np.median(values)
        
        # Calculate dispersion metric (based on entropy)
        # For perfectly stable clustering, we want all values to be 0 or 1
        values_adj = np.maximum(values, 1e-10)  # Avoid log(0)
        entropy = -np.mean(values_adj * np.log2(values_adj) + 
                           (1-values_adj) * np.log2(np.maximum(1-values_adj, 1e-10)))
        
        results.append({
            'k': k,
            'PAC': pac,
            'Mean_Consensus': mean_consensus,
            'Median_Consensus': median_consensus,
            'Entropy': entropy
        })
    
    return pd.DataFrame(results)


# Function to plot a consensus matrix with ordered clusters
def plot_consensus_matrix(consensus_matrix, k, figsize=(10, 8)):
    """
    Plot a consensus matrix with hierarchical clustering to order samples.
    
    Parameters:
    consensus_matrix: The consensus matrix to plot
    k: The number of clusters this consensus matrix represents
    figsize: Size of the figure (width, height)
    """
    # Hierarchical clustering to order the matrix
    
    # Calculate distances and linkage
    row_linkage = sch.linkage(consensus_matrix, method='ward')
    
    # Get cluster assignments for k clusters
    cluster_labels = sch.fcluster(row_linkage, k, criterion='maxclust')
    
    # Get indices for sorting samples by cluster
    idx_sort = np.argsort(cluster_labels)
    
    # Reorder the consensus matrix
    ordered_matrix = consensus_matrix[idx_sort, :][:, idx_sort]
    
    # Create figure
    plt.figure(figsize=figsize)
    
    # Create a custom colormap from white to dark blue
    colors = [(1, 1, 1), (0.8, 0.8, 1), (0.6, 0.6, 1), (0.4, 0.4, 1), (0, 0, 0.8)]
    cmap_name = 'white_blue'
    cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=100)
    # Plot heatmap
    sns.heatmap(ordered_matrix, cmap=cm, vmin=0, vmax=1, 
                cbar_kws={'label': 'Consensus Value'})
    
    # Add cluster boundaries if k > 1
    if k > 1:
        # Find boundaries between clusters in sorted matrix
        boundaries = []
        prev_label = cluster_labels[idx_sort][0]
        
        for i, label in enumerate(cluster_labels[idx_sort][1:], 1):
            if label != prev_label:
                boundaries.append(i)
                prev_label = label
        
        # Add lines to mark cluster boundaries
        for b in boundaries:
            plt.axhline(y=b, color='red', linestyle='-', linewidth=1)
            plt.axvline(x=b, color='red', linestyle='-', linewidth=1)
    
    plt.title(f'Consensus Matrix for k={k}', fontsize=16)
    plt.tight_layout()
    return plt

# Function to read cluster assignments from CSV files
def read_cluster_files(directory):
    results = {}
    if not os.path.exists(directory):
        print(f"Directory {directory} does not exist!")
        return results
        
    for filename in os.listdir(directory):
        if filename.endswith('.csv'):
            # Extract k value from filename
            if 'k' in filename:
                k = int(filename.split('k')[1].split('.')[0])
                file_path = os.path.join(directory, filename)
                df = pd.read_csv(file_path)
                results[k] = df
    return results

# Function to compute overlap between two cluster assignments
def compute_overlap(df1, df2, df1_name="Dataset1", df2_name="Dataset2", output_file="cluster_comparison.csv"):
    # Merge the two dataframes on protein ID
    merged = pd.merge(df1, df2, on='Protein', suffixes=('_1', '_2'))
    
    # Rename columns for clarity in the output file
    merged = merged.rename(columns={
        'Cluster_1': f'{df1_name}_Cluster',
        'Cluster_2': f'{df2_name}_Cluster'
    })
    
    # Create contingency table
    contingency = pd.crosstab(merged[f'{df1_name}_Cluster'], merged[f'{df2_name}_Cluster'])
    
    # Get total number of proteins
    total_proteins = len(merged)
    
    # Calculate cluster overlap counts
    overlap_count = 0
    for i in range(len(contingency)):
        max_overlap = contingency.iloc[i].max()
        overlap_count += max_overlap
    
    # Calculate overlap percentage
    overlap_percentage = (overlap_count / total_proteins) * 100
    
    # Sort the merged dataframe:
    # First proteins with matching clusters from cluster 1
    # Then proteins with matching clusters from cluster 2
    # Finally proteins with mismatching clusters
    
    # Create a column to identify if clusters match
    merged['Clusters_Match'] = merged[f'{df1_name}_Cluster'] == merged[f'{df2_name}_Cluster']
    
    # Sort the dataframe
    sorted_merged = pd.concat([
        # Matching with cluster 1
        merged[(merged['Clusters_Match']) & (merged[f'{df1_name}_Cluster'] == 1)],
        # Matching with cluster 2
        merged[(merged['Clusters_Match']) & (merged[f'{df1_name}_Cluster'] == 2)],
        # Mismatching clusters
        merged[~merged['Clusters_Match']]
    ])
    
    # Save to file (drop the match column as it's just for sorting)
    sorted_merged.drop('Clusters_Match', axis=1).to_csv(output_file, index=False)
    
    # Return both percentage and actual count of overlapping proteins
    return overlap_count, overlap_percentage, contingency

# Function to compute retained genes per threshold
def compute_gene_retention(df_cluster, percentiles):
    vals = df_cluster.values.flatten()
    vals = vals[~np.isnan(vals)]
    abs_vals = np.abs(vals)
    thresholds = {p: np.percentile(abs_vals, p) for p in percentiles}

    retained = {}
    for p, thr in thresholds.items():
        # For each gene (column), check if **any** protein has correlation above threshold
        retain_mask = (np.abs(df_cluster) >= thr).any(axis=0)
        retained[p] = retain_mask.sum()
    return thresholds, retained

# Function to compute consensus across three cluster assignments
def compute_consensus(df1, df2, df3):
    # Merge all three dataframes on protein ID
    merged = pd.merge(df1, df2, on='Protein', suffixes=('_1', '_2'))
    merged = pd.merge(merged, df3, on='Protein')
    merged.rename(columns={'Cluster': 'Cluster_3'}, inplace=True)
    
    # Count how many proteins have the same cluster assignment in all methods
    all_match = (merged['Cluster_1'] == merged['Cluster_2']) & (merged['Cluster_2'] == merged['Cluster_3'])
    total_consensus = all_match.sum()
    
    # Count pairwise agreements
    pair1_match = (merged['Cluster_1'] == merged['Cluster_2']).sum()
    pair2_match = (merged['Cluster_2'] == merged['Cluster_3']).sum()
    pair3_match = (merged['Cluster_1'] == merged['Cluster_3']).sum()
    
    total_proteins = len(merged)
    
    # Calculate percentages
    all_consensus_pct = (total_consensus / total_proteins) * 100
    pair1_pct = (pair1_match / total_proteins) * 100
    pair2_pct = (pair2_match / total_proteins) * 100
    pair3_pct = (pair3_match / total_proteins) * 100
    
    return {
        'all_consensus': all_consensus_pct,
        'pg_pp_overlap': pair1_pct,
        'pp_pe_overlap': pair2_pct,
        'pg_pe_overlap': pair3_pct,
        'total_proteins': total_proteins,
        'full_consensus_count': total_consensus
    }