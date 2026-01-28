"""
Create detailed visualization with histograms for ElasticNet predictions.
Shows distributions across cells and proteins.
"""
import joblib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


# Load elasticnet predictions from HDF5
predictions_list = []
ground_truth_list = []

#load the elasticnet results, make sure the file is in your folder!
with pd.HDFStore('elasticnet_results.h5', mode='r') as store:
    protein_keys = store.keys()
    print(f"Number of proteins: {len(protein_keys)}")
    
    for key in protein_keys:
        df = store[key]
        protein_name = key.strip('/')
        
        for _, row in df.iterrows():
            sample = row['LOO sample']
            true_val = row['true']
            pred_val = row['predicted']
            
            predictions_list.append({
                'sample': sample,
                'protein': protein_name,
                'predicted': pred_val
            })
            ground_truth_list.append({
                'sample': sample,
                'protein': protein_name,
                'true': true_val
            })

print("Creating DataFrames...")

pred_df = pd.DataFrame(predictions_list)
truth_df = pd.DataFrame(ground_truth_list)

predictions = pred_df.pivot(index='sample', columns='protein', values='predicted')
ground_truth = truth_df.pivot(index='sample', columns='protein', values='true')

print(f"Shape: {predictions.shape}")

# Define metric functions
def cosine_similarity(y_true, y_pred, eps=1e-8):
    dot_product = np.dot(y_true, y_pred)
    norm_true = np.linalg.norm(y_true)
    norm_pred = np.linalg.norm(y_pred)
    return dot_product / max(norm_true * norm_pred, eps)

def pearson_correlation(y_true, y_pred):
    if np.std(y_true) == 0 or np.std(y_pred) == 0:
        return np.nan
    return np.corrcoef(y_true, y_pred)[0, 1]

def mse(y_true, y_pred):
    return np.mean((y_pred - y_true) ** 2)

def mae(y_true, y_pred):
    return np.mean(np.abs(y_true - y_pred))

print("\nCalculating per-cell metrics...")

# Calculate per-cell metrics (across all proteins for each cell)
per_cell_metrics = {
    'sample': [],
    'cosine': [],
    'pearson': [],
    'mse': [],
    'mae': []
}

samples = predictions.index.tolist()
proteins = predictions.columns.tolist()

for sample in samples:
    y_true = ground_truth.loc[sample].values.astype(float)
    y_pred = predictions.loc[sample].values.astype(float)
    
    mask = ~(np.isnan(y_true) | np.isnan(y_pred))
    y_true_clean = y_true[mask]
    y_pred_clean = y_pred[mask]
    
    if len(y_true_clean) > 0:
        per_cell_metrics['sample'].append(sample)
        per_cell_metrics['cosine'].append(cosine_similarity(y_true_clean, y_pred_clean))
        per_cell_metrics['pearson'].append(pearson_correlation(y_true_clean, y_pred_clean))
        per_cell_metrics['mse'].append(mse(y_true_clean, y_pred_clean))
        per_cell_metrics['mae'].append(mae(y_true_clean, y_pred_clean))

per_cell_df = pd.DataFrame(per_cell_metrics)

print("\nCalculating per-protein metrics...")

# Calculate per-protein metrics (across all cells for each protein)
per_protein_metrics = {
    'protein': [],
    'pearson': []
}

for protein in proteins:
    y_true = ground_truth[protein].values.astype(float)
    y_pred = predictions[protein].values.astype(float)
    
    mask = ~(np.isnan(y_true) | np.isnan(y_pred))
    y_true_clean = y_true[mask]
    y_pred_clean = y_pred[mask]
    
    if len(y_true_clean) > 1:
        pcc = pearson_correlation(y_true_clean, y_pred_clean)
        per_protein_metrics['protein'].append(protein)
        per_protein_metrics['pearson'].append(pcc)

per_protein_df = pd.DataFrame(per_protein_metrics)

print(f"Per-cell metrics: {len(per_cell_df)} cells")
print(f"Per-protein metrics: {len(per_protein_df)} proteins")

# Get top 10 and bottom 10 proteins by PCC
per_protein_df_sorted = per_protein_df.sort_values('pearson', ascending=False)
top_10 = per_protein_df_sorted.head(10)
bottom_10 = per_protein_df_sorted.tail(10)

print("\nCreating visualization...")

# Create the figure
fig = plt.figure(figsize=(16, 10))

# Define grid for subplots
gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1], height_ratios=[1, 1], 
                       hspace=0.3, wspace=0.3)

# Colors
cell_cosine_color = '#5DADE2'  # Blue
cell_pcc_color = '#F5B041'     # Orange
cell_mse_color = '#58D68D'     # Green
protein_pcc_color = '#AF7AC5'  # Purple

# 1. Per-Cell Cosine Similarity (top left)
ax1 = fig.add_subplot(gs[0, 0])
ax1.hist(per_cell_df['cosine'], bins=30, color=cell_cosine_color, edgecolor='black', alpha=0.8)
mean_cosine = per_cell_df['cosine'].mean()
ax1.axvline(mean_cosine, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_cosine:.3f}')
ax1.set_xlabel('Cosine Similarity', fontsize=10)
ax1.set_ylabel('Count', fontsize=10)
ax1.set_title('Per-Cell Cosine Similarity', fontsize=12, fontweight='bold')
ax1.legend(loc='upper right')

# 2. Per-Cell PCC (top middle)
ax2 = fig.add_subplot(gs[0, 1])
ax2.hist(per_cell_df['pearson'], bins=30, color=cell_pcc_color, edgecolor='black', alpha=0.8)
mean_pcc = per_cell_df['pearson'].mean()
ax2.axvline(mean_pcc, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_pcc:.3f}')
ax2.set_xlabel('Pearson Correlation', fontsize=10)
ax2.set_ylabel('Count', fontsize=10)
ax2.set_title('Per-Cell PCC', fontsize=12, fontweight='bold')
ax2.legend(loc='upper right')

# 3. Per-Cell MSE (top right)
ax3 = fig.add_subplot(gs[0, 2])
ax3.hist(per_cell_df['mse'], bins=30, color=cell_mse_color, edgecolor='black', alpha=0.8)
mean_mse = per_cell_df['mse'].mean()
ax3.axvline(mean_mse, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_mse:.4f}')
ax3.set_xlabel('MSE', fontsize=10)
ax3.set_ylabel('Count', fontsize=10)
ax3.set_title('Per-Cell MSE', fontsize=12, fontweight='bold')
ax3.legend(loc='upper right')

# 4. Per-Protein PCC (bottom left)
ax4 = fig.add_subplot(gs[1, 0])
ax4.hist(per_protein_df['pearson'].dropna(), bins=30, color=protein_pcc_color, edgecolor='black', alpha=0.8)
mean_protein_pcc = per_protein_df['pearson'].mean()
ax4.axvline(mean_protein_pcc, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_protein_pcc:.3f}')
ax4.set_xlabel('Pearson Correlation', fontsize=10)
ax4.set_ylabel('Count', fontsize=10)
ax4.set_title('Per-Protein PCC (across cells)', fontsize=12, fontweight='bold')
ax4.legend(loc='upper right')

# 5. Top 10 / Bottom 10 Proteins by PCC (bottom middle)
ax5 = fig.add_subplot(gs[1, 1])

top_10_sorted = top_10.sort_values('pearson', ascending=True)  # best (IGF2BP1) will be at very bottom
# Bottom 10 sorted: worst at very top of chart  
bottom_10_sorted = bottom_10.sort_values('pearson', ascending=False)  # worst (DDX6) will be at very top
# Concatenate: top_10 first (bottom of chart), bottom_10 last (top of chart)
combined = pd.concat([top_10_sorted, bottom_10_sorted])
colors = ['#58D68D'] * 10 + ['#E74C3C'] * 10  # Green for top 10 (at bottom), Red for bottom 10 (at top)
y_pos = np.arange(len(combined))
ax5.barh(y_pos, combined['pearson'].values, color=colors, edgecolor='black', alpha=0.8)
ax5.set_yticks(y_pos)
ax5.set_yticklabels(combined['protein'].values, fontsize=8)
ax5.axvline(0, color='red', linestyle='-', linewidth=1)
ax5.set_xlabel('PCC', fontsize=10)
ax5.set_title('Top 10 / Bottom 10 Proteins by PCC', fontsize=12, fontweight='bold')

# 6. Summary Statistics (bottom right)
ax6 = fig.add_subplot(gs[1, 2])
ax6.axis('off')

# Calculate summary statistics
n_cells = len(per_cell_df)
n_proteins = len(per_protein_df)
positive_pcc_count = (per_protein_df['pearson'] > 0).sum()
positive_pcc_pct = positive_pcc_count / n_proteins * 100
high_pcc_03 = (per_protein_df['pearson'] > 0.3).sum()
high_pcc_05 = (per_protein_df['pearson'] > 0.5).sum()

summary_text = f"""Summary Statistics (ElasticNet)
{'='*40}

Per-Cell Metrics (n={n_cells} cells):

Cosine Similarity:
  Mean: {per_cell_df['cosine'].mean():.4f} +/- {per_cell_df['cosine'].std():.4f}
  Min: {per_cell_df['cosine'].min():.4f}
  Max: {per_cell_df['cosine'].max():.4f}

Pearson Correlation:
  Mean: {per_cell_df['pearson'].mean():.4f} +/- {per_cell_df['pearson'].std():.4f}
  Min: {per_cell_df['pearson'].min():.4f}
  Max: {per_cell_df['pearson'].max():.4f}

MSE:
  Mean: {per_cell_df['mse'].mean():.4f} +/- {per_cell_df['mse'].std():.4f}

MAE:
  Mean: {per_cell_df['mae'].mean():.4f} +/- {per_cell_df['mae'].std():.4f}

Per-Protein PCC (n={n_proteins} proteins):
  Mean: {per_protein_df['pearson'].mean():.4f} +/- {per_protein_df['pearson'].std():.4f}
  Positive: {positive_pcc_count} ({positive_pcc_pct:.1f}%)
  > 0.3:    {high_pcc_03} ({high_pcc_03/n_proteins*100:.1f}%)
  > 0.5:    {high_pcc_05} ({high_pcc_05/n_proteins*100:.1f}%)
"""

ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes, fontsize=9,
         verticalalignment='top', fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Main title
fig.suptitle('ElasticNet Benchmark Results', fontsize=16, fontweight='bold', y=0.98)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('elasticnet_detailed_metrics.png', dpi=150, bbox_inches='tight')
print("\nPlot saved as 'elasticnet_detailed_metrics.png'")
