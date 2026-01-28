"""
Calculate metrics found in scTranslator and create bar chart for elasticnet predictions vs proteomics ground truth.
Metrics: MSE, RMSE, MAE, Pearson Correlation, Spearman Correlation, and Cosine Similarity
"""
import joblib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

print("Loading elasticnet results...")

# Load elasticnet predictions from HDF5 - structure is:
# Each protein has columns: ['protein', 'LOO sample', 'true', 'predicted', 'params_dict']

predictions_list = []
ground_truth_list = []

#make sure elasticnet_results.h5 is in your folder!
with pd.HDFStore('elasticnet_results.h5', mode='r') as store:
    protein_keys = store.keys()
    print(f"Number of proteins: {len(protein_keys)}")
    
    # Load all proteins
    for key in protein_keys:
        df = store[key]
        protein_name = key.strip('/')
        
        # Extract predictions and ground truth for this protein
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
#debug
print("Creating DataFrames...")

# Create DataFrames
pred_df = pd.DataFrame(predictions_list)
truth_df = pd.DataFrame(ground_truth_list)

# Get the predictions vs ground truth protein abundance
predictions = pred_df.pivot(index='sample', columns='protein', values='predicted')
ground_truth = truth_df.pivot(index='sample', columns='protein', values='true')

#make sure shape is correct
print(f"\nPredictions shape: {predictions.shape}")
print(f"Ground truth shape: {ground_truth.shape}")


# just to make sure they are correctly aligned
assert list(predictions.index) == list(ground_truth.index), "Sample indices don't match!"
assert list(predictions.columns) == list(ground_truth.columns), "Protein columns don't match!"


print("Calculating metrics per cell line...")

# Define metric functions
def cosine_similarity(y_true, y_pred, eps=1e-8):
    """Cosine similarity between two vectors"""
    dot_product = np.dot(y_true, y_pred)
    norm_true = np.linalg.norm(y_true)
    norm_pred = np.linalg.norm(y_pred)
    # this is how they do the metric for cosine similarity in scTranslator below
    return dot_product / max(norm_true * norm_pred, eps)

def pearson_correlation(y_true, y_pred):
    """Pearson correlation coefficient as defined in scTranslator"""
    if np.std(y_true) == 0 or np.std(y_pred) == 0:
        return np.nan
    return np.corrcoef(y_true, y_pred)[0, 1]

def spearman_correlation(y_true, y_pred):
    """Spearman correlation coefficient"""
    return stats.spearmanr(y_true, y_pred)[0]

def mse(y_true, y_pred):
    """Mean Squared Error"""
    return np.mean((y_pred - y_true) ** 2)

def rmse(y_true, y_pred):
    """Root Mean Squared Error"""
    return np.sqrt(mse(y_true, y_pred))

def mae(y_true, y_pred):
    """Mean Absolute Error"""
    return np.mean(np.abs(y_true - y_pred))

# Calculate metrics for each cell line
metrics_per_sample = {
    'MSE': [],
    'RMSE': [],
    'MAE': [],
    'Pearson': [],
    'Spearman': [],
    'Cosine': []
}

samples = predictions.index.tolist()
for sample in samples:
    y_true = ground_truth.loc[sample].values.astype(float)
    y_pred = predictions.loc[sample].values.astype(float)
    
    # Handle any NaN values
    mask = ~(np.isnan(y_true) | np.isnan(y_pred))
    y_true_clean = y_true[mask]
    y_pred_clean = y_pred[mask]
    
    if len(y_true_clean) > 0:
        metrics_per_sample['MSE'].append(mse(y_true_clean, y_pred_clean))
        metrics_per_sample['RMSE'].append(rmse(y_true_clean, y_pred_clean))
        metrics_per_sample['MAE'].append(mae(y_true_clean, y_pred_clean))
        metrics_per_sample['Pearson'].append(pearson_correlation(y_true_clean, y_pred_clean))
        metrics_per_sample['Spearman'].append(spearman_correlation(y_true_clean, y_pred_clean))
        metrics_per_sample['Cosine'].append(cosine_similarity(y_true_clean, y_pred_clean))

# Calculate mean metrics
mean_metrics = {k: np.nanmean(v) for k, v in metrics_per_sample.items()}

print("\nMean Metrics across all cell lines:")
for metric, value in mean_metrics.items():
    print(f"  {metric}: {value:.4f}")


plt.style.use('seaborn-v0_8-whitegrid')

# Create the bar chart 
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.patch.set_facecolor('#F8F9FA')

# Main title with styling
fig.suptitle('ElasticNet Benchmark: Overall Performance Metrics', 
             fontsize=20, fontweight='bold', color='#2C3E50', y=0.98)


colors = [
    '#E74C3C',  # MSE - Red
    '#E67E22',  # RMSE - Orange  
    '#F39C12',  # MAE - Yellow/Gold
    '#27AE60',  # Pearson - Green
    '#2980B9',  # Spearman - Blue
    '#9B59B6',  # Cosine - Purple
]

# Plot configurations
plots = [
    ('Mean Squared Error', 'MSE', mean_metrics['MSE'], colors[0], ""),
    ('Root MSE', 'RMSE', mean_metrics['RMSE'], colors[1], ""),
    ('Mean Absolute Error', 'MAE', mean_metrics['MAE'], colors[2], ""),
    ('Pearson Correlation', 'PCC', mean_metrics['Pearson'], colors[3], ""),
    ('Spearman Correlation', 'SCC', mean_metrics['Spearman'], colors[4], ""),
    ('Cosine Similarity', 'COS', mean_metrics['Cosine'], colors[5], "")
]

for idx, (title, ylabel, value, color, note) in enumerate(plots):
    row = idx // 3
    col = idx % 3
    ax = axes[row, col]
    
    # Set background
    ax.set_facecolor('#FFFFFF')
    
    # Create bar with gradient effect using edge color
    bar = ax.bar(['ElasticNet'], [value], color=color, width=0.6, 
                  edgecolor='#2C3E50', linewidth=2, alpha=0.9)
    
    # Add shadow effect
    ax.bar(['ElasticNet'], [value], color='gray', width=0.6, alpha=0.1, 
           bottom=-value*0.02 if value > 0 else value*0.02)
    
    # Add value label on top of bar with better styling
    for rect in bar:
        height = rect.get_height()
        ax.annotate(f'{value:.4f}',
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 8),
                    textcoords="offset points",
                    ha='center', va='bottom',
                    fontsize=14, fontweight='bold', color='#2C3E50',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                             edgecolor=color, alpha=0.8))
    
    # Set title with color accent
    ax.set_title(title, fontsize=14, fontweight='bold', color='#2C3E50', pad=15)
    
    # Style the ylabel - black color
    ax.set_ylabel(ylabel, fontsize=12, fontweight='bold', color='#2C3E50')
    
    # Adjust y-axis limits for better visualization
    if value >= 0:
        ax.set_ylim(0, value * 1.35 if value > 0 else 1)
    else:
        ax.set_ylim(value * 1.35, abs(value) * 0.35)
    
    # Style axes
    ax.tick_params(axis='x', labelsize=12, colors='#2C3E50')
    ax.tick_params(axis='y', labelsize=10, colors='#7F8C8D')
    
    # Add horizontal grid lines
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)
    ax.set_axisbelow(True)
    
    # Style spines
    for spine in ax.spines.values():
        spine.set_color('#BDC3C7')
        spine.set_linewidth(1.5)
    
    # Remove x-axis spine
    ax.spines['bottom'].set_visible(False)

plt.tight_layout(rect=[0, 0.02, 1, 0.95])
plt.savefig('elasticnet_metrics.png', dpi=200, bbox_inches='tight', facecolor=fig.get_facecolor())
print("\nPlot saved as 'elasticnet_metrics.png'")
