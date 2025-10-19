# correlations


## Files
This folder contains three analysis notebooks. Run them **in order** (01 → 02 → 03).  
Filenames are kept as-is (typos in names are intentional to match the originals).

### 01. pearson_anlaysis.ipynb
- **What it does:** Computes pairwise **Pearson correlations** and basic stats.
- **Typical outputs:** Correlation matrix (CSV/Parquet), ranked pairs, and a heatmap (optionally clustered).
- **Use cases:** Spot multicollinearity, identify strongly related features/samples.

### 02. cluster_anlaysis.ipynb
- **What it does:** Performs **unsupervised clustering** (e.g., hierarchical / k-means) on the selected feature/correlation space.
- **Typical outputs:** Cluster labels, dendrogram or elbow/silhouette diagnostics, and 2D projections (PCA/UMAP) for visualization.
- **Use cases:** Choose a sensible number of clusters and inspect cluster structure/quality.

### 03. prot_clustering.ipynb
- **What it does:** **Protein-level clustering** and grouping based on similarity/correlation (may include graph-based views).
- **Typical outputs:** Protein clusters/groups as tables, summary plots (heatmaps/network), and export of selected clusters.
---

### Quick start
1. Install dependencies (from repo root):
   ```bash
   pip install -r requirements.txt
   ```
2. Open each notebook and set input/output paths in the first configuration cell (defaults usually point to a `data/` and `figures/` folder).
3. Run all cells top-to-bottom.

### Expected layout
```
.
├─ 01. pearson_anlaysis.ipynb
├─ 02. cluster_anlaysis.ipynb
├─ 03. prot_clustering.ipynb
└─ README.md
```
