# correlations


## Files
This folder contains a computation folder and three analysis notebooks. In order to get data to use in notebook please run computation/corr_script.py or use the run_corr.sh shell script in an external cluster; corr_test.ipynb is a controllable test script for small data.

Once full_correlations are computed and saved to disk, the three notebooks can be used to analyze the results.
Run notebooks **in order** (01 → 02 → 03).  

### 01. pearson_anlaysis.ipynb
- **What it does:** Analyzes pairwise **Pearson correlations** between each protein and each feature, 
along with basic stats, clustering of the correlation matrix.
- **Outputs:** Correlation matrix CSV, ranked pairs, clusters and a heatmap.

### 02. cluster_anlaysis.ipynb
- **What it does:** Analyzes **unsupervised clustering** hierarchical and k-means on the selected feature/correlation space.
- **Outputs:** Cluster labels, dendrogram or elbow/silhouette diagnostics, and 2D projections PCA for visualization.

### 03. prot_clustering.ipynb
- **What it does:** **Protein-level clustering** and grouping based on similarity/correlation and includes graph-based views.
- **Outputs:** Protein clusters/groups as tables, summary plots (heatmaps/network), and export of selected clusters.
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
├─ 02. prot_clustering.ipynb
├─ 03. cluster_anlaysis.ipynb
└─ README.md
```
