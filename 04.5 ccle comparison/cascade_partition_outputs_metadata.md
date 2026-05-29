# Metadata: Cascade Protein Partition Output Files


## Output folder

All files are saved inside:

```text
protein_names_by_cascade_partition/
```

This folder is created inside the notebook variable `EXPORT_TABLE_DIR`.

For `procedural_ccle.ipynb`, this usually corresponds to the CCLE procedural table output folder.

For `9.procedural_cptac.ipynb`, this usually corresponds to the CPTAC procedural table output folder.

---

## Main files created

The export creates three main metadata/output tables:

```text
cascade_partition_proteins_long.csv
cascade_partition_proteins_grouped.csv
cascade_partition_proteins_files_manifest.csv
```

It also creates one subfolder:

```text
protein_lists_one_file_per_partition/
```

---

## 1. `cascade_partition_proteins_long.csv`

### Purpose

This is the main machine-readable protein assignment table.

Each row is one protein assigned to one specific partition.

A partition is defined by:

```text
dataset + analysis + threshold + stage/model
```

### Row interpretation

Each row means:

```text
This protein belongs to this dataset, this analysis, this threshold, and this cascade stage/model.
```

Example interpretation:

```text
CCLE | pearson_threshold_sweep | Pearson >= 0.50 | ENet | AAK1
```

This means that, in CCLE, at Pearson threshold 0.50, protein `AAK1` was assigned to the ENet stage of the cascade.

### Main columns

| Column | Description |
|---|---|
| `dataset` | Dataset name. Usually `CCLE`, `iTRAQ`, `TMT`, or `LabelFree`. |
| `analysis` | Name of the analysis/plot that produced the partition. |
| `source_variable` | Notebook variable used to create the table. Useful for reproducibility/debugging. |
| `threshold_label` | Human-readable threshold description, such as `Pearson >= 0.50`, `MSE <= 0.5`, or `All proteins / lowest MSE winner / no fixed threshold`. |
| `stage` | Cascade stage or model group, such as `Proxy Pearson`, `Proxy Spearman`, `ENet`, `XGB`, or `Still not good`. |
| `protein` | Protein name assigned to the partition. |

### Optional metric columns

Some columns appear only if they existed in the original assignment table.

| Column | Description |
|---|---|
| `threshold` | Numeric threshold used in the plot, when available. |
| `prediction_threshold` | Prediction-performance threshold, usually a Pearson threshold. |
| `corr_threshold` | Proxy-correlation eligibility threshold, if used. |
| `mse_percentile_threshold` | MSE percentile cutoff, if used. |
| `metric` | Generic metric value used for assignment. |
| `test_pearson_corr` | Pearson correlation on validation/test samples. |
| `test_mse` | Mean squared error on validation/test samples. |
| `test_cosine_similarity` | Cosine similarity on validation/test samples. |
| `best_mse` | Best MSE among compared models for that protein. |
| `proxy_mse` | Proxy model MSE. |
| `enet_mse` | ElasticNet MSE. |
| `xgb_mse` | XGBoost MSE. |
| `same_gene_corr_full` | Same-gene proxy correlation, if available. |
| `proxy_corr` | Proxy correlation value, if available. |
| `run_label` | Original model/run label from the notebook. |
| `model_name` | Model name, if available. |
| `proxy_type` | Proxy type, if available. |

### Best use

Use this file when filtering proteins in Python or R.

Example use cases:

- Find all proteins rescued by ENet at `Pearson >= 0.50`.
- Find all proteins still not good at a specific threshold.
- Compare protein partitions across `iTRAQ`, `TMT`, and `LabelFree`.
- Extract proteins belonging to one model/stage across all thresholds.

---

## 2. `cascade_partition_proteins_grouped.csv`

### Purpose

This is the human-readable summary table.

Each row corresponds to one graph segment or partition.

Instead of storing one protein per row, it groups all proteins for the same dataset, analysis, threshold, and stage into one row.

### Row interpretation

Each row means:

```text
For this dataset, analysis, threshold, and stage, these are the proteins in that partition.
```

Example interpretation:

```text
TMT | pearson_threshold_sweep | Pearson >= 0.50 | XGB | 116 proteins | protein1; protein2; ...
```

This means that, in TMT, at Pearson threshold 0.50, 116 proteins were assigned to the XGB stage.

### Columns

| Column | Description |
|---|---|
| `dataset` | Dataset name. |
| `analysis` | Analysis/plot name. |
| `source_variable` | Notebook variable used to produce the assignment. |
| `threshold_label` | Human-readable threshold description. |
| `stage` | Cascade stage/model group. |
| `n_proteins` | Number of unique proteins in the partition. |
| `protein_names` | Full semicolon-separated list of protein names. |
| `protein_names_preview` | Short preview of the protein list for notebook display. |

### Best use

Use this file for inspection, reporting, and checking the contents of each bar segment in the plots.

It is less convenient than the long table for filtering, but easier to read manually.

---

## 3. `cascade_partition_proteins_files_manifest.csv`

### Purpose

This file is an index of all separate protein-list CSV files saved in the folder:

```text
protein_lists_one_file_per_partition/
```

Each row tells us which file corresponds to one exact dataset, analysis, threshold, and stage.

### Row interpretation

Each row means:

```text
This partition was saved as this separate CSV file.
```

### Columns

| Column | Description |
|---|---|
| `dataset` | Dataset name. |
| `analysis` | Analysis/plot name. |
| `source_variable` | Notebook variable used. |
| `threshold_label` | Threshold description. |
| `stage` | Cascade stage/model group. |
| `n_proteins` | Number of proteins in the partition file. |
| `file_name` | Name of the separate CSV file. |
| `relative_path` | Path to the file relative to `protein_names_by_cascade_partition/`. |

### Best use

Use this file to quickly locate the correct separate CSV file for a given partition.

For example, it can tell us where the file is for:

```text
TMT + pearson_threshold_sweep + Pearson >= 0.50 + ENet
```

---

## 4. `protein_lists_one_file_per_partition/`

### Purpose

This folder contains one CSV file per exact partition.

Each file corresponds to one combination of:

```text
dataset + analysis + threshold + stage/model
```

### Example file names

```text
CCLE__pearson_threshold_sweep__Pearson_ge_0p50__ENet.csv
CCLE__pearson_threshold_sweep__Pearson_ge_0p50__Still_not_good.csv
TMT__mse_threshold_sweep__MSE_le_0p5__XGB.csv
LabelFree__mse_winner_partition__All_proteins_lowest_MSE_winner_no_fixed_threshold__ENet_lowest_MSE.csv
```

### Row interpretation

Each row is one protein belonging to that exact partition.

### Typical columns

| Column | Description |
|---|---|
| `dataset` | Dataset name. |
| `analysis` | Analysis/plot name. |
| `source_variable` | Notebook variable used. |
| `threshold_label` | Threshold description. |
| `stage` | Cascade stage/model group. |
| `protein` | Protein name. |
| Metric columns | Extra metric columns are included when available, such as `test_mse`, `test_pearson_corr`, or `test_cosine_similarity`. |

### Best use

Use this folder when we want a clean, ready-made protein list for one exact graph segment.

---

## Analysis names

The `analysis` column identifies which plot or partitioning logic generated the protein assignment.

### Pearson-based analyses

| Analysis name | Meaning |
|---|---|
| `main_pearson_cascade` | Main cascade at selected Pearson thresholds, usually 0.30, 0.50, and 0.70. |
| `pearson_threshold_sweep` | Full Pearson threshold sweep across many Pearson cutoffs. |
| `pearson_cascade_with_sctranslator` | Cascade including scTranslator after proxy, ENet, and XGB. |
| `sctranslator_rescue_only` | Proteins rescued only by scTranslator among previous failures. |

### MSE-based analyses

| Analysis name | Meaning |
|---|---|
| `mse_threshold_sweep` | Cascade using absolute MSE thresholds. A protein is good when `MSE <= threshold`. |
| `mse_percentile_threshold_sweep` | Cascade using MSE percentile thresholds. |
| `mse_winner_partition` | No fixed threshold. Assigns each protein to the model with the lowest MSE. |
| `mse_procedural_improvement_partition` | No fixed threshold. Assigns proteins according to whether ENet improves over proxy, then whether XGB improves over the previous best. |

### Other analyses

| Analysis name | Meaning |
|---|---|
| `proxy_correlation_threshold_sweep` | Analysis varying the same-gene proxy-correlation eligibility threshold. |
| `cosine_threshold_sweep` | Cascade using cosine similarity thresholds. A protein is good when `cosine >= threshold`. |

---

## Stage/model names

The `stage` column identifies the cascade group or model responsible for the assignment.

| Stage name | Meaning |
|---|---|
| `Proxy Pearson` | Protein reached the threshold using the same-gene Pearson proxy. |
| `Proxy Spearman` | Protein reached the threshold using the same-gene Spearman proxy. |
| `ENet` | Protein was not already assigned to proxy, but reached the threshold with ElasticNet. |
| `XGB` | Protein was not already assigned to proxy or ENet, but reached the threshold with XGBoost. |
| `Still not good` | Protein did not reach the threshold with any model in that cascade. |
| `Still not good after scTranslator` | Protein remained below threshold even after adding scTranslator. |
| `ENet lowest MSE` | ENet had the lowest MSE for that protein. |
| `XGB lowest MSE` | XGB had the lowest MSE for that protein. |
| `ENet improves over proxy` | ENet improved MSE compared with the proxy. |
| `XGB improves over previous best` | XGB improved MSE compared with the previous best method. |
| `ENet first valid MSE` | No valid proxy MSE was available, so ENet was the first valid MSE comparison. |

---

## Threshold interpretation

For Pearson and cosine analyses, higher values are better.

Therefore:

```text
Pearson >= threshold
Cosine >= threshold
```

For MSE analyses, lower values are better.

Therefore:

```text
MSE <= threshold
```

For MSE winner and MSE procedural improvement analyses, there is no fixed MSE threshold. These analyses partition all available proteins by relative MSE behavior.

Therefore, they are labeled as:

```text
All proteins / lowest MSE winner / no fixed threshold
All proteins / procedural MSE improvement / no fixed threshold
```

---

## Recommended usage

Use `cascade_partition_proteins_long.csv` for filtering and downstream analysis.

Use `cascade_partition_proteins_grouped.csv` for reading, reporting, and checking which proteins belong to each bar segment.

Use `cascade_partition_proteins_files_manifest.csv` to locate a specific separate protein-list file.

Use the files in `protein_lists_one_file_per_partition/` when a clean CSV is needed for one exact threshold/model group.

