"""Reusable validation utilities for the CPTAC/Collection A notebooks.

This module collects the long helper functions that were previously embedded in
multiple notebooks. The notebooks keep configuration, orchestration, and plots;
long reusable functions live here.
"""
from __future__ import annotations

import ast
import glob
import json
import lzma
import os
import pickle
import re
import sys
import textwrap
from pathlib import Path
from collections import Counter

import h5py
import joblib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from IPython.display import display
from scipy.stats import pearsonr, spearmanr, zscore
from scipy.spatial.distance import cosine
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from sklearn.metrics.pairwise import cosine_similarity as sklearn_cossim
from sklearn.preprocessing import StandardScaler, RobustScaler, QuantileTransformer
from tqdm.auto import tqdm

try:
    import mygene
except Exception:  # mygene is only needed by standardize_to_ccle_format.
    mygene = None

# Defaults used by functions that originated in notebooks. Notebooks may call
# configure_runtime_globals(...) to override any of these at runtime.
GROUPS = ["TMT", "iTRAQ", "LabelFree"]
TYPES = ["Protein", "RNA"]
OUTPUT_DIR = Path("../results_comparison_all_models")
PREDICTIONS_PARQUET = OUTPUT_DIR / "eval_predictions_long.parquet"
PREDICTIONS_CSV = OUTPUT_DIR / "eval_predictions_long.csv"
METRICS_PROTEIN_CSV = OUTPUT_DIR / "eval_metrics_by_protein.csv"
METRICS_TISSUE_CSV = OUTPUT_DIR / "eval_metrics_by_tissue.csv"
MODEL_DIRS = {}
EVALUATION_ORDER = []
FORCE_REEVALUATE = True
VAL_DATA_DIR = Path("../val_data")
datasets = {}
model_inventory = pd.DataFrame()
metrics_partition = pd.DataFrame()
SAVE_OUTPUTS = False
FIGURE_DIR = OUTPUT_DIR / "procedural_cptac_figures"
EXPORT_TABLE_DIR = OUTPUT_DIR / "procedural_cptac_tables"
TARGET_ENET_RUN_LABEL = "res_enet_zbig | zscore | fast | k5000"
METRIC_NAMES = {
    "test_pearson_corr": "Pearson",
    "test_cosine_similarity": "cosine similarity",
    "test_mse": "MSE",
}
PREDICTION_COLUMNS = [
    "dataset", "sample", "tissue", "protein", "y_true", "y_pred", "source_folder",
    "model_family", "normalization", "grid_preset", "select_k", "run_label", "h5_path",
]
METRIC_COLUMNS = [
    "dataset", "protein", "source_folder", "model_family", "normalization", "grid_preset",
    "select_k", "run_label", "h5_path", "test_mse", "test_pearson_corr",
    "test_cosine_similarity", "n_test_used", "selected_feature_count",
    "missing_feature_count", "missing_feature_fraction",
]
COVERAGE_COLUMNS = [
    "dataset", "protein", "source_folder", "run_label", "h5_path", "selected_feature_count",
    "missing_feature_count", "missing_feature_fraction",
]
BASE_STAGES = ["Proxy Pearson", "Proxy Spearman", "Proxy baseline", "ENet", "XGB"]
PROXY_STAGES = ["Proxy Pearson", "Proxy Spearman", "Proxy baseline"]
SC_STAGES = ["scTranslator 2M", "scTranslator 160K", "scTranslator 10K", "scTranslator"]
COLOR_MAP = {}
PROXY_CORR_CANDIDATES = []
USE_ABS_PROXY_CORR_FOR_ELIGIBILITY = True
PEARSON_GOOD_THRESHOLDS = [0.30, 0.50, 0.70]
PROXY_CORR_THRESHOLDS = [round(x, 2) for x in np.arange(0.00, 0.95, 0.05)]


# -----------------------------------------------------------------------------
# Consolidated CPTAC validation helpers used by 3.results_cptac.ipynb
# and 4.results_norms.ipynb
# -----------------------------------------------------------------------------

def configure_runtime_globals(**kwargs):
    """Inject notebook-level runtime objects into this module.

    Several original helper functions were authored in notebooks and refer to
    notebook globals such as OUTPUT_DIR, model_inventory, or metrics_partition.
    Refactored notebooks call this function immediately before those helpers are
    used, preserving the original behavior while keeping long code out of cells.
    """
    globals().update(kwargs)
    return kwargs


def set_module_globals(namespace=None, **kwargs):
    """Synchronize notebook variables into this module's global namespace.

    Newer cleaned notebooks call ``uv.set_module_globals(globals())`` after
    loading data or changing configuration. Older versions of ``utils_val.py``
    only exposed ``configure_runtime_globals``; this wrapper keeps both notebook
    styles compatible.

    Parameters
    ----------
    namespace : Mapping, optional
        Usually the notebook ``globals()`` dictionary. Non-private names are
        copied into this module so extracted helper functions can resolve paths,
        dataframes, and configuration values exactly as they did in-notebook.
    **kwargs
        Additional values to inject or override.
    """
    clean = {}
    if namespace is not None:
        try:
            items = dict(namespace).items()
        except Exception:
            items = []
        for name, value in items:
            if str(name).startswith("__"):
                continue
            clean[name] = value
    clean.update(kwargs)
    globals().update(clean)
    return None


def sync_module_globals(namespace=None, **kwargs):
    """Alias for set_module_globals, kept for readability in notebooks."""
    return set_module_globals(namespace, **kwargs)


def ensure_directory(path):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def _read_result_table(path):
    path = Path(path)
    return pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)


def _result_files(results_dir, pattern):
    """Files matching <pattern>.parquet or <pattern>.csv; parquet wins when both exist for the same stem."""
    results_dir = Path(results_dir)
    found = {}
    for suffix in (".csv", ".parquet"):
        for path in sorted(results_dir.glob(pattern + suffix)):
            found[path.stem] = path
    return [found[k] for k in sorted(found)]


def _load_result_tables(results_dir, pattern):
    files = _result_files(results_dir, pattern)
    rows = []
    for path in files:
        try:
            df = _read_result_table(path)
            df["source_file"] = path.name
            rows.append(df)
        except Exception as exc:
            print(f"Could not read {path}: {type(exc).__name__}: {exc}")
    if not rows:
        return pd.DataFrame(), files
    out = pd.concat(rows, ignore_index=True)
    if "protein" in out.columns:
        out = out.sort_values("protein").reset_index(drop=True)
    return out, files


def load_best_result_tables(results_dir):
    """Load all per-protein best models (<P>_best.parquet or <P>_best.csv) from a LOO result directory."""
    return _load_result_tables(results_dir, "*_best")


def load_batch_result_tables(results_dir):
    """Load all batch files (batch_*.parquet or batch_*.csv) from a LOO result directory."""
    return _load_result_tables(results_dir, "*batch*")


def compare_model_result_tables(best_df, batch_df):
    """Compare protein coverage between individual best files and batch files."""
    if best_df is None or batch_df is None or best_df.empty or batch_df.empty:
        return {
            "can_compare": False,
            "best_rows": 0 if best_df is None else len(best_df),
            "batch_rows": 0 if batch_df is None else len(batch_df),
            "message": "One or both result tables are empty.",
        }
    best_proteins = set(best_df.get("protein", pd.Series(dtype=str)).astype(str).unique())
    batch_proteins = set(batch_df.get("protein", pd.Series(dtype=str)).astype(str).unique())
    shape_equal = best_df.shape == batch_df.shape
    exact_equal = False
    if shape_equal:
        try:
            pd.testing.assert_frame_equal(
                best_df.drop(columns=["source_file"], errors="ignore").reset_index(drop=True),
                batch_df.drop(columns=["source_file"], errors="ignore").reset_index(drop=True),
                check_like=True,
            )
            exact_equal = True
        except AssertionError:
            exact_equal = False
    return {
        "can_compare": True,
        "best_rows": len(best_df),
        "batch_rows": len(batch_df),
        "best_proteins": len(best_proteins),
        "batch_proteins": len(batch_proteins),
        "missing_in_batch": sorted(best_proteins - batch_proteins),
        "missing_in_best": sorted(batch_proteins - best_proteins),
        "shape_equal": shape_equal,
        "exact_equal_after_sort": exact_equal,
    }


def load_validation_dataset(rna_path, prot_path):
    """Load a pair of validation RNA/protein joblib matrices and align samples."""
    rna_path = Path(rna_path)
    prot_path = Path(prot_path)
    if not rna_path.exists() or not prot_path.exists():
        raise FileNotFoundError(f"Missing validation files: {rna_path} or {prot_path}")
    rna_df = joblib.load(rna_path)
    prot_df = joblib.load(prot_path)
    if not isinstance(rna_df, pd.DataFrame) or not isinstance(prot_df, pd.DataFrame):
        raise TypeError("Expected both validation objects to be pandas DataFrames.")
    rna_df = rna_df.copy()
    prot_df = prot_df.copy()
    rna_df.index = rna_df.index.astype(str)
    prot_df.index = prot_df.index.astype(str)
    rna_df.columns = rna_df.columns.astype(str)
    prot_df.columns = prot_df.columns.astype(str)
    common = rna_df.index.intersection(prot_df.index)
    if len(common) > 0:
        rna_df = rna_df.loc[common].copy()
        prot_df = prot_df.loc[common].copy()
    return rna_df, prot_df


def percent_zeros_by_tissue(rna_df, tissue_func=None):
    if rna_df is None or len(rna_df) == 0:
        return pd.DataFrame(columns=["Tissue", "percent_zeros"])
    tissue_func = tissue_func or get_tissue
    tmp = rna_df.copy()
    tmp["Tissue"] = [tissue_func(x) for x in tmp.index]
    rows = []
    for tissue, sub in tmp.groupby("Tissue"):
        values = sub.drop(columns="Tissue")
        pct = 100.0 * (values == 0).sum().sum() / max(1, values.size)
        rows.append({"Tissue": tissue, "percent_zeros": pct, "samples": len(sub)})
    return pd.DataFrame(rows).sort_values("percent_zeros", ascending=False).reset_index(drop=True)


def parse_coefficients(value):
    """Coefficient dict from a dict (parquet), a JSON string or a python-literal string (csv)."""
    if isinstance(value, dict):
        parsed = value
    elif value is None or (isinstance(value, float) and np.isnan(value)):
        return {}
    else:
        text = str(value)
        try:
            parsed = json.loads(text)
        except Exception:
            try:
                parsed = ast.literal_eval(text)
            except Exception:
                return {}
    if not isinstance(parsed, dict):
        return {}
    return {str(k): float(v) for k, v in parsed.items() if v is not None and np.isfinite(float(v))}


def feature_coverage_for_models(model_df, rna_columns, n_preview=10):
    """Summarize selected-feature coverage for each best-model row."""
    if model_df is None or model_df.empty:
        return pd.DataFrame()
    available = set(map(str, rna_columns))
    rows = []
    for _, row in model_df.iterrows():
        coeffs = parse_coefficients(row.get("feat_coeffs_dict", {}))
        features = list(map(str, coeffs.keys()))
        missing = [x for x in features if x not in available]
        rows.append({
            "protein": row.get("protein"),
            "n_features": len(features),
            "n_present": len(features) - len(missing),
            "n_missing": len(missing),
            "missing_fraction": len(missing) / max(1, len(features)),
            "missing_preview": missing[:n_preview],
        })
    return pd.DataFrame(rows)


def predict_from_best_models(model_df, rna_df, prot_df, impute_missing_features=False):
    """Predict validation proteins from coefficient dictionaries and compute metrics.

    Parameters
    ----------
    model_df : DataFrame with columns protein and feat_coeffs_dict.
    rna_df : samples x genes validation transcriptomics.
    prot_df : samples x proteins validation proteomics.
    impute_missing_features : if False, models with missing RNA features are skipped;
        if True, missing features contribute zero.
    """
    if model_df is None or model_df.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    X_all = rna_df.fillna(0).copy()
    available = set(X_all.columns.astype(str))
    result_rows = []
    error_vectors = {}
    pred_vectors = {}
    true_vectors = {}
    coverage_rows = []

    for _, row in model_df.iterrows():
        protein = str(row.get("protein"))
        coeffs = parse_coefficients(row.get("feat_coeffs_dict", {}))
        needed = list(map(str, coeffs.keys()))
        if not needed:
            continue
        missing = [g for g in needed if g not in available]
        present = [g for g in needed if g in available]
        coverage_rows.append({
            "Protein": protein,
            "n_features": len(needed),
            "n_present": len(present),
            "n_missing": len(missing),
            "missing_fraction": len(missing) / max(1, len(needed)),
            "Target_Present": protein in prot_df.columns,
        })
        if missing and not impute_missing_features:
            continue
        if impute_missing_features:
            weights = np.asarray([float(coeffs[g]) for g in present], dtype=float)
            y_pred = X_all[present].dot(weights) if present else pd.Series(0.0, index=X_all.index)
        else:
            weights = np.asarray([float(coeffs[g]) for g in needed], dtype=float)
            y_pred = X_all[needed].dot(weights)

        target_present = protein in prot_df.columns
        y_true = prot_df[protein] if target_present else pd.Series(index=X_all.index, dtype=float)
        common = y_pred.index.intersection(y_true.index)
        mse_val = np.nan
        corr_val = np.nan
        cos_val = np.nan
        n_used = 0
        if target_present and len(common) > 0:
            yp = y_pred.loc[common].astype(float)
            yt = y_true.loc[common].astype(float)
            mask = np.isfinite(yp.to_numpy()) & np.isfinite(yt.to_numpy())
            n_used = int(mask.sum())
            if n_used > 0:
                yp_v = yp.to_numpy()[mask]
                yt_v = yt.to_numpy()[mask]
                mse_val = float(np.mean((yp_v - yt_v) ** 2))
                if n_used >= 2 and np.std(yp_v) > 0 and np.std(yt_v) > 0:
                    corr_val = float(np.corrcoef(yp_v, yt_v)[0, 1])
                denom = float(np.linalg.norm(yp_v) * np.linalg.norm(yt_v))
                cos_val = float(np.dot(yp_v, yt_v) / denom) if denom > 0 else np.nan
            error_vectors[protein] = yp - yt
            pred_vectors[protein] = yp
            true_vectors[protein] = yt
        result_rows.append({
            "Protein": protein,
            "MSE": mse_val,
            "Correlation": corr_val,
            "Cosine_Similarity": cos_val,
            "Target_Present": target_present,
            "n_used": n_used,
            "n_missing_features": len(missing),
            "missing_feature_fraction": len(missing) / max(1, len(needed)),
        })

    results_df = pd.DataFrame(result_rows)
    error_table = pd.DataFrame(error_vectors).T if error_vectors else pd.DataFrame()
    predictions_df = pd.DataFrame(pred_vectors) if pred_vectors else pd.DataFrame(index=X_all.index)
    truth_df = pd.DataFrame(true_vectors) if true_vectors else pd.DataFrame(index=X_all.index)
    coverage_df = pd.DataFrame(coverage_rows)
    return results_df, error_table, predictions_df, truth_df, coverage_df


def tissue_mse_from_error_table(error_table, tissue_func=None):
    if error_table is None or error_table.empty:
        return pd.DataFrame(columns=["Sample", "MSE", "Tissue"])
    tissue_func = tissue_func or get_tissue
    mse_per_sample = (error_table.astype(float) ** 2).mean(axis=0)
    return pd.DataFrame({
        "Sample": mse_per_sample.index.astype(str),
        "MSE": mse_per_sample.values,
        "Tissue": [tissue_func(x) for x in mse_per_sample.index.astype(str)],
    })


def summarize_cptac_results(results_by_scope):
    rows = []
    for key, bundle in results_by_scope.items():
        dataset, scope = key
        df = bundle.get("results_df", pd.DataFrame())
        coverage = bundle.get("coverage_df", pd.DataFrame())
        rna_df = bundle.get("rna_df")
        prot_df = bundle.get("prot_df")
        rows.append({
            "Dataset": dataset,
            "Scope": scope,
            "Samples": 0 if rna_df is None else len(rna_df),
            "Proteins_evaluated": int(df["Protein"].nunique()) if len(df) and "Protein" in df else 0,
            "Targets_present": int(df["Target_Present"].sum()) if len(df) and "Target_Present" in df else 0,
            "Mean_MSE": df["MSE"].mean() if len(df) and "MSE" in df else np.nan,
            "Median_MSE": df["MSE"].median() if len(df) and "MSE" in df else np.nan,
            "Mean_Correlation": df["Correlation"].mean() if len(df) and "Correlation" in df else np.nan,
            "Median_Correlation": df["Correlation"].median() if len(df) and "Correlation" in df else np.nan,
            "Mean_Cosine": df["Cosine_Similarity"].mean() if len(df) and "Cosine_Similarity" in df else np.nan,
            "Median_missing_feature_fraction": coverage["missing_fraction"].median() if len(coverage) and "missing_fraction" in coverage else np.nan,
        })
    return pd.DataFrame(rows)


def plot_cptac_summary(summary_df):
    if summary_df is None or summary_df.empty:
        print("No summary data to plot.")
        return None
    plot_df = summary_df.melt(
        id_vars=["Dataset", "Scope"],
        value_vars=["Mean_Correlation", "Mean_Cosine", "Mean_MSE"],
        var_name="Metric",
        value_name="Value",
    )
    g = sns.catplot(
        data=plot_df,
        kind="bar",
        x="Dataset",
        y="Value",
        hue="Scope",
        col="Metric",
        sharey=False,
        height=4,
        aspect=1.05,
        errorbar=None,
    )
    g.set_titles("{col_name}")
    g.set_axis_labels("", "Value")
    for ax in g.axes.flat:
        ax.tick_params(axis="x", rotation=20)
        ax.grid(axis="y", alpha=0.25)
    g.fig.suptitle("CPTAC validation performance summary", y=1.05, fontsize=15, fontweight="bold")
    return g.fig


def plot_cptac_metric_distributions(long_results_df):
    if long_results_df is None or long_results_df.empty:
        print("No protein-level result rows to plot.")
        return None
    metrics = ["Correlation", "Cosine_Similarity", "MSE"]
    for metric in metrics:
        data = long_results_df.dropna(subset=[metric]).copy()
        if data.empty:
            continue
        plt.figure(figsize=(11, 4.5))
        sns.boxplot(data=data, x="Dataset", y=metric, hue="Scope", showfliers=False)
        sns.stripplot(
            data=data.sample(min(len(data), 1500), random_state=42),
            x="Dataset", y=metric, hue="Scope", dodge=True, size=1.5, alpha=0.18, legend=False,
        )
        plt.title(f"Protein-level {metric} distribution", weight="bold")
        plt.xlabel("")
        plt.grid(axis="y", alpha=0.2)
        if metric == "MSE":
            plt.yscale("log")
        plt.tight_layout()
        plt.show()


def plot_cptac_tissue_mse(tissue_df):
    if tissue_df is None or tissue_df.empty:
        print("No tissue MSE rows to plot.")
        return None
    plot_df = (
        tissue_df.groupby(["Dataset", "Scope", "Tissue"], as_index=False)["MSE"]
        .mean()
        .sort_values("MSE", ascending=False)
    )
    for dataset in plot_df["Dataset"].dropna().unique():
        sub = plot_df[plot_df["Dataset"].eq(dataset)].copy()
        keep = sub.groupby("Tissue")["MSE"].max().sort_values(ascending=False).head(20).index
        sub = sub[sub["Tissue"].isin(keep)]
        plt.figure(figsize=(10, max(4, 0.35 * len(keep))))
        sns.barplot(data=sub, y="Tissue", x="MSE", hue="Scope", errorbar=None)
        plt.title(f"{dataset}: tissue-level mean MSE", weight="bold")
        plt.xlabel("Mean MSE across proteins")
        plt.ylabel("")
        plt.grid(axis="x", alpha=0.2)
        plt.tight_layout()
        plt.show()


def evaluate_cptac_validation_datasets(dataset_configs, best_models, excluded_tissues=None, impute_missing_features=False):
    """Evaluate coefficient-dictionary models across CPTAC validation datasets."""
    excluded_tissues = excluded_tissues or []
    results = {}
    for dataset_name, cfg in dataset_configs.items():
        print("\n" + "=" * 90)
        print(f"Dataset: {dataset_name}")
        print("=" * 90)
        rna_df, prot_df = load_validation_dataset(cfg["rna_path"], cfg["prot_path"])
        for scope, filter_tissues in [("all_tissues", []), ("filtered_tissues", excluded_tissues)]:
            rna_use = rna_df.fillna(0).copy()
            prot_use = prot_df.copy()
            if filter_tissues:
                drop_samples = [s for s in rna_use.index if get_tissue(s) in set(filter_tissues)]
                rna_use = rna_use.drop(index=drop_samples, errors="ignore")
                prot_use = prot_use.drop(index=drop_samples, errors="ignore")
            pred_results, error_table, pred_df, true_df, coverage_df = predict_from_best_models(
                best_models,
                rna_use,
                prot_use,
                impute_missing_features=impute_missing_features,
            )
            tissue_df = tissue_mse_from_error_table(error_table)
            tissue_df["Dataset"] = dataset_name
            tissue_df["Scope"] = scope
            if len(pred_results):
                print(f"{scope}: proteins={len(pred_results):,}, samples={len(rna_use):,}, mean corr={pred_results['Correlation'].mean():.4f}")
            results[(dataset_name, scope)] = {
                "rna_df": rna_use,
                "prot_df": prot_use,
                "results_df": pred_results,
                "error_table": error_table,
                "predictions_df": pred_df,
                "truth_df": true_df,
                "coverage_df": coverage_df,
                "tissue_df": tissue_df,
                "zero_stats": percent_zeros_by_tissue(rna_use),
            }
    return results


def evaluate_norm_config(norm_label, results_dir, val_data_dir, datasets=None, impute_missing_features=True):
    """Evaluate one normalization result directory on all validation datasets."""
    datasets = datasets or ["iTRAQ", "LabelFree", "TMT"]
    best_models, files = load_best_result_tables(results_dir)
    if best_models.empty:
        print(f"No best-model files found for {norm_label} in {results_dir}")
        return {"norm": norm_label, "best_models": best_models, "results": {}, "summary": pd.DataFrame(), "tissue": pd.DataFrame(), "protein": pd.DataFrame()}
    dataset_configs = {
        ds: {
            "rna_path": Path(val_data_dir) / f"{ds}_transcriptomics.joblib",
            "prot_path": Path(val_data_dir) / f"{ds}_proteomics.joblib",
        }
        for ds in datasets
    }
    raw_results = {}
    protein_rows = []
    tissue_rows = []
    summary_rows = []
    for ds, cfg in dataset_configs.items():
        try:
            rna_df, prot_df = load_validation_dataset(cfg["rna_path"], cfg["prot_path"])
        except Exception as exc:
            print(f"Skipping {ds} for {norm_label}: {type(exc).__name__}: {exc}")
            continue
        results_df, error_table, pred_df, truth_df, coverage_df = predict_from_best_models(
            best_models,
            rna_df.fillna(0),
            prot_df,
            impute_missing_features=impute_missing_features,
        )
        prot_metrics = compute_detailed_tables(truth_df, pred_df)[0] if not pred_df.empty and not truth_df.empty else pd.DataFrame()
        if len(prot_metrics):
            prot_metrics["Dataset"] = ds
            prot_metrics["Normalization"] = norm_label
            protein_rows.append(prot_metrics)
        tissue_df = tissue_mse_from_error_table(error_table)
        if len(tissue_df):
            tissue_df["Dataset"] = ds
            tissue_df["Normalization"] = norm_label
            tissue_rows.append(tissue_df)
        sample_metrics = {}
        if len(error_table):
            sample_mse = (error_table.astype(float) ** 2).mean(axis=0)
            sample_metrics["Mean_MSE_sample"] = sample_mse.mean()
        summary_rows.append({
            "Normalization": norm_label,
            "Dataset": ds,
            "Samples": len(rna_df),
            "Proteins_evaluated": int(results_df["Protein"].nunique()) if len(results_df) else 0,
            "Mean_PCC_protein": results_df["Correlation"].mean() if len(results_df) else np.nan,
            "Median_PCC_protein": results_df["Correlation"].median() if len(results_df) else np.nan,
            "Mean_Cosine_protein": results_df["Cosine_Similarity"].mean() if len(results_df) else np.nan,
            "Mean_MSE_protein": results_df["MSE"].mean() if len(results_df) else np.nan,
            "Median_MSE_protein": results_df["MSE"].median() if len(results_df) else np.nan,
            **sample_metrics,
            "Median_missing_feature_fraction": coverage_df["missing_fraction"].median() if len(coverage_df) else np.nan,
        })
        raw_results[ds] = {
            "results_df": results_df,
            "error_table": error_table,
            "predictions_df": pred_df,
            "truth_df": truth_df,
            "coverage_df": coverage_df,
            "protein_metrics": prot_metrics,
            "tissue_df": tissue_df,
        }
    return {
        "norm": norm_label,
        "best_models": best_models,
        "best_files": files,
        "results": raw_results,
        "summary": pd.DataFrame(summary_rows),
        "tissue": pd.concat(tissue_rows, ignore_index=True) if tissue_rows else pd.DataFrame(),
        "protein": pd.concat(protein_rows, ignore_index=True) if protein_rows else pd.DataFrame(),
    }


def plot_norm_summary(summary_df):
    if summary_df is None or summary_df.empty:
        print("No normalization summary rows to plot.")
        return None
    plot_df = summary_df.melt(
        id_vars=["Normalization", "Dataset"],
        value_vars=["Mean_PCC_protein", "Mean_Cosine_protein", "Mean_MSE_protein"],
        var_name="Metric",
        value_name="Value",
    )
    g = sns.catplot(
        data=plot_df,
        kind="bar",
        x="Dataset",
        y="Value",
        hue="Normalization",
        col="Metric",
        sharey=False,
        height=4,
        aspect=1.05,
        errorbar=None,
    )
    for ax in g.axes.flat:
        ax.grid(axis="y", alpha=0.25)
        ax.tick_params(axis="x", rotation=20)
    g.set_axis_labels("", "Value")
    g.set_titles("{col_name}")
    g.fig.suptitle("Normalization comparison across CPTAC datasets", y=1.06, fontsize=15, fontweight="bold")
    return g.fig


def plot_norm_protein_distributions(protein_df):
    if protein_df is None or protein_df.empty:
        print("No protein metric rows to plot.")
        return None
    metric_map = {
        "pearson": "Pearson correlation",
        "cosine": "Cosine similarity",
        "mse": "MSE",
    }
    for metric, label in metric_map.items():
        if metric not in protein_df.columns:
            continue
        data = protein_df.dropna(subset=[metric]).copy()
        if data.empty:
            continue
        plt.figure(figsize=(12, 4.8))
        sns.boxplot(data=data, x="Dataset", y=metric, hue="Normalization", showfliers=False)
        plt.title(f"Protein-wise {label} by normalization", weight="bold")
        plt.xlabel("")
        plt.ylabel(label)
        if metric == "mse":
            plt.yscale("log")
        plt.grid(axis="y", alpha=0.2)
        plt.tight_layout()
        plt.show()


def build_top_protein_table(protein_df, top_n=20):
    if protein_df is None or protein_df.empty:
        return pd.DataFrame()
    rows = []
    specs = [("pearson", False), ("cosine", False), ("mse", True)]
    for (norm, ds), sub in protein_df.groupby(["Normalization", "Dataset"]):
        for metric, ascending in specs:
            if metric not in sub.columns:
                continue
            tmp = sub.dropna(subset=[metric]).sort_values(metric, ascending=ascending).head(top_n)
            for rank, (_, row) in enumerate(tmp.iterrows(), start=1):
                rows.append({
                    "Normalization": norm,
                    "Dataset": ds,
                    "Metric": metric,
                    "Rank": rank,
                    "Protein": row.get("protein"),
                    "Value": row.get(metric),
                })
    return pd.DataFrame(rows)


def high_performance_counts(protein_df, threshold=0.5):
    if protein_df is None or protein_df.empty:
        return pd.DataFrame()
    rows = []
    for (norm, ds), sub in protein_df.groupby(["Normalization", "Dataset"]):
        pcc_valid = sub["pearson"].dropna() if "pearson" in sub else pd.Series(dtype=float)
        cos_valid = sub["cosine"].dropna() if "cosine" in sub else pd.Series(dtype=float)
        rows.append({
            "Normalization": norm,
            "Dataset": ds,
            "Valid_Pearson": len(pcc_valid),
            f"Pearson_gt_{threshold}": int((pcc_valid > threshold).sum()),
            "Pearson_share": float((pcc_valid > threshold).mean()) if len(pcc_valid) else np.nan,
            "Valid_Cosine": len(cos_valid),
            f"Cosine_gt_{threshold}": int((cos_valid > threshold).sum()),
            "Cosine_share": float((cos_valid > threshold).mean()) if len(cos_valid) else np.nan,
        })
    return pd.DataFrame(rows)


def plot_norm_tissue_heatmap(tissue_df, metric="MSE"):
    if tissue_df is None or tissue_df.empty:
        print("No tissue rows to plot.")
        return None
    summary = tissue_df.groupby(["Dataset", "Normalization", "Tissue"], as_index=False)[metric].mean()
    for ds in summary["Dataset"].dropna().unique():
        sub = summary[summary["Dataset"].eq(ds)].copy()
        keep = sub.groupby("Tissue")[metric].max().sort_values(ascending=False).head(25).index
        pivot = sub[sub["Tissue"].isin(keep)].pivot_table(index="Tissue", columns="Normalization", values=metric, aggfunc="mean")
        pivot = pivot.loc[pivot.mean(axis=1).sort_values(ascending=False).index]
        plt.figure(figsize=(8, max(4, 0.35 * len(pivot))))
        vals = pivot.to_numpy().ravel()
        vals = vals[np.isfinite(vals)]
        vmax = np.nanpercentile(vals, 90) if len(vals) else None
        sns.heatmap(pivot, annot=True, fmt=".3g", linewidths=0.4, cbar_kws={"label": metric}, vmax=vmax)
        plt.title(f"{ds}: tissue-level {metric} by normalization", weight="bold")
        plt.xlabel("Normalization")
        plt.ylabel("Tissue")
        plt.tight_layout()
        plt.show()


# -----------------------------------------------------------------------------
# Extracted notebook helper functions
# -----------------------------------------------------------------------------

def _as_finite_array(values, default=0.0):
    arr = np.asarray(values, dtype=np.float32)
    arr[~np.isfinite(arr)] = default
    return arr

def _display(t):
    return (t.replace('_', ' ').replace('HEADNECK', 'Head & Neck').replace('PANCREA', 'Pancreas').replace('PEDBRAIN', 'Ped. Brain').title())

def _fmt_float(x, digits=4):
    if pd.isna(x):
        return "NA"
    try:
        return f"{float(x):.{digits}g}"
    except Exception:
        return str(x)

def _get_global(name):
    return globals().get(name, None)

def _index_by_protein(df):
    if df is None or len(df) == 0:
        return pd.DataFrame().rename_axis("protein")
    out = df.copy()
    out["protein"] = out["protein"].astype(str)
    return out.drop_duplicates("protein", keep="first").set_index("protein")

def _kde_smooth(
    ax,
    vals_my,
    vals_ref,
    lab_my,
    lab_ref,
    palette,
    n=20000,
    seed=0,
    bw_adjust=1.35,      # ↑ makes peaks smoother
    gridsize=1024,       # ↑ more points -> smoother curve
    cut=0,               # don't extend beyond data range
    jitter_if_discrete=True
):
    rng = np.random.default_rng(seed)

    v1 = np.asarray(vals_my)
    v2 = np.asarray(vals_ref)
    v1 = v1[~np.isnan(v1)]
    v2 = v2[~np.isnan(v2)]

    if v1.size > n:
        v1 = rng.choice(v1, size=n, replace=False)
    if v2.size > n:
        v2 = rng.choice(v2, size=n, replace=False)

    # If there are MANY identical values (common after imputation/clipping),
    # KDE can look pointy/jagged. Add tiny jitter for PLOTTING ONLY.
    if jitter_if_discrete:
        for arr_name, arr in [("val", v1), ("ref", v2)]:
            if arr.size > 0:
                uniq_ratio = np.unique(arr).size / arr.size
                if uniq_ratio < 0.02:  # very discrete
                    eps = max(1e-6, 1e-3 * np.std(arr))
                    if arr_name == "val":
                        v1 = arr + rng.normal(0, eps, size=arr.size)
                    else:
                        v2 = arr + rng.normal(0, eps, size=arr.size)

    sns.kdeplot(
        v1, ax=ax, fill=True,
        color=palette["Val"], alpha=0.4, label=lab_my,
        bw_adjust=bw_adjust, gridsize=gridsize, cut=cut, warn_singular=False
    )
    sns.kdeplot(
        v2, ax=ax, fill=True,
        color=palette["Ref"], alpha=0.2, label=lab_ref,
        linestyle="--",
        bw_adjust=bw_adjust, gridsize=gridsize, cut=cut, warn_singular=False
    )

def _safe_name(x, max_len=120):
    x = str(x)
    x = x.replace("≥", "ge").replace("≤", "le")
    x = x.replace(">", "gt").replace("<", "lt").replace("=", "eq")
    x = x.replace(".", "p").replace("/", "_")
    x = re.sub(r"[^A-Za-z0-9_+-]+", "_", x)
    x = re.sub(r"_+", "_", x).strip("_")

    if len(x) > max_len:
        h = hashlib.md5(x.encode("utf-8")).hexdigest()[:8]
        x = f"{x[:max_len]}_{h}"

    return x

def _sample_center_df(X, how=None):
    if how is None:
        return X
    if how == "median":
        return X.sub(X.median(axis=1), axis=0)
    if how == "mean":
        return X.sub(X.mean(axis=1), axis=0)
    raise ValueError(how)

def _standardize_assignment_table(df, analysis_name, source_variable, default_dataset):
    """
    Convert notebook-specific assignment tables into one common schema.

    Crucial fix:
    We initialize the output DataFrame with the same index as df, so scalar
    metadata columns are repeated correctly and do not become NaN.
    """

    df = df.copy()

    protein_col = None
    for c in ["protein", "protein_name", "target_protein", "target", "gene", "gene_symbol"]:
        if c in df.columns:
            protein_col = c
            break

    stage_col = None
    for c in ["stage", "model", "family", "assigned_stage", "winner_stage", "winner"]:
        if c in df.columns:
            stage_col = c
            break

    if protein_col is None:
        print(f"Skipped {source_variable}: no protein column found.")
        return None

    if stage_col is None:
        print(f"Skipped {source_variable}: no stage/model column found.")
        return None

    out = pd.DataFrame(index=df.index)

    out["analysis"] = analysis_name
    out["source_variable"] = source_variable

    if "dataset_display" in df.columns:
        out["dataset"] = df["dataset_display"].astype(str).values
    elif "dataset" in df.columns:
        out["dataset"] = df["dataset"].astype(str).values
    elif "dataset_key" in df.columns:
        out["dataset"] = df["dataset_key"].astype(str).values
    else:
        out["dataset"] = default_dataset

    out["threshold_label"] = df.apply(
        lambda row: _threshold_label(row, analysis_name),
        axis=1,
    ).values

    out["stage"] = df[stage_col].astype(str).values
    out["protein"] = df[protein_col].astype(str).values

    # Keep useful metric columns if present.
    useful_extra_cols = [
        "threshold",
        "prediction_threshold",
        "corr_threshold",
        "mse_percentile_threshold",
        "mse_percentile",
        "percentile",
        "metric",
        "test_pearson_corr",
        "test_mse",
        "test_cosine_similarity",
        "best_mse",
        "proxy_mse",
        "enet_mse",
        "xgb_mse",
        "same_gene_corr_full",
        "proxy_corr",
        "proxy_corr_column",
        "run_label",
        "model_name",
        "family",
        "proxy_type",
    ]

    for c in useful_extra_cols:
        if c in df.columns:
            out[c] = df[c].values

    out = out.dropna(subset=["protein", "stage"])
    out = out[out["protein"].astype(str).str.len() > 0]
    out = out.drop_duplicates().reset_index(drop=True)

    return out

def _threshold_label(row, analysis_name):
    """
    Creates meaningful labels instead of blindly using the plotting threshold.
    Important because MSE winner/improvement plots use threshold = 0 only as
    an internal plotting placeholder.
    """

    a = str(analysis_name).lower()

    # No-threshold MSE partitions.
    if "mse_winner" in a:
        return "All proteins / lowest MSE winner / no fixed threshold"

    if "mse_procedural_improvement" in a or "mse_improvement" in a:
        return "All proteins / procedural MSE improvement / no fixed threshold"

    # Proxy-correlation eligibility sweep.
    if "proxy_corr" in a or "proxy_correlation" in a:
        parts = []

        if "prediction_threshold" in row.index and pd.notna(row["prediction_threshold"]):
            parts.append(f"validation Pearson >= {float(row['prediction_threshold']):.2f}")

        if "corr_threshold" in row.index and pd.notna(row["corr_threshold"]):
            parts.append(f"proxy corr >= {float(row['corr_threshold']):.2f}")
        elif "threshold" in row.index and pd.notna(row["threshold"]):
            parts.append(f"proxy corr >= {float(row['threshold']):.2f}")

        return " | ".join(parts) if parts else "Proxy-correlation threshold"

    # Percentile MSE threshold sweep.
    if "mse_percentile" in a:
        pct = None
        for c in ["mse_percentile_threshold", "mse_percentile", "percentile"]:
            if c in row.index and pd.notna(row[c]):
                pct = row[c]
                break

        if "threshold" in row.index and pd.notna(row["threshold"]):
            mse_thr = row["threshold"]
            if pct is not None:
                return f"MSE <= p{int(float(pct))} ({_fmt_float(mse_thr)})"
            return f"MSE <= {_fmt_float(mse_thr)}"

        if pct is not None:
            return f"MSE <= p{int(float(pct))}"

        return "MSE percentile threshold"

    # Absolute MSE threshold sweep.
    if "mse_threshold" in a:
        if "threshold" in row.index and pd.notna(row["threshold"]):
            return f"MSE <= {_fmt_float(row['threshold'])}"
        return "MSE threshold"

    # Pearson cascades and Pearson sweeps.
    if (
        "pearson" in a
        or "sctranslator" in a
        or "sc_rescue" in a
        or "final_sc" in a
    ):
        if "threshold" in row.index and pd.notna(row["threshold"]):
            return f"Pearson >= {float(row['threshold']):.2f}"
        return "Pearson threshold"

    # Cosine sweeps.
    if "cosine" in a:
        if "threshold" in row.index and pd.notna(row["threshold"]):
            return f"Cosine >= {float(row['threshold']):.2f}"
        return "Cosine threshold"

    # Fallback.
    if "threshold" in row.index and pd.notna(row["threshold"]):
        return f"threshold = {_fmt_float(row['threshold'])}"

    return "No fixed threshold"

def _valid_df(x):
    return isinstance(x, pd.DataFrame) and len(x) > 0

def align_sctranslator_predictions(prediction_array, metadata, enet_bundle):
    arr = np.asarray(prediction_array, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"Expected a 2D scTranslator prediction matrix, got shape {arr.shape}")

    n_rows = len(enet_bundle["row_order"])
    n_proteins = len(enet_bundle["protein_order"])

    protein_labels = extract_labels_from_metadata(
        metadata,
        [
            "proteins", "protein_names", "target_proteins", "protein_ids", "columns",
            "colnames", "var_names", "genes", "gene_names",
        ],
    )
    sample_labels = extract_labels_from_metadata(
        metadata,
        [
            "samples", "sample_ids", "cell_lines", "cell_line_ids", "obs_names",
            "index", "rownames",
        ],
    )

    if protein_labels is not None and len(protein_labels) in arr.shape:
        if arr.shape[0] == len(protein_labels) and arr.shape[1] != len(protein_labels):
            arr = arr.T

        available_proteins = pd.Index([str(x) for x in protein_labels])
        wanted_proteins = pd.Index([str(x) for x in enet_bundle["protein_order"]])
        common_proteins = wanted_proteins.intersection(available_proteins)

        if len(common_proteins) == 0:
            raise ValueError(
                "Protein labels were found in scTranslator metadata, but none matched ENet proteins."
            )

        pred_col_positions = available_proteins.get_indexer(common_proteins)
        enet_col_positions = wanted_proteins.get_indexer(common_proteins)

        arr = arr[:, pred_col_positions]
        y_true = enet_bundle["y_true"][:, enet_col_positions]
        enet_pred = enet_bundle["enet_pred"][:, enet_col_positions]
        proteins = common_proteins.tolist()

    else:
        if arr.shape == (n_rows, n_proteins):
            y_true = enet_bundle["y_true"]
            enet_pred = enet_bundle["enet_pred"]
            proteins = enet_bundle["protein_order"]
        elif arr.shape == (n_proteins, n_rows):
            arr = arr.T
            y_true = enet_bundle["y_true"]
            enet_pred = enet_bundle["enet_pred"]
            proteins = enet_bundle["protein_order"]
        elif arr.shape[0] >= n_rows and arr.shape[1] >= n_proteins:
            print(
                "Warning: scTranslator matrix is larger than the ENet validation matrix. "
                "Using first rows and first proteins in ENet order."
            )
            arr = arr[:n_rows, :n_proteins]
            y_true = enet_bundle["y_true"]
            enet_pred = enet_bundle["enet_pred"]
            proteins = enet_bundle["protein_order"]
        elif arr.shape[1] >= n_rows and arr.shape[0] >= n_proteins:
            print(
                "Warning: scTranslator matrix appears transposed and larger than ENet matrix. "
                "Transposing, then using first rows/proteins in ENet order."
            )
            arr = arr.T[:n_rows, :n_proteins]
            y_true = enet_bundle["y_true"]
            enet_pred = enet_bundle["enet_pred"]
            proteins = enet_bundle["protein_order"]
        else:
            raise ValueError(
                "Could not align scTranslator predictions to ENet validation data.\n"
                f"scTranslator shape: {arr.shape}\n"
                f"Expected either {(n_rows, n_proteins)} or {(n_proteins, n_rows)}."
            )

    if sample_labels is not None and len(sample_labels) == arr.shape[0]:
        available_rows = pd.Index([str(x) for x in sample_labels])
        wanted_rows = pd.Index([str(x) for x in enet_bundle["row_order"]])
        common_rows = wanted_rows.intersection(available_rows)

        if len(common_rows) > 0:
            pred_row_positions = available_rows.get_indexer(common_rows)
            enet_row_positions = wanted_rows.get_indexer(common_rows)
            arr = arr[pred_row_positions, :]
            y_true = y_true[enet_row_positions, :]
            enet_pred = enet_pred[enet_row_positions, :]

    return y_true, enet_pred, arr, proteins

def apply_artifact_normalizer(X_new: pd.DataFrame, artifact: dict) -> tuple:
    """Return transformed design matrix in the exact selected-feature order plus coverage stats."""
    selected = list(map(str, artifact.get("selected_features", [])))
    if not selected:
        raise ValueError("Artifact has no selected_features.")

    normalization = str(artifact.get("normalization", "none")).lower()
    normalizer = artifact.get("normalizer", None) or {}
    method = str(normalizer.get("method", normalization)).lower() if isinstance(normalizer, dict) else normalization

    missing_features = [g for g in selected if g not in X_new.columns]
    missing_fraction = len(missing_features) / max(1, len(selected))

    # XGB qbin/rank/zscore style normalizers.
    if isinstance(normalizer, dict) and method == "zscore" and "mean" in normalizer and "std" in normalizer:
        cols = list(map(str, normalizer.get("columns", selected)))
        Xraw = X_new.reindex(columns=cols)
        mean = np.asarray(normalizer["mean"], dtype=np.float32)
        std = np.asarray(normalizer["std"], dtype=np.float32)
        std[~np.isfinite(std) | (std == 0)] = 1.0
        arr = Xraw.to_numpy(dtype=np.float32, copy=True)
        # fill missing values with train mean so transformed value is 0
        for j in range(arr.shape[1]):
            bad = ~np.isfinite(arr[:, j])
            if bad.any():
                arr[bad, j] = mean[j] if j < len(mean) else 0.0
        out = ((arr - mean) / std).astype(np.float32)
        return out, {"selected_feature_count": len(cols), "missing_feature_count": len(missing_features), "missing_feature_fraction": missing_fraction}

    if isinstance(normalizer, dict) and method in {"rank", "ranking"} and "sorted_values" in normalizer:
        cols = list(map(str, normalizer.get("columns", selected)))
        Xraw = X_new.reindex(columns=cols)
        out = np.zeros((Xraw.shape[0], len(cols)), dtype=np.float32)
        for j, col in enumerate(cols):
            sorted_col = np.asarray(normalizer.get("sorted_values", {}).get(col, []), dtype=np.float32)
            vals = Xraw.iloc[:, j].to_numpy(dtype=np.float32, copy=False)
            mask = np.isfinite(vals)
            if sorted_col.size > 0:
                out[mask, j] = np.searchsorted(sorted_col, vals[mask], side="right") / float(sorted_col.size)
        return out, {"selected_feature_count": len(cols), "missing_feature_count": len(missing_features), "missing_feature_fraction": missing_fraction}

    if isinstance(normalizer, dict) and (method.startswith("qbin") or "edges" in normalizer):
        cols = list(map(str, normalizer.get("columns", selected)))
        Xraw = X_new.reindex(columns=cols)
        n_bins = int(normalizer.get("n_bins", re.sub("[^0-9]", "", method) or 50))
        denom = max(1, n_bins - 1)
        out = np.zeros((Xraw.shape[0], len(cols)), dtype=np.float32)
        for j, col in enumerate(cols):
            edges = np.asarray(normalizer.get("edges", {}).get(col, []), dtype=np.float32)
            vals = Xraw.iloc[:, j].to_numpy(dtype=np.float32, copy=False)
            mask = np.isfinite(vals)
            if edges.size > 0:
                out[mask, j] = np.searchsorted(edges, vals[mask], side="right") / float(denom)
        return out, {"selected_feature_count": len(cols), "missing_feature_count": len(missing_features), "missing_feature_fraction": missing_fraction}

    # ENet row-wise relative representations.
    if normalization in {"binning", "relative_binning"} or method in {"binning", "relative_binning"}:
        n_bins = int(artifact.get("n_bins", normalizer.get("n_bins", 50) if isinstance(normalizer, dict) else 50))
        Xtrans = apply_row_binning(X_new, n_bins=n_bins).reindex(columns=selected, fill_value=0.0)
        return Xtrans.to_numpy(dtype=np.float32, copy=False), {"selected_feature_count": len(selected), "missing_feature_count": len(missing_features), "missing_feature_fraction": missing_fraction}

    if normalization in {"ranking", "rank"} or method in {"ranking", "rank"}:
        Xtrans = apply_row_ranking(X_new).reindex(columns=selected, fill_value=0.0)
        return Xtrans.to_numpy(dtype=np.float32, copy=False), {"selected_feature_count": len(selected), "missing_feature_count": len(missing_features), "missing_feature_fraction": missing_fraction}

    # No normalizer or simple none. Missing genes become zeros.
    Xmat = X_new.reindex(columns=selected, fill_value=0.0)
    arr = Xmat.to_numpy(dtype=np.float32, copy=True)
    arr[~np.isfinite(arr)] = 0.0
    return arr, {"selected_feature_count": len(selected), "missing_feature_count": len(missing_features), "missing_feature_fraction": missing_fraction}

def apply_row_binning(X_full: pd.DataFrame, n_bins: int = 50) -> pd.DataFrame:
    """Row-wise relative representation binning. Used for ENet binning when no fitted edges are stored."""
    x = X_full.to_numpy(dtype=np.float64, copy=True)
    x[x == 0] = np.nan
    qs = np.linspace(0.0, 1.0, int(n_bins), dtype=np.float64)
    q = np.nanquantile(x, qs, axis=1).T
    out = np.zeros_like(x, dtype=np.float32)
    for i in range(x.shape[0]):
        row = x[i]
        m = np.isfinite(row)
        if not m.any() or np.isnan(q[i]).all():
            continue
        idx = np.searchsorted(q[i], row[m], side="left")
        out[i, m] = float(n_bins) - idx.astype(np.float32)
    return pd.DataFrame(out, index=X_full.index, columns=X_full.columns)

def apply_row_ranking(X_full: pd.DataFrame) -> pd.DataFrame:
    x = X_full.to_numpy(dtype=np.float64, copy=False)
    ranked = x.argsort(axis=1).argsort(axis=1).astype(np.float32, copy=False)
    return pd.DataFrame(ranked, index=X_full.index, columns=X_full.columns)

def best_by_metric(df, metric, higher=True):
    data = df.dropna(subset=[metric]).copy()
    if data.empty:
        return data

    return (
        data.sort_values(
            ["protein", metric, "run_label"],
            ascending=[True, not higher, True],
        )
        .groupby("protein", as_index=False)
        .head(1)
        .reset_index(drop=True)
    )

def best_per_protein(df, metric, higher=True):
    if df is None or len(df) == 0 or metric not in df.columns:
        return pd.DataFrame()

    data = df.dropna(subset=[metric]).copy()
    if data.empty:
        return data

    ascending_metric = not higher
    sort_cols = ["protein", metric, "run_label"]
    return (
        data.sort_values(sort_cols, ascending=[True, ascending_metric, True])
        .groupby("protein", as_index=False)
        .head(1)
        .reset_index(drop=True)
    )

def build_enet_matrices_for_dataset(eval_df, dataset_name):
    dataset_norm = normalize_name(dataset_name)

    df = eval_df.copy()
    data_col = first_existing_column(df, ["dataset"])
    run_col = first_existing_column(df, ["run_label", "model", "model_label"])
    protein_col = first_existing_column(df, ["protein", "protein_id", "target", "gene"])
    true_col = first_existing_column(
        df,
        [
            "y_true", "true", "truth", "target_value", "observed", "observed_value",
            "value_true", "test_y", "y",
        ],
    )
    pred_col = first_existing_column(
        df,
        [
            "y_pred", "pred", "prediction", "predicted", "predicted_value",
            "value_pred", "test_pred", "y_hat",
        ],
    )
    sample_col = first_existing_column(
        df,
        [
            "sample", "sample_id", "cell_line", "cell_line_id", "depmap_id",
            "DepMap_ID", "model_id", "sample_name",
        ],
        required=False,
    )

    df["_dataset_norm"] = df[data_col].map(normalize_name)
    ds = df[df["_dataset_norm"].eq(dataset_norm)].copy()

    if ds.empty:
        raise ValueError(f"No validation prediction rows found for dataset {dataset_name}")

    enet_rows = ds[ds[run_col].map(is_enet_zbig_zscore_k5000_label)].copy()
    if enet_rows.empty:
        display(ds[[run_col]].drop_duplicates().sort_values(run_col).reset_index(drop=True))
        raise ValueError(
            f"No ENet zbig z-score k5000 run found for dataset {dataset_name}. "
            f"Expected something like: {TARGET_ENET_RUN_LABEL}"
        )

    target_compact = normalize_name(TARGET_ENET_RUN_LABEL)
    enet_rows["_run_compact"] = enet_rows[run_col].astype(str).map(normalize_name)
    exact_runs = (
        enet_rows.loc[enet_rows["_run_compact"].eq(target_compact), run_col]
        .drop_duplicates()
        .tolist()
    )

    if exact_runs:
        selected_run = exact_runs[0]
    else:
        run_counts = enet_rows.groupby(run_col)[protein_col].nunique().sort_values(ascending=False)
        selected_run = run_counts.index[0]
        if len(run_counts) > 1:
            print(f"Multiple ENet zbig z-score k5000-like runs found for {dataset_name}.")
            print("Using the one with the most proteins:", selected_run)
            display(run_counts.rename("n_proteins").reset_index())

    enet_rows = enet_rows[enet_rows[run_col].eq(selected_run)].copy()

    if sample_col is not None:
        enet_rows["_row_key"] = enet_rows[sample_col].astype(str)
    else:
        enet_rows = enet_rows.reset_index(drop=False).rename(columns={"index": "_original_order"})
        enet_rows = enet_rows.sort_values([protein_col, "_original_order"], kind="stable").copy()
        enet_rows["_row_key"] = (
            enet_rows.groupby(protein_col, sort=False).cumcount().astype(str)
        )

    true_wide = (
        enet_rows.pivot_table(
            index="_row_key",
            columns=protein_col,
            values=true_col,
            aggfunc="first",
        )
        .sort_index()
    )

    enet_pred_wide = (
        enet_rows.pivot_table(
            index="_row_key",
            columns=protein_col,
            values=pred_col,
            aggfunc="first",
        )
        .reindex(index=true_wide.index, columns=true_wide.columns)
    )

    proteins = list(true_wide.columns)
    rows = list(true_wide.index)

    return {
        "dataset": dataset_name,
        "selected_run": selected_run,
        "row_order": rows,
        "protein_order": proteins,
        "y_true": true_wide.to_numpy(dtype=float),
        "enet_pred": enet_pred_wide.to_numpy(dtype=float),
        "true_wide": true_wide,
        "enet_pred_wide": enet_pred_wide,
    }

def clean_label_list(values):
    values = np.asarray(values)
    if values.ndim == 0:
        return None
    values = values.ravel().tolist()
    cleaned = []
    for x in values:
        if isinstance(x, bytes):
            x = x.decode("utf-8")
        cleaned.append(str(x))
    return cleaned

def compact_label(s, max_len=34):
    s = str(s)
    return s if len(s) <= max_len else s[:max_len-1] + "…"

def compute_detailed_tables(ground_truth, predictions):
    # Generates per-protein stats for all metrics
    # Now storing EVERYTHING per protein
    per_protein = {
        "protein": [], "cosine": [], "pearson": [], 
        "spearman": [], "mse": [], "mae": []
    }
    
    proteins = predictions.columns.tolist()

    for protein in proteins:
        y_true = ground_truth[protein].values.astype(float)
        y_pred = predictions[protein].values.astype(float)

        mask = ~(np.isnan(y_true) | np.isnan(y_pred))
        y_true_clean = y_true[mask]
        y_pred_clean = y_pred[mask]

        if len(y_true_clean) > 1:
            per_protein["protein"].append(protein)
            per_protein["cosine"].append(cosine_similarity(y_true_clean, y_pred_clean))
            per_protein["pearson"].append(pearson_correlation(y_true_clean, y_pred_clean))
            per_protein["spearman"].append(spearman_correlation(y_true_clean, y_pred_clean))
            per_protein["mse"].append(mse(y_true_clean, y_pred_clean))
            per_protein["mae"].append(mae(y_true_clean, y_pred_clean))

    per_protein_df = pd.DataFrame(per_protein)
    
    # Sort by Pearson for the Top/Bottom analysis
    per_protein_df_sorted = per_protein_df.sort_values("pearson", ascending=False)
    top_10 = per_protein_df_sorted.head(10)
    bottom_10 = per_protein_df_sorted.tail(10)

    return per_protein_df, top_10, bottom_10

def compute_group_metrics(df, group_cols):
    rows = []
    if df is None or df.empty:
        return pd.DataFrame()
    for key, sub in tqdm(df.groupby(group_cols), desc="Computing grouped metrics"):
        if not isinstance(key, tuple):
            key = (key,)
        m = evaluate_metric_vector(sub["y_true"].to_numpy(), sub["y_pred"].to_numpy())
        row = dict(zip(group_cols, key))
        row.update({
            "test_mse": m["mse"],
            "test_pearson_corr": m["pearson"],
            "test_cosine_similarity": m["cosine_similarity"],
            "n_prediction_pairs": m["n_used"],
            "n_samples": sub["sample"].nunique(),
            "n_proteins": sub["protein"].nunique(),
        })
        rows.append(row)
    return pd.DataFrame(rows)

def compute_mean_metrics(ground_truth, predictions):
    # Calculates mean metrics across PROTEINS (instead of samples)
    metrics_per_protein = {
        "MSE": [], "RMSE": [], "MAE": [],
        "Pearson": [], "Spearman": [], "Cosine": []
    }

    # Iterate over columns (Proteins)
    for protein in predictions.columns:
        y_true = ground_truth[protein].values.astype(float)
        y_pred = predictions[protein].values.astype(float)

        mask = ~(np.isnan(y_true) | np.isnan(y_pred))
        y_true_clean = y_true[mask]
        y_pred_clean = y_pred[mask]

        if len(y_true_clean) > 1: # Pearson/Spearman need at least 2 points
            metrics_per_protein["MSE"].append(mse(y_true_clean, y_pred_clean))
            metrics_per_protein["RMSE"].append(rmse(y_true_clean, y_pred_clean))
            metrics_per_protein["MAE"].append(mae(y_true_clean, y_pred_clean))
            metrics_per_protein["Pearson"].append(pearson_correlation(y_true_clean, y_pred_clean))
            metrics_per_protein["Spearman"].append(spearman_correlation(y_true_clean, y_pred_clean))
            metrics_per_protein["Cosine"].append(cosine_similarity(y_true_clean, y_pred_clean))

    return {k: np.nanmean(v) for k, v in metrics_per_protein.items()}

def compute_protein_wise_metrics_all(y_true, y_pred):
    n_prots = y_true.shape[1]
    results = {
        'PCC': np.zeros(n_prots),
        'MSE': np.zeros(n_prots),
        'MAE': np.zeros(n_prots),
        'R2': np.zeros(n_prots),
        'Cosine': np.zeros(n_prots)
    }
    
    for i in range(n_prots):
        yt = y_true[:, i]
        yp = y_pred[:, i]
        mask = (~np.isnan(yt)) & (~np.isnan(yp))
        
        if mask.sum() > 2:
            r, _ = pearsonr(yt[mask], yp[mask])
            results['PCC'][i] = r
            results['MSE'][i] = mean_squared_error(yt[mask], yp[mask])
            results['MAE'][i] = mean_absolute_error(yt[mask], yp[mask])
            results['R2'][i] = r2_score(yt[mask], yp[mask])
            denom = (np.linalg.norm(yp[mask]) * np.linalg.norm(yt[mask]))
            results['Cosine'][i] = np.dot(yp[mask], yt[mask]) / denom if denom > 0 else 0.0
        else:
            results['PCC'][i] = 0.0
            results['MSE'][i] = np.nan
            results['MAE'][i] = np.nan
            results['R2'][i] = np.nan
            results['Cosine'][i] = 0.0
            
    return results

def cosine_similarity(y_true, y_pred): 
     denom = (np.linalg.norm(y_true) * np.linalg.norm(y_pred))
     return np.dot(y_true, y_pred) / denom if denom > 0 else 0.0

def dataset_display_name(dataset_key):
    mapping = {
        "itraq": "iTRAQ",
        "tmt": "TMT",
        "labelfree": "LabelFree",
    }
    return mapping.get(normalize_name(dataset_key), str(dataset_key))

def evaluate_and_save_source(source_folder: str, force: bool = None) -> dict:
    """Evaluate one folder, save its results immediately, and return a result bundle."""
    if force is None:
        force = FORCE_REEVALUATE
    if source_folder not in MODEL_DIRS:
        raise KeyError(f"Unknown source_folder {source_folder!r}. Expected one of: {list(MODEL_DIRS)}")

    paths = source_output_paths(source_folder)
    can_reload = (
        (paths["predictions_parquet"].exists() or paths["predictions_csv"].exists())
        and paths["metrics_protein_csv"].exists()
        and paths["feature_coverage_csv"].exists()
    )

    if (not force) and can_reload:
        print(f"Loading existing saved results for {source_folder}")
        predictions = read_saved_predictions(paths["predictions_parquet"], paths["predictions_csv"])
        metrics = read_csv_if_exists(paths["metrics_protein_csv"], METRIC_COLUMNS)
        coverage = read_csv_if_exists(paths["feature_coverage_csv"], COVERAGE_COLUMNS)
    else:
        sub_inventory = model_inventory[model_inventory["source_folder"].eq(source_folder)].copy()
        print(f"\nEvaluating {source_folder}: {len(sub_inventory)} model files")
        predictions, metrics, coverage = evaluate_model_inventory_subset(
            sub_inventory,
            desc=f"Evaluating {source_folder}",
        )

        print(f"\nSaving per-folder results for {source_folder}...")
        write_predictions_with_fallback(predictions, paths["predictions_parquet"], paths["predictions_csv"])
        metrics.to_csv(paths["metrics_protein_csv"], index=False)
        coverage.to_csv(paths["feature_coverage_csv"], index=False)
        print("Saved:", paths["metrics_protein_csv"])
        print("Saved:", paths["feature_coverage_csv"])

    print(f"\nFinished {source_folder}")
    print("Prediction rows:", len(predictions))
    print("Protein metric rows:", len(metrics))
    if len(metrics) and "dataset" in metrics:
        print("Datasets:", sorted(metrics["dataset"].dropna().unique()))

    return {
        "source_folder": source_folder,
        "predictions": predictions,
        "metrics_by_protein": metrics,
        "feature_coverage": coverage,
        "paths": paths,
    }

def evaluate_metric_vector(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    mask = np.isfinite(y_true) & np.isfinite(y_pred)

    out = {
        "mse": np.nan,
        "pearson": np.nan,
        "cosine_similarity": np.nan,
        "n_used": int(mask.sum()),
    }
    if mask.sum() == 0:
        return out

    yt = y_true[mask]
    yp = y_pred[mask]
    out["mse"] = float(np.mean((yt - yp) ** 2))
    if len(yt) >= 2 and np.std(yt) > 0 and np.std(yp) > 0:
        out["pearson"] = float(np.corrcoef(yt, yp)[0, 1])
    denom = float(np.linalg.norm(yt) * np.linalg.norm(yp))
    if denom > 0:
        out["cosine_similarity"] = float(np.dot(yt, yp) / denom)
    return out

def evaluate_model_inventory_subset(sub_inventory: pd.DataFrame, desc: str):
    """Evaluate a subset of model_inventory and return predictions, protein metrics, and coverage."""
    prediction_chunks = []
    protein_metric_rows = []
    coverage_rows = []

    sub_inventory = sub_inventory.copy().reset_index(drop=True)
    if sub_inventory.empty:
        print(f"No models to evaluate for: {desc}")
        return pd.DataFrame(columns=PREDICTION_COLUMNS), pd.DataFrame(columns=METRIC_COLUMNS), pd.DataFrame(columns=COVERAGE_COLUMNS)

    iterator = sub_inventory.itertuples(index=False)
    for row in tqdm(list(iterator), total=len(sub_inventory), desc=desc):
        h5_path = Path(row.h5_path)
        try:
            artifact = load_artifact_from_h5(h5_path)
            selected = list(map(str, artifact.get("selected_features", [])))
            if not selected:
                raise ValueError("No selected_features in artifact")
        except Exception as e:
            print(f"Skipping {row.file}: could not load artifact: {type(e).__name__}: {e}")
            continue

        protein = str(getattr(row, "protein"))
        model = artifact.get("model", None)
        if model is None:
            print(f"Skipping {row.file}: artifact has no model")
            continue

        for dataset_name, (X_ds, Y_ds) in datasets.items():
            if protein not in Y_ds.columns:
                continue
            try:
                Xmat, coverage = apply_artifact_normalizer(X_ds, artifact)
                y_pred = model.predict(Xmat)
            except Exception:
                # Try dataframe fallback for rare sklearn/xgb behavior.
                try:
                    X_df = X_ds.reindex(columns=selected, fill_value=0.0)
                    y_pred = model.predict(X_df)
                    coverage = {
                        "selected_feature_count": len(selected),
                        "missing_feature_count": len([g for g in selected if g not in X_ds.columns]),
                        "missing_feature_fraction": np.nan,
                    }
                except Exception as e:
                    print(f"Skipping {row.file} on {dataset_name}: prediction failed: {type(e).__name__}: {e}")
                    continue

            y_true = Y_ds[protein].to_numpy(dtype=float)
            sample_ids = Y_ds.index.astype(str).to_numpy()
            tissues = [infer_tissue(s) for s in sample_ids]

            pred_df = pd.DataFrame({
                "dataset": dataset_name,
                "sample": sample_ids,
                "tissue": tissues,
                "protein": protein,
                "y_true": y_true,
                "y_pred": np.asarray(y_pred, dtype=float),
                "source_folder": row.source_folder,
                "model_family": row.model_family,
                "normalization": row.normalization,
                "grid_preset": row.grid_preset,
                "select_k": row.select_k,
                "run_label": row.run_label,
                "h5_path": str(h5_path),
            })
            prediction_chunks.append(pred_df)

            m = evaluate_metric_vector(y_true, y_pred)
            metric_row = {
                "dataset": dataset_name,
                "protein": protein,
                "source_folder": row.source_folder,
                "model_family": row.model_family,
                "normalization": row.normalization,
                "grid_preset": row.grid_preset,
                "select_k": row.select_k,
                "run_label": row.run_label,
                "h5_path": str(h5_path),
                "test_mse": m["mse"],
                "test_pearson_corr": m["pearson"],
                "test_cosine_similarity": m["cosine_similarity"],
                "n_test_used": m["n_used"],
                **coverage,
            }
            protein_metric_rows.append(metric_row)
            coverage_rows.append({
                "dataset": dataset_name,
                "protein": protein,
                "source_folder": row.source_folder,
                "run_label": row.run_label,
                "h5_path": str(h5_path),
                **coverage,
            })

    predictions = pd.concat(prediction_chunks, ignore_index=True) if prediction_chunks else pd.DataFrame(columns=PREDICTION_COLUMNS)
    metrics = pd.DataFrame(protein_metric_rows, columns=METRIC_COLUMNS)
    coverage = pd.DataFrame(coverage_rows, columns=COVERAGE_COLUMNS)
    return predictions, metrics, coverage

def extract_labels_from_metadata(metadata, candidates):
    for key in candidates:
        if key in metadata:
            labels = clean_label_list(metadata[key])
            if labels is not None:
                return labels
    return None

def extract_prediction_array_from_dict(d):
    pred_keys = [
        k for k in d
        if "pred" in str(k).lower()
        and hasattr(d[k], "ndim")
        and np.asarray(d[k]).ndim >= 2
    ]
    if pred_keys:
        return np.asarray(d[pred_keys[0]])

    numeric_2d_keys = [
        k for k in d
        if hasattr(d[k], "ndim")
        and np.asarray(d[k]).ndim >= 2
        and np.issubdtype(np.asarray(d[k]).dtype, np.number)
    ]
    if numeric_2d_keys:
        return np.asarray(d[numeric_2d_keys[0]])

    raise ValueError(f"No usable prediction array found in dictionary keys: {list(d)}")

def filter_nans(df, threshold=0.5):
    if df.empty: return df
    # 1. Rows (Features)
    df_clean = df.dropna(thresh=int(df.shape[1] * threshold), axis=0)
    rows_dropped = df.shape[0] - df_clean.shape[0]
    # 2. Cols (Samples)
    df_clean = df_clean.dropna(thresh=int(df_clean.shape[0] * threshold), axis=1)
    cols_dropped = df.shape[1] - df_clean.shape[1]
    print(f"   -> Dropped {rows_dropped} Features & {cols_dropped} Samples")
    return df_clean

def find_proxy_corr_column(df):
    normalized_cols = {normalize_name(c): c for c in df.columns}
    for candidate in PROXY_CORR_CANDIDATES:
        key = normalize_name(candidate)
        if key in normalized_cols:
            return normalized_cols[key]
    return None

def find_sctranslator_prediction_file(dataset_folder, variant="10K"):
    search_roots = [
        Path("."),
        Path(".."),
        Path("../notebooks"),
        OUTPUT_DIR,
        Path("../../notebooks"),
    ]

    candidate_names = [
        f"predictions_{variant}.npz",
        f"predictions_{variant}.npy",
    ]

    checked = []
    for root in search_roots:
        for name in candidate_names:
            path = root / dataset_folder / name
            checked.append(path)
            if path.exists():
                return path

    raise FileNotFoundError(
        "Could not find scTranslator prediction file. Checked:\n"
        + "\n".join(str(p) for p in checked)
    )

def first_existing_column(df, candidates, required=True):
    for col in candidates:
        if col in df.columns:
            return col
    if required:
        raise ValueError(
            "Could not find any of these columns: "
            + ", ".join(candidates)
            + f"\nAvailable columns: {list(df.columns)}"
        )
    return None

def fit_on_ref_apply(
    df_target, df_ref,
    method="robust_z",
    transform=None,              # None | "log2p1" | "asinh"
    clip_nonneg=False,
    sample_center=None,          # None | "median" | "mean"
    clip_to_ref_quantiles=None,  # e.g. (0.001, 0.999)
    n_quantiles=1000,
    random_state=0
):
    # 1) Align features strictly to training (Ref) feature space
    common = df_target.columns.intersection(df_ref.columns)
    if len(common) == 0:
        return None, None, common

    X_t = df_target[common].copy()
    X_r = df_ref[common].copy()

    # 2) Optional per-sample centering (apply to BOTH so pipeline is consistent)
    X_t = _sample_center_df(X_t, sample_center)
    X_r = _sample_center_df(X_r, sample_center)

    # 3) Transform (explicit; avoid "auto" here)
    if transform == "log2p1":
        if clip_nonneg:
            X_t = X_t.clip(lower=0)
            X_r = X_r.clip(lower=0)
        X_t = np.log2(X_t + 1.0)
        X_r = np.log2(X_r + 1.0)
    elif transform == "asinh":
        X_t = np.arcsinh(X_t)
        X_r = np.arcsinh(X_r)
    elif transform is None:
        pass
    else:
        raise ValueError(transform)

    # 4) Impute using REF medians only
    ref_med = X_r.median(axis=0)
    X_r = X_r.fillna(ref_med)
    X_t = X_t.fillna(ref_med)

    # 5) Optional: clip target to REF feature-wise quantiles (prevents huge tails)
    if clip_to_ref_quantiles is not None:
        qlo, qhi = clip_to_ref_quantiles
        lo = X_r.quantile(qlo, axis=0)
        hi = X_r.quantile(qhi, axis=0)
        X_t = X_t.clip(lower=lo, upper=hi, axis=1)

    # 6) Fit on REF, apply to TARGET
    if method == "standard_z":
        scaler = StandardScaler()
        Z_r = scaler.fit_transform(X_r)
        Z_t = scaler.transform(X_t)

    elif method == "robust_z":
        scaler = RobustScaler(with_centering=True, with_scaling=True, quantile_range=(25, 75))
        Z_r = scaler.fit_transform(X_r)
        Z_t = scaler.transform(X_t)

    elif method == "quantile_normal":
        nq = int(min(n_quantiles, max(10, X_r.shape[0])))
        qt = QuantileTransformer(
            n_quantiles=nq,
            output_distribution="normal",
            subsample=int(1e9),
            random_state=random_state
        )
        Q_r = qt.fit_transform(X_r)
        Q_t = qt.transform(X_t)
        scaler = StandardScaler()
        Z_r = scaler.fit_transform(Q_r)
        Z_t = scaler.transform(Q_t)

    else:
        raise ValueError(method)

    # 7) PEAK ALIGNMENT (The "Final Fix")
    # This forces the validation features to have the exact same mean as the training features.
    # It slides the "hill" of the validation data to overlap with the training peak at zero.
    Z_t = Z_t - (np.mean(Z_t, axis=0) - np.mean(Z_r, axis=0))

    # Return as DataFrames
    Z_t_df = pd.DataFrame(Z_t, index=X_t.index, columns=common)
    Z_r_df = pd.DataFrame(Z_r, index=X_r.index, columns=common)
    return Z_t_df, Z_r_df, common

def fmt_num(x, digits=3):
    if pd.isna(x) or not np.isfinite(float(x)):
        return "NA"
    return f"{float(x):.{digits}g}"

def get_dataset_metrics(dataset_key, include_sc=False):
    key = normalize_name(dataset_key)
    stages = BASE_STAGES + (SC_STAGES if include_sc else [])
    out = metrics_partition[
        metrics_partition["dataset_key"].eq(key)
        & metrics_partition["stage"].isin(stages)
    ].copy()
    return out

def get_stats(df, name):
    vals = df.values.flatten()
    return pd.Series({
        "Mean": np.nanmean(vals),
        "Median": np.nanmedian(vals),
        "Std": np.nanstd(vals),
        "Min": np.nanmin(vals),
        "Max": np.nanmax(vals)
    }, name=name)

def get_stats_from_array(vals, name):
    vals = np.asarray(vals)
    vals = vals[~np.isnan(vals)]
    return pd.Series({
        "Mean": np.mean(vals) if vals.size else np.nan,
        "Median": np.median(vals) if vals.size else np.nan,
        "Std": np.std(vals) if vals.size else np.nan,
        "Min": np.min(vals) if vals.size else np.nan,
        "Max": np.max(vals) if vals.size else np.nan
    }, name=name)

def get_tissue(sample_name):
    parts = sample_name.rsplit('_', 1)
    if len(parts) == 2 and parts[1].isdigit(): return parts[0]
    return sample_name

def get_top_100(proteins, values, metric_name='PCC'):
    # Filter NaNs
    valid_mask = ~np.isnan(values)
    prots = np.array(proteins)[valid_mask]
    vals = values[valid_mask]
    
    # Sort
    if metric_name in ['MSE', 'MAE']:
        # Ascending sort: smallest is best
        sorted_idx = np.argsort(vals)
    else:
        # Descending sort: largest is best
        sorted_idx = np.argsort(vals)[::-1]
        
    top_100_idx = sorted_idx[:100]
    return prots[top_100_idx], vals[top_100_idx]

def infer_stage_from_run_label(label):
    s = str(label).lower()
    compact = normalize_name(s)

    if is_sctranslator_label(label):
        if "2m" in compact:
            return "scTranslator 2M"
        if "160k" in compact:
            return "scTranslator 160K"
        if "10k" in compact:
            return "scTranslator 10K"
        return "scTranslator"

    if is_proxy_label(label):
        if "spearman" in compact:
            return "Proxy Spearman"
        if "pearson" in compact:
            return "Proxy Pearson"
        return "Proxy baseline"

    if is_enet_label(label):
        return "ENet"

    if is_xgb_label(label):
        return "XGB"

    return "Other"

def infer_tissue(sample_id: str) -> str:
    """Simple tissue inference from sample IDs. Customize this if you have a better tissue map."""
    s = str(sample_id)
    if "_" in s:
        tail = s.split("_")[-1]
        if 2 <= len(tail) <= 30:
            return tail
    if "-" in s:
        return s.split("-")[0]
    return "Unknown"

def is_baseline_label(label: str) -> bool:
    label = str(label).lower()
    return "same_gene" in label or "baseline" in label

def is_enet_label(label):
    s = str(label).lower()
    return "enet" in s or "elasticnet" in s or "elastic_net" in s

def is_enet_zbig_zscore_k5000_label(label):
    s = str(label).lower()
    compact = normalize_name(s)
    target_compact = normalize_name(TARGET_ENET_RUN_LABEL)

    if compact == target_compact:
        return True

    has_enet = "enet" in compact
    has_zbig = "zbig" in compact or "zgbig" in compact
    has_zscore = "zscore" in compact or "zscored" in compact or "z_score" in s or "z-score" in s
    has_k5000 = (
        "k5000" in compact
        or "top5000" in compact
        or "5000features" in compact
        or re.search(r"(?<!\d)5000(?!\d)", s) is not None
    )

    return bool(has_enet and has_zbig and has_zscore and has_k5000)

def is_proxy_label(label):
    s = str(label).lower()
    compact = normalize_name(s)
    return (
        "same_gene" in s
        or "samegene" in compact
        or "baseline" in s
        or "proxy" in s
        or "rnagene" in compact
        or "samegenebaseline" in compact
    )

def is_sctranslator_label(label):
    s = str(label).lower()
    compact = normalize_name(s)
    return (
        "sctranslator" in compact
        or "sc_translator" in s
        or "predictions_10k" in s
        or "predictions_160k" in s
        or "predictions_2m" in s
    )

def is_xgb_label(label):
    s = str(label).lower()
    return "xgb" in s or "xgboost" in s

def label_family(label: str) -> str:
    # Baseline first in case labels contain overlapping words.
    if is_baseline_label(label):
        return "same-gene baseline"
    if is_xgb_label(label):
        return "XGB"
    if is_enet_label(label):
        return "ENet"
    return "Other"

def load_artifact_from_h5(path):
    with h5py.File(path, "r") as h5:
        payload = h5["model/artifact_pickle_xz"][:]
    return pickle.loads(lzma.decompress(np.asarray(payload, dtype=np.uint8).tobytes()))

def load_eval_predictions_long():
    if PREDICTIONS_PARQUET.exists():
        print("Loading validation predictions from:", PREDICTIONS_PARQUET)
        return pd.read_parquet(PREDICTIONS_PARQUET)

    if PREDICTIONS_CSV.exists():
        print("Loading validation predictions from:", PREDICTIONS_CSV)
        return pd.read_csv(PREDICTIONS_CSV)

    raise FileNotFoundError(
        f"Missing {PREDICTIONS_PARQUET} and {PREDICTIONS_CSV}. "
        "Run the validation notebook first."
    )

def load_feature_list(path):
    if not os.path.exists(path):
        print(f"Warning: {path} not found.")
        return set()
    with open(path, 'r') as f:
        # Assuming one gene per line or comma separated? 
        # Usually these .txt files are one per line based on previous context.
        content = f.read().strip().splitlines()
    return set(x.strip() for x in content if x.strip())

def load_prediction_array(path):
    path = Path(path)

    if path.suffix == ".npy":
        obj = np.load(path, allow_pickle=True)
        if isinstance(obj, np.ndarray) and obj.shape == () and isinstance(obj.item(), dict):
            obj = obj.item()
            return extract_prediction_array_from_dict(obj), obj
        return np.asarray(obj), {}

    if path.suffix == ".npz":
        try:
            import scipy.sparse as sp
            sparse_obj = sp.load_npz(path)
            return sparse_obj.toarray(), {}
        except Exception:
            pass

        z = np.load(path, allow_pickle=True)
        keys = list(z.files)
        arrays = {k: z[k] for k in keys}

        pred_key_candidates = [
            k for k in keys
            if "pred" in k.lower()
            and hasattr(arrays[k], "ndim")
            and arrays[k].ndim >= 2
        ]
        if pred_key_candidates:
            return np.asarray(arrays[pred_key_candidates[0]]), arrays

        numeric_2d_keys = [
            k for k in keys
            if hasattr(arrays[k], "ndim")
            and arrays[k].ndim >= 2
            and np.issubdtype(arrays[k].dtype, np.number)
        ]
        if numeric_2d_keys:
            return np.asarray(arrays[numeric_2d_keys[0]]), arrays

        raise ValueError(
            f"No 2D numeric prediction array found in {path}. Available keys: {keys}"
        )

    raise ValueError(f"Unsupported prediction file format: {path}")

def load_val_data(prefix: str):
    rna_path = VAL_DATA_DIR / f"{prefix}_transcriptomics.joblib"
    prot_path = VAL_DATA_DIR / f"{prefix}_proteomics.joblib"

    if not rna_path.exists() or not prot_path.exists():
        print(f"Skipping {prefix}: missing {rna_path.name} or {prot_path.name}")
        return None

    X = joblib.load(rna_path)
    Y = joblib.load(prot_path)

    X = X.copy()
    Y = Y.copy()
    X.index = X.index.map(str)
    X.columns = X.columns.map(str)
    Y.index = Y.index.map(str)
    Y.columns = Y.columns.map(str)

    # Proteomics should be rows=samples, columns=proteins. If needed, transpose.
    overlap_rows = len(set(X.index).intersection(set(Y.index)))
    overlap_cols = len(set(X.index).intersection(set(Y.columns)))
    if overlap_cols > overlap_rows:
        Y = Y.T.copy()
        Y.index = Y.index.map(str)
        Y.columns = Y.columns.map(str)

    common = X.index.intersection(Y.index)
    X = X.loc[common].copy()
    Y = Y.loc[common].copy()

    print(f"Loaded {prefix:10s}: X={X.shape}, Y={Y.shape}, common samples={len(common)}")
    return X, Y

def mae(y_true, y_pred): return mean_absolute_error(y_true, y_pred)

def make_family_best_tables_for_dataset(dataset_key, source_df=None):
    key = normalize_name(dataset_key)
    df = metrics_partition if source_df is None else source_df
    ds = df[df["dataset_key"].eq(key)].copy()

    base = ds[ds["stage"].isin(BASE_STAGES)].copy()
    proxy = ds[ds["stage"].isin(PROXY_STAGES)].copy()
    enet = ds[ds["stage"].eq("ENet")].copy()
    xgb = ds[ds["stage"].eq("XGB")].copy()
    sc = ds[ds["stage"].isin(SC_STAGES)].copy()

    universe_source = base if len(base) else ds
    protein_universe = pd.Index(sorted(universe_source["protein"].dropna().astype(str).unique()))

    tables = {
        "dataset_key": key,
        "dataset_display": dataset_display_name(key),
        "protein_universe": protein_universe,
        "proxy_best_pearson": best_per_protein(proxy, "test_pearson_corr", higher=True),
        "enet_best_pearson": best_per_protein(enet, "test_pearson_corr", higher=True),
        "xgb_best_pearson": best_per_protein(xgb, "test_pearson_corr", higher=True),
        "sc_best_pearson": best_per_protein(sc, "test_pearson_corr", higher=True),
        "proxy_best_cosine": best_per_protein(proxy, "test_cosine_similarity", higher=True),
        "enet_best_cosine": best_per_protein(enet, "test_cosine_similarity", higher=True),
        "xgb_best_cosine": best_per_protein(xgb, "test_cosine_similarity", higher=True),
        "sc_best_cosine": best_per_protein(sc, "test_cosine_similarity", higher=True),
        "proxy_best_mse": best_per_protein(proxy, "test_mse", higher=False),
        "enet_best_mse": best_per_protein(enet, "test_mse", higher=False),
        "xgb_best_mse": best_per_protein(xgb, "test_mse", higher=False),
        "sc_best_mse": best_per_protein(sc, "test_mse", higher=False),
    }

    return tables

def make_partition_table(
    best_df,
    metric,
    threshold,
    higher=True,
    min_proteins_to_show=200,
):
    columns = [
        "partition_label",
        "proteins",
        "median_metric",
        "share",
        "is_bad",
        "is_other",
    ]

    df = best_df.dropna(subset=[metric]).copy()

    if df.empty or threshold is None or not np.isfinite(float(threshold)):
        empty = pd.DataFrame(columns=columns)
        return empty, empty, "Bad"

    threshold = float(threshold)
    pretty_metric = metric_name(metric)

    if higher:
        df["is_good"] = df[metric] >= threshold
        bad_label = f"Bad / {pretty_metric} < {fmt_num(threshold)}"
    else:
        df["is_good"] = df[metric] <= threshold
        bad_label = f"Bad / {pretty_metric} > {fmt_num(threshold)}"

    df["run_label"] = df["run_label"].astype(str)
    df["partition_label_raw"] = np.where(
        df["is_good"],
        df["run_label"],
        bad_label,
    )

    raw_part = (
        df.groupby("partition_label_raw", dropna=False)
        .agg(
            proteins=("protein", "nunique"),
            median_metric=(metric, "median"),
        )
        .reset_index()
        .rename(columns={"partition_label_raw": "partition_label"})
    )

    raw_part["is_bad"] = raw_part["partition_label"].eq(bad_label)
    raw_part["is_other"] = False

    small_good_labels = raw_part.loc[
        (~raw_part["is_bad"])
        & (raw_part["proteins"] < int(min_proteins_to_show)),
        "partition_label",
    ].tolist()

    df["partition_label"] = df["partition_label_raw"]

    if small_good_labels:
        df.loc[
            df["partition_label_raw"].isin(small_good_labels),
            "partition_label",
        ] = f"Other good models (<{min_proteins_to_show} proteins each)"

    part = (
        df.groupby("partition_label", dropna=False)
        .agg(
            proteins=("protein", "nunique"),
            median_metric=(metric, "median"),
        )
        .reset_index()
    )

    part["is_bad"] = part["partition_label"].eq(bad_label)
    part["is_other"] = part["partition_label"].str.startswith("Other good models")

    total = part["proteins"].sum()
    part["share"] = part["proteins"] / total if total else np.nan
    raw_part["share"] = raw_part["proteins"] / total if total else np.nan

    part["sort_key"] = np.select(
        [part["is_bad"], part["is_other"]],
        [2, 1],
        default=0,
    )

    part = (
        part.sort_values(["sort_key", "proteins"], ascending=[True, False])
        .drop(columns="sort_key")
        .reset_index(drop=True)
    )

    raw_part = (
        raw_part.sort_values(["is_bad", "proteins"], ascending=[True, False])
        .reset_index(drop=True)
    )

    return part, raw_part, bad_label

def metric_name(metric):
    return METRIC_NAMES.get(metric, metric.replace("test_", "").replace("_", " "))

def mse(y_true, y_pred): return mean_squared_error(y_true, y_pred)

def mse_procedural_improvement_assignments(tables):
    dataset = tables["dataset_display"]
    protein_universe = pd.Index(tables["protein_universe"]).astype(str)

    proxy_idx = _index_by_protein(tables["proxy_best_mse"])
    enet_idx = _index_by_protein(tables["enet_best_mse"])
    xgb_idx = _index_by_protein(tables["xgb_best_mse"])

    rows = []
    for protein in protein_universe:
        proxy_mse = proxy_idx.loc[protein, "test_mse"] if protein in proxy_idx.index else np.nan
        enet_mse = enet_idx.loc[protein, "test_mse"] if protein in enet_idx.index else np.nan
        xgb_mse = xgb_idx.loc[protein, "test_mse"] if protein in xgb_idx.index else np.nan

        valid_proxy = np.isfinite(proxy_mse)
        valid_enet = np.isfinite(enet_mse)
        valid_xgb = np.isfinite(xgb_mse)

        if valid_proxy:
            current_best = proxy_mse
            stage = "Proxy retained / no lower MSE found"
            if valid_enet and enet_mse < current_best:
                current_best = enet_mse
                stage = "ENet improves over proxy"
            if valid_xgb and xgb_mse < current_best:
                current_best = xgb_mse
                stage = "XGB improves over previous best"
        elif valid_enet:
            current_best = enet_mse
            stage = "ENet first valid MSE"
            if valid_xgb and xgb_mse < current_best:
                current_best = xgb_mse
                stage = "XGB improves over previous best"
        elif valid_xgb:
            current_best = xgb_mse
            stage = "XGB first valid MSE"
        else:
            current_best = np.nan
            stage = "No valid comparison"

        rows.append({
            "dataset": dataset,
            "protein": protein,
            "threshold": 0.0,
            "stage": stage,
            "best_mse": current_best,
            "proxy_mse": proxy_mse,
            "enet_mse": enet_mse,
            "xgb_mse": xgb_mse,
        })

    return pd.DataFrame(rows)

def mse_winner_assignments(tables):
    dataset = tables["dataset_display"]
    protein_universe = pd.Index(tables["protein_universe"]).astype(str)

    proxy_idx = _index_by_protein(tables["proxy_best_mse"])
    enet_idx = _index_by_protein(tables["enet_best_mse"])
    xgb_idx = _index_by_protein(tables["xgb_best_mse"])

    rows = []
    for protein in protein_universe:
        candidates = []
        for family, idx in [("proxy", proxy_idx), ("enet", enet_idx), ("xgb", xgb_idx)]:
            if protein in idx.index:
                row = idx.loc[protein]
                mse = safe_metric(row, "test_mse")
                if np.isfinite(mse):
                    candidates.append((family, mse, row))

        if not candidates:
            stage = "No valid MSE"
            value = np.nan
            run_label = ""
        else:
            family, value, row = min(candidates, key=lambda x: x[1])
            if family == "proxy":
                stage = f"{row.get('stage', 'Proxy baseline')} lowest MSE"
            elif family == "enet":
                stage = "ENet lowest MSE"
            else:
                stage = "XGB lowest MSE"
            run_label = row.get("run_label", "")

        rows.append({
            "dataset": dataset,
            "protein": protein,
            "threshold": 0.0,
            "stage": stage,
            "test_mse": value,
            "run_label": run_label,
        })

    return pd.DataFrame(rows)

def normalize_dataset_name(x):
    return re.sub(r"[^a-z0-9]", "", str(x).lower())

def normalize_name(x):
    return re.sub(r"[^a-z0-9]", "", str(x).lower())

def parse_model_filename(path: Path):
    stem = path.stem
    protein = stem
    model_family = "unknown"
    normalization = "unknown"
    grid_preset = "unknown"
    select_k = np.nan

    if "_xgb_8020_" in stem:
        protein = stem.split("_xgb_8020_")[0]
        model_family = "xgb"
        rest = stem.split("_xgb_8020_")[1]
    elif "_enet_8020_" in stem:
        protein = stem.split("_enet_8020_")[0]
        model_family = "enet"
        rest = stem.split("_enet_8020_")[1]
    else:
        rest = stem

    toks = rest.split("_")
    known_norms = {"none", "zscore", "rank", "ranking", "binning", "qbin5", "qbin10", "qbin20", "qbin50", "qbin100"}
    for tok in toks:
        if tok in known_norms or re.fullmatch(r"qbin\d+", tok):
            normalization = tok
            break
    for tok in toks:
        if tok in {"tiny", "fast", "current"}:
            grid_preset = tok
            break
    for tok in toks:
        m = re.fullmatch(r"k(\d+)", tok)
        if m:
            select_k = int(m.group(1))
            break

    return protein, model_family, normalization, grid_preset, select_k

def pearson_correlation(y_true, y_pred): 
     if len(y_true) < 2: return 0.0
     val, _ = pearsonr(y_true, y_pred)
     return val

def per_protein_metrics(y_true_matrix, y_pred_matrix, proteins):
    rows = []
    y_true_matrix = np.asarray(y_true_matrix, dtype=float)
    y_pred_matrix = np.asarray(y_pred_matrix, dtype=float)

    for j, protein in enumerate(proteins):
        yt = y_true_matrix[:, j]
        yp = y_pred_matrix[:, j]
        mask = np.isfinite(yt) & np.isfinite(yp)
        n = int(mask.sum())

        if n == 0:
            pearson = np.nan
            cosine = np.nan
            mse = np.nan
        else:
            yt_valid = yt[mask]
            yp_valid = yp[mask]
            mse = float(np.mean((yt_valid - yp_valid) ** 2))

            if n >= 2 and np.std(yt_valid) > 0 and np.std(yp_valid) > 0:
                pearson = float(np.corrcoef(yt_valid, yp_valid)[0, 1])
            else:
                pearson = np.nan

            denom = np.linalg.norm(yt_valid) * np.linalg.norm(yp_valid)
            cosine = float(np.dot(yt_valid, yp_valid) / denom) if denom > 0 else np.nan

        rows.append({
            "protein": protein,
            "n_samples": n,
            "test_pearson_corr": pearson,
            "test_cosine_similarity": cosine,
            "test_mse": mse,
        })

    return pd.DataFrame(rows)

def plot_detailed(per_protein_df, top_10, bottom_10, title, out_path):
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)

    # Reusing your style with the new "Per-Protein" logic
    metrics_to_plot = [
        ("pearson", "#F5B041", "Pearson Correlation"),
        ("cosine", "#5DADE2", "Cosine Similarity"),
        ("mse", "#58D68D", "Mean Squared Error"),
        ("spearman", "#AF7AC5", "Spearman Correlation")
    ]

    for i, (col, color, label) in enumerate(metrics_to_plot):
        row, cp = divmod(i, 2) # Arrange in a grid
        ax = fig.add_subplot(gs[row, cp])
        ax.hist(per_protein_df[col].dropna(), bins=30, color=color, edgecolor="black", alpha=0.8)
        mean_val = per_protein_df[col].mean()
        ax.axvline(mean_val, color="red", linestyle="--", label=f"Mean: {mean_val:.3f}")
        ax.set_title(f"Per-Protein {label}", fontweight="bold")
        ax.legend()

    # Top/Bottom 10 (Using Pearson as the anchor)
    ax5 = fig.add_subplot(gs[0, 2])
    combined = pd.concat([top_10.sort_values("pearson"), bottom_10.sort_values("pearson", ascending=False)])
    colors = ["#58D68D"] * 10 + ["#E74C3C"] * 10
    ax5.barh(range(len(combined)), combined["pearson"], color=colors, edgecolor="black")
    ax5.set_yticks(range(len(combined)))
    ax5.set_yticklabels(combined["protein"], fontsize=8)
    ax5.set_title("Top/Bottom Proteins (Pearson)", fontweight="bold")

    # Stats Summary Box
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.axis("off")
    stats = per_protein_df.describe().loc[['mean', 'std', 'min', 'max']]
    
    summary_text = "Global Protein Stats\n" + "="*20 + "\n"
    summary_text += f"Total Proteins: {len(per_protein_df)}\n\n"
    summary_text += f"Mean Pearson: {stats.loc['mean', 'pearson']:.4f}\n"
    summary_text += f"Mean Cosine:  {stats.loc['mean', 'cosine']:.4f}\n"
    summary_text += f"Mean MSE:     {stats.loc['mean', 'mse']:.4f}\n"
    
    ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes, fontsize=10, 
             fontfamily="monospace", verticalalignment="top", 
             bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    fig.suptitle(title, fontsize=16, fontweight="bold", y=0.98)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.show()

def plot_normalized_distributions(data_store, ccle_prot, ccle_rna):
    sns.set_style("white")
    fig, axes = plt.subplots(3, 2, figsize=(14, 12), dpi=160)  # higher DPI helps smooth rendering
    fig.suptitle("Final Alignment: Normalized Collection A vs CCLE",
                 fontsize=20, fontweight='bold', y=0.98)

    palette = {"Val": "#3498db", "Ref": "#e67e22"}
    GROUPS = ['TMT', 'iTRAQ', 'LabelFree']
    TYPES = ['Protein', 'RNA']

    # ---- Set explicitly (MOST IMPORTANT) ----
    # If CCLE RNA is already log2(TPM+1): RNA_TRANSFORM = None
    # If CCLE RNA is TPM/CPM (not logged): RNA_TRANSFORM = "log2p1"
    RNA_TRANSFORM = None
    RNA_METHOD = "robust_z"
    RNA_CLIP_Q = (0.001, 0.999)  # reduce crazy tails that create spikes/angles

    # Protein
    PROT_TRANSFORM = "asinh"         # safe when ref has negatives/ratios
    PROT_METHOD_TMT_ITRAQ = "standard_z"
    PROT_METHOD_LABELFREE = "robust_z"
    LABELFREE_SAMPLE_CENTER = "median"

    for r, group in enumerate(GROUPS):
        for c, dtype in enumerate(TYPES):
            ax = axes[r][c]
            df_my = data_store[group][dtype]
            df_ref = ccle_prot if dtype == "Protein" else ccle_rna

            if df_my is None:
                ax.text(0.5, 0.5, "No Data", ha='center', color='gray')
                ax.set_title(f"{group} - {dtype}", fontsize=14, fontweight='semibold')
                continue

            # Normalize
            if dtype == "Protein":
                method = PROT_METHOD_LABELFREE if group == "LabelFree" else PROT_METHOD_TMT_ITRAQ
                sample_center = LABELFREE_SAMPLE_CENTER if group == "LabelFree" else None

                Z_t_df, Z_r_df, common = fit_on_ref_apply(
                    df_target=df_my,
                    df_ref=df_ref,
                    method=method,
                    transform=PROT_TRANSFORM,
                    clip_nonneg=False,
                    sample_center=sample_center,
                    clip_to_ref_quantiles=None
                )
            else:
                Z_t_df, Z_r_df, common = fit_on_ref_apply(
                    df_target=df_my,
                    df_ref=df_ref,
                    method=RNA_METHOD,
                    transform=RNA_TRANSFORM,
                    clip_nonneg=(RNA_TRANSFORM == "log2p1"),
                    sample_center=None,
                    clip_to_ref_quantiles=RNA_CLIP_Q
                )

            if Z_t_df is None:
                ax.text(0.5, 0.5, "No common features", ha='center', color='gray')
                ax.set_title(f"{group} - {dtype}", fontsize=14, fontweight='semibold')
                continue

            vals_my = Z_t_df.values.flatten()
            vals_ref = Z_r_df.values.flatten()

            # 1) Visualization (smooth KDE)
            _kde_smooth(
                ax,
                vals_my=vals_my,
                vals_ref=vals_ref,
                lab_my=f"{group} (Normalized)",
                lab_ref="CCLE Ref (fit)",
                palette=palette,
                n=20000,
                seed=0,
                bw_adjust=1.35,
                gridsize=1024,
                cut=0,
                jitter_if_discrete=True
            )

            ax.set_title(f"{group} - {dtype}", fontsize=14, fontweight='semibold')
            ax.legend(frameon=True, facecolor='white', loc='upper right')

            # Dynamic X-limits like your RAW plot (but use random sample to avoid bias)
            rng = np.random.default_rng(0)
            v1 = vals_my[~np.isnan(vals_my)]
            v2 = vals_ref[~np.isnan(vals_ref)]
            if v1.size and v2.size:
                s1 = rng.choice(v1, size=min(5000, v1.size), replace=False)
                s2 = rng.choice(v2, size=min(5000, v2.size), replace=False)
                all_vals = np.concatenate([s1, s2])
                ax.set_xlim(np.nanpercentile(all_vals, 0.5), np.nanpercentile(all_vals, 99.5))

            # 2) Stats & correlation logging (RAW style)
            print(f"\n>>> NORMALIZED STATS: {group} - {dtype}")
            s1 = get_stats_from_array(vals_my, "Validation_Norm")
            s2 = get_stats_from_array(vals_ref, "CCLE_Ref_Norm")
            print(pd.concat([s1, s2], axis=1))

            # Feature-mean correlation after normalization
            if len(common) > 10:
                m1 = Z_t_df.mean(axis=0)
                m2 = Z_r_df.mean(axis=0)
                corr, _ = pearsonr(m1.values, m2.values)
                print(f"Feature Mean Correlation (Pearson): {corr:.4f}")

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

def plot_overall_bar(mean_metrics, title, out_path):
    plt.style.use("seaborn-v0_8-whitegrid")

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.patch.set_facecolor("#F8F9FA")

    fig.suptitle(title, fontsize=20, fontweight="bold", color="#2C3E50", y=0.98)

    colors = ["#E74C3C", "#E67E22", "#F39C12", "#27AE60", "#2980B9", "#9B59B6"]
    # Keys from the user's new compute_mean_metrics
    plots = [
        ("Mean Squared Error", "MSE", mean_metrics.get("MSE", 0), colors[0]),
        ("Root MSE", "RMSE", mean_metrics.get("RMSE", 0), colors[1]),
        ("Mean Absolute Error", "MAE", mean_metrics.get("MAE", 0), colors[2]),
        ("Pearson Correlation", "Pearson", mean_metrics.get("Pearson", 0), colors[3]),
        ("Spearman Correlation", "Spearman", mean_metrics.get("Spearman", 0), colors[4]),
        ("Cosine Similarity", "Cosine", mean_metrics.get("Cosine", 0), colors[5]),
    ]

    for idx, (t, ylabel, value, color) in enumerate(plots):
        row, col = divmod(idx, 3)
        ax = axes[row, col]
        ax.set_facecolor("#FFFFFF")

        if np.isnan(value): value = 0.0

        bar = ax.bar(["Model"], [value], color=color, width=0.6,
                     edgecolor="#2C3E50", linewidth=2, alpha=0.9)

        ax.bar(["Model"], [value], color="gray", width=0.6, alpha=0.1,
               bottom=-value * 0.02 if value > 0 else value * 0.02)

        for rect in bar:
            height = rect.get_height()
            ax.annotate(
                f"{value:.4f}",
                xy=(rect.get_x() + rect.get_width() / 2, height),
                xytext=(0, 8),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=14,
                fontweight="bold",
                color="#2C3E50",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=color, alpha=0.8),
            )

        ax.set_title(t, fontsize=14, fontweight="bold", color="#2C3E50", pad=15)
        ax.set_ylabel(ylabel, fontsize=12, fontweight="bold", color="#2C3E50")

        # Dynamic ylim
        if abs(value) < 1e-9:
             ax.set_ylim(0, 1)
        elif value >= 0:
            ax.set_ylim(0, value * 1.35)
        else:
            ax.set_ylim(value * 1.35, abs(value) * 0.35)

        ax.tick_params(axis="x", labelsize=12, colors="#2C3E50")
        ax.tick_params(axis="y", labelsize=10, colors="#7F8C8D")
        ax.yaxis.grid(True, linestyle="--", alpha=0.3)
        ax.set_axisbelow(True)

        for spine in ax.spines.values():
            spine.set_color("#BDC3C7")
            spine.set_linewidth(1.5)
        ax.spines["bottom"].set_visible(False)

    plt.tight_layout(rect=[0, 0.02, 1, 0.95])
    plt.savefig(out_path, dpi=200, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(f"Saved: {out_path}")
    plt.show()

def plot_partition_bar(
    best_df,
    metric,
    threshold,
    higher=True,
    title=None,
    min_proteins_to_show=200,
    figsize=(14, 3.4),
    show_tables=True,
):
    part, raw_part, bad_label = make_partition_table(
        best_df,
        metric,
        threshold,
        higher=higher,
        min_proteins_to_show=min_proteins_to_show,
    )

    total = int(part["proteins"].sum()) if len(part) else 0

    if total == 0:
        print(f"No data for {metric}")
        return part

    bad_n = int(part.loc[part["partition_label"].eq(bad_label), "proteins"].sum())
    good_n = total - bad_n

    pretty_metric = metric_name(metric)
    threshold_text = fmt_num(threshold)

    if higher:
        good_rule = f"{pretty_metric} ≥ {threshold_text}"
    else:
        good_rule = f"{pretty_metric} ≤ {threshold_text}"

    print_family_summary(best_df, metric, threshold, higher=higher)

    if show_tables:
        visible = part.copy()
        visible["share_percent"] = (100 * visible["share"]).round(1)

        print("\nVisible bar segments:")
        display(
            visible[
                ["partition_label", "proteins", "share_percent", "median_metric"]
            ]
        )

        hidden = raw_part[
            (~raw_part["is_bad"])
            & (raw_part["proteins"] < min_proteins_to_show)
        ].copy()

        if len(hidden):
            hidden["share_percent"] = (100 * hidden["share"]).round(1)
            print(f"\nHidden individual good segments with < {min_proteins_to_show} proteins:")
            display(
                hidden[
                    ["partition_label", "proteins", "share_percent", "median_metric"]
                ]
            )

    base_colors = [
        "#3B6EA8", "#4E9F8A", "#6BA84F", "#F2A65A", "#9D79BC",
        "#7A6F9B", "#D6A77A", "#8FB8DE", "#A1C181",
    ]

    bad_color = "#D1495B"
    other_color = "#9E9E9E"

    color_map = {}
    color_i = 0

    for label in part["partition_label"]:
        if label == bad_label:
            color_map[label] = bad_color
        elif str(label).startswith("Other good models"):
            color_map[label] = other_color
        else:
            color_map[label] = base_colors[color_i % len(base_colors)]
            color_i += 1

    fig, ax = plt.subplots(figsize=figsize)

    subtitle = (
        f"Good rule: {good_rule}   |   "
        f"Good: {good_n:,} ({100 * good_n / total:.1f}%)   |   "
        f"Bad: {bad_n:,} ({100 * bad_n / total:.1f}%)   |   "
        f"Total proteins: {total:,}"
    )

    fig.suptitle(
        title or f"Best model per protein by {pretty_metric}",
        y=0.98,
        fontsize=14,
        fontweight="bold",
    )

    fig.text(
        0.5,
        0.84,
        subtitle,
        ha="center",
        va="center",
        fontsize=10,
        color="0.25",
    )

    left = 0
    handles = []
    legend_labels = []

    for _, row in part.iterrows():
        width = int(row["proteins"])
        label = row["partition_label"]
        pct = 100 * width / total
        color = color_map[label]

        bar = ax.barh(
            y=0,
            width=width,
            left=left,
            height=0.48,
            color=color,
            edgecolor="white",
            linewidth=1.6,
        )

        if pct >= 5:
            ax.text(
                left + width / 2,
                0,
                f"{width:,}\n{pct:.1f}%",
                ha="center",
                va="center",
                color="white",
                fontsize=10,
                fontweight="bold",
            )

        med = row["median_metric"]
        short_label = textwrap.shorten(str(label), width=70, placeholder="...")

        legend_label = f"{short_label} — {width:,} ({pct:.1f}%)"
        if pd.notna(med) and np.isfinite(float(med)):
            legend_label += f", median={fmt_num(med)}"

        handles.append(bar[0])
        legend_labels.append(legend_label)

        left += width

    ax.set_xlim(0, total)
    ax.set_ylim(-0.55, 0.55)
    ax.set_yticks([])
    ax.set_xlabel(f"Proteins, total = {total:,}", fontsize=10)
    ax.grid(axis="x", alpha=0.18)
    ax.grid(axis="y", visible=False)
    ax.set_axisbelow(True)

    ax.spines[["top", "right", "left"]].set_visible(False)

    ax.legend(
        handles,
        legend_labels,
        title="Segments",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=9,
        title_fontsize=10,
    )

    fig.subplots_adjust(top=0.72, bottom=0.24, right=0.68)
    plt.show()

    return part

def plot_raw_distributions(data_store, ccle_prot, ccle_rna):
    sns.set_style("white")
    fig, axes = plt.subplots(3, 2, figsize=(14, 12), dpi=100)
    fig.suptitle("Raw Distribution: Collection A vs CCLE", fontsize=20, fontweight='bold', y=0.98)

    palette = {"Val": "#3498db", "Ref": "#e67e22"}

    for r, group in enumerate(GROUPS):
        for c, dtype in enumerate(TYPES):
            ax = axes[r][c]
            df_my = data_store[group][dtype]
            df_ref = ccle_prot if dtype == "Protein" else ccle_rna
            
            if df_my is not None:
                vals_my = df_my.values.flatten()
                vals_ref = df_ref.values.flatten()
                
                # 1. Visualization
                sns.kdeplot(vals_my[~np.isnan(vals_my)][:20000], ax=ax, 
                            fill=True, color=palette["Val"], alpha=0.4, label=f"{group} (Raw)")
                sns.kdeplot(vals_ref[~np.isnan(vals_ref)][:20000], ax=ax, 
                            fill=True, color=palette["Ref"], alpha=0.2, label="CCLE Ref", ls='--')
                
                ax.set_title(f"{group} - {dtype}", fontsize=14, fontweight='semibold')
                ax.legend(frameon=True, facecolor='white', loc='upper right')
                
                # Dynamic X-limits to ignore extreme outliers
                all_vals = np.concatenate([vals_my[:5000], vals_ref[:5000]])
                ax.set_xlim(np.nanpercentile(all_vals, 0.5), np.nanpercentile(all_vals, 99.5))

                # 2. Stats & Correlation Logging
                print(f"\n>>> RAW STATS: {group} - {dtype}")
                s1 = get_stats(df_my, "Validation")
                s2 = get_stats(df_ref, "CCLE_Ref")
                print(pd.concat([s1, s2], axis=1))
                
                common = df_my.columns.intersection(df_ref.columns)
                if len(common) > 10:
                    m1 = df_my[common].mean()
                    m2 = df_ref[common].mean()
                    corr, _ = pearsonr(m1, m2)
                    print(f"Gene Mean Correlation (Pearson): {corr:.4f}")
            else:
                ax.text(0.5, 0.5, "No Data", ha='center', color='gray')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

def plot_single_partition(assignments, title, order, save_name=None):
    tmp = assignments.copy()
    tmp["threshold"] = 0.0
    summary = summarize_assignments(tmp, order=order)

    dataset = str(tmp["dataset"].iloc[0]) if len(tmp) else "dataset"
    plot_stacked_partition_rows(
        summary,
        title=title,
        order=order,
        xlabel="Proteins",
        save_name=save_name,
        figsize=(14, 3.2),
        row_label_func=lambda _: "All proteins",
    )

    display(
        summary[["dataset", "stage", "proteins", "percent"]]
        .sort_values(["dataset", "stage"])
        .reset_index(drop=True)
    )

    return summary

def plot_stacked_partition_rows(
    summary,
    title,
    order,
    xlabel="Proteins",
    save_name=None,
    figsize=None,
    row_label_func=None,
):
    if summary is None or len(summary) == 0:
        print("No summary to plot.")
        return None

    plot_df = summary.copy()
    thresholds = list(plot_df["threshold"].drop_duplicates())
    if figsize is None:
        figsize = (14, max(3.0, 0.62 * len(thresholds) + 1.8))

    fig, ax = plt.subplots(figsize=figsize)
    y_positions = np.arange(len(thresholds))
    max_total = int(plot_df.groupby("threshold")["total_proteins"].first().max())

    used_labels = set()
    for y, thr in zip(y_positions, thresholds):
        sub = plot_df[plot_df["threshold"].eq(thr)].copy()
        total = int(sub["total_proteins"].iloc[0]) if len(sub) else 0
        left = 0

        for stage in order:
            row = sub[sub["stage"].astype(str).eq(stage)]
            width = int(row["proteins"].iloc[0]) if len(row) else 0
            if width <= 0:
                continue

            pct = 100 * width / total if total else 0
            color = COLOR_MAP.get(stage, None)
            label = stage if stage not in used_labels else None

            ax.barh(
                y,
                width,
                left=left,
                color=color,
                edgecolor="white",
                linewidth=1.2,
                label=label,
            )
            used_labels.add(stage)

            if pct >= 6:
                ax.text(
                    left + width / 2,
                    y,
                    f"{width}\n{pct:.0f}%",
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=9,
                    weight="bold",
                )
            elif width > 0 and total > 0:
                ax.text(
                    left + width / 2,
                    y + 0.27,
                    f"{width}",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    color=color,
                )

            left += width

    ax.set_yticks(y_positions)
    if row_label_func is None:
        row_labels = [f"Threshold ≥ {thr:.2f}" for thr in thresholds]
    else:
        row_labels = [row_label_func(thr) for thr in thresholds]
    ax.set_yticklabels(row_labels)
    ax.set_xlim(0, max_total)
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=14, weight="bold", pad=12)
    ax.grid(axis="x", alpha=0.18)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.invert_yaxis()
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=3,
        frameon=False,
        fontsize=9,
    )
    plt.tight_layout()

    if save_name and SAVE_OUTPUTS:
        fig.savefig(FIGURE_DIR / save_name, dpi=180, bbox_inches="tight")

    plt.show()
    return fig

def print_family_summary(best_df, metric, threshold, higher=True):
    df = best_df.dropna(subset=[metric]).copy()

    if df.empty or threshold is None or not np.isfinite(float(threshold)):
        print(f"No valid values for {metric}")
        return

    threshold = float(threshold)
    df["is_good"] = df[metric] >= threshold if higher else df[metric] <= threshold
    df["family"] = df["run_label"].map(label_family)

    total = df["protein"].nunique()
    good = df[df["is_good"]].copy()
    bad = df[~df["is_good"]].copy()

    family_counts = (
        good.groupby("family")["protein"]
        .nunique()
        .reindex(["XGB", "ENet", "same-gene baseline", "Other"], fill_value=0)
    )

    print(f"\n{metric_name(metric)} family summary")
    print("-" * 80)
    print(f"Total proteins with valid metric: {total:,}")
    print(f"Good proteins: {good['protein'].nunique():,} ({100 * good['protein'].nunique() / total:.1f}%)")
    print(f"Bad proteins:  {bad['protein'].nunique():,} ({100 * bad['protein'].nunique() / total:.1f}%)")
    print(f"Good proteins best predicted by XGB:                {int(family_counts['XGB']):,}")
    print(f"Good proteins best predicted by ENet:               {int(family_counts['ENet']):,}")
    print(f"Good proteins best predicted by same-gene baseline: {int(family_counts['same-gene baseline']):,}")

    if int(family_counts["Other"]) > 0:
        print(f"Good proteins best predicted by other labels:       {int(family_counts['Other']):,}")

def procedural_metric_cascade(
    tables,
    threshold,
    metric="test_pearson_corr",
    higher=True,
    include_sc=False,
):
    dataset = tables["dataset_display"]
    protein_universe = pd.Index(tables["protein_universe"]).astype(str)
    remaining = set(protein_universe)
    rows = []

    suffix = {
        "test_pearson_corr": "pearson",
        "test_cosine_similarity": "cosine",
        "test_mse": "mse",
    }[metric]

    proxy_idx = _index_by_protein(tables[f"proxy_best_{suffix}"])
    enet_idx = _index_by_protein(tables[f"enet_best_{suffix}"])
    xgb_idx = _index_by_protein(tables[f"xgb_best_{suffix}"])
    sc_idx = _index_by_protein(tables[f"sc_best_{suffix}"]) if include_sc else pd.DataFrame().rename_axis("protein")

    def passes(value):
        if not np.isfinite(value):
            return False
        return value >= threshold if higher else value <= threshold

    # 1) Proxy/baseline first.
    for protein in list(remaining):
        if protein in proxy_idx.index:
            row = proxy_idx.loc[protein]
            value = safe_metric(row, metric)
            if passes(value):
                rows.append({
                    "dataset": dataset,
                    "protein": protein,
                    "threshold": threshold,
                    "stage": row.get("stage", "Proxy baseline"),
                    "metric": value,
                    "run_label": row.get("run_label", ""),
                })
                remaining.remove(protein)

    # 2) ENet rescues only remaining proteins.
    for protein in list(remaining):
        if protein in enet_idx.index:
            row = enet_idx.loc[protein]
            value = safe_metric(row, metric)
            if passes(value):
                rows.append({
                    "dataset": dataset,
                    "protein": protein,
                    "threshold": threshold,
                    "stage": "ENet",
                    "metric": value,
                    "run_label": row.get("run_label", ""),
                })
                remaining.remove(protein)

    # 3) XGB rescues only remaining proteins.
    for protein in list(remaining):
        if protein in xgb_idx.index:
            row = xgb_idx.loc[protein]
            value = safe_metric(row, metric)
            if passes(value):
                rows.append({
                    "dataset": dataset,
                    "protein": protein,
                    "threshold": threshold,
                    "stage": "XGB",
                    "metric": value,
                    "run_label": row.get("run_label", ""),
                })
                remaining.remove(protein)

    # 4) Optional scTranslator final rescue.
    if include_sc:
        for protein in list(remaining):
            if protein in sc_idx.index:
                row = sc_idx.loc[protein]
                value = safe_metric(row, metric)
                if passes(value):
                    rows.append({
                        "dataset": dataset,
                        "protein": protein,
                        "threshold": threshold,
                        "stage": row.get("stage", "scTranslator"),
                        "metric": value,
                        "run_label": row.get("run_label", ""),
                    })
                    remaining.remove(protein)

    final_label = "Still not good after scTranslator" if include_sc else "Still not good"
    for protein in sorted(remaining):
        rows.append({
            "dataset": dataset,
            "protein": protein,
            "threshold": threshold,
            "stage": final_label,
            "metric": np.nan,
            "run_label": "",
        })

    return pd.DataFrame(rows)

def proxy_corr_threshold_assignments_for_dataset(
    tables,
    source_df,
    pred_thresholds=PEARSON_GOOD_THRESHOLDS,
    corr_thresholds=PROXY_CORR_THRESHOLDS,
):
    dataset_key = tables["dataset_key"]
    dataset_name = tables["dataset_display"]
    universe = pd.Index(tables["protein_universe"]).astype(str)

    proxy_df = source_df[
        source_df["dataset_key"].eq(dataset_key)
        & source_df["stage"].isin(PROXY_STAGES)
    ].copy()

    if proxy_df.empty:
        print(f"No proxy/baseline rows for {dataset_name}.")
        return pd.DataFrame(), None

    corr_col = find_proxy_corr_column(proxy_df)
    if corr_col is None:
        print(
            f"No proxy-correlation metadata column found for {dataset_name}. "
            "Skipping proxy-correlation threshold inspection."
        )
        return pd.DataFrame(), None

    proxy_df[corr_col] = pd.to_numeric(proxy_df[corr_col], errors="coerce")
    proxy_best = best_per_protein(proxy_df, "test_pearson_corr", higher=True)
    proxy_idx = _index_by_protein(proxy_best)

    rows = []
    for pred_thr in pred_thresholds:
        for corr_thr in corr_thresholds:
            for protein in universe:
                if protein not in proxy_idx.index:
                    stage = "No eligible proxy"
                    metric = np.nan
                    corr_value = np.nan
                    run_label = ""
                else:
                    row = proxy_idx.loc[protein]
                    corr_value = row.get(corr_col, np.nan)
                    corr_for_rule = abs(corr_value) if USE_ABS_PROXY_CORR_FOR_ELIGIBILITY else corr_value
                    metric = row.get("test_pearson_corr", np.nan)
                    run_label = row.get("run_label", "")

                    eligible = np.isfinite(corr_for_rule) and corr_for_rule >= corr_thr
                    good = np.isfinite(metric) and metric >= pred_thr

                    if eligible and good:
                        stage = "Eligible and good"
                    elif eligible and not good:
                        stage = "Eligible but not good"
                    else:
                        stage = "No eligible proxy"

                rows.append({
                    "dataset": dataset_name,
                    "protein": protein,
                    "threshold": corr_thr,
                    "prediction_threshold": pred_thr,
                    "stage": stage,
                    "metric": metric,
                    "proxy_corr": corr_value,
                    "proxy_corr_column": corr_col,
                    "run_label": run_label,
                })

    return pd.DataFrame(rows), corr_col

def read_csv_if_exists(path: Path, columns: list) -> pd.DataFrame:
    """Read a CSV if it exists; return an empty dataframe with columns otherwise."""
    if not path.exists():
        return pd.DataFrame(columns=columns)
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=columns)

def read_saved_predictions(parquet_path: Path, csv_path: Path) -> pd.DataFrame:
    """Load saved predictions from parquet if available, with CSV fallback."""
    if parquet_path.exists():
        try:
            return pd.read_parquet(parquet_path)
        except Exception as e:
            print(f"Could not read {parquet_path.name}; trying CSV fallback. Error: {e}")
    if csv_path.exists():
        try:
            return pd.read_csv(csv_path)
        except pd.errors.EmptyDataError:
            return pd.DataFrame(columns=PREDICTION_COLUMNS)
    return pd.DataFrame(columns=PREDICTION_COLUMNS)

def refresh_combined_results(source_folders=None):
    """Rebuild combined notebook-8 files from all per-folder results that currently exist."""
    if source_folders is None:
        source_folders = EVALUATION_ORDER

    prediction_parts = []
    metric_parts = []
    coverage_parts = []
    loaded_sources = []

    for source_folder in source_folders:
        paths = source_output_paths(source_folder)
        preds = read_saved_predictions(paths["predictions_parquet"], paths["predictions_csv"])
        metrics_path = paths["metrics_protein_csv"]
        coverage_path = paths["feature_coverage_csv"]

        if preds.empty and not metrics_path.exists():
            continue

        if not preds.empty:
            prediction_parts.append(preds)
        if metrics_path.exists():
            metric_parts.append(read_csv_if_exists(metrics_path, METRIC_COLUMNS))
        if coverage_path.exists():
            coverage_parts.append(read_csv_if_exists(coverage_path, COVERAGE_COLUMNS))
        loaded_sources.append(source_folder)

    predictions = pd.concat(prediction_parts, ignore_index=True) if prediction_parts else pd.DataFrame(columns=PREDICTION_COLUMNS)
    metrics = pd.concat(metric_parts, ignore_index=True) if metric_parts else pd.DataFrame(columns=METRIC_COLUMNS)
    coverage = pd.concat(coverage_parts, ignore_index=True) if coverage_parts else pd.DataFrame(columns=COVERAGE_COLUMNS)

    print("\nRefreshing combined outputs from completed sources:", loaded_sources)
    write_predictions_with_fallback(predictions, PREDICTIONS_PARQUET, PREDICTIONS_CSV)
    metrics.to_csv(METRICS_PROTEIN_CSV, index=False)
    coverage.to_csv(FEATURE_COVERAGE_CSV, index=False)
    print("Saved:", METRICS_PROTEIN_CSV)
    print("Saved:", FEATURE_COVERAGE_CSV)

    return predictions, metrics, coverage

def require_columns(df, columns, df_name="dataframe"):
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"{df_name} is missing columns: {missing}\nAvailable: {list(df.columns)}")

def rmse(y_true, y_pred): return np.sqrt(mean_squared_error(y_true, y_pred))

def safe_metric(row, metric):
    try:
        return row.get(metric, np.nan)
    except Exception:
        return np.nan

def save_list(gene_set, filename):
    filepath = os.path.join(output_dir, filename)
    with open(filepath, 'w') as f:
        for gene in sorted(list(gene_set)):
            f.write(f"{gene}\n")
    print(f"Saved {len(gene_set)} genes to {filepath}")

def save_list_to_txt(items, filepath):
    with open(filepath, 'w') as f:
        for item in sorted(list(items)):
            f.write(f"{item}\n")
    print(f"Saved {len(items)} features to {filepath}")

def source_output_paths(source_folder: str) -> dict:
    """Return all per-folder output paths for a model source."""
    stem = str(source_folder)
    return {
        "predictions_parquet": OUTPUT_DIR / f"{stem}__eval_predictions_long.parquet",
        "predictions_csv": OUTPUT_DIR / f"{stem}__eval_predictions_long.csv",
        "metrics_protein_csv": OUTPUT_DIR / f"{stem}__eval_metrics_by_protein.csv",
        "feature_coverage_csv": OUTPUT_DIR / f"{stem}__feature_coverage.csv",
    }

def spearman_correlation(y_true, y_pred): 
     if len(y_true) < 2: return 0.0
     val, _ = spearmanr(y_true, y_pred)
     return val

def standardize_to_ccle_format(df_input):
    df = df_input.copy()
    # Force numeric
    df = df.apply(pd.to_numeric, errors='coerce')
    
    # Map ENSG -> Symbol using MyGene
    
    if mygene is None:
        # Without mygene, still coerce numeric values and collapse duplicate names.
        pass
    mg = mygene.MyGeneInfo() if mygene is not None else None
    ensembl_ids = df.columns.tolist()
    
    if mg is not None and any(str(col).startswith('ENSG') for col in ensembl_ids[:5]):
        try:
            results = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human', as_dataframe=True)
            mapping_dict = results['symbol'].dropna().to_dict()
            df = df.rename(columns=mapping_dict)
        except Exception:
            pass
    
    # Handle duplicates
    df = df.T.groupby(level=0).mean().T
    
    # Standardize Index
    df.index = df.index.astype(str).str.upper().str.strip()
    return df

def summarize_assignments(assignments, order=None):
    if assignments is None or len(assignments) == 0:
        return pd.DataFrame()

    group_cols = ["dataset", "threshold", "stage"]
    counts = (
        assignments
        .groupby(group_cols, dropna=False)["protein"]
        .nunique()
        .reset_index(name="proteins")
    )

    totals = (
        assignments
        .groupby(["dataset", "threshold"])["protein"]
        .nunique()
        .rename("total_proteins")
        .reset_index()
    )

    counts = counts.merge(totals, on=["dataset", "threshold"], how="left")
    counts["percent"] = 100 * counts["proteins"] / counts["total_proteins"]

    if order is not None:
        counts["stage"] = pd.Categorical(counts["stage"], categories=order, ordered=True)
        counts = counts.sort_values(["dataset", "threshold", "stage"]).reset_index(drop=True)

    return counts

def write_predictions_with_fallback(df: pd.DataFrame, parquet_path: Path, csv_path: Path) -> Path:
    """Save predictions as parquet when possible, otherwise CSV."""
    if df is None:
        df = pd.DataFrame()
    try:
        df.to_parquet(parquet_path, index=False)
        if csv_path.exists():
            # Avoid confusion from stale CSV fallbacks when parquet succeeds.
            csv_path.unlink()
        print("Saved:", parquet_path)
        return parquet_path
    except Exception as e:
        df.to_csv(csv_path, index=False)
        print("Parquet failed; saved CSV:", csv_path, "error:", e)
        return csv_path
