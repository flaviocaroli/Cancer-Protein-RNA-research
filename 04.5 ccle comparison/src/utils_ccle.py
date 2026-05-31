"""
Shared CCLE validation utilities.

This module contains the long/reused helpers extracted from:
  1. comparison_models.ipynb
  2. final_comparison_ccle.ipynb
  3. procedural_ccle.ipynb

The notebooks call set_module_globals(globals()) so utility functions that were
originally written as notebook functions can still access notebook variables
such as X, Y, test_samples, all_results, thresholds, and output folders.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple
import functools
import hashlib
import inspect
import json
import lzma
import pickle
import re
import warnings

import h5py
import joblib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display

warnings.filterwarnings("ignore")


def set_module_globals(namespace: dict | None = None) -> None:
    """Copy notebook globals into this module for moved notebook helpers."""
    if namespace is None:
        return
    protected = {
        "inspect", "functools", "set_module_globals", "sync_utils_globals",
        "_sync_from_caller", "_auto_sync", "_install_auto_sync", "__all__",
    }
    g = globals()
    for key, value in namespace.items():
        if key.startswith("__") or key in protected:
            continue
        g[key] = value


def sync_utils_globals(namespace: dict | None = None) -> None:
    """Alias for set_module_globals, kept for notebook readability."""
    set_module_globals(namespace)


_LAST_SYNC_SIGNATURE = None

def _sync_from_caller() -> None:
    """Sync only when a function is called directly from a notebook/user frame."""
    global _LAST_SYNC_SIGNATURE
    frame = inspect.currentframe()
    if frame is None:
        return
    caller = frame.f_back.f_back if frame.f_back is not None and frame.f_back.f_back is not None else None
    if caller is None:
        return
    caller_name = str(caller.f_globals.get("__name__", ""))
    # Jupyter/IPython user code normally has __name__ == "__main__".
    # Avoid copying large third-party module namespaces when utilities are called
    # indirectly by pandas/matplotlib internals. Cache per calling line so loops
    # over thousands of H5 files do not repeatedly copy the same namespace.
    if caller_name == "__main__":
        signature = (id(caller.f_globals), caller.f_code.co_filename, caller.f_lineno)
        if signature != _LAST_SYNC_SIGNATURE:
            set_module_globals(caller.f_globals)
            _LAST_SYNC_SIGNATURE = signature


def _auto_sync(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        _sync_from_caller()
        return func(*args, **kwargs)
    return wrapper


def _add_panel_label(ax, label):
    ax.text(
        0.01, 0.99, label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
        ha="left",
        bbox=dict(
            boxstyle="round,pad=0.2",
            facecolor="white",
            edgecolor="black",
            alpha=0.9,
        ),
    )


def _build_summary_text_from_metric_df(df, title_label="Per-Protein", count_label="Proteins"):
    n_items = len(df)

    summary_text = (
        f"{title_label} Summary\n"
        f"{'=' * 28}\n"
        f"{count_label} (n={n_items})\n\n"
        f"{_metric_summary_line(df.get('pearson', pd.Series(dtype=float)), 'Pearson')}\n"
        f"{_metric_summary_line(df.get('mse', pd.Series(dtype=float)), 'MSE')}\n"
        f"{_metric_summary_line(df.get('cosine', pd.Series(dtype=float)), 'Cosine Similarity')}\n"
        f"{_metric_summary_line(df.get('r2', pd.Series(dtype=float)), 'R²')}"
    )

    return summary_text


def _decode_h5_string_array(x):
    arr = np.asarray(x)

    if arr.shape == ():
        arr = arr.reshape(1)

    out = []
    for v in arr.ravel():
        if isinstance(v, bytes):
            out.append(v.decode("utf-8"))
        elif isinstance(v, np.bytes_):
            out.append(bytes(v).decode("utf-8"))
        else:
            out.append(str(v))

    return pd.Index(out)


def _decode_json_dataset(ds):
    """
    Decode JSON stored as:
      - uint8 array
      - bytes scalar
      - string scalar
      - array of byte strings
    """
    data = ds[()]

    if isinstance(data, bytes):
        raw = data
    elif isinstance(data, str):
        raw = data.encode("utf-8")
    else:
        arr = np.asarray(data)
        if arr.dtype.kind in {"S", "U", "O"}:
            raw = b"".join(
                x if isinstance(x, bytes) else str(x).encode("utf-8")
                for x in arr.ravel()
            )
        else:
            raw = bytes(arr.astype(np.uint8).ravel())

    return json.loads(raw.decode("utf-8"))


def _decode_json_from_h5_dataset(ds):
    data = ds[()]

    if isinstance(data, bytes):
        raw = data
    elif isinstance(data, str):
        raw = data.encode("utf-8")
    else:
        arr = np.asarray(data)

        if arr.dtype.kind in {"S", "U", "O"}:
            raw = b"".join(
                x if isinstance(x, bytes) else str(x).encode("utf-8")
                for x in arr.ravel()
            )
        else:
            raw = bytes(arr.astype(np.uint8).ravel())

    return json.loads(raw.decode("utf-8"))


def _estimate_r2_for_row(row):
    """
    Estimate R2 from stored MSE and the true Y values.

    R2 = 1 - SS_res / SS_tot
    If MSE is mean squared error on the test set, SS_res = MSE * n.

    This is only an estimate because it assumes the stored MSE used the same
    finite test values we recover here.
    """
    if "Y" not in globals():
        return np.nan

    protein = str(row.get("protein", ""))
    mse_val = row.get("mse", row.get("test_mse", np.nan))

    if protein not in Y.columns:
        return np.nan

    if not np.isfinite(mse_val):
        return np.nan

    h5_path = str(row.get("h5_path", ""))

    samples = None

    if h5_path and h5_path.lower() != "nan":
        if h5_path not in _TEST_SAMPLE_CACHE:
            _TEST_SAMPLE_CACHE[h5_path] = _read_h5_string_list(h5_path)
        samples = _TEST_SAMPLE_CACHE[h5_path]

    # Baselines in this notebook use the global test_samples variable.
    if samples is None and "test_samples" in globals():
        samples = pd.Index(test_samples)

    if samples is None:
        return np.nan

    samples = pd.Index(samples).intersection(Y.index)

    if len(samples) <= 1:
        return np.nan

    y_true = pd.to_numeric(Y.loc[samples, protein], errors="coerce").to_numpy(dtype=float)
    y_true = y_true[np.isfinite(y_true)]

    if len(y_true) <= 1:
        return np.nan

    ss_tot = float(np.sum((y_true - np.mean(y_true)) ** 2))
    if ss_tot <= 0:
        return np.nan

    ss_res = float(mse_val) * len(y_true)
    return 1.0 - ss_res / ss_tot


def _find_key_recursive(obj, target_keys):
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k in target_keys:
                return v

        for v in obj.values():
            found = _find_key_recursive(v, target_keys)
            if found is not None:
                return found

    elif isinstance(obj, list):
        for v in obj:
            found = _find_key_recursive(v, target_keys)
            if found is not None:
                return found

    return None


def _fmt_float(x, digits=4):
    if pd.isna(x):
        return "NA"
    try:
        return f"{float(x):.{digits}g}"
    except Exception:
        return str(x)


def _get_attr_float(group, name, default=np.nan):
    try:
        if group is not None and name in group.attrs:
            return _safe_float(group.attrs.get(name), default)
    except Exception:
        pass
    return default


def _get_global(name):
    return globals().get(name, None)


def _get_group_colors(group="proteins"):
    if group == "proteins":
        return {
            "pearson": "#F5B041",
            "mse": "#E74C3C",
            "cosine": "#5DADE2",
            "r2": "#AF7AC5",
        }
    elif group == "samples":
        return {
            "pearson": "#F8C471",
            "mse": "#EC7063",
            "cosine": "#85C1E9",
            "r2": "#C39BD3",
        }
    else:
        raise ValueError("group must be either 'samples' or 'proteins'")


def _get_h5_attr_float(h5, name, default=np.nan):
    try:
        if name in h5.attrs:
            return _safe_float(h5.attrs.get(name), default)
    except Exception:
        pass
    return default


def _index_by_protein(df):
    if df is None or len(df) == 0:
        return pd.DataFrame().rename_axis("protein")
    out = df.copy()
    out["protein"] = out["protein"].astype(str)
    return out.drop_duplicates("protein", keep="first").set_index("protein")


def _is_enet_label(label: str) -> bool:
    label = str(label).lower()
    return "enet" in label


def _is_xgb_label(label: str) -> bool:
    label = str(label).lower()
    return "xgb" in label


def _json_from_any_h5_dataset(ds):
    """
    Decode a JSON dataset stored either as uint8 bytes, bytes string, or scalar bytes.
    """
    data = ds[()]

    if isinstance(data, bytes):
        raw = data
    elif isinstance(data, str):
        raw = data.encode("utf-8")
    else:
        arr = np.asarray(data)
        if arr.dtype.kind in {"S", "U", "O"}:
            # array of byte strings / strings
            raw = b"".join(
                x if isinstance(x, bytes) else str(x).encode("utf-8")
                for x in arr.ravel()
            )
        else:
            # uint8 array
            raw = bytes(arr.astype(np.uint8).ravel())

    return json.loads(raw.decode("utf-8"))


def _metric_summary_line(series, name):
    s = pd.Series(series, dtype="float64").replace([np.inf, -np.inf], np.nan).dropna()

    if len(s) == 0:
        return (
            f"{name}:\n"
            f"  Mean: not available\n"
            f"  Min:  not available\n"
            f"  Max:  not available\n"
        )

    return (
        f"{name}:\n"
        f"  Mean: {s.mean():.4f} +/- {s.std():.4f}\n"
        f"  Min:  {s.min():.4f}\n"
        f"  Max:  {s.max():.4f}\n"
    )


def _plot_metric_hist(ax, values, title, xlabel, color, panel_label=None, bins=30):
    values = pd.Series(values, dtype="float64").replace([np.inf, -np.inf], np.nan).dropna()

    ax.set_facecolor("#FFFFFF")

    if len(values) > 0:
        ax.hist(values, bins=bins, color=color, edgecolor="black", alpha=0.8)

        mean_val = values.mean()
        median_val = values.median()

        if np.isfinite(mean_val):
            ax.axvline(mean_val, color="red", linestyle="--", linewidth=2, label=f"Mean = {mean_val:.3f}")

        if np.isfinite(median_val):
            ax.axvline(median_val, color="black", linestyle=":", linewidth=2, label=f"Median = {median_val:.3f}")

        ax.legend(fontsize=8, frameon=False)
    else:
        ax.text(
            0.5, 0.5,
            "No valid values",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=11,
        )

    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel("Count", fontsize=10)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(axis="y", alpha=0.25)

    for spine in ax.spines.values():
        spine.set_color("#BDC3C7")
        spine.set_linewidth(1.0)

    if panel_label is not None:
        _add_panel_label(ax, panel_label)


def _proxy_label(row):
    proxy_type = str(row.get("proxy_type", "")).lower()
    run_label = str(row.get("run_label", "")).lower()
    if "spearman" in proxy_type or "spearman" in run_label:
        return "Proxy Spearman"
    return "Proxy Pearson"


def _read_best_params(h5, summary):
    """
    Read best params from summary if present, otherwise model JSON datasets.
    """
    for key in ["best_params", "params", "model_params"]:
        if key in summary:
            return summary.get(key)

    for key in ["model/best_params_json", "model/params_json"]:
        if key in h5:
            try:
                return _decode_json_dataset(h5[key])
            except Exception:
                pass

    return {}


def _read_h5_string_list(path, candidates=None):
    """
    Tries to read test sample names from the H5 files.
    Used only if DETAILED_ESTIMATE_R2=True.
    """
    if h5py is None:
        return None

    if candidates is None:
        candidates = [
            "samples/test_samples",
            "test_samples",
            "data/test_samples",
            "split/test_samples",
        ]

    path = Path(path)
    if not path.exists():
        return None

    try:
        with h5py.File(path, "r") as h5:
            for key in candidates:
                if key not in h5:
                    continue

                arr = h5[key][()]
                arr = np.asarray(arr).ravel()

                out = []
                for x in arr:
                    if isinstance(x, bytes):
                        out.append(x.decode("utf-8"))
                    else:
                        out.append(str(x))

                return pd.Index(out)
    except Exception:
        return None

    return None


def _safe_float(x, default=np.nan):
    try:
        if x is None:
            return default
        return float(x)
    except Exception:
        return default


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


def _sanitize_filename(text, max_len=160):
    text = str(text)
    text = re.sub(r"[^\w\-\.]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text[:max_len] if len(text) > max_len else text


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


def add_family_group(df):
    df = df.copy()
    labels = df.get("run_label", pd.Series(index=df.index, dtype=object)).astype(str).str.lower()
    fam = df.get("model_family", pd.Series(index=df.index, dtype=object)).astype(str).str.lower()

    out = np.full(len(df), "other", dtype=object)
    out[(fam.str.contains("proxy") | labels.str.contains("proxy") | labels.str.contains("same_gene")).to_numpy()] = "proxy"
    out[(fam.str.contains("enet") | labels.str.contains("enet")).to_numpy()] = "enet"
    out[(fam.str.contains("xgb") | labels.str.contains("xgb")).to_numpy()] = "xgb"
    out[(fam.str.contains("sctranslator") | labels.str.contains("sctranslator")).to_numpy()] = "sctranslator"

    df["family_group"] = out
    return df


def as_dataframe_or_none(obj):
    if isinstance(obj, pd.DataFrame):
        df = obj.copy()
        df.index = df.index.map(str)
        df.columns = df.columns.map(str)
        return df
    return None


def best_by_metric(df, metric, higher=True):
    tmp = df.dropna(subset=[metric]).copy()
    if tmp.empty:
        return tmp
    tmp = tmp.sort_values(metric, ascending=not higher)
    return tmp.drop_duplicates("protein", keep="first").reset_index(drop=True)


def best_family(metric, family, higher=True):
    return best_per_protein(all_results[all_results["family_group"].eq(family)], metric=metric, higher=higher)


def best_per_protein(df, metric, higher=True):
    if df is None or len(df) == 0 or metric not in df.columns:
        return pd.DataFrame()
    use = df.dropna(subset=["protein", metric]).copy()
    if use.empty:
        return use
    use["protein"] = use["protein"].astype(str)
    use = use.sort_values(metric, ascending=not higher)
    return use.drop_duplicates("protein", keep="first").reset_index(drop=True)


def boxplot_metric(df, metric, title, higher_is_better=True):
    data = [df.loc[df["run_label"] == r, metric].dropna().values for r in RUNS_ORDER]
    labels = RUNS_ORDER
    fig_h = max(4, 0.35 * len(labels))
    fig, ax = plt.subplots(figsize=(12, fig_h))
    bp = ax.boxplot(data, vert=False, patch_artist=True, labels=labels, showfliers=False)
    for patch in bp["boxes"]:
        patch.set_facecolor("#4C78A8")
        patch.set_alpha(0.65)
    ax.set_title(title, fontsize=14, weight="bold")
    ax.set_xlabel(metric)
    ax.axvline(0, lw=1, alpha=0.35)
    ax.text(0.995, 1.02, "higher is better" if higher_is_better else "lower is better", ha="right", va="bottom", transform=ax.transAxes, alpha=0.7)
    plt.tight_layout()
    return fig, ax


def build_same_gene_baseline_from_lookup(corr_lookup, baseline_name, corr_threshold):
    rows = []

    possible = sorted(set(Y.columns.astype(str)).intersection(set(X.columns.astype(str))))

    candidate_count = 0
    eligible_count = 0

    print(f"\nBuilding {baseline_name}")
    print("- Same-name protein/gene pairs possible:", len(possible))
    print("- Correlation lookup entries:", len(corr_lookup))
    print("- Test samples used:", len(test_samples))
    print("- Train samples used:", len(train_samples))

    for protein in possible:
        c = corr_lookup.get(protein, np.nan)

        if not np.isfinite(c):
            continue

        candidate_count += 1

        if c < corr_threshold:
            continue

        eligible_count += 1

        y_true = Y.loc[test_samples, protein].to_numpy(dtype=float)
        x_test = X.loc[test_samples, protein].to_numpy(dtype=float)

        if SCALE_GENE_BASELINE_WITH_TRAIN_LINEAR_MODEL:
            x_train = X.loc[train_samples, protein].to_numpy(dtype=float)
            y_train = Y.loc[train_samples, protein].to_numpy(dtype=float)
            y_pred = linear_scale_prediction(x_train, y_train, x_test)
            variant = "scaled_train_linear"
        else:
            y_pred = x_test
            variant = "raw_gene_value"

        metrics = compute_metrics(y_true, y_pred)

        rows.append(
            {
                "protein": protein,
                "run_label": baseline_name,
                "model_family": "same_gene_baseline",
                "source_folder": baseline_name,
                "normalization": variant,
                "h5_path": "",
                "same_gene_corr_full": float(c),
                "test_mse": metrics["test_mse"],
                "test_pearson_corr": metrics["test_pearson_corr"],
                "test_cosine_similarity": metrics["test_cosine_similarity"],
                "n_test_used": metrics.get("n_test_used", np.nan),
                "selected_feature_count": 1,
                "train_n": len(train_samples),
                "test_n": len(test_samples),
                "best_params": {
                    "baseline": "same_name_gene_value",
                    "full_corr_threshold": corr_threshold,
                    "scaled": SCALE_GENE_BASELINE_WITH_TRAIN_LINEAR_MODEL,
                    "test_split_source": "read_from_xgb_enet_h5",
                },
            }
        )

    return pd.DataFrame(rows), candidate_count, eligible_count


def build_same_gene_proxy_all_from_lookup(corr_lookup, proxy_type):
    rows = []
    possible = sorted(set(Y.columns.astype(str)).intersection(set(X.columns.astype(str))))

    print(f"\nBuilding same-gene {proxy_type} proxy")
    print("- Same-name protein/gene pairs possible:", len(possible))
    print("- Correlation lookup entries:", len(corr_lookup))
    print("- Test samples used:", len(test_samples))
    print("- Train samples used:", len(train_samples))

    for protein in possible:
        c = corr_lookup.get(protein, np.nan)
        if not np.isfinite(c):
            continue

        y_true = Y.loc[test_samples, protein].to_numpy(dtype=float)
        x_test = X.loc[test_samples, protein].to_numpy(dtype=float)

        if SCALE_GENE_PROXY_WITH_TRAIN_LINEAR_MODEL:
            x_train = X.loc[train_samples, protein].to_numpy(dtype=float)
            y_train = Y.loc[train_samples, protein].to_numpy(dtype=float)
            y_pred = linear_scale_prediction(x_train, y_train, x_test)
            variant = "scaled_train_linear"
        else:
            y_pred = x_test
            variant = "raw_gene_value"

        metrics = compute_metrics(y_true, y_pred)

        rows.append({
            "protein": protein,
            "run_label": f"proxy_same_gene_{proxy_type}",
            "model_family": "same_gene_proxy",
            "family_group": "proxy",
            "source_folder": f"proxy_same_gene_{proxy_type}",
            "normalization": variant,
            "h5_path": "",
            "proxy_type": proxy_type,
            "same_gene_corr_full": float(c),
            "test_mse": metrics["test_mse"],
            "test_pearson_corr": metrics["test_pearson_corr"],
            "test_cosine_similarity": metrics["test_cosine_similarity"],
            "n_test_used": metrics.get("n_test_used", np.nan),
            "selected_feature_count": 1,
            "train_n": len(train_samples),
            "test_n": len(test_samples),
            "split_hash": model_split_hash,
            "best_params": {
                "baseline": "same_name_gene_value",
                "proxy_type": proxy_type,
                "scaled": SCALE_GENE_PROXY_WITH_TRAIN_LINEAR_MODEL,
                "test_split_source": "read_from_xgb_enet_h5",
            },
        })

    return pd.DataFrame(rows)


def build_same_name_corr_lookup_fast(corr_obj, possible_proteins):
    possible_proteins = [str(p) for p in possible_proteins]
    possible_set = set(possible_proteins)
    lookup = {}

    if isinstance(corr_obj, pd.DataFrame):
        df = corr_obj

        index_map = {str(v): i for i, v in enumerate(df.index)}
        column_map = {str(v): i for i, v in enumerate(df.columns)}

        common = possible_set.intersection(index_map).intersection(column_map)

        if len(common) > 0:
            for p in common:
                try:
                    lookup[p] = float(df.iat[index_map[p], column_map[p]])
                except Exception:
                    pass
            return lookup

        lower_cols = {str(c).lower(): c for c in df.columns}

        prot_candidates = ["protein", "protein_name", "target", "target_protein"]
        gene_candidates = ["gene", "feature", "rna", "transcript", "transcript_gene"]
        corr_candidates = [
            c for c in df.columns
            if "corr" in str(c).lower()
            or "pearson" in str(c).lower()
            or "spearman" in str(c).lower()
        ]

        prot_col = next((lower_cols[c] for c in prot_candidates if c in lower_cols), None)
        gene_col = next((lower_cols[c] for c in gene_candidates if c in lower_cols), None)
        corr_col = corr_candidates[0] if corr_candidates else None

        if prot_col and gene_col and corr_col:
            prot_vals = df[prot_col].astype(str)
            gene_vals = df[gene_col].astype(str)

            mask = (
                prot_vals.isin(possible_set)
                & gene_vals.isin(possible_set)
                & (prot_vals == gene_vals)
            )

            sub = df.loc[mask, [prot_col, corr_col]]

            for _, row in sub.iterrows():
                try:
                    lookup[str(row[prot_col])] = float(row[corr_col])
                except Exception:
                    pass

        return lookup

    if isinstance(corr_obj, dict):
        for p in possible_proteins:
            val = None

            if p in corr_obj:
                obj = corr_obj[p]

                if isinstance(obj, pd.Series):
                    idx_map = {str(v): i for i, v in enumerate(obj.index)}
                    if p in idx_map:
                        val = obj.iloc[idx_map[p]]

                elif isinstance(obj, dict):
                    if p in obj:
                        val = obj[p]

                else:
                    val = obj

            elif (p, p) in corr_obj:
                val = corr_obj[(p, p)]

            if val is not None:
                try:
                    lookup[p] = float(val)
                except Exception:
                    pass

        return lookup

    return lookup


def compute_metrics(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    mask = np.isfinite(y_true) & np.isfinite(y_pred)

    out = {
        "test_mse": np.nan,
        "test_pearson_corr": np.nan,
        "test_cosine_similarity": np.nan,
        "n_test_used": int(mask.sum()),
    }
    if mask.sum() == 0:
        return out

    yt = y_true[mask]
    yp = y_pred[mask]
    out["test_mse"] = float(np.mean((yt - yp) ** 2))

    if len(yt) >= 2 and np.std(yt) > 0 and np.std(yp) > 0:
        out["test_pearson_corr"] = float(np.corrcoef(yt, yp)[0, 1])

    denom = float(np.linalg.norm(yt) * np.linalg.norm(yp))
    if denom > 0:
        out["test_cosine_similarity"] = float(np.dot(yt, yp) / denom)

    return out


def decode_bytes_array(arr):
    out = []
    for x in arr:
        if isinstance(x, bytes):
            out.append(x.decode("utf-8"))
        else:
            out.append(str(x))
    return out


def get_attr_float(group_or_h5, key):
    try:
        if group_or_h5 is not None and key in group_or_h5.attrs:
            return safe_float(group_or_h5.attrs.get(key))
    except Exception:
        pass
    return np.nan


def get_attr_int(group_or_h5, key):
    try:
        if group_or_h5 is not None and key in group_or_h5.attrs:
            return safe_int(group_or_h5.attrs.get(key))
    except Exception:
        pass
    return np.nan


def get_same_name_corr(corr_obj, protein: str):
    """Robustly extract corr(protein, same-name gene) from several common formats."""
    protein = str(protein)

    # Matrix DataFrame: rows=proteins, columns=genes or vice versa.
    df = as_dataframe_or_none(corr_obj)
    if df is not None:
        if protein in df.index and protein in df.columns:
            try:
                return float(df.loc[protein, protein])
            except Exception:
                pass

        # Long-format table: columns like protein, gene, correlation.
        lower_cols = {c.lower(): c for c in df.columns}
        prot_candidates = ["protein", "protein_name", "target", "target_protein"]
        gene_candidates = ["gene", "feature", "rna", "transcript", "transcript_gene"]
        corr_candidates = [
            c for c in df.columns
            if "corr" in c.lower() or "pearson" in c.lower() or "spearman" in c.lower()
        ]

        prot_col = next((lower_cols[c] for c in prot_candidates if c in lower_cols), None)
        gene_col = next((lower_cols[c] for c in gene_candidates if c in lower_cols), None)
        corr_col = corr_candidates[0] if corr_candidates else None

        if prot_col and gene_col and corr_col:
            sub = df[
                (df[prot_col].astype(str) == protein) &
                (df[gene_col].astype(str) == protein)
            ]
            if len(sub):
                return float(sub.iloc[0][corr_col])

    # Dict: protein -> Series/dict containing same-name gene.
    if isinstance(corr_obj, dict):
        if protein in corr_obj:
            val = corr_obj[protein]
            if isinstance(val, pd.Series):
                val.index = val.index.map(str)
                if protein in val.index:
                    return float(val.loc[protein])
            if isinstance(val, dict) and protein in val:
                return float(val[protein])
        # dict with tuple key
        if (protein, protein) in corr_obj:
            return float(corr_obj[(protein, protein)])

    return np.nan


def is_baseline_label(label: str) -> bool:
    label = str(label).lower()
    return "same_gene" in label or "baseline" in label


def is_enet_label(label: str) -> bool:
    return "enet" in str(label).lower()


def is_xgb_label(label: str) -> bool:
    return "xgb" in str(label).lower()


def linear_scale_prediction(x_train, y_train, x_test):
    x_train = np.asarray(x_train, dtype=float)
    y_train = np.asarray(y_train, dtype=float)
    x_test = np.asarray(x_test, dtype=float)

    mask = np.isfinite(x_train) & np.isfinite(y_train)

    if mask.sum() < 2 or np.std(x_train[mask]) == 0:
        return x_test

    a, b = np.polyfit(x_train[mask], y_train[mask], deg=1)
    return a * x_test + b


def load_artifact_from_h5(path: Path):
    """
    Optional helper. Do not use during full scan unless needed,
    because loading thousands of pickled models is slow.
    """
    path = Path(path)
    with h5py.File(path, "r") as h5:
        if "model/artifact_pickle_xz" not in h5:
            return None
        payload = h5["model/artifact_pickle_xz"][:]

    return pickle.loads(lzma.decompress(np.asarray(payload, dtype=np.uint8).tobytes()))


def load_h5_summary(path: Path, source_label: str):
    """
    Robust loader for both XGB and ENet H5 formats.

    Important fix:
    XGB 80/20 files often store metrics as /metrics attrs or /metrics/summary_json,
    not only /summary/json or /metrics/json.
    """
    path = Path(path)
    info = parse_filename_info(path)

    row = {
        "protein": info["protein"],
        "run_label": None,
        "model_family": info["model_family"],
        "source_folder": source_label,
        "normalization": info["normalization"],
        "grid_preset": info["grid_preset"],
        "select_k": info["select_k"],
        "h5_path": str(path),
        "file": path.name,
        "file_size_mb": path.stat().st_size / (1024 ** 2),
        "test_mse": np.nan,
        "test_pearson_corr": np.nan,
        "test_cosine_similarity": np.nan,
        "selected_feature_count": np.nan,
        "train_n": np.nan,
        "test_n": np.nan,
        "best_params": {},
        "load_error": None,
    }

    try:
        with h5py.File(path, "r") as h5:
            root_attrs = dict(h5.attrs)
            metrics_group = h5.get("metrics", None)
            model_group = h5.get("model", None)
            samples_group = h5.get("samples", None)

            summary = read_json_dataset(
                h5,
                [
                    "summary/json",
                    "metrics/json",
                    "metrics/summary_json",
                    "metrics/summary",
                ],
            )

            # Identifiers
            row["protein"] = str(summary.get("protein") or root_attrs.get("protein") or row["protein"])
            row["normalization"] = str(summary.get("normalization") or root_attrs.get("normalization") or row["normalization"])
            row["model_family"] = str(summary.get("model_family") or row["model_family"])
            row["grid_preset"] = str(
                summary.get("grid_preset")
                or summary.get("grid")
                or root_attrs.get("grid_preset", row["grid_preset"])
                or row["grid_preset"]
            )
            row["select_k"] = safe_int(
                summary.get("select_k")
                or root_attrs.get("select_k", row["select_k"])
                or row["select_k"]
            )

            # Metrics: summary first, then /metrics attrs, then root attrs.
            row["test_mse"] = safe_float(summary.get("test_mse", summary.get("mse", np.nan)))
            if not np.isfinite(row["test_mse"]):
                row["test_mse"] = get_attr_float(metrics_group, "test_mse")
            if not np.isfinite(row["test_mse"]):
                row["test_mse"] = get_attr_float(metrics_group, "mse")
            if not np.isfinite(row["test_mse"]):
                row["test_mse"] = get_attr_float(h5, "test_mse")

            row["test_pearson_corr"] = safe_float(
                summary.get("test_pearson_corr", summary.get("pearson_corr", np.nan))
            )
            if not np.isfinite(row["test_pearson_corr"]):
                row["test_pearson_corr"] = get_attr_float(metrics_group, "test_pearson_corr")
            if not np.isfinite(row["test_pearson_corr"]):
                row["test_pearson_corr"] = get_attr_float(metrics_group, "pearson_corr")
            if not np.isfinite(row["test_pearson_corr"]):
                row["test_pearson_corr"] = get_attr_float(h5, "test_pearson_corr")

            row["test_cosine_similarity"] = safe_float(
                summary.get(
                    "test_cosine_similarity",
                    summary.get(
                        "cosine_similarity",
                        summary.get("vector_cosine_similarity_all_loo", np.nan),
                    ),
                )
            )
            if not np.isfinite(row["test_cosine_similarity"]):
                row["test_cosine_similarity"] = get_attr_float(metrics_group, "test_cosine_similarity")
            if not np.isfinite(row["test_cosine_similarity"]):
                row["test_cosine_similarity"] = get_attr_float(metrics_group, "cosine_similarity")
            if not np.isfinite(row["test_cosine_similarity"]):
                row["test_cosine_similarity"] = get_attr_float(metrics_group, "vector_cosine_similarity_all_loo")
            if not np.isfinite(row["test_cosine_similarity"]):
                row["test_cosine_similarity"] = get_attr_float(h5, "test_cosine_similarity")

            # Train/test sample counts
            train_samples = summary.get("train_samples", [])
            test_samples = summary.get("test_samples", [])
            if isinstance(train_samples, list):
                row["train_n"] = len(train_samples)
            if isinstance(test_samples, list):
                row["test_n"] = len(test_samples)

            if samples_group is not None:
                if "train_samples" in samples_group:
                    row["train_n"] = len(samples_group["train_samples"])
                if "test_samples" in samples_group:
                    row["test_n"] = len(samples_group["test_samples"])

            # Selected features
            selected_features = summary.get("selected_features", None)
            if isinstance(selected_features, list):
                row["selected_feature_count"] = len(selected_features)

            if model_group is not None:
                if "selected_feature_count" in model_group.attrs:
                    row["selected_feature_count"] = safe_int(model_group.attrs.get("selected_feature_count"))
                elif "selected_features" in model_group:
                    row["selected_feature_count"] = len(model_group["selected_features"])

            row["best_params"] = _read_best_params(h5, summary)

    except Exception as e:
        row["load_error"] = f"{type(e).__name__}: {e}"

    # Clean model family from source label if filename was ambiguous
    if row["model_family"] == "unknown":
        if "xgb" in str(source_label).lower():
            row["model_family"] = "xgb"
        elif "enet" in str(source_label).lower():
            row["model_family"] = "enet"

    # Keep labels concise but distinct.
    label_bits = [str(source_label), str(row["normalization"])]
    if row["grid_preset"] not in {"", "unknown", "None", "nan"}:
        label_bits.append(str(row["grid_preset"]))
    if pd.notna(row["select_k"]):
        label_bits.append(f"k{int(row['select_k'])}")
    row["run_label"] = " | ".join(label_bits)

    return row


def make_per_protein_metric_df(
    results_df,
    run_label=None,
    model_family=None,
    source_folder=None,
    estimate_r2=False,
    sort_by="pearson",
):
    """
    Converts the notebook's all_results table into the dataframe expected
    by the detailed plot.

    Required all_results columns:
      - protein
      - test_pearson_corr
      - test_mse
      - test_cosine_similarity
      - run_label
    """
    df = results_df.copy()

    if run_label is not None:
        if isinstance(run_label, (list, tuple, set, pd.Index, np.ndarray)):
            df = df[df["run_label"].isin(list(run_label))].copy()
        else:
            df = df[df["run_label"].astype(str).eq(str(run_label))].copy()

    if model_family is not None and "model_family" in df.columns:
        if isinstance(model_family, (list, tuple, set, pd.Index, np.ndarray)):
            df = df[df["model_family"].isin(list(model_family))].copy()
        else:
            df = df[df["model_family"].astype(str).eq(str(model_family))].copy()

    if source_folder is not None and "source_folder" in df.columns:
        if isinstance(source_folder, (list, tuple, set, pd.Index, np.ndarray)):
            df = df[df["source_folder"].isin(list(source_folder))].copy()
        else:
            df = df[df["source_folder"].astype(str).eq(str(source_folder))].copy()

    if df.empty:
        return pd.DataFrame(columns=["protein", "pearson", "mse", "cosine", "r2", "run_label"])

    out = pd.DataFrame({
        "protein": df["protein"].astype(str).values,
        "pearson": pd.to_numeric(df.get("test_pearson_corr", np.nan), errors="coerce").values,
        "mse": pd.to_numeric(df.get("test_mse", np.nan), errors="coerce").values,
        "cosine": pd.to_numeric(df.get("test_cosine_similarity", np.nan), errors="coerce").values,
        "run_label": df.get("run_label", pd.Series(["unknown"] * len(df))).astype(str).values,
    })

    # Keep useful metadata if present.
    for col in ["model_family", "source_folder", "normalization", "grid_preset", "select_k", "h5_path", "file"]:
        if col in df.columns:
            out[col] = df[col].values

    # R2 is not directly stored in the notebook's all_results.
    if "test_r2" in df.columns:
        out["r2"] = pd.to_numeric(df["test_r2"], errors="coerce").values
    elif "r2" in df.columns:
        out["r2"] = pd.to_numeric(df["r2"], errors="coerce").values
    else:
        out["r2"] = np.nan

    if estimate_r2:
        out["r2"] = out.apply(_estimate_r2_for_row, axis=1)

    # If a filtered set accidentally contains duplicate proteins,
    # keep the better one by the requested metric.
    if sort_by == "mse":
        out = (
            out.sort_values("mse", ascending=True, na_position="last")
               .drop_duplicates("protein", keep="first")
               .reset_index(drop=True)
        )
    elif sort_by == "cosine":
        out = (
            out.sort_values("cosine", ascending=False, na_position="last")
               .drop_duplicates("protein", keep="first")
               .reset_index(drop=True)
        )
    else:
        out = (
            out.sort_values("pearson", ascending=False, na_position="last")
               .drop_duplicates("protein", keep="first")
               .reset_index(drop=True)
        )

    return out


def make_run_label(row):
    bits = [str(row.get("source_folder", "unknown"))]

    for col in ["normalization", "grid_preset"]:
        v = row.get(col)
        if pd.notna(v) and str(v) not in ["", "None", "nan"]:
            bits.append(str(v))

    k = row.get("select_k")
    if pd.notna(k):
        bits.append(f"k{int(k)}")

    return " | ".join(bits)


def metric_summary_table(df, group_col="run_label"):
    if df.empty:
        return pd.DataFrame()
    return (
        df.groupby(group_col)
        .agg(
            proteins=("protein", "nunique"),
            pearson_median=("test_pearson_corr", "median"),
            pearson_mean=("test_pearson_corr", "mean"),
            cosine_median=("test_cosine_similarity", "median"),
            mse_median=("test_mse", "median"),
            mse_mean=("test_mse", "mean"),
        )
        .sort_values(["pearson_median", "cosine_median"], ascending=False)
    )


def mse_procedural_improvement_assignments():
    proxy_idx = _index_by_protein(proxy_best_mse)
    enet_idx = _index_by_protein(enet_best_mse)
    xgb_idx = _index_by_protein(xgb_best_mse)
    rows = []

    for protein in protein_universe.astype(str):
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
            "protein": protein,
            "threshold": 0.0,
            "stage": stage,
            "best_mse": current_best,
            "proxy_mse": proxy_mse,
            "enet_mse": enet_mse,
            "xgb_mse": xgb_mse,
        })

    return pd.DataFrame(rows)


def mse_winner_assignments():
    proxy_idx = _index_by_protein(proxy_best_mse)
    enet_idx = _index_by_protein(enet_best_mse)
    xgb_idx = _index_by_protein(xgb_best_mse)
    rows = []

    for protein in protein_universe.astype(str):
        candidates = []
        for fam, idx in [("proxy", proxy_idx), ("enet", enet_idx), ("xgb", xgb_idx)]:
            if protein in idx.index:
                row = idx.loc[protein]
                mse = row.get("test_mse", np.nan)
                if np.isfinite(mse):
                    candidates.append((fam, mse, row))

        if not candidates:
            stage = "No valid MSE"
            value = np.nan
            run_label = ""
        else:
            fam, value, row = min(candidates, key=lambda x: x[1])
            if fam == "proxy":
                base = _proxy_label(row).replace("Proxy", "Proxy")
                stage = f"{base} lowest MSE"
            elif fam == "enet":
                stage = "ENet lowest MSE"
            else:
                stage = "XGB lowest MSE"
            run_label = row.get("run_label", "")

        rows.append({
            "protein": protein,
            "threshold": 0.0,
            "stage": stage,
            "test_mse": value,
            "run_label": run_label,
        })

    return pd.DataFrame(rows)


def parse_filename(path):
    """
    Robust parser for names like:
      A2M_xgb_8020_zscore_fast_k2000.h5
      A2M_enet_8020_zscore_fast_k5000.h5
      A2M_xgb_8020_qbin50_tiny_k50.h5
    """
    stem = Path(path).stem

    known_norms = {
        "none", "zscore", "rank", "ranking", "binning",
        "qbin5", "qbin10", "qbin20", "qbin50", "qbin100",
    }

    tokens = stem.split("_")

    model_family = None
    if "_xgb_" in f"_{stem}_":
        model_family = "xgb"
    elif "_enet_" in f"_{stem}_":
        model_family = "enet"

    normalization = None
    grid_preset = None
    select_k = None

    norm_idx = None
    for i, tok in enumerate(tokens):
        if tok in known_norms or re.fullmatch(r"qbin\d+", tok):
            normalization = tok
            norm_idx = i
            break

    for tok in tokens:
        if tok in {"tiny", "fast", "current"}:
            grid_preset = tok
            break

    for tok in tokens:
        m = re.fullmatch(r"k(\d+)", tok)
        if m:
            select_k = int(m.group(1))
            break

    # Protein is everything before xgb/enet if present; otherwise before normalization.
    if "xgb" in tokens:
        model_idx = tokens.index("xgb")
        protein = "_".join(tokens[:model_idx])
    elif "enet" in tokens:
        model_idx = tokens.index("enet")
        protein = "_".join(tokens[:model_idx])
    elif norm_idx is not None:
        protein = "_".join(tokens[:norm_idx])
    else:
        protein = tokens[0] if tokens else stem

    return protein, model_family, normalization, grid_preset, select_k


def parse_filename_info(path: Path):
    """
    Parse filenames such as:
      A2M_xgb_8020_zscore_fast_k2000.h5
      A2M_xgb_8020_qbin50_tiny_k50.h5
      A2M_enet_8020_zscore_fast_k5000.h5
    """
    path = Path(path)
    stem = path.stem
    tokens = stem.split("_")

    model_family = "unknown"
    if "_xgb_8020_" in stem:
        model_family = "xgb"
    elif "_enet_8020_" in stem:
        model_family = "enet"

    protein = protein_from_filename(path)

    normalization = "unknown"
    grid_preset = "unknown"
    select_k = np.nan

    known_norms = {
        "none", "zscore", "rank", "ranking", "binning",
        "qbin5", "qbin10", "qbin20", "qbin50", "qbin100",
    }

    for tok in tokens:
        if tok in known_norms or re.fullmatch(r"qbin\d+", tok):
            normalization = tok
            break

    for tok in tokens:
        if tok in {"tiny", "fast", "current"}:
            grid_preset = tok
            break

    for tok in tokens:
        m = re.fullmatch(r"k(\d+)", tok)
        if m:
            select_k = int(m.group(1))
            break

    return {
        "protein": protein,
        "model_family": model_family,
        "normalization": normalization,
        "grid_preset": grid_preset,
        "select_k": select_k,
    }


def plot_detailed(
    per_sample_df=None,
    per_protein_df=None,
    *legacy_args,
    group="proteins",
    out_path=None,
    title=None,
    show=True,
):
    """
    Detailed plot function compatible with the style of the function shared earlier.

    In this notebook, the natural use is:
        plot_detailed(per_protein_df=my_df, group="proteins")

    Expected per_protein_df columns:
        protein, pearson, mse, cosine, r2

    R2 may be all NaN; the R2 panel will show "No valid values".
    """
    if group not in {"samples", "proteins"}:
        raise ValueError("group must be either 'samples' or 'proteins'")

    colors = _get_group_colors(group)

    if group == "samples":
        if per_sample_df is None:
            raise ValueError("per_sample_df is required when group='samples'")
        df = per_sample_df.copy()
        item_col = "sample"
        title_prefix = "Per-Sample"
        count_label = "Samples"
    else:
        if per_protein_df is None:
            raise ValueError("per_protein_df is required when group='proteins'")
        df = per_protein_df.copy()
        item_col = "protein"
        title_prefix = "Per-Protein"
        count_label = "Proteins"

    for needed_col in [item_col, "pearson", "mse", "cosine"]:
        if needed_col not in df.columns:
            raise ValueError(f"Missing required column: {needed_col}")

    if "r2" not in df.columns:
        df["r2"] = np.nan

    for col in ["pearson", "mse", "cosine", "r2"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    fig = plt.figure(figsize=(16, 10))
    fig.patch.set_facecolor("#F8F9FA")
    gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.32)

    if title is not None:
        fig.suptitle(title, fontsize=16, fontweight="bold", y=0.98)

    ax1 = fig.add_subplot(gs[0, 0])
    _plot_metric_hist(
        ax1,
        df["pearson"],
        f"{title_prefix} PCC",
        "Pearson Correlation",
        colors["pearson"],
        panel_label="a)",
    )

    ax2 = fig.add_subplot(gs[0, 1])
    _plot_metric_hist(
        ax2,
        df["mse"],
        f"{title_prefix} MSE",
        "MSE",
        colors["mse"],
        panel_label="b)",
    )

    ax3 = fig.add_subplot(gs[0, 2])
    _plot_metric_hist(
        ax3,
        df["cosine"],
        f"{title_prefix} Cosine Similarity",
        "Cosine Similarity",
        colors["cosine"],
        panel_label="c)",
    )

    ax4 = fig.add_subplot(gs[1, 0])
    _plot_metric_hist(
        ax4,
        df["r2"],
        f"{title_prefix} R²",
        "R² Score",
        colors["r2"],
        panel_label="d)",
    )

    ax5 = fig.add_subplot(gs[1, 1])

    valid_rank_df = df.dropna(subset=["pearson"]).copy()

    if len(valid_rank_df) > 0:
        best_10 = (
            valid_rank_df.sort_values("pearson", ascending=False, na_position="last")
            .head(10)
            .sort_values("pearson", ascending=True)
        )

        worst_10 = (
            valid_rank_df.sort_values("pearson", ascending=True, na_position="last")
            .head(10)
            .sort_values("pearson", ascending=False)
        )

        combined = pd.concat([best_10, worst_10], axis=0).reset_index(drop=True)
        bar_colors = (["#58D68D"] * len(best_10)) + (["#E74C3C"] * len(worst_10))
        y_pos = np.arange(len(combined))

        ax5.barh(
            y_pos,
            combined["pearson"].values,
            color=bar_colors,
            edgecolor="black",
            alpha=0.8,
        )

        ax5.set_yticks(y_pos)
        ax5.set_yticklabels(combined[item_col].astype(str).values, fontsize=8)
        ax5.axvline(0, color="red", linestyle="-", linewidth=1)
        ax5.set_xlabel("PCC", fontsize=10)
        ax5.set_title(f"Best 10 / Worst 10 {count_label} by PCC", fontsize=12, fontweight="bold")
        ax5.grid(axis="x", alpha=0.25)
    else:
        ax5.text(
            0.5, 0.5,
            "No valid Pearson values",
            ha="center",
            va="center",
            transform=ax5.transAxes,
        )
        ax5.set_title(f"Best 10 / Worst 10 {count_label} by PCC", fontsize=12, fontweight="bold")

    _add_panel_label(ax5, "e)")

    ax6 = fig.add_subplot(gs[1, 2])
    ax6.axis("off")
    ax6.set_title("Summary Statistics", fontsize=12, fontweight="bold")
    _add_panel_label(ax6, "f)")

    summary_text = _build_summary_text_from_metric_df(
        df=df,
        title_label=title_prefix,
        count_label=count_label,
    )

    ax6.text(
        0.5,
        0.92,
        summary_text,
        transform=ax6.transAxes,
        fontsize=9,
        va="top",
        ha="center",
        multialignment="left",
        fontfamily="monospace",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.45),
    )

    if title is not None:
        plt.tight_layout(rect=[0, 0, 1, 0.96])
    else:
        plt.tight_layout()

    if out_path is not None:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_path, dpi=DETAILED_DPI, bbox_inches="tight", facecolor=fig.get_facecolor())
        print(f"Saved: {out_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig


def plot_detailed_for_all_run_labels(
    results_df=None,
    min_proteins=None,
    include_baselines=None,
    estimate_r2=None,
    save=None,
    show=None,
    out_dir=None,
):
    """
    Plot every run_label/model present in all_results.
    """
    if results_df is None:
        results_df = all_results

    if min_proteins is None:
        min_proteins = DETAILED_MIN_PROTEINS_PER_RUN

    if include_baselines is None:
        include_baselines = DETAILED_INCLUDE_BASELINES

    counts = (
        results_df.groupby("run_label")["protein"]
        .nunique()
        .sort_values(ascending=False)
    )

    labels = counts[counts >= int(min_proteins)].index.tolist()

    if not include_baselines:
        labels = [
            lab for lab in labels
            if "baseline" not in str(lab).lower() and "same_gene" not in str(lab).lower()
        ]

    print(f"Plotting {len(labels)} run labels")
    print("-" * 80)
    for lab in labels:
        print(f"{lab}  |  proteins={int(counts.loc[lab])}")

    outputs = {}

    for lab in labels:
        print("\n" + "=" * 100)
        print(f"Plotting: {lab}")
        outputs[lab] = plot_detailed_for_run_label(
            lab,
            results_df=results_df,
            estimate_r2=estimate_r2,
            save=save,
            show=show,
            out_dir=out_dir,
        )

    return outputs


def plot_detailed_for_model_family(
    model_family,
    results_df=None,
    estimate_r2=None,
    save=None,
    show=None,
    out_dir=None,
):
    """
    Optional helper: combine all rows for one model_family.
    Example:
        plot_detailed_for_model_family("xgb")
        plot_detailed_for_model_family("enet")
        plot_detailed_for_model_family("same_gene_baseline")
    """
    if results_df is None:
        results_df = all_results

    if estimate_r2 is None:
        estimate_r2 = DETAILED_ESTIMATE_R2

    if save is None:
        save = DETAILED_SAVE_FIGURES

    if show is None:
        show = DETAILED_SHOW_FIGURES

    if out_dir is None:
        out_dir = DETAILED_PLOT_DIR

    per_protein_df = make_per_protein_metric_df(
        results_df,
        model_family=model_family,
        estimate_r2=estimate_r2,
        sort_by="pearson",
    )

    if per_protein_df.empty:
        print(f"No rows found for model_family={model_family!r}")
        return None

    safe_name = _sanitize_filename(f"family_{model_family}")
    out_path = Path(out_dir) / f"detailed_{safe_name}.png" if save else None

    title = f"Detailed per-protein metrics: model_family = {model_family}"

    fig = plot_detailed(
        per_protein_df=per_protein_df,
        group="proteins",
        out_path=out_path,
        title=title,
        show=show,
    )

    return per_protein_df


def plot_detailed_for_run_label(
    run_label,
    results_df=None,
    estimate_r2=None,
    save=None,
    show=None,
    out_dir=None,
):
    """
    Plot one model/run from all_results.
    """
    if results_df is None:
        results_df = all_results

    if estimate_r2 is None:
        estimate_r2 = DETAILED_ESTIMATE_R2

    if save is None:
        save = DETAILED_SAVE_FIGURES

    if show is None:
        show = DETAILED_SHOW_FIGURES

    if out_dir is None:
        out_dir = DETAILED_PLOT_DIR

    per_protein_df = make_per_protein_metric_df(
        results_df,
        run_label=run_label,
        estimate_r2=estimate_r2,
        sort_by="pearson",
    )

    if per_protein_df.empty:
        print(f"No rows found for run_label={run_label!r}")
        return None

    safe_name = _sanitize_filename(run_label)
    out_path = Path(out_dir) / f"detailed_{safe_name}.png" if save else None

    title = f"Detailed per-protein metrics: {run_label}"

    fig = plot_detailed(
        per_protein_df=per_protein_df,
        group="proteins",
        out_path=out_path,
        title=title,
        show=show,
    )

    return per_protein_df


def plot_metric_distributions(df):
    if df.empty:
        print("No data to plot.")
        return

    plot_df = df.copy()
    # Keep run labels with enough rows visible; group tiny labels for readability.
    counts = plot_df.groupby("run_label")["protein"].nunique()
    keep = counts[counts >= MIN_PROTEINS_TO_SHOW_AS_OWN_SEGMENT].index
    plot_df["plot_label"] = np.where(plot_df["run_label"].isin(keep), plot_df["run_label"], "Other small runs")

    for metric, title, xlabel in [
        ("test_pearson_corr", "Distribution of test Pearson by run", "test Pearson"),
        ("test_cosine_similarity", "Distribution of test cosine similarity by run", "test cosine similarity"),
        ("test_mse", "Distribution of test MSE by run", "test MSE"),
    ]:
        use = plot_df.dropna(subset=[metric]).copy()
        if use.empty:
            print(f"No valid values for {metric}")
            continue

        labels = list(use.groupby("plot_label")[metric].median().sort_values(ascending=(metric == "test_mse")).index)
        data = [use.loc[use["plot_label"] == lab, metric].values for lab in labels]

        fig, ax = plt.subplots(figsize=(12, max(4, 0.45 * len(labels))))
        bp = ax.boxplot(data, vert=False, labels=labels, patch_artist=True, showfliers=False)
        for patch in bp["boxes"]:
            patch.set_alpha(0.72)
        ax.set_title(title, fontsize=14, weight="bold")
        ax.set_xlabel(xlabel)
        ax.grid(axis="x", alpha=0.18)
        ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout()
        plt.show()


def plot_single_partition(assignments, title, order, save_name=None):
    tmp = assignments.copy()
    tmp["threshold"] = 0.0
    summary = summarize_assignments(tmp, order=order)
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
        summary[["stage", "proteins", "percent"]]
        .sort_values("stage")
        .reset_index(drop=True)
    )
    return summary


def plot_stacked_partition_rows(summary, title, order, xlabel="Proteins", save_name=None, figsize=None, row_label_func=None):
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
            ax.barh(y, width, left=left, color=color, edgecolor="white", linewidth=1.2, label=label)
            used_labels.add(stage)

            if pct >= 6:
                ax.text(left + width / 2, y, f"{width}\n{pct:.0f}%", ha="center", va="center", color="white", fontsize=9, weight="bold")
            elif width > 0 and total > 0:
                ax.text(left + width / 2, y + 0.27, f"{width}", ha="center", va="bottom", fontsize=8, color=color)

            left += width

    ax.set_yticks(y_positions)
    if row_label_func is None:
        row_labels = [f"Pearson ≥ {thr:.2f}" for thr in thresholds]
    else:
        row_labels = [row_label_func(thr) for thr in thresholds]
    ax.set_yticklabels(row_labels)
    ax.set_xlim(0, max_total)
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=14, weight="bold", pad=12)
    ax.grid(axis="x", alpha=0.18)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.invert_yaxis()
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3, frameon=False, fontsize=9)
    plt.tight_layout()

    if save_name and SAVE_OUTPUTS:
        fig.savefig(FIGURE_DIR / save_name, dpi=180, bbox_inches="tight")
    plt.show()
    return fig


def print_family_summary(best_df, metric, threshold, higher=True):
    df = best_df.dropna(subset=[metric]).copy()
    if df.empty:
        print(f"No valid values for {metric}")
        return

    df["is_good"] = df[metric] >= threshold if higher else df[metric] <= threshold
    good = df[df["is_good"]]
    bad = df[~df["is_good"]]

    xgb_good = int(good["run_label"].map(is_xgb_label).sum())
    enet_good = int(good["run_label"].map(is_enet_label).sum())
    baseline_good = int(good["run_label"].map(is_baseline_label).sum())
    other_good = len(good) - xgb_good - enet_good - baseline_good

    print(f"\n{metric} family summary")
    print("-" * 80)
    print(f"Total proteins with valid metric: {len(df)}")
    print(f"Good proteins: {len(good)} ({100 * len(good) / len(df):.1f}%)")
    print(f"Bad proteins:  {len(bad)} ({100 * len(bad) / len(df):.1f}%)")
    print(f"Good proteins best predicted by XGB:                {xgb_good}")
    print(f"Good proteins best predicted by ENet:               {enet_good}")
    print(f"Good proteins best predicted by same-gene baseline: {baseline_good}")
    if other_good > 0:
        print(f"Good proteins best predicted by other labels:       {other_good}")


def print_model_family_summary(best_df: pd.DataFrame, metric: str, threshold: float, higher: bool = True):
    """
    Print concise family-level counts:
      - how many proteins are good overall
      - how many good proteins are best predicted by XGB
      - how many good proteins are best predicted by ENet
    """
    df = best_df.dropna(subset=[metric]).copy()

    if df.empty:
        print(f"No valid rows for {metric}")
        return

    if higher:
        df["is_good"] = df[metric] >= threshold
    else:
        df["is_good"] = df[metric] <= threshold

    good = df[df["is_good"]].copy()
    bad = df[~df["is_good"]].copy()

    xgb_good = good["run_label"].map(_is_xgb_label).sum()
    enet_good = good["run_label"].map(_is_enet_label).sum()
    other_good = len(good) - int(xgb_good) - int(enet_good)

    print(f"\n{metric} family summary")
    print("-" * 70)
    print(f"Total proteins with valid metric: {len(df)}")
    print(f"Good proteins: {len(good)} ({100 * len(good) / len(df):.1f}%)")
    print(f"Bad proteins:  {len(bad)} ({100 * len(bad) / len(df):.1f}%)")
    print(f"Good proteins best predicted by XGB:  {int(xgb_good)}")
    print(f"Good proteins best predicted by ENet: {int(enet_good)}")
    if other_good > 0:
        print(f"Good proteins best predicted by other labels: {int(other_good)}")


def procedural_mse_cascade(threshold, universe=None):
    universe = pd.Index(protein_universe if universe is None else universe).astype(str)
    remaining = set(universe)
    rows = []

    proxy_idx = _index_by_protein(proxy_best_mse)
    enet_idx = _index_by_protein(enet_best_mse)
    xgb_idx = _index_by_protein(xgb_best_mse)

    # 1) Proxy first: best same-gene proxy by MSE.
    for protein in list(remaining):
        if protein in proxy_idx.index:
            row = proxy_idx.loc[protein]
            value = row.get("test_mse", np.nan)
            if np.isfinite(value) and value <= threshold:
                rows.append({
                    "protein": protein,
                    "threshold": threshold,
                    "stage": _proxy_label(row),
                    "metric": value,
                    "run_label": row.get("run_label", ""),
                })
                remaining.remove(protein)

    # 2) ENet rescues only remaining proteins.
    for protein in list(remaining):
        if protein in enet_idx.index:
            row = enet_idx.loc[protein]
            value = row.get("test_mse", np.nan)
            if np.isfinite(value) and value <= threshold:
                rows.append({
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
            value = row.get("test_mse", np.nan)
            if np.isfinite(value) and value <= threshold:
                rows.append({
                    "protein": protein,
                    "threshold": threshold,
                    "stage": "XGB",
                    "metric": value,
                    "run_label": row.get("run_label", ""),
                })
                remaining.remove(protein)

    # 4) Anything not reaching the MSE threshold.
    for protein in sorted(remaining):
        rows.append({
            "protein": protein,
            "threshold": threshold,
            "stage": "Still not good",
            "metric": np.nan,
            "run_label": "",
        })

    return pd.DataFrame(rows)


def procedural_pearson_cascade(threshold, include_sc=False, sctranslator_best=None, universe=None):
    universe = pd.Index(protein_universe if universe is None else universe).astype(str)
    remaining = set(universe)
    rows = []

    proxy_idx = _index_by_protein(proxy_best_pearson)
    enet_idx = _index_by_protein(enet_best_pearson)
    xgb_idx = _index_by_protein(xgb_best_pearson)
    sc_idx = _index_by_protein(sctranslator_best) if sctranslator_best is not None else pd.DataFrame().rename_axis("protein")

    # 1) Proxy first: best of Pearson and Spearman same-gene proxies.
    for protein in list(remaining):
        if protein in proxy_idx.index:
            row = proxy_idx.loc[protein]
            value = row.get("test_pearson_corr", np.nan)
            if np.isfinite(value) and value >= threshold:
                rows.append({
                    "protein": protein,
                    "threshold": threshold,
                    "stage": _proxy_label(row),
                    "metric": value,
                    "run_label": row.get("run_label", ""),
                })
                remaining.remove(protein)

    # 2) ENet rescues only remaining proteins.
    for protein in list(remaining):
        if protein in enet_idx.index:
            row = enet_idx.loc[protein]
            value = row.get("test_pearson_corr", np.nan)
            if np.isfinite(value) and value >= threshold:
                rows.append({
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
            value = row.get("test_pearson_corr", np.nan)
            if np.isfinite(value) and value >= threshold:
                rows.append({
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
                value = row.get("test_pearson_corr", np.nan)
                if np.isfinite(value) and value >= threshold:
                    rows.append({
                        "protein": protein,
                        "threshold": threshold,
                        "stage": f"scTranslator {row.get('sc_size', row.get('run_label', '')).replace('scTranslator ', '')}",
                        "metric": value,
                        "run_label": row.get("run_label", ""),
                    })
                    remaining.remove(protein)

    final_label = "Still not good after scTranslator" if include_sc else "Still not good"
    for protein in sorted(remaining):
        rows.append({
            "protein": protein,
            "threshold": threshold,
            "stage": final_label,
            "metric": np.nan,
            "run_label": "",
        })

    out = pd.DataFrame(rows)
    return out


def protein_from_filename(path: Path):
    name = path.stem
    for token in ["_xgb_8020_", "_enet_8020_"]:
        if token in name:
            return name.split(token)[0]
    return name.split("_")[0]


def read_h5_summary(path, folder_label):
    path = Path(path)
    protein_file, family_file, norm_file, grid_file, k_file = parse_filename(path)

    row = {
        "source_folder": folder_label,
        "h5_path": str(path),
        "file": path.name,
        "file_size_mb": path.stat().st_size / (1024 ** 2),
        "protein_from_file": protein_file,
        "model_family_from_file": family_file,
        "normalization_from_file": norm_file,
        "grid_preset_from_file": grid_file,
        "select_k_from_file": k_file,
    }

    try:
        with h5py.File(path, "r") as h5:
            attrs = dict(h5.attrs)

            # Support all known schemas.
            summary = read_json_dataset(
                h5,
                [
                    "summary/json",
                    "metrics/json",
                    "metrics/summary_json",
                    "metrics/summary",
                ],
            )

            metrics_group = h5.get("metrics", None)
            samples_group = h5.get("samples", None)
            model_group = h5.get("model", None)

            # Basic identifiers
            row["protein"] = (
                summary.get("protein")
                or attrs.get("protein")
                or protein_file
            )

            row["normalization"] = (
                summary.get("normalization")
                or attrs.get("normalization")
                or norm_file
            )

            row["model_family"] = (
                summary.get("model_family")
                or family_file
            )

            row["grid_preset"] = (
                summary.get("grid_preset")
                or summary.get("grid")
                or grid_file
            )

            row["select_k"] = (
                summary.get("select_k")
                or attrs.get("select_k")
                or k_file
            )

            # Metrics: summary first, then /metrics attrs, then root attrs.
            row["test_mse"] = _safe_float(
                summary.get("test_mse", summary.get("mse", np.nan))
            )
            if not np.isfinite(row["test_mse"]):
                row["test_mse"] = _get_attr_float(metrics_group, "test_mse", np.nan)
            if not np.isfinite(row["test_mse"]):
                row["test_mse"] = _get_attr_float(metrics_group, "mse", np.nan)
            if not np.isfinite(row["test_mse"]):
                row["test_mse"] = _get_h5_attr_float(h5, "test_mse", np.nan)

            row["test_pearson_corr"] = _safe_float(
                summary.get("test_pearson_corr", summary.get("pearson_corr", np.nan))
            )
            if not np.isfinite(row["test_pearson_corr"]):
                row["test_pearson_corr"] = _get_attr_float(metrics_group, "test_pearson_corr", np.nan)
            if not np.isfinite(row["test_pearson_corr"]):
                row["test_pearson_corr"] = _get_attr_float(metrics_group, "pearson_corr", np.nan)
            if not np.isfinite(row["test_pearson_corr"]):
                row["test_pearson_corr"] = _get_h5_attr_float(h5, "test_pearson_corr", np.nan)

            row["test_cosine_similarity"] = _safe_float(
                summary.get(
                    "test_cosine_similarity",
                    summary.get(
                        "vector_cosine_similarity_all_loo",
                        summary.get("cosine_similarity", np.nan),
                    ),
                )
            )
            if not np.isfinite(row["test_cosine_similarity"]):
                row["test_cosine_similarity"] = _get_attr_float(metrics_group, "test_cosine_similarity", np.nan)
            if not np.isfinite(row["test_cosine_similarity"]):
                row["test_cosine_similarity"] = _get_attr_float(metrics_group, "cosine_similarity", np.nan)
            if not np.isfinite(row["test_cosine_similarity"]):
                row["test_cosine_similarity"] = _get_h5_attr_float(h5, "test_cosine_similarity", np.nan)

            # Train/test sample counts
            train_samples = summary.get("train_samples", [])
            test_samples = summary.get("test_samples", [])

            if isinstance(train_samples, list):
                row["train_n"] = len(train_samples)
            else:
                row["train_n"] = np.nan

            if isinstance(test_samples, list):
                row["test_n"] = len(test_samples)
            else:
                row["test_n"] = np.nan

            if samples_group is not None:
                if "train_samples" in samples_group:
                    row["train_n"] = len(samples_group["train_samples"])
                if "test_samples" in samples_group:
                    row["test_n"] = len(samples_group["test_samples"])

            # Selected features
            selected_features = summary.get("selected_features", [])
            if isinstance(selected_features, list) and len(selected_features) > 0:
                row["selected_feature_count"] = len(selected_features)
            else:
                row["selected_feature_count"] = np.nan

            if model_group is not None:
                if "selected_feature_count" in model_group.attrs:
                    row["selected_feature_count"] = int(model_group.attrs["selected_feature_count"])
                elif "selected_features" in model_group:
                    row["selected_feature_count"] = len(model_group["selected_features"])

                row["artifact_kind"] = str(model_group.attrs.get("artifact_kind", ""))

            row["best_params"] = _read_best_params(h5, summary)

            row["h5_format"] = str(attrs.get("format", ""))
            row["load_error"] = None

    except Exception as e:
        row["load_error"] = f"{type(e).__name__}: {e}"

    return row


def read_json_dataset(h5, candidates):
    """
    Try all possible JSON paths used by different XGB/ENet H5 pipeline versions.
    """
    for key in candidates:
        if key in h5:
            try:
                return _decode_json_dataset(h5[key])
            except Exception:
                pass
    return {}


def read_test_samples_from_h5(path):
    path = str(path)

    direct_keys = [
        "samples/test_samples",
        "test_samples",
        "data/test_samples",
        "split/test_samples",
        "splits/test_samples",
        "sample_names/test_samples",
    ]

    json_keys = [
        "summary/json",
        "metrics/json",
        "metrics/summary_json",
        "metrics/summary",
    ]

    with h5py.File(path, "r") as h5:
        for key in direct_keys:
            if key in h5:
                return _decode_h5_string_array(h5[key][()])

        for key in json_keys:
            if key in h5:
                try:
                    obj = _decode_json_from_h5_dataset(h5[key])
                    found = _find_key_recursive(
                        obj,
                        target_keys={
                            "test_samples",
                            "test_sample_names",
                            "samples_test",
                        },
                    )

                    if found is not None:
                        return pd.Index([str(x) for x in found])

                except Exception:
                    pass

    return None


def require_file(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return path


def safe_float(x):
    try:
        if x is None:
            return np.nan
        v = float(x)
        return v if np.isfinite(v) else np.nan
    except Exception:
        return np.nan


def safe_int(x):
    try:
        if x is None:
            return np.nan
        return int(x)
    except Exception:
        return np.nan


def scan_all_results(result_folders):
    rows = []

    for label, folder in result_folders.items():
        folder = Path(folder)

        if not folder.exists():
            print(f"{label:10s}: missing -> {folder}")
            continue

        files = sorted(folder.glob("*.h5"))
        print(f"{label:10s}: {len(files):6d} H5 files")

        for p in files:
            rows.append(read_h5_summary(p, label))

    return pd.DataFrame(rows)


def split_hash(samples):
    samples = pd.Index(samples.astype(str))
    key = "\n".join(sorted(samples))
    return hashlib.sha1(key.encode("utf-8")).hexdigest()[:12]


def summarize_assignments(assignments, order=None):
    if assignments is None or len(assignments) == 0:
        return pd.DataFrame()
    counts = (
        assignments
        .groupby(["threshold", "stage"], dropna=False)["protein"]
        .nunique()
        .reset_index(name="proteins")
    )
    totals = assignments.groupby("threshold")["protein"].nunique().rename("total_proteins")
    counts = counts.merge(totals, on="threshold", how="left")
    counts["percent"] = 100 * counts["proteins"] / counts["total_proteins"]
    if order is not None:
        counts["stage"] = pd.Categorical(counts["stage"], categories=order, ordered=True)
        counts = counts.sort_values(["threshold", "stage"]).reset_index(drop=True)
    return counts


def summarize_best(best_df, metric, higher=True, threshold=None):
    out = (
        best_df.groupby("run_label")
        .agg(proteins=("protein", "nunique"), median_metric=(metric, "median"), mean_metric=(metric, "mean"))
        .reset_index()
        .sort_values("proteins", ascending=False)
    )
    if threshold is not None:
        good_mask = best_df[metric] >= threshold if higher else best_df[metric] <= threshold
        good_counts = best_df.loc[good_mask].groupby("run_label")["protein"].nunique().rename("good_proteins")
        out = out.merge(good_counts, on="run_label", how="left")
        out["good_proteins"] = out["good_proteins"].fillna(0).astype(int)
    return out

# -----------------------------------------------------------------------------
# Unified plotting/partition helpers. These replace near-duplicate notebook
# versions and are compatible with both comparison notebooks.
# -----------------------------------------------------------------------------

def make_partition_table(
    best_df: pd.DataFrame,
    metric: str,
    threshold: float,
    higher: bool = True,
    min_proteins_to_show: int = 200,
) -> Tuple[pd.DataFrame, pd.DataFrame, str]:
    """
    Partition best-per-protein rows into good model segments plus one bad segment.

    Small good segments are collapsed into "Other good models" for compact bars,
    while raw_part preserves the uncollapsed counts for inspection.
    """
    df = best_df.copy() if best_df is not None else pd.DataFrame()
    needed = {"protein", "run_label", metric}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in best_df: {missing}")

    df = df.dropna(subset=[metric]).copy()
    if df.empty:
        empty = pd.DataFrame(
            columns=["partition_label", "proteins", "median_metric", "is_bad", "is_other"]
        )
        return empty, empty, "Bad"

    if higher:
        df["is_good"] = df[metric] >= threshold
        bad_label = f"Bad / {metric} < {threshold:g}"
    else:
        df["is_good"] = df[metric] <= threshold
        bad_label = f"Bad / {metric} > {threshold:g}"

    df["partition_label_raw"] = np.where(df["is_good"], df["run_label"], bad_label)

    raw_part = (
        df.groupby("partition_label_raw", dropna=False)
        .agg(proteins=("protein", "nunique"), median_metric=(metric, "median"))
        .reset_index()
        .rename(columns={"partition_label_raw": "partition_label"})
    )
    raw_part["is_bad"] = raw_part["partition_label"].eq(bad_label)
    raw_part["is_other"] = False

    small_good_labels = raw_part.loc[
        (~raw_part["is_bad"]) & (raw_part["proteins"] < int(min_proteins_to_show)),
        "partition_label",
    ].tolist()

    df["partition_label"] = df["partition_label_raw"]
    if small_good_labels:
        df.loc[df["partition_label_raw"].isin(small_good_labels), "partition_label"] = (
            f"Other good models (<{min_proteins_to_show} proteins each)"
        )

    part = (
        df.groupby("partition_label", dropna=False)
        .agg(proteins=("protein", "nunique"), median_metric=(metric, "median"))
        .reset_index()
    )
    part["is_bad"] = part["partition_label"].eq(bad_label)
    part["is_other"] = part["partition_label"].astype(str).str.startswith("Other good models")
    part = part.sort_values(
        ["is_bad", "is_other", "proteins"], ascending=[True, True, False]
    ).reset_index(drop=True)

    raw_part = raw_part.sort_values(
        ["is_bad", "proteins"], ascending=[True, False]
    ).reset_index(drop=True)

    return part, raw_part, bad_label


def _family_counts_for_partition(best_df, metric, threshold, higher=True):
    df = best_df.dropna(subset=[metric]).copy()
    if df.empty:
        return df, df, {}
    df["is_good"] = df[metric] >= threshold if higher else df[metric] <= threshold
    good = df[df["is_good"]].copy()
    bad = df[~df["is_good"]].copy()
    counts = {
        "xgb": int(good["run_label"].map(is_xgb_label).sum()) if "run_label" in good else 0,
        "enet": int(good["run_label"].map(is_enet_label).sum()) if "run_label" in good else 0,
        "baseline": int(good["run_label"].map(is_baseline_label).sum()) if "run_label" in good else 0,
    }
    counts["other"] = len(good) - counts["xgb"] - counts["enet"] - counts["baseline"]
    return good, bad, counts


def plot_partition_bar(
    best_df: pd.DataFrame,
    metric: str,
    threshold: float,
    higher: bool = True,
    title: Optional[str] = None,
    filename: Optional[str] = None,
    min_proteins_to_show: int = 200,
) -> pd.DataFrame:
    """
    Compact horizontal stacked partition plot.

    Compatible with both original signatures:
    - compare_models used filename=...
    - final_comparison_ccle omitted filename
    """
    part, raw_part, bad_label = make_partition_table(
        best_df,
        metric,
        threshold,
        higher=higher,
        min_proteins_to_show=min_proteins_to_show,
    )

    total = int(part["proteins"].sum()) if len(part) else 0
    if total == 0:
        print(f"No valid data for {metric}")
        return part

    bad_n = int(part.loc[part["partition_label"].eq(bad_label), "proteins"].sum())
    good_n = total - bad_n

    good, bad, counts = _family_counts_for_partition(best_df, metric, threshold, higher=higher)
    print(f"\n{metric} family summary")
    print("-" * 80)
    print(f"Total proteins with valid metric: {total}")
    print(f"Good proteins: {len(good)} ({100 * len(good) / total:.1f}%)")
    print(f"Bad proteins:  {len(bad)} ({100 * len(bad) / total:.1f}%)")
    print(f"Good proteins best predicted by XGB:                {counts.get('xgb', 0)}")
    print(f"Good proteins best predicted by ENet:               {counts.get('enet', 0)}")
    if counts.get("baseline", 0) > 0:
        print(f"Good proteins best predicted by same-gene baseline: {counts.get('baseline', 0)}")
    if counts.get("other", 0) > 0:
        print(f"Good proteins best predicted by other labels:       {counts.get('other', 0)}")

    print("\nVisible bar segments:")
    display(part[["partition_label", "proteins", "median_metric"]])

    hidden = raw_part[(~raw_part["is_bad"]) & (raw_part["proteins"] < int(min_proteins_to_show))].copy()
    if len(hidden):
        print(f"\nHidden individual good segments with < {min_proteins_to_show} proteins:")
        display(hidden[["partition_label", "proteins", "median_metric"]])

    base_colors = [
        "#3B6EA8", "#4E9F8A", "#6BA84F", "#F2A65A", "#9D79BC",
        "#7A6F9B", "#D6A77A", "#8FB8DE", "#A1C181", "#B8B8B8",
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

    fig, ax = plt.subplots(figsize=(14, 3.2))
    left = 0
    handles, labels = [], []
    for _, row in part.iterrows():
        width = int(row["proteins"])
        label = row["partition_label"]
        pct = 100 * width / total if total else 0
        color = color_map[label]
        bar = ax.barh(
            y=[0], width=[width], left=left, height=0.48,
            color=color, edgecolor="white", linewidth=1.6,
        )
        if pct >= 5:
            ax.text(left + width / 2, 0, f"{width}\n{pct:.1f}%", ha="center", va="center", color="white", fontsize=10, weight="bold")
        elif width > 0:
            ax.text(left + width / 2, 0.38, f"{width}", ha="center", va="bottom", color=color, fontsize=9, weight="bold")

        med = row["median_metric"]
        label_text = f"{label} — {width} ({pct:.1f}%)"
        if np.isfinite(med):
            label_text += f", median={med:.3g}"
        labels.append(label_text)
        handles.append(bar[0])
        left += width

    ax.set_xlim(0, total)
    ax.set_yticks([])
    ax.set_xlabel(f"Proteins, total = {total}", fontsize=10)
    ax.set_title(title or f"Best model per protein by {metric}", fontsize=14, weight="bold", pad=14)
    ax.text(
        0.5, -0.42,
        f"Good: {good_n} ({100 * good_n / total:.1f}%)   |   Bad: {bad_n} ({100 * bad_n / total:.1f}%)",
        ha="center", va="center", transform=ax.transAxes, fontsize=10, alpha=0.85,
    )
    ax.grid(axis="x", alpha=0.15)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.72), ncol=1, frameon=False, fontsize=9)
    plt.tight_layout()

    if filename:
        out_dir = Path(globals().get("OUT_DIR", "."))
        out_dir.mkdir(parents=True, exist_ok=True)
        output_path = out_dir / filename
        fig.savefig(output_path, dpi=180, bbox_inches="tight")
        print(f"Saved: {output_path}")

    plt.show()
    return part


def proxy_corr_threshold_assignments(proxy_type, pred_thresholds=None, corr_thresholds=None):
    if pred_thresholds is None:
        pred_thresholds = globals().get("PEARSON_GOOD_THRESHOLDS", [0.30, 0.50, 0.70])
    if corr_thresholds is None:
        corr_thresholds = globals().get("PROXY_CORR_THRESHOLDS", [round(x, 2) for x in np.arange(0.00, 0.85, 0.10)])

    proxy_df = globals()["proxy_results"][globals()["proxy_results"]["proxy_type"].astype(str).eq(proxy_type)].copy()
    proxy_idx = _index_by_protein(proxy_df)
    protein_universe_local = globals()["protein_universe"]
    rows = []

    for pred_thr in pred_thresholds:
        for corr_thr in corr_thresholds:
            for protein in protein_universe_local.astype(str):
                label = "No eligible proxy"
                metric = np.nan
                corr_val = np.nan
                if protein in proxy_idx.index:
                    row = proxy_idx.loc[protein]
                    corr_val = row.get("same_gene_corr_full", np.nan)
                    metric = row.get("test_pearson_corr", np.nan)
                    if np.isfinite(corr_val) and corr_val >= corr_thr:
                        label = "Eligible and good" if np.isfinite(metric) and metric >= pred_thr else "Eligible but not good"
                rows.append({
                    "proxy_type": proxy_type,
                    "prediction_threshold": pred_thr,
                    "corr_threshold": corr_thr,
                    "threshold": corr_thr,
                    "stage": label,
                    "protein": protein,
                    "test_pearson_corr": metric,
                    "same_gene_corr_full": corr_val,
                })
    return pd.DataFrame(rows)


def plot_proxy_corr_grid(proxy_type, pred_threshold):
    assign = proxy_corr_threshold_assignments(
        proxy_type=proxy_type,
        pred_thresholds=[pred_threshold],
        corr_thresholds=globals().get("PROXY_CORR_THRESHOLDS"),
    )
    summary = summarize_assignments(assign, order=["Eligible and good", "Eligible but not good", "No eligible proxy"])
    if summary.empty:
        print("No proxy-correlation summary to plot.")
        return assign, summary

    order = ["Eligible and good", "Eligible but not good", "No eligible proxy"]
    thresholds = list(summary["threshold"].drop_duplicates())
    fig, ax = plt.subplots(figsize=(14, max(3.2, 0.45 * len(thresholds) + 1.8)))
    max_total = int(summary.groupby("threshold")["total_proteins"].first().max())
    used_labels = set()
    for y, thr in enumerate(thresholds):
        sub = summary[summary["threshold"].eq(thr)]
        total = int(sub["total_proteins"].iloc[0])
        left = 0
        for stage in order:
            row = sub[sub["stage"].astype(str).eq(stage)]
            width = int(row["proteins"].iloc[0]) if len(row) else 0
            if width <= 0:
                continue
            pct = 100 * width / total if total else 0
            color = globals().get("COLOR_MAP", {}).get(stage, None)
            label = stage if stage not in used_labels else None
            ax.barh(y, width, left=left, color=color, edgecolor="white", linewidth=1.2, label=label)
            used_labels.add(stage)
            if pct >= 6:
                ax.text(left + width / 2, y, f"{width}\n{pct:.0f}%", ha="center", va="center", color="white", fontsize=9, weight="bold")
            left += width
    ax.set_yticks(np.arange(len(thresholds)))
    ax.set_yticklabels([f"Proxy corr >= {thr:.2f}" for thr in thresholds])
    ax.set_xlim(0, max_total)
    ax.set_xlabel("Proteins")
    ax.set_title(f"{proxy_type.title()} same-gene proxy eligibility, Pearson >= {pred_threshold:.2f}", fontsize=14, weight="bold")
    ax.grid(axis="x", alpha=0.18)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.invert_yaxis()
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3, frameon=False, fontsize=9)
    plt.tight_layout()
    if globals().get("SAVE_OUTPUTS", False):
        fig_dir = Path(globals().get("FIGURE_DIR", "."))
        fig_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(fig_dir / f"proxy_corr_grid_{proxy_type}_pearson_{pred_threshold:.2f}.png", dpi=180, bbox_inches="tight")
    plt.show()
    return assign, summary


def load_sctranslator_prediction_tables(bundle_paths=None, gt=None):
    if bundle_paths is None:
        bundle_paths = globals().get("SCTRANSLATOR_BUNDLE_PATHS", {})
    if gt is None:
        gt = globals().get("GT_scaled", None)
    if gt is None:
        raise ValueError("GT_scaled/gt is required to align scTranslator predictions.")

    pred_dfs = {}
    missing = []
    for label, path in bundle_paths.items():
        path = Path(path)
        if not path.exists():
            missing.append(str(path))
            continue
        bundle = np.load(path, allow_pickle=True)
        required = {"predictions", "protein_names", "cell_names"}
        missing_keys = required.difference(set(bundle.files))
        if missing_keys:
            raise KeyError(f"{path.name} is missing keys: {sorted(missing_keys)}")
        df = pd.DataFrame(bundle["predictions"], columns=[str(x) for x in bundle["protein_names"]])
        df.index = [str(x) for x in bundle["cell_names"]]
        df = df.reindex(index=gt.index.astype(str), columns=gt.columns.astype(str))
        pred_dfs[label] = df

    if missing:
        print("Missing scTranslator files; those labels will be skipped:")
        for m in missing:
            print("-", m)
    return pred_dfs


def evaluate_sctranslator_predictions(pred_dfs, gt=None, eval_samples=None):
    if gt is None:
        gt = globals().get("GT_scaled", None)
    if gt is None:
        raise ValueError("GT_scaled/gt is required to evaluate scTranslator predictions.")
    if eval_samples is None:
        eval_samples = gt.index
    eval_samples = pd.Index(eval_samples.astype(str)).intersection(gt.index.astype(str))
    model_split_hash_local = globals().get("model_split_hash", "")

    rows = []
    for label, pred_df in pred_dfs.items():
        common_cols = gt.columns.astype(str).intersection(pred_df.columns.astype(str))
        for protein in common_cols:
            y_true = gt.loc[eval_samples, protein].to_numpy(dtype=float)
            y_pred = pred_df.loc[eval_samples, protein].to_numpy(dtype=float)
            metrics = compute_metrics(y_true, y_pred)
            rows.append({
                "protein": str(protein),
                "run_label": f"scTranslator {label}",
                "model_family": "sctranslator",
                "family_group": "sctranslator",
                "sc_size": str(label),
                "source_folder": "scTranslator npz",
                "test_mse": metrics["test_mse"],
                "test_pearson_corr": metrics["test_pearson_corr"],
                "test_cosine_similarity": metrics["test_cosine_similarity"],
                "n_test_used": metrics.get("n_test_used", np.nan),
                "train_n": np.nan,
                "test_n": len(eval_samples),
                "split_hash": model_split_hash_local,
                "h5_path": "",
            })
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Install lightweight auto-sync wrappers after all functions are defined.
# -----------------------------------------------------------------------------

def _install_auto_sync() -> None:
    skip = {
        "set_module_globals", "sync_utils_globals", "_sync_from_caller",
        "_auto_sync", "_install_auto_sync",
    }
    g = globals()
    for name, obj in list(g.items()):
        if name in skip:
            continue
        if inspect.isfunction(obj) and getattr(obj, "__module__", None) == __name__:
            g[name] = _auto_sync(obj)


_install_auto_sync()

__all__ = [
    name for name, obj in globals().items()
    if inspect.isfunction(obj) and getattr(obj, "__module__", None) == __name__
]
