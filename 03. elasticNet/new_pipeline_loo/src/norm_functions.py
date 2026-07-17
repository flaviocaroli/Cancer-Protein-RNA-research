# norm_functions_optimized.py
# Compatible with Python 3.8+
#
# Drop-in replacement for your norm_functions.py, but much faster by:
#   1) Computing BINNING/RANKING transforms ONCE per run_protein_group_norm call (not once per protein)
#   2) Removing the gigantic (n_bins, n_samples, n_genes) broadcast in binning
#   3) Removing DataFrame .T / .drop / .loc in hot loops (use NumPy indexing)
#   4) Avoiding sklearn SelectPercentile object creation in hot loops (direct f_regression + percentile threshold)
#   5) Avoiding pandas mean_squared_error overhead (use NumPy MSE)

import itertools
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

import numpy as np
import pandas as pd
from sklearn.feature_selection import f_regression
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import KFold

# Fix seed globally
np.random.seed(42)

# ============================================================
# RELATIVE REPRESENTATION HELPER FUNCTIONS
# ============================================================

def _apply_binning_matrix(X: pd.DataFrame, n_bins: int = 50) -> pd.DataFrame:
    """
    Fast binning normalization row-wise (per sample).
    Matches the original logic:
        x_transformed = (x <= nanquantile(x, qs, axis=1)[..., None]).sum(axis=0)

    Key idea:
      For a row i with thresholds t (ascending), the original returns:
        bin(x) = count_k [ x <= t_k ] = n_bins - count_k [ t_k < x ]
      which is:
        bin(x) = n_bins - searchsorted(t, x, side="left")

    Zeros are assigned to bin 0 (we convert zeros -> NaN so they become 0 after masking).
    """
    if n_bins <= 0:
        raise ValueError("n_bins must be >= 1")

    # Copy to avoid in-place modifications to caller's DataFrame
    x = X.to_numpy(dtype=np.float64, copy=True)
    n_samples, n_genes = x.shape

    # zeros -> NaN (ignored in nanquantile; later mapped to 0)
    x[x == 0] = np.nan

    qs = np.linspace(0.0, 1.0, int(n_bins), dtype=np.float64)

    # Compute quantiles once (shape: (n_bins, n_samples))
    q = np.nanquantile(x, qs, axis=1)
    qT = q.T  # (n_samples, n_bins)

    out = np.zeros((n_samples, n_genes), dtype=np.float64)

    # Per-row searchsorted (avoids creating the huge 3D boolean tensor)
    for i in range(n_samples):
        t = qT[i]
        if np.isnan(t).all():
            continue  # all zeros/NaNs in this row => stays 0

        row = x[i]
        m = ~np.isnan(row)
        if not m.any():
            continue

        # idx = number of thresholds < x  (because side="left")
        idx = np.searchsorted(t, row[m], side="left")
        out[i, m] = float(n_bins) - idx.astype(np.float64)

    return pd.DataFrame(out, index=X.index, columns=X.columns)


def _apply_ranking_matrix(X: pd.DataFrame) -> pd.DataFrame:
    """
    Ranking normalization row-wise (per sample).
    Uses argsort(argsort()).
    """
    x = X.to_numpy(dtype=np.float64, copy=False)
    ranked = x.argsort(axis=1).argsort(axis=1).astype(np.float64, copy=False)
    return pd.DataFrame(ranked, index=X.index, columns=X.columns)


# ============================================================
# Fast feature selection (percentile via f_regression)
# ============================================================

def _select_feature_indices_percentile(
    X_train: np.ndarray,   # (n_samples, n_features)
    y_train: np.ndarray,   # (n_samples,)
    percentile: int,
) -> np.ndarray:
    """
    Select features using ANOVA f_regression percentile.
    Mimics SelectPercentile behavior more closely by using a percentile threshold
    (ties may return > exact %), but is much faster than DataFrame transpose + sklearn selector objects.
    """
    p = int(percentile)
    if p <= 0:
        return np.empty(0, dtype=np.int64)
    if p >= 100:
        return np.arange(X_train.shape[1], dtype=np.int64)

    scores, _ = f_regression(X_train, y_train)
    # Replace NaNs (can happen with constant columns or constant y) with -inf so they won't be selected
    scores = np.asarray(scores, dtype=np.float64)
    np.nan_to_num(scores, copy=False, nan=-np.inf, posinf=-np.inf, neginf=-np.inf)

    # If everything is -inf (degenerate), fall back to one arbitrary feature
    if not np.isfinite(scores).any():
        return np.array([0], dtype=np.int64)

    thresh = np.percentile(scores, 100.0 - float(p))
    idx = np.flatnonzero(scores >= thresh).astype(np.int64, copy=False)

    # Safety: ensure non-empty
    if idx.size == 0:
        idx = np.array([int(np.argmax(scores))], dtype=np.int64)
    return idx


# ============================================================
# IO helpers
# ============================================================

def _safe_to_parquet_or_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        df.to_parquet(path, compression="snappy")
    except Exception:
        df.to_csv(path.with_suffix(".csv"), index=False)


# ============================================================
# Best-only summary (fast path, no need to build a big raw DF)
# ============================================================

def _best_only_row_from_tracking(
    protein_name: str,
    sample_ids: np.ndarray,
    best_i: int,
    best_abs_error: float,
    total_sse: float,
    best_sq_error: float,
    params_dict: Dict[str, Any],
    coef_dict: Dict[str, float],
    *,
    fit_intercept: bool,
    output_dir: Optional[Union[str, Path]] = None,
    normalization: str = "RelativeRepresentation",
) -> pd.DataFrame:
    n = int(sample_ids.shape[0])
    if n > 1:
        mse_all_others = (total_sse - best_sq_error) / float(n - 1)
    else:
        mse_all_others = np.nan

    out = pd.DataFrame(
        [{
            "protein": protein_name,
            "best_sample": str(sample_ids[best_i]),
            "best_abs_error": float(best_abs_error),
            "mse_all_others": float(mse_all_others) if np.isfinite(mse_all_others) else np.nan,
            "params_dict": params_dict,
            "feat_coeffs_dict": coef_dict,
            "normalization": normalization,
        }]
    )

    if output_dir is not None:
        out_path = Path(output_dir) / "{}_best.parquet".format(protein_name)
        _safe_to_parquet_or_csv(out, out_path)

    return out


# ============================================================
# Core optimized CV (expects X already transformed)
# ============================================================

def _cv_pipeline_fixed_percentile_pretransformed(
    param_combinations: List[Dict[str, Any]],
    X_transformed_df: pd.DataFrame,
    y: pd.DataFrame,
    protein_name: str,
    *,
    feature_selector: str = "anova",  # kept for API compat; only "anova" implemented here
    n_splits: int = 3,
    percentile: int = 10,
    enet_max_iter: int = 5000,
    enet_tol: float = 1e-3,
    enet_selection: str = "random",
    fit_intercept: bool = False,
    save_only_best: bool = True,
    output_dir: Optional[Union[str, Path]] = None,
    normalization: str = "RelativeRepresentation",
) -> pd.DataFrame:
    if protein_name not in y.columns:
        raise ValueError("Protein {} not found in y columns.".format(protein_name))

    # Align y to X rows
    y_protein = y.loc[X_transformed_df.index, protein_name].to_numpy(dtype=np.float64, copy=False)
    X_vals = X_transformed_df.to_numpy(dtype=np.float64, copy=False)
    sample_ids = np.asarray(X_transformed_df.index, dtype=object)
    gene_names = np.asarray(X_transformed_df.columns, dtype=object)

    n_samples, n_features = X_vals.shape
    if n_samples < 2:
        raise ValueError("Need at least 2 samples for LOO CV.")

    if n_splits < 2:
        raise ValueError("n_splits must be >= 2.")
    if (n_samples - 1) < n_splits:
        raise ValueError("n_splits={} is too large for LOO training size {}.".format(n_splits, n_samples - 1))

    n_params = len(param_combinations)
    if n_params == 0:
        param_combinations = [{}]
        n_params = 1

    # Precompute inner splits once (train size is always n_samples-1 in LOO)
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    inner_splits = list(kf.split(np.arange(n_samples - 1, dtype=np.int64)))

    # Fast best-only tracking
    best_abs_err = np.inf
    best_sq_err = np.inf
    best_i = -1
    best_params_dict: Dict[str, Any] = {}
    best_coef_dict: Dict[str, float] = {}

    total_sse = 0.0

    # If user wants full raw results
    raw_rows = [] if not save_only_best else None

    all_idx = np.arange(n_samples, dtype=np.int64)

    # Outer LOO
    for i in range(n_samples):
        train_mask = np.ones(n_samples, dtype=bool)
        train_mask[i] = False

        X_train = X_vals[train_mask]
        y_train = y_protein[train_mask]

        X_test = X_vals[i:i+1, :]
        y_true = float(y_protein[i])

        # Inner CV hyperparam tuning
        avg_mse = np.zeros(n_params, dtype=np.float64)
        folds_used = 0

        for fit_idx, val_idx in inner_splits:
            X_fit = X_train[fit_idx]
            y_fit = y_train[fit_idx]
            X_val = X_train[val_idx]
            y_val = y_train[val_idx]

            if feature_selector != "anova":
                # For now, fall back to the same scoring method
                pass

            sel = _select_feature_indices_percentile(X_fit, y_fit, percentile)
            if sel.size == 0:
                continue

            folds_used += 1
            X_fit_sel = X_fit[:, sel]
            X_val_sel = X_val[:, sel]

            # Evaluate all params for this fold
            for p_idx, params in enumerate(param_combinations):
                est = ElasticNet(
                    max_iter=enet_max_iter,
                    tol=enet_tol,
                    selection=enet_selection,
                    random_state=42,
                    fit_intercept=fit_intercept,
                    **params,
                )
                est.fit(X_fit_sel, y_fit)
                pred_val = est.predict(X_val_sel)
                err = y_val - pred_val
                avg_mse[p_idx] += float(np.mean(err * err))

        if folds_used > 0:
            avg_mse /= float(folds_used)
        else:
            avg_mse[:] = np.inf

        best_params = param_combinations[int(np.argmin(avg_mse))]

        # Final selection on full LOO-train
        final_sel = _select_feature_indices_percentile(X_train, y_train, percentile)
        if final_sel.size == 0:
            # ultra-degenerate fallback
            final_sel = np.array([0], dtype=np.int64)

        X_train_final = X_train[:, final_sel]
        X_test_final = X_test[:, final_sel]

        # Fit final model and predict left-out
        est = ElasticNet(
            max_iter=enet_max_iter,
            tol=enet_tol,
            selection=enet_selection,
            random_state=42,
            fit_intercept=fit_intercept,
            **best_params,
        )
        est.fit(X_train_final, y_train)
        y_pred = float(est.predict(X_test_final)[0])

        err = y_true - y_pred
        abs_err = abs(err)
        sq_err = err * err

        total_sse += sq_err

        # Build outputs for this LOO split (only if requested)
        params_out = dict(best_params)
        if fit_intercept:
            params_out["__intercept__"] = float(est.intercept_)

        if not save_only_best:
            sel_names = gene_names[final_sel]
            coef_dict = {str(g): float(c) for g, c in zip(sel_names, est.coef_)}
            raw_rows.append(
                {
                    "LOO sample": sample_ids[i],
                    "true": y_true,
                    "predicted": y_pred,
                    "params_dict": params_out,
                    "features_coef_dict": coef_dict,
                }
            )

        # Best-only tracking
        if abs_err < best_abs_err:
            best_abs_err = float(abs_err)
            best_sq_err = float(sq_err)
            best_i = int(i)

            sel_names = gene_names[final_sel]
            best_coef_dict = {str(g): float(c) for g, c in zip(sel_names, est.coef_)}
            best_params_dict = params_out

    if save_only_best:
        return _best_only_row_from_tracking(
            protein_name=protein_name,
            sample_ids=sample_ids,
            best_i=best_i,
            best_abs_error=best_abs_err,
            total_sse=total_sse,
            best_sq_error=best_sq_err,
            params_dict=best_params_dict,
            coef_dict=best_coef_dict,
            fit_intercept=fit_intercept,
            output_dir=output_dir,
            normalization=normalization,
        )

    raw = pd.DataFrame(raw_rows)
    # If full raw results, preserve your original post-processing if you want best-only after
    # (caller can still request save_only_best=False)
    return raw


# ============================================================
# Pipeline 1: BINNING (optimized; accepts precomputed transform optionally)
# ============================================================

def cv_pipeline_fixed_percentile_binning(
    param_grid: Dict[str, Sequence[Any]],
    X: pd.DataFrame,
    y: pd.DataFrame,
    protein_name: str,
    *,
    n_bins: int = 50,
    feature_selector: str = "anova",
    n_splits: int = 3,
    percentile: int = 10,
    enet_max_iter: int = 5000,
    enet_tol: float = 1e-3,
    enet_selection: str = "random",
    fit_intercept: bool = False,
    save_only_best: bool = True,
    output_dir: Optional[Union[str, Path]] = None,
    X_binned: Optional[pd.DataFrame] = None,  # NEW: allow caller to pass precomputed
) -> pd.DataFrame:
    if protein_name not in y.columns:
        raise ValueError("Protein {} not found in y columns.".format(protein_name))

    keys, values = zip(*param_grid.items()) if len(param_grid) else ([], [])
    param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)] if len(keys) else [{}]

    X_b = X_binned if X_binned is not None else _apply_binning_matrix(X, n_bins=n_bins)

    return _cv_pipeline_fixed_percentile_pretransformed(
        param_combinations=param_combinations,
        X_transformed_df=X_b,
        y=y,
        protein_name=protein_name,
        feature_selector=feature_selector,
        n_splits=n_splits,
        percentile=percentile,
        enet_max_iter=enet_max_iter,
        enet_tol=enet_tol,
        enet_selection=enet_selection,
        fit_intercept=fit_intercept,
        save_only_best=save_only_best,
        output_dir=output_dir,
        normalization="RelativeRepresentation",
    )


# ============================================================
# Pipeline 2: RANKING (optimized; accepts precomputed transform optionally)
# ============================================================

def cv_pipeline_fixed_percentile_ranking(
    param_grid: Dict[str, Sequence[Any]],
    X: pd.DataFrame,
    y: pd.DataFrame,
    protein_name: str,
    *,
    feature_selector: str = "anova",
    n_splits: int = 3,
    percentile: int = 10,
    enet_max_iter: int = 5000,
    enet_tol: float = 1e-3,
    enet_selection: str = "random",
    fit_intercept: bool = False,
    save_only_best: bool = True,
    output_dir: Optional[Union[str, Path]] = None,
    X_ranked: Optional[pd.DataFrame] = None,  # NEW: allow caller to pass precomputed
) -> pd.DataFrame:
    if protein_name not in y.columns:
        raise ValueError("Protein {} not found in y columns.".format(protein_name))

    keys, values = zip(*param_grid.items()) if len(param_grid) else ([], [])
    param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)] if len(keys) else [{}]

    X_r = X_ranked if X_ranked is not None else _apply_ranking_matrix(X)

    return _cv_pipeline_fixed_percentile_pretransformed(
        param_combinations=param_combinations,
        X_transformed_df=X_r,
        y=y,
        protein_name=protein_name,
        feature_selector=feature_selector,
        n_splits=n_splits,
        percentile=percentile,
        enet_max_iter=enet_max_iter,
        enet_tol=enet_tol,
        enet_selection=enet_selection,
        fit_intercept=fit_intercept,
        save_only_best=save_only_best,
        output_dir=output_dir,
        normalization="RelativeRepresentation",
    )


# ============================================================
# Wrapper Logic (major speedup: transform once for all proteins)
# ============================================================

def run_protein_group_norm(
    proteins: Sequence[str],
    pipeline: str,
    X: pd.DataFrame,
    Y: pd.DataFrame,
    *,
    param_grid: Dict[str, Sequence[Any]],
    feature_selector: str = "anova",
    percentile: int = 5,
    n_splits: int = 3,
    n_bins: int = 50,
    enet_max_iter: int = 5000,
    enet_tol: float = 1e-3,
    enet_selection: str = "random",
    fit_intercept: bool = False,
    save_only_best: bool = True,
    output_dir: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """
    Run the chosen Relative Representation pipeline (Binning or Ranking).
    OPTIMIZED: compute the transform ONCE and reuse for all proteins.
    """
    pipeline_l = str(pipeline).lower()

    keys, values = zip(*param_grid.items()) if len(param_grid) else ([], [])
    param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)] if len(keys) else [{}]

    # Precompute transformation once
    if pipeline_l == "binning":
        X_trans = _apply_binning_matrix(X, n_bins=n_bins)
    elif pipeline_l == "ranking":
        X_trans = _apply_ranking_matrix(X)
    else:
        raise ValueError("Unknown pipeline: {}. Use 'binning' or 'ranking'.".format(pipeline))

    outs = []
    for prot in proteins:
        if prot not in Y.columns:
            outs.append(pd.DataFrame([{"protein": prot, "error": "missing_in_Y"}]))
            continue

        dfp = _cv_pipeline_fixed_percentile_pretransformed(
            param_combinations=param_combinations,
            X_transformed_df=X_trans,
            y=Y,
            protein_name=prot,
            feature_selector=feature_selector,
            n_splits=n_splits,
            percentile=percentile,
            enet_max_iter=enet_max_iter,
            enet_tol=enet_tol,
            enet_selection=enet_selection,
            fit_intercept=fit_intercept,
            save_only_best=save_only_best,
            output_dir=output_dir,
            normalization="RelativeRepresentation",
        )

        # Ensure protein column exists (best-only already has it; raw has LOO rows)
        if "protein" not in dfp.columns:
            dfp = dfp.copy()
            dfp.insert(0, "protein", prot)

        outs.append(dfp)

    if not outs:
        return pd.DataFrame()

    return pd.concat(outs, ignore_index=True)
