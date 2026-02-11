# norm_functions.py
# Compatible with Python 3.8+

import itertools
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

import numpy as np
import pandas as pd
from sklearn.feature_selection import f_regression, SelectPercentile
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import KFold

# Fix seed globally
np.random.seed(42)

# ============================================================
# RELATIVE REPRESENTATION HELPER FUNCTIONS
# Source: Stahlschmidt et al., Page 2 [cite: 5, 75, 76]
# ============================================================

def _apply_binning_matrix(X: pd.DataFrame, n_bins: int = 50) -> pd.DataFrame:
    """
    Applies binning normalization row-wise (per sample).
    Zero values are assigned to Bin 0. Non-zero values are binned into n_bins quantiles.
    
    Ref: "Group transcripts with nonzero values into b quantiles... Zero values are assigned to bin 0." 
    """
    # Work on numpy array for speed
    # We copy to avoid modifying the original dataframe in place
    x_vals = X.values.astype(float).copy()
    
    # -------- BINNING PROCEDURE (User Snippet) ----------
    # x is a (NUM_CELLS, GENES) numpy array
    # Replace zeros with nans so we can use nanquantile function
    x_vals[x_vals == 0] = np.nan
    
    # Create quantiles evenly spaced
    # Note: linspace(0, 1, BINS) creates BINS points. 
    qs = np.linspace(0, 1, int(n_bins))
    
    # nanquantile gives (BINS, NUM_CELLS). We reshape to (BINS, NUM_CELLS, 1) for broadcasting
    # We compare x_vals (NUM_CELLS, GENES) against thresholds
    # The sum counts how many quantiles the value is greater than/equal to.
    
    # Optimization: doing this on the full matrix can be memory intensive. 
    # However, it ensures strict adherence to the provided vectorized logic.
    quantiles = np.nanquantile(x_vals, qs, axis=1) # Shape: (n_bins, n_samples)
    
    # Broadcasting: (n_bins, n_samples, 1) <= (n_samples, n_genes) -> Alignment requires swap
    # x_vals is (n_samples, n_genes). quantiles is (n_bins, n_samples).
    # We need to compare x_vals[i, :] with quantiles[:, i]
    
    # Transpose quantiles to (n_samples, n_bins) and reshape to (n_samples, n_bins, 1)
    # But the snippet used: (x <= np.nanquantile(x, qs, axis=1)[..., None])
    # Let's trace the snippet's dimension logic assuming x is (Cells, Genes):
    # np.nanquantile(x, qs, axis=1) -> (n_bins, n_cells).
    # [..., None] -> (n_bins, n_cells, 1).
    # x (n_cells, n_genes) broadcasts against (n_bins, n_cells, 1).
    # This results in (n_bins, n_cells, n_genes).
    # Summing axis 0 -> (n_cells, n_genes).
    
    # Implementing snippet exactly:
    x_transformed = (x_vals <= np.nanquantile(x_vals, qs, axis=1)[..., None]).sum(axis=0)
    
    # Restore dataframe structure
    return pd.DataFrame(x_transformed, index=X.index, columns=X.columns)


def _apply_ranking_matrix(X: pd.DataFrame) -> pd.DataFrame:
    """
    Applies ranking normalization row-wise (per sample).
    
    Ref: "For ranked representations, gene expression values are ranked in an ascending order." 
    """
    x_vals = X.values.astype(float)
    
    # -------- RANKING PROCEDURE (User Snippet) ----------
    # x.argsort(axis=1).argsort(axis=1) returns 0-based rank indices
    x_transformed = x_vals.argsort(axis=1).argsort(axis=1)
    
    return pd.DataFrame(x_transformed, index=X.index, columns=X.columns)


# ============================================================
# Feature selection (Generic)
# ============================================================
def _select_features_percentile(
    X_feat_by_sample: pd.DataFrame, 
    y: pd.Series,
    feature_selector: str,
    percentile: int,
) -> List[str]:
    """Select features using ANOVA (f_regression) percentile ranking."""
    if feature_selector == "anova":
        selector = SelectPercentile(score_func=f_regression, percentile=int(percentile))
        selector.fit(X_feat_by_sample.T, y.values)
        return list(X_feat_by_sample.index[selector.get_support()])

    # Fallback
    scores, _ = f_regression(X_feat_by_sample.T, y.values)
    s = pd.Series(scores, index=X_feat_by_sample.index).fillna(0.0)
    thresh = np.percentile(s.values, 100 - percentile)
    return s.index[s.values >= thresh].tolist()


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
# Best-only summary 
# ============================================================
def best_only_row_for_protein(
    raw_results: pd.DataFrame,
    X_original: pd.DataFrame, # Pass original for reference, or transformed? 
    # Note: In norm pipelines, the coeffs correspond to the TRANSFORMED data. 
    # But for interpretation we usually track the gene names.
    y_protein: pd.Series,
    protein_name: str,
    *,
    fit_intercept: bool,
    output_dir: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """
    Build & optionally save EXACTLY ONE ROW. 
    Note: zscore params are removed as this script uses Binning/Ranking.
    """
    df = raw_results.copy()
    df["error"] = df["true"].astype(float) - df["predicted"].astype(float)
    df["abs_error"] = df["error"].abs()

    best_idx = df["abs_error"].idxmin()
    best_row = df.loc[best_idx]

    best_sample = str(best_row["LOO sample"])
    best_abs_error = float(best_row["abs_error"])
    
    params_dict = best_row["params_dict"] if isinstance(best_row["params_dict"], dict) else {}
    coef_dict = best_row["features_coef_dict"] if isinstance(best_row["features_coef_dict"], dict) else {}
    intercept = float(params_dict.get("__intercept__", 0.0)) if fit_intercept else 0.0

    # Recalculate predictions for "all others" using the coefficients on the transformed data?
    # Since transformation is sample-independent, we can't just apply coeffs to raw X. 
    # We rely on the pipeline functions to have handled transformation cleanly.
    # For this simplified reporter, we calculate mse_all_others based on the loop results if available, 
    # or skip complex "predict all" logic here to avoid re-transforming the whole matrix inside this helper.
    # We will set mse_all_others to NaN for now to keep it safe, or rely on raw_results statistics.
    
    # Calculate MSE of the other folds (simple average of the errors recorded in raw_results)
    other_rows = df.drop(index=best_idx)
    if len(other_rows) > 0:
        mse_all_others = float(np.mean(other_rows["error"] ** 2))
    else:
        mse_all_others = np.nan

    out = pd.DataFrame(
        [{
            "protein": protein_name,
            "best_sample": best_sample,
            "best_abs_error": best_abs_error,
            "mse_all_others": mse_all_others,
            "params_dict": params_dict,
            "feat_coeffs_dict": coef_dict,
            "normalization": "RelativeRepresentation" 
        }]
    )

    if output_dir is not None:
        out_path = Path(output_dir) / "{}_best.parquet".format(protein_name)
        _safe_to_parquet_or_csv(out, out_path)

    return out


# ============================================================
# Pipeline 1: BINNING
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
) -> pd.DataFrame:
    """
    LOO pipeline using BINNING normalization instead of Z-score.
    Optimization: Since Binning is an SSP (Single Sample Predictor) method[cite: 40],
    we transform the data first. This is valid as per Page 2 of the paper.
    """
    if protein_name not in y.columns:
        raise ValueError("Protein {} not found in y columns.".format(protein_name))

    y_protein = y[protein_name]
    sample_list = X.index

    # 1. Apply Binning Transformation globally (Row-wise, so no leakage)
    X_binned = _apply_binning_matrix(X, n_bins=n_bins)

    keys, values = zip(*param_grid.items()) if len(param_grid) else ([], [])
    param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)] if len(keys) else [{}]

    results = []

    for sample in sample_list:
        # LOO Split on ALREADY TRANSFORMED data
        X_train = X_binned.drop(index=[sample])
        y_train = y_protein.drop(index=[sample])
        X_test = X_binned.loc[[sample]]

        kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
        avg_mse = np.zeros(len(param_combinations), dtype=float)
        folds_used = 0

        # Inner CV for hyperparam tuning
        for fit_idx, val_idx in kf.split(X_train):
            X_fit = X_train.iloc[fit_idx]
            X_val = X_train.iloc[val_idx]
            y_fit = y_train.iloc[fit_idx]
            y_val = y_train.iloc[val_idx]

            # Feature selection on BINNED data (Post-transformation selection)
            # Paper Eq 1 implies Transformed -> HVG (Selection) [cite: 81]
            sel = _select_features_percentile(X_fit.T, y_fit, feature_selector, percentile)
            if not sel:
                continue
            
            folds_used += 1

            X_fit_sel = X_fit[sel].to_numpy(dtype=float, copy=False)
            X_val_sel = X_val[sel].to_numpy(dtype=float, copy=False)

            for i, params in enumerate(param_combinations):
                est = ElasticNet(
                    max_iter=enet_max_iter, tol=enet_tol, selection=enet_selection,
                    random_state=42, fit_intercept=fit_intercept, **params
                )
                est.fit(X_fit_sel, y_fit.values)
                pred_val = est.predict(X_val_sel)
                avg_mse[i] += mean_squared_error(y_val.values, pred_val)

        if folds_used > 0:
            avg_mse /= float(folds_used)
        else:
            avg_mse[:] = np.inf

        # Refit on LOO train
        best_params = param_combinations[int(np.argmin(avg_mse))]
        final_sel = _select_features_percentile(X_train.T, y_train, feature_selector, percentile)
        
        X_train_final = X_train[final_sel].to_numpy(dtype=float)
        X_test_final = X_test[final_sel].to_numpy(dtype=float)

        est = ElasticNet(
            max_iter=enet_max_iter, tol=enet_tol, selection=enet_selection,
            random_state=42, fit_intercept=fit_intercept, **best_params
        )
        est.fit(X_train_final, y_train.values)
        
        params_out = dict(best_params)
        if fit_intercept:
            params_out["__intercept__"] = float(est.intercept_)

        results.append({
            "LOO sample": sample,
            "true": float(y_protein.loc[sample]),
            "predicted": float(est.predict(X_test_final)[0]),
            "params_dict": params_out,
            "features_coef_dict": dict(zip(final_sel, est.coef_.tolist())),
        })

    raw = pd.DataFrame(results)

    if save_only_best:
        return best_only_row_for_protein(
            raw, X_binned, y_protein, protein_name, 
            fit_intercept=fit_intercept, output_dir=output_dir
        )

    return raw


# ============================================================
# Pipeline 2: RANKING
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
) -> pd.DataFrame:
    """
    LOO pipeline using RANKING normalization instead of Z-score.
    Applies row-wise ranking (argsort().argsort()).
    """
    if protein_name not in y.columns:
        raise ValueError("Protein {} not found in y columns.".format(protein_name))

    y_protein = y[protein_name]
    sample_list = X.index

    # 1. Apply Ranking Transformation globally (Row-wise, so no leakage)
    X_ranked = _apply_ranking_matrix(X)

    keys, values = zip(*param_grid.items()) if len(param_grid) else ([], [])
    param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)] if len(keys) else [{}]

    results = []

    for sample in sample_list:
        X_train = X_ranked.drop(index=[sample])
        y_train = y_protein.drop(index=[sample])
        X_test = X_ranked.loc[[sample]]

        kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
        avg_mse = np.zeros(len(param_combinations), dtype=float)
        folds_used = 0

        for fit_idx, val_idx in kf.split(X_train):
            X_fit = X_train.iloc[fit_idx]
            X_val = X_train.iloc[val_idx]
            y_fit = y_train.iloc[fit_idx]
            y_val = y_train.iloc[val_idx]

            sel = _select_features_percentile(X_fit.T, y_fit, feature_selector, percentile)
            if not sel:
                continue
            
            folds_used += 1
            X_fit_sel = X_fit[sel].to_numpy(dtype=float, copy=False)
            X_val_sel = X_val[sel].to_numpy(dtype=float, copy=False)

            for i, params in enumerate(param_combinations):
                est = ElasticNet(
                    max_iter=enet_max_iter, tol=enet_tol, selection=enet_selection,
                    random_state=42, fit_intercept=fit_intercept, **params
                )
                est.fit(X_fit_sel, y_fit.values)
                pred_val = est.predict(X_val_sel)
                avg_mse[i] += mean_squared_error(y_val.values, pred_val)

        if folds_used > 0:
            avg_mse /= float(folds_used)
        else:
            avg_mse[:] = np.inf

        best_params = param_combinations[int(np.argmin(avg_mse))]
        final_sel = _select_features_percentile(X_train.T, y_train, feature_selector, percentile)
        
        X_train_final = X_train[final_sel].to_numpy(dtype=float)
        X_test_final = X_test[final_sel].to_numpy(dtype=float)

        est = ElasticNet(
            max_iter=enet_max_iter, tol=enet_tol, selection=enet_selection,
            random_state=42, fit_intercept=fit_intercept, **best_params
        )
        est.fit(X_train_final, y_train.values)
        
        params_out = dict(best_params)
        if fit_intercept:
            params_out["__intercept__"] = float(est.intercept_)

        results.append({
            "LOO sample": sample,
            "true": float(y_protein.loc[sample]),
            "predicted": float(est.predict(X_test_final)[0]),
            "params_dict": params_out,
            "features_coef_dict": dict(zip(final_sel, est.coef_.tolist())),
        })

    raw = pd.DataFrame(results)

    if save_only_best:
        return best_only_row_for_protein(
            raw, X_ranked, y_protein, protein_name, 
            fit_intercept=fit_intercept, output_dir=output_dir
        )

    return raw


# ============================================================
# Wrapper Logic
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
    n_bins: int = 50, # New param for binning
    enet_max_iter: int = 5000,
    enet_tol: float = 1e-3,
    enet_selection: str = "random",
    fit_intercept: bool = False,
    save_only_best: bool = True,
    output_dir: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """
    Run the chosen Relative Representation pipeline (Binning or Ranking).
    """
    pipeline = str(pipeline).lower()

    outs = []
    for prot in proteins:
        if prot not in Y.columns:
            outs.append(pd.DataFrame([{"protein": prot, "error": "missing_in_Y"}]))
            continue

        if pipeline == "binning":
            dfp = cv_pipeline_fixed_percentile_binning(
                param_grid=param_grid, X=X, y=Y, protein_name=prot,
                n_bins=n_bins,
                feature_selector=feature_selector, n_splits=n_splits, percentile=percentile,
                enet_max_iter=enet_max_iter, enet_tol=enet_tol, enet_selection=enet_selection,
                fit_intercept=fit_intercept, save_only_best=save_only_best, output_dir=output_dir,
            )
        elif pipeline == "ranking":
            dfp = cv_pipeline_fixed_percentile_ranking(
                param_grid=param_grid, X=X, y=Y, protein_name=prot,
                feature_selector=feature_selector, n_splits=n_splits, percentile=percentile,
                enet_max_iter=enet_max_iter, enet_tol=enet_tol, enet_selection=enet_selection,
                fit_intercept=fit_intercept, save_only_best=save_only_best, output_dir=output_dir,
            )
        else:
            raise ValueError("Unknown pipeline for norm_functions: {}. Use 'binning' or 'ranking'.".format(pipeline))

        if "protein" not in dfp.columns:
            dfp = dfp.copy()
            dfp.insert(0, "protein", prot)

        outs.append(dfp)

    if len(outs) == 0:
        return pd.DataFrame()

    return pd.concat(outs, ignore_index=True)