# new_functions.py
# Compatible with Python 3.8+ (no PEP604 | type unions)

import itertools
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, SelectPercentile, f_regression
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import KFold

# Fix seed globally
np.random.seed(42)


# ============================================================
# Naive cluster-based LOO feature selection (NEW)
# ============================================================
MAX_FEATURES: Dict[int, int] = {
    1: 8600,
    2: 10400,
    3: 6700,
}


def _load_clusters_map(clusters_csv_path: Union[str, Path]) -> pd.Series:
    """
    Load a clusters CSV into a Series mapping: protein_name -> cluster_id.

    Tries to be robust to column naming:
      - protein column: 'protein' or 'protein_name' else first column
      - cluster column: 'cluster_id' or 'cluster' else the first non-protein column
    """
    p = Path(clusters_csv_path)
    if not p.exists():
        raise FileNotFoundError("clusters_csv_path not found: {}".format(str(p)))

    df = pd.read_csv(p)

    if df.shape[1] < 2:
        raise ValueError("clusters_csv_path must have at least 2 columns (protein, cluster_id).")

    # protein column
    if "protein" in df.columns:
        prot_col = "protein"
    elif "protein_name" in df.columns:
        prot_col = "protein_name"
    else:
        prot_col = df.columns[0]

    # cluster column
    cluster_candidates = ["cluster_id", "cluster", "Cluster", "clusterId", "clusterID"]
    cluster_col = None
    for c in cluster_candidates:
        if c in df.columns and c != prot_col:
            cluster_col = c
            break
    if cluster_col is None:
        non_prot_cols = [c for c in df.columns if c != prot_col]
        cluster_col = non_prot_cols[0]

    s = df.set_index(prot_col)[cluster_col]
    # normalize index to str for reliable lookup
    s.index = s.index.astype(str)
    return s


def naive_cluster_feature_selection_indices(
    X: pd.DataFrame,
    y: pd.DataFrame,
    protein_name: str,
    *,
    clusters_csv_path: Union[str, Path],
    n_samples: Optional[int] = None,
) -> dict:
    """
    Highly optimized LOO feature selection using vectorized ANOVA (f_regression),
    keeping top-N features where N depends on the protein's cluster id.

    Returns a CSR-like structure:
      - samples: array of sample ids (row labels)
      - indptr: offsets into indices per sample (len = n_samples + 1)
      - indices: concatenated feature indices (column positions into X)
      - n_genes: total number of genes/features in X
    """
    if protein_name not in y.columns:
        raise ValueError("Protein {} not found in y columns.".format(protein_name))

    # 1) Resolve cluster and N features
    clusters_map = _load_clusters_map(clusters_csv_path)
    if str(protein_name) not in clusters_map.index:
        raise ValueError("Protein {} not present in clusters file.".format(protein_name))

    cluster_id = int(clusters_map.loc[str(protein_name)])
    if cluster_id not in MAX_FEATURES:
        raise ValueError("Cluster id {} not in MAX_FEATURES.".format(cluster_id))

    keep_top = int(MAX_FEATURES[cluster_id])

    # 2) Prepare Data (Slicing)
    sample_idx = X.index if n_samples is None else X.index[: int(n_samples)]
    X_subset = X.loc[sample_idx]
    y_subset = y.loc[sample_idx, protein_name]

    X_values = X_subset.values  # (n_samples, n_genes)
    y_values = y_subset.values  # (n_samples,)

    n_genes = int(X_values.shape[1])
    samples: Sequence[str] = list(sample_idx)

    # 3) Output buffers
    indptr = [0]
    indices: List[int] = []

    # 4) Optimized LOO Loop
    n_rows = int(X_values.shape[0])

    for i in range(n_rows):
        mask = np.ones(n_rows, dtype=bool)
        mask[i] = False

        X_train = X_values[mask]
        y_train = y_values[mask]

        f_scores, _ = f_regression(X_train, y_train)
        np.nan_to_num(f_scores, copy=False, nan=-1.0)

        k = min(int(keep_top), n_genes)
        if k <= 0:
            chosen_idx = np.empty(0, dtype=np.int64)
        else:
            part_indices = np.argpartition(f_scores, n_genes - k)[n_genes - k :]
            top_scores = f_scores[part_indices]
            sorted_order = np.argsort(top_scores)[::-1]
            chosen_idx = part_indices[sorted_order]

        indices.extend(chosen_idx.tolist())
        indptr.append(len(indices))

    return {
        "samples": np.asarray(samples),
        "indptr": np.asarray(indptr, dtype=np.uint32),
        "indices": np.asarray(indices, dtype=np.uint32),
        "n_genes": n_genes,
    }


# ============================================================
# Feature selection
# ============================================================
def _select_features_percentile(
    X_feat_by_sample: pd.DataFrame,  # (features x samples)
    y: pd.Series,
    feature_selector: str,
    percentile: int,
) -> List[str]:
    """Select features using ANOVA (f_regression) percentile ranking."""
    if feature_selector == "anova":
        selector = SelectPercentile(score_func=f_regression, percentile=int(percentile))
        selector.fit(X_feat_by_sample.T, y.values)
        return list(X_feat_by_sample.index[selector.get_support()])

    # Fallback: rank by f_regression scores and keep top percentile
    scores, _ = f_regression(X_feat_by_sample.T, y.values)
    s = pd.Series(scores, index=X_feat_by_sample.index).fillna(0.0)
    thresh = np.percentile(s.values, 100 - percentile)
    return s.index[s.values >= thresh].tolist()


def _select_features_pvalue(
    X_feat_by_sample: pd.DataFrame,  # (features x samples)
    y: pd.Series,
    feature_selector: str,
    p_value: float,
) -> List[str]:
    """Select features using ANOVA p-value threshold."""
    selector = SelectKBest(score_func=f_regression, k="all")
    selector.fit(X_feat_by_sample.T, y.values)
    pv = pd.Series(selector.pvalues_, index=X_feat_by_sample.index).fillna(1.0)
    selected = pv.index[pv.values <= float(p_value)].tolist()
    if not selected:
        selected = [pv.idxmin()]
    return selected


# ============================================================
# Z-score (NO leakage) helpers
# ============================================================
def _fit_zscore_params(
    X_train_df: pd.DataFrame,
    cols: Sequence[str],
    eps: float = 1e-8,
) -> Tuple[np.ndarray, np.ndarray]:
    Xt = X_train_df[list(cols)].to_numpy(dtype=float, copy=False)
    mu = Xt.mean(axis=0)
    sigma = Xt.std(axis=0)
    sigma = np.where(sigma < eps, 1.0, sigma)
    return mu, sigma


def _apply_zscore(
    X_df: pd.DataFrame,
    cols: Sequence[str],
    mu: np.ndarray,
    sigma: np.ndarray,
    eps: float = 1e-8,
) -> np.ndarray:
    Xa = X_df[list(cols)].to_numpy(dtype=float, copy=False)
    return (Xa - mu) / (sigma + eps)


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
# Best-only summary (1 row per protein)
# ============================================================
def _predict_with_coef_dict(
    X: pd.DataFrame,
    coef_dict: Dict[str, float],
    *,
    fit_intercept: bool,
    intercept: float,
    zscore: bool,
    mu: Optional[np.ndarray],
    sigma: Optional[np.ndarray],
    eps: float,
) -> np.ndarray:
    """
    Predict using coef_dict (+ optional intercept), with optional z-scoring.

    - If zscore=False:
        pred = intercept + X[genes] @ w
    - If zscore=True:
        pred = intercept + zscore(X[genes]; mu,sigma) @ w
        (mu/sigma MUST be fit on the training split that produced coef_dict)
    """
    if not isinstance(coef_dict, dict) or len(coef_dict) == 0:
        raise ValueError("coef_dict is empty/invalid.")

    genes = [g for g in coef_dict.keys() if g in X.columns]
    if len(genes) == 0:
        raise ValueError("No coef_dict genes found in X.columns.")

    w = np.asarray([float(coef_dict[g]) for g in genes], dtype=float)

    if zscore:
        if mu is None or sigma is None:
            raise ValueError("zscore=True requires mu and sigma.")
        Xmat = _apply_zscore(X, genes, mu, sigma, eps=eps)
    else:
        Xmat = X[genes].to_numpy(dtype=float, copy=False)

    pred = Xmat @ w
    if fit_intercept:
        pred = pred + float(intercept)
    return pred


def best_only_row_for_protein(
    raw_results: pd.DataFrame,
    X: pd.DataFrame,
    y_protein: pd.Series,
    protein_name: str,
    *,
    zscore: bool,
    zscore_eps: float,
    fit_intercept: bool,
    output_dir: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """
    Build & optionally save EXACTLY ONE ROW with columns:
      - protein
      - best_sample
      - best_abs_error
      - mse_all_others
      - params_dict
      - feat_coeffs_dict
      - zscore_genes
      - zscore_mu
      - zscore_sigma
      - zscore_eps
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

    # --- NEW: compute & store z-score params for THIS best-fold model ---
    zscore_genes = None
    zscore_mu = None
    zscore_sigma = None
    zscore_eps_out = None

    if zscore:
        genes = [g for g in coef_dict.keys() if g in X.columns]
        if len(genes) == 0:
            raise ValueError("Best fold has no usable genes in X.columns for z-score prediction.")
        if best_sample not in X.index:
            raise ValueError("Best sample '{}' not found in X.index.".format(best_sample))

        X_train_best = X.drop(index=[best_sample])
        mu_arr, sigma_arr = _fit_zscore_params(X_train_best, genes, eps=zscore_eps)

        # store as python lists (portable: parquet/csv)
        zscore_genes = list(genes)
        zscore_mu = [float(v) for v in mu_arr.tolist()]
        zscore_sigma = [float(v) for v in sigma_arr.tolist()]
        zscore_eps_out = float(zscore_eps)

        # use them for "predict all"
        mu, sigma = mu_arr, sigma_arr
    else:
        mu, sigma = None, None

    pred_all = _predict_with_coef_dict(
        X,
        coef_dict,
        fit_intercept=fit_intercept,
        intercept=intercept,
        zscore=zscore,
        mu=mu,
        sigma=sigma,
        eps=zscore_eps,
    )

    y_all = y_protein.loc[X.index].astype(float).to_numpy()
    mask_others = np.array([s != best_sample for s in X.index], dtype=bool)

    if mask_others.any():
        err_others = y_all[mask_others] - pred_all[mask_others]
        mse_all_others = float(np.mean(err_others ** 2))
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
            # NEW columns:
            "zscore_genes": zscore_genes,
            "zscore_mu": zscore_mu,
            "zscore_sigma": zscore_sigma,
            "zscore_eps": zscore_eps_out,
        }]
    )

    if output_dir is not None:
        out_path = Path(output_dir) / "{}_best.parquet".format(protein_name)
        _safe_to_parquet_or_csv(out, out_path)

    return out


# ============================================================
# Core per-protein pipelines (LOO)
# ============================================================
def cv_pipeline_fixed_percentile(
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
    zscore: bool = True,
    zscore_eps: float = 1e-8,
    save_only_best: bool = True,
    output_dir: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """
    LOO outer loop + inner CV hyperparam selection.
    """
    if protein_name not in y.columns:
        raise ValueError("Protein {} not found in y columns.".format(protein_name))

    y_protein = y[protein_name]
    sample_list = X.index

    keys, values = zip(*param_grid.items()) if len(param_grid) else ([], [])
    param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)] if len(keys) else [{}]

    results = []

    for sample in sample_list:
        X_train = X.drop(index=[sample])
        y_train = y_protein.drop(index=[sample])
        X_test = X.loc[[sample]]

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

            if zscore:
                mu, sigma = _fit_zscore_params(X_fit, sel, eps=zscore_eps)
                X_fit_mat = _apply_zscore(X_fit, sel, mu, sigma, eps=zscore_eps)
                X_val_mat = _apply_zscore(X_val, sel, mu, sigma, eps=zscore_eps)
            else:
                X_fit_mat = X_fit[sel].to_numpy(dtype=float, copy=False)
                X_val_mat = X_val[sel].to_numpy(dtype=float, copy=False)

            for i, params in enumerate(param_combinations):
                est = ElasticNet(
                    max_iter=enet_max_iter,
                    tol=enet_tol,
                    selection=enet_selection,
                    random_state=42,
                    fit_intercept=fit_intercept,
                    **params,
                )
                est.fit(X_fit_mat, y_fit.values)
                pred_val = est.predict(X_val_mat)
                avg_mse[i] += mean_squared_error(y_val.values, pred_val)

        if folds_used > 0:
            avg_mse /= float(folds_used)
        else:
            avg_mse[:] = np.inf

        best_params = param_combinations[int(np.argmin(avg_mse))]

        final_sel = _select_features_percentile(X_train.T, y_train, feature_selector, percentile)
        if not final_sel:
            final_sel = _select_features_pvalue(X_train.T, y_train, feature_selector, p_value=1.0)

        if zscore:
            mu, sigma = _fit_zscore_params(X_train, final_sel, eps=zscore_eps)
            X_train_mat = _apply_zscore(X_train, final_sel, mu, sigma, eps=zscore_eps)
            X_test_mat = _apply_zscore(X_test, final_sel, mu, sigma, eps=zscore_eps)
        else:
            X_train_mat = X_train[final_sel].to_numpy(dtype=float, copy=False)
            X_test_mat = X_test[final_sel].to_numpy(dtype=float, copy=False)

        est = ElasticNet(
            max_iter=enet_max_iter,
            tol=enet_tol,
            selection=enet_selection,
            random_state=42,
            fit_intercept=fit_intercept,
            **best_params,
        )
        est.fit(X_train_mat, y_train.values)

        params_out = dict(best_params)
        if fit_intercept:
            params_out["__intercept__"] = float(est.intercept_)

        results.append(
            {
                "LOO sample": sample,
                "true": float(y_protein.loc[sample]),
                "predicted": float(est.predict(X_test_mat)[0]),
                "params_dict": params_out,
                "features_coef_dict": dict(zip(final_sel, est.coef_.tolist())),
            }
        )

    raw = pd.DataFrame(results)

    if save_only_best:
        return best_only_row_for_protein(
            raw,
            X,
            y_protein,
            protein_name,
            zscore=zscore,
            zscore_eps=zscore_eps,
            fit_intercept=fit_intercept,
            output_dir=output_dir,
        )

    return raw


def cv_pipeline_fixed_p(
    param_grid: Dict[str, Sequence[Any]],
    X: pd.DataFrame,
    y: pd.DataFrame,
    protein_name: str,
    *,
    feature_selector: str = "anova",
    n_splits: int = 3,
    p_value: float = 0.05,
    enet_max_iter: int = 5000,
    enet_tol: float = 1e-3,
    enet_selection: str = "random",
    fit_intercept: bool = False,
    zscore: bool = True,
    zscore_eps: float = 1e-8,
    save_only_best: bool = True,
    output_dir: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """Same behavior as cv_pipeline_fixed_percentile, but feature selection uses a p-value threshold."""
    if protein_name not in y.columns:
        raise ValueError("Protein {} not found in y columns.".format(protein_name))

    y_protein = y[protein_name]
    sample_list = X.index

    keys, values = zip(*param_grid.items()) if len(param_grid) else ([], [])
    param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)] if len(keys) else [{}]

    results = []

    for sample in sample_list:
        X_train = X.drop(index=[sample])
        y_train = y_protein.drop(index=[sample])
        X_test = X.loc[[sample]]

        kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
        avg_mse = np.zeros(len(param_combinations), dtype=float)
        folds_used = 0

        for fit_idx, val_idx in kf.split(X_train):
            X_fit = X_train.iloc[fit_idx]
            X_val = X_train.iloc[val_idx]
            y_fit = y_train.iloc[fit_idx]
            y_val = y_train.iloc[val_idx]

            sel = _select_features_pvalue(X_fit.T, y_fit, feature_selector, p_value)
            if not sel:
                continue

            folds_used += 1

            if zscore:
                mu, sigma = _fit_zscore_params(X_fit, sel, eps=zscore_eps)
                X_fit_mat = _apply_zscore(X_fit, sel, mu, sigma, eps=zscore_eps)
                X_val_mat = _apply_zscore(X_val, sel, mu, sigma, eps=zscore_eps)
            else:
                X_fit_mat = X_fit[sel].to_numpy(dtype=float, copy=False)
                X_val_mat = X_val[sel].to_numpy(dtype=float, copy=False)

            for i, params in enumerate(param_combinations):
                est = ElasticNet(
                    max_iter=enet_max_iter,
                    tol=enet_tol,
                    selection=enet_selection,
                    random_state=42,
                    fit_intercept=fit_intercept,
                    **params,
                )
                est.fit(X_fit_mat, y_fit.values)
                pred_val = est.predict(X_val_mat)
                avg_mse[i] += mean_squared_error(y_val.values, pred_val)

        if folds_used > 0:
            avg_mse /= float(folds_used)
        else:
            avg_mse[:] = np.inf

        best_params = param_combinations[int(np.argmin(avg_mse))]

        final_sel = _select_features_pvalue(X_train.T, y_train, feature_selector, p_value)

        if zscore:
            mu, sigma = _fit_zscore_params(X_train, final_sel, eps=zscore_eps)
            X_train_mat = _apply_zscore(X_train, final_sel, mu, sigma, eps=zscore_eps)
            X_test_mat = _apply_zscore(X_test, final_sel, mu, sigma, eps=zscore_eps)
        else:
            X_train_mat = X_train[final_sel].to_numpy(dtype=float, copy=False)
            X_test_mat = X_test[final_sel].to_numpy(dtype=float, copy=False)

        est = ElasticNet(
            max_iter=enet_max_iter,
            tol=enet_tol,
            selection=enet_selection,
            random_state=42,
            fit_intercept=fit_intercept,
            **best_params,
        )
        est.fit(X_train_mat, y_train.values)

        params_out = dict(best_params)
        if fit_intercept:
            params_out["__intercept__"] = float(est.intercept_)

        results.append(
            {
                "LOO sample": sample,
                "true": float(y_protein.loc[sample]),
                "predicted": float(est.predict(X_test_mat)[0]),
                "params_dict": params_out,
                "features_coef_dict": dict(zip(final_sel, est.coef_.tolist())),
            }
        )

    raw = pd.DataFrame(results)

    if save_only_best:
        return best_only_row_for_protein(
            raw,
            X,
            y_protein,
            protein_name,
            zscore=zscore,
            zscore_eps=zscore_eps,
            fit_intercept=fit_intercept,
            output_dir=output_dir,
        )

    return raw


def cv_pipeline_fixed_features(
    param_grid: Dict[str, Sequence[Any]],
    X: pd.DataFrame,
    y: pd.DataFrame,
    protein_name: str,
    features: pd.DataFrame,
    *,
    n_splits: int = 3,
    enet_max_iter: int = 5000,
    enet_tol: float = 1e-3,
    enet_selection: str = "random",
    fit_intercept: bool = False,
    zscore: bool = True,
    zscore_eps: float = 1e-8,
    save_only_best: bool = True,
    output_dir: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """Fixed feature list LOO pipeline."""
    if protein_name not in y.columns:
        raise ValueError("Protein {} not found in y columns.".format(protein_name))
    if "features" not in features.columns:
        raise ValueError("features DataFrame must contain a 'features' column.")

    y_protein = y[protein_name]
    sample_list = X.index

    keys, values = zip(*param_grid.items()) if len(param_grid) else ([], [])
    param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)] if len(keys) else [{}]

    results = []

    for sample in sample_list:
        X_train = X.drop(index=[sample])
        X_test = X.loc[[sample]]
        y_train = y_protein.drop(index=[sample])

        selected_features = features.loc[sample, "features"]
        if not isinstance(selected_features, (list, tuple, np.ndarray)):
            raise ValueError("Each entry in features['features'] must be a list/tuple/array.")

        selected_features = [f for f in selected_features if f in X_train.columns]
        if len(selected_features) == 0:
            results.append(
                {
                    "LOO sample": sample,
                    "true": float(y_protein.loc[sample]),
                    "predicted": np.nan,
                    "params_dict": {},
                    "features_coef_dict": {},
                }
            )
            continue

        kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
        avg_mse = np.zeros(len(param_combinations), dtype=float)
        folds_used = 0

        for fit_idx, val_idx in kf.split(X_train):
            X_fit = X_train.iloc[fit_idx]
            X_val = X_train.iloc[val_idx]
            y_fit = y_train.iloc[fit_idx]
            y_val = y_train.iloc[val_idx]

            folds_used += 1

            if zscore:
                mu, sigma = _fit_zscore_params(X_fit, selected_features, eps=zscore_eps)
                X_fit_mat = _apply_zscore(X_fit, selected_features, mu, sigma, eps=zscore_eps)
                X_val_mat = _apply_zscore(X_val, selected_features, mu, sigma, eps=zscore_eps)
            else:
                X_fit_mat = X_fit[selected_features].to_numpy(dtype=float, copy=False)
                X_val_mat = X_val[selected_features].to_numpy(dtype=float, copy=False)

            for i, params in enumerate(param_combinations):
                est = ElasticNet(
                    max_iter=enet_max_iter,
                    tol=enet_tol,
                    selection=enet_selection,
                    random_state=42,
                    fit_intercept=fit_intercept,
                    **params,
                )
                try:
                    est.fit(X_fit_mat, y_fit.values)
                    pred_val = est.predict(X_val_mat)
                    avg_mse[i] += mean_squared_error(y_val.values, pred_val)
                except Exception:
                    avg_mse[i] += 1e6

        if folds_used > 0:
            avg_mse /= float(folds_used)
        else:
            avg_mse[:] = np.inf

        best_params = param_combinations[int(np.argmin(avg_mse))]

        if zscore:
            mu, sigma = _fit_zscore_params(X_train, selected_features, eps=zscore_eps)
            X_train_mat = _apply_zscore(X_train, selected_features, mu, sigma, eps=zscore_eps)
            X_test_mat = _apply_zscore(X_test, selected_features, mu, sigma, eps=zscore_eps)
        else:
            X_train_mat = X_train[selected_features].to_numpy(dtype=float, copy=False)
            X_test_mat = X_test[selected_features].to_numpy(dtype=float, copy=False)

        est = ElasticNet(
            max_iter=enet_max_iter,
            tol=enet_tol,
            selection=enet_selection,
            random_state=42,
            fit_intercept=fit_intercept,
            **best_params,
        )
        try:
            est.fit(X_train_mat, y_train.values)
            y_pred = float(est.predict(X_test_mat)[0])
            coef_dict = dict(zip(selected_features, est.coef_.tolist()))
            params_out = dict(best_params)
            if fit_intercept:
                params_out["__intercept__"] = float(est.intercept_)
        except Exception:
            y_pred = np.nan
            coef_dict = {}
            params_out = dict(best_params)

        results.append(
            {
                "LOO sample": sample,
                "true": float(y_protein.loc[sample]),
                "predicted": y_pred,
                "params_dict": params_out,
                "features_coef_dict": coef_dict,
            }
        )

    raw = pd.DataFrame(results)

    if save_only_best:
        return best_only_row_for_protein(
            raw,
            X,
            y_protein,
            protein_name,
            zscore=zscore,
            zscore_eps=zscore_eps,
            fit_intercept=fit_intercept,
            output_dir=output_dir,
        )

    return raw


def cv_union(
    param_grid: Dict[str, Sequence[Any]],
    X: pd.DataFrame,
    y: pd.DataFrame,
    protein_name: str,
    *,
    feature_selector: str = "anova",
    percentile: int = 5,
    n_splits_stage1: int = 3,
    n_splits_stage2: int = 3,
    enet_max_iter: int = 5000,
    enet_tol: float = 1e-3,
    enet_selection: str = "random",
    fit_intercept: bool = False,
    zscore: bool = True,
    zscore_eps: float = 1e-8,
    save_only_best: bool = True,
    output_dir: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """2-stage union pipeline."""
    raw_stage1 = cv_pipeline_fixed_percentile(
        param_grid=param_grid,
        X=X,
        y=y,
        protein_name=protein_name,
        feature_selector=feature_selector,
        n_splits=n_splits_stage1,
        percentile=percentile,
        enet_max_iter=enet_max_iter,
        enet_tol=enet_tol,
        enet_selection=enet_selection,
        fit_intercept=fit_intercept,
        zscore=zscore,
        zscore_eps=zscore_eps,
        save_only_best=False,
        output_dir=None,
    )

    def _keys_or_empty(d):
        return list(d.keys()) if isinstance(d, dict) else []

    union_feats = sorted(set().union(*raw_stage1["features_coef_dict"].apply(_keys_or_empty)))
    union_feats = [f for f in union_feats if f in X.columns]

    features_df = pd.DataFrame(index=X.index)
    features_df["features"] = [union_feats] * len(features_df)

    raw_stage2 = cv_pipeline_fixed_features(
        param_grid=param_grid,
        X=X,
        y=y,
        protein_name=protein_name,
        features=features_df,
        n_splits=n_splits_stage2,
        enet_max_iter=enet_max_iter,
        enet_tol=enet_tol,
        enet_selection=enet_selection,
        fit_intercept=fit_intercept,
        zscore=zscore,
        zscore_eps=zscore_eps,
        save_only_best=False,
        output_dir=None,
    )

    if save_only_best:
        return best_only_row_for_protein(
            raw_stage2,
            X,
            y[protein_name],
            protein_name,
            zscore=zscore,
            zscore_eps=zscore_eps,
            fit_intercept=fit_intercept,
            output_dir=output_dir,
        )

    return raw_stage2


# ============================================================
# NEW pipeline: naive cluster LOO feature selection + fixed-features CV
# ============================================================
def cv_pipeline_naive_cluster(
    param_grid: Dict[str, Sequence[Any]],
    X: pd.DataFrame,
    y: pd.DataFrame,
    protein_name: str,
    *,
    clusters_csv_path: Union[str, Path],
    n_samples: Optional[int] = None,
    n_splits: int = 3,
    enet_max_iter: int = 5000,
    enet_tol: float = 1e-3,
    enet_selection: str = "random",
    fit_intercept: bool = False,
    zscore: bool = True,
    zscore_eps: float = 1e-8,
    save_only_best: bool = True,
    output_dir: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """
    Pipeline that uses naive_cluster_feature_selection_indices() to produce a per-LOO-sample
    feature list (top-N depends on the protein's cluster), then runs the standard fixed-features
    ElasticNet LOO pipeline (with inner CV hyperparam selection).
    """
    # Optional slicing for quick experiments
    if n_samples is not None:
        idx = X.index[: int(n_samples)]
        X = X.loc[idx]
        y = y.loc[idx]

    fs = naive_cluster_feature_selection_indices(
        X=X,
        y=y,
        protein_name=protein_name,
        clusters_csv_path=clusters_csv_path,
        n_samples=None,  # already sliced above if requested
    )

    samples = fs["samples"]
    indptr = fs["indptr"]
    indices = fs["indices"]

    # Build features DataFrame aligned to X.index
    features_df = pd.DataFrame(index=X.index)
    feats_per_sample: List[List[str]] = []

    # Map: sample -> selected genes (by column indices)
    colnames = list(X.columns)
    for i, s in enumerate(samples):
        start = int(indptr[i])
        end = int(indptr[i + 1])
        idxs = indices[start:end].astype(int, copy=False)
        genes = [colnames[j] for j in idxs if 0 <= int(j) < len(colnames)]
        feats_per_sample.append(genes)

    # Ensure alignment to X.index order (samples is derived from X.index[:...], so should match)
    features_df["features"] = feats_per_sample

    return cv_pipeline_fixed_features(
        param_grid=param_grid,
        X=X,
        y=y,
        protein_name=protein_name,
        features=features_df,
        n_splits=n_splits,
        enet_max_iter=enet_max_iter,
        enet_tol=enet_tol,
        enet_selection=enet_selection,
        fit_intercept=fit_intercept,
        zscore=zscore,
        zscore_eps=zscore_eps,
        save_only_best=save_only_best,
        output_dir=output_dir,
    )


def run_protein_group(
    proteins: Sequence[str],
    pipeline: str,
    X: pd.DataFrame,
    Y: pd.DataFrame,
    *,
    param_grid: Dict[str, Sequence[Any]],
    feature_selector: str = "anova",
    percentile: int = 5,
    p_value: float = 0.05,
    n_splits: int = 3,
    n_splits_stage1: int = 3,
    n_splits_stage2: int = 3,
    enet_max_iter: int = 5000,
    enet_tol: float = 1e-3,
    enet_selection: str = "random",
    fit_intercept: bool = False,
    zscore: bool = True,
    zscore_eps: float = 1e-8,
    save_only_best: bool = True,
    output_dir: Optional[Union[str, Path]] = None,
    # NEW: only needed for pipeline="naive_cluster"
    clusters_csv_path: Optional[Union[str, Path]] = None,
    naive_n_samples: Optional[int] = None,
) -> pd.DataFrame:
    """Run the chosen pipeline for a list of proteins and concatenate outputs."""
    pipeline = str(pipeline)

    outs = []
    for prot in proteins:
        if prot not in Y.columns:
            outs.append(pd.DataFrame([{"protein": prot, "error": "missing_in_Y"}]))
            continue

        if pipeline == "fixed_percentile":
            dfp = cv_pipeline_fixed_percentile(
                param_grid=param_grid,
                X=X,
                y=Y,
                protein_name=prot,
                feature_selector=feature_selector,
                n_splits=n_splits,
                percentile=percentile,
                enet_max_iter=enet_max_iter,
                enet_tol=enet_tol,
                enet_selection=enet_selection,
                fit_intercept=fit_intercept,
                zscore=zscore,
                zscore_eps=zscore_eps,
                save_only_best=save_only_best,
                output_dir=output_dir,
            )
        elif pipeline == "fixed_p":
            dfp = cv_pipeline_fixed_p(
                param_grid=param_grid,
                X=X,
                y=Y,
                protein_name=prot,
                feature_selector=feature_selector,
                n_splits=n_splits,
                p_value=p_value,
                enet_max_iter=enet_max_iter,
                enet_tol=enet_tol,
                enet_selection=enet_selection,
                fit_intercept=fit_intercept,
                zscore=zscore,
                zscore_eps=zscore_eps,
                save_only_best=save_only_best,
                output_dir=output_dir,
            )
        elif pipeline == "union":
            dfp = cv_union(
                param_grid=param_grid,
                X=X,
                y=Y,
                protein_name=prot,
                feature_selector=feature_selector,
                percentile=percentile,
                n_splits_stage1=n_splits_stage1,
                n_splits_stage2=n_splits_stage2,
                enet_max_iter=enet_max_iter,
                enet_tol=enet_tol,
                enet_selection=enet_selection,
                fit_intercept=fit_intercept,
                zscore=zscore,
                zscore_eps=zscore_eps,
                save_only_best=save_only_best,
                output_dir=output_dir,
            )
        elif pipeline == "naive_cluster":
            if clusters_csv_path is None:
                raise ValueError("pipeline='naive_cluster' requires clusters_csv_path.")
            dfp = cv_pipeline_naive_cluster(
                param_grid=param_grid,
                X=X,
                y=Y,
                protein_name=prot,
                clusters_csv_path=clusters_csv_path,
                n_samples=naive_n_samples,
                n_splits=n_splits,
                enet_max_iter=enet_max_iter,
                enet_tol=enet_tol,
                enet_selection=enet_selection,
                fit_intercept=fit_intercept,
                zscore=zscore,
                zscore_eps=zscore_eps,
                save_only_best=save_only_best,
                output_dir=output_dir,
            )
        else:
            raise ValueError("Unknown pipeline: {}".format(pipeline))

        if "protein" not in dfp.columns:
            dfp = dfp.copy()
            dfp.insert(0, "protein", prot)

        outs.append(dfp)

    if len(outs) == 0:
        return pd.DataFrame()

    return pd.concat(outs, ignore_index=True)
