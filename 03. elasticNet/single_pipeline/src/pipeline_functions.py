# pipeline_functions.py

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Union, Dict, Any, List

from sklearn.model_selection import KFold
from sklearn.feature_selection import SelectPercentile, f_regression, SelectKBest
from sklearn.metrics import mean_squared_error
from sklearn.linear_model import ElasticNet
import itertools

# Fix seed globally
np.random.seed(42)

# ============================================================
# Feature selection
# ============================================================
def _select_features_percentile(
    X_feat_by_sample,  # pd.DataFrame (features x samples)
    y,                 # pd.Series
    feature_selector,  # str
    percentile,        # int
):
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
    X_feat_by_sample,  # pd.DataFrame (features x samples)
    y,                 # pd.Series
    feature_selector,  # str
    p_value,           # float
):
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
def _fit_zscore_params(X_train_df, cols, eps=1e-8):
    Xt = X_train_df[cols].to_numpy(dtype=float, copy=False)
    mu = Xt.mean(axis=0)
    sigma = Xt.std(axis=0)
    sigma = np.where(sigma < eps, 1.0, sigma)
    return mu, sigma


def _apply_zscore(X_df, cols, mu, sigma, eps=1e-8):
    Xa = X_df[cols].to_numpy(dtype=float, copy=False)
    return (Xa - mu) / (sigma + eps)


# ============================================================
# Best-only saving helpers (1-row output)
# ============================================================
def _safe_to_parquet_or_csv(df, path):
    """
    Save parquet (snappy) if possible, else fallback to CSV next to it.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        df.to_parquet(path, compression="snappy")
    except Exception:
        df.to_csv(path.with_suffix(".csv"), index=False)


def _predict_with_coef_dict(
    X,                  # pd.DataFrame
    coef_dict,          # dict
    zscore=True,        # bool
    mu=None,            # Optional[np.ndarray]
    sigma=None,         # Optional[np.ndarray]
    eps=1e-8,
):
    """
    Prediction using coef_dict.

    If zscore=False:
      pred = X[genes] @ w

    If zscore=True:
      pred = zscore(X[genes]; mu,sigma) @ w
      (mu/sigma MUST be from the training split used to fit the coefficients)
    """
    if not isinstance(coef_dict, dict) or len(coef_dict) == 0:
        raise ValueError("coef_dict is empty/invalid.")

    genes = [g for g in coef_dict.keys() if g in X.columns]
    if len(genes) == 0:
        raise ValueError("No coef_dict genes found in X.columns.")

    w = np.asarray([float(coef_dict[g]) for g in genes], dtype=float)

    if not zscore:
        Xmat = X[genes].astype(float).to_numpy()
        return Xmat @ w

    if mu is None or sigma is None:
        raise ValueError("zscore=True requires mu and sigma (fit on training split).")

    Xz = _apply_zscore(X, genes, mu, sigma, eps=eps)
    return Xz @ w


def _best_only_row(
    raw_results,        # pd.DataFrame (full LOO table)
    X,                  # pd.DataFrame
    y_protein,          # pd.Series
    protein_name,       # str
    zscore=True,        # bool
    eps=1e-8,
    output_dir=None,    # Optional[Union[str,Path]]
):
    """
    Returns (and optionally saves) EXACTLY ONE ROW with columns:
      - best_sample
      - best_abs_error
      - mse_all_others
      - params_dict
      - feat_coeffs_dict

    mse_all_others:
      Use best sample's coef dict to predict ALL other samples (no refit),
      compute mean squared error on the "others" only (exclude the best sample).

    If zscore=True:
      coefficients are assumed to have been fit on z-scored features.
      We reconstruct the fold's training split: X_train = X.drop(best_sample),
      fit (mu,sigma) on that training split only, and use it to z-score ALL samples.
    """
    df = raw_results.copy()
    df["error"] = df["true"].astype(float) - df["predicted"].astype(float)
    df["abs_error"] = df["error"].abs()

    best_idx = df["abs_error"].idxmin()
    best_row = df.loc[best_idx]

    best_sample = str(best_row["LOO sample"])
    best_abs_error = float(best_row["abs_error"])
    params_dict = best_row["params_dict"] if isinstance(best_row["params_dict"], dict) else {}
    feat_coeffs_dict = best_row["features_coef_dict"] if isinstance(best_row["features_coef_dict"], dict) else {}

    # Fit z-score params on the training split of that best fold (NO leakage)
    if zscore:
        genes = [g for g in feat_coeffs_dict.keys() if g in X.columns]
        if len(genes) == 0:
            raise ValueError("Best fold has no usable genes in X.columns for z-score prediction.")

        if best_sample not in X.index:
            raise ValueError("Best sample '{}' not found in X.index.".format(best_sample))

        X_train_best = X.drop(index=[best_sample])
        mu, sigma = _fit_zscore_params(X_train_best, genes, eps=eps)
    else:
        mu, sigma = None, None

    # Predict all samples with the best coef dict (no refit)
    pred_all = _predict_with_coef_dict(X, feat_coeffs_dict, zscore=zscore, mu=mu, sigma=sigma, eps=eps)
    y_all = y_protein.loc[X.index].astype(float).to_numpy()

    # Exclude best sample from mse_all_others
    mask_others = np.array([s != best_sample for s in X.index], dtype=bool)
    if mask_others.any():
        err_others = y_all[mask_others] - pred_all[mask_others]
        mse_all_others = float(np.mean(err_others ** 2))
    else:
        mse_all_others = np.nan

    out = pd.DataFrame([{
        "best_sample": best_sample,
        "best_abs_error": best_abs_error,
        "mse_all_others": mse_all_others,
        "params_dict": params_dict,
        "feat_coeffs_dict": feat_coeffs_dict,
    }])

    if output_dir is not None:
        out_path = Path(output_dir) / "{}_best.parquet".format(protein_name)
        _safe_to_parquet_or_csv(out, out_path)

    return out


# ============================================================
# Pipelines
# ============================================================
def cv_pipeline_fixed_percentile(
    model,
    param_grid,
    X,                          # pd.DataFrame
    y,                          # pd.DataFrame
    protein_name,               # str
    feature_selector="anova",
    n_splits=3,
    percentile=10,
    enet_max_iter=5000,
    enet_tol=1e-3,
    enet_selection="random",
    fit_intercept=False,
    zscore=True,                # default True
    zscore_eps=1e-8,
    save_only_best=True,
    output_dir=None,            # Optional[Union[str,Path]]
):
    """
    LOO outer loop + inner CV hyperparam selection.

    If zscore=True (default):
      - z-score fit on X_fit only; applied to X_val (inner CV)  => NO leakage
      - final model fit: z-score fit on X_train only; applied to X_test => NO leakage

    save_only_best=True returns/saves 1-row summary (best sample + mse_all_others).
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
                X_fit_mat = X_fit[sel].astype(float).to_numpy()
                X_val_mat = X_val[sel].astype(float).to_numpy()

            for i, params in enumerate(param_combinations):
                est = ElasticNet(
                    max_iter=enet_max_iter,
                    tol=enet_tol,
                    selection=enet_selection,
                    random_state=42,
                    fit_intercept=fit_intercept,
                    **params
                )
                est.fit(X_fit_mat, y_fit.values)
                pred_val = est.predict(X_val_mat)
                avg_mse[i] += mean_squared_error(y_val.values, pred_val)

        if folds_used > 0:
            avg_mse /= folds_used
        else:
            avg_mse[:] = np.inf

        best_params = param_combinations[int(np.argmin(avg_mse))]

        # Final train for this LOO split
        final_sel = _select_features_percentile(X_train.T, y_train, feature_selector, percentile)
        if not final_sel:
            final_sel = _select_features_pvalue(X_train.T, y_train, feature_selector, p_value=1.0)

        if zscore:
            mu, sigma = _fit_zscore_params(X_train, final_sel, eps=zscore_eps)
            X_train_mat = _apply_zscore(X_train, final_sel, mu, sigma, eps=zscore_eps)
            X_test_mat = _apply_zscore(X_test, final_sel, mu, sigma, eps=zscore_eps)
        else:
            X_train_mat = X_train[final_sel].astype(float).to_numpy()
            X_test_mat = X_test[final_sel].astype(float).to_numpy()

        est = ElasticNet(
            max_iter=enet_max_iter,
            tol=enet_tol,
            selection=enet_selection,
            random_state=42,
            fit_intercept=fit_intercept,
            **best_params
        )
        est.fit(X_train_mat, y_train.values)

        results.append({
            "LOO sample": sample,
            "true": float(y_protein.loc[sample]),
            "predicted": float(est.predict(X_test_mat)[0]),
            "params_dict": dict(best_params),
            "features_coef_dict": dict(zip(final_sel, est.coef_.tolist())),
        })

    raw = pd.DataFrame(results)

    if save_only_best:
        return _best_only_row(raw, X, y_protein, protein_name, zscore=zscore, eps=zscore_eps, output_dir=output_dir)

    return raw


def cv_pipeline_fixed_p(
    model,
    param_grid,
    X,
    y,
    protein_name,
    feature_selector="anova",
    n_splits=3,
    p_value=0.05,
    enet_max_iter=5000,
    enet_tol=1e-3,
    enet_selection="random",
    fit_intercept=False,
    zscore=True,
    zscore_eps=1e-8,
    save_only_best=False,
    output_dir=None,
):
    """Same output behavior as cv_pipeline_fixed_percentile."""
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
                X_fit_mat = X_fit[sel].astype(float).to_numpy()
                X_val_mat = X_val[sel].astype(float).to_numpy()

            for i, params in enumerate(param_combinations):
                est = ElasticNet(
                    max_iter=enet_max_iter,
                    tol=enet_tol,
                    selection=enet_selection,
                    random_state=42,
                    fit_intercept=fit_intercept,
                    **params
                )
                est.fit(X_fit_mat, y_fit.values)
                pred_val = est.predict(X_val_mat)
                avg_mse[i] += mean_squared_error(y_val.values, pred_val)

        if folds_used > 0:
            avg_mse /= folds_used
        else:
            avg_mse[:] = np.inf

        best_params = param_combinations[int(np.argmin(avg_mse))]

        final_sel = _select_features_pvalue(X_train.T, y_train, feature_selector, p_value)

        if zscore:
            mu, sigma = _fit_zscore_params(X_train, final_sel, eps=zscore_eps)
            X_train_mat = _apply_zscore(X_train, final_sel, mu, sigma, eps=zscore_eps)
            X_test_mat = _apply_zscore(X_test, final_sel, mu, sigma, eps=zscore_eps)
        else:
            X_train_mat = X_train[final_sel].astype(float).to_numpy()
            X_test_mat = X_test[final_sel].astype(float).to_numpy()

        est = ElasticNet(
            max_iter=enet_max_iter,
            tol=enet_tol,
            selection=enet_selection,
            random_state=42,
            fit_intercept=fit_intercept,
            **best_params
        )
        est.fit(X_train_mat, y_train.values)

        results.append({
            "LOO sample": sample,
            "true": float(y_protein.loc[sample]),
            "predicted": float(est.predict(X_test_mat)[0]),
            "params_dict": dict(best_params),
            "features_coef_dict": dict(zip(final_sel, est.coef_.tolist())),
        })

    raw = pd.DataFrame(results)

    if save_only_best:
        return _best_only_row(raw, X, y_protein, protein_name, zscore=zscore, eps=zscore_eps, output_dir=output_dir)

    return raw


def cv_pipeline_fixed_features(
    model,
    param_grid,
    X,
    y,
    protein_name,
    features,                    # pd.DataFrame with column "features"
    n_splits=3,
    n_samples=None,
    enet_max_iter=5000,
    enet_tol=1e-3,
    enet_selection="random",
    fit_intercept=False,
    zscore=True,
    zscore_eps=1e-8,
    save_only_best=False,
    output_dir=None,
):
    """
    LOO CV for a single protein with a fixed set of features (provided per sample).
    If zscore=True (default), z-score fit on X_fit (inner) / X_train (final) only => NO leakage.
    """
    if protein_name not in y.columns:
        raise ValueError("Protein {} not found in y columns.".format(protein_name))
    if "features" not in features.columns:
        raise ValueError("`features` DataFrame must contain a 'features' column with lists of feature names per sample.")

    y_protein = y[protein_name]
    sample_list = X.index if n_samples is None else X.index[:n_samples]
    X = X.loc[sample_list]
    y_protein = y_protein.loc[sample_list]

    keys, values = zip(*param_grid.items()) if len(param_grid) else ([], [])
    param_combinations = [dict(zip(keys, v)) for v in itertools.product(*values)] if len(keys) else [{}]

    results = []

    for sample in sample_list:
        X_train = X.drop(index=[sample])
        X_test = X.loc[[sample]]
        y_train = y_protein.drop(index=[sample])
        y_test = y_protein.loc[[sample]]

        selected_features = features.loc[sample, "features"]
        if not isinstance(selected_features, (list, tuple, np.ndarray)):
            raise ValueError("Each entry in `features['features']` must be a list/tuple/array of feature names.")

        selected_features = [f for f in selected_features if f in X_train.columns]
        if len(selected_features) == 0:
            results.append({
                "LOO sample": sample,
                "true": float(y_test.values[0]),
                "predicted": np.nan,
                "params_dict": {},
                "features_coef_dict": {},
            })
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
                X_fit_mat = X_fit[selected_features].astype(float).to_numpy()
                X_val_mat = X_val[selected_features].astype(float).to_numpy()

            for i, params in enumerate(param_combinations):
                est = ElasticNet(
                    max_iter=enet_max_iter,
                    tol=enet_tol,
                    selection=enet_selection,
                    random_state=42,
                    fit_intercept=fit_intercept,
                    **params
                )
                try:
                    est.fit(X_fit_mat, y_fit.values)
                    pred_val = est.predict(X_val_mat)
                    avg_mse[i] += mean_squared_error(y_val.values, pred_val)
                except Exception:
                    avg_mse[i] += 1e6

        if folds_used > 0:
            avg_mse /= folds_used
        else:
            avg_mse[:] = np.inf

        best_params = param_combinations[int(np.argmin(avg_mse))]

        if zscore:
            mu, sigma = _fit_zscore_params(X_train, selected_features, eps=zscore_eps)
            X_train_mat = _apply_zscore(X_train, selected_features, mu, sigma, eps=zscore_eps)
            X_test_mat = _apply_zscore(X_test, selected_features, mu, sigma, eps=zscore_eps)
        else:
            X_train_mat = X_train[selected_features].astype(float).to_numpy()
            X_test_mat = X_test[selected_features].astype(float).to_numpy()

        est = ElasticNet(
            max_iter=enet_max_iter,
            tol=enet_tol,
            selection=enet_selection,
            random_state=42,
            fit_intercept=fit_intercept,
            **best_params
        )
        try:
            est.fit(X_train_mat, y_train.values)
            y_pred = float(est.predict(X_test_mat)[0])
            coef_dict = dict(zip(selected_features, est.coef_.tolist()))
        except Exception:
            y_pred = np.nan
            coef_dict = {}

        results.append({
            "LOO sample": sample,
            "true": float(y_test.values[0]),
            "predicted": y_pred,
            "params_dict": dict(best_params),
            "features_coef_dict": coef_dict,
        })

    raw = pd.DataFrame(results)

    if save_only_best:
        return _best_only_row(raw, X, y_protein, protein_name, zscore=zscore, eps=zscore_eps, output_dir=output_dir)

    return raw


# ============================================================
# New: cv_union (2-stage LOO)
# ============================================================
def cv_union(
    model,
    param_grid,
    X,
    y,
    protein_name,
    feature_selector="anova",
    percentile=5,
    n_splits_stage1=3,
    n_splits_stage2=3,
    enet_max_iter=5000,
    enet_tol=1e-3,
    enet_selection="random",
    fit_intercept=False,
    zscore=True,
    zscore_eps=1e-8,
    save_only_best=True,
    output_dir=None,
):
    """
    2-stage LOO:
      Stage 1: cv_pipeline_fixed_percentile (save_only_best=False) to get per-fold feature sets
      Stage 2: cv_pipeline_fixed_features using the UNION of Stage-1 selected features

    If save_only_best=True:
      outputs ONE ROW (best_sample, best_abs_error, mse_all_others, params_dict, feat_coeffs_dict)
      and saves "{protein}_best.parquet" if output_dir is provided.
    """
    raw_stage1 = cv_pipeline_fixed_percentile(
        model=model,
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
        model=model,
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
        return _best_only_row(raw_stage2, X, y[protein_name], protein_name, zscore=zscore, eps=zscore_eps, output_dir=output_dir)

    return raw_stage2
