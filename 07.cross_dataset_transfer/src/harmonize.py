"""
harmonize.py -- cross-dataset harmonization, reference-anchored.

Stage A of the transfer pipeline. Brings every dataset (CCLE, TMT, iTRAQ,
LabelFree) onto a common scale by fitting the transform ONLY on an anchor
reference that is always part of training (CCLE), then applying the frozen
transform to every other dataset -- train and test alike. This is the same
fit-on-ref / apply-to-target logic used in `05. validation` (fit_on_ref_apply),
generalized so the reference is the training anchor and the "target" can be any
train-union member or the held-out test dataset.

Leakage note: the anchor (CCLE) is always in the training union, so estimating
the scaler on it never touches held-out labels. Optional peak-alignment (mean
re-centering on the target) uses only unsupervised target *feature* means; it is
exposed as a flag and OFF by default for the strict, no-test-statistics setting.
Turn it on to reproduce the 05.validation behavior.

Learned parameters are captured in a plain dict so the notebook can save them as
"global" harmonization params and the HPC pipeline can re-apply them unchanged.
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler, RobustScaler, QuantileTransformer


def _sample_center_df(df: pd.DataFrame, how: Optional[str]) -> pd.DataFrame:
    if how is None:
        return df
    if how == "median":
        return df.sub(df.median(axis=1), axis=0)
    if how == "mean":
        return df.sub(df.mean(axis=1), axis=0)
    raise ValueError(how)


def fit_reference_normalizer(
    df_ref: pd.DataFrame,
    method: str = "robust_z",
    transform: Optional[str] = None,     # None | "log2p1" | "asinh"
    clip_nonneg: bool = False,
    sample_center: Optional[str] = None,
    clip_to_ref_quantiles: Optional[Tuple[float, float]] = None,
    n_quantiles: int = 1000,
    random_state: int = 0,
) -> Dict:
    """
    Fit a harmonization transform on the reference matrix (samples x features)
    and return a serializable parameter bundle. Nothing about any target is used.
    """
    if df_ref.shape[1] == 0:
        raise ValueError(
            "fit_reference_normalizer got a reference matrix with 0 features "
            f"(shape {df_ref.shape}). This almost always means the common "
            "gene/protein intersection across the datasets in this split is "
            "empty -- e.g. protein IDs do not match between CCLE and the CPTAC "
            "platforms. Run diagnose_feature_overlap() in the notebook to see "
            "where the intersection collapses, and fix ID formatting (symbol vs "
            "Ensembl/UniProt, case, version suffixes) before harmonizing."
        )
    if df_ref.shape[0] == 0:
        raise ValueError(
            f"fit_reference_normalizer got 0 samples (shape {df_ref.shape}); "
            "the train-union concatenation is empty for this split."
        )
    X_r = df_ref.copy()
    X_r = _sample_center_df(X_r, sample_center)

    if transform == "log2p1":
        if clip_nonneg:
            X_r = X_r.clip(lower=0)
        X_r = np.log2(X_r + 1.0)
    elif transform == "asinh":
        X_r = np.arcsinh(X_r)
    elif transform is not None:
        raise ValueError(transform)

    ref_med = X_r.median(axis=0)
    X_r = X_r.fillna(ref_med)

    lo = hi = None
    if clip_to_ref_quantiles is not None:
        qlo, qhi = clip_to_ref_quantiles
        lo = X_r.quantile(qlo, axis=0)
        hi = X_r.quantile(qhi, axis=0)
        X_r = X_r.clip(lower=lo, upper=hi, axis=1)

    if method == "standard_z":
        scaler = StandardScaler().fit(X_r)
    elif method == "robust_z":
        scaler = RobustScaler(with_centering=True, with_scaling=True,
                              quantile_range=(25, 75)).fit(X_r)
    elif method == "quantile_normal":
        nq = int(min(n_quantiles, max(10, X_r.shape[0])))
        qt = QuantileTransformer(n_quantiles=nq, output_distribution="normal",
                                 subsample=int(1e9), random_state=random_state).fit(X_r)
        Q_r = qt.transform(X_r)
        scaler = StandardScaler().fit(Q_r)
        return {
            "method": method, "transform": transform, "clip_nonneg": clip_nonneg,
            "sample_center": sample_center, "columns": list(map(str, X_r.columns)),
            "ref_median": ref_med.to_dict(),
            "clip_lo": None if lo is None else lo.to_dict(),
            "clip_hi": None if hi is None else hi.to_dict(),
            "quantile_transformer": qt, "scaler": scaler,
            "ref_mean_after": np.asarray(scaler.transform(Q_r)).mean(axis=0),
        }
    else:
        raise ValueError(method)

    Z_r = scaler.transform(X_r)
    return {
        "method": method, "transform": transform, "clip_nonneg": clip_nonneg,
        "sample_center": sample_center, "columns": list(map(str, X_r.columns)),
        "ref_median": ref_med.to_dict(),
        "clip_lo": None if lo is None else lo.to_dict(),
        "clip_hi": None if hi is None else hi.to_dict(),
        "scaler": scaler, "quantile_transformer": None,
        "ref_mean_after": np.asarray(Z_r).mean(axis=0),
    }


def apply_reference_normalizer(
    df_target: pd.DataFrame,
    params: Dict,
    peak_align: bool = False,
) -> pd.DataFrame:
    """
    Apply a frozen harmonization bundle (from fit_reference_normalizer) to any
    target matrix. Features are aligned to the reference column order; missing
    reference features are imputed with the reference median.
    """
    cols = list(map(str, params["columns"]))
    X_t = df_target.copy()
    X_t.columns = X_t.columns.map(str)
    X_t = X_t.reindex(columns=cols)

    X_t = _sample_center_df(X_t, params.get("sample_center"))

    transform = params.get("transform")
    if transform == "log2p1":
        if params.get("clip_nonneg"):
            X_t = X_t.clip(lower=0)
        X_t = np.log2(X_t + 1.0)
    elif transform == "asinh":
        X_t = np.arcsinh(X_t)

    ref_med = pd.Series(params["ref_median"]).reindex(cols)
    X_t = X_t.fillna(ref_med)

    if params.get("clip_lo") is not None:
        lo = pd.Series(params["clip_lo"]).reindex(cols)
        hi = pd.Series(params["clip_hi"]).reindex(cols)
        X_t = X_t.clip(lower=lo, upper=hi, axis=1)

    qt = params.get("quantile_transformer")
    if qt is not None:
        X_t = pd.DataFrame(qt.transform(X_t), index=X_t.index, columns=cols)

    Z_t = params["scaler"].transform(X_t)

    if peak_align:
        # 05.validation behaviour: slide target peak onto reference peak.
        Z_t = Z_t - (np.nanmean(Z_t, axis=0) - np.asarray(params["ref_mean_after"]))

    return pd.DataFrame(Z_t, index=df_target.index, columns=cols)


# --------------------------------------------------------------------------- #
# Diagnostics: why did the common feature space collapse?
# --------------------------------------------------------------------------- #
def diagnose_feature_overlap(raw: Dict[str, Dict[str, pd.DataFrame]],
                             datasets, modality: str = "Protein", n_show: int = 8):
    """
    Print, for the datasets in a split, how the shared feature set shrinks as
    each dataset is intersected in, plus a sample of each dataset's feature IDs
    so mismatches (symbol vs Ensembl/UniProt, case, version suffixes) are visible.

    raw: {dataset: {"RNA": df, "Protein": df}}  (samples x features)
    Returns the final common feature list.
    """
    print(f"=== feature overlap diagnosis: modality={modality} ===")
    running = None
    for d in datasets:
        cols = set(map(str, raw[d][modality].columns))
        print(f"{d:10s} {modality} features: {len(cols):6d} | e.g. "
              f"{sorted(list(cols))[:n_show]}")
        running = cols if running is None else (running & cols)
        print(f"   running intersection after {d}: {len(running)}")
    running = sorted(running or [])
    if not running:
        print("!! EMPTY intersection. IDs are formatted differently across "
              "datasets. Harmonize IDs (uppercase gene symbols; strip Ensembl "
              "version suffixes like .3; map UniProt->symbol) before matching.")
    else:
        print(f"final common {modality}: {len(running)} (e.g. {running[:n_show]})")
    return running


def normalize_symbols(cols) -> list:
    """
    Best-effort ID cleanup for matching: uppercase, strip whitespace, drop
    Ensembl version suffixes (e.g. trailing dot-number). Apply the SAME
    function to every dataset's columns before intersecting.
    """
    out = []
    for c in map(str, cols):
        c = c.strip().upper()
        if c.startswith("ENSG") and "." in c:
            c = c.split(".")[0]
        out.append(c)
    return out


def map_ensembl_columns_to_symbol(df: pd.DataFrame, id2sym: Dict[str, str],
                                  collapse: str = "mean") -> pd.DataFrame:
    """
    Map DataFrame COLUMNS from Ensembl gene IDs (ENSG...) to gene symbols using
    id2sym {ensembl_gene_id -> gene_symbol}. Version suffixes (ENSG....3) are
    stripped before lookup. Columns already symbols (or unmapped) are kept as-is.
    Duplicate symbols created by the mapping are collapsed (default: mean).
    Use to bring CPTAC protein tables (Ensembl) into CCLE's symbol namespace.
    """
    def _one(c):
        c = str(c)
        key = c.split(".")[0] if c.upper().startswith("ENSG") else c
        return id2sym.get(key, id2sym.get(c, c))
    out = df.copy()
    out.columns = [_one(c) for c in out.columns]
    if out.columns.duplicated().any():
        if collapse == "mean":
            out = out.T.groupby(level=0).mean().T
        else:
            out = out.loc[:, ~out.columns.duplicated(keep="first")]
    return out
