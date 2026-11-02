"""Protein-vs-feature correlations for the chunked correlation pipeline.

Features are either RNA genes (target="rna", -> full_pearson100 / full_spearman100)
or other proteins (target="protein", -> full_prot).

Significance rule (unchanged from the original implementation):
  * bootstrap=True : resample the protein's valid samples n_bootstraps times; keep the
    bootstrap median if the 95% percentile CI excludes 0, else 0. A feature needs at least
    max(10, n_bootstraps//10) usable bootstrap draws.
  * bootstrap=False: keep r if the parametric p-value < p_threshold, else 0.
Computation is vectorised over all features (same results distributionally as the old
per-gene loop; random streams differ, so bootstrap medians are not bit-identical).
"""
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import numpy as np
import pandas as pd
from scipy.stats import rankdata, t as student_t

METHODS = ("pearson", "spearman")


def _sanitize_filename(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_")


def _safe_to_parquet_or_csv(df: pd.DataFrame, path: Path, *, index: bool = False) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        df.to_parquet(path, compression="snappy", index=index)
        return path
    except Exception:
        csv_path = path.with_suffix(".csv")
        df.to_csv(csv_path, index=index)
        return csv_path


def _rank_cols(a: np.ndarray) -> np.ndarray:
    """Average ranks per column; NaNs stay NaN."""
    out = np.full(a.shape, np.nan)
    for j in range(a.shape[1]):
        col = a[:, j]
        ok = ~np.isnan(col)
        if ok.any():
            out[ok, j] = rankdata(col[ok])
    return out


def _masked_corr(x: np.ndarray, Y: np.ndarray, method: str):
    """Correlation of vector x with every column of Y using pairwise-complete rows.
    Returns (r, n_used). r is NaN where undefined (n<2 or zero variance)."""
    W = (~np.isnan(Y)) & (~np.isnan(x))[:, None]
    n = W.sum(0).astype(float)
    if method == "spearman":
        r = np.full(Y.shape[1], np.nan)
        full = W.all(0)
        if full.any():  # fast path: columns fully observed
            xr = rankdata(x)
            Yr = np.apply_along_axis(rankdata, 0, Y[:, full])
            r[full] = _pearson_dense(xr, Yr)
        for j in np.where(~full)[0]:
            m = W[:, j]
            if m.sum() > 1:
                r[j] = _pearson_dense(rankdata(x[m]), rankdata(Y[m, j])[:, None])[0]
        return r, n
    Xm = np.where(W, x[:, None], 0.0)
    Ym = np.where(W, Y, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        sx, sy = Xm.sum(0), Ym.sum(0)
        sxx, syy, sxy = (Xm * Xm).sum(0), (Ym * Ym).sum(0), (Xm * Ym).sum(0)
        cov = sxy - sx * sy / n
        vx = sxx - sx * sx / n
        vy = syy - sy * sy / n
        r = cov / np.sqrt(vx * vy)
    r[(n < 2) | ~(vx > 1e-12) | ~(vy > 1e-12)] = np.nan
    return np.clip(r, -1.0, 1.0), n


def _pearson_dense(x: np.ndarray, Y: np.ndarray) -> np.ndarray:
    xc = x - x.mean()
    Yc = Y - Y.mean(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        den = np.sqrt((xc ** 2).sum() * (Yc ** 2).sum(0))
        r = (xc @ Yc) / den
    r[~(den > 1e-12)] = np.nan
    return np.clip(r, -1.0, 1.0)


def compute_correlation_for_protein(
    transcriptomics: pd.DataFrame,   # samples x features (genes, or proteins for target=protein)
    proteomics: pd.DataFrame,        # samples x proteins
    protein_name: str,
    p_threshold: float = 0.05,
    bootstrap: bool = True,
    n_bootstraps: int = 100,
    random_state: Optional[int] = None,
    method: str = "pearson",
) -> pd.Series:
    """Correlate one protein with every column of `transcriptomics` (same API as before,
    plus `method`). Returns pd.Series(index=features) with r or 0 if not significant."""
    if method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}")
    if protein_name not in proteomics.columns:
        raise ValueError(f"Protein '{protein_name}' not found in proteomics.")

    common_idx = proteomics.index.intersection(transcriptomics.index)
    y = proteomics.loc[common_idx, protein_name]
    feats = transcriptomics.columns
    out = pd.Series(0.0, index=feats, dtype=float)
    y = y.dropna()
    if len(y) <= 1 or np.std(y.to_numpy(float)) == 0:
        return out
    x = y.to_numpy(float)
    Y = transcriptomics.loc[y.index].to_numpy(float)

    if not bootstrap:
        r, n = _masked_corr(x, Y, method)
        with np.errstate(invalid="ignore", divide="ignore"):
            df = n - 2
            tt = r * np.sqrt(df / np.maximum(1e-300, 1.0 - r * r))
            p = 2 * student_t.sf(np.abs(tt), df)
        p[np.abs(r) >= 1.0] = 0.0
        keep = np.isfinite(r) & (df > 0) & (p < p_threshold)
        out.iloc[np.where(keep)[0]] = r[keep]
        return out

    rng = np.random.default_rng(random_state)
    m = len(x)
    boots = np.full((int(n_bootstraps), Y.shape[1]), np.nan)
    for b in range(int(n_bootstraps)):
        idx = rng.integers(0, m, size=m)
        xb = x[idx]
        if np.std(xb) == 0.0:
            continue
        boots[b], _ = _masked_corr(xb, Y[idx], method)
    n_ok = np.isfinite(boots).sum(0)
    with np.errstate(all="ignore"):
        lo = np.nanpercentile(boots, 2.5, axis=0)
        hi = np.nanpercentile(boots, 97.5, axis=0)
        med = np.nanmedian(boots, axis=0)
    keep = (n_ok >= max(10, int(n_bootstraps) // 10)) & ((lo > 0.0) | (hi < 0.0))
    out.iloc[np.where(keep)[0]] = med[keep]
    return out


def run_correlation_group(
    proteins: Sequence[str],
    transcriptomics: pd.DataFrame,
    proteomics: pd.DataFrame,
    *,
    p_threshold: float = 0.05,
    bootstrap: bool = True,
    n_bootstraps: int = 100,
    random_state: Optional[int] = 42,
    output_dir: Optional[Union[str, Path]] = None,
    per_protein_subdir: str = "per_protein",
    per_protein_format: str = "parquet",
    method: str = "pearson",
) -> pd.DataFrame:
    """Loop over proteins; optionally write <output_dir>/<per_protein_subdir>/<protein>_corr.*
    (1 row x features). Returns compact batch table: protein, n_sig, corr_json[, error]."""
    out_rows: List[Dict[str, object]] = []
    out_path = Path(output_dir) if output_dir is not None else None
    if out_path is not None:
        (out_path / per_protein_subdir).mkdir(parents=True, exist_ok=True)

    for i, prot in enumerate(proteins):
        if prot not in proteomics.columns:
            out_rows.append({"protein": prot, "n_sig": 0, "corr_json": json.dumps({}),
                             "error": "missing_in_proteomics"})
            continue
        seed = None if random_state is None else int(random_state) + int(i)
        s = compute_correlation_for_protein(
            transcriptomics, proteomics, prot, p_threshold=p_threshold, bootstrap=bootstrap,
            n_bootstraps=n_bootstraps, random_state=seed, method=method)
        nz = s[s != 0.0]
        out_rows.append({"protein": prot, "n_sig": int(len(nz)),
                         "corr_json": json.dumps({str(g): float(v) for g, v in nz.items()})})
        if out_path is not None:
            wide = s.to_frame().T
            wide.index = pd.Index([prot], name="protein")
            wide.columns = wide.columns.astype(str)
            safe = _sanitize_filename(prot)
            if per_protein_format == "csv":
                wide.to_csv(out_path / per_protein_subdir / f"{safe}_corr.csv", index=True)
            else:
                _safe_to_parquet_or_csv(wide, out_path / per_protein_subdir / f"{safe}_corr.parquet", index=True)
    return pd.DataFrame(out_rows)


def assemble_full_matrix(per_protein_dir: Union[str, Path], protein_order: Optional[Sequence[str]] = None,
                         square: bool = False) -> pd.DataFrame:
    """Stack per-protein wide files into a proteins x features DataFrame (missing -> 0).
    square=True (protein-protein run) forces a symmetric-labelled protein x protein matrix
    with diagonal 1."""
    per_protein_dir = Path(per_protein_dir)
    files = sorted(list(per_protein_dir.glob("*_corr.parquet")) + list(per_protein_dir.glob("*_corr.csv")))
    if not files:
        raise FileNotFoundError(f"No *_corr.parquet/csv in {per_protein_dir}")
    rows = []
    for f in files:
        df = pd.read_parquet(f) if f.suffix == ".parquet" else pd.read_csv(f, index_col=0)
        if "protein" in df.columns:
            df = df.set_index("protein")
        rows.append(df)
    full = pd.concat(rows, axis=0).astype("float32")
    full = full[~full.index.duplicated(keep="last")]
    if protein_order is not None:
        order = [p for p in protein_order if p in full.index]
        full = full.loc[order]
    full.index.name = None
    full.columns.name = None
    if square:
        prots = list(full.index)
        arr = full.reindex(columns=prots).fillna(0.0).to_numpy(dtype="float32", copy=True)
        np.fill_diagonal(arr, 1.0)
        full = pd.DataFrame(arr, index=pd.Index(prots), columns=pd.Index(prots))
    return full.fillna(0.0)
