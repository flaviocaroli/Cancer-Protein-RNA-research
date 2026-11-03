"""paper_validation – shared module for the paper external-validation notebooks (06/notebooks 6–10).

Implements the analysis described in the paper (Methods: "Input transformation and external
validation"; Supplementary S5):

* CCLE reference = the 295 development samples of the frozen 80/20 split shared by every primary model.
* External RNA is aligned to the CCLE RNA schema; missing genes/values are completed with CCLE-training
  gene medians (Eq. 7). Z-score models then apply their *saved* CCLE normalizer (Eq. 8). No refitting.
* External proteins stay on their cohort-native scale (identity target).
* Primary models: z-score ElasticNet and XGBoost (k = 5000; folders res_enet_zbig / res_xgb_zbig) and the
  same-gene proxies (RawProxy, ZscoreProxy, LinearProxy) whose parameters are fitted on CCLE training only.
* Sensitivity scenarios (never primary): RNA clipped to CCLE-training quantiles, and per-sample RNA scale
  alignment. TMT1/TMT2 (native protein scale groups) is a reporting stratum only.

Data contracts (consumed by the notebooks):
  cohort bundle  : dict(platform, cohort, xu_shortname, tmt_scale_sensitivity, rna_raw, protein_raw, samples)
  metric table   : one row per cohort x scenario x model x protein (see METRIC_COLUMNS)
  prediction npz : samples, proteins, y, pred_<model> for the primary scenario
  status json    : one per cohort (manifest_hash, metric/prediction paths and sha256, counts)
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import lzma
import os
import pickle
import platform as _platform
import re
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import joblib
import numpy as np
import pandas as pd

PIPELINE_VERSION = "pRNA-external-validation-1.0"
PRIMARY_SCENARIO = "primary_raw"   # notebook 10 filters on this exact label
SCENARIOS = ("primary_raw", "clipped", "scale_aligned")
SCENARIO_TO_VARIANT = {"primary_raw": "raw", "clipped": "clipped", "scale_aligned": "scale_aligned"}
PRIMARY_MODELS = ("LinearProxy", "ElasticNet", "XGBoost")
PROXY_MODELS = ("RawProxy", "ZscoreProxy", "LinearProxy")
MULTIGENE_MODELS = ("ElasticNet", "XGBoost")
ALL_MODELS = ("RawProxy", "ZscoreProxy", "LinearProxy", "ElasticNet", "XGBoost")
MIN_PAIRS = 3

# Primary CCLE-trained models (paper: z-score, k = 5000, fast grid). Folders are looked up under PROJECT_ROOT.
MODEL_SPECS = {
    "ElasticNet": {"directory": "data/res_enet_zbig", "pattern": "{protein}_enet_8020_zscore_fast_k5000.h5"},
    "XGBoost": {"directory": "data/res_xgb_zbig", "pattern": "{protein}_xgb_8020_zscore_fast_k5000.h5"},
}

# The 20 external cohorts (Xu et al. 2022; BioStudies S-BSST733). Cohort labels are the upper-case
# Collection_A folder names that notebook 1 writes as sample prefixes (<COHORT>_<i>).
COHORT_TABLE = pd.DataFrame(
    [
        ("BRAIN_2021", "TMT", "BN_tmt", "brain"),
        ("BREAST2020", "TMT", "BR_tmt", "breast"),
        ("COLON_2019", "TMT", "CO_tmt", "colon"),
        ("ENDOMET_VAL", "TMT", "EC_tmt", "endometrium"),
        ("LIVER_VAL", "TMT", "LV_tmt", "liver"),
        ("LUNG_VAL", "TMT", "LU_1_tmt", "lung"),
        ("LUNG_2020_SMOKE", "TMT", "LU_2_tmt", "lung"),
        ("LUNG2021", "TMT", "LU_3_tmt", "lung"),
        ("PEDBRAIN_2020", "TMT", "PBN_tmt", "brain"),
        ("HEADNECK_2021", "TMT", "HN_tmt", "head_and_neck"),
        ("PANCREA2021", "TMT", "PA_tmt", "pancreas"),
        ("RENAL_2019", "TMT", "RC_tmt", "kidney"),
        ("BREAST_2016", "iTRAQ", "BR_itraq", "breast"),
        ("OVARY_2016", "iTRAQ", "OV_itraq", "ovary"),
        ("STOMACH_2019", "iTRAQ", "GC_itraq", "stomach"),
        ("COLON", "LabelFree", "CO_labelfree", "colon"),
        ("BRAIN", "LabelFree", "BN_labelfree", "brain"),
        ("PROSTATE", "LabelFree", "PR_labelfree", "prostate"),
        ("LIVER", "LabelFree", "LV_labelfree", "liver"),
        ("LUNG", "LabelFree", "LU_labelfree", "lung"),
    ],
    columns=["cohort", "platform", "xu_shortname", "tissue"],
)
EXPECTED_COHORTS_BY_PLATFORM = {
    p: COHORT_TABLE.loc[COHORT_TABLE["platform"].eq(p), "cohort"].tolist() for p in ["TMT", "iTRAQ", "LabelFree"]
}
# TMT cohorts whose native protein values are absolute-abundance-like (median > threshold) form the
# "TMT2" sensitivity stratum; log-ratio cohorts form "TMT1". Reporting stratum only.
TMT2_MEDIAN_THRESHOLD = 5.0

METRIC_COLUMNS = [
    "manifest_hash", "cohort", "platform", "xu_shortname", "tmt_scale_sensitivity", "scenario", "model",
    "protein", "n", "pearson", "spearman", "mse", "rmse", "normalized_rmse", "cosine",
    "calibration_slope", "calibration_intercept",
]


# =============================================================================
# Generic I/O, hashing, provenance
# =============================================================================

def _json_default(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return None if not np.isfinite(obj) else float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (Path,)):
        return str(obj)
    if isinstance(obj, (pd.Timestamp, _dt.datetime, _dt.date)):
        return obj.isoformat()
    if isinstance(obj, (set, tuple)):
        return list(obj)
    return str(obj)


def stable_hash(obj: Any) -> str:
    """sha256 of a canonical JSON serialisation (sorted keys)."""
    text = json.dumps(obj, sort_keys=True, default=_json_default, ensure_ascii=True, separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def stable_int(text: str, modulo: int = 2 ** 31 - 1) -> int:
    return int(hashlib.sha256(str(text).encode("utf-8")).hexdigest()[:12], 16) % modulo


def split_hash(values: Iterable[Any], sorted_values: bool = True) -> str:
    vals = [str(v) for v in values]
    if sorted_values:
        vals = sorted(vals)
    return hashlib.sha256("\n".join(vals).encode("utf-8")).hexdigest()


def file_sha256(path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def safe_name(text: str, max_len: int = 120) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(text)).strip("_")[:max_len] or "item"


def _atomic_replace(tmp: Path, path: Path) -> Path:
    os.replace(tmp, path)
    return path


def atomic_joblib_dump(obj: Any, path, compress: int = 3) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix=path.name, suffix=".tmp")
    os.close(fd)
    joblib.dump(obj, tmp, compress=compress)
    return _atomic_replace(Path(tmp), path)


def atomic_json_dump(obj: Any, path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix=path.name, suffix=".tmp")
    with os.fdopen(fd, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=True, default=_json_default)
    return _atomic_replace(Path(tmp), path)


def write_table(df: pd.DataFrame, base, index: bool = False) -> Path:
    """Write <base>.parquet (fallback <base>.csv). Returns the written path."""
    base = Path(base)
    base.parent.mkdir(parents=True, exist_ok=True)
    stem = base.with_suffix("") if base.suffix in (".parquet", ".csv") else base
    out = df.copy()
    for col in out.columns:  # object columns holding lists/dicts -> JSON for portable parquet/csv
        if out[col].dtype == object and out[col].map(lambda v: isinstance(v, (dict, list, tuple, set, np.ndarray))).any():
            out[col] = out[col].map(lambda v: json.dumps(v, default=_json_default)
                                    if isinstance(v, (dict, list, tuple, set, np.ndarray)) else v)
    try:
        path = stem.with_suffix(".parquet")
        fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix=path.name, suffix=".tmp")
        os.close(fd)
        out.to_parquet(tmp, index=index)
        stale = stem.with_suffix(".csv")
        if stale.exists():
            stale.unlink()
        return _atomic_replace(Path(tmp), path)
    except Exception:
        try:
            os.unlink(tmp)
        except Exception:
            pass
        path = stem.with_suffix(".csv")
        out.to_csv(path, index=index)
        return path


def read_table(path_or_base) -> pd.DataFrame:
    p = Path(path_or_base)
    candidates = [p] if p.suffix in (".parquet", ".csv") else [p.with_suffix(".parquet"), p.with_suffix(".csv")]
    for c in candidates:
        if c.exists():
            return pd.read_parquet(c) if c.suffix == ".parquet" else pd.read_csv(c, low_memory=False)
    raise FileNotFoundError(f"No table at {p} (.parquet/.csv)")


def save_manifest(core: Dict[str, Any], path) -> Dict[str, Any]:
    manifest = dict(core)
    manifest.pop("manifest_hash", None)
    manifest["manifest_hash"] = stable_hash(manifest)
    manifest.setdefault("created_at", _dt.datetime.now().isoformat(timespec="seconds"))
    atomic_json_dump(manifest, path)
    return manifest


def load_manifest(path) -> Dict[str, Any]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Manifest not found: {path}. Run notebook 6 first.")
    return json.loads(path.read_text(encoding="utf-8"))


def notebook_provenance(notebook: str, manifest_hash: str, inputs: Dict[str, Any], outputs: Dict[str, Any]) -> Dict[str, Any]:
    versions = {}
    for mod in ("numpy", "pandas", "scipy", "sklearn", "xgboost", "joblib", "matplotlib"):
        try:
            versions[mod] = __import__(mod).__version__
        except Exception:
            versions[mod] = None
    return {
        "notebook": notebook,
        "pipeline_version": PIPELINE_VERSION,
        "manifest_hash": manifest_hash,
        "created_at": _dt.datetime.now().isoformat(timespec="seconds"),
        "python": sys.version.split()[0],
        "platform": _platform.platform(),
        "package_versions": versions,
        "inputs": inputs,
        "outputs": outputs,
    }


def assert_manifest_column(df: pd.DataFrame, manifest_hash: str, label: str = "table") -> None:
    if "manifest_hash" not in df.columns:
        raise RuntimeError(f"{label}: no manifest_hash column")
    observed = sorted(df["manifest_hash"].dropna().astype(str).unique().tolist())
    if observed != [str(manifest_hash)]:
        raise RuntimeError(f"{label}: manifest_hash mismatch, expected {manifest_hash}, found {observed}")


# =============================================================================
# Project layout
# =============================================================================

def detect_project_root(start: Optional[Path] = None) -> Path:
    """'05. validation' folder: env PAPER_VALIDATION_PROJECT_ROOT, else the first ancestor of the cwd
    holding data/Collection_A or data/processed."""
    env = os.environ.get("PAPER_VALIDATION_PROJECT_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    start = Path(start or Path.cwd()).resolve()
    for base in [start, *start.parents]:
        if (base / "data" / "Collection_A").is_dir() or (base / "data" / "processed").is_dir():
            return base
    for base in [start, *start.parents]:
        if base.name == "05. validation":
            return base
    return start.parent if start.name == "notebooks" else start


def make_output_layout(run_root) -> Dict[str, Path]:
    run_root = Path(run_root)
    layout = {
        "root": run_root,
        "contract": run_root / "contract",
        "cohorts": run_root / "contract" / "cohorts",
        "evaluation": run_root / "evaluation",
        "evaluation_status": run_root / "evaluation" / "status",
        "metric_partitions": run_root / "evaluation" / "metrics",
        "prediction_checkpoints": run_root / "evaluation" / "predictions",
        "tables": run_root / "tables",
        "procedural": run_root / "procedural",
        "xu": run_root / "xu",
    }
    for p in layout.values():
        p.mkdir(parents=True, exist_ok=True)
    return layout


def infer_cohort(sample_id: str, platform: Optional[str] = None) -> str:
    """Cohort label from a notebook-1 sample id '<COHORT>_<i>'."""
    s = str(sample_id).strip()
    s = re.sub(r"_\d+$", "", s)
    return s.upper()


# =============================================================================
# Data loading
# =============================================================================

def _as_float32(df: pd.DataFrame) -> pd.DataFrame:
    out = df.apply(pd.to_numeric, errors="coerce").astype(np.float32)
    out.index = out.index.astype(str)
    out.columns = out.columns.astype(str)
    return out


def load_ccle_matrices(processed_dir) -> Tuple[pd.DataFrame, pd.DataFrame]:
    processed_dir = Path(processed_dir)
    rna = joblib.load(processed_dir / "transcriptomics_filtered.joblib")
    protein = joblib.load(processed_dir / "proteomics_matched.joblib")
    rna, protein = _as_float32(rna), _as_float32(protein)
    rna = rna.loc[:, ~rna.columns.duplicated()]
    protein = protein.loc[:, ~protein.columns.duplicated()]
    common = rna.index.intersection(protein.index, sort=False)
    if len(common) == 0:
        raise RuntimeError("CCLE RNA and protein matrices share no sample identifiers.")
    return rna.loc[common], protein.loc[common]


def load_external_platform(processed_dir, platform: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Notebook-1 output: Collection_A_<platform>_RNA.csv / _Protein.csv (samples x features)."""
    processed_dir = Path(processed_dir)
    if (processed_dir / "processed").is_dir():          # 05. validation/data/processed (notebook 1)
        processed_dir = processed_dir / "processed"
    rna = pd.read_csv(processed_dir / f"Collection_A_{platform}_RNA.csv", index_col=0, low_memory=False)
    protein = pd.read_csv(processed_dir / f"Collection_A_{platform}_Protein.csv", index_col=0, low_memory=False)
    rna, protein = _as_float32(rna), _as_float32(protein)
    rna = rna.loc[~rna.index.duplicated(), ~rna.columns.duplicated()]
    protein = protein.loc[~protein.index.duplicated(), ~protein.columns.duplicated()]
    return rna, protein


# =============================================================================
# Saved CCLE models (80/20 H5 files)
# =============================================================================

def _decode_list(arr) -> List[str]:
    return [x.decode("utf-8") if isinstance(x, (bytes, np.bytes_)) else str(x) for x in arr[:]]


def read_h5_split(path) -> Tuple[List[str], List[str]]:
    import h5py
    with h5py.File(path, "r") as h5:
        group = "split" if "split" in h5 else "samples"
        return _decode_list(h5[f"{group}/train_samples"]), _decode_list(h5[f"{group}/test_samples"])


def load_h5_artifact(path) -> Dict[str, Any]:
    """Model artifact (dict with model, normalizer, selected_features, ...) from either 80/20 H5 format."""
    import h5py
    with h5py.File(path, "r") as h5:
        payload = np.asarray(h5["model/artifact_pickle_xz"][:], dtype=np.uint8)
    artifact = pickle.loads(lzma.decompress(payload.tobytes()))
    artifact.setdefault("normalizer", {"method": "none", "columns": list(map(str, artifact["selected_features"]))})
    return artifact


def apply_saved_normalizer(X: pd.DataFrame, normalizer: Dict[str, Any]) -> pd.DataFrame:
    """Apply a saved 80/20 normalizer (none / zscore) on its own column schema."""
    method = str(normalizer.get("method", "none"))
    cols = [str(c) for c in normalizer.get("columns", X.columns)]
    Xa = X.reindex(columns=cols)
    arr = Xa.to_numpy(dtype=np.float32, copy=True)
    if method == "none":
        pass
    elif method == "zscore":
        arr = (arr - np.asarray(normalizer["mean"], dtype=np.float32)) / np.asarray(normalizer["std"], dtype=np.float32)
    else:
        raise ValueError(f"Normalizer '{method}' is not used by the primary external models (z-score / none only).")
    return pd.DataFrame(arr.astype(np.float32, copy=False), index=Xa.index, columns=cols)


def predict_artifact(artifact: Dict[str, Any], X: pd.DataFrame) -> np.ndarray:
    """ENet and XGB 80/20 artifacts: normalise on the normalizer schema, then take selected features."""
    features = [str(f) for f in artifact["selected_features"]]
    Xn = apply_saved_normalizer(X, artifact["normalizer"])
    Xin = Xn.reindex(columns=features).to_numpy(dtype=np.float32, copy=False)
    return np.asarray(artifact["model"].predict(Xin), dtype=np.float32)


def build_primary_model_inventory(project_root, protein_columns, *, expected_sample_ids=None,
                                  compute_file_hashes: bool = True, verify_all_splits: bool = True,
                                  progress_every: int = 500, model_specs: Optional[Dict] = None):
    """One row per (model family, protein) and the single frozen 80/20 split shared by all of them."""
    project_root = Path(project_root)
    specs = model_specs or MODEL_SPECS
    proteins = [str(p) for p in protein_columns]
    allow_missing = os.environ.get("PAPER_VALIDATION_ALLOW_MISSING_MODELS", "0") == "1"
    rows, splits = [], {}
    for family, spec in specs.items():
        folder = Path(os.environ.get(f"PAPER_VALIDATION_{family.upper()}_DIR", project_root / spec["directory"]))
        if not folder.is_dir():
            raise FileNotFoundError(f"{family} model folder not found: {folder}")
        missing = []
        for i, prot in enumerate(proteins):
            path = folder / spec["pattern"].format(protein=prot)
            if not path.exists():
                missing.append(prot)
                continue
            row = {"display_model": family, "protein": prot, "path": str(path.resolve()),
                   "file_name": path.name, "size_bytes": int(path.stat().st_size)}
            if compute_file_hashes:
                row["sha256"] = file_sha256(path)
            if verify_all_splits or not splits:
                tr, te = read_h5_split(path)
                key = (split_hash(tr), split_hash(te))
                splits.setdefault(key, {"train": tr, "test": te, "path": str(path.resolve()), "n": 0})
                splits[key]["n"] += 1
            rows.append(row)
            if progress_every and (i + 1) % progress_every == 0:
                print(f"  {family}: {i + 1}/{len(proteins)} model files checked")
        if missing and not allow_missing:
            raise FileNotFoundError(f"{family}: {len(missing)} proteins have no model file in {folder} "
                                    f"(first: {missing[:5]}). Set PAPER_VALIDATION_ALLOW_MISSING_MODELS=1 to continue.")
    if len(splits) != 1:
        raise RuntimeError(f"Primary model files disagree on the 80/20 split ({len(splits)} distinct splits).")
    (train_h, test_h), info = next(iter(splits.items()))
    if expected_sample_ids is not None:
        exp = set(map(str, expected_sample_ids))
        extra = (set(info["train"]) | set(info["test"])) - exp
        if extra:
            raise RuntimeError(f"Model split contains {len(extra)} samples absent from the CCLE matrices.")
    inventory = pd.DataFrame(rows).sort_values(["display_model", "protein"]).reset_index(drop=True)
    enet = inventory[inventory["display_model"].eq("ElasticNet")]
    canonical = enet["path"].iloc[0] if len(enet) else inventory["path"].iloc[0]
    split_info = {
        "test_samples": list(info["test"]), "train_samples": list(info["train"]),
        "test_sample_count": len(info["test"]), "train_sample_count": len(info["train"]),
        "test_split_hash_sorted": test_h, "train_split_hash_sorted": train_h,
        "test_split_location": f"{info['path']}::test_samples",
        "train_split_location": f"{info['path']}::train_samples",
        "canonical_model_path": canonical,
        "inventory_hash": stable_hash(inventory.sort_values(["display_model", "protein"]).to_dict("records")),
    }
    return inventory, split_info


# =============================================================================
# CCLE reference and same-gene proxies
# =============================================================================

def fit_ccle_reference(ccle_rna: pd.DataFrame, ccle_protein: pd.DataFrame, test_samples,
                       clip_quantiles=(0.005, 0.995)) -> Dict[str, Any]:
    test = [s for s in map(str, test_samples) if s in ccle_rna.index]
    if len(test) != len(set(map(str, test_samples))):
        raise RuntimeError("Some frozen test samples are missing from the CCLE matrices.")
    train = [s for s in ccle_rna.index.astype(str) if s not in set(test)]
    rna_train = ccle_rna.loc[train].astype(np.float32)
    protein_train = ccle_protein.loc[train].astype(np.float32)
    arr = rna_train.to_numpy(np.float32)
    lo, hi = np.nanquantile(arr, clip_quantiles, axis=0)
    sample_medians = np.nanmedian(arr, axis=1)
    return {
        "train_samples": train, "test_samples": test,
        "test_split_hash_sorted": split_hash(test), "train_split_hash_sorted": split_hash(train),
        "rna_columns": [str(c) for c in ccle_rna.columns], "protein_columns": [str(c) for c in ccle_protein.columns],
        "rna_train": rna_train, "protein_train": protein_train,
        "rna_median": np.nanmedian(arr, axis=0).astype(np.float32),
        "rna_mean": np.nanmean(arr, axis=0).astype(np.float32),
        "rna_std": np.nanstd(arr, axis=0).astype(np.float32),
        "rna_clip_low": lo.astype(np.float32), "rna_clip_high": hi.astype(np.float32),
        "clip_quantiles": tuple(float(q) for q in clip_quantiles),
        "rna_sample_median_reference": float(np.nanmean(sample_medians)),
        "protein_median": np.nanmedian(protein_train.to_numpy(np.float32), axis=0).astype(np.float32),
    }


def fit_linear_proxy_parameters(reference: Dict[str, Any]) -> pd.DataFrame:
    """Per protein, on CCLE training samples only: cognate-RNA mean/SD, protein mean/SD, OLS slope/intercept."""
    rna, prot = reference["rna_train"], reference["protein_train"]
    rna_cols = set(rna.columns)
    rows = []
    for p in reference["protein_columns"]:
        row = {"protein": p, "same_gene_available": False, "n_train": 0, "x_mean": np.nan, "x_std": np.nan,
               "y_mean": np.nan, "y_std": np.nan, "slope": np.nan, "intercept": np.nan}
        if p in rna_cols:
            x = rna[p].to_numpy(np.float64)
            y = prot[p].to_numpy(np.float64)
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() >= MIN_PAIRS and np.std(x[ok]) > 0:
                xm, ym = x[ok].mean(), y[ok].mean()
                xs, ys = x[ok].std(), y[ok].std()
                slope = float(np.sum((x[ok] - xm) * (y[ok] - ym)) / np.sum((x[ok] - xm) ** 2))
                row.update(same_gene_available=True, n_train=int(ok.sum()), x_mean=xm, x_std=xs, y_mean=ym,
                           y_std=ys, slope=slope, intercept=float(ym - slope * xm))
        rows.append(row)
    return pd.DataFrame(rows)


def _same_gene_matrix(rna_aligned: pd.DataFrame, protein_order: Sequence[str], params: pd.DataFrame) -> np.ndarray:
    """Cognate RNA (already median-completed) for every protein; NaN where no same-gene transcript."""
    eligible = params.set_index("protein").reindex(protein_order)["same_gene_available"].fillna(False).to_numpy(bool)
    x = rna_aligned.reindex(columns=list(protein_order)).to_numpy(np.float32, copy=True)
    x[:, ~eligible] = np.nan
    return x


def proxy_predictions(x_same: np.ndarray, protein_order, params: pd.DataFrame) -> Dict[str, np.ndarray]:
    """RawProxy (cognate RNA), ZscoreProxy (y_mean + z * y_std), LinearProxy (slope * x + intercept)."""
    p = params.set_index("protein").reindex(list(protein_order))
    get = lambda c: pd.to_numeric(p[c], errors="coerce").to_numpy(np.float32)
    eligible = p["same_gene_available"].fillna(False).to_numpy(bool)
    x_mean, x_std, y_mean, y_std, slope, intercept = map(get, ["x_mean", "x_std", "y_mean", "y_std", "slope", "intercept"])
    safe = eligible & np.isfinite(x_std) & (x_std > 0) & np.isfinite(y_mean) & np.isfinite(y_std)
    z = np.full_like(x_same, np.nan)
    z[:, safe] = y_mean[safe] + (x_same[:, safe] - x_mean[safe]) / x_std[safe] * y_std[safe]
    lin = x_same * slope[None, :] + intercept[None, :]
    lin[:, ~(eligible & np.isfinite(slope) & np.isfinite(intercept))] = np.nan
    raw = x_same.copy()
    raw[:, ~eligible] = np.nan
    return {"RawProxy": raw, "ZscoreProxy": z.astype(np.float32), "LinearProxy": lin.astype(np.float32)}


# =============================================================================
# External cohort contracts
# =============================================================================

def build_cohort_bundle(platform: str, cohort: str, rna_platform: pd.DataFrame, protein_platform: pd.DataFrame,
                        reference: Dict[str, Any], *, max_sample_missing: float = 0.5,
                        max_rna_feature_missing: float = 0.5, max_protein_feature_missing: float = 0.5):
    """Paper preprocessing per cohort: matched samples; drop samples with >50% missing in either matrix,
    then RNA genes / proteins with >50% missing. Values stay on their native scales."""
    rna_c = rna_platform.loc[[s for s in rna_platform.index if infer_cohort(s, platform) == cohort]]
    prot_c = protein_platform.loc[[s for s in protein_platform.index if infer_cohort(s, platform) == cohort]]
    samples_input = int(len(set(rna_c.index) | set(prot_c.index)))
    common = rna_c.index.intersection(prot_c.index, sort=False)
    rna_c, prot_c = rna_c.loc[common], prot_c.loc[common]
    rna_c = rna_c.dropna(axis=1, how="all")
    prot_c = prot_c.dropna(axis=1, how="all")
    keep = (rna_c.isna().mean(axis=1) <= max_sample_missing) & (prot_c.isna().mean(axis=1) <= max_sample_missing)
    rna_c, prot_c = rna_c.loc[keep], prot_c.loc[keep]
    rna_c = rna_c.loc[:, rna_c.isna().mean(axis=0) <= max_rna_feature_missing]
    prot_c = prot_c.loc[:, prot_c.isna().mean(axis=0) <= max_protein_feature_missing]
    if platform == "TMT":
        med = float(np.nanmedian(prot_c.to_numpy(np.float32))) if prot_c.size else np.nan
        scale_group = "TMT2" if np.isfinite(med) and med > TMT2_MEDIAN_THRESHOLD else "TMT1"
    else:
        scale_group = "not_applicable"
    meta = COHORT_TABLE.set_index("cohort").loc[cohort]
    bundle = {
        "platform": platform, "cohort": cohort, "xu_shortname": str(meta["xu_shortname"]),
        "tissue": str(meta["tissue"]), "tmt_scale_sensitivity": scale_group,
        "samples": [str(s) for s in rna_c.index], "rna_raw": rna_c.astype(np.float32),
        "protein_raw": prot_c.astype(np.float32), "pipeline_version": PIPELINE_VERSION,
    }
    rna_overlap = len(set(rna_c.columns) & set(reference["rna_columns"]))
    prot_overlap = len(set(prot_c.columns) & set(reference["protein_columns"]))
    audit = {
        "platform": platform, "cohort": cohort, "xu_shortname": bundle["xu_shortname"],
        "tmt_scale_sensitivity": scale_group, "samples_input": samples_input,
        "samples_matched": int(len(common)), "samples_retained": int(len(rna_c)),
        "rna_features_retained": int(rna_c.shape[1]), "proteins_retained": int(prot_c.shape[1]),
        "rna_overlap_ccle": int(rna_overlap), "protein_overlap_ccle": int(prot_overlap),
        "rna_missing_fraction": float(rna_c.isna().to_numpy().mean()) if rna_c.size else np.nan,
        "protein_missing_fraction": float(prot_c.isna().to_numpy().mean()) if prot_c.size else np.nan,
        "protein_native_median": float(np.nanmedian(prot_c.to_numpy(np.float32))) if prot_c.size else np.nan,
    }
    return bundle, audit


def prepare_rna_variants(bundle: Dict[str, Any], reference: Dict[str, Any]):
    """'raw' = primary model input before the saved normalizer: CCLE RNA schema, missing completed with
    CCLE-training medians (Eq. 7). 'clipped' and 'scale_aligned' are sensitivity variants."""
    cols = reference["rna_columns"]
    med = np.asarray(reference["rna_median"], dtype=np.float32)
    rna = bundle["rna_raw"].reindex(columns=cols)
    arr = rna.to_numpy(np.float32, copy=True)
    missing = ~np.isfinite(arr)
    arr = np.where(missing, med[None, :], arr).astype(np.float32)
    clipped = np.clip(arr, reference["rna_clip_low"][None, :], reference["rna_clip_high"][None, :]).astype(np.float32)
    shift = reference["rna_sample_median_reference"] - np.nanmedian(arr, axis=1)
    aligned = (arr + shift[:, None]).astype(np.float32)
    idx = rna.index
    variants = {
        "raw": pd.DataFrame(arr, index=idx, columns=cols),
        "clipped": pd.DataFrame(clipped, index=idx, columns=cols),
        "scale_aligned": pd.DataFrame(aligned, index=idx, columns=cols),
    }
    observed_genes = int(len(set(bundle["rna_raw"].columns) & set(cols)))
    obs = ~missing
    n_obs = max(1, int(obs.sum()))
    audit = pd.DataFrame([{
        "cohort": bundle["cohort"], "platform": bundle["platform"],
        "tmt_scale_sensitivity": bundle.get("tmt_scale_sensitivity", "not_applicable"),
        "samples": int(len(idx)), "ccle_rna_features": int(len(cols)), "observed_ccle_genes": observed_genes,
        "genes_completed_by_median": int(len(cols) - observed_genes),
        "completed_values_pct": float(100.0 * missing.mean()) if missing.size else np.nan,
        "clipped_observed_values_pct": float(100.0 * np.sum(obs & ~np.isclose(clipped, arr)) / n_obs),
        "scale_aligned_observed_values_changed_pct": float(100.0 * np.sum(obs & ~np.isclose(aligned, arr)) / n_obs),
        "median_sample_shift": float(np.nanmedian(shift)) if len(shift) else np.nan,
    }])
    return variants, audit


def prepare_protein_targets(bundle: Dict[str, Any], reference: Dict[str, Any]) -> Dict[str, pd.DataFrame]:
    """'raw' = primary identity target (cohort-native scale, CCLE protein schema, NaN kept).
    'harmonized' = per-protein cohort-median centering (diagnostic only, never scored as primary)."""
    raw = bundle["protein_raw"].reindex(columns=reference["protein_columns"]).astype(np.float32)
    harmonized = raw - raw.median(axis=0, skipna=True)
    return {"raw": raw, "harmonized": harmonized.astype(np.float32)}


# =============================================================================
# Metrics
# =============================================================================

def metric_frame(y, pred, protein_order, compute_spearman: bool = True) -> pd.DataFrame:
    """Protein-wise metrics over finite (y, pred) pairs. y, pred: samples x proteins arrays."""
    y = np.asarray(y, dtype=np.float64)
    p = np.asarray(pred, dtype=np.float64)
    ok = np.isfinite(y) & np.isfinite(p)
    n = ok.sum(0).astype(float)
    yz, pz = np.where(ok, y, 0.0), np.where(ok, p, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        ym, pm = yz.sum(0) / n, pz.sum(0) / n
        yc, pc = np.where(ok, y - ym, 0.0), np.where(ok, p - pm, 0.0)
        syy, spp, syp = (yc ** 2).sum(0), (pc ** 2).sum(0), (yc * pc).sum(0)
        pearson = syp / np.sqrt(syy * spp)
        mse = ((yz - pz) ** 2).sum(0) / n
        cosine = (yz * pz).sum(0) / np.sqrt((yz ** 2).sum(0) * (pz ** 2).sum(0))
        slope = syp / spp
        intercept = ym - slope * pm
        sd_y = np.sqrt(syy / n)
        rmse = np.sqrt(mse)
        nrmse = rmse / sd_y
    bad = n < MIN_PAIRS
    pearson[bad | ~(syy > 1e-12) | ~(spp > 1e-12)] = np.nan
    slope[bad | ~(spp > 1e-12)] = np.nan
    intercept[~np.isfinite(slope)] = np.nan
    for arr in (mse, rmse, cosine, nrmse):
        arr[bad] = np.nan
    nrmse[~(sd_y > 1e-12)] = np.nan
    spearman = np.full(len(n), np.nan)
    if compute_spearman:
        for j in np.where(~bad)[0]:
            m = ok[:, j]
            ry = pd.Series(y[m, j]).rank().to_numpy()
            rp = pd.Series(p[m, j]).rank().to_numpy()
            if ry.std() > 0 and rp.std() > 0:
                spearman[j] = float(np.corrcoef(ry, rp)[0, 1])
    return pd.DataFrame({
        "protein": [str(x) for x in protein_order], "n": n.astype(int),
        "pearson": np.clip(pearson, -1, 1), "spearman": spearman, "mse": mse, "rmse": rmse,
        "normalized_rmse": nrmse, "cosine": np.clip(cosine, -1, 1),
        "calibration_slope": slope, "calibration_intercept": intercept,
    })


def bootstrap_stat_ci(values, statistic: str = "median", repetitions: int = 2000, seed: int = 0, alpha: float = 0.05):
    v = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return np.nan, np.nan, np.nan
    fn = np.median if statistic == "median" else np.mean
    point = float(fn(v))
    if len(v) == 1 or repetitions <= 0:
        return point, point, point
    rng = np.random.default_rng(int(seed))
    idx = rng.integers(0, len(v), size=(int(repetitions), len(v)))
    boots = fn(v[idx], axis=1)
    return point, float(np.quantile(boots, alpha / 2)), float(np.quantile(boots, 1 - alpha / 2))


def matched_universe(primary: pd.DataFrame, scenario: str = PRIMARY_SCENARIO) -> pd.DataFrame:
    """Per cohort x protein: measurable, proxy-/multigene-eligible and matched across PRIMARY_MODELS."""
    df = primary[primary["scenario"].eq(scenario)].copy()
    df["valid"] = pd.to_numeric(df["pearson"], errors="coerce").notna()
    keys = ["cohort", "platform", "xu_shortname", "tmt_scale_sensitivity", "protein"]
    wide = df.pivot_table(index=keys, columns="model", values="valid", aggfunc="max", fill_value=False).reset_index()
    for m in ALL_MODELS:
        if m not in wide:
            wide[m] = False
    meas = df.groupby(keys, observed=True)["n"].max().reset_index().rename(columns={"n": "max_pairs"})
    wide = wide.merge(meas, on=keys, how="left")
    wide["measurable_external"] = pd.to_numeric(wide["max_pairs"], errors="coerce").fillna(0) >= MIN_PAIRS
    wide["proxy_eligible"] = wide["LinearProxy"].astype(bool)
    wide["multigene_eligible"] = wide["ElasticNet"].astype(bool) & wide["XGBoost"].astype(bool)
    wide["matched_primary"] = wide[list(PRIMARY_MODELS)].astype(bool).all(axis=1)
    wide["scenario"] = scenario
    wide.columns.name = None
    return wide[keys + ["scenario", "max_pairs", "measurable_external", "proxy_eligible",
                        "multigene_eligible", "matched_primary"]]


def summarize_cohort_metrics(primary: pd.DataFrame, universe: pd.DataFrame, scenario: str = PRIMARY_SCENARIO,
                             repetitions: int = 2000, seed: int = 2026) -> pd.DataFrame:
    matched = set(universe.loc[universe["matched_primary"], "protein"].astype(str))
    df = primary[primary["scenario"].eq(scenario) & primary["protein"].astype(str).isin(matched)]
    rows = []
    for model, g in df.groupby("model", observed=True):
        row = {k: str(g[k].iloc[0]) for k in ["cohort", "platform", "xu_shortname", "tmt_scale_sensitivity"]}
        row.update(scenario=scenario, model=str(model), matched_proteins=int(len(matched)))
        for metric in ["pearson", "spearman", "cosine", "mse", "normalized_rmse"]:
            vals = pd.to_numeric(g[metric], errors="coerce")
            point, lo, hi = bootstrap_stat_ci(vals, "median", repetitions,
                                              seed + stable_int(f"{row['cohort']}|{model}|{metric}"))
            row[f"mean_{metric}"] = float(vals.mean())
            row[f"median_{metric}"] = point
            row[f"median_{metric}_ci_low"] = lo
            row[f"median_{metric}_ci_high"] = hi
            row[f"finite_{metric}_proteins"] = int(vals.notna().sum())
        rows.append(row)
    return pd.DataFrame(rows)


# =============================================================================
# Notebook 7 – memory-safe inference and evaluation
# =============================================================================

def _cohort_meta(bundle):
    return {k: bundle[k] for k in ["cohort", "platform", "xu_shortname", "tmt_scale_sensitivity"]}


def _work_array(work_dir: Path, family: str, cohort: str, scenario: str, shape, create: bool):
    from numpy.lib.format import open_memmap
    path = work_dir / f"{family}__{safe_name(cohort)}__{scenario}.npy"
    if path.exists():
        return open_memmap(path, mode="r+")
    if not create:
        return None
    arr = open_memmap(path, mode="w+", dtype=np.float32, shape=shape)
    arr[:] = np.nan
    return arr


def _verify_normalizer(artifact, reference, rtol=1e-3, atol=1e-4) -> bool:
    norm = artifact.get("normalizer", {})
    if norm.get("method") != "zscore":
        return True
    cols = [str(c) for c in norm["columns"]]
    pos = pd.Index(reference["rna_columns"]).get_indexer(cols)
    if (pos < 0).any():
        return False
    return bool(np.allclose(np.asarray(norm["mean"]), reference["rna_mean"][pos], rtol=rtol, atol=atol)
                and np.allclose(np.asarray(norm["std"]), np.where(reference["rna_std"][pos] > 0,
                                                                   reference["rna_std"][pos], 1.0), rtol=rtol, atol=atol))


def evaluate_all_cohorts_memory_safe(cohort_bundle_paths: Dict[str, str], reference: Dict[str, Any],
                                     proxy_parameters: pd.DataFrame, model_inventory: pd.DataFrame, *,
                                     manifest_hash_value: str, prediction_dir, metric_dir, status_dir, work_dir,
                                     compute_spearman_primary: bool = True, compressed_npz: bool = True,
                                     verify_every_normalizer: bool = True, sample_batch_size: int = 512,
                                     checkpoint_every: int = 50, force: bool = False,
                                     cleanup_work_on_success: bool = False, model_families=("ElasticNet", "XGBoost"),
                                     show_progress: bool = True, finalize: bool = True) -> pd.DataFrame:
    """Phase 1 (per family): each saved model is loaded once and predicts every cohort and scenario;
    predictions go to resumable on-disk arrays. Phase 2 (finalize=True): proxies + metrics per cohort,
    prediction checkpoint (primary scenario) and a status JSON per cohort."""
    import shutil
    prediction_dir, metric_dir, status_dir, work_dir = map(Path, (prediction_dir, metric_dir, status_dir, work_dir))
    if force and work_dir.exists():
        shutil.rmtree(work_dir)
    for d in (prediction_dir, metric_dir, status_dir, work_dir):
        d.mkdir(parents=True, exist_ok=True)
    state_path = work_dir / "state.json"
    state = json.loads(state_path.read_text()) if state_path.exists() else {}
    if state.get("manifest_hash") not in (None, manifest_hash_value):
        raise RuntimeError("Work directory belongs to a different manifest; use force=True.")
    state.setdefault("manifest_hash", manifest_hash_value)
    state.setdefault("done", {})
    state.setdefault("progress", {})
    state.setdefault("failures", {})
    protein_order = list(reference["protein_columns"])
    ppos = {p: i for i, p in enumerate(protein_order)}
    cohorts = list(cohort_bundle_paths)

    # Pre-compute normalised-input variants once per cohort (memory: one cohort at a time).
    def _variants(cohort):
        bundle = joblib.load(cohort_bundle_paths[cohort])
        variants, _ = prepare_rna_variants(bundle, reference)
        return bundle, variants

    try:
        from tqdm.auto import tqdm
    except Exception:  # pragma: no cover
        tqdm = None

    for family in model_families:
        inv = model_inventory[model_inventory["display_model"].eq(family)].sort_values("protein")
        done = set(state["done"].get(family, []))
        todo = [r for r in inv.itertuples(index=False) if r.protein not in done]
        if not todo:
            continue
        cohort_inputs = {c: _variants(c) for c in cohorts}
        arrays = {(c, s): _work_array(work_dir, family, c, s, (len(cohort_inputs[c][1]["raw"]), len(protein_order)), True)
                  for c in cohorts for s in SCENARIOS}
        it = tqdm(todo, desc=family, unit="protein") if (tqdm and show_progress) else todo
        since = 0
        failures = state["failures"].setdefault(family, {})
        for row in it:
            try:
                art = load_h5_artifact(row.path)
                if verify_every_normalizer and not _verify_normalizer(art, reference):
                    raise RuntimeError("saved normalizer does not match CCLE-training statistics")
                j = ppos[row.protein]
                for c in cohorts:
                    variants = cohort_inputs[c][1]
                    for s in SCENARIOS:
                        X = variants[SCENARIO_TO_VARIANT[s]]
                        out = arrays[(c, s)]
                        for start in range(0, len(X), max(1, int(sample_batch_size))):
                            out[start:start + sample_batch_size, j] = predict_artifact(art, X.iloc[start:start + sample_batch_size])
            except Exception as exc:
                failures[row.protein] = f"{type(exc).__name__}: {exc}"
            done.add(row.protein)
            since += 1
            if since >= max(1, int(checkpoint_every)):
                for a in arrays.values():
                    a.flush()
                state["done"][family] = sorted(done)
                state["progress"][family] = len(done)
                atomic_json_dump(state, state_path)
                since = 0
        for a in arrays.values():
            a.flush()
        state["done"][family] = sorted(done)
        state["progress"][family] = len(done)
        atomic_json_dump(state, state_path)
        del arrays, cohort_inputs

    progress = pd.DataFrame([{
        "display_model": f, "completed_proteins": int(state["progress"].get(f, 0)),
        "total_proteins": int(model_inventory["display_model"].eq(f).sum()),
        "failures": int(len(state["failures"].get(f, {}))),
    } for f in MULTIGENE_MODELS])
    if not finalize:
        return progress

    for f in MULTIGENE_MODELS:
        total = int(model_inventory["display_model"].eq(f).sum())
        if int(state["progress"].get(f, 0)) < total:
            raise RuntimeError(f"{f} predictions incomplete ({state['progress'].get(f, 0)}/{total}); run its cell first.")

    rows = []
    for c in cohorts:
        bundle, variants = _variants(c)
        meta = _cohort_meta(bundle)
        y = prepare_protein_targets(bundle, reference)["raw"].to_numpy(np.float32)
        frames, primary_preds = [], {}
        for s in SCENARIOS:
            X = variants[SCENARIO_TO_VARIANT[s]]
            preds = proxy_predictions(_same_gene_matrix(X, protein_order, proxy_parameters), protein_order, proxy_parameters)
            for f in MULTIGENE_MODELS:
                arr = _work_array(work_dir, f, c, s, None, False)
                preds[f] = np.asarray(arr, dtype=np.float32) if arr is not None else np.full_like(y, np.nan)
            for model, pred in preds.items():
                fr = metric_frame(y, pred, protein_order,
                                  compute_spearman=compute_spearman_primary and s == PRIMARY_SCENARIO)
                fr.insert(0, "model", model)
                fr.insert(0, "scenario", s)
                for k, v in reversed(list(meta.items())):
                    fr.insert(0, k, v)
                fr.insert(0, "manifest_hash", manifest_hash_value)
                frames.append(fr[METRIC_COLUMNS])
            if s == PRIMARY_SCENARIO:
                primary_preds = preds
        metrics = pd.concat(frames, ignore_index=True)
        metric_path = write_table(metrics, metric_dir / safe_name(c), index=False)
        pred_path = prediction_dir / f"{safe_name(c)}.npz"
        saver = np.savez_compressed if compressed_npz else np.savez
        saver(pred_path, samples=np.asarray(bundle["samples"], dtype=str), proteins=np.asarray(protein_order, dtype=str),
              y=y, **{f"pred_{m}": v.astype(np.float32) for m, v in primary_preds.items()})
        fails = sum(1 for f in MULTIGENE_MODELS for p in state["failures"].get(f, {}))
        status = {
            **meta, "manifest_hash": manifest_hash_value, "pipeline_version": PIPELINE_VERSION,
            "samples": int(len(bundle["samples"])), "proteins": int(len(protein_order)),
            "models": list(ALL_MODELS), "scenarios": list(SCENARIOS), "model_failures": int(fails),
            "metric_path": str(Path(metric_path).resolve()), "metric_sha256": file_sha256(metric_path),
            "prediction_path": str(pred_path.resolve()), "prediction_sha256": file_sha256(pred_path),
            "completed_at": _dt.datetime.now().isoformat(timespec="seconds"),
        }
        atomic_json_dump(status, status_dir / f"{safe_name(c)}.json")
        rows.append(status)
        print(f"  finalized {c}: {status['samples']} samples")
    if cleanup_work_on_success:
        shutil.rmtree(work_dir, ignore_errors=True)
    return validate_evaluation_completion(status_dir, manifest_hash_value, cohorts)


def validate_evaluation_completion(status_dir, manifest_hash: str, cohorts: Sequence[str]) -> pd.DataFrame:
    status_dir = Path(status_dir)
    rows, problems = [], []
    for c in cohorts:
        p = status_dir / f"{safe_name(c)}.json"
        if not p.exists():
            problems.append(f"{c}: no status file")
            continue
        st = json.loads(p.read_text(encoding="utf-8"))
        if str(st.get("manifest_hash")) != str(manifest_hash):
            problems.append(f"{c}: manifest {st.get('manifest_hash')} != {manifest_hash}")
        for kind in ("metric", "prediction"):
            fp = Path(st.get(f"{kind}_path", ""))
            if not fp.exists():
                problems.append(f"{c}: missing {kind} file {fp}")
            elif file_sha256(fp) != st.get(f"{kind}_sha256"):
                problems.append(f"{c}: {kind} file changed since notebook 7 ({fp})")
        row = dict(st)
        row["models"] = json.dumps(st.get("models", []))
        row["scenarios"] = json.dumps(st.get("scenarios", []))
        rows.append(row)
    if problems:
        raise RuntimeError("Evaluation incomplete or inconsistent:\n  " + "\n  ".join(problems))
    return pd.DataFrame(rows)


def pooled_platform_metrics(prediction_paths: Dict[str, str], cohort_table: pd.DataFrame, manifest_hash: str,
                            scenario: str = PRIMARY_SCENARIO) -> pd.DataFrame:
    """Secondary analysis: protein-wise metrics after pooling samples of all cohorts of a platform."""
    platform_of = cohort_table.set_index("cohort")["platform"].astype(str).to_dict()
    frames = []
    for platform in ["TMT", "iTRAQ", "LabelFree"]:
        ys, preds, proteins = [], {m: [] for m in PRIMARY_MODELS}, None
        for c, path in prediction_paths.items():
            if platform_of.get(c) != platform:
                continue
            z = np.load(path, allow_pickle=False)
            proteins = z["proteins"].tolist()
            ys.append(z["y"])
            for m in PRIMARY_MODELS:
                preds[m].append(z[f"pred_{m}"])
        if not ys:
            continue
        y = np.concatenate(ys, axis=0)
        for m in PRIMARY_MODELS:
            fr = metric_frame(y, np.concatenate(preds[m], axis=0), proteins, compute_spearman=False)
            fr.insert(0, "model", m)
            fr.insert(0, "scenario", scenario)
            fr.insert(0, "platform", platform)
            fr.insert(0, "manifest_hash", manifest_hash)
            frames.append(fr)
    return pd.concat(frames, ignore_index=True)


def proxy_metrics_by_tissue_protein(metric_dir, manifest_hash: str, scenario: str = PRIMARY_SCENARIO) -> pd.DataFrame:
    """Compact same-gene proxy table read by notebook 8 (RawProxy / ZscoreProxy / LinearProxy, one row per
    cohort x model x protein, CCLE-trained parameters, primary scenario)."""
    frames = []
    for path in sorted(Path(metric_dir).glob("*.parquet")) + sorted(Path(metric_dir).glob("*.csv")):
        m = read_table(path)
        m = m[m["scenario"].eq(scenario) & m["model"].isin(PROXY_MODELS)]
        frames.append(m)
    if not frames:
        raise FileNotFoundError(f"No metric partitions in {metric_dir}; run notebook 7 first.")
    df = pd.concat(frames, ignore_index=True)
    out = pd.DataFrame({
        "dataset": df["platform"], "tissue": df["cohort"], "ShortName": df["xu_shortname"],
        "display_model": df["model"], "proxy_variant": df["model"], "protein": df["protein"],
        "pearson": df["pearson"], "cosine": df["cosine"], "mse": df["mse"],
        "run_label": "ccle_trained_external_primary", "proxy_pipeline_version": PIPELINE_VERSION,
    })
    out["manifest_hash"] = manifest_hash
    return out
