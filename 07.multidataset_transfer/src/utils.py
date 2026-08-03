from __future__ import annotations

import hashlib
import itertools
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import joblib
import numpy as np
import pandas as pd
import yaml
from scipy.stats import pearsonr, spearmanr, rankdata
from sklearn.feature_selection import f_regression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import StratifiedKFold, GroupKFold
from sklearn.decomposition import PCA

EPS_DEFAULT = 1e-8


# ----------------------------- general I/O -----------------------------

def load_yaml(path: str | Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def dump_json(obj: Any, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, default=json_default)


def json_default(obj: Any) -> Any:
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, pd.Timestamp):
        return obj.isoformat()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def safe_name(value: Any, max_len: int = 120) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")
    if len(text) <= max_len:
        return text
    digest = hashlib.md5(text.encode("utf-8")).hexdigest()[:10]
    return f"{text[:max_len-11]}_{digest}"


def dataframe_hash(df: pd.DataFrame) -> str:
    h = hashlib.sha256()
    h.update("\n".join(map(str, df.index)).encode())
    h.update("\n".join(map(str, df.columns)).encode())
    h.update(np.ascontiguousarray(df.to_numpy(dtype=np.float32, copy=False)).tobytes())
    return h.hexdigest()


def write_table(df: pd.DataFrame, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".parquet":
        try:
            df.to_parquet(path, index=False)
            return path
        except Exception:
            fallback = path.with_suffix(".csv.gz")
            df.to_csv(fallback, index=False)
            return fallback
    df.to_csv(path, index=False)
    return path


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.exists():
        if path.suffix == ".parquet":
            return pd.read_parquet(path)
        return pd.read_csv(path)
    alt = path.with_suffix(".csv.gz")
    if alt.exists():
        return pd.read_csv(alt)
    raise FileNotFoundError(path)


def load_dataframe(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith(".joblib"):
        obj = joblib.load(path)
    elif suffixes.endswith(".parquet"):
        obj = pd.read_parquet(path)
    elif suffixes.endswith(".csv") or suffixes.endswith(".csv.gz"):
        obj = pd.read_csv(path, index_col=0)
    elif suffixes.endswith(".tsv") or suffixes.endswith(".tsv.gz"):
        obj = pd.read_csv(path, sep="\t", index_col=0)
    else:
        raise ValueError(f"Unsupported table format: {path}")
    if not isinstance(obj, pd.DataFrame):
        raise TypeError(f"Expected DataFrame in {path}, got {type(obj)}")
    return obj


# ----------------------------- Collection A loading -----------------------------

DEFAULT_COLLECTION_A_GROUPS = {
    "TMT": [
        "colon_2019", "renal_2019", "Endomet_val", "breast2020", "pedbrain_2020",
        "lung_2020_smoke", "lung_val", "lung2021", "headneck_2021", "pancrea2021",
        "brain_2021", "liver_val",
    ],
    "iTRAQ": ["breast_2016", "ovary_2016", "stomach_2019"],
    "LabelFree": ["brain", "colon", "prostate", "liver", "lung"],
}


def strip_ensembl_version(value: Any) -> str:
    text = str(value).strip()
    return text.split(".", 1)[0] if text.upper().startswith("ENSG") else text


def looks_like_ensembl(values: Iterable[Any], min_fraction: float = 0.20) -> bool:
    vals = [strip_ensembl_version(v).upper() for v in list(values)[:500]]
    if not vals:
        return False
    return sum(v.startswith("ENSG") for v in vals) / len(vals) >= min_fraction


def read_gene_mapping(path: str | Path | None) -> dict[str, str]:
    """Read a local Ensembl-to-symbol map.

    Accepted Ensembl columns include ensembl_gene_id, ensembl, query, and gene_id.
    Accepted symbol columns include gene_symbol, symbol, external_gene_name, and hgnc_symbol.
    """
    if path is None or not Path(path).exists():
        return {}
    path = Path(path)
    sep = "\t" if ".tsv" in "".join(path.suffixes).lower() else ","
    table = pd.read_csv(path, sep=sep)
    ens_col = next((c for c in ["ensembl_gene_id", "ensembl", "query", "gene_id"] if c in table.columns), None)
    sym_col = next((c for c in ["gene_symbol", "symbol", "external_gene_name", "hgnc_symbol"] if c in table.columns), None)
    if ens_col is None or sym_col is None:
        raise ValueError(
            f"Could not identify Ensembl/symbol columns in {path}. Columns: {list(table.columns)}"
        )
    table = table[[ens_col, sym_col]].dropna()
    mapping = {
        strip_ensembl_version(e).upper(): str(s).strip().upper()
        for e, s in table.itertuples(index=False)
        if str(s).strip()
    }
    return mapping


def query_mygene_mapping(ensembl_ids: Iterable[Any], cache_path: str | Path | None = None) -> dict[str, str]:
    """Query MyGene only when a local mapping is unavailable; save a reusable cache."""
    ids = sorted({strip_ensembl_version(x).upper() for x in ensembl_ids if str(x).upper().startswith("ENSG")})
    if not ids:
        return {}
    try:
        import mygene
    except ImportError as exc:
        raise RuntimeError(
            "Ensembl identifiers were detected but no local mapping was supplied and mygene is not installed. "
            "Provide data/reference/ensembl_to_symbol.tsv or install mygene."
        ) from exc
    mg = mygene.MyGeneInfo()
    result = mg.querymany(
        ids,
        scopes="ensembl.gene",
        fields="symbol",
        species="human",
        as_dataframe=False,
        returnall=False,
        verbose=False,
    )
    mapping: dict[str, str] = {}
    for row in result:
        query = strip_ensembl_version(row.get("query", "")).upper()
        symbol = row.get("symbol")
        if query and symbol and not row.get("notfound", False):
            mapping[query] = str(symbol).strip().upper()
    if cache_path is not None:
        cache_path = Path(cache_path)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            {"ensembl_gene_id": list(mapping.keys()), "gene_symbol": list(mapping.values())}
        ).to_csv(cache_path, sep="\t", index=False)
    return mapping


def map_feature_index_to_symbols(
    matrix_features_by_samples: pd.DataFrame,
    mapping: dict[str, str],
    feature_case: str = "upper",
    drop_unmapped_ensembl: bool = True,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Map a features-by-samples matrix to gene symbols and collapse duplicates by mean."""
    df = matrix_features_by_samples.copy()
    raw_ids = pd.Index([strip_ensembl_version(x) for x in df.index])
    mapped = []
    n_mapped = 0
    n_unmapped_ensembl = 0
    for raw in raw_ids:
        key = str(raw).upper()
        if key.startswith("ENSG"):
            symbol = mapping.get(key)
            if symbol:
                mapped.append(symbol)
                n_mapped += 1
            else:
                mapped.append(np.nan if drop_unmapped_ensembl else key)
                n_unmapped_ensembl += 1
        else:
            mapped.append(str(raw))
    df.index = pd.Index(mapped)
    df = df.loc[df.index.notna()].copy()
    df.index = _normalize_labels(df.index, feature_case)
    before = len(df)
    if df.index.duplicated().any():
        df = df.groupby(level=0, sort=False).mean(numeric_only=True)
    stats = {
        "n_input_features": len(raw_ids),
        "n_mapped_ensembl": n_mapped,
        "n_unmapped_ensembl": n_unmapped_ensembl,
        "n_duplicate_rows_collapsed": before - len(df),
        "n_output_features": len(df),
    }
    return df, stats


def filter_missing_features_samples(
    df_features_by_samples: pd.DataFrame,
    feature_nonmissing_fraction: float = 0.50,
    sample_nonmissing_fraction: float = 0.50,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Match the Collection A notebook rule: retain features/samples with >= threshold observed."""
    if df_features_by_samples.empty:
        return df_features_by_samples, {
            "features_dropped": 0, "samples_dropped": 0,
        }
    feature_keep = df_features_by_samples.notna().mean(axis=1) >= float(feature_nonmissing_fraction)
    after_features = df_features_by_samples.loc[feature_keep].copy()
    sample_keep = after_features.notna().mean(axis=0) >= float(sample_nonmissing_fraction)
    out = after_features.loc[:, sample_keep].copy()
    return out, {
        "features_dropped": int((~feature_keep).sum()),
        "samples_dropped": int((~sample_keep).sum()),
    }


def build_collection_a_platforms(
    collection_dir: str | Path,
    output_dir: str | Path,
    cohort_groups: dict[str, list[str]] | None = None,
    mapping_path: str | Path | None = None,
    allow_mygene_query: bool = False,
    feature_nonmissing_fraction: float = 0.50,
    sample_nonmissing_fraction: float = 0.50,
    map_protein_ids: bool = True,
    feature_case: str = "upper",
    force: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load Collection A exactly once and merge cohorts into TMT/iTRAQ/LabelFree matrices.

    Input cohort folders must contain protein_matrix_*.txt and RNA_matrix_*.txt with
    features on rows and samples on columns. Outputs are samples x features joblibs plus
    sample metadata, matching the downstream transfer pipeline.
    """
    collection_dir = Path(collection_dir)
    output_dir = Path(output_dir)
    cohort_groups = cohort_groups or DEFAULT_COLLECTION_A_GROUPS
    output_dir.mkdir(parents=True, exist_ok=True)
    if not collection_dir.exists():
        raise FileNotFoundError(collection_dir)

    cohort_to_platform = {
        str(cohort).lower(): platform
        for platform, cohorts in cohort_groups.items()
        for cohort in cohorts
    }
    storage = {
        platform: {"rna": [], "protein": [], "metadata": []}
        for platform in cohort_groups
    }
    cohort_rows: list[dict[str, Any]] = []
    all_feature_ids: list[str] = []

    subdirs = sorted([p for p in collection_dir.iterdir() if p.is_dir()])
    for subdir in subdirs:
        cohort = subdir.name
        platform = cohort_to_platform.get(cohort.lower())
        if platform is None:
            cohort_rows.append({"cohort": cohort, "platform": "UNASSIGNED", "status": "skipped_not_in_groups"})
            continue
        prot_files = sorted(subdir.glob("protein_matrix_*.txt"))
        rna_files = sorted(subdir.glob("RNA_matrix_*.txt"))
        if not prot_files or not rna_files:
            cohort_rows.append({"cohort": cohort, "platform": platform, "status": "missing_matrix_file"})
            continue
        prot = pd.read_csv(prot_files[0], sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce")
        rna = pd.read_csv(rna_files[0], sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce")
        prot.columns = prot.columns.astype(str)
        rna.columns = rna.columns.astype(str)
        common = rna.columns.intersection(prot.columns)
        if len(common) == 0:
            cohort_rows.append({"cohort": cohort, "platform": platform, "status": "no_shared_samples"})
            continue
        rna = rna.loc[:, common].copy()
        prot = prot.loc[:, common].copy()
        all_feature_ids.extend(map(str, rna.index))
        if map_protein_ids:
            all_feature_ids.extend(map(str, prot.index))

        safe_cohort = re.sub(r"[^A-Za-z0-9]+", "_", cohort).strip("_").upper()
        new_ids = [f"{safe_cohort}_{i+1:04d}" for i in range(len(common))]
        rna.columns = new_ids
        prot.columns = new_ids
        meta = pd.DataFrame({
            "sample_id": new_ids,
            "original_sample_id": list(common.astype(str)),
            "cohort": cohort,
            "platform": platform,
            "dataset": platform,
        }).set_index("sample_id")
        storage[platform]["rna"].append(rna)
        storage[platform]["protein"].append(prot)
        storage[platform]["metadata"].append(meta)
        cohort_rows.append({
            "cohort": cohort,
            "platform": platform,
            "status": "loaded",
            "n_shared_samples": len(common),
            "n_rna_features": rna.shape[0],
            "n_protein_features": prot.shape[0],
            "rna_file": str(rna_files[0]),
            "protein_file": str(prot_files[0]),
        })

    mapping = read_gene_mapping(mapping_path)
    if looks_like_ensembl(all_feature_ids) and not mapping:
        cache_path = output_dir.parent.parent / "reference" / "ensembl_to_symbol_cache.tsv"
        if cache_path.exists():
            mapping = read_gene_mapping(cache_path)
        elif allow_mygene_query:
            mapping = query_mygene_mapping(all_feature_ids, cache_path=cache_path)
        else:
            raise RuntimeError(
                "Collection A contains Ensembl identifiers, but no mapping is available. "
                "Set MAPPING_PATH to a local TSV or run with ALLOW_MYGENE_QUERY=1 once."
            )

    platform_rows: list[dict[str, Any]] = []
    for platform in cohort_groups:
        platform_dir = output_dir / platform
        platform_dir.mkdir(parents=True, exist_ok=True)
        out_rna = platform_dir / "transcriptomics.joblib"
        out_prot = platform_dir / "proteomics.joblib"
        if out_rna.exists() and out_prot.exists() and not force:
            platform_rows.append({"platform": platform, "status": "skipped_existing"})
            continue
        if not storage[platform]["rna"]:
            platform_rows.append({"platform": platform, "status": "no_loaded_cohorts"})
            continue

        rna_fs = pd.concat(storage[platform]["rna"], axis=1, join="outer")
        prot_fs = pd.concat(storage[platform]["protein"], axis=1, join="outer")
        rna_fs, rna_map_stats = map_feature_index_to_symbols(rna_fs, mapping, feature_case=feature_case)
        prot_fs, prot_map_stats = map_feature_index_to_symbols(
            prot_fs, mapping if map_protein_ids else {}, feature_case=feature_case,
            drop_unmapped_ensembl=map_protein_ids,
        )
        rna_fs, rna_filter = filter_missing_features_samples(
            rna_fs, feature_nonmissing_fraction, sample_nonmissing_fraction
        )
        prot_fs, prot_filter = filter_missing_features_samples(
            prot_fs, feature_nonmissing_fraction, sample_nonmissing_fraction
        )
        common_samples = rna_fs.columns.intersection(prot_fs.columns)
        rna = rna_fs.loc[:, common_samples].T.astype(np.float32)
        prot = prot_fs.loc[:, common_samples].T.astype(np.float32)
        meta = pd.concat(storage[platform]["metadata"], axis=0).reindex(common_samples)

        joblib.dump(rna, out_rna, compress=3)
        joblib.dump(prot, out_prot, compress=3)
        meta.to_csv(platform_dir / "metadata.tsv", sep="\t")
        dump_json(
            {
                "platform": platform,
                "collection_dir": str(collection_dir.resolve()),
                "rna_shape": list(rna.shape),
                "protein_shape": list(prot.shape),
                "rna_mapping": rna_map_stats,
                "protein_mapping": prot_map_stats,
                "rna_filtering": rna_filter,
                "protein_filtering": prot_filter,
                "feature_nonmissing_fraction": feature_nonmissing_fraction,
                "sample_nonmissing_fraction": sample_nonmissing_fraction,
                "mapping_path": str(mapping_path) if mapping_path else None,
            },
            platform_dir / "loader_manifest.json",
        )
        platform_rows.append({
            "platform": platform,
            "status": "saved",
            "n_samples": len(common_samples),
            "n_rna_features": rna.shape[1],
            "n_proteins": prot.shape[1],
            "rna_missing_fraction": float(rna.isna().mean().mean()),
            "protein_missing_fraction": float(prot.isna().mean().mean()),
        })

    cohort_summary = pd.DataFrame(cohort_rows)
    platform_summary = pd.DataFrame(platform_rows)
    cohort_summary.to_csv(output_dir / "collection_a_cohort_summary.tsv", sep="\t", index=False)
    platform_summary.to_csv(output_dir / "collection_a_platform_summary.tsv", sep="\t", index=False)
    return cohort_summary, platform_summary


# ----------------------------- matrix matching -----------------------------

def _normalize_labels(values: Iterable[Any], case: str = "upper") -> pd.Index:
    idx = pd.Index([str(v).strip() for v in values])
    if case == "upper":
        idx = idx.str.upper()
    elif case == "lower":
        idx = idx.str.lower()
    return idx


def collapse_duplicate_columns(df: pd.DataFrame) -> pd.DataFrame:
    if not df.columns.duplicated().any():
        return df
    return df.T.groupby(level=0, sort=False).mean(numeric_only=True).T


def orient_rna_protein(rna: pd.DataFrame, protein: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Choose the pair of orientations that maximizes exact sample-ID overlap."""
    candidates = []
    for r_t in (False, True):
        r = rna.T if r_t else rna
        for p_t in (False, True):
            p = protein.T if p_t else protein
            overlap = len(pd.Index(r.index.astype(str)).intersection(pd.Index(p.index.astype(str))))
            candidates.append((overlap, r_t, p_t, r, p))
    candidates.sort(key=lambda x: x[0], reverse=True)
    overlap, r_t, p_t, rna_o, protein_o = candidates[0]
    if overlap == 0:
        details = [(x[0], x[1], x[2]) for x in candidates]
        raise ValueError(f"No shared sample IDs between RNA and protein axes. Orientation overlaps: {details}")
    info = {
        "rna_transposed": r_t,
        "protein_transposed": p_t,
        "sample_overlap": overlap,
        "candidate_overlaps": [
            {"overlap": x[0], "rna_transposed": x[1], "protein_transposed": x[2]}
            for x in candidates
        ],
    }
    return rna_o.copy(), protein_o.copy(), info


def load_and_align_dataset(
    dataset_name: str,
    rna_path: str | Path,
    protein_path: str | Path,
    metadata_path: str | Path | None = None,
    feature_case: str = "upper",
    prefix_sample_ids: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    rna_raw = load_dataframe(rna_path)
    protein_raw = load_dataframe(protein_path)
    rna, protein, orientation = orient_rna_protein(rna_raw, protein_raw)

    rna.index = rna.index.astype(str)
    protein.index = protein.index.astype(str)
    common_samples = rna.index.intersection(protein.index)
    rna = rna.loc[common_samples].copy()
    protein = protein.loc[common_samples].copy()

    rna.columns = _normalize_labels(rna.columns, feature_case)
    protein.columns = _normalize_labels(protein.columns, feature_case)
    rna = collapse_duplicate_columns(rna)
    protein = collapse_duplicate_columns(protein)

    rna = rna.apply(pd.to_numeric, errors="coerce").astype(np.float32)
    protein = protein.apply(pd.to_numeric, errors="coerce").astype(np.float32)

    original_ids = rna.index.astype(str)
    if prefix_sample_ids:
        new_ids = pd.Index([f"{dataset_name}::{x}" for x in original_ids])
        rna.index = new_ids
        protein.index = new_ids
    else:
        new_ids = original_ids

    meta = pd.DataFrame(index=new_ids)
    meta.index.name = "sample_id"
    meta["dataset"] = dataset_name
    meta["original_sample_id"] = list(original_ids)
    meta["cohort"] = [infer_cohort(x) for x in original_ids]

    if metadata_path:
        metadata_path = Path(metadata_path)
        if metadata_path.exists():
            sep = "\t" if ".tsv" in "".join(metadata_path.suffixes).lower() else ","
            extra = pd.read_csv(metadata_path, sep=sep)
            sample_col = next((c for c in ["sample_id", "sample", "Sample", "DepMap_ID"] if c in extra.columns), None)
            if sample_col is not None:
                extra[sample_col] = extra[sample_col].astype(str)
                extra = extra.set_index(sample_col)
                extra = extra.reindex(original_ids)
                extra.index = new_ids
                for col in extra.columns:
                    if col not in meta.columns:
                        meta[col] = extra[col].values
                for candidate in ["cohort", "lineage", "tissue", "primary_disease"]:
                    if candidate in extra.columns:
                        meta[candidate] = extra[candidate].values

    info = {
        "dataset": dataset_name,
        "rna_path": str(rna_path),
        "protein_path": str(protein_path),
        "rna_raw_shape": list(rna_raw.shape),
        "protein_raw_shape": list(protein_raw.shape),
        "rna_aligned_shape": list(rna.shape),
        "protein_aligned_shape": list(protein.shape),
        "orientation": orientation,
        "n_common_samples": len(common_samples),
        "rna_missing_fraction": float(rna.isna().mean().mean()),
        "protein_missing_fraction": float(protein.isna().mean().mean()),
    }
    return rna, protein, meta, info


def infer_cohort(sample_id: str) -> str:
    text = str(sample_id)
    if "::" in text:
        text = text.split("::", 1)[1]
    # Collection A samples are commonly COHORT_1, COHORT_2, etc.
    return re.sub(r"[_-]?\d+$", "", text)


def feature_intersection(frames: dict[str, pd.DataFrame], datasets: list[str]) -> pd.Index:
    common: pd.Index | None = None
    for ds in datasets:
        cols = pd.Index(frames[ds].columns.astype(str))
        common = cols if common is None else common.intersection(cols)
    return pd.Index([]) if common is None else common.sort_values()


def filter_protein_universe(
    protein_frames: dict[str, pd.DataFrame],
    train_datasets: list[str],
    test_datasets: list[str],
    min_train_fraction: float,
    min_test_fraction: float,
    mode: str = "all_domains_intersection",
) -> pd.Index:
    datasets = train_datasets + test_datasets if mode == "all_domains_intersection" else train_datasets
    common = feature_intersection(protein_frames, datasets)
    keep = []
    for protein in common:
        train_ok = all(protein_frames[d][protein].notna().mean() >= min_train_fraction for d in train_datasets)
        test_ok = all(protein_frames[d][protein].notna().mean() >= min_test_fraction for d in test_datasets)
        if train_ok and test_ok:
            keep.append(protein)
    return pd.Index(keep)


def concatenate_domains(frames: dict[str, pd.DataFrame], datasets: list[str], columns: pd.Index) -> pd.DataFrame:
    parts = [frames[d].reindex(columns=columns) for d in datasets]
    return pd.concat(parts, axis=0)


# ----------------------------- normalization -----------------------------

def domain_balanced_sample_weights(domains: pd.Series | np.ndarray) -> np.ndarray:
    domains = pd.Series(domains).astype(str)
    counts = domains.value_counts()
    weights = domains.map({d: 1.0 / n for d, n in counts.items()}).to_numpy(dtype=float)
    weights *= len(weights) / weights.sum()
    return weights


def weighted_nan_mean_var(X: np.ndarray, weights: np.ndarray, eps: float = EPS_DEFAULT) -> tuple[np.ndarray, np.ndarray]:
    X = np.asarray(X, dtype=float)
    w = np.asarray(weights, dtype=float)[:, None]
    mask = np.isfinite(X)
    denom = np.sum(w * mask, axis=0)
    denom = np.where(denom < eps, 1.0, denom)
    mean = np.nansum(np.where(mask, X, 0.0) * w, axis=0) / denom
    var = np.nansum(np.where(mask, (X - mean) ** 2, 0.0) * w, axis=0) / denom
    return mean.astype(np.float32), var.astype(np.float32)


def balanced_quantile_training_view(
    X: pd.DataFrame,
    domains: pd.Series,
    max_samples_per_domain: int,
    seed: int,
) -> pd.DataFrame:
    """Return an equal-size sample from every training domain for quantile fitting."""
    rng = np.random.default_rng(seed)
    rows = []
    domains = domains.reindex(X.index).astype(str)
    groups = {str(domain): np.asarray(list(idx)) for domain, idx in domains.groupby(domains).groups.items()}
    if not groups:
        return X
    target_n = min(int(max_samples_per_domain), min(len(idx) for idx in groups.values()))
    target_n = max(1, target_n)
    for domain in sorted(groups):
        idx = groups[domain]
        if len(idx) > target_n:
            idx = rng.choice(idx, size=target_n, replace=False)
        rows.append(X.loc[idx])
    return pd.concat(rows, axis=0)


def fit_x_normalizer(
    X: pd.DataFrame,
    domains: pd.Series,
    method: str,
    n_bins: int = 50,
    eps: float = EPS_DEFAULT,
    domain_balanced: bool = True,
    quantile_fit_max_samples_per_domain: int = 1000,
    seed: int = 42,
) -> dict[str, Any]:
    method = method.lower()
    bundle: dict[str, Any] = {
        "method": method,
        "features": list(map(str, X.columns)),
        "n_bins": int(n_bins),
        "eps": float(eps),
        "domain_balanced": bool(domain_balanced),
    }
    if method == "zscore":
        weights = domain_balanced_sample_weights(domains) if domain_balanced else np.ones(len(X), dtype=float)
        mean, var = weighted_nan_mean_var(X.to_numpy(), weights, eps)
        sd = np.sqrt(var)
        sd = np.where(sd < eps, 1.0, sd).astype(np.float32)
        bundle.update({"mean": mean, "sd": sd})
    elif method == "rank":
        pass
    elif method in {"binning", "quantile_bin"}:
        fit_view = balanced_quantile_training_view(X, domains, quantile_fit_max_samples_per_domain, seed) if domain_balanced else X
        quantiles = np.linspace(0, 1, n_bins + 1)[1:-1]
        fit_arr = fit_view.to_numpy(dtype=float)
        edges = np.nanquantile(fit_arr, quantiles, axis=0).astype(np.float32)
        median = np.nanmedian(fit_arr, axis=0).astype(np.float32)
        median = np.where(np.isfinite(median), median, 0.0).astype(np.float32)
        bundle.update({"quantiles": quantiles.astype(np.float32), "edges": edges, "median": median})
    elif method == "none":
        med = np.nanmedian(X.to_numpy(dtype=float), axis=0).astype(np.float32)
        bundle.update({"median": med})
    else:
        raise ValueError(f"Unknown X normalization: {method}")
    return bundle


def transform_x(X: pd.DataFrame, bundle: dict[str, Any]) -> pd.DataFrame:
    features = pd.Index(bundle["features"])
    X = X.reindex(columns=features)
    method = bundle["method"]
    arr = X.to_numpy(dtype=np.float32)
    if method == "zscore":
        mean = np.asarray(bundle["mean"], dtype=np.float32)
        sd = np.asarray(bundle["sd"], dtype=np.float32)
        arr = np.where(np.isfinite(arr), arr, mean)
        arr = (arr - mean) / sd
    elif method == "rank":
        out = np.empty_like(arr, dtype=np.float32)
        for i in range(arr.shape[0]):
            row = arr[i]
            finite = np.isfinite(row)
            out[i] = 0.5
            if finite.any():
                out[i, finite] = rankdata(row[finite], method="average").astype(np.float32) / float(finite.sum())
        arr = out
    elif method in {"binning", "quantile_bin"}:
        edges = np.asarray(bundle["edges"], dtype=np.float32)
        median = np.asarray(bundle["median"], dtype=np.float32)
        n_bins = int(bundle["n_bins"])
        out = np.empty_like(arr, dtype=np.float32)
        for j in range(arr.shape[1]):
            col = arr[:, j]
            col = np.where(np.isfinite(col), col, median[j])
            out[:, j] = np.searchsorted(edges[:, j], col, side="right") / float(max(1, n_bins - 1))
        arr = out
    elif method == "none":
        med = np.asarray(bundle["median"], dtype=np.float32)
        arr = np.where(np.isfinite(arr), arr, med)
    else:
        raise ValueError(method)
    return pd.DataFrame(arr, index=X.index, columns=features, dtype=np.float32)


def fit_y_transformer(
    Y: pd.DataFrame,
    domains: pd.Series,
    method: str,
    eps: float = EPS_DEFAULT,
) -> dict[str, Any]:
    method = method.lower()
    domains = domains.reindex(Y.index).astype(str)
    bundle: dict[str, Any] = {"method": method, "proteins": list(map(str, Y.columns)), "eps": eps, "params": {}}
    if method == "domain_zscore":
        for domain, idx in domains.groupby(domains).groups.items():
            sub = Y.loc[list(idx)]
            mean = sub.mean(axis=0, skipna=True).to_numpy(dtype=np.float32)
            sd = sub.std(axis=0, ddof=0, skipna=True).to_numpy(dtype=np.float32)
            sd = np.where(np.isfinite(sd) & (sd >= eps), sd, 1.0).astype(np.float32)
            mean = np.where(np.isfinite(mean), mean, 0.0).astype(np.float32)
            bundle["params"][str(domain)] = {"mean": mean, "sd": sd}
    elif method == "pooled_train_zscore":
        mean = Y.mean(axis=0, skipna=True).to_numpy(dtype=np.float32)
        sd = Y.std(axis=0, ddof=0, skipna=True).to_numpy(dtype=np.float32)
        sd = np.where(np.isfinite(sd) & (sd >= eps), sd, 1.0).astype(np.float32)
        mean = np.where(np.isfinite(mean), mean, 0.0).astype(np.float32)
        bundle["params"]["__pooled__"] = {"mean": mean, "sd": sd}
    elif method == "none":
        pass
    else:
        raise ValueError(f"Unknown Y normalization: {method}")
    return bundle


def transform_y(Y: pd.DataFrame, domains: pd.Series, bundle: dict[str, Any], strict_unknown: bool = True) -> pd.DataFrame:
    proteins = pd.Index(bundle["proteins"])
    Y = Y.reindex(columns=proteins)
    method = bundle["method"]
    if method == "none":
        return Y.astype(np.float32)
    domains = domains.reindex(Y.index).astype(str)
    raw = Y.to_numpy(dtype=np.float32)
    if method == "pooled_train_zscore":
        p = bundle["params"]["__pooled__"]
        arr = (raw - np.asarray(p["mean"], dtype=np.float32)) / np.asarray(p["sd"], dtype=np.float32)
        return pd.DataFrame(arr.astype(np.float32), index=Y.index, columns=proteins)
    out = np.full(raw.shape, np.nan, dtype=np.float32)
    positions = pd.Series(np.arange(len(Y)), index=Y.index)
    for domain, idx in domains.groupby(domains).groups.items():
        pos = positions.loc[list(idx)].to_numpy(dtype=int)
        if str(domain) not in bundle["params"]:
            if strict_unknown:
                raise KeyError(f"No Y transform parameters for unseen domain {domain}")
            out[pos] = raw[pos]
            continue
        p = bundle["params"][str(domain)]
        out[pos] = (raw[pos] - np.asarray(p["mean"], dtype=np.float32)) / np.asarray(p["sd"], dtype=np.float32)
    return pd.DataFrame(out, index=Y.index, columns=proteins)


def fit_eval_y_zscore(Y: pd.DataFrame, eps: float = EPS_DEFAULT) -> dict[str, Any]:
    mean = Y.mean(axis=0, skipna=True).to_numpy(dtype=np.float32)
    sd = Y.std(axis=0, ddof=0, skipna=True).to_numpy(dtype=np.float32)
    sd = np.where(np.isfinite(sd) & (sd >= eps), sd, 1.0).astype(np.float32)
    mean = np.where(np.isfinite(mean), mean, 0.0).astype(np.float32)
    return {"method": "evaluation_only_test_zscore", "proteins": list(Y.columns), "mean": mean, "sd": sd, "eps": eps}


def apply_eval_y_zscore(Y: pd.DataFrame, bundle: dict[str, Any]) -> pd.DataFrame:
    Y = Y.reindex(columns=bundle["proteins"])
    return pd.DataFrame((Y.to_numpy(dtype=float) - bundle["mean"]) / bundle["sd"], index=Y.index, columns=Y.columns, dtype=np.float32)


# ----------------------------- feature selection -----------------------------

def percentile_ranks(scores: np.ndarray) -> np.ndarray:
    scores = np.asarray(scores, dtype=float)
    scores = np.where(np.isfinite(scores), scores, -np.inf)
    order = np.argsort(np.argsort(scores))
    if len(scores) <= 1:
        return np.ones_like(scores, dtype=float)
    return order / float(len(scores) - 1)


def f_scores(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    try:
        scores, _ = f_regression(X, y, force_finite=True)
    except TypeError:
        scores, _ = f_regression(X, y)
    return np.where(np.isfinite(scores), scores, 0.0)


def spearman_scores(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    y_rank = rankdata(y, method="average")
    y_rank = y_rank - y_rank.mean()
    y_norm = np.sqrt(np.sum(y_rank ** 2))
    scores = np.zeros(X.shape[1], dtype=float)
    if y_norm == 0:
        return scores
    for j in range(X.shape[1]):
        x = X[:, j]
        xr = rankdata(x, method="average")
        xr = xr - xr.mean()
        denom = np.sqrt(np.sum(xr ** 2)) * y_norm
        scores[j] = abs(float(np.dot(xr, y_rank) / denom)) if denom > 0 else 0.0
    return scores


def select_feature_indices(
    X: np.ndarray,
    y: np.ndarray,
    domains: np.ndarray,
    model_type: str,
    k: int,
    score_mode: str = "domain_rank_mean",
) -> tuple[np.ndarray, np.ndarray]:
    model_type = model_type.lower()
    scorer = f_scores if model_type == "enet" else spearman_scores
    k = min(int(k), X.shape[1])
    if score_mode == "pooled":
        scores = scorer(X, y)
    elif score_mode == "domain_rank_mean":
        rank_scores = []
        for domain in pd.unique(domains):
            mask = domains == domain
            if mask.sum() < 4 or np.nanstd(y[mask]) == 0:
                continue
            rank_scores.append(percentile_ranks(scorer(X[mask], y[mask])))
        scores = np.mean(rank_scores, axis=0) if rank_scores else scorer(X, y)
    else:
        raise ValueError(f"Unknown feature_score_mode: {score_mode}")
    idx = np.argpartition(scores, -k)[-k:]
    idx = idx[np.argsort(scores[idx])[::-1]]
    return idx.astype(int), scores[idx].astype(np.float32)


# ----------------------------- CV and metrics -----------------------------

def make_cv_splits(
    domains: np.ndarray,
    cohorts: np.ndarray | None,
    n_splits: int,
    mode: str,
    seed: int,
) -> list[tuple[np.ndarray, np.ndarray]]:
    n_splits = int(n_splits)
    if mode == "group_cohort" and cohorts is not None and len(pd.unique(cohorts)) >= n_splits:
        cv = GroupKFold(n_splits=n_splits)
        return list(cv.split(np.zeros(len(domains)), groups=cohorts))
    counts = pd.Series(domains).value_counts()
    effective = min(n_splits, int(counts.min()))
    if effective < 2:
        raise ValueError("Not enough samples per domain for stratified CV")
    cv = StratifiedKFold(n_splits=effective, shuffle=True, random_state=seed)
    return list(cv.split(np.zeros(len(domains)), domains))


def domain_mean_pearson(y_true: np.ndarray, y_pred: np.ndarray, domains: np.ndarray) -> float:
    vals = []
    for domain in pd.unique(domains):
        mask = (domains == domain) & np.isfinite(y_true) & np.isfinite(y_pred)
        if mask.sum() >= 3 and np.std(y_true[mask]) > 0 and np.std(y_pred[mask]) > 0:
            vals.append(float(pearsonr(y_true[mask], y_pred[mask]).statistic))
    if vals:
        return float(np.mean(vals))
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    if mask.sum() >= 3 and np.std(y_true[mask]) > 0 and np.std(y_pred[mask]) > 0:
        return float(pearsonr(y_true[mask], y_pred[mask]).statistic)
    return -1.0


def score_predictions(y_true: np.ndarray, y_pred: np.ndarray, domains: np.ndarray, metric: str) -> float:
    metric = metric.lower()
    if metric == "domain_mean_pearson":
        return domain_mean_pearson(y_true, y_pred, domains)
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    if mask.sum() == 0:
        return -np.inf
    if metric in {"neg_mse", "domain_mean_mse"}:
        if metric == "domain_mean_mse":
            vals = []
            for domain in pd.unique(domains):
                dm = mask & (domains == domain)
                if dm.sum():
                    vals.append(mean_squared_error(y_true[dm], y_pred[dm]))
            return -float(np.mean(vals)) if vals else -np.inf
        return -float(mean_squared_error(y_true[mask], y_pred[mask]))
    raise ValueError(metric)


def compute_metrics(y_true_raw: np.ndarray, y_true_eval: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    raw_mask = np.isfinite(y_true_raw) & np.isfinite(y_pred)
    eval_mask = np.isfinite(y_true_eval) & np.isfinite(y_pred)
    out: dict[str, float] = {"n_test": int(raw_mask.sum())}
    if raw_mask.sum() >= 3 and np.std(y_true_raw[raw_mask]) > 0 and np.std(y_pred[raw_mask]) > 0:
        out["pearson"] = float(pearsonr(y_true_raw[raw_mask], y_pred[raw_mask]).statistic)
        out["spearman"] = float(spearmanr(y_true_raw[raw_mask], y_pred[raw_mask]).statistic)
    else:
        out["pearson"] = np.nan
        out["spearman"] = np.nan
    if eval_mask.sum() > 0:
        yt, yp = y_true_eval[eval_mask], y_pred[eval_mask]
        mse = float(mean_squared_error(yt, yp))
        out.update({
            "mse_eval_scaled": mse,
            "rmse_eval_scaled": math.sqrt(mse),
            "mae_eval_scaled": float(mean_absolute_error(yt, yp)),
            "r2_eval_scaled": float(r2_score(yt, yp)) if len(yt) >= 2 else np.nan,
            "observed_sd": float(np.std(yt, ddof=0)),
            "predicted_sd": float(np.std(yp, ddof=0)),
        })
        out["sd_ratio"] = out["predicted_sd"] / out["observed_sd"] if out["observed_sd"] > 0 else np.nan
        if len(yt) >= 2 and np.std(yp) > 0:
            slope, intercept = np.polyfit(yp, yt, 1)
            out["calibration_slope"] = float(slope)
            out["calibration_intercept"] = float(intercept)
        else:
            out["calibration_slope"] = np.nan
            out["calibration_intercept"] = np.nan
    else:
        for key in ["mse_eval_scaled", "rmse_eval_scaled", "mae_eval_scaled", "r2_eval_scaled", "observed_sd", "predicted_sd", "sd_ratio", "calibration_slope", "calibration_intercept"]:
            out[key] = np.nan
    return out


def parameter_grid(grid: dict[str, list[Any]], max_combos: int, seed: int) -> list[dict[str, Any]]:
    fixed = {k: v for k, v in grid.items() if not isinstance(v, list)}
    variable = {k: v for k, v in grid.items() if isinstance(v, list)}
    keys = list(variable)
    combos = [dict(zip(keys, vals), **fixed) for vals in itertools.product(*(variable[k] for k in keys))]
    if max_combos and len(combos) > max_combos:
        rng = np.random.default_rng(seed)
        keep = np.sort(rng.choice(len(combos), size=max_combos, replace=False))
        combos = [combos[i] for i in keep]
    return combos


# ----------------------------- diagnostics plotting -----------------------------

def setup_paper_style() -> None:
    import matplotlib as mpl
    import seaborn as sns
    sns.set_theme(context="paper", style="ticks", font_scale=1.15)
    mpl.rcParams.update({
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.titleweight": "bold",
        "font.family": "DejaVu Sans",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def save_figure(fig, base_path: str | Path) -> None:
    base = Path(base_path)
    base.parent.mkdir(parents=True, exist_ok=True)
    for ext in ["pdf", "svg", "png"]:
        fig.savefig(base.with_suffix(f".{ext}"), bbox_inches="tight", dpi=300)


def sample_values_by_domain(df: pd.DataFrame, meta: pd.DataFrame, max_values: int = 50000, seed: int = 42) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    domains = meta.reindex(df.index)["dataset"]
    for domain in pd.unique(domains):
        arr = df.loc[domains.eq(domain)].to_numpy(dtype=float).ravel()
        arr = arr[np.isfinite(arr)]
        if len(arr) > max_values:
            arr = rng.choice(arr, max_values, replace=False)
        rows.append(pd.DataFrame({"value": arr, "dataset": domain}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["value", "dataset"])


def pca_frame(df: pd.DataFrame, meta: pd.DataFrame, max_features: int = 5000, seed: int = 42) -> pd.DataFrame:
    X = df.copy()
    if X.shape[1] > max_features:
        variances = X.var(axis=0).sort_values(ascending=False)
        X = X.loc[:, variances.head(max_features).index]
    arr = X.to_numpy(dtype=float)
    med = np.nanmedian(arr, axis=0)
    arr = np.where(np.isfinite(arr), arr, med)
    arr = arr - arr.mean(axis=0)
    sd = arr.std(axis=0)
    sd = np.where(sd < 1e-8, 1.0, sd)
    arr = arr / sd
    pcs = PCA(n_components=2, random_state=seed).fit_transform(arr)
    out = meta.reindex(X.index).copy()
    out["PC1"] = pcs[:, 0]
    out["PC2"] = pcs[:, 1]
    return out.reset_index()


def create_preparation_figures(
    raw_train_x: pd.DataFrame,
    normalized_train_x: pd.DataFrame,
    raw_train_y: pd.DataFrame,
    normalized_train_y: pd.DataFrame,
    train_meta: pd.DataFrame,
    dataset_summary: pd.DataFrame,
    output_dir: str | Path,
    scenario: str,
    x_norm: str,
) -> list[str]:
    import matplotlib.pyplot as plt
    import seaborn as sns

    setup_paper_style()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    palette = dict(zip(sorted(train_meta["dataset"].unique()), sns.color_palette("colorblind", n_colors=train_meta["dataset"].nunique())))
    saved = []

    # Figure 1: matched dimensions and missingness
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    ds = dataset_summary.sort_values("dataset")
    axes[0].barh(ds["dataset"], ds["n_samples"], color=[palette.get(x, "0.5") for x in ds["dataset"]])
    axes[0].set_xlabel("Matched samples")
    axes[0].set_ylabel("")
    axes[0].set_title("Matched sample coverage")
    width = 0.36
    ypos = np.arange(len(ds))
    axes[1].barh(ypos - width/2, 100*ds["rna_missing_fraction"], height=width, label="RNA")
    axes[1].barh(ypos + width/2, 100*ds["protein_missing_fraction"], height=width, label="Protein")
    axes[1].set_yticks(ypos, ds["dataset"])
    axes[1].set_xlabel("Missing values (%)")
    axes[1].set_title("Missingness before matching")
    axes[1].legend(frameon=False)
    fig.suptitle(f"{scenario}: data substrate", y=1.02, fontweight="bold")
    fig.tight_layout()
    base = output_dir / "fig1_data_substrate"
    save_figure(fig, base); saved.append(str(base)); plt.close(fig)

    # Figure 2: distributions before / after
    raw_long = sample_values_by_domain(raw_train_x, train_meta)
    norm_long = sample_values_by_domain(normalized_train_x, train_meta)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    for domain in sorted(raw_long["dataset"].unique()):
        sns.ecdfplot(raw_long[raw_long.dataset.eq(domain)], x="value", ax=axes[0], label=domain, color=palette[domain])
        sns.ecdfplot(norm_long[norm_long.dataset.eq(domain)], x="value", ax=axes[1], label=domain, color=palette[domain])
    axes[0].set_title("RNA before normalization")
    axes[1].set_title(f"RNA after {x_norm}")
    for ax in axes:
        ax.set_ylabel("Empirical cumulative probability")
        ax.legend(frameon=False, title="Training domain")
    fig.tight_layout()
    base = output_dir / "fig2_rna_ecdf_before_after"
    save_figure(fig, base); saved.append(str(base)); plt.close(fig)

    # Figure 3: PCA before / after
    raw_pca = pca_frame(raw_train_x, train_meta)
    norm_pca = pca_frame(normalized_train_x, train_meta)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))
    for domain in sorted(train_meta["dataset"].unique()):
        for ax, frame, title in [(axes[0], raw_pca, "Before normalization"), (axes[1], norm_pca, f"After {x_norm}")]:
            sub = frame[frame["dataset"].eq(domain)]
            ax.scatter(sub["PC1"], sub["PC2"], s=12, alpha=0.55, label=domain, color=palette[domain], linewidths=0)
            ax.set_title(title)
            ax.set_xlabel("PC1")
            ax.set_ylabel("PC2")
    axes[1].legend(frameon=False, title="Training domain", bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    base = output_dir / "fig3_rna_pca_before_after"
    save_figure(fig, base); saved.append(str(base)); plt.close(fig)

    # Figure 4: protein target distributions
    raw_y_long = sample_values_by_domain(raw_train_y, train_meta)
    norm_y_long = sample_values_by_domain(normalized_train_y, train_meta)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    for domain in sorted(raw_y_long["dataset"].unique()):
        sns.ecdfplot(raw_y_long[raw_y_long.dataset.eq(domain)], x="value", ax=axes[0], label=domain, color=palette[domain])
        sns.ecdfplot(norm_y_long[norm_y_long.dataset.eq(domain)], x="value", ax=axes[1], label=domain, color=palette[domain])
    axes[0].set_title("Protein abundance before harmonization")
    axes[1].set_title("Protein target after training-domain harmonization")
    for ax in axes:
        ax.set_ylabel("Empirical cumulative probability")
        ax.legend(frameon=False, title="Training domain")
    fig.tight_layout()
    base = output_dir / "fig4_protein_ecdf_before_after"
    save_figure(fig, base); saved.append(str(base)); plt.close(fig)

    return saved
