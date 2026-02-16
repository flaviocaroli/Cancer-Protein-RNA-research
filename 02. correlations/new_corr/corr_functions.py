import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import numpy as np
import pandas as pd
from scipy.stats import pearsonr


def _sanitize_filename(name: str) -> str:
    """Make a string safe for filenames across most filesystems."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_")


def _safe_to_parquet_or_csv(df: pd.DataFrame, path: Path, *, index: bool = False) -> None:
    """Try parquet first; fallback to csv if parquet engine is unavailable."""
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        df.to_parquet(path, compression="snappy", index=index)
    except Exception:
        df.to_csv(path.with_suffix(".csv"), index=index)


def compute_correlation_for_protein(
    transcriptomics: pd.DataFrame,   # samples x genes
    proteomics: pd.DataFrame,        # samples x proteins
    protein_name: str,
    p_threshold: float = 0.05,
    bootstrap: bool = True,
    n_bootstraps: int = 100,
    random_state: Optional[int] = None,
) -> pd.Series:
    """
    Compute Pearson correlations between one protein and all genes.

    For each gene:
      - bootstrap=False: pearsonr(protein, gene), keep corr if p < p_threshold else 0.
      - bootstrap=True : bootstrap 95% CI + median; keep median if CI excludes 0 else 0.

    Returns:
      pd.Series index=genes, value=correlation (or 0 if not significant)
    """
    if protein_name not in proteomics.columns:
        raise ValueError(f"Protein '{protein_name}' not found in proteomics.")

    rng = np.random.default_rng(random_state)
    correlations = pd.Series(index=transcriptomics.columns, dtype=float)

    prot_series = proteomics[protein_name]

    for gene in transcriptomics.columns:
        gene_series = transcriptomics[gene]
        common_samples = prot_series.dropna().index.intersection(gene_series.dropna().index)

        if len(common_samples) <= 1:
            correlations[gene] = 0.0
            continue

        data_prot = prot_series.loc[common_samples].to_numpy(dtype=float, copy=False)
        data_gene = gene_series.loc[common_samples].to_numpy(dtype=float, copy=False)

        # pearsonr is undefined for constant vectors
        if np.std(data_prot) == 0.0 or np.std(data_gene) == 0.0:
            correlations[gene] = 0.0
            continue

        if not bootstrap:
            corr_val, p_val = pearsonr(data_prot, data_gene)
            if np.isnan(corr_val) or np.isnan(p_val):
                correlations[gene] = 0.0
            else:
                correlations[gene] = float(corr_val) if float(p_val) < float(p_threshold) else 0.0
            continue

        # bootstrap mode
        boot_corrs: List[float] = []
        n = len(common_samples)

        for _ in range(int(n_bootstraps)):
            idx = rng.integers(0, n, size=n)
            xp = data_prot[idx]
            xg = data_gene[idx]

            if np.std(xp) == 0.0 or np.std(xg) == 0.0:
                continue

            c, _ = pearsonr(xp, xg)
            if not np.isnan(c):
                boot_corrs.append(float(c))

        if len(boot_corrs) < max(10, int(n_bootstraps) // 10):
            correlations[gene] = 0.0
            continue

        ci_lower, ci_upper = np.percentile(boot_corrs, [2.5, 97.5])
        median_boot = float(np.median(boot_corrs))

        correlations[gene] = median_boot if (ci_lower > 0.0 or ci_upper < 0.0) else 0.0

    return correlations


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
    per_protein_format: str = "parquet",  # parquet|csv
) -> pd.DataFrame:
    """
    Run correlations for a list of proteins and (optionally) save per-protein wide files.

    If output_dir is provided, writes:
      - per-protein wide row (genes as columns) to:
          <output_dir>/<per_protein_subdir>/<protein>_corr.(parquet|csv)

    Returns:
      compact batch DataFrame (1 row per protein):
        - protein
        - n_sig
        - corr_json  (JSON dict gene->corr for non-zero entries)
        - error (optional)
    """
    out_rows: List[Dict[str, object]] = []

    out_path = Path(output_dir) if output_dir is not None else None
    if out_path is not None:
        (out_path / per_protein_subdir).mkdir(parents=True, exist_ok=True)

    for i, prot in enumerate(proteins):
        if prot not in proteomics.columns:
            out_rows.append(
                {"protein": prot, "n_sig": 0, "corr_json": json.dumps({}), "error": "missing_in_proteomics"}
            )
            continue

        seed = None if random_state is None else int(random_state) + int(i)

        s = compute_correlation_for_protein(
            transcriptomics=transcriptomics,
            proteomics=proteomics,
            protein_name=prot,
            p_threshold=p_threshold,
            bootstrap=bootstrap,
            n_bootstraps=n_bootstraps,
            random_state=seed,
        )

        nz = s[s != 0.0]
        corr_dict = {str(g): float(v) for g, v in nz.items()}
        out_rows.append({"protein": prot, "n_sig": int(len(corr_dict)), "corr_json": json.dumps(corr_dict)})

        if out_path is not None:
            safe = _sanitize_filename(prot)
            wide = s.to_frame().T
            wide.index = [prot]
            wide.index.name = "protein"

            if per_protein_format == "csv":
                wide.to_csv(out_path / per_protein_subdir / f"{safe}_corr.csv", index=True)
            else:
                _safe_to_parquet_or_csv(wide, out_path / per_protein_subdir / f"{safe}_corr.parquet", index=True)

    return pd.DataFrame(out_rows)
