#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import numpy as np
import pandas as pd


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--project-root", default=str(Path(__file__).resolve().parents[1]))
    p.add_argument("--genes", type=int, default=80)
    p.add_argument("--proteins", type=int, default=12)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()
    root = Path(args.project_root).resolve()
    rng = np.random.default_rng(args.seed)

    symbols = [f"G{i:04d}" for i in range(args.genes)]
    ensembl = [f"ENSG{i:011d}" for i in range(1, args.genes + 1)]
    proteins = symbols[:args.proteins]
    protein_ensembl = ensembl[:args.proteins]
    mapping = pd.DataFrame({"ensembl_gene_id": ensembl, "gene_symbol": symbols})
    mapping_dir = root / "data" / "reference"
    mapping_dir.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(mapping_dir / "ensembl_to_symbol.tsv", sep="\t", index=False)

    W = rng.normal(0, 0.30, size=(args.genes, args.proteins))
    W[rng.random(W.shape) < 0.90] = 0

    # CCLE arrives already as samples x gene symbols.
    n_ccle = 40
    X = rng.normal(size=(n_ccle, args.genes))
    Y = X @ W + rng.normal(scale=0.5, size=(n_ccle, args.proteins))
    ccle_ids = [f"CCLE_{i:03d}" for i in range(n_ccle)]
    ccle_dir = root / "data" / "processed" / "CCLE"
    ccle_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(pd.DataFrame(X, index=ccle_ids, columns=symbols), ccle_dir / "transcriptomics.joblib")
    joblib.dump(pd.DataFrame(Y, index=ccle_ids, columns=proteins), ccle_dir / "proteomics.joblib")
    pd.DataFrame({"sample_id": ccle_ids, "cohort": "CCLE", "dataset": "CCLE"}).to_csv(
        ccle_dir / "metadata.tsv", sep="\t", index=False
    )

    # Minimal Collection A raw folder: one cohort per platform, features x samples.
    cohorts = [
        ("colon_2019", "TMT", 32, 0.4),
        ("breast_2016", "iTRAQ", 28, 0.8),
        ("brain", "LabelFree", 26, 1.2),
    ]
    collection_dir = root / "data" / "Collection_A"
    for cohort, platform, n, shift in cohorts:
        X = rng.normal(loc=shift, scale=1.1, size=(n, args.genes))
        Y = X @ W + rng.normal(scale=0.6, size=(n, args.proteins)) + shift
        X[rng.random(X.shape) < 0.01] = np.nan
        Y[rng.random(Y.shape) < 0.02] = np.nan
        original_ids = [f"{cohort}_sample_{i:03d}" for i in range(n)]
        cohort_dir = collection_dir / cohort
        cohort_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(X.T, index=ensembl, columns=original_ids).to_csv(
            cohort_dir / f"RNA_matrix_{cohort}.txt", sep="\t"
        )
        pd.DataFrame(Y.T, index=protein_ensembl, columns=original_ids).to_csv(
            cohort_dir / f"protein_matrix_{cohort}.txt", sep="\t"
        )

    print("[DONE] Synthetic CCLE and Collection A data written")


if __name__ == "__main__":
    main()
