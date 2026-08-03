#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import joblib
import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from utils import load_dataframe, transform_x, write_table


def main() -> None:
    p = argparse.ArgumentParser(description="Apply one saved per-protein model to a new samples-by-genes RNA matrix.")
    p.add_argument("--model", required=True, help="Per-protein .joblib model artifact")
    p.add_argument("--normalizer", default="", help="Optional x_normalizer.joblib; defaults to the path stored in the artifact")
    p.add_argument("--rna", required=True, help="Samples x genes RNA DataFrame/joblib")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    artifact = joblib.load(args.model)
    normalizer_path = args.normalizer or artifact.get("shared_normalization_bundle")
    if not normalizer_path:
        raise ValueError("No normalizer supplied and the model artifact contains no shared_normalization_bundle path")
    normalizer = joblib.load(normalizer_path)
    X = load_dataframe(args.rna)
    X.index = X.index.astype(str)
    X.columns = X.columns.astype(str).str.upper()
    X_norm = transform_x(X, normalizer)
    features = artifact["selected_features"]
    pred = artifact["model"].predict(X_norm[features].to_numpy())
    out = pd.DataFrame({
        "sample_id": X_norm.index,
        "protein": artifact["protein"],
        "prediction_relative_protein_abundance": pred,
        "model_type": artifact["model_type"],
        "scenario": artifact["scenario"],
    })
    write_table(out, args.out)
    print(f"[DONE] {args.out}")


if __name__ == "__main__":
    main()
