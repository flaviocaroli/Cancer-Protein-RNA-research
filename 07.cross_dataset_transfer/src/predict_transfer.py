#!/usr/bin/env python
"""
predict_transfer.py -- apply trained per-protein artifacts to the held-out test
dataset(s) of a split and write per-protein metrics + long predictions.

Test matrices are already harmonized (notebook 02) with train-anchored params;
each artifact re-applies its own train-fitted feature normalizer, so no test
statistics leak into the transform.
"""
import argparse
from pathlib import Path

import joblib
import pandas as pd

import io_utils as IO
from models import predict_with_artifact
from metrics import protein_metrics


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--split_dir", required=True)
    ap.add_argument("--models_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--save_predictions", type=int, default=1)
    args = ap.parse_args()

    cfg = IO.load_split_config(args.split_dir)
    tests = IO.load_tests(args.split_dir)
    models_dir = Path(args.models_dir)
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    model_files = sorted(models_dir.glob("*.joblib"))
    print(f"{len(model_files)} models; test datasets={list(tests.keys())}")

    metric_rows, pred_rows = [], []
    for mf in model_files:
        art = joblib.load(mf); protein = art["protein"]
        for dset, (X_test, Y_test) in tests.items():
            if protein not in Y_test.columns:
                continue
            y_true = Y_test[protein].to_numpy(float)
            y_pred = predict_with_artifact(art, X_test)
            m = protein_metrics(y_true, y_pred)
            metric_rows.append({
                "split": Path(args.split_dir).name, "test_dataset": dset, "protein": protein,
                "model_family": art["model_family"], "normalization": art["normalization"],
                "select_k": art["select_k"], "n_train": art["n_train"],
                "cv_pearson_mean": art.get("cv_pearson_mean"),
                "n_nonzero_coef": art.get("n_nonzero_coef"), **m})
            if args.save_predictions:
                for s, yt, yp in zip(Y_test.index, y_true, y_pred):
                    pred_rows.append((Path(args.split_dir).name, dset, protein, s, float(yt), float(yp)))

    met = pd.DataFrame(metric_rows)
    met.to_csv(out_dir / "metrics_by_protein.csv", index=False)
    print(f"Wrote {out_dir/'metrics_by_protein.csv'} ({len(met)} rows)")
    if len(met):
        print(met.groupby("test_dataset")["test_pearson_corr"].agg(["median","mean","count"]).to_string())
    if args.save_predictions and pred_rows:
        preds = pd.DataFrame(pred_rows, columns=["split","test_dataset","protein","sample","y_true","y_pred"])
        try: preds.to_parquet(out_dir / "predictions_long.parquet", index=False)
        except Exception: preds.to_csv(out_dir / "predictions_long.csv.gz", index=False)


if __name__ == "__main__":
    main()
