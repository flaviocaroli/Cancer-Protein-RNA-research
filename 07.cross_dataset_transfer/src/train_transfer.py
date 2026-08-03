#!/usr/bin/env python
"""
train_transfer.py -- train per-protein models for ONE split + model family over
a PROTEIN CHUNK, saving one self-contained artifact per protein.

Reads a split directory produced by notebook 02 (already matched + harmonized).
Applies model-level normalization (zscore/binning/ranking) INSIDE the pipeline,
fit on X_train only, so each artifact carries its own normalizer + coefficients
and can predict on any dataset standalone.

Chunking for the 20-job HPC cap: the launcher passes --protein_start/--protein_end;
this script loops over exactly the proteins in that chunk.

Example:
  python src/train_transfer.py --split_dir data/splits/ccle_tmt_itraq__labelfree \
    --model enet --normalization zscore --select_k 5000 --cv_folds 5 \
    --protein_start 0 --protein_end 256 \
    --out_dir results/ccle_tmt_itraq__labelfree/enet_zscore_k5000
"""
import argparse
from pathlib import Path

import joblib

import io_utils as IO
from models import train_protein_model


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--split_dir", required=True)
    ap.add_argument("--model", choices=["enet", "xgb"], required=True)
    ap.add_argument("--normalization", default=None,
                    choices=["zscore", "binning", "ranking", "none", None],
                    help="Override; default = split_config['model_norm'].")
    ap.add_argument("--select_k", type=int, default=5000)
    ap.add_argument("--cv_folds", type=int, default=5)
    ap.add_argument("--n_bins", type=int, default=50)
    ap.add_argument("--max_combos", type=int, default=8)
    ap.add_argument("--protein_start", type=int, default=0)
    ap.add_argument("--protein_end", type=int, default=-1)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--skip_existing", type=int, default=1)
    args = ap.parse_args()

    cfg = IO.load_split_config(args.split_dir)
    norm = args.normalization or cfg.get("model_norm", "zscore")
    X_train, Y_train = IO.load_train(args.split_dir)
    proteins = list(Y_train.columns)
    end = len(proteins) if args.protein_end < 0 else min(args.protein_end, len(proteins))
    chunk = proteins[args.protein_start:end]

    out_dir = Path(args.out_dir); (out_dir / "models").mkdir(parents=True, exist_ok=True)
    print(f"Split={Path(args.split_dir).name} train={cfg['train']} -> test={cfg['test']}")
    print(f"X_train={X_train.shape} | model={args.model} norm={norm} k={args.select_k}")
    print(f"Protein chunk [{args.protein_start}:{end}] -> {len(chunk)} proteins")

    done = 0
    for protein in chunk:
        mpath = out_dir / "models" / f"{protein}.joblib"
        if args.skip_existing and mpath.exists():
            done += 1; continue
        if protein not in Y_train.columns:
            continue
        art = train_protein_model(
            X_train, Y_train[protein].to_numpy(), protein,
            model_family=args.model, normalization=norm,
            select_k=args.select_k, cv_folds=args.cv_folds, n_bins=args.n_bins,
            max_combos=args.max_combos, seed=args.seed, scenario=Path(args.split_dir).name,
            train_datasets=cfg["train"], test_datasets=cfg["test"],
        )
        if art is not None:
            joblib.dump(art, mpath, compress=3)
            done += 1
    print(f"Trained/kept {done}/{len(chunk)} models -> {out_dir/'models'}")


if __name__ == "__main__":
    main()
