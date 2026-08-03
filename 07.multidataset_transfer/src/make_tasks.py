#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--prepared-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--protein-chunks", type=int, default=20)
    p.add_argument("--protein-subset", default="")
    args = p.parse_args()

    prepared = Path(args.prepared_dir)
    proteins = pd.read_csv(prepared / "protein_universe.tsv", sep="\t")["protein"].astype(str).tolist()
    if args.protein_subset:
        sub = pd.read_csv(args.protein_subset, sep=None, engine="python")
        col = "protein" if "protein" in sub.columns else sub.columns[0]
        wanted = set(sub[col].astype(str))
        proteins = [p for p in proteins if p in wanted]
    chunks = np.array_split(np.arange(len(proteins)), min(args.protein_chunks, max(1, len(proteins))))
    rows = []
    for task_id, idx in enumerate(chunks, start=1):
        if len(idx) == 0:
            continue
        rows.append({
            "task_id": task_id,
            "protein_start": int(idx[0]),
            "protein_end": int(idx[-1]) + 1,
            "n_proteins": len(idx),
        })
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"[OK] wrote {len(rows)} tasks to {out}; proteins={len(proteins)}")
    if args.protein_subset:
        pd.Series(proteins, name="protein").to_csv(out.parent / "protein_subset_resolved.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
