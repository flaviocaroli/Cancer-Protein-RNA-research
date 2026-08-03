"""
io_utils.py -- load the per-split matrices written by notebook 02.

Notebook 02 (match + harmonize) writes, for a chosen split, one directory:

    <split_dir>/
        X_train.joblib          (train-union samples x genes, harmonized)
        Y_train.joblib          (train-union samples x proteins, harmonized)
        X_test__<D>.joblib      one per held-out test dataset D
        Y_test__<D>.joblib
        common_proteins.joblib  list[str]
        common_genes.joblib     list[str]
        split_config.json       {train, test, harmonization, model_norm, ...}
        harmonization_params.joblib

The HPC pipeline reads a split_dir directly; it does NOT re-assemble unions or
re-harmonize. Model-level normalization (zscore/binning/ranking) is applied
inside the pipeline (fit on X_train, applied to X_test) so each per-protein
artifact stays self-contained.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Tuple

import joblib
import pandas as pd


def load_joblib(path):
    return joblib.load(Path(path))


def load_split_config(split_dir) -> Dict:
    return json.loads((Path(split_dir) / "split_config.json").read_text())


def load_train(split_dir) -> Tuple[pd.DataFrame, pd.DataFrame]:
    d = Path(split_dir)
    X = load_joblib(d / "X_train.joblib"); Y = load_joblib(d / "Y_train.joblib")
    X.columns = X.columns.map(str); Y.columns = Y.columns.map(str)
    return X, Y


def load_tests(split_dir) -> Dict[str, Tuple[pd.DataFrame, pd.DataFrame]]:
    d = Path(split_dir)
    cfg = load_split_config(split_dir)
    tests = {}
    for ds in cfg["test"]:
        X = load_joblib(d / f"X_test__{ds}.joblib")
        Y = load_joblib(d / f"Y_test__{ds}.joblib")
        X.columns = X.columns.map(str); Y.columns = Y.columns.map(str)
        tests[ds] = (X, Y)
    return tests


def protein_list(split_dir) -> List[str]:
    return [str(p) for p in load_joblib(Path(split_dir) / "common_proteins.joblib")]
