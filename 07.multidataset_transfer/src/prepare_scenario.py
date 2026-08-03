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

from utils import (
    apply_eval_y_zscore,
    concatenate_domains,
    create_preparation_figures,
    dump_json,
    feature_intersection,
    filter_protein_universe,
    fit_eval_y_zscore,
    fit_x_normalizer,
    fit_y_transformer,
    load_and_align_dataset,
    load_yaml,
    transform_x,
    transform_y,
    write_table,
)


def resolve(root: Path, value: str | None) -> Path | None:
    if value is None:
        return None
    p = Path(value)
    return p if p.is_absolute() else root / p


def prepare_scenario(
    project_root: Path,
    scenario_name: str,
    x_norm: str,
    y_norm: str,
    n_bins: int,
    eps: float,
    domain_balanced_normalizer: bool,
    quantile_fit_max_samples_per_domain: int,
    make_figures: bool,
    force: bool,
    train_datasets_override: list[str] | None = None,
    test_datasets_override: list[str] | None = None,
    output_name: str | None = None,
) -> Path:
    datasets_cfg = load_yaml(project_root / "configs/datasets.yaml")
    scenarios_cfg = load_yaml(project_root / "configs/scenarios.yaml")
    if scenario_name not in scenarios_cfg["scenarios"] and not (train_datasets_override and test_datasets_override):
        raise KeyError(f"Unknown scenario {scenario_name!r}; provide explicit train/test overrides or add it to configs/scenarios.yaml")
    scenario = scenarios_cfg["scenarios"].get(scenario_name, {"description": "Custom notebook-defined scenario"})
    train_datasets = list(train_datasets_override or scenario.get("train", []))
    test_datasets = list(test_datasets_override or scenario.get("test", []))
    if not train_datasets or not test_datasets:
        raise ValueError("Both train_datasets and test_datasets must be non-empty")
    overlap = set(train_datasets).intersection(test_datasets)
    if overlap:
        raise ValueError(f"Datasets cannot be both train and test: {sorted(overlap)}")
    all_datasets = train_datasets + test_datasets
    matching = datasets_cfg.get("matching", {})

    scenario_output_name = output_name or scenario_name
    out_dir = project_root / "prepared" / scenario_output_name / f"x-{x_norm}__y-{y_norm}"
    manifest_path = out_dir / "manifest.json"
    if manifest_path.exists() and not force:
        print(f"[SKIP] Prepared bundle already exists: {out_dir}")
        return out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "tests").mkdir(exist_ok=True)
    (out_dir / "normalization").mkdir(exist_ok=True)

    rna_frames: dict[str, pd.DataFrame] = {}
    protein_frames: dict[str, pd.DataFrame] = {}
    metadata_frames: dict[str, pd.DataFrame] = {}
    dataset_info = []

    for dataset_name in all_datasets:
        cfg = datasets_cfg["datasets"][dataset_name]
        rna_path = resolve(project_root, cfg["rna"])
        protein_path = resolve(project_root, cfg["protein"])
        metadata_path = resolve(project_root, cfg.get("metadata"))
        rna, protein, meta, info = load_and_align_dataset(
            dataset_name=dataset_name,
            rna_path=rna_path,
            protein_path=protein_path,
            metadata_path=metadata_path,
            feature_case=matching.get("feature_case", "upper"),
            prefix_sample_ids=bool(matching.get("prefix_sample_ids", True)),
        )
        rna_frames[dataset_name] = rna
        protein_frames[dataset_name] = protein
        metadata_frames[dataset_name] = meta
        dataset_info.append(info)
        print(f"[LOAD] {dataset_name}: RNA={rna.shape}, protein={protein.shape}")

    if matching.get("gene_universe", "all_domains_intersection") == "all_domains_intersection":
        gene_universe = feature_intersection(rna_frames, all_datasets)
    else:
        gene_universe = feature_intersection(rna_frames, train_datasets)

    protein_universe = filter_protein_universe(
        protein_frames=protein_frames,
        train_datasets=train_datasets,
        test_datasets=test_datasets,
        min_train_fraction=float(matching.get("min_train_nonmissing_fraction", 0.50)),
        min_test_fraction=float(matching.get("min_test_nonmissing_fraction", 0.20)),
        mode=matching.get("protein_universe", "all_domains_intersection"),
    )

    if len(gene_universe) == 0:
        raise RuntimeError("The matched gene universe is empty.")
    if len(protein_universe) == 0:
        raise RuntimeError("The matched protein universe is empty.")

    print(f"[MATCH] genes={len(gene_universe):,}; proteins={len(protein_universe):,}")

    train_X_raw = concatenate_domains(rna_frames, train_datasets, gene_universe)
    train_Y_raw = concatenate_domains(protein_frames, train_datasets, protein_universe)
    train_meta = pd.concat([metadata_frames[d] for d in train_datasets], axis=0).reindex(train_X_raw.index)

    x_bundle = fit_x_normalizer(
        train_X_raw,
        train_meta["dataset"],
        method=x_norm,
        n_bins=n_bins,
        eps=eps,
        domain_balanced=domain_balanced_normalizer,
        quantile_fit_max_samples_per_domain=quantile_fit_max_samples_per_domain,
        seed=42,
    )
    train_X_norm = transform_x(train_X_raw, x_bundle)

    y_bundle = fit_y_transformer(train_Y_raw, train_meta["dataset"], method=y_norm, eps=eps)
    train_Y_norm = transform_y(train_Y_raw, train_meta["dataset"], y_bundle)

    joblib.dump(train_X_raw, out_dir / "train_X_raw.joblib", compress=3)
    joblib.dump(train_X_norm, out_dir / "train_X_normalized.joblib", compress=3)
    joblib.dump(train_Y_raw, out_dir / "train_Y_raw.joblib", compress=3)
    joblib.dump(train_Y_norm, out_dir / "train_Y_normalized.joblib", compress=3)
    train_meta.to_csv(out_dir / "train_metadata.tsv", sep="\t")
    pd.Series(gene_universe, name="gene").to_csv(out_dir / "gene_universe.tsv", sep="\t", index=False)
    pd.Series(protein_universe, name="protein").to_csv(out_dir / "protein_universe.tsv", sep="\t", index=False)
    joblib.dump(x_bundle, out_dir / "normalization/x_normalizer.joblib", compress=3)
    joblib.dump(y_bundle, out_dir / "normalization/y_train_transformer.joblib", compress=3)

    # Human-readable/auditable copies of the learned parameters.
    x_parameter_files = []
    if x_bundle["method"] == "zscore":
        x_params = pd.DataFrame({
            "gene": x_bundle["features"],
            "training_mean": x_bundle["mean"],
            "training_sd": x_bundle["sd"],
        })
        x_parameter_files.append(str(write_table(x_params, out_dir / "normalization/x_zscore_parameters.parquet")))
    elif x_bundle["method"] in {"binning", "quantile_bin"}:
        edge_table = pd.DataFrame(
            x_bundle["edges"].T,
            columns=[f"q_{float(q):.6f}" for q in x_bundle["quantiles"]],
        )
        edge_table.insert(0, "gene", x_bundle["features"])
        edge_table["training_median"] = x_bundle["median"]
        x_parameter_files.append(str(write_table(edge_table, out_dir / "normalization/x_binning_edges.parquet")))
    elif x_bundle["method"] == "none":
        x_params = pd.DataFrame({"gene": x_bundle["features"], "training_median": x_bundle["median"]})
        x_parameter_files.append(str(write_table(x_params, out_dir / "normalization/x_imputation_medians.parquet")))
    else:
        rank_note = {
            "method": "sample_wise_percentile_rank",
            "fitted_parameters": False,
            "missing_value_rule": "missing genes receive percentile 0.5",
            "features": len(x_bundle["features"]),
        }
        dump_json(rank_note, out_dir / "normalization/x_rank_definition.json")
        x_parameter_files.append(str(out_dir / "normalization/x_rank_definition.json"))

    y_parameter_rows = []
    if y_bundle["method"] == "domain_zscore":
        for domain, params in y_bundle["params"].items():
            for protein, mean, sd in zip(y_bundle["proteins"], params["mean"], params["sd"]):
                y_parameter_rows.append({"training_domain": domain, "protein": protein, "training_mean": mean, "training_sd": sd})
    elif y_bundle["method"] == "pooled_train_zscore":
        params = y_bundle["params"]["__pooled__"]
        for protein, mean, sd in zip(y_bundle["proteins"], params["mean"], params["sd"]):
            y_parameter_rows.append({"training_domain": "__pooled__", "protein": protein, "training_mean": mean, "training_sd": sd})
    y_parameter_file = None
    if y_parameter_rows:
        y_parameter_file = str(write_table(pd.DataFrame(y_parameter_rows), out_dir / "normalization/y_training_parameters.parquet"))

    test_manifests = {}
    for test_name in test_datasets:
        test_dir = out_dir / "tests" / test_name
        test_dir.mkdir(parents=True, exist_ok=True)
        X_raw = rna_frames[test_name].reindex(columns=gene_universe)
        Y_raw = protein_frames[test_name].reindex(columns=protein_universe)
        meta = metadata_frames[test_name].reindex(X_raw.index)
        X_norm = transform_x(X_raw, x_bundle)

        # Correlation metrics use raw Y. Error/calibration metrics use a clearly labeled,
        # evaluation-only within-test z-score. The latter is never used to fit or predict.
        if y_norm == "pooled_train_zscore":
            Y_eval = transform_y(Y_raw, meta["dataset"], y_bundle, strict_unknown=False)
            y_eval_bundle = {"method": "pooled_train_zscore_applied_to_test"}
        else:
            y_eval_bundle = fit_eval_y_zscore(Y_raw, eps=eps)
            Y_eval = apply_eval_y_zscore(Y_raw, y_eval_bundle)

        joblib.dump(X_raw, test_dir / "X_raw.joblib", compress=3)
        joblib.dump(X_norm, test_dir / "X_normalized.joblib", compress=3)
        joblib.dump(Y_raw, test_dir / "Y_raw.joblib", compress=3)
        joblib.dump(Y_eval, test_dir / "Y_evaluation_scaled.joblib", compress=3)
        joblib.dump(y_eval_bundle, test_dir / "y_evaluation_scaler.joblib", compress=3)
        meta.to_csv(test_dir / "metadata.tsv", sep="\t")
        test_manifests[test_name] = {
            "n_samples": len(X_raw),
            "n_genes": X_raw.shape[1],
            "n_proteins": Y_raw.shape[1],
            "rna_missing_fraction": float(X_raw.isna().mean().mean()),
            "protein_missing_fraction": float(Y_raw.isna().mean().mean()),
        }

    summary_rows = []
    for info in dataset_info:
        summary_rows.append({
            "dataset": info["dataset"],
            "n_samples": info["n_common_samples"],
            "rna_features_before": info["rna_aligned_shape"][1],
            "protein_features_before": info["protein_aligned_shape"][1],
            "rna_missing_fraction": info["rna_missing_fraction"],
            "protein_missing_fraction": info["protein_missing_fraction"],
            "role": "train" if info["dataset"] in train_datasets else "test",
        })
    dataset_summary = pd.DataFrame(summary_rows)
    dataset_summary.to_csv(out_dir / "dataset_summary.tsv", sep="\t", index=False)

    manifest = {
        "scenario": scenario_output_name,
        "description": scenario.get("description", ""),
        "train_datasets": train_datasets,
        "test_datasets": test_datasets,
        "x_normalization": x_norm,
        "y_normalization": y_norm,
        "n_bins": n_bins,
        "eps": eps,
        "domain_balanced_normalizer": domain_balanced_normalizer,
        "x_parameter_files": x_parameter_files,
        "y_parameter_file": y_parameter_file,
        "gene_universe_size": len(gene_universe),
        "protein_universe_size": len(protein_universe),
        "train_samples": len(train_X_raw),
        "test_manifests": test_manifests,
        "dataset_info": dataset_info,
        "important_target_note": (
            "domain_zscore trains a dimensionless within-training-domain relative-abundance target. "
            "Held-out-test Pearson/Spearman use raw test protein values. Test-scaled MSE/R2/calibration "
            "use an evaluation-only test z-score and must not be described as deployment calibration."
        ),
    }
    dump_json(manifest, manifest_path)

    if make_figures:
        figure_dir = project_root / "figures" / "preparation" / scenario_output_name / f"x-{x_norm}__y-{y_norm}"
        create_preparation_figures(
            raw_train_x=train_X_raw,
            normalized_train_x=train_X_norm,
            raw_train_y=train_Y_raw,
            normalized_train_y=train_Y_norm,
            train_meta=train_meta,
            dataset_summary=dataset_summary,
            output_dir=figure_dir,
            scenario=scenario_output_name,
            x_norm=x_norm,
        )
        print(f"[FIGURES] {figure_dir}")

    print(f"[DONE] Prepared bundle: {out_dir}")
    return out_dir


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Match and normalize a multi-dataset transfer scenario.")
    p.add_argument("--project-root", default=str(Path(__file__).resolve().parents[1]))
    p.add_argument("--scenario", required=True)
    p.add_argument("--train-datasets", default="", help="Optional comma-separated override, e.g. CCLE,TMT")
    p.add_argument("--test-datasets", default="", help="Optional comma-separated override, e.g. iTRAQ,LabelFree")
    p.add_argument("--output-name", default="", help="Optional prepared-folder name for a custom split")
    p.add_argument("--x-norm", choices=["zscore", "rank", "binning", "none"], default="zscore")
    p.add_argument("--y-norm", choices=["domain_zscore", "pooled_train_zscore", "none"], default="domain_zscore")
    p.add_argument("--n-bins", type=int, default=50)
    p.add_argument("--eps", type=float, default=1e-8)
    p.add_argument("--domain-balanced-normalizer", type=int, choices=[0, 1], default=1)
    p.add_argument("--quantile-fit-max-samples-per-domain", type=int, default=1000)
    p.add_argument("--make-figures", type=int, choices=[0, 1], default=1)
    p.add_argument("--force", action="store_true")
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    prepare_scenario(
        project_root=Path(args.project_root).resolve(),
        scenario_name=args.scenario,
        x_norm=args.x_norm,
        y_norm=args.y_norm,
        n_bins=args.n_bins,
        eps=args.eps,
        domain_balanced_normalizer=bool(args.domain_balanced_normalizer),
        quantile_fit_max_samples_per_domain=args.quantile_fit_max_samples_per_domain,
        make_figures=bool(args.make_figures),
        force=args.force,
        train_datasets_override=[x.strip() for x in args.train_datasets.split(",") if x.strip()] or None,
        test_datasets_override=[x.strip() for x in args.test_datasets.split(",") if x.strip()] or None,
        output_name=args.output_name or None,
    )
