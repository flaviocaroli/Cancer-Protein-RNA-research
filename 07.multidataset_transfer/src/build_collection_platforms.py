#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from utils import build_collection_a_platforms, load_yaml


def main() -> None:
    p = argparse.ArgumentParser(
        description="Load Collection A cohort folders and merge them into TMT, iTRAQ, and LabelFree matrices."
    )
    p.add_argument("--project-root", default=str(Path(__file__).resolve().parents[1]))
    p.add_argument("--collection-dir", default="data/Collection_A")
    p.add_argument("--output-dir", default="data/processed/platforms")
    p.add_argument("--groups-config", default="configs/collection_a_groups.yaml")
    p.add_argument("--mapping-path", default="")
    p.add_argument("--allow-mygene-query", type=int, choices=[0, 1], default=0)
    p.add_argument("--feature-nonmissing-fraction", type=float, default=0.50)
    p.add_argument("--sample-nonmissing-fraction", type=float, default=0.50)
    p.add_argument("--map-protein-ids", type=int, choices=[0, 1], default=1)
    p.add_argument("--force", action="store_true")
    args = p.parse_args()

    root = Path(args.project_root).resolve()
    collection_dir = Path(args.collection_dir)
    if not collection_dir.is_absolute():
        collection_dir = root / collection_dir
    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = root / output_dir
    groups_path = Path(args.groups_config)
    if not groups_path.is_absolute():
        groups_path = root / groups_path
    mapping_path = Path(args.mapping_path) if args.mapping_path else None
    if mapping_path is not None and not mapping_path.is_absolute():
        mapping_path = root / mapping_path

    cohort_groups = load_yaml(groups_path)["collection_a"]
    cohort_summary, platform_summary = build_collection_a_platforms(
        collection_dir=collection_dir,
        output_dir=output_dir,
        cohort_groups=cohort_groups,
        mapping_path=mapping_path,
        allow_mygene_query=bool(args.allow_mygene_query),
        feature_nonmissing_fraction=args.feature_nonmissing_fraction,
        sample_nonmissing_fraction=args.sample_nonmissing_fraction,
        map_protein_ids=bool(args.map_protein_ids),
        force=args.force,
    )
    print("\nCollection A cohorts")
    print(cohort_summary.to_string(index=False))
    print("\nMerged platform matrices")
    print(platform_summary.to_string(index=False))


if __name__ == "__main__":
    main()
