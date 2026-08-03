#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PROJECT_DIR}"

COLLECTION_DIR="${COLLECTION_DIR:-${PROJECT_DIR}/data/Collection_A}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_DIR}/data/processed/platforms}"
GROUPS_CONFIG="${GROUPS_CONFIG:-${PROJECT_DIR}/configs/collection_a_groups.yaml}"
MAPPING_PATH="${MAPPING_PATH:-}"
ALLOW_MYGENE_QUERY="${ALLOW_MYGENE_QUERY:-0}"
FEATURE_NONMISSING_FRACTION="${FEATURE_NONMISSING_FRACTION:-0.50}"
SAMPLE_NONMISSING_FRACTION="${SAMPLE_NONMISSING_FRACTION:-0.50}"
MAP_PROTEIN_IDS="${MAP_PROTEIN_IDS:-1}"
FORCE="${FORCE:-0}"
PYTHON_EXECUTABLE="${PYTHON_EXECUTABLE:-$(command -v python)}"

ARGS=(
  --project-root "${PROJECT_DIR}"
  --collection-dir "${COLLECTION_DIR}"
  --output-dir "${OUTPUT_DIR}"
  --groups-config "${GROUPS_CONFIG}"
  --allow-mygene-query "${ALLOW_MYGENE_QUERY}"
  --feature-nonmissing-fraction "${FEATURE_NONMISSING_FRACTION}"
  --sample-nonmissing-fraction "${SAMPLE_NONMISSING_FRACTION}"
  --map-protein-ids "${MAP_PROTEIN_IDS}"
)
if [[ -n "${MAPPING_PATH}" ]]; then
  ARGS+=(--mapping-path "${MAPPING_PATH}")
fi
if [[ "${FORCE}" == "1" ]]; then
  ARGS+=(--force)
fi

"${PYTHON_EXECUTABLE}" src/build_collection_platforms.py "${ARGS[@]}"
