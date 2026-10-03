#!/usr/bin/env python3
"""Validate frozen model-input artifacts for leakage and row-count consistency."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cardio_pipeline.features import (
    BASE_CATEGORICAL_FEATURES,
    BASE_NUMERIC_FEATURES,
    CLASS_TO_NUMERIC,
    ENRICHED_CATEGORICAL_FEATURES,
    ENRICHED_NUMERIC_FEATURES,
    LEAKAGE_COLUMNS,
)


def read_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader), list(reader.fieldnames or [])


def count_jsonl(path: Path) -> int:
    with path.open(encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default="outputs/model_inputs")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    manifest_path = input_dir / "model_input_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    failures: list[str] = []

    baseline_rows, baseline_cols = read_rows(input_dir / "ml_baseline_features.csv")
    enriched_rows, enriched_cols = read_rows(input_dir / "ml_enriched_features.csv")
    modeling_rows, modeling_cols = read_rows(input_dir / "modeling_cohort.csv")
    variant_rows, _ = read_rows(input_dir / "variant_only_cohort.csv")
    criterion_rows, criterion_cols = read_rows(input_dir / "criterion_reference.csv")

    expected_primary = manifest.get("primary_model_rows_with_phenotype")
    expected_preferred = manifest.get("preferred_labeled_rows")

    if len(baseline_rows) != expected_primary:
        failures.append(f"Baseline feature rows {len(baseline_rows)} != manifest primary rows {expected_primary}")
    if len(enriched_rows) != expected_primary:
        failures.append(f"Enriched feature rows {len(enriched_rows)} != manifest primary rows {expected_primary}")
    if len(modeling_rows) != expected_primary:
        failures.append(f"Modeling cohort rows {len(modeling_rows)} != manifest primary rows {expected_primary}")
    if len(variant_rows) != expected_preferred:
        failures.append(f"Variant-only rows {len(variant_rows)} != manifest preferred rows {expected_preferred}")
    if len(criterion_rows) != expected_preferred:
        failures.append(f"Criterion reference rows {len(criterion_rows)} != manifest preferred rows {expected_preferred}")

    if count_jsonl(input_dir / "llm_baseline_evidence.jsonl") != expected_primary:
        failures.append("LLM baseline JSONL count does not match primary modeling row count")
    if count_jsonl(input_dir / "llm_enriched_evidence.jsonl") != expected_primary:
        failures.append("LLM enriched JSONL count does not match primary modeling row count")

    baseline_required = set(BASE_NUMERIC_FEATURES + BASE_CATEGORICAL_FEATURES + ["target_3class", "target_numeric"])
    enriched_required = baseline_required | set(ENRICHED_NUMERIC_FEATURES + ENRICHED_CATEGORICAL_FEATURES)
    missing_baseline = sorted(baseline_required - set(baseline_cols))
    missing_enriched = sorted(enriched_required - set(enriched_cols))
    if missing_baseline:
        failures.append(f"Baseline feature file missing columns: {missing_baseline}")
    if missing_enriched:
        failures.append(f"Enriched feature file missing columns: {missing_enriched}")

    leakage_present = sorted(set(LEAKAGE_COLUMNS) & set(baseline_cols + enriched_cols))
    if leakage_present:
        failures.append(f"Leakage/sensitivity columns present in ML feature files: {leakage_present}")

    for name, rows in [("baseline", baseline_rows), ("enriched", enriched_rows)]:
        for i, row in enumerate(rows, start=2):
            if row.get("target_3class") not in CLASS_TO_NUMERIC:
                failures.append(f"{name} line {i}: invalid target_3class {row.get('target_3class')!r}")
            if str(row.get("target_numeric")) not in {"0", "1", "2"}:
                failures.append(f"{name} line {i}: invalid target_numeric {row.get('target_numeric')!r}")

    for required in ["criteria_codes", "ref_PM2", "ref_PP3", "ref_PVS1", "ref_BA1", "ref_BS1"]:
        if required not in criterion_cols:
            failures.append(f"Criterion reference missing {required}")

    if failures:
        print("VALIDATION FAILED")
        for failure in failures:
            print(f"- {failure}")
        raise SystemExit(1)

    print("VALIDATION PASSED")
    print(f"- primary modeling rows: {expected_primary}")
    print(f"- preferred/labeled rows: {expected_preferred}")
    print(f"- baseline ML columns: {len(baseline_cols)}")
    print(f"- enriched ML columns: {len(enriched_cols)}")


if __name__ == "__main__":
    main()
