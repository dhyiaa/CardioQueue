#!/usr/bin/env python3
"""Validate the normalized April 2026 new-data outputs."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def read_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader), list(reader.fieldnames or [])


def require(condition: bool, message: str, failures: list[str]) -> None:
    if not condition:
        failures.append(message)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="outputs/new_data")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    cohort, cohort_cols = read_rows(output_dir / "analysis_cohort.csv")
    variants, _ = read_rows(output_dir / "normalized_variants.csv")
    patients, patient_cols = read_rows(output_dir / "normalized_patients.csv")
    dictionary_rows, _ = read_rows(output_dir / "variable_dictionary_used.csv")
    report = json.loads((output_dir / "preprocessing_report.json").read_text(encoding="utf-8"))

    failures: list[str] = []

    require(len(variants) == 512, f"Expected 512 variant rows, found {len(variants)}", failures)
    require(len(cohort) == 512, f"Expected 512 cohort rows, found {len(cohort)}", failures)
    require(len(patients) == 234, f"Expected 234 patient rows, found {len(patients)}", failures)
    require(dictionary_rows, "Expected variable_dictionary_used.csv to contain rows", failures)

    sensitive_tokens = ("dob", "birth", "diagnosis_date")
    exported_cols = [*cohort_cols, *patient_cols]
    leaking_cols = [col for col in exported_cols if any(token in col.lower() for token in sensitive_tokens)]
    require(not leaking_cols, f"Date/birth columns exported unexpectedly: {leaking_cols}", failures)

    allowed_classes = {"Pathogenic", "VUS", "Benign", "Missing"}
    observed_classes = {row["classification_3class"] for row in cohort}
    require(
        observed_classes.issubset(allowed_classes),
        f"Unexpected 3-class labels: {sorted(observed_classes - allowed_classes)}",
        failures,
    )

    preferred_count = sum(row.get("preferred_for_training") == "true" for row in cohort)
    require(preferred_count == 511, f"Expected 511 preferred rows, found {preferred_count}", failures)

    match_count = sum(row.get("patient_match") == "true" for row in cohort)
    require(match_count == 483, f"Expected 483 patient-matched rows, found {match_count}", failures)

    metadata = report.get("metadata", {})
    dict_meta = metadata.get("hiro_variable_dictionary_file", {})
    require(
        dict_meta.get("definitions_loaded", 0) > 0,
        "HiRO dictionary metadata missing or empty",
        failures,
    )

    if failures:
        print("VALIDATION FAILED")
        for failure in failures:
            print(f"- {failure}")
        raise SystemExit(1)

    print("VALIDATION PASSED")
    print(f"- variants: {len(variants)}")
    print(f"- patients: {len(patients)}")
    print(f"- cohort rows: {len(cohort)}")
    print(f"- patient-matched variants: {match_count}")
    print(f"- preferred training rows: {preferred_count}")


if __name__ == "__main__":
    main()
