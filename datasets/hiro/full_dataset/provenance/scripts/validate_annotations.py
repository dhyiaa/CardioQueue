#!/usr/bin/env python3
"""Validate external annotation output files."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cardio_pipeline.annotation import ANNOTATION_COLUMNS


VALID_CLINVAR = {"matched", "no_match", "skipped", "error", ""}
VALID_GNOMAD = {"matched", "absent", "skipped", "error", ""}
VALID_MYVARIANT = {"matched", "no_scores", "no_match", "skipped", "error", ""}
VALID_OVERALL = {"complete", "partial", "partial_with_error", "error", "no_external_match", "planned", ""}


def read_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader), list(reader.fieldnames or [])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation-dir", default="outputs/annotations")
    args = parser.parse_args()

    annotation_dir = Path(args.annotation_dir)
    rows, columns = read_rows(annotation_dir / "variant_annotations.csv")
    report = json.loads((annotation_dir / "annotation_report.json").read_text(encoding="utf-8"))

    failures: list[str] = []
    missing_cols = [col for col in ANNOTATION_COLUMNS if col not in columns]
    if missing_cols:
        failures.append(f"Missing annotation columns: {missing_cols}")

    for i, row in enumerate(rows, start=2):
        if row.get("clinvar_status") not in VALID_CLINVAR:
            failures.append(f"Line {i}: invalid ClinVar status {row.get('clinvar_status')!r}")
        if row.get("gnomad_status") not in VALID_GNOMAD:
            failures.append(f"Line {i}: invalid gnomAD status {row.get('gnomad_status')!r}")
        if row.get("myvariant_status") not in VALID_MYVARIANT:
            failures.append(f"Line {i}: invalid MyVariant status {row.get('myvariant_status')!r}")
        if row.get("annotation_status") not in VALID_OVERALL:
            failures.append(f"Line {i}: invalid overall status {row.get('annotation_status')!r}")
        if row.get("gnomad_status") == "absent" and row.get("gnomad_max_af") not in {"0", "0.0"}:
            failures.append(f"Line {i}: absent gnomAD variant should have AF 0.0")
        if row.get("clinvar_status") == "error" and not row.get("clinvar_error"):
            failures.append(f"Line {i}: ClinVar error status without error text")
        if row.get("gnomad_status") == "error" and not row.get("gnomad_error"):
            failures.append(f"Line {i}: gnomAD error status without error text")
        if row.get("myvariant_status") == "error" and not row.get("myvariant_error"):
            failures.append(f"Line {i}: MyVariant error status without error text")

    if report.get("rows_requested") != len(rows):
        failures.append(f"Report rows_requested={report.get('rows_requested')} but CSV has {len(rows)} rows")

    if failures:
        print("VALIDATION FAILED")
        for failure in failures:
            print(f"- {failure}")
        raise SystemExit(1)

    print("VALIDATION PASSED")
    print(f"- annotation rows: {len(rows)}")
    print(f"- status counts: {report.get('status_counts', {})}")
    print(f"- cache hits: {report.get('cache_hits')}")
    print(f"- new annotations: {report.get('new_annotations')}")


if __name__ == "__main__":
    main()
