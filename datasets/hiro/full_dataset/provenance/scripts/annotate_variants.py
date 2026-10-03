#!/usr/bin/env python3
"""Annotate normalized variants with ClinVar, gnomAD, and MyVariant/dbNSFP."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cardio_pipeline.annotation import annotate_rows, read_csv_rows, write_annotation_outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="outputs/new_data/analysis_cohort.csv")
    parser.add_argument("--output-dir", default="outputs/annotations")
    parser.add_argument("--limit", type=int, default=None, help="Annotate only the first N rows.")
    parser.add_argument("--refresh", action="store_true", help="Ignore existing cache and fetch again.")
    parser.add_argument("--dry-run", action="store_true", help="Parse/query-plan only; do not call APIs.")
    parser.add_argument("--gnomad-dataset", default="gnomad_r4")
    parser.add_argument("--progress-every", type=int, default=25, help="Print progress every N rows.")
    parser.add_argument(
        "--preferred-labeled-only",
        action="store_true",
        help="Restrict to rows preferred for training with non-missing labels.",
    )
    args = parser.parse_args()

    rows, _ = read_csv_rows(args.input)
    if args.preferred_labeled_only:
        rows = [
            row
            for row in rows
            if row.get("preferred_for_training") == "true" and row.get("classification_3class") != "Missing"
        ]

    annotations, summary = annotate_rows(
        rows,
        output_dir=args.output_dir,
        limit=args.limit,
        refresh=args.refresh,
        dry_run=args.dry_run,
        dataset=args.gnomad_dataset,
        progress_every=args.progress_every,
    )
    paths = write_annotation_outputs(rows, annotations, summary, args.output_dir)

    print(json.dumps(summary, indent=2, sort_keys=True))
    print("\nWrote:")
    for name, path in paths.items():
        print(f"  {name}: {path}")


if __name__ == "__main__":
    main()
