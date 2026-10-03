#!/usr/bin/env python3
"""Merge in silico feature tables by variant_id."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_table(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "variant_id" not in reader.fieldnames:
            raise ValueError(f"{path} must be a TSV with a variant_id column")
        rows = {}
        for row in reader:
            variant_id = row["variant_id"]
            if variant_id in rows:
                raise ValueError(f"Duplicate variant_id in {path}: {variant_id}")
            rows[variant_id] = row
    return reader.fieldnames, rows


def merge(args: argparse.Namespace) -> int:
    alpha_fields, alpha_rows = read_table(args.alphamissense)
    cadd_fields, cadd_rows = read_table(args.cadd)
    variant_ids = sorted(set(alpha_rows) | set(cadd_rows))
    fields = ["variant_id"]
    for source_fields in [alpha_fields, cadd_fields]:
        for field in source_fields:
            if field != "variant_id" and field not in fields:
                fields.append(field)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for variant_id in variant_ids:
            row = {"variant_id": variant_id}
            row.update(alpha_rows.get(variant_id, {}))
            row.update(cadd_rows.get(variant_id, {}))
            writer.writerow(row)
    print(
        f"merged in silico features: variants={len(variant_ids)} "
        f"alphamissense={len(alpha_rows)} cadd={len(cadd_rows)}",
        file=None,
    )
    return 0


def validate(args: argparse.Namespace) -> int:
    with args.input.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    required = ["variant_id", "alphamissense_status", "cadd_status"]
    missing = [field for field in required if field not in (reader.fieldnames or [])]
    if missing:
        raise ValueError(f"{args.input} missing required columns: {missing}")
    if not rows:
        raise ValueError(f"{args.input} has no data rows")
    print(f"validated merged in silico rows={len(rows)}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(required=True)

    merge_parser = subparsers.add_parser("merge")
    merge_parser.add_argument("--alphamissense", required=True, type=Path)
    merge_parser.add_argument("--cadd", required=True, type=Path)
    merge_parser.add_argument("--output", required=True, type=Path)
    merge_parser.set_defaults(func=merge)

    validate_parser = subparsers.add_parser("validate")
    validate_parser.add_argument("--input", required=True, type=Path)
    validate_parser.set_defaults(func=validate)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())

