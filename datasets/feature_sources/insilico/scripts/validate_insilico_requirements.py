#!/usr/bin/env python3
"""Validate in silico feature tables against required feature coverage."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


IMPLEMENTED_COLUMN_MAP = {
    "AlphaMissense": ["am_pathogenicity", "am_class", "alphamissense_status"],
    "CADD_PHRED": ["cadd_phred", "cadd_status"],
    "REVEL": ["revel", "revel_status"],
    "SIFT": ["sift", "sift_status"],
    "PolyPhen2_HDIV": ["polyphen2_hdiv", "polyphen2_hdiv_status"],
    "MetaLR": ["metalr", "metalr_status"],
    "FATHMM_XF": ["fathmm_xf", "fathmm_xf_status"],
}


def load_requirements(path: Path) -> dict:
    return json.loads(path.read_text())


def read_header(path: Path) -> list[str]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            return next(reader)
        except StopIteration as exc:
            raise ValueError(f"{path} is empty") from exc


def validate(args: argparse.Namespace) -> int:
    requirements = load_requirements(args.requirements)
    header = set(read_header(args.input))
    required = [item["feature"] for item in requirements["minimum_from_previous_si"]]
    missing = {}
    for feature in required:
        expected_columns = IMPLEMENTED_COLUMN_MAP.get(feature, [])
        absent = [column for column in expected_columns if column not in header]
        if absent:
            missing[feature] = absent
    if missing and args.strict:
        raise ValueError(f"Missing required in silico feature columns: {missing}")
    print(f"checked old-SI features={len(required)} missing_groups={len(missing)}")
    if missing:
        print(json.dumps(missing, indent=2, sort_keys=True))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument(
        "--requirements",
        default=Path("metadata/insilico_feature_requirements.json"),
        type=Path,
    )
    parser.add_argument("--strict", action="store_true")
    return parser


def main() -> int:
    return validate(build_parser().parse_args())


if __name__ == "__main__":
    raise SystemExit(main())

