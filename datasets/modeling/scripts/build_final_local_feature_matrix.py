#!/usr/bin/env python3
"""Freeze the current local feature blocks into one modeling matrix.

This intentionally excludes feature blocks still running/downloading, such as
gnomAD browser chunks, VEP/SpliceAI, and the full FoldX batch.
"""

from __future__ import annotations

import csv
import json
import argparse
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]

BASE = ROOT / "datasets/modeling/interim/final_modeling_table_with_clingen_protein.tsv"
DBNSFP = ROOT / "datasets/modeling/interim/dbnsfp_selected_features.tsv"
ALPHAMISSENSE = ROOT / "datasets/modeling/interim/alphamissense_direct_features.tsv"

OUT = ROOT / "datasets/modeling/interim/final_modeling_table_local_features_frozen.tsv"
SUMMARY = ROOT / "datasets/modeling/interim/final_modeling_table_local_features_frozen.summary.json"

KEY_FIELDS = {"variant_id", "genome_build", "chrom", "pos", "ref", "alt", "primary_gene", "genes"}
REPLACE_FROM_JOIN = {
    "alphamissense_direct_status",
    "alphamissense_direct_missing_reason",
}


def rel(path: Path) -> str:
    path = path.resolve()
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def read_feature_table(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = {}
        for row in reader:
            rows[row["variant_id"]] = row
    return fields, rows


def ordered_add_fields(base_fields: list[str], source_fields: list[str]) -> list[str]:
    add_fields = []
    for field in source_fields:
        if field in KEY_FIELDS:
            continue
        if field in base_fields and field not in REPLACE_FROM_JOIN:
            continue
        if field not in add_fields:
            add_fields.append(field)
    return add_fields


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=BASE)
    parser.add_argument("--dbnsfp", type=Path, default=DBNSFP)
    parser.add_argument("--alphamissense", type=Path, default=ALPHAMISSENSE)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--summary", type=Path, default=SUMMARY)
    args = parser.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)

    dbnsfp_fields, dbnsfp_rows = read_feature_table(args.dbnsfp)
    am_fields, am_rows = read_feature_table(args.alphamissense)

    counts: Counter[str] = Counter()
    with args.base.open(newline="") as in_handle, args.out.open("w", newline="") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        base_fields = list(reader.fieldnames or [])

        dbnsfp_add_fields = ordered_add_fields(base_fields, dbnsfp_fields)
        am_add_fields = ordered_add_fields(base_fields + dbnsfp_add_fields, am_fields)

        fieldnames = []
        for field in base_fields + dbnsfp_add_fields + am_add_fields:
            if field not in fieldnames:
                fieldnames.append(field)

        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()

        seen = set()
        for row in reader:
            counts["base_rows"] += 1
            variant_id = row["variant_id"]
            seen.add(variant_id)

            db_row = dbnsfp_rows.get(variant_id)
            if db_row:
                counts["dbnsfp_joined"] += 1
                for field in dbnsfp_add_fields:
                    row[field] = db_row.get(field, "")
            else:
                counts["dbnsfp_missing_payload"] += 1
                for field in dbnsfp_add_fields:
                    row[field] = ""

            am_row = am_rows.get(variant_id)
            if am_row:
                counts["alphamissense_direct_joined"] += 1
                for field in am_add_fields + sorted(REPLACE_FROM_JOIN):
                    if field in fieldnames:
                        row[field] = am_row.get(field, "")
                if am_row.get("alphamissense_direct_status") == "ok":
                    counts["alphamissense_direct_ok"] += 1
            else:
                counts["alphamissense_direct_missing_payload"] += 1
                for field in am_add_fields:
                    row[field] = ""

            writer.writerow(row)

    summary = {
        "inputs": {
            "base": rel(args.base),
            "dbnsfp_selected": rel(args.dbnsfp),
            "alphamissense_direct": rel(args.alphamissense),
        },
        "output": rel(args.out),
        "output_rows": counts["base_rows"],
        "output_columns": len(fieldnames),
        "added_columns": {
            "dbnsfp_selected": dbnsfp_add_fields,
            "alphamissense_direct": am_add_fields,
            "replaced_placeholders": sorted(REPLACE_FROM_JOIN),
        },
        "counts": dict(sorted(counts.items())),
        "notes": [
            "ClinGen and protein/structure features were already present in the base table.",
            "dbNSFP selected features were appended without overwriting existing coordinate/provenance columns.",
            "Direct AlphaMissense status/missing-reason placeholders were replaced with the completed local join.",
            "gnomAD browser chunks, VEP/SpliceAI, and full FoldX DDG are intentionally not included yet.",
        ],
    }
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
