#!/usr/bin/env python3
"""Create a model-ready direct AlphaMissense feature table.

This reads the exact GRCh38 AlphaMissense join output and writes a compact
modeling table. It keeps direct AlphaMissense separate from dbNSFP-derived
AlphaMissense fields so the two sources can be compared or reconciled later.
"""

from __future__ import annotations

import csv
import json
import argparse
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]

INPUT = ROOT / "datasets/feature_sources/insilico/interim/alphamissense_final_modeling_skeleton.tsv"
OUTPUT = ROOT / "datasets/modeling/interim/alphamissense_direct_features.tsv"
SUMMARY = ROOT / "datasets/modeling/interim/alphamissense_direct_features.summary.json"


def rel(path: Path) -> str:
    path = path.resolve()
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def clean(value: object) -> str:
    text = str(value if value is not None else "").strip()
    return "" if text in {"", ".", "NA", "N/A", "nan", "None"} else text


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--summary", type=Path, default=SUMMARY)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    counts: Counter = Counter()
    fieldnames = [
        "variant_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "alphamissense_direct_status",
        "alphamissense_direct_missing_reason",
        "alphamissense_direct_score",
        "alphamissense_direct_class",
        "alphamissense_direct_uniprot_id",
        "alphamissense_direct_transcript_id",
        "alphamissense_direct_protein_variant",
        "alphamissense_direct_score_is_missing",
        "alphamissense_direct_class_is_missing",
    ]

    with args.input.open(newline="") as in_handle, args.output.open("w", newline="") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in reader:
            counts["rows"] += 1
            status = clean(row.get("alphamissense_status"))
            score = clean(row.get("am_pathogenicity"))
            am_class = clean(row.get("am_class"))
            out = {
                "variant_id": clean(row.get("variant_id")),
                "chrom": clean(row.get("chrom")),
                "pos": clean(row.get("pos")),
                "ref": clean(row.get("ref")),
                "alt": clean(row.get("alt")),
                "alphamissense_direct_status": status,
                "alphamissense_direct_missing_reason": clean(row.get("alphamissense_missing_reason")),
                "alphamissense_direct_score": score,
                "alphamissense_direct_class": am_class,
                "alphamissense_direct_uniprot_id": clean(row.get("am_uniprot_id")),
                "alphamissense_direct_transcript_id": clean(row.get("am_transcript_id")),
                "alphamissense_direct_protein_variant": clean(row.get("am_protein_variant")),
                "alphamissense_direct_score_is_missing": "true" if not score else "false",
                "alphamissense_direct_class_is_missing": "true" if not am_class else "false",
            }
            counts[f"alphamissense_direct_status_{status or 'blank'}"] += 1
            counts[f"alphamissense_direct_class_{am_class or 'blank'}"] += 1
            writer.writerow(out)

    summary = {
        "inputs": {
            "exact_join": rel(args.input),
            "raw_alphamissense": "datasets/feature_sources/insilico/raw/alphamissense/AlphaMissense_hg38.tsv.gz",
        },
        "outputs": {
            "features": rel(args.output),
            "summary": rel(args.summary),
        },
        "counts": dict(sorted(counts.items())),
        "notes": [
            "Direct exact GRCh38 chrom-pos-ref-alt AlphaMissense join.",
            "This is independent of dbNSFP AlphaMissense columns.",
            "Missing status includes non-missense variants and variants absent from AlphaMissense.",
        ],
    }
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
