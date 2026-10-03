#!/usr/bin/env python3
"""Audit dbNSFP coverage in the final modeling skeleton."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def variant_len_type(row: dict[str, str]) -> str:
    ref = (row.get("ref") or "").strip()
    alt = (row.get("alt") or "").strip()
    if not ref or not alt:
        return "missing"
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    return "indel_or_MNV"


def nested_count(rows: list[dict[str, str]], group_col: str) -> dict[str, dict[str, int]]:
    out: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        group = row.get(group_col) or "blank"
        out[group][row.get("dbnsfp_status") or "blank"] += 1
    return {group: dict(counter) for group, counter in sorted(out.items())}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="datasets/modeling/interim/final_modeling_table_skeleton.tsv")
    parser.add_argument("--output", default="datasets/feature_sources/insilico/interim/dbnsfp_modeling_skeleton_coverage.summary.json")
    args = parser.parse_args()

    rows = read_rows(Path(args.input))
    for row in rows:
        row["variant_len_type"] = variant_len_type(row)
        row["is_protein_mapped"] = "true" if row.get("protein_feature_status") == "ok" else "false"

    summary = {
        "input": args.input,
        "rows": len(rows),
        "dbnsfp_status_counts": dict(Counter(row.get("dbnsfp_status") or "blank" for row in rows)),
        "by_variant_len_type": nested_count(rows, "variant_len_type"),
        "by_model_label_3class": nested_count(rows, "model_label_3class"),
        "by_sources": nested_count(rows, "sources"),
        "by_protein_variant_type": nested_count(rows, "protein_variant_type"),
        "by_protein_mapped": nested_count(rows, "is_protein_mapped"),
        "interpretation": [
            "The local dbNSFP GRCh38 file is valid and exact chrom-pos-ref-alt joins work.",
            "All dbNSFP hits in this skeleton are SNVs; indel/MNV variants are expected to be not_found in this dbNSFP SNV-oriented join.",
            "Many remaining not_found SNVs are ClinVar variants outside dbNSFP's functional/transcript-scored scope, especially noncoding/UTR/synonymous or otherwise not protein-mapped rows.",
            "For non-missense/indel/splice variants, use VEP/SpliceAI, ClinVar consequence fields, and protein/LoF-specific annotations rather than expecting dbNSFP missense scores.",
        ],
    }
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
