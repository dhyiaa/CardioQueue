#!/usr/bin/env python3
"""Extract immediate gnomAD-like population features from local dbNSFP joins."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
INPUT = ROOT / "datasets/variant_registry/interim/clean_combined_local_feature_registry.tsv"
OUT_DIR = ROOT / "datasets/feature_sources/gnomad_v4/interim"
OUTPUT = OUT_DIR / "gnomad_from_dbnsfp.clean_combined.tsv"
SUMMARY = OUT_DIR / "gnomad_from_dbnsfp.clean_combined.summary.json"

FIELDS = [
    "variant_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "genes",
    "sources",
    "labels_3class",
    "gnomad_proxy_status",
    "gnomad_proxy_missing_reason",
    "gnomad_joint_af",
    "gnomad_joint_nhomalt",
    "gnomad_joint_popmax_af",
    "gnomad_joint_popmax_nhomalt",
    "gnomad_dbnsfp_popmax_af",
    "gnomad_dbnsfp_popmax_pop",
    "gnomad_proxy_source",
]


def has_value(value: str | None) -> bool:
    return bool((value or "").strip())


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    counts: Counter = Counter()
    with INPUT.open(newline="") as source, OUTPUT.open("w", newline="") as out:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(out, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        for row in reader:
            counts["rows"] += 1
            joint_af = row.get("dbnsfp_gnomAD4.1_joint_AF", "")
            popmax_af = row.get("dbnsfp_gnomAD4.1_joint_POPMAX_AF", "")
            dbnsfp_popmax_af = row.get("dbnsfp_dbNSFP_POPMAX_AF", "")
            has_any = any(has_value(x) for x in [joint_af, popmax_af, dbnsfp_popmax_af])
            status = "ok" if has_any else "not_found"
            missing_reason = "" if has_any else "no gnomAD population-frequency fields in local dbNSFP row"
            counts[f"gnomad_proxy_{status}"] += 1
            if has_value(joint_af):
                counts["nonmissing_gnomad_joint_af"] += 1
            if has_value(popmax_af):
                counts["nonmissing_gnomad_joint_popmax_af"] += 1
            if has_value(row.get("dbnsfp_gnomAD4.1_joint_nhomalt", "")):
                counts["nonmissing_gnomad_joint_nhomalt"] += 1
            if has_value(dbnsfp_popmax_af):
                counts["nonmissing_dbnsfp_popmax_af"] += 1
            writer.writerow(
                {
                    "variant_id": row["variant_id"],
                    "chrom": row["chrom"],
                    "pos": row["pos"],
                    "ref": row["ref"],
                    "alt": row["alt"],
                    "genes": row["genes"],
                    "sources": row["sources"],
                    "labels_3class": row["labels_3class"],
                    "gnomad_proxy_status": status,
                    "gnomad_proxy_missing_reason": missing_reason,
                    "gnomad_joint_af": joint_af,
                    "gnomad_joint_nhomalt": row.get("dbnsfp_gnomAD4.1_joint_nhomalt", ""),
                    "gnomad_joint_popmax_af": popmax_af,
                    "gnomad_joint_popmax_nhomalt": row.get("dbnsfp_gnomAD4.1_joint_POPMAX_nhomalt", ""),
                    "gnomad_dbnsfp_popmax_af": dbnsfp_popmax_af,
                    "gnomad_dbnsfp_popmax_pop": row.get("dbnsfp_dbNSFP_POPMAX_POP", ""),
                    "gnomad_proxy_source": "dbNSFP5.3.1a_GRCh38",
                }
            )

    summary = {
        "input": str(INPUT.relative_to(ROOT)),
        "output": str(OUTPUT.relative_to(ROOT)),
        "counts": dict(sorted(counts.items())),
        "notes": [
            "This table is the immediate local population-frequency feature layer.",
            "It uses gnomAD-related fields already present in dbNSFP v5.3.1a.",
            "It is not a full replacement for direct gnomAD browser/Hail/VCF extraction if filters, transcript consequences, or exhaustive benign sampling are needed.",
        ],
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
