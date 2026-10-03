#!/usr/bin/env python3
"""Audit ClinVar overlap without substituting ClinVar for external source labels."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FINAL = ROOT / "datasets/external_validation_candidates/final_model_heldout_2026"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def clinvar_binary_category(significance: str) -> str:
    value = significance.lower()
    if "conflict" in value or "uncertain" in value or value in {"-", "not provided"}:
        return "nonbinary_or_conflicting"
    if "pathogenic" in value and "benign" not in value:
        return "Pathogenic"
    if "benign" in value and "pathogenic" not in value:
        return "Benign"
    return "nonbinary_or_conflicting"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observations", type=Path, default=FINAL / "external_source_observations.tsv")
    parser.add_argument("--manifest", type=Path, default=FINAL / "external_unique_variant_manifest.tsv")
    parser.add_argument("--clinvar", type=Path, default=ROOT / "datasets/clinvar/variant_summary.txt.gz")
    parser.add_argument("--share", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=FINAL / "label_provenance_audit.json")
    args = parser.parse_args()

    target_sources = {"SHaRe_HCM_2026Q1", "Fernandez_Falgueras_PLOS_2024"}
    observations = [row for row in read_tsv(args.observations) if row["source"] in target_sources]
    target_ids = {row["variant_id"] for row in read_tsv(args.manifest)}

    clinvar_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    with gzip.open(args.clinvar, "rt", newline="", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("Assembly") != "GRCh38":
                continue
            chromosome = (row.get("Chromosome") or "").removeprefix("chr")
            position = row.get("PositionVCF") or row.get("Start")
            variant_id = f"{chromosome}-{position}-{row.get('ReferenceAlleleVCF')}-{row.get('AlternateAlleleVCF')}"
            if variant_id in target_ids:
                clinvar_rows[variant_id].append(row)

    summary = {
        "final_benchmark_current_clinvar_presence": {
            "total": len(target_ids),
            "exact_current_clinvar_records": len(set(clinvar_rows)),
            "absent_from_current_clinvar": len(target_ids - set(clinvar_rows)),
        }
    }
    for source in sorted(target_sources):
        source_rows = [row for row in observations if row["source"] == source]
        source_labels = {row["variant_id"]: row["external_label"] for row in source_rows}
        comparison = Counter()
        significance = Counter()
        for variant_id, source_label in source_labels.items():
            if variant_id not in clinvar_rows:
                comparison["absent_from_current_clinvar"] += 1
                continue
            categories = set()
            for row in clinvar_rows[variant_id]:
                significance[row["ClinicalSignificance"]] += 1
                categories.add(clinvar_binary_category(row["ClinicalSignificance"]))
            if categories == {source_label}:
                comparison["concordant_binary_clinvar"] += 1
            elif source_label in categories and "nonbinary_or_conflicting" in categories:
                comparison["concordant_plus_nonbinary_clinvar"] += 1
            else:
                comparison["discordant_or_nonbinary_clinvar"] += 1
        summary[source] = {
            "eligible_observations": len(source_rows),
            "exact_current_clinvar_records": sum(variant_id in clinvar_rows for variant_id in source_labels),
            "comparison_to_current_clinvar": dict(comparison),
            "current_clinvar_significance": dict(significance),
        }

    share_ids = {row["variant_id"] for row in observations if row["source"] == "SHaRe_HCM_2026Q1"}
    with args.share.open(newline="", encoding="utf-8-sig") as handle:
        share_rows = [row for row in csv.DictReader(handle) if row["VariantID"] in share_ids]
    summary["SHaRe_browser_internal_fields"] = {
        "matched_eligible_rows": len(share_rows),
        "varclass_status": dict(Counter((row.get("VarClass_status") or "missing") for row in share_rows)),
        "share_clinvar_conflict_flag": dict(Counter((row.get("Flag_conflicting_SHaRe_ClinVar_VarClass") or "missing") for row in share_rows)),
        "clinvar_browser_class": dict(Counter((row.get("ClinVar_VarClass") or "missing") for row in share_rows)),
        "interpretation": (
            "SHaRe stores SHaRe and ClinVar classifications in separate fields and explicitly flags disagreement; "
            "this supports a distinct SHaRe curation layer but not complete evidence independence."
        ),
    }
    summary["interpretation_limits"] = [
        "Current ClinVar presence is not equivalent to use in CardioQueue training.",
        "A distinct source label can still use literature or database evidence shared with ClinVar.",
        "Only PLOS/SHaRE alleles absent from the original train and validation partitions enter the final benchmark.",
        "ClinVar coordinate resolution for PLOS does not determine the PLOS outcome label.",
    ]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
