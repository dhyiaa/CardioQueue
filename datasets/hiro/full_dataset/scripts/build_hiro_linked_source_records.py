#!/usr/bin/env python3
"""Build a HiRO patient/source-record table linked to variant-level registry rows.

The modeling matrix is variant-level, but HiRO evidence is source-record level:
the same genotype may appear in multiple participants with different phenotype
fields and adjudicated classifications. This table preserves every HiRO row and
adds the coordinate/registry link used for feature joins.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]

HIRO = ROOT / "datasets/hiro/full_dataset/model_inputs/ml_baseline_features.csv"
RESCUES = ROOT / "datasets/hiro/full_dataset/interim/hiro_coordinate_rescue_candidates.tsv"
REGISTRY = ROOT / "datasets/variant_registry/interim/clean_combined_variant_registry_with_hiro_rescues.tsv"

OUT = ROOT / "datasets/hiro/full_dataset/interim/hiro_linked_source_records.tsv"
SUMMARY = ROOT / "datasets/hiro/full_dataset/interim/hiro_linked_source_records.summary.json"

MISSING = {"", "Missing", "missing", "NA", "N/A", "nan", "None", "-", "."}


def clean(value: object) -> str:
    return str(value if value is not None else "").strip()


def present(value: object) -> bool:
    return clean(value) not in MISSING


def norm_chrom(value: object) -> str:
    text = clean(value)
    if text.lower().startswith("chr"):
        text = text[3:]
    return "MT" if text.upper() == "M" else text.upper() if text.upper() in {"X", "Y", "MT"} else text


def make_variant_id(chrom: object, pos: object, ref: object, alt: object) -> str:
    if not all(present(x) for x in [chrom, pos, ref, alt]):
        return ""
    return f"{norm_chrom(chrom)}-{clean(pos)}-{clean(ref).upper()}-{clean(alt).upper()}"


def load_rescues() -> dict[str, dict[str, str]]:
    by_uid = {}
    if not RESCUES.exists():
        return by_uid
    with RESCUES.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_uid[row["variant_uid"]] = row
    return by_uid


def load_registry_ids() -> set[str]:
    with REGISTRY.open(newline="") as handle:
        return {row["variant_id"] for row in csv.DictReader(handle, delimiter="\t")}


def resolved_link(row: dict[str, str], rescue: dict[str, str] | None, registry_ids: set[str]) -> dict[str, str]:
    original_id = make_variant_id(row.get("chrom"), row.get("pos"), row.get("ref"), row.get("alt"))
    rescue_id = ""
    if rescue and rescue.get("rescue_confidence") == "high":
        rescue_id = make_variant_id(
            rescue.get("rescued_chrom"),
            rescue.get("rescued_pos"),
            rescue.get("rescued_ref"),
            rescue.get("rescued_alt"),
        )
    resolved_id = original_id or rescue_id
    if resolved_id and resolved_id in registry_ids:
        status = "linked_to_registry"
    elif rescue_id:
        status = "rescued_but_not_in_promoted_registry"
    elif original_id:
        status = "coordinate_complete_but_not_in_registry"
    else:
        status = "unresolved_no_coordinates"
    return {
        "resolved_variant_id": resolved_id,
        "original_variant_id": original_id,
        "rescued_variant_id": rescue_id,
        "hiro_registry_link_status": status,
        "coordinate_rescue_status": rescue.get("rescue_status", "") if rescue else "",
        "coordinate_rescue_confidence": rescue.get("rescue_confidence", "") if rescue else "",
        "coordinate_rescue_source": rescue.get("rescue_source", "") if rescue else "",
        "coordinate_rescue_clinvar_variation_id": rescue.get("rescue_clinvar_variation_id", "") if rescue else "",
        "coordinate_rescue_note": rescue.get("rescue_note", "") if rescue else "",
    }


def main() -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    rescues = load_rescues()
    registry_ids = load_registry_ids()
    counts: Counter[str] = Counter()
    by_resolved_variant: defaultdict[str, list[dict[str, str]]] = defaultdict(list)

    with HIRO.open(newline="") as in_handle, OUT.open("w", newline="") as out_handle:
        reader = csv.DictReader(in_handle)
        added_fields = [
            "resolved_variant_id",
            "original_variant_id",
            "rescued_variant_id",
            "hiro_registry_link_status",
            "coordinate_rescue_status",
            "coordinate_rescue_confidence",
            "coordinate_rescue_source",
            "coordinate_rescue_clinvar_variation_id",
            "coordinate_rescue_note",
            "is_distinct_hiro_source_record",
        ]
        fieldnames = added_fields + list(reader.fieldnames or [])
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            counts["hiro_raw_rows"] += 1
            link = resolved_link(row, rescues.get(row["variant_uid"]), registry_ids)
            out_row = {**row, **link, "is_distinct_hiro_source_record": "true"}
            writer.writerow(out_row)
            counts[f"link_status_{link['hiro_registry_link_status']}"] += 1
            counts[f"label_{row.get('target_3class') or 'blank'}"] += 1
            if link["resolved_variant_id"]:
                by_resolved_variant[link["resolved_variant_id"]].append(out_row)

    duplicate_variant_groups = {k: v for k, v in by_resolved_variant.items() if len(v) > 1}
    discordant_label_groups = 0
    phenotype_variable_groups = 0
    for rows in duplicate_variant_groups.values():
        if len({r.get("target_3class", "") for r in rows}) > 1:
            discordant_label_groups += 1
        phenotype_signature = {
            (
                r.get("patient_sym_presyncope", ""),
                r.get("patient_sym_syncope", ""),
                r.get("patient_sym_palpitations", ""),
                r.get("patient_sym_chest_pain", ""),
                r.get("patient_sym_cardiac_arrest", ""),
                r.get("patient_sym_death", ""),
                r.get("patient_qt_latest", ""),
                r.get("patient_lvef_latest", ""),
            )
            for r in rows
        }
        if len(phenotype_signature) > 1:
            phenotype_variable_groups += 1

    summary = {
        "inputs": {
            "hiro": str(HIRO.relative_to(ROOT)),
            "rescue_candidates": str(RESCUES.relative_to(ROOT)),
            "promoted_registry": str(REGISTRY.relative_to(ROOT)),
        },
        "output": str(OUT.relative_to(ROOT)),
        "counts": dict(sorted(counts.items())),
        "unique_resolved_variant_ids": len(by_resolved_variant),
        "duplicate_resolved_variant_groups": len(duplicate_variant_groups),
        "duplicate_resolved_source_records": sum(len(v) for v in duplicate_variant_groups.values()),
        "duplicate_groups_with_discordant_3class_labels": discordant_label_groups,
        "duplicate_groups_with_variable_phenotype_signature": phenotype_variable_groups,
        "notes": [
            "Every HiRO raw row is preserved as a distinct source record.",
            "The variant-level modeling matrix may collapse records by resolved_variant_id for feature learning/leakage control.",
            "Use this table for patient-level phenotype, ACMG/adjudication review, and source-record stratified analyses.",
        ],
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
