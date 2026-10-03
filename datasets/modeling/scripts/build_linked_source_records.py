#!/usr/bin/env python3
"""Build source-record preservation tables for non-HiRO local datasets."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]

DATASETS = {
    "emerge": {
        "input": ROOT / "datasets/emerge/full_arrhythmia_gene_dataset/data/emerge_arrhythmia_ml_baseline_features.csv",
        "out": ROOT / "datasets/emerge/full_arrhythmia_gene_dataset/data/emerge_linked_source_records.tsv",
        "summary": ROOT / "datasets/emerge/full_arrhythmia_gene_dataset/provenance/emerge_linked_source_records.summary.json",
    },
    "cardioboost": {
        "input": ROOT
        / "datasets/cardioboost/public_dataset/processed_for_cardiogenetics/model_inputs/reannotated_ml_baseline_features.csv",
        "out": ROOT
        / "datasets/cardioboost/public_dataset/processed_for_cardiogenetics/model_inputs/cardioboost_linked_source_records.tsv",
        "summary": ROOT
        / "datasets/cardioboost/public_dataset/processed_for_cardiogenetics/provenance/cardioboost_linked_source_records.summary.json",
    },
}

REGISTRY = ROOT / "datasets/variant_registry/interim/clean_combined_variant_registry_with_hiro_rescues.tsv"

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


def make_variant_id(row: dict[str, str]) -> str:
    if not all(present(row.get(field)) for field in ["chrom", "pos", "ref", "alt"]):
        return ""
    return f"{norm_chrom(row.get('chrom'))}-{clean(row.get('pos'))}-{clean(row.get('ref')).upper()}-{clean(row.get('alt')).upper()}"


def load_registry_ids() -> set[str]:
    with REGISTRY.open(newline="") as handle:
        return {row["variant_id"] for row in csv.DictReader(handle, delimiter="\t")}


def phenotype_signature(row: dict[str, str]) -> tuple[str, ...]:
    fields = [
        "patient_sym_presyncope",
        "patient_sym_syncope",
        "patient_sym_palpitations",
        "patient_sym_chest_pain",
        "patient_sym_cardiac_arrest",
        "patient_sym_death",
        "patient_qt_latest",
        "patient_lvef_latest",
    ]
    return tuple(row.get(field, "") for field in fields)


def build_one(name: str, paths: dict[str, Path], registry_ids: set[str]) -> dict[str, object]:
    paths["out"].parent.mkdir(parents=True, exist_ok=True)
    paths["summary"].parent.mkdir(parents=True, exist_ok=True)
    counts: Counter[str] = Counter()
    by_variant: defaultdict[str, list[dict[str, str]]] = defaultdict(list)

    with paths["input"].open(newline="") as in_handle, paths["out"].open("w", newline="") as out_handle:
        reader = csv.DictReader(in_handle)
        added = [
            "resolved_variant_id",
            "source_registry_link_status",
            "is_distinct_source_record",
        ]
        fieldnames = added + list(reader.fieldnames or [])
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            counts["raw_rows"] += 1
            variant_id = make_variant_id(row)
            if not variant_id:
                status = "unresolved_no_coordinates"
            elif variant_id in registry_ids:
                status = "linked_to_registry"
            else:
                status = "coordinate_complete_but_not_in_registry"
            out_row = {
                **row,
                "resolved_variant_id": variant_id,
                "source_registry_link_status": status,
                "is_distinct_source_record": "true",
            }
            writer.writerow(out_row)
            counts[f"link_status_{status}"] += 1
            counts[f"label_{row.get('target_3class') or 'blank'}"] += 1
            if variant_id:
                by_variant[variant_id].append(out_row)

    duplicate_groups = {key: rows for key, rows in by_variant.items() if len(rows) > 1}
    discordant_labels = sum(
        1 for rows in duplicate_groups.values() if len({row.get("target_3class", "") for row in rows}) > 1
    )
    variable_phenotypes = sum(
        1 for rows in duplicate_groups.values() if len({phenotype_signature(row) for row in rows}) > 1
    )
    summary = {
        "dataset": name,
        "input": str(paths["input"].relative_to(ROOT)),
        "output": str(paths["out"].relative_to(ROOT)),
        "counts": dict(sorted(counts.items())),
        "unique_resolved_variant_ids": len(by_variant),
        "duplicate_resolved_variant_groups": len(duplicate_groups),
        "duplicate_resolved_source_records": sum(len(rows) for rows in duplicate_groups.values()),
        "duplicate_groups_with_discordant_3class_labels": discordant_labels,
        "duplicate_groups_with_variable_phenotype_signature": variable_phenotypes,
        "notes": [
            "Every input row is preserved as a distinct source record.",
            "The variant-level modeling matrix may still collapse records by resolved_variant_id for leakage control.",
        ],
    }
    paths["summary"].write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> int:
    registry_ids = load_registry_ids()
    summaries = {name: build_one(name, paths, registry_ids) for name, paths in DATASETS.items()}
    print(json.dumps(summaries, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
