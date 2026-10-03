#!/usr/bin/env python3
"""Create clean registry outputs with high-confidence HiRO rescues promoted."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]

REGISTRY_IN = ROOT / "datasets/variant_registry/interim/clean_combined_variant_registry.tsv"
SOURCES_IN = ROOT / "datasets/variant_registry/interim/clean_combined_variant_registry_sources.tsv"
RESCUES = ROOT / "datasets/hiro/full_dataset/interim/hiro_coordinate_rescue_candidates.tsv"
HIRO_RAW = ROOT / "datasets/hiro/full_dataset/model_inputs/ml_baseline_features.csv"
GENE_DECISIONS = ROOT / "datasets/gene_panels/interim/gene_panel_decisions.tsv"

OUT_DIR = ROOT / "datasets/variant_registry/interim"
REGISTRY_OUT = OUT_DIR / "clean_combined_variant_registry_with_hiro_rescues.tsv"
SOURCES_OUT = OUT_DIR / "clean_combined_variant_registry_sources_with_hiro_rescues.tsv"
PROMOTED_OUT = ROOT / "datasets/hiro/full_dataset/interim/hiro_coordinate_rescues_promoted.tsv"
SUMMARY_OUT = OUT_DIR / "clean_combined_variant_registry_with_hiro_rescues.summary.json"

SOURCES = ("clinvar", "hiro", "emerge", "cardioboost")


def split_pipe(value: str | None) -> set[str]:
    return {piece for piece in str(value or "").split("|") if piece}


def joined(values: set[str]) -> str:
    return "|".join(sorted(v for v in values if v))


def norm_chrom(value: str) -> str:
    text = str(value or "").strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    return "MT" if text.upper() == "M" else text.upper() if text.upper() in {"X", "Y", "MT"} else text


def variant_id(chrom: str, pos: str, ref: str, alt: str) -> str:
    return f"{norm_chrom(chrom)}-{str(pos).strip()}-{str(ref).strip().upper()}-{str(alt).strip().upper()}"


def sortable_key(key: str) -> tuple[tuple[int, str], int, str]:
    chrom, pos, *_ = key.split("-")
    order = {"X": 23, "Y": 24, "MT": 25}
    chrom_sort = (int(chrom), chrom) if chrom.isdigit() else (order.get(chrom, 99), chrom)
    return chrom_sort, int(pos), key


def load_gene_decisions() -> dict[str, str]:
    decisions = {}
    with GENE_DECISIONS.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            decisions[row["gene"]] = row.get("decision", "")
    return decisions


def gene_decision_for(gene: str, decisions: dict[str, str]) -> str:
    genes = [piece.strip().upper() for piece in gene.replace(",", ";").replace("|", ";").split(";")]
    return joined({decisions.get(g, "") for g in genes if g})


def load_registry() -> tuple[list[str], dict[str, dict[str, object]]]:
    records: dict[str, dict[str, object]] = {}
    with REGISTRY_IN.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        for row in reader:
            records[row["variant_id"]] = {
                "variant_id": row["variant_id"],
                "chrom": row["chrom"],
                "pos": row["pos"],
                "ref": row["ref"],
                "alt": row["alt"],
                "genes": split_pipe(row.get("genes")),
                "sources": split_pipe(row.get("sources")),
                "source_row_count": int(row.get("source_row_count") or 0),
                "phenotype_available": row.get("phenotype_available") == "true",
                "private_or_patient_linked": row.get("private_or_patient_linked") == "true",
                "labels_3class": split_pipe(row.get("labels_3class")),
                "raw_labels": split_pipe(row.get("raw_labels")),
                "review_stars": split_pipe(row.get("max_review_stars")),
                "review_statuses": split_pipe(row.get("review_statuses")),
                "confidence_tiers": split_pipe(row.get("confidence_tiers")),
                "gene_panel_decisions": split_pipe(row.get("gene_panel_decisions")),
            }
    return fields, records


def load_hiro_raw() -> dict[str, dict[str, str]]:
    by_uid = {}
    with HIRO_RAW.open(newline="") as handle:
        for row in csv.DictReader(handle):
            by_uid[row["variant_uid"]] = row
    return by_uid


def add_member_source_fields(fields: list[str]) -> list[str]:
    extra = [
        "coordinate_rescue_status",
        "coordinate_rescue_confidence",
        "coordinate_rescue_source",
        "coordinate_rescue_clinvar_variation_id",
        "coordinate_rescue_note",
    ]
    return fields + [field for field in extra if field not in fields]


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    registry_fields, records = load_registry()
    decisions = load_gene_decisions()
    hiro_raw = load_hiro_raw()
    counts: Counter[str] = Counter()

    with SOURCES_IN.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        source_fields = add_member_source_fields(list(reader.fieldnames or []))
        members = [{field: row.get(field, "") for field in source_fields} for row in reader]

    existing_hiro_uids = {
        row.get("source_variant_uid", "")
        for row in members
        if row.get("source") == "hiro" and row.get("source_variant_uid")
    }

    promoted_rows: list[dict[str, str]] = []
    with RESCUES.open(newline="") as handle:
        for rescue in csv.DictReader(handle, delimiter="\t"):
            counts["rescue_candidates_seen"] += 1
            if rescue.get("rescue_confidence") != "high":
                counts["rescue_candidates_not_promoted_non_high_confidence"] += 1
                continue

            uid = rescue.get("variant_uid", "")
            if uid in existing_hiro_uids:
                counts["rescue_candidates_skipped_existing_hiro_uid"] += 1
                continue

            key = variant_id(
                rescue.get("rescued_chrom", ""),
                rescue.get("rescued_pos", ""),
                rescue.get("rescued_ref", ""),
                rescue.get("rescued_alt", ""),
            )
            raw = hiro_raw.get(uid, {})
            gene = rescue.get("gene", "") or raw.get("gene", "")
            label = raw.get("target_3class", "")
            raw_label = raw.get("target_classification_raw", "")
            decision = gene_decision_for(gene, decisions)

            if key not in records:
                chrom, pos, ref, alt = key.split("-", 3)
                records[key] = {
                    "variant_id": key,
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "genes": set(),
                    "sources": set(),
                    "source_row_count": 0,
                    "phenotype_available": False,
                    "private_or_patient_linked": False,
                    "labels_3class": set(),
                    "raw_labels": set(),
                    "review_stars": set(),
                    "review_statuses": set(),
                    "confidence_tiers": set(),
                    "gene_panel_decisions": set(),
                }
                counts["new_variant_rows_added"] += 1
            else:
                counts["existing_variant_rows_enriched"] += 1

            rec = records[key]
            rec["genes"].update(split_pipe(gene.replace(";", "|")))  # type: ignore[union-attr]
            rec["sources"].add("hiro")  # type: ignore[union-attr]
            rec["source_row_count"] = int(rec["source_row_count"]) + 1
            rec["phenotype_available"] = True
            rec["private_or_patient_linked"] = True
            if label:
                rec["labels_3class"].add(label)  # type: ignore[union-attr]
            if raw_label:
                rec["raw_labels"].add(raw_label)  # type: ignore[union-attr]
            if decision:
                rec["gene_panel_decisions"].update(split_pipe(decision.replace(";", "|")))  # type: ignore[union-attr]

            member = {field: "" for field in source_fields}
            member.update(
                {
                    "variant_id": key,
                    "source": "hiro",
                    "source_row_id": rescue.get("source_row", ""),
                    "source_variant_uid": uid,
                    "source_variant_key": rescue.get("variant_key", ""),
                    "gene": gene,
                    "label_3class": label,
                    "raw_label": raw_label,
                    "phenotype_available": "true",
                    "private_or_patient_linked": "true",
                    "gene_panel_decision": decision,
                    "coordinate_rescue_status": rescue.get("rescue_status", ""),
                    "coordinate_rescue_confidence": rescue.get("rescue_confidence", ""),
                    "coordinate_rescue_source": rescue.get("rescue_source", ""),
                    "coordinate_rescue_clinvar_variation_id": rescue.get("rescue_clinvar_variation_id", ""),
                    "coordinate_rescue_note": rescue.get("rescue_note", ""),
                }
            )
            members.append(member)

            promoted = dict(rescue)
            promoted["promoted_variant_id"] = key
            promoted["promoted_label_3class"] = label
            promoted["promoted_raw_label"] = raw_label
            promoted["promoted_gene_panel_decision"] = decision
            promoted_rows.append(promoted)
            counts["high_confidence_rescues_promoted"] += 1
            counts[f"promoted_label_{label or 'blank'}"] += 1

    with REGISTRY_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=registry_fields, delimiter="\t")
        writer.writeheader()
        for key in sorted(records, key=sortable_key):
            rec = records[key]
            sources = rec["sources"]
            stars = [int(x) for x in rec["review_stars"] if str(x).isdigit()]  # type: ignore[index]
            row = {
                "variant_id": key,
                "chrom": rec["chrom"],
                "pos": rec["pos"],
                "ref": rec["ref"],
                "alt": rec["alt"],
                "genes": joined(rec["genes"]),  # type: ignore[arg-type]
                "sources": joined(sources),  # type: ignore[arg-type]
                "source_count": str(len(sources)),  # type: ignore[arg-type]
                "source_row_count": str(rec["source_row_count"]),
                "phenotype_available": "true" if rec["phenotype_available"] else "false",
                "private_or_patient_linked": "true" if rec["private_or_patient_linked"] else "false",
                "labels_3class": joined(rec["labels_3class"]),  # type: ignore[arg-type]
                "raw_labels": joined(rec["raw_labels"]),  # type: ignore[arg-type]
                "max_review_stars": str(max(stars)) if stars else "",
                "review_statuses": joined(rec["review_statuses"]),  # type: ignore[arg-type]
                "confidence_tiers": joined(rec["confidence_tiers"]),  # type: ignore[arg-type]
                "gene_panel_decisions": joined(rec["gene_panel_decisions"]),  # type: ignore[arg-type]
            }
            for source in SOURCES:
                row[f"in_{source}"] = "true" if source in sources else "false"  # type: ignore[operator]
            writer.writerow(row)

    with SOURCES_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=source_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(members)

    promoted_fields = []
    for row in promoted_rows:
        for field in row:
            if field not in promoted_fields:
                promoted_fields.append(field)
    with PROMOTED_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=promoted_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(promoted_rows)

    summary = {
        "inputs": {
            "registry": str(REGISTRY_IN.relative_to(ROOT)),
            "sources": str(SOURCES_IN.relative_to(ROOT)),
            "rescue_candidates": str(RESCUES.relative_to(ROOT)),
            "hiro_raw": str(HIRO_RAW.relative_to(ROOT)),
        },
        "outputs": {
            "registry_with_hiro_rescues": str(REGISTRY_OUT.relative_to(ROOT)),
            "sources_with_hiro_rescues": str(SOURCES_OUT.relative_to(ROOT)),
            "promoted_rescues": str(PROMOTED_OUT.relative_to(ROOT)),
            "summary": str(SUMMARY_OUT.relative_to(ROOT)),
        },
        "original_registry_rows": sum(1 for _ in REGISTRY_IN.open()) - 1,
        "promoted_registry_rows": len(records),
        "original_source_rows": sum(1 for _ in SOURCES_IN.open()) - 1,
        "promoted_source_rows": len(members),
        "counts": dict(sorted(counts.items())),
        "notes": [
            "Only rescue_confidence=high rows are promoted.",
            "The original clean registry and source-membership files are unchanged.",
            "Promoted HiRO rows are marked phenotype_available=true and private_or_patient_linked=true.",
            "Coordinate provenance is preserved in the promoted source-membership file and promoted-rescues audit table.",
        ],
    }
    SUMMARY_OUT.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
