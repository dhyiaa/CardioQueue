#!/usr/bin/env python3
"""Promote validated eMERGE liftover coordinates into rescue-aware registry files."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]

REGISTRY_IN = ROOT / "datasets/variant_registry/interim/clean_combined_variant_registry_with_hiro_rescues.tsv"
SOURCES_IN = ROOT / "datasets/variant_registry/interim/clean_combined_variant_registry_sources_with_hiro_rescues.tsv"
EMERGE_RESCUE = ROOT / "datasets/emerge/full_arrhythmia_gene_dataset/interim/emerge_coordinate_rescue_liftover.tsv"

OUT_DIR = ROOT / "datasets/variant_registry/interim"
REGISTRY_OUT = OUT_DIR / "clean_combined_variant_registry_with_hiro_emerge_rescues.tsv"
SOURCES_OUT = OUT_DIR / "clean_combined_variant_registry_sources_with_hiro_emerge_rescues.tsv"
SUMMARY_OUT = OUT_DIR / "clean_combined_variant_registry_with_hiro_emerge_rescues.summary.json"

SOURCES = ("clinvar", "hiro", "emerge", "cardioboost")


def split_pipe(value: str | None) -> set[str]:
    return {piece for piece in str(value or "").split("|") if piece}


def joined(values: set[str]) -> str:
    return "|".join(sorted(v for v in values if v))


def sortable_key(key: str) -> tuple[tuple[int, str], int, str]:
    chrom, pos, *_ = key.split("-")
    order = {"X": 23, "Y": 24, "MT": 25}
    chrom_sort = (int(chrom), chrom) if chrom.isdigit() else (order.get(chrom, 99), chrom)
    try:
        pos_i = int(pos)
    except ValueError:
        pos_i = 10**12
    return chrom_sort, pos_i, key


def load_rescue() -> dict[str, dict[str, str]]:
    with EMERGE_RESCUE.open(newline="") as handle:
        return {row["variant_uid"]: row for row in csv.DictReader(handle, delimiter="\t")}


def aggregate_registry(members: list[dict[str, str]], registry_fields: list[str]) -> dict[str, dict[str, str]]:
    grouped: dict[str, dict[str, object]] = defaultdict(
        lambda: {
            "genes": set(),
            "sources": set(),
            "labels_3class": set(),
            "raw_labels": set(),
            "review_stars": set(),
            "review_statuses": set(),
            "confidence_tiers": set(),
            "gene_panel_decisions": set(),
            "phenotype_available": False,
            "private_or_patient_linked": False,
            "source_row_count": 0,
        }
    )
    for member in members:
        vid = member["variant_id"]
        if not vid:
            continue
        chrom, pos, ref, alt = vid.split("-", 3)
        rec = grouped[vid]
        rec.update({"variant_id": vid, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt})
        for field, target in [
            ("gene", "genes"),
            ("source", "sources"),
            ("label_3class", "labels_3class"),
            ("raw_label", "raw_labels"),
            ("review_stars", "review_stars"),
            ("review_status", "review_statuses"),
            ("confidence_tier", "confidence_tiers"),
            ("gene_panel_decision", "gene_panel_decisions"),
        ]:
            text = member.get(field, "")
            if text:
                for piece in text.replace(";", "|").split("|"):
                    piece = piece.strip()
                    if piece:
                        rec[target].add(piece)  # type: ignore[index]
        if member.get("phenotype_available") == "true":
            rec["phenotype_available"] = True
        if member.get("private_or_patient_linked") == "true":
            rec["private_or_patient_linked"] = True
        rec["source_row_count"] = int(rec["source_row_count"]) + 1

    out: dict[str, dict[str, str]] = {}
    for vid, rec in grouped.items():
        row = {field: "" for field in registry_fields}
        row.update(
            {
                "variant_id": vid,
                "chrom": str(rec["chrom"]),
                "pos": str(rec["pos"]),
                "ref": str(rec["ref"]),
                "alt": str(rec["alt"]),
                "genes": joined(rec["genes"]),  # type: ignore[arg-type]
                "sources": joined(rec["sources"]),  # type: ignore[arg-type]
                "source_count": str(len(rec["sources"])),  # type: ignore[arg-type]
                "source_row_count": str(rec["source_row_count"]),
                "phenotype_available": "true" if rec["phenotype_available"] else "false",
                "private_or_patient_linked": "true" if rec["private_or_patient_linked"] else "false",
                "labels_3class": joined(rec["labels_3class"]),  # type: ignore[arg-type]
                "raw_labels": joined(rec["raw_labels"]),  # type: ignore[arg-type]
                "review_statuses": joined(rec["review_statuses"]),  # type: ignore[arg-type]
                "confidence_tiers": joined(rec["confidence_tiers"]),  # type: ignore[arg-type]
                "gene_panel_decisions": joined(rec["gene_panel_decisions"]),  # type: ignore[arg-type]
            }
        )
        stars = [int(x) for x in rec["review_stars"] if str(x).isdigit()]  # type: ignore[index]
        row["max_review_stars"] = str(max(stars)) if stars else ""
        for source in SOURCES:
            row[f"in_{source}"] = "true" if source in rec["sources"] else "false"  # type: ignore[operator]
        out[vid] = row
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rescue = load_rescue()
    counts: Counter[str] = Counter()

    with REGISTRY_IN.open(newline="") as handle:
        registry_fields = list(csv.DictReader(handle, delimiter="\t").fieldnames or [])

    with SOURCES_IN.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        source_fields = list(reader.fieldnames or [])
        extra_fields = [
            "pre_emerge_liftover_variant_id",
            "emerge_liftover_status",
            "emerge_liftover_chain_strand",
            "emerge_grch38_ref_at_pos",
            "emerge_allele_ref_validation",
            "emerge_coordinate_qc_flag",
            "emerge_source_hg19_variant_id",
        ]
        for field in extra_fields:
            if field not in source_fields:
                source_fields.append(field)

        members: list[dict[str, str]] = []
        for row in reader:
            row = {field: row.get(field, "") for field in source_fields}
            counts["source_rows_seen"] += 1
            if row.get("source") == "emerge":
                counts["emerge_source_rows_seen"] += 1
                rec = rescue.get(row.get("source_variant_uid", ""))
                row["pre_emerge_liftover_variant_id"] = row.get("variant_id", "")
                if rec and rec.get("coordinate_qc_flag") == "pass":
                    row["variant_id"] = rec["grch38_variant_id"]
                    row["emerge_liftover_status"] = rec.get("liftover_status", "")
                    row["emerge_liftover_chain_strand"] = rec.get("liftover_chain_strand", "")
                    row["emerge_grch38_ref_at_pos"] = rec.get("grch38_ref_at_pos", "")
                    row["emerge_allele_ref_validation"] = rec.get("allele_ref_validation", "")
                    row["emerge_coordinate_qc_flag"] = rec.get("coordinate_qc_flag", "")
                    row["emerge_source_hg19_variant_id"] = rec.get("source_hg19_variant_id", "")
                    counts["emerge_source_rows_promoted"] += 1
                else:
                    row["emerge_coordinate_qc_flag"] = rec.get("coordinate_qc_flag", "missing_rescue") if rec else "missing_rescue"
                    counts["emerge_source_rows_not_promoted"] += 1
            members.append(row)

    registry = aggregate_registry(members, registry_fields)

    with REGISTRY_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=registry_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for vid in sorted(registry, key=sortable_key):
            writer.writerow(registry[vid])

    with SOURCES_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=source_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(members)

    summary = {
        "inputs": {
            "registry": str(REGISTRY_IN.relative_to(ROOT)),
            "sources": str(SOURCES_IN.relative_to(ROOT)),
            "emerge_coordinate_rescue": str(EMERGE_RESCUE.relative_to(ROOT)),
        },
        "outputs": {
            "registry": str(REGISTRY_OUT.relative_to(ROOT)),
            "sources": str(SOURCES_OUT.relative_to(ROOT)),
            "summary": str(SUMMARY_OUT.relative_to(ROOT)),
        },
        "input_registry_rows": sum(1 for _ in REGISTRY_IN.open()) - 1,
        "output_registry_rows": len(registry),
        "input_source_rows": counts["source_rows_seen"],
        "output_source_rows": len(members),
        "counts": dict(sorted(counts.items())),
        "notes": [
            "All eMERGE source membership rows are remapped to validated GRCh38 liftover IDs when coordinate_qc_flag == pass.",
            "Registry rows are rebuilt from the source-membership table so corrected eMERGE variants can merge with existing ClinVar/HiRO/CardioBoost rows if coordinates overlap.",
            "Old eMERGE variant IDs are preserved in pre_emerge_liftover_variant_id and emerge_source_hg19_variant_id.",
        ],
    }
    SUMMARY_OUT.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
