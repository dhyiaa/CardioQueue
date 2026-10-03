#!/usr/bin/env python3
"""Build a coordinate-level master variant registry.

The registry is a deduplicated GRCh38-oriented variant spine. It keeps one row
per chrom-pos-ref-alt and a companion long table that records every source row
contributing to each key.
"""

from __future__ import annotations

import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[3]

CORE_SOURCES = {
    "hiro": ROOT / "datasets/hiro/full_dataset/model_inputs/ml_baseline_features.csv",
    "emerge": ROOT / "datasets/emerge/full_arrhythmia_gene_dataset/data/emerge_arrhythmia_ml_baseline_features.csv",
    "cardioboost": ROOT
    / "datasets/cardioboost/public_dataset/processed_for_cardiogenetics/model_inputs/reannotated_ml_baseline_features.csv",
}
CLINVAR = ROOT / "datasets/clinvar/variant_summary.txt.gz"

OUT_DIR = ROOT / "datasets/variant_registry/interim"
REGISTRY_OUT = OUT_DIR / "master_variant_registry.tsv"
SOURCES_OUT = OUT_DIR / "master_variant_registry_sources.tsv"
SUMMARY_OUT = OUT_DIR / "master_variant_registry.summary.json"
GENE_PANEL_OUT = OUT_DIR / "master_variant_registry_gene_panel.txt"
SOURCE_INVENTORY_OUT = OUT_DIR / "master_variant_registry_source_inventory.tsv"

MISSING = {"", "Missing", "missing", "NA", "N/A", "na", "nan", "None", "-", "."}


def clean(value: object) -> str:
    return str(value if value is not None else "").strip()


def is_present(value: object) -> bool:
    return clean(value) not in MISSING


def norm_chrom(chrom: object) -> str:
    value = clean(chrom)
    if value.lower().startswith("chr"):
        value = value[3:]
    return value.upper() if value.upper() in {"X", "Y", "MT", "M"} else value


def norm_allele(allele: object) -> str:
    return clean(allele).upper()


def variant_key(chrom: object, pos: object, ref: object, alt: object) -> str | None:
    c = norm_chrom(chrom)
    p = clean(pos)
    r = norm_allele(ref)
    a = norm_allele(alt)
    if not all(is_present(x) for x in [c, p, r, a]):
        return None
    try:
        if int(p) <= 0:
            return None
    except ValueError:
        return None
    return f"{c}-{p}-{r}-{a}"


def split_genes(raw: object) -> list[str]:
    text = clean(raw).replace(";", ",").replace("|", ",")
    genes = []
    for piece in text.split(","):
        gene = piece.strip().upper()
        if gene and gene not in MISSING:
            genes.append(gene)
    return genes


def label_3class(raw: object) -> str:
    text = clean(raw).lower()
    if not text:
        return ""
    if "conflict" in text or "drug response" in text or "association" in text:
        return "Conflict_or_other"
    has_path = "pathogenic" in text
    has_benign = "benign" in text
    has_vus = "uncertain" in text or "vus" in text
    if has_path and has_benign:
        return "Conflict_or_other"
    if has_path:
        return "Pathogenic"
    if has_benign:
        return "Benign"
    if has_vus:
        return "VUS"
    return "Other"


def clinvar_stars(review_status: object) -> str:
    text = clean(review_status).lower()
    if "practice guideline" in text:
        return "4"
    if "reviewed by expert panel" in text:
        return "3"
    if "multiple submitters" in text and "no conflicts" in text:
        return "2"
    if "criteria provided" in text:
        return "1"
    if "no assertion" in text or "no classification" in text:
        return "0"
    return ""


def load_gene_panel() -> set[str]:
    genes: set[str] = set()
    panel_files = [
        ROOT / "datasets/feature_sources/gnomad_v4/external/gene_panel.example.txt",
        ROOT / "datasets/cardioboost/public_dataset/raw_public_cardioBoost_repo/data/arrhythmia/arm_gene.txt",
        ROOT / "datasets/cardioboost/public_dataset/raw_public_cardioBoost_repo/data/cardiomyopathy/cm_gene.txt",
    ]
    for path in panel_files:
        if not path.exists():
            continue
        with path.open(newline="") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                gene = line.split("\t")[0].split(",")[0].strip().upper()
                if gene and gene != "GENE_SYMBOL":
                    genes.add(gene)

    for path in CORE_SOURCES.values():
        if not path.exists():
            continue
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                genes.update(split_genes(row.get("gene", "")))

    # Remove artifacts from compound gene strings that are not standard symbols
    # for this project's main panel.
    return {g for g in genes if g and not g.startswith("FPGT-")}


def add_member(
    records: dict[str, dict[str, object]],
    members: list[dict[str, str]],
    *,
    key: str,
    source: str,
    source_row_id: str,
    gene: str,
    label: str,
    raw_label: str,
    variant_uid: str,
    variant_key_source: str,
    assembly: str,
    review_status: str = "",
    review_stars: str = "",
    phenotype_available: str = "false",
) -> None:
    chrom, pos, ref, alt = key.split("-", 3)
    rec = records.setdefault(
        key,
        {
            "variant_registry_key": key,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "assembly": assembly,
            "genes": set(),
            "sources": set(),
            "labels_3class": set(),
            "raw_labels": set(),
            "review_statuses": set(),
            "max_clinvar_stars": "",
            "phenotype_available": False,
            "source_variant_count": 0,
            "source_row_count": 0,
        },
    )
    rec["genes"].update(split_genes(gene))  # type: ignore[union-attr]
    rec["sources"].add(source)  # type: ignore[union-attr]
    if label:
        rec["labels_3class"].add(label)  # type: ignore[union-attr]
    if raw_label:
        rec["raw_labels"].add(raw_label)  # type: ignore[union-attr]
    if review_status:
        rec["review_statuses"].add(review_status)  # type: ignore[union-attr]
    if phenotype_available == "true":
        rec["phenotype_available"] = True
    if review_stars:
        current = rec.get("max_clinvar_stars", "")
        if not current or int(review_stars) > int(str(current)):
            rec["max_clinvar_stars"] = review_stars
    rec["source_row_count"] = int(rec["source_row_count"]) + 1

    members.append(
        {
            "variant_registry_key": key,
            "source": source,
            "source_row_id": source_row_id,
            "source_variant_uid": variant_uid,
            "source_variant_key": variant_key_source,
            "gene": gene,
            "label_3class": label,
            "raw_label": raw_label,
            "assembly": assembly,
            "review_status": review_status,
            "review_stars": review_stars,
            "phenotype_available": phenotype_available,
        }
    )


def iter_core_rows(source: str, path: Path) -> Iterable[dict[str, str]]:
    with path.open(newline="") as handle:
        yield from csv.DictReader(handle)


def add_core_sources(records: dict[str, dict[str, object]], members: list[dict[str, str]]) -> Counter:
    counts: Counter = Counter()
    for source, path in CORE_SOURCES.items():
        for idx, row in enumerate(iter_core_rows(source, path), start=1):
            counts[f"{source}_rows_seen"] += 1
            key = variant_key(row.get("chrom"), row.get("pos"), row.get("ref"), row.get("alt"))
            if not key:
                counts[f"{source}_rows_without_coordinates"] += 1
                continue
            counts[f"{source}_rows_with_coordinates"] += 1
            label = clean(row.get("target_3class", ""))
            add_member(
                records,
                members,
                key=key,
                source=source,
                source_row_id=str(idx),
                gene=clean(row.get("gene", "")),
                label=label,
                raw_label=clean(row.get("target_classification_raw", "")),
                variant_uid=clean(row.get("variant_uid", "")),
                variant_key_source=clean(row.get("variant_key", "")),
                assembly="GRCh38_assumed_from_processed_table",
                phenotype_available="true" if source == "hiro" else "false",
            )
    return counts


def write_source_inventory(gene_panel: set[str]) -> list[dict[str, str]]:
    inventory: list[dict[str, str]] = []
    for source, path in CORE_SOURCES.items():
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        coordinate_keys = []
        rows_with_coordinates = 0
        rows_without_coordinates = 0
        label_counts = Counter()
        for row in rows:
            label_counts[clean(row.get("target_3class", "")) or "blank"] += 1
            key = variant_key(row.get("chrom"), row.get("pos"), row.get("ref"), row.get("alt"))
            if key:
                rows_with_coordinates += 1
                coordinate_keys.append(key)
            else:
                rows_without_coordinates += 1
        inventory.append(
            {
                "source": source,
                "source_table": str(path.relative_to(ROOT)),
                "source_rows_total": str(len(rows)),
                "rows_with_complete_coordinate_key": str(rows_with_coordinates),
                "rows_without_complete_coordinate_key": str(rows_without_coordinates),
                "unique_coordinate_keys": str(len(set(coordinate_keys))),
                "private_or_patient_linked": "true" if source == "hiro" else "false",
                "phenotype_available": "true" if source == "hiro" else "false",
                "label_counts": json.dumps(dict(sorted(label_counts.items())), sort_keys=True),
                "note": (
                    "HiRO source rows are private patient-linked observations; coordinate registry counts are deduplicated variants, not patient/row counts."
                    if source == "hiro"
                    else "Coordinate registry counts are deduplicated variants, not necessarily source-row counts."
                ),
            }
        )

    clinvar_rows = sum(1 for _ in gzip.open(CLINVAR, "rt")) - 1
    inventory.append(
        {
            "source": "clinvar",
            "source_table": str(CLINVAR.relative_to(ROOT)),
            "source_rows_total": str(clinvar_rows),
            "rows_with_complete_coordinate_key": "",
            "rows_without_complete_coordinate_key": "",
            "unique_coordinate_keys": "",
            "private_or_patient_linked": "false",
            "phenotype_available": "false",
            "label_counts": "",
            "note": f"Registry includes GRCh38 rows intersecting first-pass {len(gene_panel)}-gene panel, not all ClinVar rows.",
        }
    )

    fields = [
        "source",
        "source_table",
        "source_rows_total",
        "rows_with_complete_coordinate_key",
        "rows_without_complete_coordinate_key",
        "unique_coordinate_keys",
        "private_or_patient_linked",
        "phenotype_available",
        "label_counts",
        "note",
    ]
    with SOURCE_INVENTORY_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(inventory)
    return inventory


def add_clinvar(
    records: dict[str, dict[str, object]],
    members: list[dict[str, str]],
    gene_panel: set[str],
) -> Counter:
    counts: Counter = Counter()
    with gzip.open(CLINVAR, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for idx, row in enumerate(reader, start=1):
            counts["clinvar_rows_seen"] += 1
            if row.get("Assembly") != "GRCh38":
                continue
            genes = split_genes(row.get("GeneSymbol", ""))
            if not any(g in gene_panel for g in genes):
                continue
            key = variant_key(
                row.get("Chromosome"),
                row.get("PositionVCF"),
                row.get("ReferenceAlleleVCF"),
                row.get("AlternateAlleleVCF"),
            )
            if not key:
                counts["clinvar_gene_panel_rows_without_coordinates"] += 1
                continue
            raw = clean(row.get("ClinicalSignificance", ""))
            label = label_3class(raw)
            counts["clinvar_gene_panel_rows_with_coordinates"] += 1
            counts[f"clinvar_label_{label or 'blank'}"] += 1
            add_member(
                records,
                members,
                key=key,
                source="clinvar",
                source_row_id=clean(row.get("VariationID", "")) or str(idx),
                gene=";".join(genes),
                label=label,
                raw_label=raw,
                variant_uid=clean(row.get("VariationID", "")),
                variant_key_source=clean(row.get("Name", "")),
                assembly="GRCh38",
                review_status=clean(row.get("ReviewStatus", "")),
                review_stars=clinvar_stars(row.get("ReviewStatus", "")),
                phenotype_available="false",
            )
    return counts


def serialize_set(value: object) -> str:
    if isinstance(value, set):
        return "|".join(sorted(str(v) for v in value if str(v)))
    return clean(value)


def write_outputs(
    records: dict[str, dict[str, object]],
    members: list[dict[str, str]],
    gene_panel: set[str],
    counts: Counter,
) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    GENE_PANEL_OUT.write_text("\n".join(sorted(gene_panel)) + "\n")
    inventory = write_source_inventory(gene_panel)

    source_names = ["clinvar", "hiro", "emerge", "cardioboost"]
    registry_fields = [
        "variant_registry_key",
        "chrom",
        "pos",
        "ref",
        "alt",
        "assembly",
        "genes",
        "sources",
        "source_count",
        "source_row_count",
        "in_clinvar",
        "in_hiro",
        "in_emerge",
        "in_cardioboost",
        "phenotype_available",
        "labels_3class",
        "raw_labels",
        "review_statuses",
        "max_clinvar_stars",
    ]
    with REGISTRY_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=registry_fields, delimiter="\t")
        writer.writeheader()
        for key in sorted(records, key=lambda k: (norm_chrom(k.split("-")[0]), int(k.split("-")[1]), k)):
            rec = records[key]
            sources = rec["sources"]
            row = {field: "" for field in registry_fields}
            for field in ["variant_registry_key", "chrom", "pos", "ref", "alt", "assembly", "source_row_count"]:
                row[field] = serialize_set(rec.get(field, ""))
            row["genes"] = serialize_set(rec["genes"])
            row["sources"] = serialize_set(sources)
            row["source_count"] = str(len(sources))  # type: ignore[arg-type]
            for source in source_names:
                row[f"in_{source}"] = "true" if source in sources else "false"  # type: ignore[operator]
            row["phenotype_available"] = "true" if rec.get("phenotype_available") else "false"
            row["labels_3class"] = serialize_set(rec["labels_3class"])
            row["raw_labels"] = serialize_set(rec["raw_labels"])
            row["review_statuses"] = serialize_set(rec["review_statuses"])
            row["max_clinvar_stars"] = serialize_set(rec.get("max_clinvar_stars", ""))
            writer.writerow(row)

    source_fields = [
        "variant_registry_key",
        "source",
        "source_row_id",
        "source_variant_uid",
        "source_variant_key",
        "gene",
        "label_3class",
        "raw_label",
        "assembly",
        "review_status",
        "review_stars",
        "phenotype_available",
    ]
    with SOURCES_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=source_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(members)

    source_combo_counts = Counter(serialize_set(rec["sources"]) for rec in records.values())
    label_counts = Counter()
    for rec in records.values():
        labels = serialize_set(rec["labels_3class"]) or "unlabeled"
        label_counts[labels] += 1
    summary = {
        "registry_rows": len(records),
        "source_membership_rows": len(members),
        "gene_panel_size": len(gene_panel),
        "outputs": {
            "registry": str(REGISTRY_OUT.relative_to(ROOT)),
            "sources": str(SOURCES_OUT.relative_to(ROOT)),
            "gene_panel": str(GENE_PANEL_OUT.relative_to(ROOT)),
            "source_inventory": str(SOURCE_INVENTORY_OUT.relative_to(ROOT)),
        },
        "source_inventory": inventory,
        "source_counts": dict(counts),
        "source_combo_counts": dict(source_combo_counts),
        "registry_label_counts": dict(label_counts),
    }
    SUMMARY_OUT.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


def main() -> None:
    gene_panel = load_gene_panel()
    records: dict[str, dict[str, object]] = {}
    members: list[dict[str, str]] = []
    counts = Counter()
    counts.update(add_core_sources(records, members))
    counts.update(add_clinvar(records, members, gene_panel))
    write_outputs(records, members, gene_panel, counts)
    print(json.dumps(json.loads(SUMMARY_OUT.read_text()), indent=2))


if __name__ == "__main__":
    main()
