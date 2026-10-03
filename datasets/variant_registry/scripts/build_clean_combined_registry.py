#!/usr/bin/env python3
"""Build the clean combined modeling registry.

This registry is narrower than master_variant_registry.tsv. It uses the
post-gene-QC ClinVar slice plus coordinate-complete HiRO, eMERGE, and
CardioBoost rows. Raw/private source rows are not deleted; rows without a
complete coordinate key remain tracked in the summary.
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]

CLINVAR_PRIMARY = ROOT / "datasets/clinvar/interim/clinvar_primary_after_gene_qc.tsv"
CORE_SOURCES = {
    "hiro": ROOT / "datasets/hiro/full_dataset/model_inputs/ml_baseline_features.csv",
    "emerge": ROOT / "datasets/emerge/full_arrhythmia_gene_dataset/data/emerge_arrhythmia_ml_baseline_features.csv",
    "cardioboost": ROOT
    / "datasets/cardioboost/public_dataset/processed_for_cardiogenetics/model_inputs/reannotated_ml_baseline_features.csv",
}
GENE_DECISIONS = ROOT / "datasets/gene_panels/interim/gene_panel_decisions.tsv"

OUT_DIR = ROOT / "datasets/variant_registry/interim"
REGISTRY_OUT = OUT_DIR / "clean_combined_variant_registry.tsv"
SOURCES_OUT = OUT_DIR / "clean_combined_variant_registry_sources.tsv"
SUMMARY_OUT = OUT_DIR / "clean_combined_variant_registry.summary.json"

MISSING = {"", "Missing", "missing", "NA", "N/A", "na", "nan", "None", "-", "."}


def clean(value: object) -> str:
    return str(value if value is not None else "").strip()


def norm_chrom(chrom: object) -> str:
    value = clean(chrom)
    if value.lower().startswith("chr"):
        value = value[3:]
    return "MT" if value.upper() == "M" else value.upper() if value.upper() in {"X", "Y", "MT"} else value


def norm_allele(allele: object) -> str:
    return clean(allele).upper()


def is_present(value: object) -> bool:
    return clean(value) not in MISSING


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
    text = clean(raw).replace("|", ";").replace(",", ";")
    genes = []
    for piece in text.split(";"):
        gene = piece.strip().upper()
        if gene and gene not in MISSING:
            genes.append(gene)
    return genes


def load_gene_decisions() -> dict[str, dict[str, str]]:
    if not GENE_DECISIONS.exists():
        return {}
    with GENE_DECISIONS.open(newline="") as handle:
        return {row["gene"]: row for row in csv.DictReader(handle, delimiter="\t")}


def sortable_key(key: str) -> tuple[tuple[int, str], int, str]:
    chrom, pos, *_ = key.split("-")
    chrom_order = {"X": 23, "Y": 24, "MT": 25}
    if chrom.isdigit():
        chrom_sort = (int(chrom), chrom)
    else:
        chrom_sort = (chrom_order.get(chrom, 99), chrom)
    return chrom_sort, int(pos), key


def add_record(
    records: dict[str, dict[str, object]],
    members: list[dict[str, str]],
    *,
    key: str,
    source: str,
    source_row_id: str,
    source_variant_uid: str,
    source_variant_key: str,
    gene: str,
    label: str,
    raw_label: str,
    review_stars: str = "",
    review_status: str = "",
    confidence_tier: str = "",
    phenotype_available: str = "false",
    private_or_patient_linked: str = "false",
    gene_panel_decision: str = "",
) -> None:
    chrom, pos, ref, alt = key.split("-", 3)
    rec = records.setdefault(
        key,
        {
            "variant_id": key,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
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
        },
    )
    rec["genes"].update(split_genes(gene))  # type: ignore[union-attr]
    rec["sources"].add(source)  # type: ignore[union-attr]
    if label:
        rec["labels_3class"].add(label)  # type: ignore[union-attr]
    if raw_label:
        rec["raw_labels"].add(raw_label)  # type: ignore[union-attr]
    if review_stars:
        rec["review_stars"].add(review_stars)  # type: ignore[union-attr]
    if review_status:
        rec["review_statuses"].add(review_status)  # type: ignore[union-attr]
    if confidence_tier:
        rec["confidence_tiers"].add(confidence_tier)  # type: ignore[union-attr]
    if gene_panel_decision:
        rec["gene_panel_decisions"].add(gene_panel_decision)  # type: ignore[union-attr]
    if phenotype_available == "true":
        rec["phenotype_available"] = True
    if private_or_patient_linked == "true":
        rec["private_or_patient_linked"] = True
    rec["source_row_count"] = int(rec["source_row_count"]) + 1

    members.append(
        {
            "variant_id": key,
            "source": source,
            "source_row_id": source_row_id,
            "source_variant_uid": source_variant_uid,
            "source_variant_key": source_variant_key,
            "gene": gene,
            "label_3class": label,
            "raw_label": raw_label,
            "review_stars": review_stars,
            "review_status": review_status,
            "confidence_tier": confidence_tier,
            "phenotype_available": phenotype_available,
            "private_or_patient_linked": private_or_patient_linked,
            "gene_panel_decision": gene_panel_decision,
        }
    )


def add_clinvar(records: dict[str, dict[str, object]], members: list[dict[str, str]]) -> Counter:
    counts: Counter = Counter()
    with CLINVAR_PRIMARY.open(newline="") as handle:
        for idx, row in enumerate(csv.DictReader(handle, delimiter="\t"), start=1):
            counts["clinvar_rows_seen"] += 1
            key = variant_key(row.get("chrom"), row.get("pos"), row.get("ref"), row.get("alt"))
            if not key:
                counts["clinvar_rows_without_coordinates"] += 1
                continue
            counts["clinvar_rows_with_coordinates"] += 1
            counts[f"clinvar_label_{clean(row.get('label_3class')) or 'blank'}"] += 1
            add_record(
                records,
                members,
                key=key,
                source="clinvar",
                source_row_id=clean(row.get("variation_id")) or str(idx),
                source_variant_uid=clean(row.get("variant_id")),
                source_variant_key=clean(row.get("name")),
                gene=clean(row.get("gene")),
                label=clean(row.get("label_3class")),
                raw_label=clean(row.get("clinical_significance")),
                review_stars=clean(row.get("review_stars")),
                review_status=clean(row.get("review_status")),
                confidence_tier=clean(row.get("confidence_tier")),
                gene_panel_decision=clean(row.get("gene_panel_decision")),
            )
    return counts


def decision_for_source_gene(gene: str, decisions: dict[str, dict[str, str]]) -> str:
    source_genes = split_genes(gene)
    found = []
    for source_gene in source_genes:
        row = decisions.get(source_gene)
        if row:
            found.append(row.get("decision", ""))
    return "|".join(sorted(set(x for x in found if x)))


def add_core_sources(
    records: dict[str, dict[str, object]],
    members: list[dict[str, str]],
    decisions: dict[str, dict[str, str]],
) -> Counter:
    counts: Counter = Counter()
    for source, path in CORE_SOURCES.items():
        with path.open(newline="") as handle:
            for idx, row in enumerate(csv.DictReader(handle), start=1):
                counts[f"{source}_rows_seen"] += 1
                key = variant_key(row.get("chrom"), row.get("pos"), row.get("ref"), row.get("alt"))
                if not key:
                    counts[f"{source}_rows_without_coordinates"] += 1
                    continue
                counts[f"{source}_rows_with_coordinates"] += 1
                label = clean(row.get("target_3class"))
                counts[f"{source}_label_{label or 'blank'}"] += 1
                add_record(
                    records,
                    members,
                    key=key,
                    source=source,
                    source_row_id=str(idx),
                    source_variant_uid=clean(row.get("variant_uid")),
                    source_variant_key=clean(row.get("variant_key")),
                    gene=clean(row.get("gene")),
                    label=label,
                    raw_label=clean(row.get("target_classification_raw")),
                    phenotype_available="true" if source == "hiro" else "false",
                    private_or_patient_linked="true" if source == "hiro" else "false",
                    gene_panel_decision=decision_for_source_gene(clean(row.get("gene")), decisions),
                )
    return counts


def serialize(value: object) -> str:
    if isinstance(value, set):
        return "|".join(sorted(str(v) for v in value if str(v)))
    return clean(value)


def write_outputs(records: dict[str, dict[str, object]], members: list[dict[str, str]], counts: Counter) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sources = ["clinvar", "hiro", "emerge", "cardioboost"]
    registry_fields = [
        "variant_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "genes",
        "sources",
        "source_count",
        "source_row_count",
        "in_clinvar",
        "in_hiro",
        "in_emerge",
        "in_cardioboost",
        "phenotype_available",
        "private_or_patient_linked",
        "labels_3class",
        "raw_labels",
        "max_review_stars",
        "review_statuses",
        "confidence_tiers",
        "gene_panel_decisions",
    ]
    with REGISTRY_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=registry_fields, delimiter="\t")
        writer.writeheader()
        for key in sorted(records, key=sortable_key):
            rec = records[key]
            rec_sources = rec["sources"]
            review_stars = [int(x) for x in rec["review_stars"] if str(x).isdigit()]  # type: ignore[index]
            row = {field: "" for field in registry_fields}
            for field in ["variant_id", "chrom", "pos", "ref", "alt", "source_row_count"]:
                row[field] = serialize(rec.get(field, ""))
            row["genes"] = serialize(rec["genes"])
            row["sources"] = serialize(rec_sources)
            row["source_count"] = str(len(rec_sources))  # type: ignore[arg-type]
            for source in sources:
                row[f"in_{source}"] = "true" if source in rec_sources else "false"  # type: ignore[operator]
            row["phenotype_available"] = "true" if rec.get("phenotype_available") else "false"
            row["private_or_patient_linked"] = "true" if rec.get("private_or_patient_linked") else "false"
            row["labels_3class"] = serialize(rec["labels_3class"])
            row["raw_labels"] = serialize(rec["raw_labels"])
            row["max_review_stars"] = str(max(review_stars)) if review_stars else ""
            row["review_statuses"] = serialize(rec["review_statuses"])
            row["confidence_tiers"] = serialize(rec["confidence_tiers"])
            row["gene_panel_decisions"] = serialize(rec["gene_panel_decisions"])
            writer.writerow(row)

    source_fields = [
        "variant_id",
        "source",
        "source_row_id",
        "source_variant_uid",
        "source_variant_key",
        "gene",
        "label_3class",
        "raw_label",
        "review_stars",
        "review_status",
        "confidence_tier",
        "phenotype_available",
        "private_or_patient_linked",
        "gene_panel_decision",
    ]
    with SOURCES_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=source_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(members)

    summary = {
        "registry_rows": len(records),
        "source_membership_rows": len(members),
        "inputs": {
            "clinvar_primary": str(CLINVAR_PRIMARY.relative_to(ROOT)),
            **{source: str(path.relative_to(ROOT)) for source, path in CORE_SOURCES.items()},
        },
        "outputs": {
            "registry": str(REGISTRY_OUT.relative_to(ROOT)),
            "sources": str(SOURCES_OUT.relative_to(ROOT)),
            "summary": str(SUMMARY_OUT.relative_to(ROOT)),
        },
        "source_counts": dict(sorted(counts.items())),
        "source_combo_counts": dict(sorted(Counter(serialize(rec["sources"]) for rec in records.values()).items())),
        "label_combo_counts": dict(sorted(Counter(serialize(rec["labels_3class"]) for rec in records.values()).items())),
        "notes": [
            "ClinVar rows come from the post-gene-QC primary ClinVar slice, not raw variant_summary.",
            "HiRO/eMERGE/CardioBoost source rows are preserved when coordinate-complete, with gene-panel decisions flagged.",
            "Rows without complete coordinates remain counted in source_counts and require HGVS/protein fallback mapping.",
        ],
    }
    SUMMARY_OUT.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


def main() -> int:
    decisions = load_gene_decisions()
    records: dict[str, dict[str, object]] = {}
    members: list[dict[str, str]] = []
    counts = Counter()
    counts.update(add_clinvar(records, members))
    counts.update(add_core_sources(records, members, decisions))
    write_outputs(records, members, counts)
    print(json.dumps(json.loads(SUMMARY_OUT.read_text()), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
