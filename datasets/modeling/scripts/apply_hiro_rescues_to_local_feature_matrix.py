#!/usr/bin/env python3
"""Apply promoted HiRO rescue registry/source flags to the local feature matrix."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]

MATRIX_IN = ROOT / "datasets/modeling/interim/final_modeling_table_local_features_frozen.tsv"
REGISTRY = ROOT / "datasets/variant_registry/interim/clean_combined_variant_registry_with_hiro_rescues.tsv"
SOURCES = ROOT / "datasets/variant_registry/interim/clean_combined_variant_registry_sources_with_hiro_rescues.tsv"
GENE_DECISIONS = ROOT / "datasets/gene_panels/interim/gene_panel_decisions.tsv"

MATRIX_OUT = ROOT / "datasets/modeling/interim/final_modeling_table_local_features_frozen_with_hiro_rescues.tsv"
SUMMARY_OUT = ROOT / "datasets/modeling/interim/final_modeling_table_local_features_frozen_with_hiro_rescues.summary.json"

SOURCE_NAMES = ("clinvar", "hiro", "emerge", "cardioboost")
LABEL_NUMERIC = {"Benign": "0", "VUS": "1", "Pathogenic": "2"}


def split_pipe(value: str | None) -> list[str]:
    return [piece for piece in str(value or "").split("|") if piece]


def joined(values: set[str]) -> str:
    return "|".join(sorted(v for v in values if v))


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def is_true(value: str | None) -> bool:
    return str(value or "").strip().lower() in {"true", "1", "yes", "y"}


def load_gene_decisions() -> dict[str, dict[str, str]]:
    with GENE_DECISIONS.open(newline="") as handle:
        return {row["gene"]: row for row in csv.DictReader(handle, delimiter="\t")}


def primary_model_inclusion(gene_decisions: list[str], decision_rows: dict[str, dict[str, str]]) -> str:
    if not gene_decisions:
        return "unknown"
    flags = []
    for decision in gene_decisions:
        flags.extend(
            row.get("primary_model_inclusion", "")
            for row in decision_rows.values()
            if row.get("decision") == decision
        )
    flags = [flag for flag in flags if flag]
    if "exclude" in flags:
        return "exclude"
    if flags and all(flag == "include" for flag in flags):
        return "include"
    if flags:
        return joined(set(flags))
    if any(decision in {"remove_contamination", "alias_or_quarantine", "sensitivity_only", "separate_ttn"} for decision in gene_decisions):
        return "exclude"
    return "include"


def modeling_scope(gene_decisions: list[str]) -> str:
    decisions = set(gene_decisions)
    if not decisions:
        return "unknown"
    if "remove_contamination" in decisions or "alias_or_quarantine" in decisions:
        return "exclude_or_quarantine"
    if "separate_ttn" in decisions:
        return "ttn_separate"
    if "sensitivity_only" in decisions or "sensitivity_or_pool" in decisions:
        return "sensitivity_only"
    if "conditional_phenocopy" in decisions:
        return "expanded_scope"
    if "primary_keep_monitor" in decisions:
        return "primary_monitor"
    if decisions == {"primary_keep"}:
        return "primary"
    return "mixed_gene_scope"


def label_conflict_type(labels: list[str]) -> str:
    label_set = set(labels)
    if len(label_set) <= 1:
        return "none"
    if label_set == {"Benign", "VUS"}:
        return "benign_vus_conflict"
    if label_set == {"Pathogenic", "VUS"}:
        return "pathogenic_vus_conflict"
    if label_set == {"Benign", "Pathogenic"}:
        return "benign_pathogenic_conflict"
    return "three_way_conflict"


def load_registry() -> dict[str, dict[str, str]]:
    with REGISTRY.open(newline="") as handle:
        return {row["variant_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def load_source_membership() -> dict[str, dict[str, object]]:
    by_variant: dict[str, dict[str, object]] = defaultdict(
        lambda: {
            "source_rows": 0,
            "labels_by_source": defaultdict(set),
            "raw_labels_by_source": defaultdict(set),
            "genes_by_source": defaultdict(set),
            "confidence_by_source": defaultdict(set),
            "review_stars_by_source": defaultdict(set),
        }
    )
    with SOURCES.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            variant_id = row["variant_id"]
            source = row["source"]
            rec = by_variant[variant_id]
            rec["source_rows"] = int(rec["source_rows"]) + 1
            if row.get("label_3class"):
                rec["labels_by_source"][source].add(row["label_3class"])  # type: ignore[index]
            if row.get("raw_label"):
                rec["raw_labels_by_source"][source].add(row["raw_label"])  # type: ignore[index]
            if row.get("gene"):
                rec["genes_by_source"][source].add(row["gene"])  # type: ignore[index]
            if row.get("confidence_tier"):
                rec["confidence_by_source"][source].add(row["confidence_tier"])  # type: ignore[index]
            if row.get("review_stars"):
                rec["review_stars_by_source"][source].add(row["review_stars"])  # type: ignore[index]
    return by_variant


def source_value(membership: dict[str, object], field: str, source: str) -> str:
    source_map = membership.get(field, {})
    return joined(source_map.get(source, set()))  # type: ignore[union-attr]


def apply_registry_fields(
    row: dict[str, str],
    registry_row: dict[str, str],
    membership: dict[str, object],
    decision_rows: dict[str, dict[str, str]],
) -> dict[str, str]:
    for field in [
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
    ]:
        if field in row:
            row[field] = registry_row.get(field, "")

    labels = sorted(set(split_pipe(registry_row.get("labels_3class"))))
    gene_decisions = split_pipe(registry_row.get("gene_panel_decisions"))
    genes = split_pipe(registry_row.get("genes"))
    conflict = label_conflict_type(labels)
    clean_label = len(labels) == 1
    model_label = labels[0] if clean_label else ""
    inclusion = primary_model_inclusion(gene_decisions, decision_rows)

    row["primary_gene"] = genes[0] if genes else ""
    row["modeling_scope"] = modeling_scope(gene_decisions)
    row["primary_model_inclusion"] = inclusion
    row["model_label_3class"] = model_label
    row["model_label_numeric"] = LABEL_NUMERIC.get(model_label, "")
    row["clean_supervised_label"] = bool_text(clean_label)
    row["label_conflict_type"] = conflict
    row["label_count"] = str(len(labels))
    row["source_label_count"] = str(
        sum(1 for source in SOURCE_NAMES if source_value(membership, "labels_by_source", source))
    )
    if not clean_label:
        row["default_split_role"] = "label_conflict_review"
    elif inclusion == "exclude":
        row["default_split_role"] = "excluded_gene_scope"
    elif is_true(registry_row.get("in_hiro")):
        row["default_split_role"] = "hiro_high_trust_stratum"
    elif is_true(registry_row.get("in_emerge")):
        row["default_split_role"] = "emerge_candidate_validation_or_source_stratum"
    elif is_true(registry_row.get("in_cardioboost")):
        row["default_split_role"] = "cardioboost_binary_source_stratum"
    else:
        row["default_split_role"] = "training_candidate_unassigned"
    row["source_leakage_group"] = row["variant_id"]

    for source in SOURCE_NAMES:
        row[f"{source}_labels"] = source_value(membership, "labels_by_source", source)
        row[f"{source}_raw_labels"] = source_value(membership, "raw_labels_by_source", source)
        row[f"{source}_genes"] = source_value(membership, "genes_by_source", source)
        row[f"{source}_confidence_tiers"] = source_value(membership, "confidence_by_source", source)
        row[f"{source}_review_stars"] = source_value(membership, "review_stars_by_source", source)

    return row


def blank_new_row(fieldnames: list[str], registry_row: dict[str, str]) -> dict[str, str]:
    row = {field: "" for field in fieldnames}
    for field in ["variant_id", "chrom", "pos", "ref", "alt"]:
        row[field] = registry_row.get(field, "")
    row["modeling_row_id"] = row["variant_id"]
    row["genome_build"] = "GRCh38"
    row["coordinate_key"] = row["variant_id"]
    placeholder_statuses = {
        "dbnsfp_status": "not_joined",
        "dbnsfp_missing_reason": "added after dbNSFP selected-feature table was built; rerun dbNSFP join for rescued rows",
        "gnomad_browser_status": "not_joined",
        "gnomad_browser_missing_reason": "gnomAD browser Hail chunks are still running",
        "vep_status": "not_joined",
        "vep_missing_reason": "local VEP cache/download not complete",
        "spliceai_status": "not_joined",
        "spliceai_missing_reason": "SpliceAI source not merged yet",
        "clingen_gene_validity_status": "not_joined",
        "clingen_variant_evidence_status": "not_joined",
        "protein_feature_status": "not_available",
        "protein_feature_missing_reason": "added from HiRO coordinate rescue; rerun protein feature parser if protein HGVS is eligible",
        "alphafold_plddt_status": "not_available",
        "dssp_status": "not_joined",
        "freesasa_status": "not_joined",
        "foldx_ddg_status": "not_joined",
        "alphamissense_direct_status": "not_joined",
        "alphamissense_direct_missing_reason": "added after direct AlphaMissense selected-feature table was built; rerun direct AlphaMissense join",
    }
    for field, value in placeholder_statuses.items():
        if field in row:
            row[field] = value
    return row


def main() -> int:
    registry = load_registry()
    memberships = load_source_membership()
    decision_rows = load_gene_decisions()
    counts: Counter[str] = Counter()

    with MATRIX_IN.open(newline="") as in_handle, MATRIX_OUT.open("w", newline="") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()

        seen = set()
        for row in reader:
            variant_id = row["variant_id"]
            seen.add(variant_id)
            registry_row = registry.get(variant_id)
            if registry_row:
                before_in_hiro = row.get("in_hiro")
                row = apply_registry_fields(row, registry_row, memberships.get(variant_id, {}), decision_rows)
                if before_in_hiro != row.get("in_hiro"):
                    counts["existing_matrix_rows_hiro_status_changed"] += 1
            else:
                counts["matrix_rows_missing_promoted_registry"] += 1
            writer.writerow(row)
            counts["matrix_input_rows"] += 1

        for variant_id in sorted(set(registry) - seen):
            registry_row = registry[variant_id]
            row = blank_new_row(fieldnames, registry_row)
            row = apply_registry_fields(row, registry_row, memberships.get(variant_id, {}), decision_rows)
            writer.writerow(row)
            counts["rescued_registry_rows_appended_to_matrix"] += 1

    summary = {
        "inputs": {
            "matrix": str(MATRIX_IN.relative_to(ROOT)),
            "registry_with_hiro_rescues": str(REGISTRY.relative_to(ROOT)),
            "sources_with_hiro_rescues": str(SOURCES.relative_to(ROOT)),
        },
        "output": str(MATRIX_OUT.relative_to(ROOT)),
        "output_rows": counts["matrix_input_rows"] + counts["rescued_registry_rows_appended_to_matrix"],
        "counts": dict(sorted(counts.items())),
        "notes": [
            "This table applies promoted HiRO registry/source flags to the frozen local-feature matrix.",
            "New rescued coordinate rows are appended with feature-block placeholders; rerun feature joins later for full annotation.",
            "The pre-rescue frozen local-feature matrix is unchanged.",
        ],
    }
    SUMMARY_OUT.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
