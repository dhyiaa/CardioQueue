#!/usr/bin/env python3
"""Build the first final modeling table skeleton.

The table is one row per cleaned coordinate variant. It starts from the current
local feature registry so dbNSFP and existing protein/AlphaFold/UniProt columns
are already available, then adds standardized model-facing label, provenance,
inclusion, and feature-block status fields.
"""

from __future__ import annotations

import csv
import json
import argparse
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]

LOCAL_FEATURE_REGISTRY = ROOT / "datasets/variant_registry/interim/clean_combined_local_feature_registry.tsv"
SOURCE_MEMBERSHIP = ROOT / "datasets/variant_registry/interim/clean_combined_variant_registry_sources.tsv"
GENE_DECISIONS = ROOT / "datasets/gene_panels/interim/gene_panel_decisions.tsv"

OUT_DIR = ROOT / "datasets/modeling/interim"
TABLE_OUT = OUT_DIR / "final_modeling_table_skeleton.tsv"
SUMMARY_OUT = OUT_DIR / "final_modeling_table_skeleton.summary.json"

LABEL_NUMERIC = {"Benign": "0", "VUS": "1", "Pathogenic": "2"}
SOURCES = ("clinvar", "hiro", "emerge", "cardioboost")
FEATURE_PLACEHOLDERS = {
    "gnomad_browser_status": "not_joined",
    "gnomad_browser_missing_reason": "gnomAD browser Hail chunks are still running",
    "vep_status": "not_joined",
    "vep_missing_reason": "local VEP cache/download not complete",
    "spliceai_status": "not_joined",
    "spliceai_missing_reason": "SpliceAI source not merged yet",
    "clingen_gene_validity_status": "not_joined",
    "clingen_variant_evidence_status": "not_joined",
    "dssp_status": "not_joined",
    "dssp_missing_reason": "DSSP batch extraction not run yet",
    "freesasa_status": "not_joined",
    "freesasa_missing_reason": "FreeSASA batch extraction not run yet",
    "foldx_ddg_status": "not_joined",
    "foldx_ddg_missing_reason": "FoldX DDG pilot/batch not run yet",
    "alphamissense_direct_status": "not_joined",
    "alphamissense_direct_missing_reason": "direct AlphaMissense join not run yet; dbNSFP AlphaMissense may already be present",
}


def split_pipe(value: str | None) -> list[str]:
    if value is None:
        return []
    text = str(value).strip()
    if not text:
        return []
    return [piece for piece in text.split("|") if piece]


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def is_true(value: str | None) -> bool:
    return str(value or "").strip().lower() in {"true", "1", "yes", "y"}


def rel(path: Path) -> str:
    path = path.resolve()
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def load_gene_decisions(path: Path = GENE_DECISIONS) -> dict[str, dict[str, str]]:
    decisions: dict[str, dict[str, str]] = {}
    if not path.exists():
        return decisions
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            decisions[row["gene"]] = row
    return decisions


def load_source_membership(path: Path) -> dict[str, dict[str, object]]:
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
    with path.open(newline="") as handle:
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


def joined(values: set[str]) -> str:
    return "|".join(sorted(v for v in values if v))


def source_value(
    membership: dict[str, object],
    field: str,
    source: str,
) -> str:
    source_map = membership.get(field, {})
    return joined(source_map.get(source, set()))  # type: ignore[union-attr]


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


def primary_model_inclusion(gene_decisions: list[str], decision_rows: dict[str, dict[str, str]]) -> str:
    if not gene_decisions:
        return "unknown"
    flags = []
    for decision in gene_decisions:
        matches = [row for row in decision_rows.values() if row.get("decision") == decision]
        if matches:
            flags.extend(row.get("primary_model_inclusion", "") for row in matches)
    flags = [flag for flag in flags if flag]
    if "exclude" in flags:
        return "exclude"
    if flags and all(flag == "include" for flag in flags):
        return "include"
    if flags:
        return "|".join(sorted(set(flags)))
    if any(decision in {"remove_contamination", "alias_or_quarantine", "sensitivity_only", "separate_ttn"} for decision in gene_decisions):
        return "exclude"
    return "include"


def modeling_scope(gene_decisions: list[str]) -> str:
    decision_set = set(gene_decisions)
    if not decision_set:
        return "unknown"
    if "remove_contamination" in decision_set or "alias_or_quarantine" in decision_set:
        return "exclude_or_quarantine"
    if "separate_ttn" in decision_set:
        return "ttn_separate"
    if "sensitivity_only" in decision_set or "sensitivity_or_pool" in decision_set:
        return "sensitivity_only"
    if "conditional_phenocopy" in decision_set:
        return "expanded_scope"
    if "primary_keep_monitor" in decision_set:
        return "primary_monitor"
    if decision_set == {"primary_keep"}:
        return "primary"
    return "mixed_gene_scope"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local-feature-registry", type=Path, default=LOCAL_FEATURE_REGISTRY)
    parser.add_argument("--source-membership", type=Path, default=SOURCE_MEMBERSHIP)
    parser.add_argument("--gene-decisions", type=Path, default=GENE_DECISIONS)
    parser.add_argument("--out", type=Path, default=TABLE_OUT)
    parser.add_argument("--summary", type=Path, default=SUMMARY_OUT)
    args = parser.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    source_membership = load_source_membership(args.source_membership)
    gene_decision_rows = load_gene_decisions(args.gene_decisions)
    counts: Counter = Counter()

    with args.local_feature_registry.open(newline="") as in_handle, args.out.open("w", newline="") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        base_fieldnames = list(reader.fieldnames or [])
        added_fields = [
            "modeling_row_id",
            "genome_build",
            "coordinate_key",
            "primary_gene",
            "modeling_scope",
            "primary_model_inclusion",
            "model_label_3class",
            "model_label_numeric",
            "clean_supervised_label",
            "label_conflict_type",
            "label_count",
            "source_label_count",
            "default_split_role",
            "source_leakage_group",
        ]
        for source in SOURCES:
            added_fields.extend(
                [
                    f"{source}_labels",
                    f"{source}_raw_labels",
                    f"{source}_genes",
                    f"{source}_confidence_tiers",
                    f"{source}_review_stars",
                ]
            )
        added_fields.extend(FEATURE_PLACEHOLDERS)

        fieldnames = []
        for col in base_fieldnames[:5] + added_fields + base_fieldnames[5:]:
            if col not in fieldnames:
                fieldnames.append(col)

        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()

        for row in reader:
            counts["rows"] += 1
            variant_id = row["variant_id"]
            labels = split_pipe(row.get("labels_3class"))
            unique_labels = sorted(set(labels))
            gene_decisions = split_pipe(row.get("gene_panel_decisions"))
            genes = split_pipe(row.get("genes"))
            membership = source_membership.get(variant_id, {})
            source_label_count = sum(
                1
                for source in SOURCES
                if source_value(membership, "labels_by_source", source)
            )

            conflict = label_conflict_type(unique_labels)
            clean_label = len(unique_labels) == 1
            model_label = unique_labels[0] if clean_label else ""
            scope = modeling_scope(gene_decisions)
            inclusion = primary_model_inclusion(gene_decisions, gene_decision_rows)

            if not clean_label:
                default_split_role = "label_conflict_review"
            elif inclusion == "exclude":
                default_split_role = "excluded_gene_scope"
            elif is_true(row.get("in_hiro")):
                default_split_role = "hiro_high_trust_stratum"
            elif is_true(row.get("in_emerge")):
                default_split_role = "emerge_candidate_validation_or_source_stratum"
            elif is_true(row.get("in_cardioboost")):
                default_split_role = "cardioboost_binary_source_stratum"
            else:
                default_split_role = "training_candidate_unassigned"

            row.update(
                {
                    "modeling_row_id": variant_id,
                    "genome_build": "GRCh38",
                    "coordinate_key": variant_id,
                    "primary_gene": genes[0] if genes else "",
                    "modeling_scope": scope,
                    "primary_model_inclusion": inclusion,
                    "model_label_3class": model_label,
                    "model_label_numeric": LABEL_NUMERIC.get(model_label, ""),
                    "clean_supervised_label": bool_text(clean_label),
                    "label_conflict_type": conflict,
                    "label_count": str(len(unique_labels)),
                    "source_label_count": str(source_label_count),
                    "default_split_role": default_split_role,
                    "source_leakage_group": variant_id,
                }
            )

            for source in SOURCES:
                row[f"{source}_labels"] = source_value(membership, "labels_by_source", source)
                row[f"{source}_raw_labels"] = source_value(membership, "raw_labels_by_source", source)
                row[f"{source}_genes"] = source_value(membership, "genes_by_source", source)
                row[f"{source}_confidence_tiers"] = source_value(membership, "confidence_by_source", source)
                row[f"{source}_review_stars"] = source_value(membership, "review_stars_by_source", source)

            row.update(FEATURE_PLACEHOLDERS)
            writer.writerow(row)

            counts[f"model_label_{model_label or 'blank'}"] += 1
            counts[f"label_conflict_{conflict}"] += 1
            counts[f"modeling_scope_{scope}"] += 1
            counts[f"primary_model_inclusion_{inclusion}"] += 1
            counts[f"default_split_role_{default_split_role}"] += 1
            counts[f"dbnsfp_{row.get('dbnsfp_status') or 'blank'}"] += 1
            counts[f"protein_{row.get('protein_feature_status') or 'blank'}"] += 1
            counts[f"alphafold_{row.get('alphafold_plddt_status') or 'blank'}"] += 1

    summary = {
        "inputs": {
            "local_feature_registry": rel(args.local_feature_registry),
            "source_membership": rel(args.source_membership),
            "gene_decisions": rel(args.gene_decisions),
        },
        "outputs": {
            "table": rel(args.out),
            "summary": rel(args.summary),
        },
        "counts": dict(sorted(counts.items())),
        "notes": [
            "Variant-level skeleton: one row per cleaned chrom-pos-ref-alt key.",
            "Mixed-label variants are retained but are not clean supervised labels.",
            "default_split_role is a placeholder and should be replaced by an explicit leakage-aware split file before modeling.",
            "Feature status placeholders mark downloads/batches that are not yet merged.",
        ],
    }
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
