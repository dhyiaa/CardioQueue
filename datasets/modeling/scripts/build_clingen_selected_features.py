#!/usr/bin/env python3
"""Create selected model-ready ClinGen features and merge them into the skeleton."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]

SKELETON = ROOT / "datasets/modeling/interim/final_modeling_table_skeleton.tsv"
CLINGEN = ROOT / "datasets/feature_sources/clingen/interim/clingen_modeling_features.tsv"
OUT_DIR = ROOT / "datasets/modeling/interim"
SELECTED_OUT = OUT_DIR / "clingen_selected_features.tsv"
MERGED_OUT = OUT_DIR / "final_modeling_table_with_clingen.tsv"
SUMMARY_OUT = OUT_DIR / "clingen_selected_features.summary.json"

MISSING_VALUES = {"", ".", "NA", "N/A", "na", "nan", "None", "Missing"}

KEY_COLUMNS = ["variant_id", "primary_gene"]

SELECTED_COLUMNS = [
    "clingen_gene_validity_status",
    "clingen_gene_validity_max_classification",
    "clingen_gene_validity_max_rank",
    "clingen_gene_validity_has_definitive_strong",
    "clingen_gene_validity_has_moderate_plus",
    "clingen_gene_validity_rows",
    "clingen_dosage_status",
    "clingen_haploinsufficiency_score",
    "clingen_haploinsufficiency_description",
    "clingen_triplosensitivity_score",
    "clingen_triplosensitivity_description",
    "clingen_variant_evidence_status",
    "clingen_variant_evidence_match_count",
    "clingen_variant_evidence_classification",
    "clingen_variant_evidence_expert_panel",
    "clingen_variant_evidence_condition",
    "clingen_variant_evidence_met_codes",
]

TRACE_COLUMNS = [
    "clingen_gene_validity_classifications",
    "clingen_gene_validity_diseases",
    "clingen_gene_validity_expert_panels",
    "clingen_dosage_date",
    "clingen_variant_evidence_mondo",
    "clingen_variant_evidence_ca_id",
    "clingen_variant_evidence_clinvar_variation_id",
    "clingen_variant_evidence_pcer_doc_id",
    "clingen_variant_evidence_approved_date",
    "clingen_variant_evidence_published_date",
    "clingen_variant_evidence_hgvs_match",
]

HELPER_COLUMNS = [
    "clingen_gene_validity_is_missing",
    "clingen_dosage_is_missing",
    "clingen_variant_evidence_is_missing",
    "clingen_variant_evidence_has_pathogenic_assertion",
    "clingen_variant_evidence_has_benign_assertion",
    "clingen_haploinsufficiency_is_sufficient",
    "clingen_triplosensitivity_is_sufficient",
]


def is_missing(value: object) -> bool:
    return str(value if value is not None else "").strip() in MISSING_VALUES


def clean(value: object) -> str:
    text = str(value if value is not None else "").strip()
    return "" if text in MISSING_VALUES else text


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def sufficient_dosage(value: str) -> bool:
    text = clean(value).lower()
    return text.startswith("3") or "sufficient evidence" in text


def selected_row(row: dict[str, str]) -> dict[str, str]:
    out = {field: clean(row.get(field, "")) for field in KEY_COLUMNS + SELECTED_COLUMNS + TRACE_COLUMNS}
    gene_missing = out["clingen_gene_validity_status"] != "ok"
    dosage_missing = out["clingen_dosage_status"] != "ok"
    variant_missing = out["clingen_variant_evidence_status"] != "ok"
    variant_class = out["clingen_variant_evidence_classification"].lower()

    out["clingen_gene_validity_is_missing"] = bool_text(gene_missing)
    out["clingen_dosage_is_missing"] = bool_text(dosage_missing)
    out["clingen_variant_evidence_is_missing"] = bool_text(variant_missing)
    out["clingen_variant_evidence_has_pathogenic_assertion"] = bool_text(
        variant_class in {"pathogenic", "likely pathogenic"}
    )
    out["clingen_variant_evidence_has_benign_assertion"] = bool_text(
        variant_class in {"benign", "likely benign"}
    )
    out["clingen_haploinsufficiency_is_sufficient"] = bool_text(
        sufficient_dosage(out["clingen_haploinsufficiency_score"])
        or sufficient_dosage(out["clingen_haploinsufficiency_description"])
    )
    out["clingen_triplosensitivity_is_sufficient"] = bool_text(
        sufficient_dosage(out["clingen_triplosensitivity_score"])
        or sufficient_dosage(out["clingen_triplosensitivity_description"])
    )
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    counts: Counter[str] = Counter()

    clingen_rows = read_tsv(CLINGEN)
    selected_rows = [selected_row(row) for row in clingen_rows]
    by_variant = {row["variant_id"]: row for row in selected_rows}

    selected_fields = KEY_COLUMNS + SELECTED_COLUMNS + HELPER_COLUMNS + TRACE_COLUMNS
    write_tsv(SELECTED_OUT, selected_rows, selected_fields)

    skeleton_rows = read_tsv(SKELETON)
    skeleton_fields = list(skeleton_rows[0].keys()) if skeleton_rows else []
    add_fields = [field for field in selected_fields if field not in skeleton_fields]
    merged_fields = skeleton_fields + add_fields

    merged_rows: list[dict[str, str]] = []
    for row in skeleton_rows:
        variant_id = row["variant_id"]
        payload = by_variant.get(variant_id)
        if not payload:
            counts["merge_missing_clingen_payload"] += 1
            payload = {field: "" for field in selected_fields}
            payload["variant_id"] = variant_id
            payload["primary_gene"] = row.get("primary_gene", "")
            payload["clingen_gene_validity_status"] = "not_joined"
            payload["clingen_dosage_status"] = "not_joined"
            payload["clingen_variant_evidence_status"] = "not_joined"

        # Replace skeleton placeholder columns with real selected ClinGen values.
        for field in selected_fields:
            if field in {"variant_id", "primary_gene"}:
                continue
            row[field] = payload.get(field, "")

        counts[f"gene_validity_{row.get('clingen_gene_validity_status') or 'blank'}"] += 1
        counts[f"dosage_{row.get('clingen_dosage_status') or 'blank'}"] += 1
        counts[f"variant_evidence_{row.get('clingen_variant_evidence_status') or 'blank'}"] += 1
        if row.get("clingen_gene_validity_has_definitive_strong") == "true":
            counts["gene_validity_definitive_strong_true"] += 1
        if row.get("clingen_variant_evidence_has_pathogenic_assertion") == "true":
            counts["variant_evidence_pathogenic_assertion_true"] += 1
        if row.get("clingen_variant_evidence_has_benign_assertion") == "true":
            counts["variant_evidence_benign_assertion_true"] += 1
        merged_rows.append(row)

    write_tsv(MERGED_OUT, merged_rows, merged_fields)

    summary = {
        "inputs": {
            "modeling_skeleton": str(SKELETON.relative_to(ROOT)),
            "clingen_modeling_features": str(CLINGEN.relative_to(ROOT)),
        },
        "outputs": {
            "selected_features": str(SELECTED_OUT.relative_to(ROOT)),
            "merged_modeling_table": str(MERGED_OUT.relative_to(ROOT)),
            "summary": str(SUMMARY_OUT.relative_to(ROOT)),
        },
        "rows": {
            "clingen_input": len(clingen_rows),
            "selected_features": len(selected_rows),
            "skeleton_input": len(skeleton_rows),
            "merged_output": len(merged_rows),
        },
        "selected_columns": SELECTED_COLUMNS,
        "helper_columns": HELPER_COLUMNS,
        "trace_columns": TRACE_COLUMNS,
        "counts": dict(sorted(counts.items())),
        "notes": [
            "Non-destructive merge: final_modeling_table_skeleton.tsv is unchanged.",
            "Merged table replaces ClinGen placeholder status columns with real ClinGen feature values.",
            "Trace columns are retained for QC but can be excluded from model training if too high-cardinality.",
        ],
    }
    SUMMARY_OUT.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
