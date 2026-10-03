#!/usr/bin/env python3
"""Select and merge protein/structure features into the modeling table."""

from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]

MODELING_INPUT = ROOT / "datasets/modeling/interim/final_modeling_table_with_clingen.tsv"
PROTEIN_FEATURES = ROOT / "datasets/feature_sources/protein_structure/interim/protein_features.tsv"
DSSP_FREESASA = ROOT / "datasets/feature_sources/protein_structure/interim/dssp_freesasa_features.tsv"
OUT_DIR = ROOT / "datasets/modeling/interim"
SELECTED_OUT = OUT_DIR / "protein_structure_selected_features.tsv"
MERGED_OUT = OUT_DIR / "final_modeling_table_with_clingen_protein.tsv"
SUMMARY_OUT = OUT_DIR / "protein_structure_selected_features.summary.json"

MISSING_VALUES = {"", ".", "NA", "N/A", "na", "nan", "None", "Missing"}
SOURCE_PRIORITY = {"hiro": 0, "cardioboost": 1, "emerge": 2}

KEY_COLUMNS = ["variant_id"]

PROTEIN_COLUMNS = [
    "protein_feature_status",
    "protein_feature_missing_reason",
    "uniprot_accession",
    "uniprot_length",
    "protein_position",
    "protein_ref_aa",
    "protein_alt_aa",
    "protein_variant_type",
    "protein_position_normalized",
    "uniprot_domain_hit",
    "uniprot_domain_names",
    "uniprot_region_hit",
    "uniprot_region_names",
    "uniprot_motif_hit",
    "uniprot_motif_names",
    "uniprot_binding_site_hit",
    "uniprot_binding_site_names",
    "uniprot_active_site_hit",
    "uniprot_active_site_names",
    "uniprot_any_functional_feature_hit",
    "alphafold_plddt_status",
    "alphafold_plddt_missing_reason",
    "alphafold_residue_plddt",
    "alphafold_plddt_bin",
    "alphafold_fragment_count",
]

DSSP_FREESASA_COLUMNS = [
    "dssp_status",
    "dssp_missing_reason",
    "dssp_secondary_structure",
    "dssp_secondary_structure_class",
    "dssp_acc",
    "dssp_phi",
    "dssp_psi",
    "freesasa_status",
    "freesasa_missing_reason",
    "freesasa_total",
    "freesasa_polar",
    "freesasa_apolar",
    "freesasa_relative",
    "freesasa_exposure_bin",
]

QC_COLUMNS = [
    "protein_selected_source_dataset",
    "protein_selected_variant_uid",
    "protein_source_row_count",
    "protein_duplicate_variant_flag",
]

HELPER_COLUMNS = [
    "protein_feature_is_missing",
    "alphafold_plddt_is_missing",
    "dssp_is_missing",
    "freesasa_is_missing",
    "protein_structure_features_available",
    "protein_functional_feature_hit_flag",
    "protein_domain_hit_flag",
    "protein_binding_or_active_site_hit_flag",
]


def clean(value: object) -> str:
    text = str(value if value is not None else "").strip()
    return "" if text in MISSING_VALUES else text


def is_missing(value: object) -> bool:
    return clean(value) == ""


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def normalize_pos(value: object) -> str:
    text = clean(value)
    if not text:
        return ""
    return re.sub(r"\.0$", "", text)


def variant_id_from_row(row: dict[str, str]) -> str:
    chrom = clean(row.get("chrom", ""))
    pos = normalize_pos(row.get("pos", ""))
    ref = clean(row.get("ref", "")).upper()
    alt = clean(row.get("alt", "")).upper()
    if not chrom or not pos or not ref or not alt:
        return ""
    return f"{chrom}-{pos}-{ref}-{alt}"


def dssp_key(row: dict[str, str]) -> tuple[str, str]:
    return (row.get("source_dataset", ""), row.get("variant_uid", ""))


def row_score(row: dict[str, str]) -> tuple[int, int, int, int, int, str]:
    protein_ok = 0 if row.get("protein_feature_status") == "ok" else 1
    alphafold_ok = 0 if row.get("alphafold_plddt_status") == "ok" else 1
    dssp_ok = 0 if row.get("dssp_status") == "ok" else 1
    freesasa_ok = 0 if row.get("freesasa_status") == "ok" else 1
    source = SOURCE_PRIORITY.get(row.get("source_dataset", ""), 9)
    return (protein_ok, alphafold_ok, dssp_ok, freesasa_ok, source, row.get("variant_uid", ""))


def add_helpers(row: dict[str, str]) -> None:
    row["protein_feature_is_missing"] = bool_text(row.get("protein_feature_status") != "ok")
    row["alphafold_plddt_is_missing"] = bool_text(row.get("alphafold_plddt_status") != "ok")
    row["dssp_is_missing"] = bool_text(row.get("dssp_status") != "ok")
    row["freesasa_is_missing"] = bool_text(row.get("freesasa_status") != "ok")
    row["protein_structure_features_available"] = bool_text(
        row.get("alphafold_plddt_status") == "ok"
        and row.get("dssp_status") == "ok"
        and row.get("freesasa_status") == "ok"
    )
    row["protein_functional_feature_hit_flag"] = bool_text(
        str(row.get("uniprot_any_functional_feature_hit", "")).lower() == "true"
    )
    row["protein_domain_hit_flag"] = bool_text(
        str(row.get("uniprot_domain_hit", "")).lower() == "true"
    )
    row["protein_binding_or_active_site_hit_flag"] = bool_text(
        str(row.get("uniprot_binding_site_hit", "")).lower() == "true"
        or str(row.get("uniprot_active_site_hit", "")).lower() == "true"
    )


def build_selected() -> tuple[list[dict[str, str]], dict[str, int]]:
    protein_rows = read_tsv(PROTEIN_FEATURES)
    dssp_rows = {dssp_key(row): row for row in read_tsv(DSSP_FREESASA)}

    by_variant: dict[str, list[dict[str, str]]] = defaultdict(list)
    input_counts = Counter()
    for row in protein_rows:
        input_counts["protein_input_rows"] += 1
        variant_id = variant_id_from_row(row)
        if not variant_id:
            input_counts["protein_rows_without_complete_coordinates"] += 1
            continue
        out = {"variant_id": variant_id}
        for field in PROTEIN_COLUMNS:
            out[field] = clean(row.get(field, ""))
        out["protein_selected_source_dataset"] = row.get("source_dataset", "")
        out["protein_selected_variant_uid"] = row.get("variant_uid", "")

        dssp = dssp_rows.get(dssp_key(row), {})
        for field in DSSP_FREESASA_COLUMNS:
            out[field] = clean(dssp.get(field, ""))
        if not dssp and out.get("alphafold_plddt_status") == "ok":
            out["dssp_status"] = "not_joined"
            out["dssp_missing_reason"] = "No DSSP row for selected source variant."
            out["freesasa_status"] = "not_joined"
            out["freesasa_missing_reason"] = "No FreeSASA row for selected source variant."
        elif not dssp:
            reason = out.get("alphafold_plddt_missing_reason") or "AlphaFold residue was not queryable/resolved."
            out["dssp_status"] = "not_queryable"
            out["dssp_missing_reason"] = reason
            out["freesasa_status"] = "not_queryable"
            out["freesasa_missing_reason"] = reason
        by_variant[variant_id].append(out)

    selected: list[dict[str, str]] = []
    for variant_id, rows in by_variant.items():
        rows.sort(key=row_score)
        best = rows[0]
        best["protein_source_row_count"] = str(len(rows))
        best["protein_duplicate_variant_flag"] = bool_text(len(rows) > 1)
        add_helpers(best)
        selected.append(best)

    selected.sort(key=lambda row: row["variant_id"])
    input_counts["protein_coordinate_variant_rows"] = sum(len(rows) for rows in by_variant.values())
    input_counts["protein_unique_variant_ids"] = len(selected)
    input_counts["protein_duplicate_variant_ids"] = sum(1 for rows in by_variant.values() if len(rows) > 1)
    return selected, dict(input_counts)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    selected_rows, input_counts = build_selected()
    selected_by_variant = {row["variant_id"]: row for row in selected_rows}

    selected_fields = KEY_COLUMNS + PROTEIN_COLUMNS + DSSP_FREESASA_COLUMNS + HELPER_COLUMNS + QC_COLUMNS
    write_tsv(SELECTED_OUT, selected_rows, selected_fields)

    modeling_rows = read_tsv(MODELING_INPUT)
    modeling_fields = list(modeling_rows[0].keys()) if modeling_rows else []
    add_fields = [field for field in selected_fields if field not in modeling_fields]
    merged_fields = modeling_fields + add_fields

    counts: Counter[str] = Counter(input_counts)
    merged_rows: list[dict[str, str]] = []
    for row in modeling_rows:
        payload = selected_by_variant.get(row["variant_id"])
        if payload:
            counts["merged_rows_with_protein_payload"] += 1
            for field in selected_fields:
                if field == "variant_id":
                    continue
                row[field] = payload.get(field, "")
        else:
            counts["merged_rows_without_protein_payload"] += 1
            # Replace the existing local-registry status with an explicit final merge status.
            row["protein_feature_status"] = "not_available"
            row["protein_feature_missing_reason"] = "No protein feature row for this variant in current HiRO/eMERGE/CardioBoost protein table."
            row["alphafold_plddt_status"] = "not_available"
            row["alphafold_plddt_missing_reason"] = "No protein feature row for this variant."
            row["dssp_status"] = "not_available"
            row["dssp_missing_reason"] = "No protein feature row for this variant."
            row["freesasa_status"] = "not_available"
            row["freesasa_missing_reason"] = "No protein feature row for this variant."
            for field in HELPER_COLUMNS:
                row[field] = "false"
            row["protein_feature_is_missing"] = "true"
            row["alphafold_plddt_is_missing"] = "true"
            row["dssp_is_missing"] = "true"
            row["freesasa_is_missing"] = "true"
            row["protein_source_row_count"] = "0"
            row["protein_duplicate_variant_flag"] = "false"

        counts[f"protein_feature_{row.get('protein_feature_status') or 'blank'}"] += 1
        counts[f"alphafold_plddt_{row.get('alphafold_plddt_status') or 'blank'}"] += 1
        counts[f"dssp_{row.get('dssp_status') or 'blank'}"] += 1
        counts[f"freesasa_{row.get('freesasa_status') or 'blank'}"] += 1
        if row.get("protein_structure_features_available") == "true":
            counts["protein_structure_features_available_true"] += 1
        if row.get("protein_functional_feature_hit_flag") == "true":
            counts["protein_functional_feature_hit_true"] += 1
        merged_rows.append(row)

    write_tsv(MERGED_OUT, merged_rows, merged_fields)

    summary = {
        "inputs": {
            "modeling_input": str(MODELING_INPUT.relative_to(ROOT)),
            "protein_features": str(PROTEIN_FEATURES.relative_to(ROOT)),
            "dssp_freesasa_features": str(DSSP_FREESASA.relative_to(ROOT)),
        },
        "outputs": {
            "selected_features": str(SELECTED_OUT.relative_to(ROOT)),
            "merged_modeling_table": str(MERGED_OUT.relative_to(ROOT)),
            "summary": str(SUMMARY_OUT.relative_to(ROOT)),
        },
        "rows": {
            "selected_features": len(selected_rows),
            "modeling_input": len(modeling_rows),
            "merged_output": len(merged_rows),
        },
        "protein_columns": PROTEIN_COLUMNS,
        "dssp_freesasa_columns": DSSP_FREESASA_COLUMNS,
        "helper_columns": HELPER_COLUMNS,
        "qc_columns": QC_COLUMNS,
        "counts": dict(sorted(counts.items())),
        "notes": [
            "Non-destructive merge: final_modeling_table_with_clingen.tsv is unchanged.",
            "Protein rows currently come from HiRO, eMERGE, and CardioBoost source tables.",
            "When duplicate source rows collapse to one variant_id, rows with protein/AlphaFold/DSSP/FreeSASA coverage are preferred.",
        ],
    }
    SUMMARY_OUT.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
