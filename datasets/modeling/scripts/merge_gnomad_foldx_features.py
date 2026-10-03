#!/usr/bin/env python3
"""Merge parsed gnomAD browser features and FoldX DDG into the current matrix."""

from __future__ import annotations

import csv
import json
import argparse
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]

MATRIX = ROOT / "datasets/modeling/interim/final_modeling_table_local_features_frozen_with_hiro_rescues.tsv"
GNOMAD = ROOT / "datasets/feature_sources/gnomad_v4/interim/gnomad_browser_v4_1_parsed_features.tsv"
FOLDX = ROOT / "datasets/feature_sources/protein_structure/interim/foldxpro_ddg_batch.tsv"

OUT = ROOT / "datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx.tsv"
SUMMARY = ROOT / "datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx.summary.json"
FOLDX_SELECTED = ROOT / "datasets/modeling/interim/foldx_ddg_variant_selected_features.tsv"

KEY_FIELDS = {"variant_id", "chrom", "pos", "ref", "alt", "genes", "sources", "labels_3class"}
MISSING_VALUES = {"", ".", "NA", "N/A", "Missing", "missing", "nan", "None", "-"}

GNOMAD_FINAL_FIELDS = [
    "gnomad_final_status",
    "gnomad_final_source",
    "gnomad_final_missing_reason",
    "gnomad_final_af",
    "gnomad_final_popmax_af",
    "gnomad_final_popmax_pop",
    "gnomad_final_homozygote_count",
    "gnomad_final_hemizygote_count",
    "gnomad_final_af_is_missing",
    "gnomad_confirmed_absent_flag",
    "gnomad_observed_flag",
    "gnomad_not_joined_flag",
]

FOLDX_FIELDS = [
    "foldx_ddg_method",
    "foldx_ddg_status",
    "foldx_ddg_missing_reason",
    "foldx_ddg_kcal_mol",
    "foldx_ddg_abs",
    "foldx_ddg_batch_mode",
    "foldx_ddg_mutation",
    "foldx_ddg_uniprot_accession",
    "foldx_ddg_protein_position",
    "foldx_ddg_ref_aa",
    "foldx_ddg_alt_aa",
    "foldx_ddg_alphafold_plddt",
    "foldx_ddg_selected_source_dataset",
    "foldx_ddg_selected_variant_uid",
    "foldx_ddg_source_row_count",
    "foldx_ddg_duplicate_variant_flag",
    "foldx_ddg_discordant_duplicate_flag",
    "foldx_ddg_is_missing",
    "foldx_ddg_ref_mismatch_flag",
]


def rel(path: Path) -> str:
    path = path.resolve()
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def numeric_or_blank(value: str | None) -> str:
    value = str(value or "").strip()
    if not value:
        return ""
    try:
        return repr(float(value))
    except ValueError:
        return ""


def is_nonzero_count(value: str | None) -> bool:
    try:
        return float(str(value or "").strip()) > 0
    except ValueError:
        return False


def first_nonblank(*values: str | None) -> str:
    for value in values:
        text = str(value or "").strip()
        if text and text not in MISSING_VALUES:
            return text
    return ""


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def load_gnomad(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    fields, rows = read_tsv(path)
    return fields, {row["variant_id"]: row for row in rows}


def dbnsfp_gnomad_values(row: dict[str, str]) -> dict[str, str]:
    return {
        "af": first_nonblank(row.get("gnomad41_joint_af_dbnsfp"), row.get("dbnsfp_gnomAD4.1_joint_AF")),
        "popmax_af": first_nonblank(
            row.get("gnomad41_joint_popmax_af_dbnsfp"),
            row.get("dbnsfp_popmax_af"),
            row.get("dbnsfp_gnomAD4.1_joint_POPMAX_AF"),
            row.get("dbnsfp_dbNSFP_POPMAX_AF"),
        ),
        "popmax_pop": first_nonblank(row.get("dbnsfp_popmax_pop"), row.get("dbnsfp_dbNSFP_POPMAX_POP")),
        "homozygote_count": first_nonblank(
            row.get("gnomad41_joint_nhomalt_dbnsfp"),
            row.get("dbnsfp_gnomAD4.1_joint_nhomalt"),
        ),
        "hemizygote_count": "",
    }


def browser_gnomad_values(row: dict[str, str] | None) -> dict[str, str]:
    if not row:
        return {"status": "", "missing_reason": ""}
    return {
        "status": row.get("gnomad_browser_status", ""),
        "missing_reason": row.get("gnomad_browser_missing_reason", ""),
        "af": row.get("gnomad_browser_max_af", ""),
        "popmax_af": row.get("gnomad_browser_max_popmax_af", ""),
        "popmax_pop": first_nonblank(
            row.get("gnomad_browser_joint_popmax_pop"),
            row.get("gnomad_browser_exome_popmax_pop"),
            row.get("gnomad_browser_genome_popmax_pop"),
        ),
        "homozygote_count": row.get("gnomad_browser_max_homozygote_count", ""),
        "hemizygote_count": row.get("gnomad_browser_max_hemizygote_count", ""),
    }


def final_gnomad(row: dict[str, str], browser: dict[str, str] | None) -> dict[str, str]:
    db = dbnsfp_gnomad_values(row)
    br = browser_gnomad_values(browser)
    has_db = any(db.values())
    browser_ok = br.get("status") == "ok"
    browser_not_found = br.get("status") == "not_found"

    out = {field: "" for field in GNOMAD_FINAL_FIELDS}
    if browser_ok:
        source = "browser_and_dbnsfp" if has_db else "browser"
        out.update(
            {
                "gnomad_final_status": "observed",
                "gnomad_final_source": source,
                "gnomad_final_af": first_nonblank(br.get("af"), db.get("af")),
                "gnomad_final_popmax_af": first_nonblank(br.get("popmax_af"), db.get("popmax_af")),
                "gnomad_final_popmax_pop": first_nonblank(br.get("popmax_pop"), db.get("popmax_pop")),
                "gnomad_final_homozygote_count": first_nonblank(
                    br.get("homozygote_count"), db.get("homozygote_count")
                ),
                "gnomad_final_hemizygote_count": first_nonblank(
                    br.get("hemizygote_count"), db.get("hemizygote_count")
                ),
                "gnomad_confirmed_absent_flag": "false",
                "gnomad_observed_flag": "true",
                "gnomad_not_joined_flag": "false",
            }
        )
    elif has_db:
        out.update(
            {
                "gnomad_final_status": "observed",
                "gnomad_final_source": "dbnsfp",
                "gnomad_final_af": db.get("af", ""),
                "gnomad_final_popmax_af": db.get("popmax_af", ""),
                "gnomad_final_popmax_pop": db.get("popmax_pop", ""),
                "gnomad_final_homozygote_count": db.get("homozygote_count", ""),
                "gnomad_final_hemizygote_count": db.get("hemizygote_count", ""),
                "gnomad_confirmed_absent_flag": "false",
                "gnomad_observed_flag": "true",
                "gnomad_not_joined_flag": "false",
            }
        )
    elif browser_not_found:
        out.update(
            {
                "gnomad_final_status": "confirmed_absent",
                "gnomad_final_source": "browser_not_found",
                "gnomad_final_missing_reason": br.get("missing_reason")
                or "exact GRCh38 locus+alleles absent from gnomAD browser v4.1 table",
                "gnomad_final_af": "0",
                "gnomad_final_popmax_af": "0",
                "gnomad_final_homozygote_count": "0",
                "gnomad_final_hemizygote_count": "0",
                "gnomad_confirmed_absent_flag": "true",
                "gnomad_observed_flag": "false",
                "gnomad_not_joined_flag": "false",
            }
        )
    else:
        out.update(
            {
                "gnomad_final_status": "not_joined",
                "gnomad_final_source": "none",
                "gnomad_final_missing_reason": "no dbNSFP gnomAD value and no completed browser query row",
                "gnomad_confirmed_absent_flag": "false",
                "gnomad_observed_flag": "false",
                "gnomad_not_joined_flag": "true",
            }
        )
    out["gnomad_final_af_is_missing"] = bool_text(out.get("gnomad_final_af", "") == "")
    return out


def foldx_rank(row: dict[str, str]) -> tuple[int, int, float]:
    status_rank = {"ok": 0, "ref_mismatch": 1}.get(row.get("ddg_status", ""), 2)
    source_rank = {"hiro": 0, "cardioboost": 1, "emerge": 2}.get(row.get("source_dataset", ""), 9)
    try:
        plddt = -float(row.get("alphafold_residue_plddt", ""))
    except ValueError:
        plddt = 0.0
    return status_rank, source_rank, plddt


def load_foldx_selected(path: Path) -> tuple[dict[str, dict[str, str]], list[dict[str, str]], Counter[str]]:
    _, rows = read_tsv(path)
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["variant_id"]].append(row)

    selected: dict[str, dict[str, str]] = {}
    audit_rows: list[dict[str, str]] = []
    counts: Counter[str] = Counter()
    for variant_id, group in grouped.items():
        group_sorted = sorted(group, key=foldx_rank)
        chosen = group_sorted[0]
        ok_values = {
            numeric_or_blank(row.get("ddg_kcal_mol"))
            for row in group
            if row.get("ddg_status") == "ok" and numeric_or_blank(row.get("ddg_kcal_mol"))
        }
        discordant = len(ok_values) > 1
        duplicate = len(group) > 1
        counts["foldx_variant_rows"] += 1
        counts[f"foldx_selected_status_{chosen.get('ddg_status') or 'blank'}"] += 1
        if duplicate:
            counts["foldx_duplicate_variant_groups"] += 1
        if discordant:
            counts["foldx_discordant_duplicate_groups"] += 1

        selected[variant_id] = {
            "foldx_ddg_method": chosen.get("ddg_method", ""),
            "foldx_ddg_status": chosen.get("ddg_status", ""),
            "foldx_ddg_missing_reason": chosen.get("ddg_missing_reason", ""),
            "foldx_ddg_kcal_mol": chosen.get("ddg_kcal_mol", ""),
            "foldx_ddg_abs": chosen.get("ddg_abs", ""),
            "foldx_ddg_batch_mode": chosen.get("foldx_batch_mode", ""),
            "foldx_ddg_mutation": chosen.get("foldx_mutation", ""),
            "foldx_ddg_uniprot_accession": chosen.get("uniprot_accession", ""),
            "foldx_ddg_protein_position": chosen.get("protein_position", ""),
            "foldx_ddg_ref_aa": chosen.get("protein_ref_aa", ""),
            "foldx_ddg_alt_aa": chosen.get("protein_alt_aa", ""),
            "foldx_ddg_alphafold_plddt": chosen.get("alphafold_residue_plddt", ""),
            "foldx_ddg_selected_source_dataset": chosen.get("source_dataset", ""),
            "foldx_ddg_selected_variant_uid": chosen.get("variant_uid", ""),
            "foldx_ddg_source_row_count": str(len(group)),
            "foldx_ddg_duplicate_variant_flag": bool_text(duplicate),
            "foldx_ddg_discordant_duplicate_flag": bool_text(discordant),
            "foldx_ddg_is_missing": bool_text(chosen.get("ddg_status") != "ok" or not chosen.get("ddg_kcal_mol")),
            "foldx_ddg_ref_mismatch_flag": bool_text(chosen.get("ddg_status") == "ref_mismatch"),
        }
        audit_rows.append({"variant_id": variant_id, **selected[variant_id]})
    return selected, audit_rows, counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=MATRIX)
    parser.add_argument("--gnomad", type=Path, default=GNOMAD)
    parser.add_argument("--foldx", type=Path, default=FOLDX)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--summary", type=Path, default=SUMMARY)
    parser.add_argument("--foldx-selected", type=Path, default=FOLDX_SELECTED)
    args = parser.parse_args()

    gnomad_fields, gnomad = load_gnomad(args.gnomad)
    foldx, foldx_audit_rows, foldx_counts = load_foldx_selected(args.foldx)
    counts: Counter[str] = Counter(foldx_counts)

    gnomad_add_fields = [f for f in gnomad_fields if f not in KEY_FIELDS]
    add_fields = []
    for field in gnomad_add_fields + GNOMAD_FINAL_FIELDS + FOLDX_FIELDS:
        if field not in add_fields:
            add_fields.append(field)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.matrix.open(newline="") as in_handle, args.out.open("w", newline="") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        matrix_fields = list(reader.fieldnames or [])
        fieldnames = []
        for field in matrix_fields + [f for f in add_fields if f not in matrix_fields]:
            if field not in fieldnames:
                fieldnames.append(field)
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()

        seen = set()
        for row in reader:
            counts["matrix_rows"] += 1
            variant_id = row["variant_id"]
            seen.add(variant_id)
            browser = gnomad.get(variant_id)

            if browser:
                counts["gnomad_browser_joined"] += 1
                for field in gnomad_add_fields:
                    row[field] = browser.get(field, "")
            else:
                counts["gnomad_browser_not_joined_to_matrix"] += 1
                for field in gnomad_add_fields:
                    row.setdefault(field, "")

            final_g = final_gnomad(row, browser)
            row.update(final_g)
            counts[f"gnomad_final_status_{final_g['gnomad_final_status']}"] += 1
            if final_g.get("gnomad_final_af"):
                counts["gnomad_final_af_nonmissing"] += 1

            fx = foldx.get(variant_id)
            if fx:
                counts["foldx_joined"] += 1
                row.update(fx)
            else:
                counts["foldx_not_available"] += 1
                row.update({field: "" for field in FOLDX_FIELDS})
                row["foldx_ddg_status"] = "not_available"
                row["foldx_ddg_missing_reason"] = "variant was not eligible for FoldX primary batch or not selected"
                row["foldx_ddg_is_missing"] = "true"
                row["foldx_ddg_ref_mismatch_flag"] = "false"
            writer.writerow(row)

    args.foldx_selected.parent.mkdir(parents=True, exist_ok=True)
    with args.foldx_selected.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["variant_id"] + FOLDX_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(foldx_audit_rows)

    summary = {
        "inputs": {
            "matrix": rel(args.matrix),
            "gnomad_browser_parsed": rel(args.gnomad),
            "foldx_batch": rel(args.foldx),
        },
        "outputs": {
            "merged_matrix": rel(args.out),
            "foldx_variant_selected": rel(args.foldx_selected),
            "summary": rel(args.summary),
        },
        "output_rows": counts["matrix_rows"],
        "counts": dict(sorted(counts.items())),
        "notes": [
            "gnomAD final fields distinguish observed, confirmed_absent, and not_joined states.",
            "Confirmed absent variants receive AF=0 only with gnomad_confirmed_absent_flag=true.",
            "FoldX source-record outputs are collapsed to one row per variant, preferring ok over ref_mismatch.",
            "Previous modeling matrices are unchanged.",
        ],
    }
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
