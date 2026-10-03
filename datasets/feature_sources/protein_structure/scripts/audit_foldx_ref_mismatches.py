#!/usr/bin/env python3
"""Audit FoldX ref-AA mismatch rows and classify likely causes."""

from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]

FOLDX_BATCH = ROOT / "datasets/feature_sources/protein_structure/interim/foldxpro_ddg_batch.tsv"
FOLDX_SELECTED = ROOT / "datasets/modeling/interim/foldx_ddg_variant_selected_features.tsv"
MATRIX = ROOT / "datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx.tsv"

OUT_DIR = ROOT / "datasets/feature_sources/protein_structure/interim"
BATCH_AUDIT = OUT_DIR / "foldx_ref_mismatch_source_row_audit.tsv"
VARIANT_AUDIT = OUT_DIR / "foldx_ref_mismatch_variant_audit.tsv"
SUMMARY = OUT_DIR / "foldx_ref_mismatch_audit.summary.json"

MISMATCH_RE = re.compile(
    r"Variant ref AA (?P<variant_ref>[A-Z]) does not match AlphaFold PDB chain (?P<chain>[^ ]+) "
    r"residue (?P<position>\d+)=(?P<structure_ref>[A-Z])"
)


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def parse_mismatch(reason: str) -> dict[str, str]:
    match = MISMATCH_RE.search(reason or "")
    if not match:
        return {
            "mismatch_variant_ref_aa": "",
            "mismatch_structure_chain": "",
            "mismatch_position": "",
            "mismatch_structure_ref_aa": "",
        }
    return {
        "mismatch_variant_ref_aa": match.group("variant_ref"),
        "mismatch_structure_chain": match.group("chain"),
        "mismatch_position": match.group("position"),
        "mismatch_structure_ref_aa": match.group("structure_ref"),
    }


def recommendation(row: dict[str, str], gene_count: Counter[str], uniprot_count: Counter[str]) -> str:
    gene = row.get("gene", "")
    uniprot = row.get("uniprot_accession", "")
    if gene == "CACNA1C" or uniprot == "Q13936":
        return "exclude_from_primary_ddg; likely CACNA1C transcript/isoform coordinate mismatch; remap with VEP/canonical transcript before rerun"
    if gene_count[gene] >= 3 or uniprot_count[uniprot] >= 3:
        return "exclude_from_primary_ddg; clustered gene/protein mismatch suggests systematic isoform or residue-numbering issue"
    return "exclude_from_primary_ddg; single-row ref-AA mismatch requires manual transcript/protein remap before rerun"


def main() -> int:
    _, batch_rows = read_tsv(FOLDX_BATCH)
    _, selected_rows = read_tsv(FOLDX_SELECTED)
    _, matrix_rows = read_tsv(MATRIX)

    batch_mismatch = [row for row in batch_rows if row.get("ddg_status") == "ref_mismatch"]
    selected_mismatch = [row for row in selected_rows if row.get("foldx_ddg_status") == "ref_mismatch"]
    gene_count = Counter(row.get("gene", "") for row in batch_mismatch)
    uniprot_count = Counter(row.get("uniprot_accession", "") for row in batch_mismatch)
    selected_variant_ids = {row["variant_id"] for row in selected_mismatch}

    matrix_by_variant = {row["variant_id"]: row for row in matrix_rows}
    counts: Counter[str] = Counter()

    source_fields = [
        "variant_id",
        "source_dataset",
        "variant_uid",
        "variant_key",
        "target_3class",
        "gene",
        "uniprot_accession",
        "protein_position",
        "protein_ref_aa",
        "protein_alt_aa",
        "alphafold_residue_plddt",
        "foldx_mutation",
        "ddg_status",
        "ddg_missing_reason",
        "mismatch_variant_ref_aa",
        "mismatch_structure_chain",
        "mismatch_position",
        "mismatch_structure_ref_aa",
        "foldx_mismatch_recommendation",
        "selected_in_variant_matrix",
    ]
    with BATCH_AUDIT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=source_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in batch_mismatch:
            parsed = parse_mismatch(row.get("ddg_missing_reason", ""))
            rec = recommendation(row, gene_count, uniprot_count)
            out = {
                **row,
                **parsed,
                "foldx_mismatch_recommendation": rec,
                "selected_in_variant_matrix": "true" if row["variant_id"] in selected_variant_ids else "false",
            }
            writer.writerow(out)
            counts["source_mismatch_rows"] += 1
            counts[f"source_gene_{row.get('gene') or 'blank'}"] += 1
            counts[f"source_dataset_{row.get('source_dataset') or 'blank'}"] += 1
            if out["selected_in_variant_matrix"] == "true":
                counts["source_rows_selected_variant_has_ref_mismatch"] += 1

    variant_fields = [
        "variant_id",
        "matrix_model_label_3class",
        "matrix_primary_gene",
        "matrix_sources",
        "foldx_ddg_selected_source_dataset",
        "foldx_ddg_selected_variant_uid",
        "foldx_ddg_uniprot_accession",
        "foldx_ddg_protein_position",
        "foldx_ddg_ref_aa",
        "foldx_ddg_alt_aa",
        "foldx_ddg_alphafold_plddt",
        "foldx_ddg_status",
        "foldx_ddg_missing_reason",
        "mismatch_variant_ref_aa",
        "mismatch_structure_chain",
        "mismatch_position",
        "mismatch_structure_ref_aa",
        "foldx_ddg_source_row_count",
        "foldx_ddg_duplicate_variant_flag",
        "foldx_ddg_discordant_duplicate_flag",
        "foldx_mismatch_recommendation",
    ]
    with VARIANT_AUDIT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=variant_fields, delimiter="\t")
        writer.writeheader()
        for row in selected_mismatch:
            matrix = matrix_by_variant.get(row["variant_id"], {})
            parsed = parse_mismatch(row.get("foldx_ddg_missing_reason", ""))
            rec_row = {
                "gene": matrix.get("primary_gene", ""),
                "uniprot_accession": row.get("foldx_ddg_uniprot_accession", ""),
            }
            out = {
                "variant_id": row["variant_id"],
                "matrix_model_label_3class": matrix.get("model_label_3class", ""),
                "matrix_primary_gene": matrix.get("primary_gene", ""),
                "matrix_sources": matrix.get("sources", ""),
                "foldx_ddg_selected_source_dataset": row.get("foldx_ddg_selected_source_dataset", ""),
                "foldx_ddg_selected_variant_uid": row.get("foldx_ddg_selected_variant_uid", ""),
                "foldx_ddg_uniprot_accession": row.get("foldx_ddg_uniprot_accession", ""),
                "foldx_ddg_protein_position": row.get("foldx_ddg_protein_position", ""),
                "foldx_ddg_ref_aa": row.get("foldx_ddg_ref_aa", ""),
                "foldx_ddg_alt_aa": row.get("foldx_ddg_alt_aa", ""),
                "foldx_ddg_alphafold_plddt": row.get("foldx_ddg_alphafold_plddt", ""),
                "foldx_ddg_status": row.get("foldx_ddg_status", ""),
                "foldx_ddg_missing_reason": row.get("foldx_ddg_missing_reason", ""),
                **parsed,
                "foldx_ddg_source_row_count": row.get("foldx_ddg_source_row_count", ""),
                "foldx_ddg_duplicate_variant_flag": row.get("foldx_ddg_duplicate_variant_flag", ""),
                "foldx_ddg_discordant_duplicate_flag": row.get("foldx_ddg_discordant_duplicate_flag", ""),
                "foldx_mismatch_recommendation": recommendation(rec_row, gene_count, uniprot_count),
            }
            writer.writerow(out)
            counts["variant_mismatch_rows"] += 1
            counts[f"variant_gene_{out['matrix_primary_gene'] or 'blank'}"] += 1
            counts[f"variant_source_{out['foldx_ddg_selected_source_dataset'] or 'blank'}"] += 1

    summary = {
        "inputs": {
            "foldx_batch": str(FOLDX_BATCH.relative_to(ROOT)),
            "foldx_selected": str(FOLDX_SELECTED.relative_to(ROOT)),
            "modeling_matrix": str(MATRIX.relative_to(ROOT)),
        },
        "outputs": {
            "source_row_audit": str(BATCH_AUDIT.relative_to(ROOT)),
            "variant_audit": str(VARIANT_AUDIT.relative_to(ROOT)),
            "summary": str(SUMMARY.relative_to(ROOT)),
        },
        "batch_rows": len(batch_rows),
        "selected_variant_rows": len(selected_rows),
        "batch_ref_mismatch_rows": len(batch_mismatch),
        "selected_ref_mismatch_variants": len(selected_mismatch),
        "top_batch_mismatch_genes": gene_count.most_common(20),
        "top_batch_mismatch_uniprots": uniprot_count.most_common(20),
        "counts": dict(sorted(counts.items())),
        "recommendation": (
            "Keep all ref_mismatch rows excluded from numeric FoldX DDG features in the primary model. "
            "Treat them as missing DDG with foldx_ddg_ref_mismatch_flag=true. "
            "Do not impute DDG. Rerun only after transcript/protein remapping, especially for CACNA1C/Q13936."
        ),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
