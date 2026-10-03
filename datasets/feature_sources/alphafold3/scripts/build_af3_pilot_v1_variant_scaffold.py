#!/usr/bin/env python3
"""Create a variant-level AF3 pilot-v1 scaffold without merging into the model matrix."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
MODEL_TSV = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
PILOT_GENE_SET_TSV = ROOT / "results/af3/pilot_v1/af3_pilot_v1_gene_set.tsv"
PILOT_COMPLEX_PLAN_TSV = ROOT / "results/af3/pilot_v1/af3_pilot_v1_complex_plan.tsv"
OUT_TSV = ROOT / "datasets/feature_sources/alphafold3/processed/pilot_v1/af3_pilot_v1_variant_feature_scaffold.tsv"
SUMMARY_JSON = ROOT / "results/af3/pilot_v1/af3_pilot_v1_variant_scaffold.summary.json"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def split_values(value: str, sep: str = ",") -> list[str]:
    return [v.strip() for v in value.split(sep) if v.strip()]


def split_sources(value: str) -> list[str]:
    value = value.replace(";", "|")
    return [v.strip() for v in value.split("|") if v.strip()]


def protein_position(row: dict[str, str]) -> str:
    for col in ["protein_position", "vep_protein_position", "dbnsfp_aapos", "foldx_ddg_protein_position"]:
        value = row.get(col, "")
        if value and value not in {".", "NA", "nan"}:
            return value
    return ""


def main() -> None:
    gene_rows = read_tsv(PILOT_GENE_SET_TSV)
    complex_rows = read_tsv(PILOT_COMPLEX_PLAN_TSV)
    gene_to_complexes: dict[str, list[str]] = defaultdict(list)
    gene_to_stages: dict[str, list[str]] = defaultdict(list)
    gene_to_statuses: dict[str, list[str]] = defaultdict(list)

    for row in complex_rows:
        for gene in split_values(row["target_genes"]):
            gene_to_complexes[gene].append(row["complex_group_id"])
            gene_to_stages[gene].append(row["pilot_stage"])
            gene_to_statuses[gene].append(row["submission_status"])

    pilot_genes = {row["gene"] for row in gene_rows}
    rows = []
    counts = Counter()
    source_counts = Counter()
    label_counts = Counter()
    with MODEL_TSV.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene = row.get("primary_gene", "")
            covered = gene in pilot_genes
            pos = protein_position(row)
            if covered and pos:
                missing_reason = "af3_server_output_not_yet_available"
            elif covered:
                missing_reason = "protein_position_missing_before_af3_mapping"
            else:
                missing_reason = "not_af3_pilot_v1_gene"

            out = {
                "variant_id": row["variant_id"],
                "primary_gene": gene,
                "model_label_3class": row.get("model_label_3class", ""),
                "primary_model_inclusion": row.get("primary_model_inclusion", ""),
                "sources": row.get("sources", ""),
                "protein_position": pos,
                "uniprot_accession": row.get("uniprot_accession", ""),
                "protein_ref_aa": row.get("protein_ref_aa", ""),
                "protein_alt_aa": row.get("protein_alt_aa", ""),
                "af3_pilot_v1_gene_covered": "1" if covered else "0",
                "af3_pilot_v1_complex_count_covering_gene": str(len(set(gene_to_complexes.get(gene, [])))),
                "af3_pilot_v1_complex_ids": ";".join(sorted(set(gene_to_complexes.get(gene, [])))),
                "af3_pilot_v1_stage": ";".join(sorted(set(gene_to_stages.get(gene, [])))),
                "af3_pilot_v1_submission_status": ";".join(sorted(set(gene_to_statuses.get(gene, [])))),
                "af3_server_output_available": "0",
                "af3_residue_mapped_to_model": "0",
                "af3_feature_missing_reason": missing_reason,
            }
            rows.append(out)
            counts["total_variants"] += 1
            counts["covered_variants"] += int(covered)
            counts["covered_with_protein_position"] += int(covered and bool(pos))
            counts["covered_without_protein_position"] += int(covered and not bool(pos))
            if covered:
                label_counts[row.get("model_label_3class", "")] += 1
                for source in split_sources(row.get("sources", "")):
                    source_counts[source] += 1

    fields = [
        "variant_id",
        "primary_gene",
        "model_label_3class",
        "primary_model_inclusion",
        "sources",
        "protein_position",
        "uniprot_accession",
        "protein_ref_aa",
        "protein_alt_aa",
        "af3_pilot_v1_gene_covered",
        "af3_pilot_v1_complex_count_covering_gene",
        "af3_pilot_v1_complex_ids",
        "af3_pilot_v1_stage",
        "af3_pilot_v1_submission_status",
        "af3_server_output_available",
        "af3_residue_mapped_to_model",
        "af3_feature_missing_reason",
    ]
    write_tsv(OUT_TSV, rows, fields)

    summary = {
        "date": date.today().isoformat(),
        "input_model_matrix": str(MODEL_TSV),
        "output": str(OUT_TSV),
        "pilot_gene_count": len(pilot_genes),
        "counts": dict(counts),
        "covered_label_counts": dict(label_counts),
        "covered_source_counts": dict(source_counts),
        "safety_note": "This is a standalone AF3 scaffold. It was not merged into the existing modeling table.",
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
