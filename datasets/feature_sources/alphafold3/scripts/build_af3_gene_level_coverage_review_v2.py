#!/usr/bin/env python3
"""Merge first-wave and second-wave AF3 gene coverage decisions."""

from __future__ import annotations

import csv
import json
from collections import Counter
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
FIRST_TSV = ROOT / "results/af3/af3_gene_level_coverage_review.tsv"
SECOND_TSV = ROOT / "results/af3/af3_second_wave_curated_partner_decisions.tsv"
OUT_TSV = ROOT / "results/af3/af3_gene_level_coverage_review_v2.tsv"
SUMMARY_JSON = ROOT / "results/af3/af3_gene_level_coverage_review_v2.summary.json"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def v2_status(first_status: str, second: dict[str, str] | None) -> str:
    if first_status == "packet_draft_ready":
        return "packet_draft_ready"
    if first_status == "assigned_hold_for_manual_review":
        return "assigned_hold_for_manual_review"
    if first_status == "assigned_domain_or_split_required":
        return "assigned_domain_or_split_required"
    if not second:
        return first_status
    manual = second["manual_review_status"]
    if manual == "ready_for_evidence_check":
        return "second_wave_packet_candidate_needs_final_evidence_check"
    if manual == "needs_ligand_ccd_review":
        return "second_wave_ligand_or_homomer_context_needs_ccd_review"
    if manual == "needs_domain_range_selection":
        return "second_wave_domain_review_required"
    if manual == "needs_cardiac_specificity_review":
        return "second_wave_needs_cardiac_specificity_review"
    if manual == "needs_partner_selection":
        return "second_wave_needs_partner_selection"
    if manual == "defer_until_rna_complex_design":
        return "defer_until_rna_complex_design"
    return "second_wave_manual_review_needed"


def main() -> None:
    first_rows = read_tsv(FIRST_TSV)
    second_rows = read_tsv(SECOND_TSV)
    second_by_gene = {r["target_gene"]: r for r in second_rows}

    out = []
    for row in first_rows:
        gene = row["gene"]
        second = second_by_gene.get(gene)
        merged_status = v2_status(row["partner_assignment_status"], second)
        merged = dict(row)
        merged["partner_assignment_status_v2"] = merged_status
        merged["second_wave_complex_group_id"] = second.get("complex_group_id", "") if second else ""
        merged["second_wave_complex_name"] = second.get("complex_name", "") if second else ""
        merged["second_wave_protein_entities"] = second.get("protein_entities", "") if second else ""
        merged["second_wave_nonprotein_entities"] = second.get("nonprotein_entities", "") if second else ""
        merged["second_wave_estimated_tokens"] = second.get("estimated_tokens", "") if second else ""
        merged["second_wave_evidence_anchor"] = second.get("evidence_anchor", "") if second else ""
        merged["second_wave_manual_review_status"] = second.get("manual_review_status", "") if second else ""
        merged["second_wave_pubmed_ids_from_uniprot_gene_records"] = second.get("pubmed_ids_from_uniprot_gene_records", "") if second else ""
        merged["second_wave_pubmed_ids_from_intact_gene_records"] = second.get("pubmed_ids_from_intact_gene_records", "") if second else ""
        merged["second_wave_reactome_ids_from_gene_records"] = second.get("reactome_ids_from_gene_records", "") if second else ""
        merged["second_wave_pdb_ids_from_gene_records"] = second.get("pdb_ids_from_gene_records", "") if second else ""
        merged["second_wave_rationale"] = second.get("rationale", "") if second else ""
        if second:
            merged["fix_needed_v2"] = {
                "second_wave_packet_candidate_needs_final_evidence_check": "Verify citations/evidence rows, then generate draft AF3 packet.",
                "second_wave_ligand_or_homomer_context_needs_ccd_review": "Verify ligand/ion CCD handling or homomer stoichiometry, then generate packet.",
                "second_wave_domain_review_required": "Select domain windows before packet generation.",
                "second_wave_needs_cardiac_specificity_review": "Confirm partner is cardiac-relevant enough before packet.",
                "second_wave_needs_partner_selection": "Choose specific direct partner or keep as no_curated_partner.",
                "defer_until_rna_complex_design": "Defer until RNA/protein spliceosome design is explicitly scoped.",
            }.get(merged_status, row["fix_needed"])
        else:
            merged["fix_needed_v2"] = row["fix_needed"]
        out.append(merged)

    fields = list(first_rows[0].keys()) + [
        "partner_assignment_status_v2",
        "second_wave_complex_group_id",
        "second_wave_complex_name",
        "second_wave_protein_entities",
        "second_wave_nonprotein_entities",
        "second_wave_estimated_tokens",
        "second_wave_evidence_anchor",
        "second_wave_manual_review_status",
        "second_wave_pubmed_ids_from_uniprot_gene_records",
        "second_wave_pubmed_ids_from_intact_gene_records",
        "second_wave_reactome_ids_from_gene_records",
        "second_wave_pdb_ids_from_gene_records",
        "second_wave_rationale",
        "fix_needed_v2",
    ]
    write_tsv(OUT_TSV, out, fields)

    summary = {
        "date": date.today().isoformat(),
        "af3_target_gene_count": len(out),
        "genes_with_reviewed_uniprot_sequence": sum(r["protein_sequence_status"] == "mapped_reviewed_uniprot" for r in out),
        "v1_status_counts": dict(Counter(r["partner_assignment_status"] for r in out)),
        "v2_status_counts": dict(Counter(r["partner_assignment_status_v2"] for r in out)),
        "genes_without_any_partner_path_v2": [r["gene"] for r in out if r["partner_assignment_status_v2"] == "no_partner_assigned"],
        "genes_with_packet_or_packet_candidate_path": [
            r["gene"] for r in out
            if r["partner_assignment_status_v2"] in {
                "packet_draft_ready",
                "second_wave_packet_candidate_needs_final_evidence_check",
                "second_wave_ligand_or_homomer_context_needs_ccd_review",
            }
        ],
        "genes_requiring_domain_or_split_review": [
            r["gene"] for r in out
            if r["partner_assignment_status_v2"] in {
                "assigned_domain_or_split_required",
                "second_wave_domain_review_required",
            }
        ],
        "output": str(OUT_TSV),
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
