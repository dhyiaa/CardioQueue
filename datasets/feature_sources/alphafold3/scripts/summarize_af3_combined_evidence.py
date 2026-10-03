#!/usr/bin/env python3
"""Create combined AF3 evidence readiness summaries."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
FINAL_PAIR_TSV = ROOT / "results/af3/af3_full_panel_partner_manifest.final_evidence_stack.tsv"
PAIR_TSV = FINAL_PAIR_TSV if FINAL_PAIR_TSV.exists() else ROOT / "results/af3/af3_full_panel_partner_manifest.external_evidence.tsv"
COMPLEX_TSV = ROOT / "results/af3/af3_complex_combined_evidence_readiness.tsv"
QUEUE_TSV = ROOT / "results/af3/af3_partner_combined_evidence_review_queue.tsv"
SUMMARY_JSON = ROOT / "results/af3/af3_combined_evidence_readiness.summary.json"

DIRECT_OR_STRUCTURAL = {
    "direct_uniprot_intact_interaction",
    "uniprot_subunit_text_mention",
    "intact_psicquic_pair_support",
    "rcsb_shared_pdb_support",
}

CONTEXT_ONLY = {
    "uniprot_calcium_context",
    "uniprot_gtp_context",
    "reactome_shared_pathway_context",
    "reactome_target_context_for_nonprotein_entity",
}

REVIEW_NEEDED = {
    "calmodulin_family_proxy_review_needed",
    "calmodulin_family_context_review_needed",
    "nonprotein_entity_review_needed",
    "not_confirmed_by_uniprot_yet",
    "not_supported_by_external_structured_sources_yet",
    "no_partner_assigned",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def combined_support(row: dict[str, str]) -> str:
    uni = row["uniprot_pair_support_level"]
    ext = row["external_structured_support_level"]
    if uni in {"direct_uniprot_intact_interaction", "uniprot_subunit_text_mention"}:
        return uni
    if ext in {"intact_psicquic_pair_support", "rcsb_shared_pdb_support"}:
        return ext
    if uni in {"uniprot_calcium_context", "uniprot_gtp_context"}:
        return uni
    if ext in {"reactome_shared_pathway_context", "reactome_target_context_for_nonprotein_entity"}:
        return ext
    if uni.startswith("calmodulin_family"):
        return uni
    if row["partner_type"] == "none_found":
        return "no_partner_assigned"
    return "needs_manual_literature_or_domain_review"


def readiness_for_complex(rows: list[dict[str, str]]) -> str:
    supports = [r["combined_support_level"] for r in rows]
    n_direct = sum(s in DIRECT_OR_STRUCTURAL for s in supports)
    n_context = sum(s in CONTEXT_ONLY for s in supports)
    n_unresolved = sum(s in REVIEW_NEEDED or s == "needs_manual_literature_or_domain_review" for s in supports)
    protein_rows = [r for r in rows if r["partner_type"] == "protein"]
    protein_direct = sum(r["combined_support_level"] in DIRECT_OR_STRUCTURAL for r in protein_rows)
    if n_unresolved == 0 and n_direct > 0:
        return "structured_evidence_ready_for_manual_af3_design"
    if protein_direct > 0 and n_context > 0:
        return "mostly_supported_needs_domain_and_literature_review"
    if protein_direct > 0:
        return "partially_supported_needs_context_or_domain_review"
    if n_context > 0:
        return "context_only_needs_direct_pair_or_literature_evidence"
    return "not_ready_needs_partner_or_evidence"


def next_action(row: dict[str, str]) -> str:
    s = row["combined_support_level"]
    if s in DIRECT_OR_STRUCTURAL:
        return "manual_check_publication_and_choose_full_length_or_domain_ranges"
    if s in CONTEXT_ONLY:
        return "keep_as_context_only_then_seek_direct_pair_or_structure_evidence"
    if s.startswith("calmodulin_family"):
        return "accept_single_CALM_sequence_only_after_manual_isoform_note"
    if s == "no_partner_assigned":
        return "review_candidate_partner_table_or assign no_curated_partner"
    return "targeted_literature_review_required"


def main() -> None:
    rows = read_tsv(PAIR_TSV)
    for row in rows:
        row["combined_support_level"] = combined_support(row)
        row["combined_next_action"] = next_action(row)

    queue = sorted(rows, key=lambda r: (
        r["combined_support_level"] in DIRECT_OR_STRUCTURAL,
        r["combined_support_level"] in CONTEXT_ONLY,
        r["complex_group_id"],
        r["target_gene"],
        r["partner_gene_or_entity"],
    ))
    queue_fields = [
        "target_gene",
        "complex_group_id",
        "complex_name",
        "partner_gene_or_entity",
        "partner_type",
        "uniprot_pair_support_level",
        "external_structured_support_level",
        "combined_support_level",
        "uniprot_pair_pubmed_ids",
        "intact_psicquic_pair_pubmed_ids",
        "reactome_shared_pathway_ids",
        "rcsb_shared_pdb_ids",
        "target_hpa_status",
        "target_hpa_subcellular_main_location",
        "target_hpa_protein_tissue_specificity",
        "target_hpa_protein_tissue_distribution",
        "partner_hpa_status",
        "partner_hpa_subcellular_main_location",
        "partner_hpa_protein_tissue_specificity",
        "partner_hpa_protein_tissue_distribution",
        "hpa_pair_localization_note",
        "combined_next_action",
    ]
    write_tsv(QUEUE_TSV, queue, queue_fields)

    by_complex: dict[str, list[dict[str, str]]] = defaultdict(list)
    no_complex = []
    for row in rows:
        if row["complex_group_id"]:
            by_complex[row["complex_group_id"]].append(row)
        else:
            no_complex.append(row)

    complex_rows = []
    for complex_id, cr in sorted(by_complex.items()):
        support_counts = Counter(r["combined_support_level"] for r in cr)
        readiness = readiness_for_complex(cr)
        complex_rows.append({
            "complex_group_id": complex_id,
            "complex_name": cr[0]["complex_name"],
            "targets_covered": ",".join(sorted({r["target_gene"] for r in cr})),
            "partners_or_entities": ",".join(sorted({r["partner_gene_or_entity"] for r in cr if r["partner_gene_or_entity"]})),
            "n_rows": str(len(cr)),
            "n_direct_or_structural_support": str(sum(v for k, v in support_counts.items() if k in DIRECT_OR_STRUCTURAL)),
            "n_context_only_support": str(sum(v for k, v in support_counts.items() if k in CONTEXT_ONLY)),
            "n_review_needed": str(sum(v for k, v in support_counts.items() if k in REVIEW_NEEDED or k == "needs_manual_literature_or_domain_review")),
            "combined_support_counts_json": json.dumps(dict(support_counts), sort_keys=True),
            "combined_readiness_status": readiness,
            "next_action": "do targeted literature/domain review before final AF3 packet",
        })

    complex_fields = [
        "complex_group_id",
        "complex_name",
        "targets_covered",
        "partners_or_entities",
        "n_rows",
        "n_direct_or_structural_support",
        "n_context_only_support",
        "n_review_needed",
        "combined_support_counts_json",
        "combined_readiness_status",
        "next_action",
    ]
    write_tsv(COMPLEX_TSV, complex_rows, complex_fields)

    summary = {
        "date": date.today().isoformat(),
        "pair_rows": len(rows),
        "complex_count": len(complex_rows),
        "no_complex_gene_rows": len(no_complex),
        "combined_pair_support_counts": dict(Counter(r["combined_support_level"] for r in rows)),
        "complex_readiness_counts": dict(Counter(r["combined_readiness_status"] for r in complex_rows)),
        "outputs": {
            "input_pair_tsv": str(PAIR_TSV),
            "combined_complex_readiness_tsv": str(COMPLEX_TSV),
            "combined_pair_review_queue_tsv": str(QUEUE_TSV),
        },
        "hpa_note": "HPA localization/expression fields are included when available, but are not counted as direct interaction evidence.",
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
