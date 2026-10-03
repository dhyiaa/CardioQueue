#!/usr/bin/env python3
"""Summarize AF3 partner evidence into review queues."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
ENRICHED_TSV = ROOT / "results/af3/af3_full_panel_partner_manifest.uniprot_evidence.tsv"
COMPLEX_TSV = ROOT / "results/af3/af3_complex_evidence_readiness.tsv"
QUEUE_TSV = ROOT / "results/af3/af3_partner_evidence_review_queue.tsv"
SUMMARY_JSON = ROOT / "results/af3/af3_evidence_review.summary.json"

STRONG_SUPPORT = {"direct_uniprot_intact_interaction", "uniprot_subunit_text_mention"}
CONTEXT_SUPPORT = {"uniprot_calcium_context", "uniprot_gtp_context"}
REVIEW_SUPPORT = {
    "calmodulin_family_proxy_review_needed",
    "calmodulin_family_context_review_needed",
    "nonprotein_entity_review_needed",
    "not_confirmed_by_uniprot_yet",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def recommend_next_source(row: dict[str, str]) -> str:
    level = row["uniprot_pair_support_level"]
    partner_type = row["partner_type"]
    complex_id = row["complex_group_id"]
    if level in STRONG_SUPPORT:
        return "manual_verify_pubmed_and_domain_ranges"
    if level in CONTEXT_SUPPORT:
        return "manual_verify_nonprotein_entity_and_af3_ccd_code"
    if level.startswith("calmodulin_family"):
        return "verify_calmodulin_isoform_proxy_and_use_single_calm_sequence"
    if partner_type == "none_found":
        return "search_uniprot_reactome_intact_rcsb_for_partner"
    if "TROPONIN" in complex_id or "THIN_FILAMENT" in complex_id:
        return "check_reactome_or_rcsb_thin_filament_complex"
    if "CASQ2_TRDN_RYR2" in complex_id:
        return "check_reactome_intact_and_literature_sr_calcium_release_complex"
    if "DESMOSOME" in complex_id:
        return "check_reactome_intact_rcsb_desmosome_pair"
    if "RAS_MAPK" in complex_id:
        return "split_into_small_ras_raf_mek_ptpn11_domain_jobs_after_literature_review"
    return "check_uniprot_reactome_intact_rcsb_literature"


def main() -> None:
    rows = read_tsv(ENRICHED_TSV)

    queue_rows = []
    for row in rows:
        level = row["uniprot_pair_support_level"]
        if level in STRONG_SUPPORT or level in CONTEXT_SUPPORT:
            priority = "review_confirmed_support"
        elif row["partner_type"] == "none_found":
            priority = "needs_partner_assignment"
        else:
            priority = "needs_external_evidence"
        out = {
            "review_priority": priority,
            "target_gene": row["target_gene"],
            "complex_group_id": row["complex_group_id"],
            "complex_name": row["complex_name"],
            "partner_gene_or_entity": row["partner_gene_or_entity"],
            "partner_type": row["partner_type"],
            "uniprot_pair_support_level": level,
            "uniprot_pair_pubmed_ids": row["uniprot_pair_pubmed_ids"],
            "uniprot_pair_intact_ids": row["uniprot_pair_intact_ids"],
            "recommended_next_source": recommend_next_source(row),
        }
        queue_rows.append(out)

    queue_rows.sort(key=lambda r: (
        {"needs_external_evidence": 0, "needs_partner_assignment": 1, "review_confirmed_support": 2}[r["review_priority"]],
        r["complex_group_id"],
        r["target_gene"],
        r["partner_gene_or_entity"],
    ))
    queue_fields = [
        "review_priority",
        "target_gene",
        "complex_group_id",
        "complex_name",
        "partner_gene_or_entity",
        "partner_type",
        "uniprot_pair_support_level",
        "uniprot_pair_pubmed_ids",
        "uniprot_pair_intact_ids",
        "recommended_next_source",
    ]
    write_tsv(QUEUE_TSV, queue_rows, queue_fields)

    by_complex: dict[str, list[dict[str, str]]] = defaultdict(list)
    no_complex_rows = []
    for row in rows:
        if row["complex_group_id"]:
            by_complex[row["complex_group_id"]].append(row)
        else:
            no_complex_rows.append(row)

    complex_rows = []
    for complex_id, complex_rows_raw in sorted(by_complex.items()):
        support = Counter(r["uniprot_pair_support_level"] for r in complex_rows_raw)
        targets = sorted({r["target_gene"] for r in complex_rows_raw})
        partners = sorted({r["partner_gene_or_entity"] for r in complex_rows_raw if r["partner_gene_or_entity"]})
        n_strong = sum(support[k] for k in STRONG_SUPPORT)
        n_context = sum(support[k] for k in CONTEXT_SUPPORT)
        n_needs = sum(support[k] for k in REVIEW_SUPPORT)
        if n_needs == 0 and n_strong > 0:
            readiness = "uniprot_supported_ready_for_manual_final_check"
        elif n_strong > 0 or n_context > 0:
            readiness = "partially_supported_needs_external_evidence"
        else:
            readiness = "needs_external_evidence_before_af3"
        complex_rows.append({
            "complex_group_id": complex_id,
            "complex_name": complex_rows_raw[0]["complex_name"],
            "targets_covered": ",".join(targets),
            "partners_or_entities": ",".join(partners),
            "n_rows": str(len(complex_rows_raw)),
            "n_direct_uniprot_intact": str(support["direct_uniprot_intact_interaction"]),
            "n_uniprot_subunit_text": str(support["uniprot_subunit_text_mention"]),
            "n_uniprot_nonprotein_context": str(n_context),
            "n_needing_external_evidence_or_manual_review": str(n_needs),
            "support_counts_json": json.dumps(dict(support), sort_keys=True),
            "readiness_status": readiness,
            "next_action": "fill Reactome/IntAct/RCSB/literature evidence and domain ranges before AF3 jobs",
        })

    complex_fields = [
        "complex_group_id",
        "complex_name",
        "targets_covered",
        "partners_or_entities",
        "n_rows",
        "n_direct_uniprot_intact",
        "n_uniprot_subunit_text",
        "n_uniprot_nonprotein_context",
        "n_needing_external_evidence_or_manual_review",
        "support_counts_json",
        "readiness_status",
        "next_action",
    ]
    write_tsv(COMPLEX_TSV, complex_rows, complex_fields)

    summary = {
        "date": date.today().isoformat(),
        "complex_count": len(complex_rows),
        "genes_without_complex_count": len(no_complex_rows),
        "queue_counts": dict(Counter(r["review_priority"] for r in queue_rows)),
        "complex_readiness_counts": dict(Counter(r["readiness_status"] for r in complex_rows)),
        "outputs": {
            "review_queue_tsv": str(QUEUE_TSV),
            "complex_readiness_tsv": str(COMPLEX_TSV),
        },
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
