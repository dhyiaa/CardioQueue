#!/usr/bin/env python3
"""Build one-row-per-gene AF3 coverage review with evidence/citation fields."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
AUDIT_TSV = ROOT / "results/af3/gene_panel_audit.tsv"
SEQ_TSV = ROOT / "datasets/feature_sources/alphafold3/interim/uniprot_gene_to_sequence.tsv"
PAIR_TSV = ROOT / "results/af3/af3_full_panel_partner_manifest.external_evidence.tsv"
COMBINED_PAIR_TSV = ROOT / "results/af3/af3_partner_combined_evidence_review_queue.tsv"
DESIGN_TSV = ROOT / "results/af3/af3_complex_design_review.tsv"
JOB_TSV = ROOT / "datasets/feature_sources/alphafold3/jobs/af3_job_manifest.tsv"
OUT_TSV = ROOT / "results/af3/af3_gene_level_coverage_review.tsv"
SUMMARY_JSON = ROOT / "results/af3/af3_gene_level_coverage_review.summary.json"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def split_csv(value: str) -> list[str]:
    return [v.strip() for v in value.split(",") if v.strip()]


def split_semicolon(value: str) -> list[str]:
    return [v.strip() for v in value.split(";") if v.strip()]


def unique_join(values: list[str]) -> str:
    return ";".join(sorted(set(v for v in values if v)))


def gene_fix_needed(status: str, design_modes: str) -> str:
    if status == "packet_draft_ready":
        return "Review draft AF3 packet manually, then submit if approved."
    if status == "assigned_hold_for_manual_review":
        return "Resolve evidence/domain notes before packet. Do not submit yet."
    if status == "assigned_domain_or_split_required":
        return "Choose domain windows or split complex into smaller biologically specific AF3 jobs."
    if status == "no_partner_assigned":
        return "Manually review candidate partners and promote only biologically defensible cardiac partners, or mark no_curated_partner."
    if status == "excluded_or_quarantine":
        return "No AF3 work unless scope changes."
    if status == "ttn_domain_only_or_defer":
        return "TTN requires separate domain-level strategy; full-length AF3 is not feasible."
    return "Manual review needed."


def main() -> None:
    audit_rows = read_tsv(AUDIT_TSV)
    seq_rows = read_tsv(SEQ_TSV)
    pair_rows = read_tsv(PAIR_TSV)
    combined_rows = read_tsv(COMBINED_PAIR_TSV)
    design_rows = read_tsv(DESIGN_TSV)
    job_rows = read_tsv(JOB_TSV) if JOB_TSV.exists() else []

    seq_by_gene = {r["gene"]: r for r in seq_rows}
    audit_by_gene = {r["gene_normalized"]: r for r in audit_rows}

    # Only AF3 review genes, not excluded/quarantine or TTN defer.
    target_genes = sorted(
        r["gene_normalized"]
        for r in audit_rows
        if r["af3_inclusion_decision"].startswith("include")
    )

    design_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in design_rows:
        for gene in split_csv(row["genes_covered"]):
            design_by_gene[gene].append(row)

    job_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in job_rows:
        for gene in split_csv(row["genes_covered"]):
            job_by_gene[gene].append(row)

    pair_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in pair_rows:
        pair_by_gene[row["target_gene"]].append(row)

    combined_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in combined_rows:
        combined_by_gene[row["target_gene"]].append(row)

    out_rows = []
    for gene in target_genes:
        seq = seq_by_gene.get(gene, {})
        designs = design_by_gene.get(gene, [])
        jobs = job_by_gene.get(gene, [])
        pairs = pair_by_gene.get(gene, [])
        combined = combined_by_gene.get(gene, [])
        audit = audit_by_gene.get(gene, {})

        if jobs:
            status = "packet_draft_ready"
        elif designs and any(d["packet_draft_status"] == "hold_for_manual_review" for d in designs):
            status = "assigned_hold_for_manual_review"
        elif designs and any(d["packet_draft_status"] == "not_ready_for_packet" for d in designs):
            status = "assigned_domain_or_split_required"
        elif pairs and all(p["partner_type"] == "none_found" for p in pairs):
            status = "no_partner_assigned"
        elif designs:
            status = "assigned_manual_review_needed"
        else:
            status = "no_partner_assigned"

        partner_entities = []
        partner_types = []
        complex_ids = []
        complex_names = []
        support_levels = []
        pubmed_ids = []
        intact_pubmed_ids = []
        reactome_ids = []
        pdb_ids = []
        intact_ids = []
        evidence_source_urls = set()

        for row in combined:
            if row["partner_gene_or_entity"]:
                partner_entities.append(row["partner_gene_or_entity"])
            if row["partner_type"]:
                partner_types.append(row["partner_type"])
            if row["complex_group_id"]:
                complex_ids.append(row["complex_group_id"])
            if row["complex_name"]:
                complex_names.append(row["complex_name"])
            support_levels.append(row["combined_support_level"])
            pubmed_ids.extend(split_semicolon(row["uniprot_pair_pubmed_ids"]))
            intact_pubmed_ids.extend(split_semicolon(row["intact_psicquic_pair_pubmed_ids"]))
            reactome_ids.extend(split_semicolon(row["reactome_shared_pathway_ids"]))
            pdb_ids.extend(split_semicolon(row["rcsb_shared_pdb_ids"]))
            if row["combined_support_level"] in {"direct_uniprot_intact_interaction", "uniprot_subunit_text_mention"}:
                evidence_source_urls.add("https://rest.uniprot.org/uniprotkb/")
            if row["external_structured_support_level"] == "intact_psicquic_pair_support":
                evidence_source_urls.add("https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/")
            if row["reactome_shared_pathway_ids"]:
                evidence_source_urls.add("https://reactome.org/ContentService/")
            if row["rcsb_shared_pdb_ids"]:
                evidence_source_urls.add("https://search.rcsb.org/rcsbsearch/v2/query")

        for row in pairs:
            intact_ids.extend(split_semicolon(row.get("uniprot_pair_intact_ids", "")))

        design_modes = unique_join([d["design_mode_recommendation"] for d in designs])
        packet_statuses = unique_join([d["packet_draft_status"] for d in designs])
        job_ids = unique_join([j["job_id"] for j in jobs])
        packet_dirs = unique_join([j["packet_dir"] for j in jobs])

        out_rows.append({
            "gene": gene,
            "panel_bucket": audit.get("panel_bucket", ""),
            "primary_model_inclusion_values": audit.get("primary_model_inclusion_values", ""),
            "n_total_variants": audit.get("n_total_variants", ""),
            "n_missense": audit.get("n_missense", ""),
            "protein_sequence_status": "mapped_reviewed_uniprot" if seq.get("uniprot_accession") else "missing",
            "uniprot_accession": seq.get("uniprot_accession", ""),
            "entry_name": seq.get("entry_name", ""),
            "protein_name": seq.get("protein_name", ""),
            "sequence_length": seq.get("sequence_length", ""),
            "sequence_source": seq.get("retrieval_source", ""),
            "uniprot_url": f"https://rest.uniprot.org/uniprotkb/{seq.get('uniprot_accession', '')}.json" if seq.get("uniprot_accession") else "",
            "partner_assignment_status": status,
            "complex_group_ids": unique_join(complex_ids),
            "complex_names": unique_join(complex_names),
            "partner_entities": unique_join(partner_entities),
            "partner_types": unique_join(partner_types),
            "combined_support_levels": unique_join(support_levels),
            "pubmed_ids_from_uniprot": unique_join(pubmed_ids),
            "pubmed_ids_from_intact": unique_join(intact_pubmed_ids),
            "intact_ids_from_uniprot": unique_join(intact_ids),
            "reactome_ids": unique_join(reactome_ids),
            "pdb_ids": unique_join(pdb_ids),
            "evidence_source_urls": unique_join(list(evidence_source_urls)),
            "design_mode_recommendation": design_modes,
            "packet_draft_status": packet_statuses,
            "draft_job_ids": job_ids,
            "draft_packet_dirs": packet_dirs,
            "fix_needed": gene_fix_needed(status, design_modes),
        })

    fields = [
        "gene", "panel_bucket", "primary_model_inclusion_values",
        "n_total_variants", "n_missense",
        "protein_sequence_status", "uniprot_accession", "entry_name",
        "protein_name", "sequence_length", "sequence_source", "uniprot_url",
        "partner_assignment_status", "complex_group_ids", "complex_names",
        "partner_entities", "partner_types", "combined_support_levels",
        "pubmed_ids_from_uniprot", "pubmed_ids_from_intact",
        "intact_ids_from_uniprot", "reactome_ids", "pdb_ids",
        "evidence_source_urls", "design_mode_recommendation",
        "packet_draft_status", "draft_job_ids", "draft_packet_dirs",
        "fix_needed",
    ]
    write_tsv(OUT_TSV, out_rows, fields)

    summary = {
        "date": date.today().isoformat(),
        "af3_target_gene_count": len(target_genes),
        "genes_with_reviewed_uniprot_sequence": sum(r["protein_sequence_status"] == "mapped_reviewed_uniprot" for r in out_rows),
        "partner_assignment_status_counts": dict(Counter(r["partner_assignment_status"] for r in out_rows)),
        "genes_packet_draft_ready": [r["gene"] for r in out_rows if r["partner_assignment_status"] == "packet_draft_ready"],
        "genes_assigned_but_not_packet_ready": [r["gene"] for r in out_rows if r["partner_assignment_status"] in {"assigned_hold_for_manual_review", "assigned_domain_or_split_required", "assigned_manual_review_needed"}],
        "genes_without_partner_assignment": [r["gene"] for r in out_rows if r["partner_assignment_status"] == "no_partner_assigned"],
        "output": str(OUT_TSV),
        "source_urls": [
            "https://rest.uniprot.org/uniprotkb/",
            "https://reactome.org/ContentService/",
            "https://search.rcsb.org/rcsbsearch/v2/query",
            "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/",
        ],
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
