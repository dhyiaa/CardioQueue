#!/usr/bin/env python3
"""
Build AF3 sequence and starter partner manifests.

This script only writes into the AlphaFold 3 branch folders. It does not touch
the existing modeling matrices or previous protein-structure outputs.
"""

from __future__ import annotations

import csv
import json
import re
from collections import defaultdict
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
AUDIT_TSV = ROOT / "results/af3/gene_panel_audit.tsv"
UNIPROT_FASTA = (
    ROOT / "datasets/feature_sources/protein_structure/raw/"
    "uniprot_human_reviewed_2026_02.fasta"
)
OUT_DIR = ROOT / "datasets/feature_sources/alphafold3/interim"
RESULTS_DIR = ROOT / "results/af3"

SEQUENCE_TSV = OUT_DIR / "uniprot_gene_to_sequence.tsv"
ENTITIES_FASTA = OUT_DIR / "af3_all_entities.fasta"
WHITELIST_TSV = OUT_DIR / "af3_manual_cardiac_complex_whitelist.tsv"
PARTNER_TSV = RESULTS_DIR / "af3_full_panel_partner_manifest.tsv"
SUMMARY_JSON = RESULTS_DIR / "af3_sequence_partner_manifest.summary.json"


MANUAL_COMPLEX_WHITELIST = [
    {
        "complex_group_id": "AF3_COMPLEX_KCNQ1_KCNE1_CALM_CA",
        "complex_name": "KCNQ1-KCNE1-calmodulin-calcium complex",
        "genes_covered": ["KCNQ1", "KCNE1", "CALM1", "CALM2", "CALM3"],
        "protein_entities": ["KCNQ1", "KCNE1", "CALM1"],
        "ion_entities": ["CA"],
        "rationale": "IKs channel auxiliary subunit and calmodulin/calcium context.",
        "evidence_status": "high_confidence_manual_whitelist_literature_review_needed",
        "recommended_af3_input_mode": "full_length_complex",
    },
    {
        "complex_group_id": "AF3_COMPLEX_CACNA1C_CACNB2_CALM_CA",
        "complex_name": "CACNA1C-CACNB2-calmodulin-calcium channel complex",
        "genes_covered": ["CACNA1C", "CACNB2", "CALM1", "CALM2", "CALM3"],
        "protein_entities": ["CACNA1C", "CACNB2", "CALM1"],
        "ion_entities": ["CA"],
        "rationale": "L-type calcium-channel alpha/beta subunit and calmodulin context.",
        "evidence_status": "high_confidence_manual_whitelist_literature_review_needed",
        "recommended_af3_input_mode": "domain_level_complex",
    },
    {
        "complex_group_id": "AF3_COMPLEX_SCN5A_CALM",
        "complex_name": "SCN5A/Nav1.5 calmodulin complex",
        "genes_covered": ["SCN5A", "CALM1", "CALM2", "CALM3"],
        "protein_entities": ["SCN5A", "CALM1"],
        "ion_entities": ["CA"],
        "rationale": "Nav1.5 intracellular calmodulin regulatory context.",
        "evidence_status": "high_confidence_manual_whitelist_literature_review_needed",
        "recommended_af3_input_mode": "domain_level_complex",
    },
    {
        "complex_group_id": "AF3_COMPLEX_RYR2_CALM_FKBP1B_CA",
        "complex_name": "RYR2 calmodulin/FKBP1B calcium-release complex",
        "genes_covered": ["RYR2", "CALM1", "CALM2", "CALM3"],
        "protein_entities": ["RYR2", "CALM1", "FKBP1B"],
        "ion_entities": ["CA"],
        "rationale": "Cardiac ryanodine receptor regulatory-partner context.",
        "evidence_status": "high_confidence_manual_whitelist_literature_review_needed",
        "recommended_af3_input_mode": "domain_level_complex",
    },
    {
        "complex_group_id": "AF3_COMPLEX_CASQ2_TRDN_RYR2_CA",
        "complex_name": "CASQ2-triadin-RYR2 calcium-release context",
        "genes_covered": ["CASQ2", "TRDN", "RYR2"],
        "protein_entities": ["CASQ2", "TRDN", "RYR2"],
        "ion_entities": ["CA"],
        "rationale": "Sarcoplasmic-reticulum calcium-release complex context.",
        "evidence_status": "high_confidence_manual_whitelist_literature_review_needed",
        "recommended_af3_input_mode": "domain_level_complex",
    },
    {
        "complex_group_id": "AF3_COMPLEX_TROPONIN_THIN_FILAMENT_CA",
        "complex_name": "Cardiac troponin/thin-filament calcium complex",
        "genes_covered": ["ACTC1", "TPM1", "TNNI3", "TNNT2"],
        "protein_entities": ["ACTC1", "TPM1", "TNNI3", "TNNT2", "TNNC1"],
        "ion_entities": ["CA"],
        "rationale": "Core cardiac thin-filament regulatory complex.",
        "evidence_status": "high_confidence_manual_whitelist_literature_review_needed",
        "recommended_af3_input_mode": "full_length_complex_or_fragmented_complex",
    },
    {
        "complex_group_id": "AF3_COMPLEX_DESMOSOME_CORE",
        "complex_name": "Cardiac desmosome core complex",
        "genes_covered": ["DSP", "PKP2", "DSG2", "DSC2", "JUP"],
        "protein_entities": ["DSP", "PKP2", "DSG2", "DSC2", "JUP"],
        "ion_entities": [],
        "rationale": "Arrhythmogenic cardiomyopathy desmosomal complex context.",
        "evidence_status": "high_confidence_manual_whitelist_literature_review_needed",
        "recommended_af3_input_mode": "domain_level_complex",
    },
    {
        "complex_group_id": "AF3_COMPLEX_LMNA_EMD",
        "complex_name": "Lamin A/C-emerin nuclear-envelope complex",
        "genes_covered": ["LMNA", "EMD"],
        "protein_entities": ["LMNA", "EMD"],
        "ion_entities": [],
        "rationale": "Nuclear-envelope cardiomyopathy complex context.",
        "evidence_status": "high_confidence_manual_whitelist_literature_review_needed",
        "recommended_af3_input_mode": "domain_level_complex",
    },
    {
        "complex_group_id": "AF3_COMPLEX_PLN_ATP2A2",
        "complex_name": "Phospholamban-SERCA2a calcium-handling complex",
        "genes_covered": ["PLN"],
        "protein_entities": ["PLN", "ATP2A2"],
        "ion_entities": ["CA"],
        "rationale": "Cardiac calcium reuptake regulation context.",
        "evidence_status": "high_confidence_manual_whitelist_literature_review_needed",
        "recommended_af3_input_mode": "full_length_complex_or_domain_level_complex",
    },
    {
        "complex_group_id": "AF3_COMPLEX_RAS_MAPK_GTP",
        "complex_name": "RASopathy pathway GTP/MAPK signaling context",
        "genes_covered": ["HRAS", "KRAS", "NRAS", "RAF1", "MAP2K2", "PTPN11", "RIT1", "SOS1"],
        "protein_entities": ["HRAS", "RAF1", "MAP2K2", "PTPN11"],
        "ion_entities": [],
        "ligand_entities": ["GTP"],
        "rationale": "Syndromic cardiomyopathy/RASopathy signaling context; direct AF3 complex design needs manual review.",
        "evidence_status": "manual_review_needed_before_submission",
        "recommended_af3_input_mode": "domain_level_complex",
    },
]


PARTNER_COLUMNS = [
    "target_gene",
    "target_uniprot",
    "target_protein_name",
    "target_length",
    "target_panel_status",
    "complex_group_id",
    "complex_name",
    "partner_gene_or_entity",
    "partner_type",
    "partner_uniprot_or_identifier",
    "partner_length_or_token_count",
    "interaction_source",
    "interaction_evidence_type",
    "pubmed_ids",
    "reactome_pathway_or_complex",
    "uniprot_interaction_note",
    "intact_or_string_id",
    "pdb_or_rcsb_id",
    "ligand_or_ion_id",
    "recommended_af3_input_mode",
    "domain_start",
    "domain_end",
    "reason_for_domain_trimming",
    "token_estimate",
    "manual_review_status",
    "final_include",
    "exclusion_reason",
]


def parse_fasta(path: Path) -> list[dict[str, str]]:
    entries: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    seq_parts: list[str] = []
    header_re = re.compile(r"^>([^|]+)\|([^|]+)\|(\S+)\s+(.+)$")

    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    current["canonical_sequence"] = "".join(seq_parts)
                    current["sequence_length"] = str(len(current["canonical_sequence"]))
                    entries.append(current)
                match = header_re.match(line)
                if not match:
                    raise ValueError(f"Could not parse FASTA header: {line}")
                db, accession, entry_name, description = match.groups()
                gene_match = re.search(r"\bGN=([^\s]+)", description)
                organism_match = re.search(r"\bOS=(.+?)\s+OX=", description)
                protein_name = description.split(" OS=")[0]
                current = {
                    "gene": gene_match.group(1).upper() if gene_match else "",
                    "uniprot_accession": accession,
                    "entry_name": entry_name,
                    "protein_name": protein_name,
                    "organism": organism_match.group(1) if organism_match else "",
                    "reviewed": "true" if db == "sp" else "false",
                    "isoform_note": "canonical_fasta_entry",
                    "retrieval_date": date.today().isoformat(),
                    "retrieval_source": str(path),
                    "manual_review_needed": "false",
                    "fasta_header": line[1:],
                }
                seq_parts = []
            else:
                seq_parts.append(line)
    if current is not None:
        current["canonical_sequence"] = "".join(seq_parts)
        current["sequence_length"] = str(len(current["canonical_sequence"]))
        entries.append(current)
    return entries


def load_audit(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def select_entry(gene: str, by_gene: dict[str, list[dict[str, str]]]) -> tuple[dict[str, str] | None, str]:
    candidates = [e for e in by_gene.get(gene, []) if e["reviewed"] == "true" and e["organism"] == "Homo sapiens"]
    if not candidates:
        return None, "missing_reviewed_human_uniprot_entry"
    candidates = sorted(candidates, key=lambda e: (len(e["canonical_sequence"]), e["uniprot_accession"]))
    if len(candidates) > 1:
        # Keep deterministic longest canonical protein, but force manual review.
        return candidates[-1], "multiple_reviewed_entries_selected_longest_needs_review"
    return candidates[0], "ok"


def row_for_sequence(gene: str, entry: dict[str, str] | None, status: str, panel_status: str) -> dict[str, str]:
    if entry is None:
        return {
            "gene": gene,
            "uniprot_accession": "",
            "entry_name": "",
            "protein_name": "",
            "organism": "",
            "reviewed": "",
            "canonical_sequence": "",
            "sequence_length": "",
            "isoform_note": "",
            "retrieval_date": date.today().isoformat(),
            "retrieval_source": str(UNIPROT_FASTA),
            "manual_review_needed": "true",
            "mapping_status": status,
            "panel_status": panel_status,
        }
    row = {
        k: entry[k]
        for k in [
            "gene",
            "uniprot_accession",
            "entry_name",
            "protein_name",
            "organism",
            "reviewed",
            "canonical_sequence",
            "sequence_length",
            "isoform_note",
            "retrieval_date",
            "retrieval_source",
        ]
    }
    row["manual_review_needed"] = "true" if "needs_review" in status else "false"
    row["mapping_status"] = status
    row["panel_status"] = panel_status
    return row


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    audit_rows = load_audit(AUDIT_TSV)
    fasta_entries = parse_fasta(UNIPROT_FASTA)
    by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for entry in fasta_entries:
        if entry["gene"]:
            by_gene[entry["gene"]].append(entry)

    af3_targets = [
        r for r in audit_rows
        if r["af3_inclusion_decision"].startswith("include")
    ]
    target_genes = sorted({r["gene_normalized"] for r in af3_targets})

    partner_genes = sorted({
        gene
        for complex_row in MANUAL_COMPLEX_WHITELIST
        for gene in complex_row["protein_entities"]
    })
    all_sequence_genes = sorted(set(target_genes) | set(partner_genes))

    audit_by_gene = {r["gene_normalized"]: r for r in audit_rows}
    sequence_rows: list[dict[str, str]] = []
    sequence_by_gene: dict[str, dict[str, str]] = {}
    for gene in all_sequence_genes:
        entry, status = select_entry(gene, by_gene)
        panel_status = (
            audit_by_gene.get(gene, {}).get("panel_bucket", "external_partner_not_in_variant_panel")
        )
        seq_row = row_for_sequence(gene, entry, status, panel_status)
        sequence_rows.append(seq_row)
        sequence_by_gene[gene] = seq_row

    seq_fields = [
        "gene",
        "uniprot_accession",
        "entry_name",
        "protein_name",
        "organism",
        "reviewed",
        "canonical_sequence",
        "sequence_length",
        "isoform_note",
        "retrieval_date",
        "retrieval_source",
        "manual_review_needed",
        "mapping_status",
        "panel_status",
    ]
    write_tsv(SEQUENCE_TSV, sequence_rows, seq_fields)

    with ENTITIES_FASTA.open("w") as handle:
        for row in sequence_rows:
            if not row["canonical_sequence"]:
                continue
            handle.write(
                f">{row['gene']}|{row['uniprot_accession']}|{row['entry_name']} "
                f"{row['protein_name']}\n"
            )
            seq = row["canonical_sequence"]
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")

    whitelist_rows = []
    for item in MANUAL_COMPLEX_WHITELIST:
        row = {k: "" for k in [
            "complex_group_id", "complex_name", "genes_covered", "protein_entities",
            "ion_entities", "ligand_entities", "rationale", "evidence_status",
            "recommended_af3_input_mode",
        ]}
        row.update(item)
        row["genes_covered"] = ",".join(item.get("genes_covered", []))
        row["protein_entities"] = ",".join(item.get("protein_entities", []))
        row["ion_entities"] = ",".join(item.get("ion_entities", []))
        row["ligand_entities"] = ",".join(item.get("ligand_entities", []))
        whitelist_rows.append(row)
    write_tsv(
        WHITELIST_TSV,
        whitelist_rows,
        [
            "complex_group_id",
            "complex_name",
            "genes_covered",
            "protein_entities",
            "ion_entities",
            "ligand_entities",
            "rationale",
            "evidence_status",
            "recommended_af3_input_mode",
        ],
    )

    partner_rows: list[dict[str, str]] = []
    covered_targets = set()
    for complex_row in MANUAL_COMPLEX_WHITELIST:
        genes_covered = set(complex_row["genes_covered"]) & set(target_genes)
        for target_gene in sorted(genes_covered):
            covered_targets.add(target_gene)
            target_seq = sequence_by_gene[target_gene]
            for partner in complex_row.get("protein_entities", []):
                if partner == target_gene:
                    continue
                partner_seq = sequence_by_gene.get(partner)
                partner_rows.append({
                    "target_gene": target_gene,
                    "target_uniprot": target_seq["uniprot_accession"],
                    "target_protein_name": target_seq["protein_name"],
                    "target_length": target_seq["sequence_length"],
                    "target_panel_status": audit_by_gene.get(target_gene, {}).get("panel_bucket", ""),
                    "complex_group_id": complex_row["complex_group_id"],
                    "complex_name": complex_row["complex_name"],
                    "partner_gene_or_entity": partner,
                    "partner_type": "protein",
                    "partner_uniprot_or_identifier": partner_seq["uniprot_accession"] if partner_seq else "",
                    "partner_length_or_token_count": partner_seq["sequence_length"] if partner_seq else "",
                    "interaction_source": "manual_cardiac_complex_whitelist",
                    "interaction_evidence_type": complex_row["evidence_status"],
                    "pubmed_ids": "",
                    "reactome_pathway_or_complex": "",
                    "uniprot_interaction_note": "",
                    "intact_or_string_id": "",
                    "pdb_or_rcsb_id": "",
                    "ligand_or_ion_id": "",
                    "recommended_af3_input_mode": complex_row["recommended_af3_input_mode"],
                    "domain_start": "",
                    "domain_end": "",
                    "reason_for_domain_trimming": (
                        "requires_manual_domain_selection_before_af3"
                        if "domain" in complex_row["recommended_af3_input_mode"] else ""
                    ),
                    "token_estimate": "",
                    "manual_review_status": "needs_literature_api_confirmation_before_submission",
                    "final_include": "review_pending",
                    "exclusion_reason": "",
                })
            for ion in complex_row.get("ion_entities", []):
                partner_rows.append({
                    "target_gene": target_gene,
                    "target_uniprot": target_seq["uniprot_accession"],
                    "target_protein_name": target_seq["protein_name"],
                    "target_length": target_seq["sequence_length"],
                    "target_panel_status": audit_by_gene.get(target_gene, {}).get("panel_bucket", ""),
                    "complex_group_id": complex_row["complex_group_id"],
                    "complex_name": complex_row["complex_name"],
                    "partner_gene_or_entity": ion,
                    "partner_type": "ion",
                    "partner_uniprot_or_identifier": ion,
                    "partner_length_or_token_count": "1",
                    "interaction_source": "manual_cardiac_complex_whitelist",
                    "interaction_evidence_type": complex_row["evidence_status"],
                    "pubmed_ids": "",
                    "reactome_pathway_or_complex": "",
                    "uniprot_interaction_note": "",
                    "intact_or_string_id": "",
                    "pdb_or_rcsb_id": "",
                    "ligand_or_ion_id": ion,
                    "recommended_af3_input_mode": complex_row["recommended_af3_input_mode"],
                    "domain_start": "",
                    "domain_end": "",
                    "reason_for_domain_trimming": "",
                    "token_estimate": "",
                    "manual_review_status": "needs_literature_api_confirmation_before_submission",
                    "final_include": "review_pending",
                    "exclusion_reason": "",
                })
            for ligand in complex_row.get("ligand_entities", []):
                partner_rows.append({
                    "target_gene": target_gene,
                    "target_uniprot": target_seq["uniprot_accession"],
                    "target_protein_name": target_seq["protein_name"],
                    "target_length": target_seq["sequence_length"],
                    "target_panel_status": audit_by_gene.get(target_gene, {}).get("panel_bucket", ""),
                    "complex_group_id": complex_row["complex_group_id"],
                    "complex_name": complex_row["complex_name"],
                    "partner_gene_or_entity": ligand,
                    "partner_type": "ligand",
                    "partner_uniprot_or_identifier": ligand,
                    "partner_length_or_token_count": "1",
                    "interaction_source": "manual_cardiac_complex_whitelist",
                    "interaction_evidence_type": complex_row["evidence_status"],
                    "pubmed_ids": "",
                    "reactome_pathway_or_complex": "",
                    "uniprot_interaction_note": "",
                    "intact_or_string_id": "",
                    "pdb_or_rcsb_id": "",
                    "ligand_or_ion_id": ligand,
                    "recommended_af3_input_mode": complex_row["recommended_af3_input_mode"],
                    "domain_start": "",
                    "domain_end": "",
                    "reason_for_domain_trimming": "",
                    "token_estimate": "",
                    "manual_review_status": "needs_literature_api_confirmation_before_submission",
                    "final_include": "review_pending",
                    "exclusion_reason": "",
                })

    for target_gene in target_genes:
        if target_gene in covered_targets:
            continue
        target_seq = sequence_by_gene[target_gene]
        partner_rows.append({
            "target_gene": target_gene,
            "target_uniprot": target_seq["uniprot_accession"],
            "target_protein_name": target_seq["protein_name"],
            "target_length": target_seq["sequence_length"],
            "target_panel_status": audit_by_gene.get(target_gene, {}).get("panel_bucket", ""),
            "complex_group_id": "",
            "complex_name": "",
            "partner_gene_or_entity": "",
            "partner_type": "none_found",
            "partner_uniprot_or_identifier": "",
            "partner_length_or_token_count": "",
            "interaction_source": "none_yet",
            "interaction_evidence_type": "manual_review_needed",
            "pubmed_ids": "",
            "reactome_pathway_or_complex": "",
            "uniprot_interaction_note": "",
            "intact_or_string_id": "",
            "pdb_or_rcsb_id": "",
            "ligand_or_ion_id": "",
            "recommended_af3_input_mode": "",
            "domain_start": "",
            "domain_end": "",
            "reason_for_domain_trimming": "",
            "token_estimate": "",
            "manual_review_status": "manual_review_needed",
            "final_include": "no",
            "exclusion_reason": "no_curated_partner_assigned_yet",
        })

    write_tsv(PARTNER_TSV, partner_rows, PARTNER_COLUMNS)

    summary = {
        "date": date.today().isoformat(),
        "audit_tsv": str(AUDIT_TSV),
        "uniprot_fasta": str(UNIPROT_FASTA),
        "af3_target_gene_count": len(target_genes),
        "af3_target_genes": target_genes,
        "sequence_gene_count": len(sequence_rows),
        "missing_sequence_genes": [
            r["gene"] for r in sequence_rows if not r["uniprot_accession"]
        ],
        "manual_review_sequence_genes": [
            r["gene"] for r in sequence_rows if r["manual_review_needed"] == "true"
        ],
        "manual_whitelist_complex_count": len(MANUAL_COMPLEX_WHITELIST),
        "partner_manifest_rows": len(partner_rows),
        "targets_with_manual_whitelist_context": sorted(covered_targets),
        "targets_manual_review_needed_no_partner_yet": sorted(set(target_genes) - covered_targets),
        "outputs": {
            "sequence_tsv": str(SEQUENCE_TSV),
            "all_entities_fasta": str(ENTITIES_FASTA),
            "manual_whitelist_tsv": str(WHITELIST_TSV),
            "partner_manifest_tsv": str(PARTNER_TSV),
        },
        "safety_note": "No existing model matrix or previous dataset file is modified by this script.",
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
