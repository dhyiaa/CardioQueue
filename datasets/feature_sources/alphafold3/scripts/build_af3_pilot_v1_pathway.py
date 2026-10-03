#!/usr/bin/env python3
"""Build the AF3 pilot-v1 pathway for packet-ready and near-ready genes.

This script is intentionally additive. It creates AF3 pilot manifests, review
packets, output templates, and feature schemas without touching the existing
modeling matrices.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
GENE_REVIEW_TSV = ROOT / "results/af3/af3_gene_level_coverage_review_v2.tsv"
FIRST_WAVE_JOB_MANIFEST = ROOT / "datasets/feature_sources/alphafold3/jobs/af3_job_manifest.tsv"
SECOND_WAVE_TSV = ROOT / "results/af3/af3_second_wave_curated_partner_decisions.tsv"
BASE_SEQ_TSV = ROOT / "datasets/feature_sources/alphafold3/interim/uniprot_gene_to_sequence.tsv"
SECOND_SEQ_TSV = ROOT / "datasets/feature_sources/alphafold3/interim/af3_second_wave_entities_to_sequence.tsv"

PILOT_RESULTS_DIR = ROOT / "results/af3/pilot_v1"
PILOT_JOBS_DIR = ROOT / "datasets/feature_sources/alphafold3/jobs/pilot_v1"
PILOT_PACKET_DIR = PILOT_JOBS_DIR / "manual_submission_packets"
PILOT_RAW_RETURN_DIR = ROOT / "datasets/feature_sources/alphafold3/raw/pilot_v1_server_outputs"
PILOT_PROCESSED_DIR = ROOT / "datasets/feature_sources/alphafold3/processed/pilot_v1"

PILOT_GENE_SET_TSV = PILOT_RESULTS_DIR / "af3_pilot_v1_gene_set.tsv"
PILOT_COMPLEX_PLAN_TSV = PILOT_RESULTS_DIR / "af3_pilot_v1_complex_plan.tsv"
PILOT_SUBMISSION_MANIFEST_TSV = PILOT_JOBS_DIR / "af3_pilot_v1_submission_manifest.tsv"
PILOT_OUTPUT_TEMPLATE_TSV = PILOT_JOBS_DIR / "af3_pilot_v1_completed_job_manifest_template.tsv"
FEATURE_SCHEMA_TSV = PILOT_PROCESSED_DIR / "af3_pilot_v1_feature_schema.tsv"
SUMMARY_JSON = PILOT_RESULTS_DIR / "af3_pilot_v1_pathway.summary.json"


PILOT_STATUSES = {
    "packet_draft_ready",
    "second_wave_packet_candidate_needs_final_evidence_check",
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


def split_values(value: str, sep: str = ",") -> list[str]:
    return [v.strip() for v in value.split(sep) if v.strip()]


def uniq_join(values: list[str], sep: str = ";") -> str:
    return sep.join(sorted({v for v in values if v}))


def wrap_sequence(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def load_sequences() -> dict[str, dict[str, str]]:
    seq_by_gene: dict[str, dict[str, str]] = {}
    for path in [BASE_SEQ_TSV, SECOND_SEQ_TSV]:
        for row in read_tsv(path):
            seq_by_gene.setdefault(row["gene"], row)
    return seq_by_gene


def first_wave_complexes(
    genes: list[dict[str, str]],
    job_manifest: list[dict[str, str]],
) -> list[dict[str, str]]:
    rows = []
    genes_by_job: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in genes:
        if row["partner_assignment_status_v2"] != "packet_draft_ready":
            continue
        for job_id in split_values(row["draft_job_ids"], ";"):
            genes_by_job[job_id].append(row)

    job_by_id = {row["job_id"]: row for row in job_manifest}
    for job_id, gene_rows in sorted(genes_by_job.items()):
        job = job_by_id.get(job_id, {})
        cid = job.get("complex_group_id", "")
        if not job or not cid:
            continue
        rows.append(
            {
                "pilot_stage": "already_draft_packet_ready",
                "complex_group_id": cid,
                "complex_name": job.get("job_name", cid),
                "target_genes": ",".join(sorted({r["gene"] for r in gene_rows})),
                "protein_entities": job.get("protein_entities", ""),
                "nonprotein_entities": ",".join([v for v in [job.get("ligands", ""), job.get("ions", "")] if v]),
                "estimated_tokens": job.get("total_tokens_estimated", ""),
                "submission_status": "draft_packet_exists_needs_user_review",
                "job_id": job_id,
                "packet_dir": job.get("packet_dir", ""),
                "evidence_status": "structured_evidence_ready_first_wave",
                "notes": "Existing draft packet from first-wave AF3 design review.",
            }
        )
    return rows


def second_wave_complexes(second_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    selected = [r for r in second_rows if r["manual_review_status"] == "ready_for_evidence_check"]
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in selected:
        grouped[row["complex_group_id"]].append(row)

    rows = []
    for cid, group in sorted(grouped.items()):
        first = group[0]
        rows.append(
            {
                "pilot_stage": "second_wave_packet_candidate",
                "complex_group_id": cid,
                "complex_name": first["complex_name"],
                "target_genes": ",".join(sorted({r["target_gene"] for r in group})),
                "protein_entities": first["protein_entities"],
                "nonprotein_entities": first["nonprotein_entities"],
                "estimated_tokens": first["estimated_tokens"],
                "submission_status": "review_packet_created_not_submission_ready",
                "job_id": f"AF3_PILOTV1_REVIEW_{cid.replace('AF3_COMPLEX_', '')}",
                "packet_dir": str(PILOT_PACKET_DIR / f"AF3_PILOTV1_REVIEW_{cid.replace('AF3_COMPLEX_', '')}"),
                "evidence_status": "needs_final_user_evidence_check",
                "notes": "Second-wave candidate. Use packet to review citations/evidence before AF3 submission.",
            }
        )
    return rows


def make_review_packet(row: dict[str, str], seq_by_gene: dict[str, dict[str, str]], evidence_rows: list[dict[str, str]]) -> None:
    out_dir = Path(row["packet_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    proteins = split_values(row["protein_entities"])
    nonprotein = split_values(row["nonprotein_entities"])
    seen_counts: dict[str, int] = defaultdict(int)
    fasta_parts = []
    chain_rows = []
    for protein in proteins:
        seen_counts[protein] += 1
        seq = seq_by_gene[protein]
        copy_suffix = f"_copy{seen_counts[protein]}" if proteins.count(protein) > 1 else ""
        fasta_parts.append(
            f">{protein}{copy_suffix}|{seq['uniprot_accession']}|{seq['entry_name']}|{seq['protein_name']}\n"
            f"{wrap_sequence(seq['canonical_sequence'])}"
        )
        chain_rows.append(
            {
                "entity_order": str(len(chain_rows) + 1),
                "entity_type": "protein",
                "entity_name": protein,
                "copy_number": str(seen_counts[protein]),
                "uniprot_accession_or_id": seq["uniprot_accession"],
                "sequence_range": f"1-{seq['sequence_length']}",
                "notes": "full_length_canonical_reviewed_uniprot_sequence",
            }
        )
    for entity in nonprotein:
        chain_rows.append(
            {
                "entity_order": str(len(chain_rows) + 1),
                "entity_type": "ligand_or_ion",
                "entity_name": entity,
                "copy_number": "1",
                "uniprot_accession_or_id": entity,
                "sequence_range": "",
                "notes": "Add manually in AlphaFold Server UI only after CCD/entity review.",
            }
        )

    (out_dir / "job_name.txt").write_text(row["job_id"].replace("AF3_PILOTV1_REVIEW_", "AF3_") + "\n")
    (out_dir / "entities.fasta").write_text("\n\n".join(fasta_parts) + "\n")
    write_tsv(
        out_dir / "expected_entity_map.tsv",
        chain_rows,
        ["entity_order", "entity_type", "entity_name", "copy_number", "uniprot_accession_or_id", "sequence_range", "notes"],
    )

    relevant_evidence = [r for r in evidence_rows if r["complex_group_id"] == row["complex_group_id"]]
    evidence_lines = [
        f"complex_group_id: {row['complex_group_id']}",
        f"complex_name: {row['complex_name']}",
        f"target_genes: {row['target_genes']}",
        f"protein_entities: {row['protein_entities']}",
        f"nonprotein_entities: {row['nonprotein_entities'] or 'none'}",
        f"estimated_tokens: {row['estimated_tokens']}",
        "",
        "User/manual checks before submission:",
        "- Confirm these entities form a plausible same-complex or pairwise biological context.",
        "- Confirm PubMed/Reactome/PDB/IntAct support is not only pathway-level noise.",
        "- Confirm full-length modeling is biologically acceptable for this complex.",
        "- Confirm duplicate chains or homomer stoichiometry in the AlphaFold Server UI.",
        "",
        "Evidence rows:",
    ]
    for ev in relevant_evidence:
        evidence_lines.append(
            f"- {ev['target_gene']}: {ev['rationale']} "
            f"Evidence anchor={ev['evidence_anchor']}; "
            f"UniProt PubMed={ev['pubmed_ids_from_uniprot_gene_records']}; "
            f"IntAct PubMed={ev['pubmed_ids_from_intact_gene_records']}; "
            f"Reactome={ev['reactome_ids_from_gene_records']}; "
            f"PDB={ev['pdb_ids_from_gene_records']}"
        )
    (out_dir / "evidence_review_checklist.txt").write_text("\n".join(evidence_lines) + "\n")
    (out_dir / "notes_for_server_ui.txt").write_text(
        "This is a second-wave review packet, not an approved submission packet.\n"
        "Submit only after the evidence_review_checklist.txt items are manually approved.\n"
        "Use the reviewed UniProt sequences in entities.fasta.\n"
        "After AF3 output is downloaded, place the zip or extracted folder in raw/pilot_v1_server_outputs/ and record it in af3_pilot_v1_completed_job_manifest_template.tsv.\n"
    )


def build_feature_schema() -> list[dict[str, str]]:
    return [
        {"feature_name": "af3_pilot_v1_gene_covered", "level": "variant", "type": "binary", "description": "Variant gene belongs to the AF3 pilot-v1 target set.", "missing_value": "0"},
        {"feature_name": "af3_pilot_v1_complex_count_covering_gene", "level": "variant", "type": "integer", "description": "Number of pilot-v1 AF3 complex jobs covering the variant gene.", "missing_value": "0"},
        {"feature_name": "af3_pilot_v1_complex_ids", "level": "variant", "type": "categorical_list", "description": "AF3 complex IDs that can annotate this gene.", "missing_value": ""},
        {"feature_name": "af3_residue_mapped_to_model", "level": "variant", "type": "binary", "description": "Variant residue maps to a modeled AF3 chain after server output parsing.", "missing_value": "0"},
        {"feature_name": "af3_local_plddt", "level": "variant", "type": "numeric", "description": "AF3 local residue confidence for the variant residue.", "missing_value": "NA"},
        {"feature_name": "af3_min_distance_to_any_partner_A", "level": "variant", "type": "numeric", "description": "Minimum heavy-atom distance in Angstroms from variant residue to any partner protein chain.", "missing_value": "NA"},
        {"feature_name": "af3_any_partner_interface_5A", "level": "variant", "type": "binary", "description": "Variant residue is within 5 Angstroms of any partner protein chain.", "missing_value": "0"},
        {"feature_name": "af3_any_partner_interface_8A", "level": "variant", "type": "binary", "description": "Variant residue is within 8 Angstroms of any partner protein chain.", "missing_value": "0"},
        {"feature_name": "af3_best_partner_entity", "level": "variant", "type": "categorical", "description": "Closest partner entity to the variant residue.", "missing_value": ""},
        {"feature_name": "af3_min_distance_to_ion_or_ligand_A", "level": "variant", "type": "numeric", "description": "Minimum distance from variant residue to modeled ion or ligand, if present.", "missing_value": "NA"},
        {"feature_name": "af3_near_ion_or_ligand_5A", "level": "variant", "type": "binary", "description": "Variant residue lies within 5 Angstroms of modeled ion/ligand.", "missing_value": "0"},
        {"feature_name": "af3_local_interface_pae", "level": "variant", "type": "numeric", "description": "Local interface PAE involving the variant chain and closest partner, if parseable from AF3 output.", "missing_value": "NA"},
        {"feature_name": "af3_interface_confidence_class", "level": "variant", "type": "categorical", "description": "Derived confidence class using pLDDT and interface PAE.", "missing_value": "missing"},
        {"feature_name": "af3_feature_missing_reason", "level": "variant", "type": "categorical", "description": "Reason AF3 feature values are unavailable.", "missing_value": "not_pilot_gene_or_output_missing"},
    ]


def main() -> None:
    PILOT_PACKET_DIR.mkdir(parents=True, exist_ok=True)
    PILOT_RAW_RETURN_DIR.mkdir(parents=True, exist_ok=True)
    PILOT_PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    gene_rows = read_tsv(GENE_REVIEW_TSV)
    job_manifest_rows = read_tsv(FIRST_WAVE_JOB_MANIFEST)
    second_rows = read_tsv(SECOND_WAVE_TSV)
    seq_by_gene = load_sequences()
    job_by_id = {r["job_id"]: r for r in job_manifest_rows}

    pilot_gene_rows = []
    for row in gene_rows:
        status = row["partner_assignment_status_v2"]
        if status not in PILOT_STATUSES:
            continue
        if status == "packet_draft_ready":
            ready_jobs = [job_by_id[j] for j in split_values(row["draft_job_ids"], ";") if j in job_by_id]
            complex_ids = ";".join(j["complex_group_id"] for j in ready_jobs)
            complex_names = ";".join(j["job_name"] for j in ready_jobs)
            protein_entities = ";".join(j["protein_entities"] for j in ready_jobs)
            nonprotein_entities = ";".join(",".join([v for v in [j["ligands"], j["ions"]] if v]) for j in ready_jobs)
        else:
            complex_ids = row["second_wave_complex_group_id"]
            complex_names = row["second_wave_complex_name"]
            protein_entities = row["second_wave_protein_entities"]
            nonprotein_entities = row["second_wave_nonprotein_entities"]
        pilot_gene_rows.append(
            {
                "gene": row["gene"],
                "pilot_stage": "already_draft_packet_ready" if status == "packet_draft_ready" else "second_wave_packet_candidate",
                "partner_assignment_status_v2": status,
                "panel_bucket": row["panel_bucket"],
                "primary_model_inclusion_values": row["primary_model_inclusion_values"],
                "n_total_variants": row["n_total_variants"],
                "n_missense": row["n_missense"],
                "uniprot_accession": row["uniprot_accession"],
                "sequence_length": row["sequence_length"],
                "complex_group_ids": complex_ids,
                "complex_names": complex_names,
                "protein_entities": protein_entities,
                "nonprotein_entities": nonprotein_entities,
                "evidence_summary": row["combined_support_levels"] or row["second_wave_evidence_anchor"],
                "fix_needed": row["fix_needed_v2"],
            }
        )

    first_complex_rows = first_wave_complexes(gene_rows, job_manifest_rows)
    second_complex_rows = second_wave_complexes(second_rows)
    complex_rows = first_complex_rows + second_complex_rows

    second_evidence_rows = [r for r in second_rows if r["manual_review_status"] == "ready_for_evidence_check"]
    for row in second_complex_rows:
        make_review_packet(row, seq_by_gene, second_evidence_rows)

    submission_rows = []
    for order, row in enumerate(complex_rows, start=1):
        submission_rows.append(
            {
                "submission_order": str(order),
                "job_id": row["job_id"],
                "complex_group_id": row["complex_group_id"],
                "complex_name": row["complex_name"],
                "pilot_stage": row["pilot_stage"],
                "target_genes": row["target_genes"],
                "protein_entities": row["protein_entities"],
                "nonprotein_entities": row["nonprotein_entities"],
                "estimated_tokens": row["estimated_tokens"],
                "submission_status": row["submission_status"],
                "packet_dir": row["packet_dir"],
                "af3_status": "not_submitted",
                "af3_server_job_name": "",
                "af3_output_path": "",
                "notes": row["notes"],
            }
        )

    completed_template_rows = [
        {
            "job_id": r["job_id"],
            "complex_group_id": r["complex_group_id"],
            "af3_status": "not_submitted",
            "af3_server_job_name": "",
            "downloaded_zip_or_folder_path": "",
            "extracted_output_dir": "",
            "ranking_json_or_summary_path": "",
            "model_cif_or_pdb_path": "",
            "confidence_json_path": "",
            "chain_mapping_confirmed": "false",
            "parse_ready": "false",
            "notes": "Fill this row after AF3 server output is downloaded.",
        }
        for r in submission_rows
    ]

    gene_fields = [
        "gene",
        "pilot_stage",
        "partner_assignment_status_v2",
        "panel_bucket",
        "primary_model_inclusion_values",
        "n_total_variants",
        "n_missense",
        "uniprot_accession",
        "sequence_length",
        "complex_group_ids",
        "complex_names",
        "protein_entities",
        "nonprotein_entities",
        "evidence_summary",
        "fix_needed",
    ]
    complex_fields = [
        "pilot_stage",
        "complex_group_id",
        "complex_name",
        "target_genes",
        "protein_entities",
        "nonprotein_entities",
        "estimated_tokens",
        "submission_status",
        "job_id",
        "packet_dir",
        "evidence_status",
        "notes",
    ]
    submission_fields = [
        "submission_order",
        "job_id",
        "complex_group_id",
        "complex_name",
        "pilot_stage",
        "target_genes",
        "protein_entities",
        "nonprotein_entities",
        "estimated_tokens",
        "submission_status",
        "packet_dir",
        "af3_status",
        "af3_server_job_name",
        "af3_output_path",
        "notes",
    ]
    output_fields = [
        "job_id",
        "complex_group_id",
        "af3_status",
        "af3_server_job_name",
        "downloaded_zip_or_folder_path",
        "extracted_output_dir",
        "ranking_json_or_summary_path",
        "model_cif_or_pdb_path",
        "confidence_json_path",
        "chain_mapping_confirmed",
        "parse_ready",
        "notes",
    ]
    write_tsv(PILOT_GENE_SET_TSV, pilot_gene_rows, gene_fields)
    write_tsv(PILOT_COMPLEX_PLAN_TSV, complex_rows, complex_fields)
    write_tsv(PILOT_SUBMISSION_MANIFEST_TSV, submission_rows, submission_fields)
    write_tsv(PILOT_OUTPUT_TEMPLATE_TSV, completed_template_rows, output_fields)
    write_tsv(FEATURE_SCHEMA_TSV, build_feature_schema(), ["feature_name", "level", "type", "description", "missing_value"])

    summary = {
        "date": date.today().isoformat(),
        "pilot_name": "af3_pilot_v1",
        "target_gene_count": len(pilot_gene_rows),
        "status_counts": dict(Counter(r["partner_assignment_status_v2"] for r in pilot_gene_rows)),
        "complex_count": len(complex_rows),
        "first_wave_complex_count": len(first_complex_rows),
        "second_wave_review_complex_count": len(second_complex_rows),
        "target_genes": [r["gene"] for r in sorted(pilot_gene_rows, key=lambda x: x["gene"])],
        "outputs": {
            "gene_set": str(PILOT_GENE_SET_TSV),
            "complex_plan": str(PILOT_COMPLEX_PLAN_TSV),
            "submission_manifest": str(PILOT_SUBMISSION_MANIFEST_TSV),
            "completed_job_manifest_template": str(PILOT_OUTPUT_TEMPLATE_TSV),
            "feature_schema": str(FEATURE_SCHEMA_TSV),
            "review_packet_dir": str(PILOT_PACKET_DIR),
            "expected_server_output_dir": str(PILOT_RAW_RETURN_DIR),
        },
        "safety_note": "Pilot-v1 files are additive only. Existing modeling matrices were not modified.",
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
