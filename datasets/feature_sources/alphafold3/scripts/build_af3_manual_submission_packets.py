#!/usr/bin/env python3
"""Create draft manual AlphaFold 3 submission packets for ready complexes."""

from __future__ import annotations

import csv
import json
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SEQ_TSV = ROOT / "datasets/feature_sources/alphafold3/interim/uniprot_gene_to_sequence.tsv"
DESIGN_TSV = ROOT / "results/af3/af3_complex_design_review.tsv"
PAIR_EVIDENCE_TSV = ROOT / "results/af3/af3_partner_combined_evidence_review_queue.tsv"
PACKET_DIR = ROOT / "datasets/feature_sources/alphafold3/jobs/manual_submission_packets"
JOB_MANIFEST_TSV = ROOT / "datasets/feature_sources/alphafold3/jobs/af3_job_manifest.tsv"
SUMMARY_JSON = ROOT / "results/af3/af3_manual_submission_packets.summary.json"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def split_csv(value: str) -> list[str]:
    return [v.strip() for v in value.split(",") if v.strip()]


def wrap_sequence(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main() -> None:
    seq_rows = read_tsv(SEQ_TSV)
    design_rows = read_tsv(DESIGN_TSV)
    evidence_rows = read_tsv(PAIR_EVIDENCE_TSV)

    seq_by_gene = {r["gene"]: r for r in seq_rows}
    evidence_by_complex = {}
    for row in evidence_rows:
        cid = row["complex_group_id"]
        if not cid:
            continue
        evidence_by_complex.setdefault(cid, []).append(row)

    ready_rows = [r for r in design_rows if r["packet_draft_status"] == "ready_for_manual_packet_draft"]
    manifest_rows = []
    PACKET_DIR.mkdir(parents=True, exist_ok=True)

    for idx, row in enumerate(ready_rows, start=1):
        cid = row["complex_group_id"]
        job_id = f"AF3_JOB_{idx:03d}_{cid.replace('AF3_COMPLEX_', '')}"
        job_name = cid.replace("AF3_COMPLEX_", "AF3_")
        out_dir = PACKET_DIR / job_id
        out_dir.mkdir(parents=True, exist_ok=True)

        proteins = split_csv(row["protein_entities"])
        ions = split_csv(row["ion_entities"])
        ligands = split_csv(row["ligand_entities"])

        fasta_parts = []
        chain_rows = []
        chain_ord = ord("A")
        for protein in proteins:
            seq = seq_by_gene[protein]["canonical_sequence"]
            accession = seq_by_gene[protein]["uniprot_accession"]
            entry_name = seq_by_gene[protein]["entry_name"]
            fasta_parts.append(
                f">{protein}|{accession}|{entry_name}|{seq_by_gene[protein]['protein_name']}\n"
                f"{wrap_sequence(seq)}"
            )
            chain_rows.append({
                "expected_chain_or_entity": chr(chain_ord),
                "entity_type": "protein",
                "entity_name": protein,
                "uniprot_accession_or_id": accession,
                "sequence_range": f"1-{len(seq)}",
                "notes": "full_length_canonical_uniprot_sequence",
            })
            chain_ord += 1
        for ion in ions:
            chain_rows.append({
                "expected_chain_or_entity": chr(chain_ord),
                "entity_type": "ion",
                "entity_name": ion,
                "uniprot_accession_or_id": ion,
                "sequence_range": "",
                "notes": "Add in AlphaFold Server ligand/ion UI; likely CCD code CA for calcium ion.",
            })
            chain_ord += 1
        for ligand in ligands:
            chain_rows.append({
                "expected_chain_or_entity": chr(chain_ord),
                "entity_type": "ligand",
                "entity_name": ligand,
                "uniprot_accession_or_id": ligand,
                "sequence_range": "",
                "notes": "Add in AlphaFold Server ligand UI; verify CCD code manually.",
            })
            chain_ord += 1

        (out_dir / "job_name.txt").write_text(job_name + "\n")
        (out_dir / "entities.fasta").write_text("\n\n".join(fasta_parts) + "\n")
        (out_dir / "ligands_and_ions.txt").write_text(
            "ions: " + (",".join(ions) if ions else "none") + "\n"
            "ligands: " + (",".join(ligands) if ligands else "none") + "\n"
            "manual_note: enter non-protein entities through the AlphaFold Server UI and verify supported CCD/entity names.\n"
        )
        write_tsv(
            out_dir / "expected_chain_map.tsv",
            chain_rows,
            ["expected_chain_or_entity", "entity_type", "entity_name", "uniprot_accession_or_id", "sequence_range", "notes"],
        )

        evidence_lines = [
            f"complex_group_id: {cid}",
            f"complex_name: {row['complex_name']}",
            f"combined_readiness_status: {row['combined_readiness_status']}",
            f"estimated_af3_tokens: {row['estimated_af3_tokens']}",
            "",
            "Evidence rows:",
        ]
        for ev in evidence_by_complex.get(cid, []):
            evidence_lines.append(
                f"- {ev['target_gene']} -> {ev['partner_gene_or_entity']} "
                f"({ev['partner_type']}): {ev['combined_support_level']}; "
                f"PubMed={ev['uniprot_pair_pubmed_ids'] or ev['intact_psicquic_pair_pubmed_ids']}; "
                f"Reactome={ev['reactome_shared_pathway_ids']}; PDB={ev['rcsb_shared_pdb_ids']}"
            )
        (out_dir / "evidence_summary.txt").write_text("\n".join(evidence_lines) + "\n")

        notes = [
            "This is a draft manual AlphaFold 3 server submission packet.",
            "Do not treat it as submitted or final until manually reviewed.",
            "Use canonical reviewed human UniProt protein sequences from entities.fasta.",
            "After AF3 returns output, verify actual chain mapping from the mmCIF/job JSON.",
            "CALM note: CALM1/CALM2/CALM3 canonical protein sequences are identical in the local UniProt FASTA; CALM1 is used as the sequence entity while variants remain mapped back to their original CALM gene.",
            f"Design note: {row['manual_domain_review_note'] or 'Full-length token-feasible draft.'}",
        ]
        (out_dir / "notes_for_server_ui.txt").write_text("\n".join(notes) + "\n")

        manifest_rows.append({
            "job_id": job_id,
            "complex_group_id": cid,
            "job_name": job_name,
            "genes_covered": row["genes_covered"],
            "protein_entities": row["protein_entities"],
            "ligands": row["ligand_entities"],
            "ions": row["ion_entities"],
            "total_tokens_estimated": row["estimated_af3_tokens"],
            "full_length_or_domain": "full_length",
            "reason_for_job": row["combined_readiness_status"],
            "evidence_summary": str(out_dir / "evidence_summary.txt"),
            "manual_ready_for_submission": "draft_needs_user_review",
            "submission_day": "",
            "submission_order": str(idx),
            "af3_status": "not_submitted",
            "af3_output_zip_path": "",
            "packet_dir": str(out_dir),
        })

    manifest_fields = [
        "job_id", "complex_group_id", "job_name", "genes_covered", "protein_entities",
        "ligands", "ions", "total_tokens_estimated", "full_length_or_domain",
        "reason_for_job", "evidence_summary", "manual_ready_for_submission",
        "submission_day", "submission_order", "af3_status", "af3_output_zip_path",
        "packet_dir",
    ]
    write_tsv(JOB_MANIFEST_TSV, manifest_rows, manifest_fields)

    summary = {
        "date": date.today().isoformat(),
        "packet_count": len(manifest_rows),
        "job_manifest": str(JOB_MANIFEST_TSV),
        "packet_dir": str(PACKET_DIR),
        "packet_job_ids": [r["job_id"] for r in manifest_rows],
        "safety_note": "Draft manual packets only; nothing was submitted to AlphaFold Server.",
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
