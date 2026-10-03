#!/usr/bin/env python3
"""Build AF3 complex design review table before submission packets."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SEQ_TSV = ROOT / "datasets/feature_sources/alphafold3/interim/uniprot_gene_to_sequence.tsv"
WHITELIST_TSV = ROOT / "datasets/feature_sources/alphafold3/interim/af3_manual_cardiac_complex_whitelist.tsv"
READINESS_TSV = ROOT / "results/af3/af3_complex_combined_evidence_readiness.tsv"
MATRIX_TSV = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
OUT_TSV = ROOT / "results/af3/af3_complex_design_review.tsv"
GROUPS_TSV = ROOT / "results/af3/af3_complex_groups.tsv"
SUMMARY_JSON = ROOT / "results/af3/af3_complex_design_review.summary.json"

FULL_LENGTH_SOFT_LIMIT = 4500
FULL_LENGTH_HARD_LIMIT = 5000

DOMAIN_REVIEW_NOTES = {
    "AF3_COMPLEX_CASQ2_TRDN_RYR2_CA": (
        "Full-length RYR2 makes the complex exceed the AF3 server token limit. "
        "Use full CASQ2/TRDN plus a literature-supported RYR2 luminal/cytosolic "
        "domain only after manual domain selection."
    ),
    "AF3_COMPLEX_DESMOSOME_CORE": (
        "Full desmosome core exceeds token limit. Split into smaller jobs such as "
        "DSG2/DSC2 extracellular cadherin context and DSP/PKP2/JUP cytoplasmic "
        "plaque context after manual domain selection."
    ),
    "AF3_COMPLEX_RYR2_CALM_FKBP1B_CA": (
        "Full RYR2 with CALM/FKBP1B exceeds token limit. Use a known CALM/FKBP1B "
        "RYR2 regulatory domain only after literature/domain review."
    ),
    "AF3_COMPLEX_CACNA1C_CACNB2_CALM_CA": (
        "Token count is feasible, but full-length CACNA1C is a large membrane "
        "protein. Prefer domain review around beta-subunit interaction and CALM/IQ "
        "regulatory regions before submission."
    ),
    "AF3_COMPLEX_KCNQ1_KCNE1_CALM_CA": (
        "Token count is feasible. Full-length KCNQ1/KCNE1/CALM may be a reasonable "
        "first AF3 job, but membrane-span context and CALM isoform proxy should be "
        "explicitly reviewed."
    ),
    "AF3_COMPLEX_RAS_MAPK_GTP": (
        "This is pathway context, not one guaranteed physical complex. Split into "
        "small biologically specific domain jobs such as RAS-RAF and PTPN11 domain "
        "contexts after literature review."
    ),
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_csv_cell(value: str) -> list[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def parse_position(value: str) -> int | None:
    if not value:
        return None
    try:
        numeric = float(value)
    except ValueError:
        return None
    if numeric.is_integer() and numeric > 0:
        return int(numeric)
    return None


def variant_position_summary() -> dict[str, dict[str, str]]:
    wanted = {
        "primary_gene", "protein_position", "protein_variant_type",
        "primary_model_inclusion",
    }
    summaries: dict[str, list[int]] = defaultdict(list)
    missense_counts = defaultdict(int)
    row_counts = defaultdict(int)
    with MATRIX_TSV.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = wanted - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing expected modeling columns: {sorted(missing)}")
        for row in reader:
            gene = (row["primary_gene"] or "").upper()
            if not gene:
                continue
            row_counts[gene] += 1
            pos = parse_position(row["protein_position"])
            if pos is not None:
                summaries[gene].append(pos)
            if row["protein_variant_type"] == "missense":
                missense_counts[gene] += 1
    out = {}
    for gene, vals in summaries.items():
        vals = sorted(vals)
        out[gene] = {
            "n_model_rows": str(row_counts[gene]),
            "n_missense_rows": str(missense_counts[gene]),
            "n_protein_position_rows": str(len(vals)),
            "protein_position_min": str(vals[0]) if vals else "",
            "protein_position_max": str(vals[-1]) if vals else "",
            "protein_position_p05": str(vals[int(0.05 * (len(vals) - 1))]) if vals else "",
            "protein_position_p95": str(vals[int(0.95 * (len(vals) - 1))]) if vals else "",
        }
    return out


def feasibility(total_tokens: int, readiness: str, complex_id: str) -> tuple[str, str, str]:
    if "RAS_MAPK" in complex_id:
        return (
            "split_into_specific_domain_jobs",
            "not_ready_for_packet",
            DOMAIN_REVIEW_NOTES[complex_id],
        )
    if total_tokens <= FULL_LENGTH_SOFT_LIMIT:
        if readiness == "structured_evidence_ready_for_manual_af3_design":
            return ("full_length_feasible", "ready_for_manual_packet_draft", DOMAIN_REVIEW_NOTES.get(complex_id, ""))
        return ("full_length_token_feasible_but_evidence_or_domain_review_needed", "hold_for_manual_review", DOMAIN_REVIEW_NOTES.get(complex_id, ""))
    if total_tokens <= FULL_LENGTH_HARD_LIMIT:
        return ("borderline_full_length_feasible", "hold_for_manual_review", DOMAIN_REVIEW_NOTES.get(complex_id, ""))
    return ("domain_or_split_required", "not_ready_for_packet", DOMAIN_REVIEW_NOTES.get(complex_id, "Full-length complex exceeds token limit; select domains/fragments first."))


def main() -> None:
    seq_rows = read_tsv(SEQ_TSV)
    whitelist_rows = read_tsv(WHITELIST_TSV)
    readiness_rows = read_tsv(READINESS_TSV)
    pos_summary = variant_position_summary()

    seq_by_gene = {r["gene"]: r for r in seq_rows}
    ready_by_complex = {r["complex_group_id"]: r for r in readiness_rows}

    design_rows = []
    group_rows = []
    for w in whitelist_rows:
        cid = w["complex_group_id"]
        protein_entities = parse_csv_cell(w["protein_entities"])
        ion_entities = parse_csv_cell(w.get("ion_entities", ""))
        ligand_entities = parse_csv_cell(w.get("ligand_entities", ""))
        genes_covered = parse_csv_cell(w["genes_covered"])
        missing_entities = [g for g in protein_entities if g not in seq_by_gene]
        lengths = [int(seq_by_gene[g]["sequence_length"]) for g in protein_entities if g in seq_by_gene]
        total_tokens = sum(lengths) + len(ion_entities) + len(ligand_entities)
        readiness = ready_by_complex.get(cid, {}).get("combined_readiness_status", "missing_readiness")
        design_mode, packet_status, note = feasibility(total_tokens, readiness, cid)

        per_gene_positions = []
        for gene in genes_covered:
            ps = pos_summary.get(gene, {})
            if ps:
                per_gene_positions.append(
                    f"{gene}:n_pos={ps['n_protein_position_rows']},range={ps['protein_position_min']}-{ps['protein_position_max']},p05_p95={ps['protein_position_p05']}-{ps['protein_position_p95']}"
                )

        row = {
            "complex_group_id": cid,
            "complex_name": w["complex_name"],
            "genes_covered": ",".join(genes_covered),
            "protein_entities": ",".join(protein_entities),
            "ion_entities": ",".join(ion_entities),
            "ligand_entities": ",".join(ligand_entities),
            "protein_entity_lengths": ",".join(f"{g}:{seq_by_gene[g]['sequence_length']}" for g in protein_entities if g in seq_by_gene),
            "estimated_af3_tokens": str(total_tokens),
            "full_length_soft_limit": str(FULL_LENGTH_SOFT_LIMIT),
            "full_length_hard_limit": str(FULL_LENGTH_HARD_LIMIT),
            "combined_readiness_status": readiness,
            "design_mode_recommendation": design_mode,
            "packet_draft_status": packet_status,
            "missing_sequence_entities": ",".join(missing_entities),
            "variant_position_summary": " | ".join(per_gene_positions),
            "manual_domain_review_note": note,
            "next_action": (
                "build_manual_submission_packet"
                if packet_status == "ready_for_manual_packet_draft"
                else "resolve_domain_ranges_or_literature_review_before_packet"
            ),
        }
        design_rows.append(row)
        group_rows.append({
            "complex_group_id": cid,
            "complex_name": w["complex_name"],
            "genes_covered": row["genes_covered"],
            "protein_entities": row["protein_entities"],
            "ligand_entities": row["ligand_entities"],
            "ion_entities": row["ion_entities"],
            "dna_entities": "",
            "rna_entities": "",
            "total_tokens_estimated": str(total_tokens),
            "full_length_feasible": "true" if total_tokens <= FULL_LENGTH_HARD_LIMIT else "false",
            "domain_trim_required": "true" if design_mode in {"domain_or_split_required", "split_into_specific_domain_jobs"} else "review",
            "manual_submission_name": cid.replace("AF3_COMPLEX_", "AF3_"),
            "expected_variant_coverage": row["variant_position_summary"],
            "notes": note,
        })

    fields = [
        "complex_group_id", "complex_name", "genes_covered", "protein_entities",
        "ion_entities", "ligand_entities", "protein_entity_lengths",
        "estimated_af3_tokens", "full_length_soft_limit", "full_length_hard_limit",
        "combined_readiness_status", "design_mode_recommendation",
        "packet_draft_status", "missing_sequence_entities",
        "variant_position_summary", "manual_domain_review_note", "next_action",
    ]
    write_tsv(OUT_TSV, design_rows, fields)

    group_fields = [
        "complex_group_id", "complex_name", "genes_covered", "protein_entities",
        "ligand_entities", "ion_entities", "dna_entities", "rna_entities",
        "total_tokens_estimated", "full_length_feasible", "domain_trim_required",
        "manual_submission_name", "expected_variant_coverage", "notes",
    ]
    write_tsv(GROUPS_TSV, group_rows, group_fields)

    summary = {
        "date": date.today().isoformat(),
        "complexes_reviewed": len(design_rows),
        "full_length_soft_limit": FULL_LENGTH_SOFT_LIMIT,
        "full_length_hard_limit": FULL_LENGTH_HARD_LIMIT,
        "packet_status_counts": dict(__import__("collections").Counter(r["packet_draft_status"] for r in design_rows)),
        "design_mode_counts": dict(__import__("collections").Counter(r["design_mode_recommendation"] for r in design_rows)),
        "outputs": {
            "design_review_tsv": str(OUT_TSV),
            "complex_groups_tsv": str(GROUPS_TSV),
        },
        "safety_note": "This is a design-review table only; no AF3 server submission packets were generated.",
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
