#!/usr/bin/env python3
"""Create AlphaFold Server JSON draft inputs for AF3 pilot-v1 jobs."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SEQ_TSVS = [
    ROOT / "datasets/feature_sources/alphafold3/interim/uniprot_gene_to_sequence.tsv",
    ROOT / "datasets/feature_sources/alphafold3/interim/af3_second_wave_entities_to_sequence.tsv",
]
PILOT_GENE_TSV = ROOT / "results/af3/pilot_v1/af3_pilot_v1_gene_set.tsv"
PILOT_COMPLEX_TSV = ROOT / "results/af3/pilot_v1/af3_pilot_v1_complex_plan.tsv"
HPA_GENE_TSV = ROOT / "datasets/feature_sources/alphafold3/interim/evidence/hpa_evidence_by_gene.tsv"

JSON_DIR = ROOT / "datasets/feature_sources/alphafold3/jobs/pilot_v1/server_json"
INDIVIDUAL_DIR = JSON_DIR / "individual_jobs"
ALL_JSON = JSON_DIR / "af3_pilot_v1_all_jobs_alphafoldserver.json"
TESTABLE_TSV = ROOT / "results/af3/pilot_v1/af3_pilot_v1_testable_genes_and_interactions.tsv"
SUMMARY_JSON = ROOT / "results/af3/pilot_v1/af3_pilot_v1_server_json.summary.json"


ION_MAP = {
    "CA": "CA",
    "Ca2+": "CA",
    "MG": "MG",
    "NA": "NA",
    "K": "K",
}

LIGAND_MAP = {
    "ATP": "CCD_ATP",
    "AMP": "CCD_AMP",
    "GTP": "CCD_GTP",
    "cAMP": "CCD_CMP",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def split_csv(value: str) -> list[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def load_sequences() -> dict[str, dict[str, str]]:
    seq = {}
    for path in SEQ_TSVS:
        for row in read_tsv(path):
            seq.setdefault(row["gene"], row)
    return seq


def load_hpa() -> dict[str, dict[str, str]]:
    return {row["gene"]: row for row in read_tsv(HPA_GENE_TSV)}


def server_job(row: dict[str, str], seq_by_gene: dict[str, dict[str, str]]) -> dict:
    sequences = []
    protein_counts = Counter(split_csv(row["protein_entities"]))
    for gene, count in protein_counts.items():
        seq = seq_by_gene[gene]["canonical_sequence"]
        sequences.append(
            {
                "proteinChain": {
                    "sequence": seq,
                    "count": count,
                }
            }
        )

    for entity in split_csv(row.get("nonprotein_entities", "")):
        entity = entity.strip()
        if not entity:
            continue
        if entity in ION_MAP:
            sequences.append({"ion": {"ion": ION_MAP[entity], "count": 1}})
        elif entity in LIGAND_MAP:
            sequences.append({"ligand": {"ligand": LIGAND_MAP[entity], "count": 1}})
        else:
            sequences.append({"ligand": {"ligand": entity, "count": 1}})

    return {
        "name": row["job_id"] or row["complex_group_id"].replace("AF3_COMPLEX_", "AF3_"),
        "modelSeeds": [],
        "sequences": sequences,
        "dialect": "alphafoldserver",
        "version": 1,
    }


def main() -> None:
    seq_by_gene = load_sequences()
    hpa_by_gene = load_hpa()
    genes = read_tsv(PILOT_GENE_TSV)
    complexes = read_tsv(PILOT_COMPLEX_TSV)

    INDIVIDUAL_DIR.mkdir(parents=True, exist_ok=True)

    jobs = []
    for row in complexes:
        missing = [g for g in split_csv(row["protein_entities"]) if g not in seq_by_gene]
        if missing:
            raise ValueError(f"Missing sequence for {row['complex_group_id']}: {missing}")
        job = server_job(row, seq_by_gene)
        jobs.append(job)
        out = INDIVIDUAL_DIR / f"{job['name']}.alphafoldserver.json"
        out.write_text(json.dumps([job], indent=2) + "\n")

    ALL_JSON.write_text(json.dumps(jobs, indent=2) + "\n")

    complexes_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in complexes:
        for gene in split_csv(row["target_genes"]):
            complexes_by_gene[gene].append(row)

    detail_rows = []
    for gene_row in genes:
        gene = gene_row["gene"]
        hpa = hpa_by_gene.get(gene, {})
        for complex_row in complexes_by_gene.get(gene, []):
            detail_rows.append(
                {
                    "gene": gene,
                    "pilot_stage": gene_row["pilot_stage"],
                    "complex_group_id": complex_row["complex_group_id"],
                    "job_id": complex_row["job_id"],
                    "complex_name": complex_row["complex_name"],
                    "protein_entities": complex_row["protein_entities"],
                    "nonprotein_entities": complex_row["nonprotein_entities"],
                    "interaction_design_type": "one_to_many_or_multichain" if len(split_csv(complex_row["protein_entities"])) > 2 else "one_to_one_or_homomer",
                    "submission_status": complex_row["submission_status"],
                    "hpa_status": hpa.get("hpa_status", ""),
                    "hpa_subcellular_main_location": hpa.get("hpa_subcellular_main_location", ""),
                    "hpa_subcellular_location": hpa.get("hpa_subcellular_location", ""),
                    "hpa_protein_tissue_specificity": hpa.get("hpa_protein_tissue_specificity", ""),
                    "hpa_protein_tissue_distribution": hpa.get("hpa_protein_tissue_distribution", ""),
                    "hpa_rna_tissue_specificity": hpa.get("hpa_rna_tissue_specificity", ""),
                    "hpa_rna_tissue_distribution": hpa.get("hpa_rna_tissue_distribution", ""),
                    "rationale_or_review_note": complex_row["notes"],
                    "server_json_file": str(INDIVIDUAL_DIR / f"{(complex_row['job_id'] or complex_row['complex_group_id'].replace('AF3_COMPLEX_', 'AF3_'))}.alphafoldserver.json"),
                }
            )

    fields = [
        "gene",
        "pilot_stage",
        "complex_group_id",
        "job_id",
        "complex_name",
        "protein_entities",
        "nonprotein_entities",
        "interaction_design_type",
        "submission_status",
        "hpa_status",
        "hpa_subcellular_main_location",
        "hpa_subcellular_location",
        "hpa_protein_tissue_specificity",
        "hpa_protein_tissue_distribution",
        "hpa_rna_tissue_specificity",
        "hpa_rna_tissue_distribution",
        "rationale_or_review_note",
        "server_json_file",
    ]
    write_tsv(TESTABLE_TSV, detail_rows, fields)

    summary = {
        "date": date.today().isoformat(),
        "json_dialect": "alphafoldserver",
        "job_count": len(jobs),
        "individual_json_count": len(jobs),
        "testable_gene_count": len({r["gene"] for r in detail_rows}),
        "testable_gene_complex_rows": len(detail_rows),
        "combined_json": str(ALL_JSON),
        "individual_json_dir": str(INDIVIDUAL_DIR),
        "testable_genes_tsv": str(TESTABLE_TSV),
        "source_note": "AlphaFold Server JSON is a top-level list of job dictionaries using dialect=alphafoldserver and version=1.",
        "review_note": "These JSON files are draft upload inputs. Review second-wave evidence packets and ligand/ion handling before submission.",
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
