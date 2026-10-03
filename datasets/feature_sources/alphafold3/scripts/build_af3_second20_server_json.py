#!/usr/bin/env python3
"""Build the next 20 AlphaFold Server jobs after pilot-v1 batch QC."""

from __future__ import annotations

import csv
import json
from collections import Counter
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SEQ_TSVS = [
    ROOT / "datasets/feature_sources/alphafold3/interim/uniprot_gene_to_sequence.tsv",
    ROOT / "datasets/feature_sources/alphafold3/interim/af3_second_wave_entities_to_sequence.tsv",
]
OUT_DIR = ROOT / "datasets/feature_sources/alphafold3/jobs/pilot_v2_second20/server_json"
INDIVIDUAL_DIR = OUT_DIR / "individual_jobs"
COMBINED_JSON = OUT_DIR / "af3_pilot_v2_second20_alphafoldserver.json"
SUMMARY_TSV = ROOT / "results/af3/pilot_v2_second20/af3_pilot_v2_second20_job_summary.tsv"
PARTNER_SUMMARY_TSV = ROOT / "results/af3/pilot_v2_second20/af3_pilot_v2_partner_source_summary.tsv"
SUMMARY_JSON = ROOT / "results/af3/pilot_v2_second20/af3_pilot_v2_second20.summary.json"


ION_MAP = {"CA": "CA"}


JOBS = [
    {
        "job_id": "AF3_PILOTV2_001_BAG3_HSPB8_HSPA8_STUB1",
        "complex_name": "BAG3 CASA chaperone complex",
        "protein_entities": ["BAG3", "HSPB8", "HSPA8", "STUB1"],
        "nonprotein_entities": [],
        "target_genes": ["BAG3"],
        "source_decision": "second_wave_ready_for_evidence_check",
        "first10_context": "new_job_not_in_first10",
        "rationale": "BAG3 has curated CASA/chaperone complex context with HSPB8, HSPA8, and STUB1.",
    },
    {
        "job_id": "AF3_PILOTV2_002_LAMP2_HSPA8",
        "complex_name": "LAMP2-HSPA8 chaperone-mediated autophagy context",
        "protein_entities": ["LAMP2", "HSPA8"],
        "nonprotein_entities": [],
        "target_genes": ["LAMP2"],
        "source_decision": "second_wave_ready_for_evidence_check",
        "first10_context": "new_job_not_in_first10",
        "rationale": "LAMP2-HSPA8 is the compact protein-protein part of the chaperone-mediated autophagy context.",
    },
    {
        "job_id": "AF3_PILOTV2_003_SCN5A_SNTA1",
        "complex_name": "SCN5A-syntrophin channel-associated context",
        "protein_entities": ["SCN5A", "SNTA1"],
        "nonprotein_entities": [],
        "target_genes": ["SNTA1", "SCN5A"],
        "source_decision": "second_wave_ready_for_evidence_check",
        "first10_context": "new_scn5a_partner_after_scn5a_calm_and_scn5a_scn1b_worked",
        "rationale": "SCN5A gave strong interfaces with CALM and SCN1B in the first 10, so SNTA1 is a reasonable next SCN5A-partner test.",
    },
    {
        "job_id": "AF3_PILOTV2_004_TRPM4_HOMOTETRAMER",
        "complex_name": "TRPM4 homotetramer channel context",
        "protein_entities": ["TRPM4", "TRPM4", "TRPM4", "TRPM4"],
        "nonprotein_entities": [],
        "target_genes": ["TRPM4"],
        "source_decision": "second_wave_ready_for_evidence_check",
        "first10_context": "new_homotetramer_following_successful_kcnj2_homotetramer",
        "rationale": "KCNJ2 homotetramer worked well in the first 10; TRPM4 is another channel homomer candidate.",
    },
    {
        "job_id": "AF3_PILOTV2_005_KCNQ1_KCNE1_CALM1_CA",
        "complex_name": "KCNQ1-KCNE1-calmodulin-calcium IKs context",
        "protein_entities": ["KCNQ1", "KCNE1", "CALM1"],
        "nonprotein_entities": ["CA"],
        "target_genes": ["KCNQ1", "KCNE1", "CALM1", "CALM2", "CALM3"],
        "source_decision": "manual_cardiac_complex_whitelist",
        "first10_context": "new_core_cardiac_channel_complex",
        "rationale": "Major IKs cardiogenetics complex; use CALM1 chain as calmodulin proxy for CALM1/2/3 variant mapping.",
    },
    {
        "job_id": "AF3_PILOTV2_006_KCNQ1_KCNE1",
        "complex_name": "KCNQ1-KCNE1 direct IKs subunit context",
        "protein_entities": ["KCNQ1", "KCNE1"],
        "nonprotein_entities": [],
        "target_genes": ["KCNQ1", "KCNE1"],
        "source_decision": "manual_cardiac_complex_whitelist_pair_simplification",
        "first10_context": "paired_simplification_to_compare_with_multichain_job",
        "rationale": "Pair-only design helps diagnose whether the KCNE1 interface is clearer without calmodulin/calcium.",
    },
    {
        "job_id": "AF3_PILOTV2_007_KCNQ1_CALM1_CA",
        "complex_name": "KCNQ1-calmodulin-calcium context",
        "protein_entities": ["KCNQ1", "CALM1"],
        "nonprotein_entities": ["CA"],
        "target_genes": ["KCNQ1", "CALM1", "CALM2", "CALM3"],
        "source_decision": "manual_cardiac_complex_whitelist_pair_simplification",
        "first10_context": "paired_simplification_to_compare_with_multichain_job",
        "rationale": "Pair-only calmodulin design helps isolate KCNQ1-CALM signal.",
    },
    {
        "job_id": "AF3_PILOTV2_008_CACNA1C_CACNB2_CALM1_CA",
        "complex_name": "CACNA1C-CACNB2-calmodulin-calcium channel context",
        "protein_entities": ["CACNA1C", "CACNB2", "CALM1"],
        "nonprotein_entities": ["CA"],
        "target_genes": ["CACNA1C", "CACNB2", "CALM1", "CALM2", "CALM3"],
        "source_decision": "manual_cardiac_complex_whitelist",
        "first10_context": "new_core_cardiac_channel_complex",
        "rationale": "L-type calcium-channel context for CACNA1C/CACNB2 and calmodulin variants.",
    },
    {
        "job_id": "AF3_PILOTV2_009_CACNA1C_CACNB2",
        "complex_name": "CACNA1C-CACNB2 direct channel-subunit context",
        "protein_entities": ["CACNA1C", "CACNB2"],
        "nonprotein_entities": [],
        "target_genes": ["CACNA1C", "CACNB2"],
        "source_decision": "manual_cardiac_complex_whitelist_pair_simplification",
        "first10_context": "paired_simplification_to_compare_with_multichain_job",
        "rationale": "Pair-only design tests whether CACNB2 interface confidence is clearer without CALM/calcium.",
    },
    {
        "job_id": "AF3_PILOTV2_010_CASQ2_TRDN",
        "complex_name": "CASQ2-triadin calcium-release context",
        "protein_entities": ["CASQ2", "TRDN"],
        "nonprotein_entities": [],
        "target_genes": ["CASQ2", "TRDN"],
        "source_decision": "manual_cardiac_complex_whitelist_pair_simplification",
        "first10_context": "avoids_full_length_ryr2_token_and_domain_problem",
        "rationale": "RYR2 full-length is too large; CASQ2-TRDN is a compact triadic calcium-release pair.",
    },
    {
        "job_id": "AF3_PILOTV2_011_GLA_HOMODIMER",
        "complex_name": "GLA homodimer protein-only context",
        "protein_entities": ["GLA", "GLA"],
        "nonprotein_entities": [],
        "target_genes": ["GLA"],
        "source_decision": "second_wave_ligand_review_protein_only_first",
        "first10_context": "new_homomer_job",
        "rationale": "Use protein-only GLA homodimer first; glycosphingolipid substrate identity remains ligand-review work.",
    },
    {
        "job_id": "AF3_PILOTV2_012_HCN4_HOMOTETRAMER",
        "complex_name": "HCN4 homotetramer protein-only context",
        "protein_entities": ["HCN4", "HCN4", "HCN4", "HCN4"],
        "nonprotein_entities": [],
        "target_genes": ["HCN4"],
        "source_decision": "second_wave_ligand_review_protein_only_first",
        "first10_context": "new_homotetramer_following_successful_kcnj2_homotetramer",
        "rationale": "Protein-only HCN4 tetramer first; cAMP can be added after CCD review.",
    },
    {
        "job_id": "AF3_PILOTV2_013_PRKAA2_PRKAB2_PRKAG2",
        "complex_name": "AMPK heterotrimer protein-only PRKAG2 context",
        "protein_entities": ["PRKAA2", "PRKAB2", "PRKAG2"],
        "nonprotein_entities": [],
        "target_genes": ["PRKAG2"],
        "source_decision": "second_wave_ligand_review_protein_only_first",
        "first10_context": "new_multichain_metabolic_complex",
        "rationale": "Protein-only AMPK heterotrimer first; AMP can be added after ligand/CCD review.",
    },
    {
        "job_id": "AF3_PILOTV2_014_SLC22A5_PDZK1",
        "complex_name": "SLC22A5-PDZK1 transporter scaffold context",
        "protein_entities": ["SLC22A5", "PDZK1"],
        "nonprotein_entities": [],
        "target_genes": ["SLC22A5"],
        "source_decision": "second_wave_ligand_review_protein_only_first",
        "first10_context": "new_transporter_partner_job",
        "rationale": "Protein-only SLC22A5-PDZK1 first; L-carnitine ligand identity remains review work.",
    },
    {
        "job_id": "AF3_PILOTV2_015_FHOD3_SQSTM1",
        "complex_name": "FHOD3-SQSTM1 sarcomere/autophagy context",
        "protein_entities": ["FHOD3", "SQSTM1"],
        "nonprotein_entities": [],
        "target_genes": ["FHOD3"],
        "source_decision": "second_wave_cardiac_specificity_review",
        "first10_context": "new_review_job",
        "rationale": "Included as review-grade because the partner exists in UniProt, but cardiac specificity still needs caution.",
    },
    {
        "job_id": "AF3_PILOTV2_016_KCNH2_HOMOTETRAMER",
        "complex_name": "KCNH2/hERG homotetramer replacement context",
        "protein_entities": ["KCNH2", "KCNH2", "KCNH2", "KCNH2"],
        "nonprotein_entities": [],
        "target_genes": ["KCNH2"],
        "source_decision": "replacement_after_weak_kcnh2_kcne2",
        "first10_context": "kcnh2_kcne2_failed_interface_confidence",
        "rationale": "The first KCNH2-KCNE2 run was not interface-usable; test the core KCNH2 tetramer instead.",
    },
    {
        "job_id": "AF3_PILOTV2_017_DES_HOMODIMER",
        "complex_name": "DES homomeric filament assembly proxy",
        "protein_entities": ["DES", "DES"],
        "nonprotein_entities": [],
        "target_genes": ["DES"],
        "source_decision": "replacement_after_weak_des_cryab",
        "first10_context": "des_cryab_failed_interface_confidence",
        "rationale": "DES-CRYAB was weak; desmin homomeric assembly is a simpler protein-protein context.",
    },
    {
        "job_id": "AF3_PILOTV2_018_DSP_JUP",
        "complex_name": "DSP-JUP desmosome linker context",
        "protein_entities": ["DSP", "JUP"],
        "nonprotein_entities": [],
        "target_genes": ["DSP", "JUP"],
        "source_decision": "manual_cardiac_complex_whitelist_pair_simplification",
        "first10_context": "new_desmosome_pair",
        "rationale": "Compact desmosome pair that avoids oversized full desmosome assembly.",
    },
    {
        "job_id": "AF3_PILOTV2_019_DSG2_DSC2",
        "complex_name": "DSG2-DSC2 desmosome cadherin context",
        "protein_entities": ["DSG2", "DSC2"],
        "nonprotein_entities": [],
        "target_genes": ["DSG2", "DSC2"],
        "source_decision": "manual_cardiac_complex_whitelist_pair_simplification",
        "first10_context": "new_desmosome_pair",
        "rationale": "Compact cadherin pair for ARVC/desmosome genes.",
    },
    {
        "job_id": "AF3_PILOTV2_020_PKP2_DSP_JUP",
        "complex_name": "PKP2-DSP-JUP desmosome plaque context",
        "protein_entities": ["PKP2", "DSP", "JUP"],
        "nonprotein_entities": [],
        "target_genes": ["PKP2", "DSP", "JUP"],
        "source_decision": "manual_cardiac_complex_whitelist_reduced_multichain",
        "first10_context": "new_desmosome_reduced_multichain",
        "rationale": "Reduced desmosome plaque context below the full-core token burden.",
    },
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_sequences() -> dict[str, dict[str, str]]:
    seq = {}
    for path in SEQ_TSVS:
        for row in read_tsv(path):
            seq.setdefault(row["gene"], row)
    return seq


def server_job(row: dict[str, object], seq_by_gene: dict[str, dict[str, str]]) -> dict:
    sequences = []
    counts = Counter(row["protein_entities"])
    for gene, count in counts.items():
        sequences.append({"proteinChain": {"sequence": seq_by_gene[gene]["canonical_sequence"], "count": count}})
    for entity in row["nonprotein_entities"]:
        if entity in ION_MAP:
            sequences.append({"ion": {"ion": ION_MAP[entity], "count": 1}})
        else:
            sequences.append({"ligand": {"ligand": entity, "count": 1}})
    return {
        "name": row["job_id"],
        "modelSeeds": [],
        "sequences": sequences,
        "dialect": "alphafoldserver",
        "version": 1,
    }


def main() -> None:
    seq_by_gene = load_sequences()
    missing = sorted({g for row in JOBS for g in row["protein_entities"] if g not in seq_by_gene})
    if missing:
        raise ValueError(f"Missing sequences for second-20 jobs: {missing}")

    INDIVIDUAL_DIR.mkdir(parents=True, exist_ok=True)
    jobs = []
    summary_rows = []
    partner_rows = []
    for index, row in enumerate(JOBS, start=1):
        job = server_job(row, seq_by_gene)
        jobs.append(job)
        (INDIVIDUAL_DIR / f"{row['job_id']}.alphafoldserver.json").write_text(json.dumps([job], indent=2) + "\n")
        protein_lengths = [f"{g}:{seq_by_gene[g].get('sequence_length') or seq_by_gene[g].get('length')}" for g in row["protein_entities"]]
        unique_targets = ",".join(row["target_genes"])
        unique_proteins = ",".join(row["protein_entities"])
        token_estimate = sum(int(seq_by_gene[g].get("sequence_length") or seq_by_gene[g].get("length")) for g in row["protein_entities"]) + len(row["nonprotein_entities"])
        summary_rows.append(
            {
                "submission_order": index,
                "job_id": row["job_id"],
                "complex_name": row["complex_name"],
                "target_genes": unique_targets,
                "protein_entities": unique_proteins,
                "nonprotein_entities": ",".join(row["nonprotein_entities"]),
                "protein_lengths": ";".join(protein_lengths),
                "estimated_tokens": token_estimate,
                "source_decision": row["source_decision"],
                "first10_context": row["first10_context"],
                "rationale": row["rationale"],
                "server_json_file": str(INDIVIDUAL_DIR / f"{row['job_id']}.alphafoldserver.json"),
            }
        )
        for gene in row["target_genes"]:
            partner_rows.append(
                {
                    "target_gene": gene,
                    "job_id": row["job_id"],
                    "complex_name": row["complex_name"],
                    "protein_entities": unique_proteins,
                    "nonprotein_entities": ",".join(row["nonprotein_entities"]),
                    "source_decision": row["source_decision"],
                    "first10_context": row["first10_context"],
                    "rationale": row["rationale"],
                }
            )

    COMBINED_JSON.write_text(json.dumps(jobs, indent=2) + "\n")
    write_tsv(
        SUMMARY_TSV,
        summary_rows,
        [
            "submission_order",
            "job_id",
            "complex_name",
            "target_genes",
            "protein_entities",
            "nonprotein_entities",
            "protein_lengths",
            "estimated_tokens",
            "source_decision",
            "first10_context",
            "rationale",
            "server_json_file",
        ],
    )
    write_tsv(
        PARTNER_SUMMARY_TSV,
        partner_rows,
        [
            "target_gene",
            "job_id",
            "complex_name",
            "protein_entities",
            "nonprotein_entities",
            "source_decision",
            "first10_context",
            "rationale",
        ],
    )
    summary = {
        "date": date.today().isoformat(),
        "job_count": len(jobs),
        "combined_json": str(COMBINED_JSON),
        "individual_json_dir": str(INDIVIDUAL_DIR),
        "job_summary_tsv": str(SUMMARY_TSV),
        "partner_summary_tsv": str(PARTNER_SUMMARY_TSV),
        "safety_note": "This is an additive AF3 pilot-v2 batch only. It does not modify existing CatBoost/modeling matrices.",
        "first10_informed_note": "Weak first-10 jobs were not repeated unchanged; replacement/simplified designs are labeled in first10_context.",
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
