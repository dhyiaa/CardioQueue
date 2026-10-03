#!/usr/bin/env python3
"""Curate second-wave AF3 partner decisions for previously unresolved genes."""

from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
BASE_SEQ_TSV = ROOT / "datasets/feature_sources/alphafold3/interim/uniprot_gene_to_sequence.tsv"
UNIPROT_FASTA = ROOT / "datasets/feature_sources/protein_structure/raw/uniprot_human_reviewed_2026_02.fasta"
GENE_REVIEW_TSV = ROOT / "results/af3/af3_gene_level_coverage_review.tsv"
UNIPROT_EVIDENCE_TSV = ROOT / "datasets/feature_sources/alphafold3/interim/evidence/uniprot_evidence_by_gene.tsv"
EXTERNAL_EVIDENCE_TSV = ROOT / "datasets/feature_sources/alphafold3/interim/evidence/external_evidence_by_gene.tsv"

OUT_TSV = ROOT / "results/af3/af3_second_wave_curated_partner_decisions.tsv"
OUT_SEQ_TSV = ROOT / "datasets/feature_sources/alphafold3/interim/af3_second_wave_entities_to_sequence.tsv"
OUT_FASTA = ROOT / "datasets/feature_sources/alphafold3/interim/af3_second_wave_entities.fasta"
SUMMARY_JSON = ROOT / "results/af3/af3_second_wave_partner_decisions.summary.json"


SECOND_WAVE_DECISIONS = [
    {
        "target_gene": "ANK2",
        "decision": "assign_domain_review",
        "complex_group_id": "AF3_COMPLEX_ANK2_SPTBN1_DOMAIN",
        "complex_name": "ANK2-spectrin cytoskeletal adaptor domain context",
        "protein_entities": ["ANK2", "SPTBN1"],
        "nonprotein_entities": [],
        "rationale": "UniProt reports ANK2 interaction/colocalization with SPTBN1 and membrane ion-transport complexes; full proteins are too large, so this requires domain windows.",
        "evidence_anchor": "UniProt_SUBUNIT",
        "manual_review_status": "needs_domain_range_selection",
    },
    {
        "target_gene": "BAG3",
        "decision": "assign_full_length_candidate",
        "complex_group_id": "AF3_COMPLEX_BAG3_HSPB8_HSPA8_STUB1",
        "complex_name": "BAG3 CASA chaperone complex",
        "protein_entities": ["BAG3", "HSPB8", "HSPA8", "STUB1"],
        "nonprotein_entities": [],
        "rationale": "UniProt reports BAG3 in chaperone-assisted selective autophagy/CASA complex with HSPA8, HSPB8, and STUB1.",
        "evidence_anchor": "UniProt_SUBUNIT",
        "manual_review_status": "ready_for_evidence_check",
    },
    {
        "target_gene": "CSRP3",
        "decision": "assign_full_length_candidate",
        "complex_group_id": "AF3_COMPLEX_CSRP3_TCAP_ACTN2",
        "complex_name": "CSRP3-telethonin-alpha-actinin Z-disc context",
        "protein_entities": ["CSRP3", "TCAP", "ACTN2"],
        "nonprotein_entities": [],
        "rationale": "UniProt reports CSRP3 interactions with TCAP and ACTN2; this is a cardiac Z-disc context.",
        "evidence_anchor": "UniProt_SUBUNIT",
        "manual_review_status": "ready_for_evidence_check",
    },
    {
        "target_gene": "DES",
        "decision": "assign_full_length_candidate",
        "complex_group_id": "AF3_COMPLEX_DES_CRYAB",
        "complex_name": "Desmin-alphaB-crystallin intermediate-filament stress context",
        "protein_entities": ["DES", "CRYAB"],
        "nonprotein_entities": [],
        "rationale": "UniProt reports DES interaction with CRYAB and desmin homomeric filament assembly.",
        "evidence_anchor": "UniProt_SUBUNIT",
        "manual_review_status": "ready_for_evidence_check",
    },
    {
        "target_gene": "FHOD3",
        "decision": "assign_full_length_candidate",
        "complex_group_id": "AF3_COMPLEX_FHOD3_SQSTM1",
        "complex_name": "FHOD3-SQSTM1 sarcomere/autophagy context",
        "protein_entities": ["FHOD3", "SQSTM1"],
        "nonprotein_entities": [],
        "rationale": "UniProt reports FHOD3 interaction with SQSTM1; keep as review-needed because cardiac specificity requires checking.",
        "evidence_anchor": "UniProt_SUBUNIT",
        "manual_review_status": "needs_cardiac_specificity_review",
    },
    {
        "target_gene": "FLNC",
        "decision": "assign_domain_review",
        "complex_group_id": "AF3_COMPLEX_FLNC_MYOT_DOMAIN",
        "complex_name": "FLNC-myotilin myofibrillar context",
        "protein_entities": ["FLNC", "MYOT"],
        "nonprotein_entities": [],
        "rationale": "UniProt/IntAct report FLNC interaction with MYOT; FLNC is large and should use domain windows.",
        "evidence_anchor": "UniProt_SUBUNIT_IntAct",
        "manual_review_status": "needs_domain_range_selection",
    },
    {
        "target_gene": "GLA",
        "decision": "assign_homomer_ligand_context",
        "complex_group_id": "AF3_COMPLEX_GLA_HOMODIMER_GLYCOSPHINGOLIPID",
        "complex_name": "GLA homodimer glycosphingolipid-catabolism context",
        "protein_entities": ["GLA", "GLA"],
        "nonprotein_entities": ["glycosphingolipid_substrate_review"],
        "rationale": "UniProt reports GLA homodimer; Reactome maps GLA to glycosphingolipid catabolism. Use homodimer/substrate context rather than weak protein partners.",
        "evidence_anchor": "UniProt_SUBUNIT_Reactome",
        "manual_review_status": "needs_ligand_ccd_review",
    },
    {
        "target_gene": "HCN4",
        "decision": "assign_homomer_ligand_context",
        "complex_group_id": "AF3_COMPLEX_HCN4_HOMOTETRAMER_CAMP",
        "complex_name": "HCN4 homotetramer cyclic-nucleotide context",
        "protein_entities": ["HCN4", "HCN4", "HCN4", "HCN4"],
        "nonprotein_entities": ["cAMP"],
        "rationale": "UniProt reports HCN4 homotetramer and cAMP-regulated channel behavior; use homotetramer/cyclic nucleotide context.",
        "evidence_anchor": "UniProt_SUBUNIT_Reactome_RCSB",
        "manual_review_status": "needs_ligand_ccd_review",
    },
    {
        "target_gene": "KCNE2",
        "decision": "assign_full_length_candidate",
        "complex_group_id": "AF3_COMPLEX_KCNH2_KCNE2",
        "complex_name": "KCNH2-KCNE2 repolarization-channel complex",
        "protein_entities": ["KCNH2", "KCNE2"],
        "nonprotein_entities": [],
        "rationale": "UniProt reports KCNE2 association with KCNH2/ERG1 and KCNH2 stable complex with KCNE2.",
        "evidence_anchor": "UniProt_SUBUNIT_Reactome",
        "manual_review_status": "ready_for_evidence_check",
    },
    {
        "target_gene": "KCNH2",
        "decision": "assign_full_length_candidate",
        "complex_group_id": "AF3_COMPLEX_KCNH2_KCNE2",
        "complex_name": "KCNH2-KCNE2 repolarization-channel complex",
        "protein_entities": ["KCNH2", "KCNE2"],
        "nonprotein_entities": [],
        "rationale": "UniProt reports KCNH2 forms a stable complex with KCNE2 and Reactome maps both to repolarization pathways.",
        "evidence_anchor": "UniProt_SUBUNIT_Reactome",
        "manual_review_status": "ready_for_evidence_check",
    },
    {
        "target_gene": "KCNJ2",
        "decision": "assign_homomer_context",
        "complex_group_id": "AF3_COMPLEX_KCNJ2_HOMOTETRAMER",
        "complex_name": "KCNJ2/Kir2.1 homotetramer channel context",
        "protein_entities": ["KCNJ2", "KCNJ2", "KCNJ2", "KCNJ2"],
        "nonprotein_entities": [],
        "rationale": "UniProt reports KCNJ2 homotetramer and Reactome maps it to classical Kir channels/resting membrane potential.",
        "evidence_anchor": "UniProt_SUBUNIT_Reactome_RCSB",
        "manual_review_status": "ready_for_evidence_check",
    },
    {
        "target_gene": "LAMP2",
        "decision": "assign_full_length_candidate",
        "complex_group_id": "AF3_COMPLEX_LAMP2_HSPA8",
        "complex_name": "LAMP2-HSPA8 chaperone-mediated autophagy context",
        "protein_entities": ["LAMP2", "HSPA8"],
        "nonprotein_entities": [],
        "rationale": "UniProt reports LAMP2 interaction with HSPA8 in chaperone-mediated autophagy; Reactome supports CMA pathway context.",
        "evidence_anchor": "UniProt_SUBUNIT_Reactome",
        "manual_review_status": "ready_for_evidence_check",
    },
    {
        "target_gene": "MIB1",
        "decision": "assign_manual_review",
        "complex_group_id": "AF3_COMPLEX_MIB1_NOTCH_CONTEXT_REVIEW",
        "complex_name": "MIB1 Notch-signaling context",
        "protein_entities": ["MIB1"],
        "nonprotein_entities": [],
        "rationale": "Reactome places MIB1 in Notch signaling, but a compact direct cardiac AF3 partner is not yet settled.",
        "evidence_anchor": "Reactome",
        "manual_review_status": "needs_partner_selection",
    },
    {
        "target_gene": "MYBPC3",
        "decision": "assign_domain_review",
        "complex_group_id": "AF3_COMPLEX_MYOSIN_THICK_FILAMENT",
        "complex_name": "Cardiac myosin thick-filament context",
        "protein_entities": ["MYH7", "MYBPC3", "MYL2", "MYL3"],
        "nonprotein_entities": [],
        "rationale": "IntAct/Reactome support MYBPC3 with ACTC1/MYH7 and striated-muscle contraction; use myosin thick-filament context with domain review.",
        "evidence_anchor": "IntAct_Reactome_RCSB",
        "manual_review_status": "needs_domain_range_selection",
    },
    {
        "target_gene": "MYH7",
        "decision": "assign_domain_review",
        "complex_group_id": "AF3_COMPLEX_MYOSIN_THICK_FILAMENT",
        "complex_name": "Cardiac myosin thick-filament context",
        "protein_entities": ["MYH7", "MYBPC3", "MYL2", "MYL3"],
        "nonprotein_entities": [],
        "rationale": "UniProt reports myosin hexamer with heavy and light chains; MYBPC3/MYH7 have IntAct/Reactome support.",
        "evidence_anchor": "UniProt_SUBUNIT_IntAct_Reactome_RCSB",
        "manual_review_status": "needs_domain_range_selection",
    },
    {
        "target_gene": "MYL2",
        "decision": "assign_domain_review",
        "complex_group_id": "AF3_COMPLEX_MYOSIN_THICK_FILAMENT",
        "complex_name": "Cardiac myosin thick-filament context",
        "protein_entities": ["MYH7", "MYBPC3", "MYL2", "MYL3"],
        "nonprotein_entities": [],
        "rationale": "UniProt reports myosin as hexamer with heavy and light chains; Reactome/PDB support striated muscle contraction context.",
        "evidence_anchor": "UniProt_SUBUNIT_Reactome_RCSB",
        "manual_review_status": "needs_domain_range_selection",
    },
    {
        "target_gene": "MYL3",
        "decision": "assign_domain_review",
        "complex_group_id": "AF3_COMPLEX_MYOSIN_THICK_FILAMENT",
        "complex_name": "Cardiac myosin thick-filament context",
        "protein_entities": ["MYH7", "MYBPC3", "MYL2", "MYL3"],
        "nonprotein_entities": [],
        "rationale": "UniProt reports myosin as hexamer with heavy and light chains; Reactome/PDB support striated muscle contraction context.",
        "evidence_anchor": "UniProt_SUBUNIT_Reactome_RCSB",
        "manual_review_status": "needs_domain_range_selection",
    },
    {
        "target_gene": "PRKAG2",
        "decision": "assign_full_length_candidate",
        "complex_group_id": "AF3_COMPLEX_PRKAA2_PRKAB2_PRKAG2_AMP",
        "complex_name": "AMPK heterotrimer PRKAG2 context",
        "protein_entities": ["PRKAA2", "PRKAB2", "PRKAG2"],
        "nonprotein_entities": ["AMP"],
        "rationale": "UniProt reports AMPK heterotrimer of alpha, beta, and gamma subunits; IntAct supports PRKAA/PRKAB partners.",
        "evidence_anchor": "UniProt_SUBUNIT_IntAct_Reactome",
        "manual_review_status": "needs_ligand_ccd_review",
    },
    {
        "target_gene": "RBM20",
        "decision": "defer_or_rna_complex_review",
        "complex_group_id": "AF3_COMPLEX_RBM20_SPLICEOSOME_REVIEW",
        "complex_name": "RBM20 spliceosome/RNA-binding context",
        "protein_entities": ["RBM20"],
        "nonprotein_entities": ["RNA_context_review"],
        "rationale": "UniProt reports association with U1/U2 snRNP complexes, but a specific AF3 protein/RNA entity set needs separate spliceosome-focused review.",
        "evidence_anchor": "UniProt_SUBUNIT",
        "manual_review_status": "defer_until_rna_complex_design",
    },
    {
        "target_gene": "SCN1B",
        "decision": "assign_full_length_candidate",
        "complex_group_id": "AF3_COMPLEX_SCN5A_SCN1B",
        "complex_name": "SCN5A/Nav1.5 beta-1 subunit context",
        "protein_entities": ["SCN5A", "SCN1B"],
        "nonprotein_entities": [],
        "rationale": "UniProt reports SCN1B interaction/regulatory subunit relationship with SCN5A/Nav1.5.",
        "evidence_anchor": "UniProt_SUBUNIT_Reactome",
        "manual_review_status": "ready_for_evidence_check",
    },
    {
        "target_gene": "SLC22A5",
        "decision": "assign_full_length_candidate",
        "complex_group_id": "AF3_COMPLEX_SLC22A5_PDZK1_CARNITINE",
        "complex_name": "SLC22A5-PDZK1 carnitine-transporter context",
        "protein_entities": ["SLC22A5", "PDZK1"],
        "nonprotein_entities": ["L-carnitine"],
        "rationale": "UniProt reports SLC22A5 interaction with PDZK1 and Reactome supports carnitine shuttle/transport context.",
        "evidence_anchor": "UniProt_SUBUNIT_Reactome_RCSB",
        "manual_review_status": "needs_ligand_ccd_review",
    },
    {
        "target_gene": "SNTA1",
        "decision": "assign_full_length_candidate",
        "complex_group_id": "AF3_COMPLEX_SCN5A_SNTA1",
        "complex_name": "SCN5A-syntrophin channel-associated context",
        "protein_entities": ["SCN5A", "SNTA1"],
        "nonprotein_entities": [],
        "rationale": "UniProt reports SNTA1 interaction with sodium channel proteins including SCN5A and DGC context.",
        "evidence_anchor": "UniProt_SUBUNIT_Reactome",
        "manual_review_status": "ready_for_evidence_check",
    },
    {
        "target_gene": "TCAP",
        "decision": "assign_full_length_candidate",
        "complex_group_id": "AF3_COMPLEX_CSRP3_TCAP_ACTN2",
        "complex_name": "CSRP3-telethonin-alpha-actinin Z-disc context",
        "protein_entities": ["CSRP3", "TCAP", "ACTN2"],
        "nonprotein_entities": [],
        "rationale": "UniProt reports TCAP interaction with CSRP3 and titin, and CSRP3 interaction with TCAP/ACTN2; use compact Z-disc context.",
        "evidence_anchor": "UniProt_SUBUNIT_Reactome_RCSB",
        "manual_review_status": "ready_for_evidence_check",
    },
    {
        "target_gene": "TMEM43",
        "decision": "assign_full_length_candidate",
        "complex_group_id": "AF3_COMPLEX_TMEM43_EMD_LMNA",
        "complex_name": "TMEM43-emerin-lamin nuclear-envelope context",
        "protein_entities": ["TMEM43", "EMD", "LMNA"],
        "nonprotein_entities": [],
        "rationale": "UniProt reports TMEM43 interaction with EMD, LMNA, LMNB2, and SUN2.",
        "evidence_anchor": "UniProt_SUBUNIT_IntAct",
        "manual_review_status": "ready_for_evidence_check",
    },
    {
        "target_gene": "TRPM4",
        "decision": "assign_homomer_context",
        "complex_group_id": "AF3_COMPLEX_TRPM4_HOMOTETRAMER",
        "complex_name": "TRPM4 homotetramer channel context",
        "protein_entities": ["TRPM4", "TRPM4", "TRPM4", "TRPM4"],
        "nonprotein_entities": [],
        "rationale": "UniProt reports TRPM4 homotetramer and Reactome maps TRPM4 to TRP channels.",
        "evidence_anchor": "UniProt_SUBUNIT_Reactome_RCSB",
        "manual_review_status": "ready_for_evidence_check",
    },
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_fasta(path: Path) -> dict[str, list[dict[str, str]]]:
    by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    header_re = re.compile(r"^>([^|]+)\|([^|]+)\|(\S+)\s+(.+)$")
    current = None
    seq_parts: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    current["canonical_sequence"] = "".join(seq_parts)
                    current["sequence_length"] = str(len(current["canonical_sequence"]))
                    by_gene[current["gene"]].append(current)
                m = header_re.match(line)
                if not m:
                    raise ValueError(line)
                db, accession, entry_name, desc = m.groups()
                gene_match = re.search(r"\bGN=([^\s]+)", desc)
                current = {
                    "gene": gene_match.group(1).upper() if gene_match else "",
                    "uniprot_accession": accession,
                    "entry_name": entry_name,
                    "protein_name": desc.split(" OS=")[0],
                    "organism": "Homo sapiens" if "OS=Homo sapiens" in desc else "",
                    "reviewed": "true" if db == "sp" else "false",
                    "canonical_sequence": "",
                    "sequence_length": "",
                    "retrieval_source": str(path),
                }
                seq_parts = []
            else:
                seq_parts.append(line)
    if current:
        current["canonical_sequence"] = "".join(seq_parts)
        current["sequence_length"] = str(len(current["canonical_sequence"]))
        by_gene[current["gene"]].append(current)
    return by_gene


def select_gene_sequence(gene: str, local_seq: dict[str, dict[str, str]], fasta_by_gene: dict[str, list[dict[str, str]]]) -> dict[str, str] | None:
    if gene in local_seq:
        return local_seq[gene]
    candidates = fasta_by_gene.get(gene, [])
    if not candidates:
        return None
    candidates = sorted(candidates, key=lambda r: (int(r["sequence_length"]), r["uniprot_accession"]))
    return candidates[-1]


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def evidence_for_gene(gene: str, uni: dict[str, dict[str, str]], ext: dict[str, dict[str, str]]) -> dict[str, str]:
    u = uni.get(gene, {})
    e = ext.get(gene, {})
    return {
        "uniprot_pubmed_ids": u.get("subunit_pubmed_ids") or u.get("function_pubmed_ids", ""),
        "intact_pubmed_ids": e.get("intact_pubmed_ids", ""),
        "reactome_ids": e.get("reactome_pathway_ids", ""),
        "pdb_ids": e.get("rcsb_pdb_ids", ""),
        "uniprot_url": u.get("uniprot_url", ""),
    }


def main() -> None:
    base_seq_rows = read_tsv(BASE_SEQ_TSV)
    local_seq = {r["gene"]: r for r in base_seq_rows}
    fasta_by_gene = parse_fasta(UNIPROT_FASTA)
    uni_evidence = {r["gene"]: r for r in read_tsv(UNIPROT_EVIDENCE_TSV)}
    ext_evidence = {r["gene"]: r for r in read_tsv(EXTERNAL_EVIDENCE_TSV)}

    all_entities = sorted({e for d in SECOND_WAVE_DECISIONS for e in d["protein_entities"]})
    seq_rows = []
    missing = []
    for gene in all_entities:
        seq = select_gene_sequence(gene, local_seq, fasta_by_gene)
        if not seq:
            missing.append(gene)
            continue
        seq_rows.append({
            "gene": gene,
            "uniprot_accession": seq["uniprot_accession"],
            "entry_name": seq["entry_name"],
            "protein_name": seq["protein_name"],
            "sequence_length": seq["sequence_length"],
            "canonical_sequence": seq["canonical_sequence"],
            "retrieval_source": seq.get("retrieval_source", ""),
            "source_status": "existing_af3_sequence" if gene in local_seq else "new_external_partner_sequence",
        })

    seq_by_gene = {r["gene"]: r for r in seq_rows}
    out_rows = []
    for d in SECOND_WAVE_DECISIONS:
        lengths = []
        accessions = []
        entity_missing = []
        pmids = []
        intact_pmids = []
        reactome = []
        pdbs = []
        uniprot_urls = []
        for gene in sorted(set(d["protein_entities"])):
            seq = seq_by_gene.get(gene)
            if seq:
                lengths.append(f"{gene}:{seq['sequence_length']}")
                accessions.append(f"{gene}:{seq['uniprot_accession']}")
            else:
                entity_missing.append(gene)
            ev = evidence_for_gene(gene, uni_evidence, ext_evidence)
            pmids.extend(x for x in ev["uniprot_pubmed_ids"].split(";") if x)
            intact_pmids.extend(x for x in ev["intact_pubmed_ids"].split(";") if x)
            reactome.extend(x for x in ev["reactome_ids"].split(";") if x)
            pdbs.extend(x for x in ev["pdb_ids"].split(";") if x)
            if ev["uniprot_url"]:
                uniprot_urls.append(ev["uniprot_url"])
        total_tokens = sum(int(x.split(":")[1]) for x in lengths)
        out = dict(d)
        out["protein_entities"] = ",".join(d["protein_entities"])
        out["nonprotein_entities"] = ",".join(d["nonprotein_entities"])
        out["protein_entity_lengths"] = ",".join(lengths)
        out["protein_entity_accessions"] = ",".join(accessions)
        out["estimated_tokens"] = str(total_tokens + len(d["nonprotein_entities"]))
        out["missing_sequence_entities"] = ",".join(entity_missing)
        out["pubmed_ids_from_uniprot_gene_records"] = ";".join(sorted(set(pmids), key=lambda x: (0, int(x)) if x.isdigit() else (1, x)))
        out["pubmed_ids_from_intact_gene_records"] = ";".join(sorted(set(intact_pmids), key=lambda x: (0, int(x)) if x.isdigit() else (1, x)))
        out["reactome_ids_from_gene_records"] = ";".join(sorted(set(reactome)))
        out["pdb_ids_from_gene_records"] = ";".join(sorted(set(pdbs)))
        out["uniprot_urls"] = ";".join(sorted(set(uniprot_urls)))
        out["evidence_source_urls"] = "https://rest.uniprot.org/uniprotkb/;https://reactome.org/ContentService/;https://search.rcsb.org/rcsbsearch/v2/query;https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/"
        out_rows.append(out)

    fields = [
        "target_gene", "decision", "complex_group_id", "complex_name",
        "protein_entities", "nonprotein_entities", "protein_entity_lengths",
        "protein_entity_accessions", "estimated_tokens", "missing_sequence_entities",
        "rationale", "evidence_anchor", "manual_review_status",
        "pubmed_ids_from_uniprot_gene_records", "pubmed_ids_from_intact_gene_records",
        "reactome_ids_from_gene_records", "pdb_ids_from_gene_records",
        "uniprot_urls", "evidence_source_urls",
    ]
    write_tsv(OUT_TSV, out_rows, fields)

    seq_fields = [
        "gene", "uniprot_accession", "entry_name", "protein_name",
        "sequence_length", "canonical_sequence", "retrieval_source", "source_status",
    ]
    write_tsv(OUT_SEQ_TSV, seq_rows, seq_fields)

    with OUT_FASTA.open("w") as handle:
        for r in seq_rows:
            handle.write(f">{r['gene']}|{r['uniprot_accession']}|{r['entry_name']} {r['protein_name']}\n")
            handle.write(wrap(r["canonical_sequence"]) + "\n")

    summary = {
        "date": date.today().isoformat(),
        "second_wave_decision_rows": len(out_rows),
        "unique_target_genes": len({r["target_gene"] for r in out_rows}),
        "unique_complex_groups": len({r["complex_group_id"] for r in out_rows}),
        "decision_counts": dict(Counter(r["decision"] for r in out_rows)),
        "manual_review_status_counts": dict(Counter(r["manual_review_status"] for r in out_rows)),
        "sequence_entities": len(seq_rows),
        "new_external_partner_sequences": [r["gene"] for r in seq_rows if r["source_status"] == "new_external_partner_sequence"],
        "missing_sequence_entities": missing,
        "outputs": {
            "second_wave_decisions": str(OUT_TSV),
            "second_wave_sequences": str(OUT_SEQ_TSV),
            "second_wave_fasta": str(OUT_FASTA),
        },
        "safety_note": "Second-wave decisions are curated review decisions only; no AF3 server jobs were submitted.",
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
