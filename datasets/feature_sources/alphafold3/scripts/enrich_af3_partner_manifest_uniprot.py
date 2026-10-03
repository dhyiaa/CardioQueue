#!/usr/bin/env python3
"""
Add UniProt-derived evidence to the AF3 starter partner manifest.

This creates cached UniProt JSON snapshots and a new evidence-enriched manifest.
It does not overwrite the starter manifest or any model matrix.
"""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SEQ_TSV = ROOT / "datasets/feature_sources/alphafold3/interim/uniprot_gene_to_sequence.tsv"
PARTNER_TSV = ROOT / "results/af3/af3_full_panel_partner_manifest.tsv"
RAW_DIR = ROOT / "datasets/feature_sources/alphafold3/raw/uniprot_evidence"
INTERIM_DIR = ROOT / "datasets/feature_sources/alphafold3/interim/evidence"
RESULTS_DIR = ROOT / "results/af3"

GENE_EVIDENCE_TSV = INTERIM_DIR / "uniprot_evidence_by_gene.tsv"
ENRICHED_TSV = RESULTS_DIR / "af3_full_panel_partner_manifest.uniprot_evidence.tsv"
SUMMARY_JSON = RESULTS_DIR / "af3_uniprot_partner_evidence.summary.json"

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fetch_uniprot(accession: str) -> dict:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = RAW_DIR / f"{accession}.json"
    if cache_path.exists():
        return json.loads(cache_path.read_text())
    url = UNIPROT_URL.format(accession=accession)
    req = urllib.request.Request(url, headers={"User-Agent": "cardiogenetics-af3-manifest/0.1"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=45) as response:
                data = json.load(response)
            cache_path.write_text(json.dumps(data, indent=2) + "\n")
            time.sleep(0.15)
            return data
        except (urllib.error.URLError, TimeoutError) as exc:
            if attempt == 2:
                raise RuntimeError(f"Failed UniProt fetch for {accession}: {exc}") from exc
            time.sleep(1 + attempt)
    raise RuntimeError(f"Failed UniProt fetch for {accession}")


def evidence_pmids_from_texts(texts: list[dict]) -> list[str]:
    pmids = []
    for text in texts:
        for ev in text.get("evidences", []):
            if ev.get("source") == "PubMed" and ev.get("id"):
                pmids.append(str(ev["id"]))
    return sorted(set(pmids), key=lambda x: int(x) if x.isdigit() else x)


def extract_gene_evidence(gene: str, seq_row: dict[str, str], data: dict) -> dict[str, str]:
    comments = data.get("comments", [])
    subunit_texts = []
    subunit_pmids = []
    function_texts = []
    function_pmids = []
    interaction_partners = []
    interaction_accessions = []
    intact_ids = []
    interaction_experiments = []

    for comment in comments:
        ctype = comment.get("commentType", "")
        if ctype == "SUBUNIT":
            texts = comment.get("texts", [])
            subunit_texts.extend(t.get("value", "") for t in texts if t.get("value"))
            subunit_pmids.extend(evidence_pmids_from_texts(texts))
        elif ctype in {"FUNCTION", "ACTIVITY REGULATION", "CATALYTIC ACTIVITY"}:
            texts = comment.get("texts", [])
            function_texts.extend(t.get("value", "") for t in texts if t.get("value"))
            function_pmids.extend(evidence_pmids_from_texts(texts))
        elif ctype == "INTERACTION":
            for inter in comment.get("interactions", []):
                two = inter.get("interactantTwo", {})
                one = inter.get("interactantOne", {})
                if two.get("geneName"):
                    interaction_partners.append(str(two["geneName"]).upper())
                if two.get("uniProtKBAccession"):
                    interaction_accessions.append(two["uniProtKBAccession"].split("-")[0])
                for participant in (one, two):
                    if participant.get("intActId"):
                        intact_ids.append(participant["intActId"])
                if "numberOfExperiments" in inter:
                    interaction_experiments.append(str(inter["numberOfExperiments"]))

    xrefs_by_db: dict[str, list[str]] = defaultdict(list)
    for xref in data.get("uniProtKBCrossReferences", []):
        db = xref.get("database", "")
        xid = xref.get("id", "")
        if db and xid:
            xrefs_by_db[db].append(xid)

    references = []
    for ref in data.get("references", []):
        citation = ref.get("citation", {})
        if citation.get("citationCrossReferences"):
            for xref in citation["citationCrossReferences"]:
                if xref.get("database") == "PubMed" and xref.get("id"):
                    references.append(str(xref["id"]))

    return {
        "gene": gene,
        "uniprot_accession": seq_row["uniprot_accession"],
        "entry_name": seq_row["entry_name"],
        "protein_name": seq_row["protein_name"],
        "sequence_length": seq_row["sequence_length"],
        "uniprot_url": f"https://rest.uniprot.org/uniprotkb/{seq_row['uniprot_accession']}.json",
        "subunit_text": " | ".join(subunit_texts),
        "subunit_pubmed_ids": ";".join(sorted(set(subunit_pmids), key=lambda x: int(x) if x.isdigit() else x)),
        "function_context_text": " | ".join(function_texts),
        "function_pubmed_ids": ";".join(sorted(set(function_pmids), key=lambda x: int(x) if x.isdigit() else x)),
        "interaction_partner_genes": ";".join(sorted(set(interaction_partners))),
        "interaction_partner_accessions": ";".join(sorted(set(interaction_accessions))),
        "interaction_experiment_counts": ";".join(interaction_experiments),
        "intact_ids": ";".join(sorted(set(intact_ids))),
        "reactome_ids": ";".join(sorted(set(xrefs_by_db.get("Reactome", [])))),
        "complexportal_ids": ";".join(sorted(set(xrefs_by_db.get("ComplexPortal", [])))),
        "pdb_ids": ";".join(sorted(set(xrefs_by_db.get("PDB", [])))),
        "corum_ids": ";".join(sorted(set(xrefs_by_db.get("CORUM", [])))),
        "biogrid_ids": ";".join(sorted(set(xrefs_by_db.get("BioGRID", [])))),
        "dip_ids": ";".join(sorted(set(xrefs_by_db.get("DIP", [])))),
        "all_reference_pubmed_ids": ";".join(sorted(set(references), key=lambda x: int(x) if x.isdigit() else x)),
    }


def contains_gene(text: str, gene: str) -> bool:
    if not text or not gene:
        return False
    return re.search(rf"(?<![A-Z0-9]){re.escape(gene.upper())}(?![A-Z0-9])", text.upper()) is not None


def classify_pair(row: dict[str, str], evidence_by_gene: dict[str, dict[str, str]]) -> tuple[str, str, str]:
    partner_type = row["partner_type"]
    target_gene = row["target_gene"].upper()
    partner = row["partner_gene_or_entity"].upper()
    target_ev = evidence_by_gene.get(target_gene, {})
    partner_ev = evidence_by_gene.get(partner, {})

    if partner_type == "protein":
        target_partners = set(filter(None, target_ev.get("interaction_partner_genes", "").split(";")))
        target_partner_accs = set(filter(None, target_ev.get("interaction_partner_accessions", "").split(";")))
        partner_acc = row.get("partner_uniprot_or_identifier", "").split("-")[0]
        if partner in target_partners or (partner_acc and partner_acc in target_partner_accs):
            return "direct_uniprot_intact_interaction", target_ev.get("intact_ids", ""), target_ev.get("subunit_pubmed_ids", "")

        target_subunit = target_ev.get("subunit_text", "")
        partner_subunit = partner_ev.get("subunit_text", "")
        if contains_gene(target_subunit, partner) or contains_gene(partner_subunit, target_gene):
            pmids = sorted(set(
                filter(None, (target_ev.get("subunit_pubmed_ids", "") + ";" + partner_ev.get("subunit_pubmed_ids", "")).split(";"))
            ), key=lambda x: int(x) if x.isdigit() else x)
            return "uniprot_subunit_text_mention", "", ";".join(pmids)

        if target_gene.startswith("CALM") and partner_ev:
            # CALM1/2/3 encode identical proteins in the local canonical FASTA.
            return "calmodulin_family_proxy_review_needed", "", target_ev.get("subunit_pubmed_ids", "")

        if partner.startswith("CALM") and target_ev:
            return "calmodulin_family_context_review_needed", "", target_ev.get("subunit_pubmed_ids", "")

        return "not_confirmed_by_uniprot_yet", "", ""

    context = " ".join([
        target_ev.get("subunit_text", ""),
        target_ev.get("function_context_text", ""),
    ]).upper()
    if partner_type == "ion" and partner in {"CA", "CA2+", "CALCIUM"}:
        if "CA(2+)" in context or "CALCIUM" in context:
            return "uniprot_calcium_context", "", target_ev.get("function_pubmed_ids", "")
        return "nonprotein_entity_review_needed", "", ""
    if partner_type == "ligand" and partner == "GTP":
        if "GTP" in context:
            return "uniprot_gtp_context", "", target_ev.get("function_pubmed_ids", "")
        return "nonprotein_entity_review_needed", "", ""
    if partner_type == "none_found":
        return "no_partner_assigned", "", ""
    return "nonprotein_entity_review_needed", "", ""


def main() -> None:
    seq_rows = read_tsv(SEQ_TSV)
    partner_rows = read_tsv(PARTNER_TSV)

    evidence_rows = []
    evidence_by_gene = {}
    for row in seq_rows:
        accession = row["uniprot_accession"]
        if not accession:
            continue
        data = fetch_uniprot(accession)
        ev = extract_gene_evidence(row["gene"], row, data)
        evidence_rows.append(ev)
        evidence_by_gene[row["gene"].upper()] = ev

    gene_fields = [
        "gene",
        "uniprot_accession",
        "entry_name",
        "protein_name",
        "sequence_length",
        "uniprot_url",
        "subunit_text",
        "subunit_pubmed_ids",
        "function_context_text",
        "function_pubmed_ids",
        "interaction_partner_genes",
        "interaction_partner_accessions",
        "interaction_experiment_counts",
        "intact_ids",
        "reactome_ids",
        "complexportal_ids",
        "pdb_ids",
        "corum_ids",
        "biogrid_ids",
        "dip_ids",
        "all_reference_pubmed_ids",
    ]
    write_tsv(GENE_EVIDENCE_TSV, evidence_rows, gene_fields)

    enriched_rows = []
    for row in partner_rows:
        support_level, intact_ids, pmids = classify_pair(row, evidence_by_gene)
        target_ev = evidence_by_gene.get(row["target_gene"].upper(), {})
        partner_ev = evidence_by_gene.get(row["partner_gene_or_entity"].upper(), {})
        enriched = dict(row)
        enriched["uniprot_pair_support_level"] = support_level
        enriched["uniprot_target_url"] = target_ev.get("uniprot_url", "")
        enriched["uniprot_partner_url"] = partner_ev.get("uniprot_url", "")
        enriched["uniprot_pair_intact_ids"] = intact_ids
        enriched["uniprot_pair_pubmed_ids"] = pmids
        enriched["uniprot_target_subunit_text"] = target_ev.get("subunit_text", "")
        enriched["uniprot_partner_subunit_text"] = partner_ev.get("subunit_text", "")
        enriched["uniprot_target_reactome_ids"] = target_ev.get("reactome_ids", "")
        enriched["uniprot_target_pdb_ids"] = target_ev.get("pdb_ids", "")
        enriched["evidence_review_date"] = date.today().isoformat()
        enriched_rows.append(enriched)

    enriched_fields = list(partner_rows[0].keys()) + [
        "uniprot_pair_support_level",
        "uniprot_target_url",
        "uniprot_partner_url",
        "uniprot_pair_intact_ids",
        "uniprot_pair_pubmed_ids",
        "uniprot_target_subunit_text",
        "uniprot_partner_subunit_text",
        "uniprot_target_reactome_ids",
        "uniprot_target_pdb_ids",
        "evidence_review_date",
    ]
    write_tsv(ENRICHED_TSV, enriched_rows, enriched_fields)

    summary = {
        "date": date.today().isoformat(),
        "sequence_entities_checked": len(seq_rows),
        "uniprot_json_cached": len(list(RAW_DIR.glob("*.json"))),
        "gene_evidence_rows": len(evidence_rows),
        "partner_rows_enriched": len(enriched_rows),
        "pair_support_counts": dict(Counter(r["uniprot_pair_support_level"] for r in enriched_rows)),
        "outputs": {
            "gene_evidence_tsv": str(GENE_EVIDENCE_TSV),
            "enriched_partner_manifest_tsv": str(ENRICHED_TSV),
            "raw_uniprot_json_dir": str(RAW_DIR),
        },
        "source": "UniProt REST API",
        "source_url": "https://rest.uniprot.org/uniprotkb/",
        "safety_note": "The starter partner manifest is preserved; this script writes an evidence-enriched copy only.",
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
