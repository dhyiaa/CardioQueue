#!/usr/bin/env python3
"""
Enrich AF3 partner evidence with Reactome, RCSB/PDB, and IntAct/PSICQUIC.

Writes a new evidence table. Does not overwrite previous manifests or model data.
"""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.parse
import urllib.error
import urllib.request
from json import JSONDecodeError
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SEQ_TSV = ROOT / "datasets/feature_sources/alphafold3/interim/uniprot_gene_to_sequence.tsv"
UNIPROT_EVIDENCE_TSV = ROOT / "datasets/feature_sources/alphafold3/interim/evidence/uniprot_evidence_by_gene.tsv"
ENRICHED_IN = ROOT / "results/af3/af3_full_panel_partner_manifest.uniprot_evidence.tsv"

RAW_REACTOME = ROOT / "datasets/feature_sources/alphafold3/raw/reactome_evidence"
RAW_RCSB = ROOT / "datasets/feature_sources/alphafold3/raw/rcsb_evidence"
RAW_INTACT = ROOT / "datasets/feature_sources/alphafold3/raw/intact_evidence"
INTERIM_DIR = ROOT / "datasets/feature_sources/alphafold3/interim/evidence"
RESULTS_DIR = ROOT / "results/af3"

GENE_EXTERNAL_TSV = INTERIM_DIR / "external_evidence_by_gene.tsv"
PAIR_EXTERNAL_TSV = RESULTS_DIR / "af3_full_panel_partner_manifest.external_evidence.tsv"
NONE_FOUND_CANDIDATES_TSV = RESULTS_DIR / "af3_none_found_candidate_partners.tsv"
SUMMARY_JSON = RESULTS_DIR / "af3_external_partner_evidence.summary.json"

HEADERS = {
    "User-Agent": "cardiogenetics-af3-evidence/0.1",
    "Accept": "application/json,text/plain,*/*",
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


def get_url(url: str, cache_path: Path, binary: bool = False) -> bytes:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    if cache_path.exists():
        return cache_path.read_bytes()
    req = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(req, timeout=60) as response:
        data = response.read()
    cache_path.write_bytes(data)
    time.sleep(0.15)
    return data


def post_json(url: str, payload: dict, cache_path: Path) -> dict:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    if cache_path.exists():
        return json.loads(cache_path.read_text())
    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode(),
        headers={
            "User-Agent": HEADERS["User-Agent"],
            "Accept": "application/json",
            "Content-Type": "application/json",
        },
    )
    with urllib.request.urlopen(req, timeout=60) as response:
        raw = response.read()
    if not raw.strip():
        data = {"result_set": [], "total_count": 0, "empty_response": True}
        cache_path.write_text(json.dumps(data, indent=2) + "\n")
        time.sleep(0.15)
        return data
    try:
        data = json.loads(raw.decode())
    except JSONDecodeError:
        data = {
            "result_set": [],
            "total_count": 0,
            "non_json_response": raw.decode("utf-8", "ignore")[:500],
        }
    cache_path.write_text(json.dumps(data, indent=2) + "\n")
    time.sleep(0.15)
    return data


def fetch_reactome(accession: str) -> list[dict]:
    url = (
        "https://reactome.org/ContentService/data/mapping/UniProt/"
        f"{urllib.parse.quote(accession)}/pathways?species=9606"
    )
    cache_path = RAW_REACTOME / f"{accession}.pathways.json"
    try:
        raw = get_url(url, cache_path)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_text("[]\n")
            return []
        raise
    return json.loads(raw.decode())


def fetch_rcsb(accession: str) -> dict:
    payload = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": accession,
            },
        },
        "return_type": "polymer_entity",
        "request_options": {"paginate": {"start": 0, "rows": 500}},
    }
    return post_json("https://search.rcsb.org/rcsbsearch/v2/query", payload, RAW_RCSB / f"{accession}.polymer_entity_search.json")


def fetch_intact_psicquic(accession: str) -> str:
    url = (
        "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/"
        f"search/query/{urllib.parse.quote(accession)}?format=tab25&firstResult=0&maxResults=500"
    )
    return get_url(url, RAW_INTACT / f"{accession}.tab25.txt").decode("utf-8", "ignore")


def ids_from_field(field: str) -> set[str]:
    values = set()
    for token in field.split("|"):
        token = token.strip()
        if not token:
            continue
        if ":" in token:
            db, val = token.split(":", 1)
            values.add(val.split("(")[0].upper())
        values.add(token.upper())
    return values


def gene_names_from_aliases(field: str) -> set[str]:
    genes = set()
    for token in field.split("|"):
        token = token.strip()
        if not token:
            continue
        match = re.search(r"[:)]([A-Za-z0-9_-]+)\(gene name", token)
        if match:
            genes.add(match.group(1).upper())
        if "gene name" in token and ":" in token:
            genes.add(token.split(":", 1)[1].split("(")[0].upper())
    return genes


def pubmed_ids_from_field(field: str) -> list[str]:
    pmids = []
    for token in field.split("|"):
        if token.startswith("pubmed:"):
            pmids.append(token.split(":", 1)[1].split("(")[0])
    return sorted(set(pmids), key=sort_token)


def sort_token(value: str) -> tuple[int, int | str]:
    return (0, int(value)) if value.isdigit() else (1, value)


def parse_tab25(text: str, query_accession: str) -> list[dict[str, str]]:
    rows = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 15:
            continue
        ids_a = ids_from_field(parts[0] + "|" + parts[2])
        ids_b = ids_from_field(parts[1] + "|" + parts[3])
        genes_a = gene_names_from_aliases(parts[4])
        genes_b = gene_names_from_aliases(parts[5])
        if query_accession.upper() in ids_a:
            partner_accessions = ids_b
            partner_genes = genes_b
        elif query_accession.upper() in ids_b:
            partner_accessions = ids_a
            partner_genes = genes_a
        else:
            partner_accessions = ids_a | ids_b
            partner_genes = genes_a | genes_b
        rows.append({
            "query_accession": query_accession,
            "partner_accessions": ";".join(sorted(partner_accessions)),
            "partner_genes": ";".join(sorted(partner_genes)),
            "detection_method": parts[6],
            "publication_ids": ";".join(pubmed_ids_from_field(parts[8])),
            "taxid_a": parts[9],
            "taxid_b": parts[10],
            "interaction_type": parts[11],
            "source_db": parts[12],
            "interaction_ids": parts[13],
            "confidence": parts[14],
        })
    return rows


def pdb_codes_from_rcsb(data: dict) -> set[str]:
    codes = set()
    for item in data.get("result_set", []) or []:
        ident = item.get("identifier", "")
        if "_" in ident:
            codes.add(ident.split("_", 1)[0].upper())
    return codes


def shared(a: str, b: str) -> set[str]:
    return set(filter(None, a.split(";"))) & set(filter(None, b.split(";")))


def support_from_external(row: dict[str, str]) -> str:
    if row["partner_type"] == "none_found":
        return "no_partner_assigned"
    if row.get("intact_psicquic_pair_rows") not in {"", "0"}:
        return "intact_psicquic_pair_support"
    if row.get("rcsb_shared_pdb_ids"):
        return "rcsb_shared_pdb_support"
    if row.get("reactome_shared_pathway_ids"):
        return "reactome_shared_pathway_context"
    if row["partner_type"] in {"ion", "ligand"} and row.get("reactome_target_pathway_ids"):
        return "reactome_target_context_for_nonprotein_entity"
    return "not_supported_by_external_structured_sources_yet"


def main() -> None:
    seq_rows = read_tsv(SEQ_TSV)
    uniprot_evidence_rows = read_tsv(UNIPROT_EVIDENCE_TSV)
    pair_rows = read_tsv(ENRICHED_IN)
    seq_by_gene = {r["gene"].upper(): r for r in seq_rows}
    gene_by_accession = {r["uniprot_accession"].split("-")[0].upper(): r["gene"].upper() for r in seq_rows}

    external_by_gene: dict[str, dict[str, str]] = {}
    intact_rows_by_accession: dict[str, list[dict[str, str]]] = {}
    for row in seq_rows:
        gene = row["gene"].upper()
        acc = row["uniprot_accession"].split("-")[0]
        reactome = fetch_reactome(acc)
        rcsb = fetch_rcsb(acc)
        intact_text = fetch_intact_psicquic(acc)
        intact_rows = parse_tab25(intact_text, acc)
        intact_rows_by_accession[acc.upper()] = intact_rows

        reactome_ids = sorted({p.get("stId", "") for p in reactome if p.get("stId")})
        reactome_names = sorted({p.get("displayName", "") for p in reactome if p.get("displayName")})
        pdb_codes = sorted(pdb_codes_from_rcsb(rcsb))
        interact_gene_counts = Counter()
        interact_acc_counts = Counter()
        pubmed_ids = set()
        for ir in intact_rows:
            for g in filter(None, ir["partner_genes"].split(";")):
                interact_gene_counts[g] += 1
            for a in filter(None, ir["partner_accessions"].split(";")):
                interact_acc_counts[a] += 1
            pubmed_ids.update(filter(None, ir["publication_ids"].split(";")))
        external_by_gene[gene] = {
            "gene": gene,
            "uniprot_accession": acc,
            "reactome_pathway_ids": ";".join(reactome_ids),
            "reactome_pathway_names": " | ".join(reactome_names),
            "n_reactome_pathways": str(len(reactome_ids)),
            "rcsb_pdb_ids": ";".join(pdb_codes),
            "n_rcsb_pdb_ids": str(len(pdb_codes)),
            "intact_partner_genes": ";".join(g for g, _ in interact_gene_counts.most_common(100)),
            "intact_partner_accessions": ";".join(a for a, _ in interact_acc_counts.most_common(100)),
            "n_intact_psicquic_rows": str(len(intact_rows)),
            "intact_pubmed_ids": ";".join(sorted(pubmed_ids, key=sort_token)),
        }

    gene_fields = [
        "gene", "uniprot_accession", "reactome_pathway_ids", "reactome_pathway_names",
        "n_reactome_pathways", "rcsb_pdb_ids", "n_rcsb_pdb_ids", "intact_partner_genes",
        "intact_partner_accessions", "n_intact_psicquic_rows", "intact_pubmed_ids",
    ]
    write_tsv(GENE_EXTERNAL_TSV, [external_by_gene[g] for g in sorted(external_by_gene)], gene_fields)

    enriched_rows = []
    for row in pair_rows:
        out = dict(row)
        target_gene = row["target_gene"].upper()
        partner = row["partner_gene_or_entity"].upper()
        target_ext = external_by_gene.get(target_gene, {})
        partner_ext = external_by_gene.get(partner, {})
        target_acc = row["target_uniprot"].split("-")[0].upper()
        partner_acc = row["partner_uniprot_or_identifier"].split("-")[0].upper()

        reactome_shared_ids = shared(target_ext.get("reactome_pathway_ids", ""), partner_ext.get("reactome_pathway_ids", ""))
        rcsb_shared_ids = shared(target_ext.get("rcsb_pdb_ids", ""), partner_ext.get("rcsb_pdb_ids", ""))

        pair_intact_rows = []
        pair_pubmed = set()
        for ir in intact_rows_by_accession.get(target_acc, []):
            partner_genes = set(filter(None, ir["partner_genes"].split(";")))
            partner_accs = set(filter(None, ir["partner_accessions"].split(";")))
            if partner in partner_genes or partner_acc in partner_accs:
                pair_intact_rows.append(ir)
                pair_pubmed.update(filter(None, ir["publication_ids"].split(";")))

        out["reactome_target_pathway_ids"] = target_ext.get("reactome_pathway_ids", "")
        out["reactome_target_pathway_names"] = target_ext.get("reactome_pathway_names", "")
        out["reactome_partner_pathway_ids"] = partner_ext.get("reactome_pathway_ids", "")
        out["reactome_shared_pathway_ids"] = ";".join(sorted(reactome_shared_ids))
        out["rcsb_target_pdb_ids"] = target_ext.get("rcsb_pdb_ids", "")
        out["rcsb_partner_pdb_ids"] = partner_ext.get("rcsb_pdb_ids", "")
        out["rcsb_shared_pdb_ids"] = ";".join(sorted(rcsb_shared_ids))
        out["intact_psicquic_pair_rows"] = str(len(pair_intact_rows))
        out["intact_psicquic_pair_pubmed_ids"] = ";".join(sorted(pair_pubmed, key=sort_token))
        out["external_structured_support_level"] = support_from_external(out)
        enriched_rows.append(out)

    pair_fields = list(pair_rows[0].keys()) + [
        "reactome_target_pathway_ids",
        "reactome_target_pathway_names",
        "reactome_partner_pathway_ids",
        "reactome_shared_pathway_ids",
        "rcsb_target_pdb_ids",
        "rcsb_partner_pdb_ids",
        "rcsb_shared_pdb_ids",
        "intact_psicquic_pair_rows",
        "intact_psicquic_pair_pubmed_ids",
        "external_structured_support_level",
    ]
    write_tsv(PAIR_EXTERNAL_TSV, enriched_rows, pair_fields)

    panel_genes = set(seq_by_gene)
    candidate_rows = []
    none_found_genes = sorted({r["target_gene"].upper() for r in pair_rows if r["partner_type"] == "none_found"})
    for gene in none_found_genes:
        acc = seq_by_gene[gene]["uniprot_accession"].split("-")[0].upper()
        ext = external_by_gene[gene]
        candidates: dict[str, dict[str, str]] = {}
        for ir in intact_rows_by_accession.get(acc, []):
            pubmed = ir["publication_ids"]
            for cand_gene in filter(None, ir["partner_genes"].split(";")):
                cand_gene = cand_gene.upper()
                if cand_gene == gene:
                    continue
                rec = candidates.setdefault(cand_gene, {
                    "target_gene": gene,
                    "candidate_partner_gene": cand_gene,
                    "candidate_partner_source": "IntAct_PSICQUIC",
                    "n_intact_rows": 0,
                    "pubmed_ids": set(),
                    "candidate_in_current_sequence_manifest": "true" if cand_gene in panel_genes else "false",
                    "candidate_reactome_shared_pathway_ids": "",
                    "candidate_rcsb_shared_pdb_ids": "",
                    "manual_review_status": "manual_review_needed",
                })
                rec["n_intact_rows"] += 1
                rec["pubmed_ids"].update(filter(None, pubmed.split(";")))
        for cand_gene, rec in candidates.items():
            cand_ext = external_by_gene.get(cand_gene, {})
            rec["candidate_reactome_shared_pathway_ids"] = ";".join(sorted(shared(ext.get("reactome_pathway_ids", ""), cand_ext.get("reactome_pathway_ids", ""))))
            rec["candidate_rcsb_shared_pdb_ids"] = ";".join(sorted(shared(ext.get("rcsb_pdb_ids", ""), cand_ext.get("rcsb_pdb_ids", ""))))
        for rec in candidates.values():
            rec["pubmed_ids"] = ";".join(sorted(rec["pubmed_ids"], key=sort_token))
        top = sorted(candidates.values(), key=lambda r: (-r["n_intact_rows"], r["candidate_in_current_sequence_manifest"] != "true", r["candidate_partner_gene"]))[:10]
        candidate_rows.extend(top)

    cand_fields = [
        "target_gene", "candidate_partner_gene", "candidate_partner_source", "n_intact_rows",
        "pubmed_ids", "candidate_in_current_sequence_manifest", "candidate_reactome_shared_pathway_ids",
        "candidate_rcsb_shared_pdb_ids", "manual_review_status",
    ]
    write_tsv(NONE_FOUND_CANDIDATES_TSV, candidate_rows, cand_fields)

    summary = {
        "date": date.today().isoformat(),
        "sequence_entities_checked": len(seq_rows),
        "pair_rows_enriched": len(enriched_rows),
        "external_support_counts": dict(Counter(r["external_structured_support_level"] for r in enriched_rows)),
        "none_found_genes": none_found_genes,
        "none_found_candidate_rows": len(candidate_rows),
        "raw_cache_dirs": {
            "reactome": str(RAW_REACTOME),
            "rcsb": str(RAW_RCSB),
            "intact": str(RAW_INTACT),
        },
        "outputs": {
            "gene_external_evidence_tsv": str(GENE_EXTERNAL_TSV),
            "pair_external_evidence_tsv": str(PAIR_EXTERNAL_TSV),
            "none_found_candidate_partners_tsv": str(NONE_FOUND_CANDIDATES_TSV),
        },
        "sources": {
            "Reactome ContentService": "https://reactome.org/ContentService/",
            "RCSB Search API": "https://search.rcsb.org/rcsbsearch/v2/query",
            "IntAct PSICQUIC": "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/",
        },
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
