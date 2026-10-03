#!/usr/bin/env python3
"""Annotate a small GRCh38 variant table with the official Ensembl VEP REST API."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import pandas as pd
import requests


CONSEQUENCE_SEVERITY = {
    "transcript_ablation": 1, "splice_acceptor_variant": 2, "splice_donor_variant": 3,
    "stop_gained": 4, "frameshift_variant": 5, "stop_lost": 6, "start_lost": 7,
    "transcript_amplification": 8, "inframe_insertion": 9, "inframe_deletion": 10,
    "missense_variant": 11, "protein_altering_variant": 12, "splice_region_variant": 13,
    "incomplete_terminal_codon_variant": 14, "start_retained_variant": 15,
    "stop_retained_variant": 16, "synonymous_variant": 17, "coding_sequence_variant": 18,
    "mature_miRNA_variant": 19, "5_prime_UTR_variant": 20, "3_prime_UTR_variant": 21,
    "non_coding_transcript_exon_variant": 22, "intron_variant": 23,
    "NMD_transcript_variant": 24, "non_coding_transcript_variant": 25,
    "upstream_gene_variant": 26, "downstream_gene_variant": 27,
    "TFBS_ablation": 28, "TFBS_amplification": 29, "TF_binding_site_variant": 30,
    "regulatory_region_ablation": 31, "regulatory_region_amplification": 32,
    "feature_elongation": 33, "regulatory_region_variant": 34,
    "feature_truncation": 35, "intergenic_variant": 36,
}
IMPACT_SEVERITY = {"HIGH": 1, "MODERATE": 2, "LOW": 3, "MODIFIER": 4}


def first(values: list[object]) -> object:
    for value in values:
        if value not in {None, "", "-"}:
            return value
    return ""


def joined(values: list[object], limit: int = 16) -> str:
    out = []
    for value in values:
        if isinstance(value, list):
            pieces = value
        else:
            pieces = str(value or "").replace("&", ",").split(",")
        for piece in pieces:
            piece = str(piece).strip()
            if piece and piece not in out:
                out.append(piece)
            if len(out) >= limit:
                return "|".join(out)
    return "|".join(out)


def representative(transcripts: list[dict]) -> dict:
    def rank(row: dict) -> tuple:
        consequences = row.get("consequence_terms", [])
        consequence_rank = min((CONSEQUENCE_SEVERITY.get(x, 998) for x in consequences), default=999)
        return (
            0 if row.get("canonical") else 1,
            0 if row.get("mane_select") else 1,
            IMPACT_SEVERITY.get(row.get("impact", ""), 999),
            consequence_rank,
        )
    return sorted(transcripts, key=rank)[0] if transcripts else {}


def collapse(result: dict, variant_id: str) -> dict[str, object]:
    transcripts = result.get("transcript_consequences", []) or []
    regulatory = result.get("regulatory_feature_consequences", []) or []
    intergenic = result.get("intergenic_consequences", []) or []
    all_rows = transcripts + regulatory + intergenic
    rep = representative(transcripts) if transcripts else (all_rows[0] if all_rows else {})
    terms = [term for row in all_rows for term in row.get("consequence_terms", [])]
    worst = min(terms, key=lambda x: CONSEQUENCE_SEVERITY.get(x, 998), default="")
    impacts = [row.get("impact", "") for row in all_rows]
    return {
        "variant_id": variant_id,
        "vep_annotation_status": "ok" if all_rows else "not_returned",
        "vep_transcript_count": len(transcripts),
        "vep_symbol": rep.get("gene_symbol", ""),
        "vep_gene_id": rep.get("gene_id", ""),
        "vep_feature": rep.get("transcript_id", rep.get("regulatory_feature_id", "")),
        "vep_feature_type": "Transcript" if rep in transcripts else "RegulatoryFeature",
        "vep_consequence": ",".join(rep.get("consequence_terms", [])),
        "vep_worst_consequence": worst,
        "vep_impact": rep.get("impact", ""),
        "vep_canonical": "YES" if rep.get("canonical") else "",
        "vep_mane_select": rep.get("mane_select", ""),
        "vep_mane_plus_clinical": rep.get("mane_plus_clinical", ""),
        "vep_hgvsc": rep.get("hgvsc", ""),
        "vep_hgvsp": rep.get("hgvsp", ""),
        "vep_protein_position": rep.get("protein_start", ""),
        "vep_amino_acids": rep.get("amino_acids", ""),
        "vep_codons": rep.get("codons", ""),
        "vep_existing_variation": joined(result.get("colocated_variants", [])),
        "vep_all_consequences": joined(terms),
        "vep_all_impacts": joined(impacts),
        "vep_any_high_impact": "HIGH" in impacts,
        "vep_any_moderate_impact": "MODERATE" in impacts,
        "vep_any_missense": "missense_variant" in terms,
        "vep_any_splice_region": any("splice" in term for term in terms),
        "vep_any_frameshift": "frameshift_variant" in terms,
        "vep_any_stop_gained": "stop_gained" in terms,
        "vep_annotation_missing_reason": "" if all_rows else "variant_not_returned_by_vep_rest",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--raw-json", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--chunk-size", type=int, default=50)
    args = parser.parse_args()

    variants = pd.read_csv(args.input, sep="\t", dtype=str)
    required = ["variant_id", "chrom", "pos", "ref", "alt"]
    if any(column not in variants for column in required):
        raise ValueError(f"Input must contain {required}")
    query_by_id = {
        row.variant_id: f"{str(row.chrom).removeprefix('chr')} {row.pos} . {row.ref} {row.alt} . . ."
        for row in variants[required].itertuples(index=False)
    }
    response_by_id: dict[str, dict] = {}
    session = requests.Session()
    endpoint = "https://rest.ensembl.org/vep/homo_sapiens/region"
    ids = list(query_by_id)
    for start in range(0, len(ids), args.chunk_size):
        chunk_ids = ids[start : start + args.chunk_size]
        payload = {
            "variants": [query_by_id[variant_id] for variant_id in chunk_ids],
            "canonical": 1, "mane": 1, "hgvs": 1, "symbol": 1, "protein": 1,
        }
        for attempt in range(5):
            response = session.post(
                endpoint,
                headers={"Content-Type": "application/json", "Accept": "application/json"},
                json=payload,
                timeout=180,
            )
            if response.status_code == 200:
                break
            if response.status_code not in {429, 500, 502, 503, 504} or attempt == 4:
                response.raise_for_status()
            time.sleep(2 ** attempt)
        returned = response.json()
        by_input = {item.get("input", ""): item for item in returned}
        for variant_id in chunk_ids:
            response_by_id[variant_id] = by_input.get(query_by_id[variant_id], {})
        time.sleep(0.2)

    rows = [collapse(response_by_id.get(variant_id, {}), variant_id) for variant_id in ids]
    output = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)
    args.raw_json.write_text(json.dumps(response_by_id, indent=2) + "\n")
    summary = {
        "endpoint": endpoint,
        "input_variants": len(ids),
        "annotated": int((output["vep_annotation_status"] == "ok").sum()),
        "not_returned": int((output["vep_annotation_status"] != "ok").sum()),
        "backend_note": "Official Ensembl VEP REST used because the local VEP 116 Perl environment segfaulted.",
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
