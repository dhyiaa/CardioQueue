#!/usr/bin/env python3
"""Add Human Protein Atlas localization/expression evidence to AF3 manifests."""

from __future__ import annotations

import csv
import json
import zipfile
from collections import Counter
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
HPA_ZIP = ROOT / "datasets/feature_sources/alphafold3/raw/hpa_evidence/proteinatlas.tsv.zip"
BASE_PAIR_TSV = ROOT / "results/af3/af3_full_panel_partner_manifest.external_evidence.tsv"
OUT_GENE_TSV = ROOT / "datasets/feature_sources/alphafold3/interim/evidence/hpa_evidence_by_gene.tsv"
OUT_PAIR_TSV = ROOT / "results/af3/af3_full_panel_partner_manifest.final_evidence_stack.tsv"
SUMMARY_JSON = ROOT / "results/af3/af3_hpa_evidence.summary.json"


HPA_FIELDS = [
    "Gene",
    "Uniprot",
    "Gene description",
    "Protein class",
    "Biological process",
    "Molecular function",
    "Disease involvement",
    "Evidence",
    "RNA tissue specificity",
    "RNA tissue distribution",
    "Protein tissue specificity",
    "Protein tissue distribution",
    "Reliability (IF)",
    "Subcellular location",
    "Secretome location",
    "Secretome function",
    "Subcellular main location",
    "Subcellular additional location",
]


def clean(value: str | None) -> str:
    if value is None:
        return ""
    value = str(value).strip()
    return "" if value in {"NA", "N/A", "None", "nan"} else value


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_hpa() -> dict[str, dict[str, str]]:
    with zipfile.ZipFile(HPA_ZIP) as zf:
        with zf.open("proteinatlas.tsv") as raw:
            reader = csv.DictReader((line.decode("utf-8") for line in raw), delimiter="\t")
            by_gene: dict[str, dict[str, str]] = {}
            for row in reader:
                gene = clean(row.get("Gene")).upper()
                if not gene:
                    continue
                payload = {field: clean(row.get(field)) for field in HPA_FIELDS}
                payload["hpa_url"] = f"https://www.proteinatlas.org/{gene}"
                payload["hpa_source_file"] = str(HPA_ZIP)
                by_gene[gene] = payload
            return by_gene


def hpa_status(payload: dict[str, str] | None) -> str:
    if not payload:
        return "not_found_in_hpa_proteinatlas"
    if payload.get("Subcellular main location") or payload.get("Subcellular location"):
        return "hpa_subcellular_and_expression_available"
    if payload.get("RNA tissue specificity") or payload.get("Protein tissue specificity"):
        return "hpa_expression_available_subcellular_missing"
    return "hpa_row_found_minimal_annotation"


def prefixed(payload: dict[str, str] | None, prefix: str) -> dict[str, str]:
    cols = {
        "hpa_status": hpa_status(payload),
        "hpa_uniprot": "",
        "hpa_gene_description": "",
        "hpa_protein_class": "",
        "hpa_biological_process": "",
        "hpa_molecular_function": "",
        "hpa_disease_involvement": "",
        "hpa_evidence": "",
        "hpa_rna_tissue_specificity": "",
        "hpa_rna_tissue_distribution": "",
        "hpa_protein_tissue_specificity": "",
        "hpa_protein_tissue_distribution": "",
        "hpa_reliability_if": "",
        "hpa_subcellular_location": "",
        "hpa_subcellular_main_location": "",
        "hpa_subcellular_additional_location": "",
        "hpa_secretome_location": "",
        "hpa_secretome_function": "",
        "hpa_url": "",
    }
    if payload:
        cols.update(
            {
                "hpa_uniprot": payload.get("Uniprot", ""),
                "hpa_gene_description": payload.get("Gene description", ""),
                "hpa_protein_class": payload.get("Protein class", ""),
                "hpa_biological_process": payload.get("Biological process", ""),
                "hpa_molecular_function": payload.get("Molecular function", ""),
                "hpa_disease_involvement": payload.get("Disease involvement", ""),
                "hpa_evidence": payload.get("Evidence", ""),
                "hpa_rna_tissue_specificity": payload.get("RNA tissue specificity", ""),
                "hpa_rna_tissue_distribution": payload.get("RNA tissue distribution", ""),
                "hpa_protein_tissue_specificity": payload.get("Protein tissue specificity", ""),
                "hpa_protein_tissue_distribution": payload.get("Protein tissue distribution", ""),
                "hpa_reliability_if": payload.get("Reliability (IF)", ""),
                "hpa_subcellular_location": payload.get("Subcellular location", ""),
                "hpa_subcellular_main_location": payload.get("Subcellular main location", ""),
                "hpa_subcellular_additional_location": payload.get("Subcellular additional location", ""),
                "hpa_secretome_location": payload.get("Secretome location", ""),
                "hpa_secretome_function": payload.get("Secretome function", ""),
                "hpa_url": payload.get("hpa_url", ""),
            }
        )
    return {prefix + key: val for key, val in cols.items()}


def main() -> None:
    hpa = load_hpa()
    pair_rows = read_tsv(BASE_PAIR_TSV)
    genes = sorted(
        {
            r["target_gene"].upper()
            for r in pair_rows
            if r.get("target_gene")
        }
        | {
            r["partner_gene_or_entity"].upper()
            for r in pair_rows
            if r.get("partner_type") == "protein" and r.get("partner_gene_or_entity")
        }
    )

    gene_rows = []
    for gene in genes:
        payload = hpa.get(gene)
        row = {"gene": gene, **prefixed(payload, "")}
        gene_rows.append(row)

    gene_fields = ["gene"] + list(prefixed(None, "").keys())
    write_tsv(OUT_GENE_TSV, gene_rows, gene_fields)

    out_rows = []
    for row in pair_rows:
        merged = dict(row)
        target_payload = hpa.get(row.get("target_gene", "").upper())
        partner_payload = hpa.get(row.get("partner_gene_or_entity", "").upper()) if row.get("partner_type") == "protein" else None
        merged.update(prefixed(target_payload, "target_"))
        merged.update(prefixed(partner_payload, "partner_"))
        merged["hpa_pair_localization_note"] = ""
        if target_payload and partner_payload:
            tloc = target_payload.get("Subcellular main location") or target_payload.get("Subcellular location")
            ploc = partner_payload.get("Subcellular main location") or partner_payload.get("Subcellular location")
            if tloc and ploc:
                merged["hpa_pair_localization_note"] = "both_have_hpa_subcellular_annotations"
            elif tloc or ploc:
                merged["hpa_pair_localization_note"] = "one_partner_has_hpa_subcellular_annotation"
            else:
                merged["hpa_pair_localization_note"] = "hpa_expression_only_no_subcellular_pair_call"
        elif target_payload:
            merged["hpa_pair_localization_note"] = "target_hpa_available_partner_nonprotein_or_missing"
        else:
            merged["hpa_pair_localization_note"] = "target_missing_hpa"
        out_rows.append(merged)

    extra_fields = list(prefixed(None, "target_").keys()) + list(prefixed(None, "partner_").keys()) + ["hpa_pair_localization_note"]
    write_tsv(OUT_PAIR_TSV, out_rows, list(pair_rows[0].keys()) + extra_fields)

    summary = {
        "date": date.today().isoformat(),
        "hpa_source": str(HPA_ZIP),
        "hpa_gene_rows_loaded": len(hpa),
        "af3_evidence_genes": len(genes),
        "af3_genes_with_hpa_row": sum(g in hpa for g in genes),
        "af3_genes_with_hpa_subcellular": sum(bool((hpa.get(g) or {}).get("Subcellular main location") or (hpa.get(g) or {}).get("Subcellular location")) for g in genes),
        "target_hpa_status_counts": dict(Counter(r["target_hpa_status"] for r in out_rows)),
        "partner_hpa_status_counts": dict(Counter(r["partner_hpa_status"] for r in out_rows)),
        "outputs": {
            "hpa_evidence_by_gene": str(OUT_GENE_TSV),
            "final_pair_evidence_stack": str(OUT_PAIR_TSV),
        },
        "note": "HPA evidence is used for localization/expression context. It is not treated as direct protein-protein interaction evidence.",
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
