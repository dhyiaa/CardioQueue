#!/usr/bin/env python3
"""Build final joinable ClinGen features for the modeling skeleton."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


VALIDITY_RANK = {
    "Definitive": 6,
    "Strong": 5,
    "Moderate": 4,
    "Limited": 3,
    "Disputed": 2,
    "Refuted": 1,
    "No Known Disease Relationship": 0,
}

VARIANT_CLASS_RANK = {
    "Pathogenic": 5,
    "Likely Pathogenic": 4,
    "Uncertain Significance": 3,
    "Likely Benign": 2,
    "Benign": 1,
}

OUTPUT_FIELDS = [
    "variant_id",
    "primary_gene",
    "clingen_gene_validity_status",
    "clingen_gene_validity_max_classification",
    "clingen_gene_validity_max_rank",
    "clingen_gene_validity_has_definitive_strong",
    "clingen_gene_validity_has_moderate_plus",
    "clingen_gene_validity_classifications",
    "clingen_gene_validity_diseases",
    "clingen_gene_validity_expert_panels",
    "clingen_gene_validity_rows",
    "clingen_dosage_status",
    "clingen_haploinsufficiency_score",
    "clingen_haploinsufficiency_description",
    "clingen_triplosensitivity_score",
    "clingen_triplosensitivity_description",
    "clingen_dosage_date",
    "clingen_variant_evidence_status",
    "clingen_variant_evidence_match_count",
    "clingen_variant_evidence_classification",
    "clingen_variant_evidence_expert_panel",
    "clingen_variant_evidence_condition",
    "clingen_variant_evidence_mondo",
    "clingen_variant_evidence_met_codes",
    "clingen_variant_evidence_ca_id",
    "clingen_variant_evidence_clinvar_variation_id",
    "clingen_variant_evidence_pcer_doc_id",
    "clingen_variant_evidence_approved_date",
    "clingen_variant_evidence_published_date",
    "clingen_variant_evidence_hgvs_match",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def uniq_join(values: list[str]) -> str:
    seen: list[str] = []
    for value in values:
        value = (value or "").strip()
        if value and value not in seen:
            seen.append(value)
    return "|".join(seen)


def aggregate_validity(rows: list[dict[str, str]]) -> dict[str, dict[str, Any]]:
    by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_gene[row["input_gene"].upper()].append(row)

    out: dict[str, dict[str, Any]] = {}
    for gene, gene_rows in by_gene.items():
        ok_rows = [row for row in gene_rows if row.get("clingen_gene_validity_status") == "ok"]
        if not ok_rows:
            out[gene] = {
                "clingen_gene_validity_status": "not_found",
                "clingen_gene_validity_max_classification": "",
                "clingen_gene_validity_max_rank": "",
                "clingen_gene_validity_has_definitive_strong": "false",
                "clingen_gene_validity_has_moderate_plus": "false",
                "clingen_gene_validity_classifications": "",
                "clingen_gene_validity_diseases": "",
                "clingen_gene_validity_expert_panels": "",
                "clingen_gene_validity_rows": "0",
            }
            continue
        best = max(ok_rows, key=lambda row: VALIDITY_RANK.get(row.get("classification", ""), -1))
        best_rank = VALIDITY_RANK.get(best.get("classification", ""), "")
        out[gene] = {
            "clingen_gene_validity_status": "ok",
            "clingen_gene_validity_max_classification": best.get("classification", ""),
            "clingen_gene_validity_max_rank": best_rank,
            "clingen_gene_validity_has_definitive_strong": str(bool(isinstance(best_rank, int) and best_rank >= 5)).lower(),
            "clingen_gene_validity_has_moderate_plus": str(bool(isinstance(best_rank, int) and best_rank >= 4)).lower(),
            "clingen_gene_validity_classifications": uniq_join([row.get("classification", "") for row in ok_rows]),
            "clingen_gene_validity_diseases": uniq_join([row.get("disease_name", "") for row in ok_rows]),
            "clingen_gene_validity_expert_panels": uniq_join([row.get("expert_panel", "") for row in ok_rows]),
            "clingen_gene_validity_rows": str(len(ok_rows)),
        }
    return out


def parse_dosage_csv(path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8", newline="", errors="replace") as handle:
        reader = csv.reader(handle)
        header = None
        for raw in reader:
            if not raw:
                continue
            if raw[0] == "GENE SYMBOL":
                header = raw
                continue
            if header is None or raw[0].startswith("+") or len(raw) < 6:
                continue
            item = dict(zip(header, raw))
            gene = item.get("GENE SYMBOL", "").upper()
            if gene:
                rows[gene] = {
                    "clingen_dosage_status": "ok",
                    "clingen_haploinsufficiency_score": item.get("HAPLOINSUFFICIENCY", ""),
                    "clingen_haploinsufficiency_description": item.get("HAPLOINSUFFICIENCY", ""),
                    "clingen_triplosensitivity_score": item.get("TRIPLOSENSITIVITY", ""),
                    "clingen_triplosensitivity_description": item.get("TRIPLOSENSITIVITY", ""),
                    "clingen_dosage_date": item.get("DATE", ""),
                }
    return rows


def parse_dosage_grch38(path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    header: list[str] | None = None
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("#Gene Symbol\t"):
                header = line[1:].split("\t")
                continue
            if line.startswith("#") or not line or header is None:
                continue
            parts = line.split("\t")
            item = dict(zip(header, parts))
            gene = item.get("Gene Symbol", "").upper()
            if gene:
                rows[gene] = {
                    "clingen_dosage_status": "ok",
                    "clingen_haploinsufficiency_score": item.get("Haploinsufficiency Score", ""),
                    "clingen_haploinsufficiency_description": item.get("Haploinsufficiency Description", ""),
                    "clingen_triplosensitivity_score": item.get("Triplosensitivity Score", ""),
                    "clingen_triplosensitivity_description": item.get("Triplosensitivity Description", ""),
                    "clingen_dosage_date": item.get("Date Last Evaluated", "") or item.get("Date Last Evaluated ", ""),
                }
    return rows


def genomic_keys_from_hgvs(hgvs: str) -> set[str]:
    keys: set[str] = set()
    for token in re.split(r"[|;,\s]+", hgvs or ""):
        token = token.strip()
        match = re.search(r"NC_0000(?P<num>\d{2})\.\d+:g\.(?P<pos>\d+)(?P<ref>[ACGT])>(?P<alt>[ACGT])", token)
        if match:
            chrom_num = str(int(match.group("num")))
            keys.add(f"{chrom_num}-{match.group('pos')}-{match.group('ref')}-{match.group('alt')}")
            continue
        match = re.search(r"NC_000023\.\d+:g\.(?P<pos>\d+)(?P<ref>[ACGT])>(?P<alt>[ACGT])", token)
        if match:
            keys.add(f"X-{match.group('pos')}-{match.group('ref')}-{match.group('alt')}")
            continue
        match = re.search(r"NC_000024\.\d+:g\.(?P<pos>\d+)(?P<ref>[ACGT])>(?P<alt>[ACGT])", token)
        if match:
            keys.add(f"Y-{match.group('pos')}-{match.group('ref')}-{match.group('alt')}")
            continue
        match = re.search(r"CM\d+\.\d+:g\.(?P<pos>\d+)(?P<ref>[ACGT])>(?P<alt>[ACGT])", token)
        # CM accessions do not encode chromosome safely enough for a key here.
    return keys


def build_variant_index(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    by_variant: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row.get("clingen_variant_pathogenicity_status") != "ok":
            continue
        if str(row.get("retracted", "")).lower() == "true":
            continue
        for key in genomic_keys_from_hgvs(row.get("hgvs", "")):
            by_variant[key].append(row)
    return by_variant


def choose_variant_evidence(matches: list[dict[str, str]]) -> dict[str, str]:
    if not matches:
        return {
            "clingen_variant_evidence_status": "not_found",
            "clingen_variant_evidence_match_count": "0",
        }
    best = max(matches, key=lambda row: VARIANT_CLASS_RANK.get(row.get("classification", ""), 0))
    return {
        "clingen_variant_evidence_status": "ok",
        "clingen_variant_evidence_match_count": str(len(matches)),
        "clingen_variant_evidence_classification": best.get("classification", ""),
        "clingen_variant_evidence_expert_panel": best.get("expert_panel", ""),
        "clingen_variant_evidence_condition": best.get("condition", ""),
        "clingen_variant_evidence_mondo": best.get("mondo_id", ""),
        "clingen_variant_evidence_met_codes": best.get("met_codes", ""),
        "clingen_variant_evidence_ca_id": best.get("ca_id", ""),
        "clingen_variant_evidence_clinvar_variation_id": best.get("clinvar_variation_id", ""),
        "clingen_variant_evidence_pcer_doc_id": best.get("pcer_doc_id", ""),
        "clingen_variant_evidence_approved_date": best.get("approved_date", ""),
        "clingen_variant_evidence_published_date": best.get("published_date", ""),
        "clingen_variant_evidence_hgvs_match": uniq_join([m.get("preferred_variant_title", "") for m in matches]),
    }


def build(args: argparse.Namespace) -> dict[str, Any]:
    skeleton = read_tsv(Path(args.modeling_skeleton))
    validity = aggregate_validity(read_tsv(Path(args.gene_validity)))
    dosage = parse_dosage_grch38(Path(args.dosage_grch38))
    if not dosage:
        dosage = parse_dosage_csv(Path(args.dosage_csv))
    variant_index = build_variant_index(read_tsv(Path(args.variant_pathogenicity)))

    output_rows: list[dict[str, Any]] = []
    counters: Counter[str] = Counter()
    for row in skeleton:
        gene = (row.get("primary_gene") or "").upper()
        out = {field: "" for field in OUTPUT_FIELDS}
        out["variant_id"] = row.get("variant_id", "")
        out["primary_gene"] = gene

        validity_payload = validity.get(gene)
        if validity_payload:
            out.update(validity_payload)
        else:
            out["clingen_gene_validity_status"] = "not_found"
            out["clingen_gene_validity_has_definitive_strong"] = "false"
            out["clingen_gene_validity_has_moderate_plus"] = "false"
            out["clingen_gene_validity_rows"] = "0"

        dosage_payload = dosage.get(gene)
        if dosage_payload:
            out.update(dosage_payload)
        else:
            out["clingen_dosage_status"] = "not_found"

        variant_payload = choose_variant_evidence(variant_index.get(row.get("variant_id", ""), []))
        out.update(variant_payload)

        counters[f"gene_validity_{out['clingen_gene_validity_status']}"] += 1
        counters[f"dosage_{out['clingen_dosage_status']}"] += 1
        counters[f"variant_evidence_{out['clingen_variant_evidence_status']}"] += 1
        output_rows.append(out)

    write_tsv(Path(args.output), output_rows, OUTPUT_FIELDS)
    summary = {
        "modeling_rows": len(skeleton),
        "output": args.output,
        "gene_validity_input_rows": sum(1 for _ in read_tsv(Path(args.gene_validity))),
        "dosage_gene_rows": len(dosage),
        "variant_pathogenicity_input_rows": sum(1 for _ in read_tsv(Path(args.variant_pathogenicity))),
        "variant_pathogenicity_genomic_keys": len(variant_index),
        "counts": dict(counters),
    }
    Path(args.summary).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modeling-skeleton", default="datasets/modeling/interim/final_modeling_table_skeleton.tsv")
    parser.add_argument("--gene-validity", default="datasets/feature_sources/clingen/interim/gene_validity.master_registry_panel.tsv")
    parser.add_argument("--dosage-csv", default="datasets/feature_sources/clingen/raw/dosage_sensitivity/gene_dosage.csv")
    parser.add_argument("--dosage-grch38", default="datasets/feature_sources/clingen/raw/dosage_sensitivity/gene_dosage_GRCh38.tsv")
    parser.add_argument("--variant-pathogenicity", default="datasets/feature_sources/clingen/interim/variant_pathogenicity.master_registry_panel.tsv")
    parser.add_argument("--output", default="datasets/feature_sources/clingen/interim/clingen_modeling_features.tsv")
    parser.add_argument("--summary", default="datasets/feature_sources/clingen/interim/clingen_modeling_features.summary.json")
    args = parser.parse_args()
    print(json.dumps(build(args), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
