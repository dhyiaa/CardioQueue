#!/usr/bin/env python3
"""ClinGen targeted annotation helpers.

This script supports ClinGen gene-disease validity curations, ClinGen Evidence
Repository variant pathogenicity summaries, and dosage-sensitivity downloads.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any


VALIDITY_URL = "https://search.clinicalgenome.org/api/validity"
GENE_LOOK_URL = "https://search.clinicalgenome.org/api/genes/look/{query}"
EREPO_CLASSIFICATIONS_URL = "https://erepo.clinicalgenome.org/evrepo/api/summary/classifications"
DOSAGE_GENE_CSV_URL = "https://search.clinicalgenome.org/kb/gene-dosage/download"
DOSAGE_GRCH38_TSV_URL = "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"
DEFAULT_TIMEOUT_SECONDS = 60


VALIDITY_FIELDS = [
    "input_gene",
    "clingen_gene_validity_status",
    "clingen_gene_validity_missing_reason",
    "symbol",
    "hgnc_id",
    "disease_name",
    "mondo",
    "moi",
    "classification",
    "expert_panel",
    "affiliate_id",
    "sop",
    "perm_id",
    "report_id",
    "released",
    "date",
    "animal_model_only",
]

VARIANT_PATHOGENICITY_FIELDS = [
    "input_gene",
    "clingen_variant_pathogenicity_status",
    "clingen_variant_pathogenicity_missing_reason",
    "gene",
    "ca_id",
    "clinvar_variation_id",
    "classification",
    "expert_panel",
    "condition",
    "mondo_id",
    "moi",
    "preferred_variant_title",
    "approved_date",
    "published_date",
    "retracted",
    "met_codes",
    "hgvs",
    "summary",
    "pcer_doc_id",
]


def fetch_json(url: str, timeout: int = DEFAULT_TIMEOUT_SECONDS) -> Any:
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "catboost-largedata-synth/0.1",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = response.read()
    except urllib.error.URLError as exc:
        raise RuntimeError(f"request failed for {url}: {exc}") from exc

    try:
        return json.loads(payload.decode("utf-8"))
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"response was not valid JSON for {url}: {exc}") from exc


def fetch_bytes(url: str, timeout: int = DEFAULT_TIMEOUT_SECONDS) -> bytes:
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "*/*",
            "User-Agent": "catboost-largedata-synth/0.1",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return response.read()
    except urllib.error.URLError as exc:
        raise RuntimeError(f"request failed for {url}: {exc}") from exc


def load_gene_panel(path: Path) -> list[str]:
    genes: list[str] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            genes.append(line.split()[0].upper())
    if not genes:
        raise ValueError(f"no genes found in {path}")
    return genes


def fetch_validity_rows(timeout: int) -> list[dict[str, Any]]:
    payload = fetch_json(VALIDITY_URL, timeout=timeout)
    if not isinstance(payload, dict):
        raise RuntimeError("ClinGen validity response was not a JSON object")
    rows = payload.get("rows")
    if not isinstance(rows, list):
        raise RuntimeError("ClinGen validity response missing list field: rows")
    for row in rows:
        if not isinstance(row, dict):
            raise RuntimeError("ClinGen validity response contains a non-object row")
    return rows


def normalize_validity_row(input_gene: str, row: dict[str, Any]) -> dict[str, Any]:
    return {
        "input_gene": input_gene,
        "clingen_gene_validity_status": "ok",
        "clingen_gene_validity_missing_reason": "",
        "symbol": row.get("symbol", ""),
        "hgnc_id": row.get("hgnc_id", ""),
        "disease_name": row.get("disease_name", ""),
        "mondo": row.get("mondo", ""),
        "moi": row.get("moi", ""),
        "classification": row.get("classification", ""),
        "expert_panel": row.get("ep", ""),
        "affiliate_id": row.get("affiliate_id", ""),
        "sop": row.get("sop", ""),
        "perm_id": row.get("perm_id", ""),
        "report_id": row.get("report_id", ""),
        "released": row.get("released", ""),
        "date": row.get("date", ""),
        "animal_model_only": row.get("animal_model_only", ""),
    }


def not_found_row(input_gene: str) -> dict[str, Any]:
    return {
        field: ""
        for field in VALIDITY_FIELDS
    } | {
        "input_gene": input_gene,
        "clingen_gene_validity_status": "not_found",
        "clingen_gene_validity_missing_reason": "No ClinGen gene-disease validity row found for the input gene symbol.",
    }


def cmd_doctor(args: argparse.Namespace) -> int:
    started = time.time()
    rows = fetch_validity_rows(timeout=args.timeout)
    gene_payload = fetch_json(GENE_LOOK_URL.format(query="MYH7"), timeout=args.timeout)
    if not isinstance(gene_payload, list):
        raise RuntimeError("ClinGen gene lookup response was not a JSON list")
    erepo_payload = fetch_json(f"{EREPO_CLASSIFICATIONS_URL}?gene=MYH7", timeout=args.timeout)
    if not isinstance(erepo_payload, dict) or not isinstance(erepo_payload.get("data"), list):
        raise RuntimeError("ClinGen ERepo classifications response missing data list")
    myh7_rows = [row for row in rows if row.get("symbol") == "MYH7"]
    print(
        json.dumps(
            {
                "status": "ok",
                "validity_rows": len(rows),
                "myh7_validity_rows": len(myh7_rows),
                "gene_lookup_hits": len(gene_payload),
                "myh7_variant_pathogenicity_rows": len(erepo_payload["data"]),
                "elapsed_seconds": round(time.time() - started, 3),
            },
            indent=2,
        )
    )
    return 0


def normalize_variant_pathogenicity_row(input_gene: str, row: dict[str, Any]) -> dict[str, Any]:
    return {
        "input_gene": input_gene,
        "clingen_variant_pathogenicity_status": "ok",
        "clingen_variant_pathogenicity_missing_reason": "",
        "gene": row.get("gene", ""),
        "ca_id": row.get("caId", ""),
        "clinvar_variation_id": row.get("cvId", ""),
        "classification": row.get("classification", ""),
        "expert_panel": row.get("ep", ""),
        "condition": row.get("condition", ""),
        "mondo_id": row.get("mondoId", ""),
        "moi": row.get("moi", ""),
        "preferred_variant_title": row.get("preferredVarTitle", ""),
        "approved_date": row.get("approvedDate", ""),
        "published_date": row.get("publishedDate", ""),
        "retracted": row.get("retracted", ""),
        "met_codes": "|".join(row.get("metCodes") or []),
        "hgvs": "|".join(row.get("hgvs") or []),
        "summary": row.get("summaryDesc", ""),
        "pcer_doc_id": row.get("PCERDocID", ""),
    }


def variant_pathogenicity_not_found(input_gene: str) -> dict[str, Any]:
    return {
        field: ""
        for field in VARIANT_PATHOGENICITY_FIELDS
    } | {
        "input_gene": input_gene,
        "clingen_variant_pathogenicity_status": "not_found",
        "clingen_variant_pathogenicity_missing_reason": "No ClinGen Evidence Repository variant pathogenicity rows found for the input gene.",
    }


def fetch_variant_pathogenicity_for_gene(gene: str, timeout: int) -> list[dict[str, Any]]:
    try:
        payload = fetch_json(f"{EREPO_CLASSIFICATIONS_URL}?gene={gene}", timeout=timeout)
    except RuntimeError as exc:
        if "HTTP Error 404" in str(exc):
            return []
        raise
    if not isinstance(payload, dict):
        raise RuntimeError(f"ERepo response for {gene} was not a JSON object")
    rows = payload.get("data")
    if not isinstance(rows, list):
        raise RuntimeError(f"ERepo response for {gene} missing data list")
    return rows


def cmd_variant_pathogenicity(args: argparse.Namespace) -> int:
    genes = load_gene_panel(Path(args.gene_panel))
    output_rows: list[dict[str, Any]] = []
    raw: dict[str, Any] = {"source_url": EREPO_CLASSIFICATIONS_URL, "genes": {}}
    for gene in genes:
        rows = fetch_variant_pathogenicity_for_gene(gene, timeout=args.timeout)
        raw["genes"][gene] = rows
        if rows:
            output_rows.extend(normalize_variant_pathogenicity_row(gene, row) for row in rows)
        else:
            output_rows.append(variant_pathogenicity_not_found(gene))
        if args.sleep:
            time.sleep(args.sleep)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=VARIANT_PATHOGENICITY_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    if args.raw_json:
        raw_path = Path(args.raw_json)
        raw_path.parent.mkdir(parents=True, exist_ok=True)
        raw_path.write_text(json.dumps(raw, indent=2))

    ok_genes = {
        row["input_gene"]
        for row in output_rows
        if row["clingen_variant_pathogenicity_status"] == "ok"
    }
    print(
        json.dumps(
            {
                "status": "ok",
                "input_genes": len(genes),
                "output_rows": len(output_rows),
                "genes_with_variant_curations": len(ok_genes),
                "genes_without_variant_curations": len(genes) - len(ok_genes),
                "output": str(output_path),
            },
            indent=2,
        )
    )
    return 0


def cmd_download_dosage(args: argparse.Namespace) -> int:
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    files = {
        "gene_dosage.csv": DOSAGE_GENE_CSV_URL,
        "gene_dosage_GRCh38.tsv": DOSAGE_GRCH38_TSV_URL,
    }
    manifest = []
    for filename, url in files.items():
        payload = fetch_bytes(url, timeout=args.timeout)
        path = output_dir / filename
        path.write_bytes(payload)
        manifest.append({"filename": filename, "url": url, "bytes": len(payload)})
    (output_dir / "download_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps({"status": "ok", "files": manifest}, indent=2))
    return 0


def cmd_gene_validity(args: argparse.Namespace) -> int:
    genes = load_gene_panel(Path(args.gene_panel))
    rows = fetch_validity_rows(timeout=args.timeout)
    rows_by_symbol: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        symbol = str(row.get("symbol", "")).upper()
        if symbol:
            rows_by_symbol.setdefault(symbol, []).append(row)

    output_rows: list[dict[str, Any]] = []
    for gene in genes:
        matches = rows_by_symbol.get(gene, [])
        if matches:
            output_rows.extend(normalize_validity_row(gene, row) for row in matches)
        else:
            output_rows.append(not_found_row(gene))

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=VALIDITY_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    if args.raw_json:
        raw_path = Path(args.raw_json)
        raw_path.parent.mkdir(parents=True, exist_ok=True)
        raw_path.write_text(json.dumps({"source_url": VALIDITY_URL, "rows": rows}, indent=2))

    print(
        json.dumps(
            {
                "status": "ok",
                "input_genes": len(genes),
                "output_rows": len(output_rows),
                "genes_with_curations": sum(1 for gene in genes if rows_by_symbol.get(gene)),
                "genes_without_curations": sum(1 for gene in genes if not rows_by_symbol.get(gene)),
                "output": str(output_path),
            },
            indent=2,
        )
    )
    return 0


def cmd_validate(args: argparse.Namespace) -> int:
    path = Path(args.input)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [field for field in VALIDITY_FIELDS if field not in (reader.fieldnames or [])]
        if missing:
            raise RuntimeError(f"missing required columns: {missing}")
        rows = list(reader)

    if not rows:
        raise RuntimeError("ClinGen gene validity table has no rows")

    bad_ok = [
        row["input_gene"]
        for row in rows
        if row["clingen_gene_validity_status"] == "ok"
        and not row["classification"]
    ]
    if bad_ok:
        raise RuntimeError(f"ok ClinGen rows missing classification: {bad_ok[:10]}")

    bad_missing = [
        row["input_gene"]
        for row in rows
        if row["clingen_gene_validity_status"] != "ok"
        and not row["clingen_gene_validity_missing_reason"]
    ]
    if bad_missing:
        raise RuntimeError(f"non-ok ClinGen rows missing reason: {bad_missing[:10]}")

    status_counts: dict[str, int] = {}
    for row in rows:
        status = row["clingen_gene_validity_status"]
        status_counts[status] = status_counts.get(status, 0) + 1
    print(json.dumps({"rows": len(rows), "status_counts": status_counts}, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    doctor = subparsers.add_parser("doctor", help="Check ClinGen API shape")
    doctor.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT_SECONDS)
    doctor.set_defaults(func=cmd_doctor)

    validity = subparsers.add_parser(
        "gene-validity",
        help="Annotate a gene panel with ClinGen gene-disease validity rows",
    )
    validity.add_argument("--gene-panel", required=True)
    validity.add_argument("--output", required=True)
    validity.add_argument("--raw-json")
    validity.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT_SECONDS)
    validity.set_defaults(func=cmd_gene_validity)

    variant_pathogenicity = subparsers.add_parser(
        "variant-pathogenicity",
        help="Annotate a gene panel with ClinGen Evidence Repository variant pathogenicity summaries",
    )
    variant_pathogenicity.add_argument("--gene-panel", required=True)
    variant_pathogenicity.add_argument("--output", required=True)
    variant_pathogenicity.add_argument("--raw-json")
    variant_pathogenicity.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT_SECONDS)
    variant_pathogenicity.add_argument("--sleep", type=float, default=0.0)
    variant_pathogenicity.set_defaults(func=cmd_variant_pathogenicity)

    dosage = subparsers.add_parser("download-dosage", help="Download ClinGen dosage-sensitivity files")
    dosage.add_argument("--output-dir", required=True)
    dosage.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT_SECONDS)
    dosage.set_defaults(func=cmd_download_dosage)

    validate = subparsers.add_parser("validate", help="Validate ClinGen output")
    validate.add_argument("--input", required=True)
    validate.set_defaults(func=cmd_validate)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        return args.func(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
