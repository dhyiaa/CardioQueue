#!/usr/bin/env python3
"""Targeted AlphaFold DB and UniProt annotation helpers.

This script deliberately avoids downloading the full human AlphaFold proteome.
It maps gene symbols to reviewed human UniProt accessions, then records
AlphaFold DB model metadata and confidence links. Residue-level features can be
added from AlphaFold confidence JSON when a UniProt accession and protein
position are available.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any


UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
ALPHAFOLD_PREDICTION_URL = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
DEFAULT_TIMEOUT_SECONDS = 60

GENE_PANEL_FIELDS = [
    "input_gene",
    "uniprot_status",
    "uniprot_missing_reason",
    "uniprot_accession",
    "uniprot_gene_names",
    "uniprot_protein_name",
    "uniprot_length",
    "alphafold_status",
    "alphafold_missing_reason",
    "alphafold_entry_id",
    "alphafold_model_entity_id",
    "alphafold_latest_version",
    "alphafold_model_created_date",
    "alphafold_sequence_version_date",
    "alphafold_sequence_start",
    "alphafold_sequence_end",
    "alphafold_global_plddt",
    "alphafold_fraction_plddt_very_low",
    "alphafold_fraction_plddt_low",
    "alphafold_fraction_plddt_confident",
    "alphafold_fraction_plddt_very_high",
    "alphafold_pdb_url",
    "alphafold_cif_url",
    "alphafold_confidence_json_url",
    "alphafold_pae_json_url",
    "alphafold_am_hg38_url",
]

RESIDUE_FIELDS = [
    "input_id",
    "uniprot_accession",
    "protein_position",
    "alphafold_residue_status",
    "alphafold_residue_missing_reason",
    "alphafold_residue_plddt",
]


def fetch_text(url: str, timeout: int = DEFAULT_TIMEOUT_SECONDS) -> str:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "catboost-largedata-synth/0.1"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8")


def fetch_json(url: str, timeout: int = DEFAULT_TIMEOUT_SECONDS) -> Any:
    return json.loads(fetch_text(url, timeout=timeout))


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


def uniprot_search_url(gene: str) -> str:
    query = f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true"
    params = {
        "query": query,
        "fields": "accession,gene_names,protein_name,length",
        "format": "tsv",
        "size": "10",
    }
    return f"{UNIPROT_SEARCH_URL}?{urllib.parse.urlencode(params)}"


def lookup_uniprot(gene: str, timeout: int) -> tuple[str, dict[str, str]]:
    text = fetch_text(uniprot_search_url(gene), timeout=timeout)
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    if not rows:
        return "not_found", {}
    exact_rows = [
        row for row in rows
        if gene in {item.upper() for item in row.get("Gene Names", "").split()}
    ]
    selected = exact_rows[0] if exact_rows else rows[0]
    return "ok", {
        "uniprot_accession": selected.get("Entry", ""),
        "uniprot_gene_names": selected.get("Gene Names", ""),
        "uniprot_protein_name": selected.get("Protein names", ""),
        "uniprot_length": selected.get("Length", ""),
    }


def select_alphafold_entry(entries: list[dict[str, Any]], accession: str) -> dict[str, Any] | None:
    canonical = [
        entry for entry in entries
        if entry.get("uniprotAccession") == accession
        and entry.get("isUniProtReviewed")
        and entry.get("isUniProtReferenceProteome")
        and str(entry.get("entryId", "")).startswith(f"AF-{accession}-F")
    ]
    if canonical:
        return canonical[0]
    accession_entries = [entry for entry in entries if entry.get("uniprotAccession") == accession]
    return accession_entries[0] if accession_entries else None


def alphafold_metadata(accession: str, timeout: int) -> tuple[str, dict[str, Any]]:
    try:
        entries = fetch_json(ALPHAFOLD_PREDICTION_URL.format(accession=accession), timeout=timeout)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return "not_found", {}
        raise
    if not isinstance(entries, list):
        raise RuntimeError(f"AlphaFold response for {accession} was not a list")
    if not entries:
        return "not_found", {}
    entry = select_alphafold_entry(entries, accession)
    if entry is None:
        return "not_found", {}
    return "ok", {
        "alphafold_entry_id": entry.get("entryId", ""),
        "alphafold_model_entity_id": entry.get("modelEntityId", ""),
        "alphafold_latest_version": entry.get("latestVersion", ""),
        "alphafold_model_created_date": entry.get("modelCreatedDate", ""),
        "alphafold_sequence_version_date": entry.get("sequenceVersionDate", ""),
        "alphafold_sequence_start": entry.get("sequenceStart", ""),
        "alphafold_sequence_end": entry.get("sequenceEnd", ""),
        "alphafold_global_plddt": entry.get("globalMetricValue", ""),
        "alphafold_fraction_plddt_very_low": entry.get("fractionPlddtVeryLow", ""),
        "alphafold_fraction_plddt_low": entry.get("fractionPlddtLow", ""),
        "alphafold_fraction_plddt_confident": entry.get("fractionPlddtConfident", ""),
        "alphafold_fraction_plddt_very_high": entry.get("fractionPlddtVeryHigh", ""),
        "alphafold_pdb_url": entry.get("pdbUrl", ""),
        "alphafold_cif_url": entry.get("cifUrl", ""),
        "alphafold_confidence_json_url": entry.get("plddtDocUrl", ""),
        "alphafold_pae_json_url": entry.get("paeDocUrl", ""),
        "alphafold_am_hg38_url": entry.get("amAnnotationsHg38Url", ""),
    }


def empty_gene_row(gene: str, status: str, reason: str) -> dict[str, Any]:
    row = {field: "" for field in GENE_PANEL_FIELDS}
    row["input_gene"] = gene
    row["uniprot_status"] = status
    row["uniprot_missing_reason"] = reason
    if status != "ok":
        row["alphafold_status"] = "not_queryable"
        row["alphafold_missing_reason"] = "UniProt accession unavailable."
    return row


def cmd_doctor(args: argparse.Namespace) -> int:
    started = time.time()
    status, uniprot = lookup_uniprot(args.gene, timeout=args.timeout)
    if status != "ok":
        raise RuntimeError(f"UniProt lookup failed for {args.gene}: {status}")
    af_status, af = alphafold_metadata(uniprot["uniprot_accession"], timeout=args.timeout)
    if af_status != "ok":
        raise RuntimeError(f"AlphaFold lookup failed for {uniprot['uniprot_accession']}: {af_status}")
    print(json.dumps({
        "status": "ok",
        "gene": args.gene,
        "uniprot_accession": uniprot["uniprot_accession"],
        "alphafold_entry_id": af["alphafold_entry_id"],
        "alphafold_latest_version": af["alphafold_latest_version"],
        "alphafold_global_plddt": af["alphafold_global_plddt"],
        "elapsed_seconds": round(time.time() - started, 3),
    }, indent=2))
    return 0


def cmd_gene_panel(args: argparse.Namespace) -> int:
    genes = load_gene_panel(Path(args.gene_panel))
    output_rows: list[dict[str, Any]] = []
    raw_records: dict[str, Any] = {}

    for gene in genes:
        row = {field: "" for field in GENE_PANEL_FIELDS}
        row["input_gene"] = gene
        try:
            uniprot_status, uniprot = lookup_uniprot(gene, timeout=args.timeout)
        except Exception as exc:
            row = empty_gene_row(gene, "error", f"UniProt lookup failed: {exc}")
            output_rows.append(row)
            continue

        if uniprot_status != "ok":
            output_rows.append(empty_gene_row(gene, uniprot_status, "No reviewed human UniProt entry found."))
            continue

        row["uniprot_status"] = "ok"
        row.update(uniprot)
        try:
            af_status, af = alphafold_metadata(uniprot["uniprot_accession"], timeout=args.timeout)
        except Exception as exc:
            row["alphafold_status"] = "error"
            row["alphafold_missing_reason"] = f"AlphaFold lookup failed: {exc}"
            output_rows.append(row)
            continue

        row["alphafold_status"] = af_status
        if af_status == "ok":
            row.update(af)
            raw_records[gene] = {"uniprot": uniprot, "alphafold": af}
        else:
            row["alphafold_missing_reason"] = "No AlphaFold DB model found for UniProt accession."
        output_rows.append(row)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=GENE_PANEL_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    if args.raw_json:
        raw_path = Path(args.raw_json)
        raw_path.parent.mkdir(parents=True, exist_ok=True)
        raw_path.write_text(json.dumps(raw_records, indent=2))

    status_counts: dict[str, int] = {}
    for row in output_rows:
        status = str(row["alphafold_status"])
        status_counts[status] = status_counts.get(status, 0) + 1
    print(json.dumps({
        "status": "ok",
        "input_genes": len(genes),
        "output_rows": len(output_rows),
        "alphafold_status_counts": status_counts,
        "output": str(output_path),
    }, indent=2))
    return 0


def cmd_validate_gene_panel(args: argparse.Namespace) -> int:
    path = Path(args.input)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [field for field in GENE_PANEL_FIELDS if field not in (reader.fieldnames or [])]
        if missing:
            raise RuntimeError(f"missing required columns: {missing}")
        rows = list(reader)
    if not rows:
        raise RuntimeError("AlphaFold gene-panel table has no rows")
    bad_ok = [
        row["input_gene"] for row in rows
        if row["alphafold_status"] == "ok" and not row["alphafold_global_plddt"]
    ]
    if bad_ok:
        raise RuntimeError(f"ok AlphaFold rows missing global pLDDT: {bad_ok[:10]}")
    bad_missing = [
        row["input_gene"] for row in rows
        if row["alphafold_status"] != "ok" and not row["alphafold_missing_reason"]
    ]
    if bad_missing:
        raise RuntimeError(f"non-ok AlphaFold rows missing reason: {bad_missing[:10]}")
    counts: dict[str, int] = {}
    for row in rows:
        status = row["alphafold_status"]
        counts[status] = counts.get(status, 0) + 1
    print(json.dumps({"rows": len(rows), "alphafold_status_counts": counts}, indent=2))
    return 0


def cmd_residue_features(args: argparse.Namespace) -> int:
    input_path = Path(args.input)
    with input_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"input_id", "uniprot_accession", "protein_position"}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise RuntimeError(f"missing required input columns: {missing}")
        input_rows = list(reader)

    confidence_cache: dict[str, dict[str, Any]] = {}
    output_rows: list[dict[str, Any]] = []

    for input_row in input_rows:
        output = {field: "" for field in RESIDUE_FIELDS}
        output["input_id"] = input_row["input_id"]
        output["uniprot_accession"] = input_row["uniprot_accession"]
        output["protein_position"] = input_row["protein_position"]
        try:
            position = int(input_row["protein_position"])
        except ValueError:
            output["alphafold_residue_status"] = "not_queryable"
            output["alphafold_residue_missing_reason"] = "protein_position is not an integer"
            output_rows.append(output)
            continue

        accession = input_row["uniprot_accession"]
        if accession not in confidence_cache:
            af_status, af = alphafold_metadata(accession, timeout=args.timeout)
            if af_status != "ok":
                confidence_cache[accession] = {"status": af_status, "reason": "No AlphaFold metadata found"}
            else:
                confidence_cache[accession] = fetch_json(af["alphafold_confidence_json_url"], timeout=args.timeout)

        confidence = confidence_cache[accession]
        if "status" in confidence:
            output["alphafold_residue_status"] = confidence["status"]
            output["alphafold_residue_missing_reason"] = confidence["reason"]
            output_rows.append(output)
            continue

        residue_numbers = confidence.get("residueNumber", [])
        plddt_values = confidence.get("confidenceScore", [])
        if position not in residue_numbers:
            output["alphafold_residue_status"] = "not_found"
            output["alphafold_residue_missing_reason"] = "position outside AlphaFold confidence JSON residue range"
        else:
            idx = residue_numbers.index(position)
            output["alphafold_residue_status"] = "ok"
            output["alphafold_residue_plddt"] = plddt_values[idx]
        output_rows.append(output)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESIDUE_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    print(json.dumps({"status": "ok", "rows": len(output_rows), "output": str(output_path)}, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    doctor = subparsers.add_parser("doctor")
    doctor.add_argument("--gene", default="MYH7")
    doctor.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT_SECONDS)
    doctor.set_defaults(func=cmd_doctor)

    gene_panel = subparsers.add_parser("gene-panel")
    gene_panel.add_argument("--gene-panel", required=True)
    gene_panel.add_argument("--output", required=True)
    gene_panel.add_argument("--raw-json")
    gene_panel.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT_SECONDS)
    gene_panel.set_defaults(func=cmd_gene_panel)

    validate = subparsers.add_parser("validate-gene-panel")
    validate.add_argument("--input", required=True)
    validate.set_defaults(func=cmd_validate_gene_panel)

    residue = subparsers.add_parser("residue-features")
    residue.add_argument("--input", required=True)
    residue.add_argument("--output", required=True)
    residue.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT_SECONDS)
    residue.set_defaults(func=cmd_residue_features)
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
