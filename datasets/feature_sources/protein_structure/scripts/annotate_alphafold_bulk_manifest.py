#!/usr/bin/env python3
"""Annotate UniProt accessions using a local AlphaFold bulk tar manifest."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


NAME_RE = re.compile(r"^AF-(?P<accession>.+)-(?P<fragment>F\d+)-model_v(?P<version>\d+)\.(?P<kind>pdb|cif)\.gz$")

OUTPUT_FIELDS = [
    "input_gene",
    "uniprot_accession",
    "alphafold_bulk_status",
    "alphafold_bulk_missing_reason",
    "alphafold_bulk_latest_version",
    "alphafold_bulk_fragment_count",
    "alphafold_bulk_has_pdb",
    "alphafold_bulk_has_cif",
    "alphafold_bulk_fragments",
]


def load_manifest(path: Path) -> dict[str, dict[str, object]]:
    index: dict[str, dict[str, object]] = {}
    with path.open() as handle:
        for raw_line in handle:
            name = raw_line.strip()
            match = NAME_RE.match(name)
            if not match:
                continue
            accession = match.group("accession")
            fragment = match.group("fragment")
            version = int(match.group("version"))
            kind = match.group("kind")
            record = index.setdefault(
                accession,
                {
                    "versions": set(),
                    "fragments": set(),
                    "kinds": set(),
                },
            )
            record["versions"].add(version)
            record["fragments"].add(fragment)
            record["kinds"].add(kind)
    return index


def normalize_record(record: dict[str, object]) -> dict[str, object]:
    versions = sorted(record["versions"])
    fragments = sorted(record["fragments"], key=lambda item: int(str(item)[1:]))
    kinds = set(record["kinds"])
    return {
        "latest_version": versions[-1] if versions else "",
        "fragment_count": len(fragments),
        "has_pdb": "pdb" in kinds,
        "has_cif": "cif" in kinds,
        "fragments": ",".join(fragments),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="TSV with input_gene and uniprot_accession")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    manifest_index = load_manifest(Path(args.manifest))
    output_rows = []
    with Path(args.input).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"input_gene", "uniprot_accession"}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise RuntimeError(f"input missing required columns: {missing}")
        for row in reader:
            output = {field: "" for field in OUTPUT_FIELDS}
            output["input_gene"] = row["input_gene"]
            output["uniprot_accession"] = row["uniprot_accession"]
            accession = row["uniprot_accession"]
            if not accession:
                output["alphafold_bulk_status"] = "not_queryable"
                output["alphafold_bulk_missing_reason"] = "uniprot_accession missing"
            elif accession not in manifest_index:
                output["alphafold_bulk_status"] = "not_found"
                output["alphafold_bulk_missing_reason"] = "UniProt accession not present in AlphaFold bulk manifest"
            else:
                normalized = normalize_record(manifest_index[accession])
                output["alphafold_bulk_status"] = "ok"
                output["alphafold_bulk_latest_version"] = normalized["latest_version"]
                output["alphafold_bulk_fragment_count"] = normalized["fragment_count"]
                output["alphafold_bulk_has_pdb"] = normalized["has_pdb"]
                output["alphafold_bulk_has_cif"] = normalized["has_cif"]
                output["alphafold_bulk_fragments"] = normalized["fragments"]
            output_rows.append(output)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    status_counts: dict[str, int] = {}
    for row in output_rows:
        status = row["alphafold_bulk_status"]
        status_counts[status] = status_counts.get(status, 0) + 1
    summary = {
        "input_rows": len(output_rows),
        "manifest_accessions": len(manifest_index),
        "status_counts": status_counts,
        "output": str(output_path),
    }
    Path(args.summary).write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
