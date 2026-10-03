#!/usr/bin/env python3
"""Join AlphaMissense hg38 precomputed scores onto a variant list."""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path


DEFAULT_ALPHAMISSENSE = Path("datasets/feature_sources/insilico/raw/alphamissense/AlphaMissense_hg38.tsv.gz")


OUTPUT_FIELDS = [
    "variant_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "alphamissense_status",
    "alphamissense_missing_reason",
    "am_pathogenicity",
    "am_class",
    "am_uniprot_id",
    "am_transcript_id",
    "am_protein_variant",
]


def normalize_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def variant_key(chrom: str, pos: str, ref: str, alt: str) -> tuple[str, str, str, str]:
    return (normalize_chrom(chrom), pos.strip(), ref.strip(), alt.strip())


def parse_variant_id(variant_id: str) -> tuple[str, str, str, str]:
    parts = variant_id.strip().split("-")
    if len(parts) != 4:
        raise ValueError(f"Invalid variant_id {variant_id!r}; expected chrom-pos-ref-alt")
    return variant_key(*parts)


def read_variants(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        except csv.Error:
            dialect = csv.excel_tab
        reader = csv.DictReader(handle, dialect=dialect)
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {path}")
        normalized = {name.strip().lower(): name for name in reader.fieldnames}
        rows = []
        for row in reader:
            if "variant_id" in normalized:
                key = parse_variant_id(row[normalized["variant_id"]])
            else:
                required = ["chrom", "pos", "ref", "alt"]
                missing = [column for column in required if column not in normalized]
                if missing:
                    raise ValueError(
                        f"{path} needs variant_id or chrom,pos,ref,alt columns; missing {missing}"
                    )
                key = variant_key(
                    row[normalized["chrom"]],
                    row[normalized["pos"]],
                    row[normalized["ref"]],
                    row[normalized["alt"]],
                )
            chrom, pos, ref, alt = key
            rows.append(
                {
                    "variant_id": f"{chrom}-{pos}-{ref}-{alt}",
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "_key": key,
                }
            )
    return rows


def load_alpha_matches(path: Path, wanted: set[tuple[str, str, str, str]]) -> dict[tuple[str, str, str, str], dict[str, str]]:
    matches = {}
    with gzip.open(path, "rt", newline="") as handle:
        header = None
        for line in handle:
            if line.startswith("##") or line.startswith("# Copyright") or line.strip() == "#":
                continue
            if line.startswith("#CHROM"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                break
        if header is None:
            raise ValueError(f"Could not find AlphaMissense header in {path}")
        reader = csv.DictReader(handle, fieldnames=header, delimiter="\t")
        for row in reader:
            key = variant_key(row["CHROM"], row["POS"], row["REF"], row["ALT"])
            if key in wanted:
                matches[key] = row
                if len(matches) == len(wanted):
                    break
    return matches


def join(args: argparse.Namespace) -> int:
    variants = read_variants(args.input)
    wanted = {row["_key"] for row in variants}
    matches = load_alpha_matches(args.alphamissense, wanted)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        for row in variants:
            match = matches.get(row["_key"])
            out = {field: row.get(field, "") for field in ["variant_id", "chrom", "pos", "ref", "alt"]}
            if match:
                out.update(
                    {
                        "alphamissense_status": "ok",
                        "alphamissense_missing_reason": "",
                        "am_pathogenicity": match.get("am_pathogenicity", ""),
                        "am_class": match.get("am_class", ""),
                        "am_uniprot_id": match.get("uniprot_id", ""),
                        "am_transcript_id": match.get("transcript_id", ""),
                        "am_protein_variant": match.get("protein_variant", ""),
                    }
                )
            else:
                out.update(
                    {
                        "alphamissense_status": "not_applicable_or_not_found",
                        "alphamissense_missing_reason": "No exact hg38 SNV/missense AlphaMissense record for chrom-pos-ref-alt",
                    }
                )
            writer.writerow(out)
    found = sum(1 for row in variants if row["_key"] in matches)
    print(f"AlphaMissense join: total={len(variants)} found={found} missing={len(variants)-found}", file=sys.stderr)
    return 0


def validate(args: argparse.Namespace) -> int:
    with args.input.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{args.input} has no data rows")
    required = ["variant_id", "alphamissense_status", "alphamissense_missing_reason"]
    missing_columns = [column for column in required if column not in (reader.fieldnames or [])]
    if missing_columns:
        raise ValueError(f"{args.input} missing required columns: {missing_columns}")
    statuses: dict[str, int] = {}
    bad_ok = []
    for i, row in enumerate(rows, start=2):
        statuses[row["alphamissense_status"]] = statuses.get(row["alphamissense_status"], 0) + 1
        if row["alphamissense_status"] == "ok" and not row.get("am_pathogenicity"):
            bad_ok.append(i)
    if bad_ok:
        raise ValueError(f"Rows marked ok without am_pathogenicity: {bad_ok[:10]}")
    print(f"validated AlphaMissense rows={len(rows)} statuses={statuses}", file=sys.stderr)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(required=True)

    join_parser = subparsers.add_parser("join")
    join_parser.add_argument("--input", required=True, type=Path)
    join_parser.add_argument("--output", required=True, type=Path)
    join_parser.add_argument("--alphamissense", default=DEFAULT_ALPHAMISSENSE, type=Path)
    join_parser.set_defaults(func=join)

    validate_parser = subparsers.add_parser("validate")
    validate_parser.add_argument("--input", required=True, type=Path)
    validate_parser.set_defaults(func=validate)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())

