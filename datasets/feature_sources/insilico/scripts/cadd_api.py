#!/usr/bin/env python3
"""Targeted CADD SNV annotation via the official CADD API."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path


API_TEMPLATE = "https://cadd.gs.washington.edu/api/v1.0/{version}/{chrom}:{pos}_{ref}_{alt}"
DEFAULT_VERSION = "GRCh38-v1.7"

OUTPUT_FIELDS = [
    "variant_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "cadd_version",
    "cadd_status",
    "cadd_missing_reason",
    "cadd_rawscore",
    "cadd_phred",
]


def normalize_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def parse_variant_id(variant_id: str) -> tuple[str, str, str, str]:
    parts = variant_id.strip().split("-")
    if len(parts) != 4:
        raise ValueError(f"Invalid variant_id {variant_id!r}; expected chrom-pos-ref-alt")
    return normalize_chrom(parts[0]), parts[1].strip(), parts[2].strip(), parts[3].strip()


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
        variants = []
        for row in reader:
            if "variant_id" in normalized:
                chrom, pos, ref, alt = parse_variant_id(row[normalized["variant_id"]])
            else:
                required = ["chrom", "pos", "ref", "alt"]
                missing = [column for column in required if column not in normalized]
                if missing:
                    raise ValueError(
                        f"{path} needs variant_id or chrom,pos,ref,alt columns; missing {missing}"
                    )
                chrom = normalize_chrom(row[normalized["chrom"]])
                pos = row[normalized["pos"]].strip()
                ref = row[normalized["ref"]].strip()
                alt = row[normalized["alt"]].strip()
            variants.append(
                {
                    "variant_id": f"{chrom}-{pos}-{ref}-{alt}",
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                }
            )
    return variants


def query_cadd(variant: dict[str, str], version: str, retries: int = 4) -> list[dict]:
    url = API_TEMPLATE.format(version=version, **variant)
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return json.loads(response.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            if attempt == retries:
                raise RuntimeError(f"CADD API request failed after {retries} attempts: {url}") from exc
            time.sleep(2 * attempt)
    raise RuntimeError("unreachable")


def flatten(variant: dict[str, str], result: list[dict], version: str) -> dict[str, str]:
    out = {field: variant[field] for field in ["variant_id", "chrom", "pos", "ref", "alt"]}
    out["cadd_version"] = version
    if len(variant["ref"]) != 1 or len(variant["alt"]) != 1:
        out.update(
            {
                "cadd_status": "not_applicable",
                "cadd_missing_reason": "CADD API endpoint used here supports SNVs only",
            }
        )
        return out
    if not result:
        out.update(
            {
                "cadd_status": "not_found",
                "cadd_missing_reason": "CADD API returned an empty list for this SNV",
            }
        )
        return out
    record = result[0]
    out.update(
        {
            "cadd_status": "ok",
            "cadd_missing_reason": "",
            "cadd_rawscore": record.get("RawScore", ""),
            "cadd_phred": record.get("PHRED", ""),
        }
    )
    return out


def annotate(args: argparse.Namespace) -> int:
    rows = []
    for i, variant in enumerate(read_variants(args.input), start=1):
        if len(variant["ref"]) == 1 and len(variant["alt"]) == 1:
            if args.sleep:
                time.sleep(args.sleep)
            result = query_cadd(variant, args.version)
        else:
            result = []
        rows.append(flatten(variant, result, args.version))
        if args.progress and i % args.batch_size == 0:
            print(f"queried {i} CADD variants", file=sys.stderr)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    found = sum(row["cadd_status"] == "ok" for row in rows)
    print(f"CADD annotation: total={len(rows)} ok={found} non_ok={len(rows)-found}", file=sys.stderr)
    return 0


def validate(args: argparse.Namespace) -> int:
    with args.input.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{args.input} has no data rows")
    required = ["variant_id", "cadd_status", "cadd_missing_reason"]
    missing_columns = [column for column in required if column not in (reader.fieldnames or [])]
    if missing_columns:
        raise ValueError(f"{args.input} missing required columns: {missing_columns}")
    statuses: dict[str, int] = {}
    bad_ok = []
    for i, row in enumerate(rows, start=2):
        statuses[row["cadd_status"]] = statuses.get(row["cadd_status"], 0) + 1
        if row["cadd_status"] == "ok" and (not row.get("cadd_rawscore") or not row.get("cadd_phred")):
            bad_ok.append(i)
    if bad_ok:
        raise ValueError(f"Rows marked ok without CADD scores: {bad_ok[:10]}")
    print(f"validated CADD rows={len(rows)} statuses={statuses}", file=sys.stderr)
    return 0


def doctor(args: argparse.Namespace) -> int:
    variant = {
        "variant_id": "14-23425386-C-T",
        "chrom": "14",
        "pos": "23425386",
        "ref": "C",
        "alt": "T",
    }
    row = flatten(variant, query_cadd(variant, args.version), args.version)
    if row["cadd_status"] != "ok":
        raise RuntimeError(f"CADD doctor failed: {row}")
    print(f"CADD doctor ok: {row['variant_id']} {row['cadd_version']} PHRED={row['cadd_phred']}", file=sys.stderr)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(required=True)

    annotate_parser = subparsers.add_parser("annotate")
    annotate_parser.add_argument("--input", required=True, type=Path)
    annotate_parser.add_argument("--output", required=True, type=Path)
    annotate_parser.add_argument("--version", default=DEFAULT_VERSION)
    annotate_parser.add_argument("--sleep", default=1.0, type=float)
    annotate_parser.add_argument("--batch-size", default=20, type=int)
    annotate_parser.add_argument("--progress", action="store_true")
    annotate_parser.set_defaults(func=annotate)

    validate_parser = subparsers.add_parser("validate")
    validate_parser.add_argument("--input", required=True, type=Path)
    validate_parser.set_defaults(func=validate)

    doctor_parser = subparsers.add_parser("doctor")
    doctor_parser.add_argument("--version", default=DEFAULT_VERSION)
    doctor_parser.set_defaults(func=doctor)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())

