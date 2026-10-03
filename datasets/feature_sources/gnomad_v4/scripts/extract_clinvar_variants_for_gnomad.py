#!/usr/bin/env python3
"""Create a gnomAD annotation input table from ClinVar variant_summary."""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path


DEFAULT_LABELS = {
    "Pathogenic",
    "Likely pathogenic",
    "Pathogenic/Likely pathogenic",
    "Uncertain significance",
    "Benign",
    "Likely benign",
    "Benign/Likely benign",
}


def normalize_header(fieldnames: list[str]) -> list[str]:
    return [name.lstrip("#") for name in fieldnames]


def make_variant_id(row: dict[str, str]) -> str | None:
    chrom = row.get("Chromosome", "").strip().removeprefix("chr")
    pos = row.get("PositionVCF", "").strip()
    ref = row.get("ReferenceAlleleVCF", "").strip()
    alt = row.get("AlternateAlleleVCF", "").strip()
    if not chrom or not pos or not ref or not alt or "-" in {chrom, pos, ref, alt}:
        return None
    return f"{chrom}-{pos}-{ref}-{alt}"


def wanted(row: dict[str, str], args: argparse.Namespace) -> bool:
    if row.get("Assembly") != args.assembly:
        return False
    if args.germline_only and row.get("OriginSimple") not in {"germline", "germline/somatic"}:
        return False
    if args.exclude_no_assertion and row.get("ReviewStatus") == "no assertion criteria provided":
        return False
    if args.labels and row.get("ClinicalSignificance") not in args.labels:
        return False
    if args.genes and row.get("GeneSymbol") not in args.genes:
        return False
    return make_variant_id(row) is not None


def read_gene_panel(path: Path | None) -> set[str] | None:
    if path is None:
        return None
    genes = set()
    for line in path.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            genes.add(line)
    return genes


def extract(args: argparse.Namespace) -> int:
    args.genes = read_gene_panel(args.gene_panel)
    args.labels = set(args.labels.split(",")) if args.labels else DEFAULT_LABELS

    args.output.parent.mkdir(parents=True, exist_ok=True)
    seen = set()
    written = 0
    fields = [
        "variant_id",
        "VariationID",
        "AlleleID",
        "GeneSymbol",
        "ClinicalSignificance",
        "ClinSigSimple",
        "ReviewStatus",
        "NumberSubmitters",
        "Assembly",
        "Chromosome",
        "PositionVCF",
        "ReferenceAlleleVCF",
        "AlternateAlleleVCF",
    ]
    with gzip.open(args.clinvar, "rt", newline="") as source, args.output.open("w", newline="") as out:
        reader = csv.DictReader(source, delimiter="\t")
        reader.fieldnames = normalize_header(reader.fieldnames or [])
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in reader:
            if not wanted(row, args):
                continue
            variant_id = make_variant_id(row)
            if args.unique and variant_id in seen:
                continue
            seen.add(variant_id)
            writer.writerow({"variant_id": variant_id, **{field: row.get(field, "") for field in fields[1:]}})
            written += 1
            if args.limit and written >= args.limit:
                break
    return written


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinvar", default=Path("datasets/clinvar/variant_summary.txt.gz"), type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--assembly", default="GRCh38")
    parser.add_argument("--gene-panel", type=Path)
    parser.add_argument("--labels", help="Comma-separated ClinVar ClinicalSignificance values.")
    parser.add_argument("--germline-only", action="store_true", default=True)
    parser.add_argument("--include-somatic", dest="germline_only", action="store_false")
    parser.add_argument("--exclude-no-assertion", action="store_true")
    parser.add_argument("--unique", action="store_true", default=True)
    parser.add_argument("--allow-duplicates", dest="unique", action="store_false")
    parser.add_argument("--limit", type=int)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    written = extract(args)
    print(f"wrote {written} variants to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

