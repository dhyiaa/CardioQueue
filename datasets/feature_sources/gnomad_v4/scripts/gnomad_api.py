#!/usr/bin/env python3
"""Targeted gnomAD v4 GraphQL queries for project feature tables."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path


API_URL = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r4"
REFERENCE_GENOME = "GRCh38"
API_VERSION_NOTE = "gnomAD browser GraphQL API, dataset gnomad_r4, GRCh38"


VARIANT_QUERY = """
query($variantId:String!,$dataset:DatasetId!){
  variant(variantId:$variantId,dataset:$dataset){
    variant_id chrom pos ref alt flags
    exome{
      ac an af homozygote_count hemizygote_count filters
      faf95{popmax popmax_population}
      faf99{popmax popmax_population}
      populations{id ac an homozygote_count hemizygote_count}
    }
    genome{
      ac an af homozygote_count hemizygote_count filters
      faf95{popmax popmax_population}
      faf99{popmax popmax_population}
      populations{id ac an homozygote_count hemizygote_count}
    }
    joint{
      ac an homozygote_count hemizygote_count filters
      faf95{popmax popmax_population}
      faf99{popmax popmax_population}
      populations{id ac an homozygote_count hemizygote_count}
    }
    transcript_consequences{
      gene_symbol gene_id transcript_id is_mane_select is_canonical
      major_consequence consequence_terms hgvsc hgvsp
      lof lof_filter lof_flags polyphen_prediction sift_prediction
    }
  }
}
"""


GENE_QUERY = """
query($gene_symbol:String!,$reference_genome:ReferenceGenomeId!){
  gene(gene_symbol:$gene_symbol,reference_genome:$reference_genome){
    symbol gene_id chrom start stop canonical_transcript_id flags
    gnomad_constraint{
      pLI pli mis_z lof_z oe_lof oe_lof_lower oe_lof_upper
      oe_lof_percentile oe_mis oe_mis_lower oe_mis_upper flags
    }
  }
}
"""


VARIANT_FIELDS = [
    "requested_variant_id",
    "api_dataset",
    "api_reference_genome",
    "query_status",
    "found",
    "missing_reason",
    "variant_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "variant_flags",
    "has_exome",
    "has_genome",
    "has_joint",
    "has_transcript_consequence",
    "exome_ac",
    "exome_an",
    "exome_af",
    "exome_homozygote_count",
    "exome_hemizygote_count",
    "exome_filters",
    "exome_popmax_af",
    "exome_popmax_population",
    "exome_faf95_popmax",
    "exome_faf95_popmax_population",
    "exome_faf99_popmax",
    "exome_faf99_popmax_population",
    "genome_ac",
    "genome_an",
    "genome_af",
    "genome_homozygote_count",
    "genome_hemizygote_count",
    "genome_filters",
    "genome_popmax_af",
    "genome_popmax_population",
    "joint_ac",
    "joint_an",
    "joint_homozygote_count",
    "joint_hemizygote_count",
    "joint_filters",
    "joint_faf95_popmax",
    "joint_faf95_popmax_population",
    "joint_faf99_popmax",
    "joint_faf99_popmax_population",
    "gene_symbol",
    "gene_id",
    "transcript_id",
    "is_mane_select",
    "is_canonical",
    "major_consequence",
    "consequence_terms",
    "hgvsc",
    "hgvsp",
    "lof",
    "lof_filter",
    "lof_flags",
    "polyphen_prediction",
    "sift_prediction",
]


GENE_FIELDS = [
    "requested_gene_symbol",
    "api_dataset",
    "api_reference_genome",
    "query_status",
    "found",
    "missing_reason",
    "gene_symbol",
    "gene_id",
    "chrom",
    "start",
    "stop",
    "canonical_transcript_id",
    "gene_flags",
    "has_constraint",
    "pLI",
    "pli",
    "mis_z",
    "lof_z",
    "oe_lof",
    "oe_lof_lower",
    "oe_lof_upper",
    "loeuf",
    "oe_lof_percentile",
    "oe_mis",
    "oe_mis_lower",
    "oe_mis_upper",
    "constraint_flags",
]


def post_graphql(
    query: str,
    variables: dict,
    retries: int = 6,
    missing_ok_key: str | None = None,
) -> dict:
    payload = json.dumps({"query": query, "variables": variables}).encode("utf-8")
    request = urllib.request.Request(
        API_URL,
        data=payload,
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                result = json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            if attempt == retries:
                raise RuntimeError(f"GraphQL request failed after {retries} attempts") from exc
            retry_after = exc.headers.get("Retry-After")
            if retry_after and retry_after.isdigit():
                delay = int(retry_after)
            elif exc.code == 429:
                delay = min(60, 5 * attempt)
            else:
                delay = 2 * attempt
            time.sleep(delay)
            continue
        except (urllib.error.URLError, TimeoutError) as exc:
            if attempt == retries:
                raise RuntimeError(f"GraphQL request failed after {retries} attempts") from exc
            time.sleep(2 * attempt)
            continue
        if result.get("errors") and missing_ok_key:
            messages = [error.get("message", "") for error in result["errors"]]
            if all("not found" in message.lower() for message in messages):
                return {missing_ok_key: None}
        if result.get("errors"):
            raise RuntimeError(json.dumps(result["errors"], indent=2))
        return result["data"]
    raise RuntimeError("unreachable")


def stable_id(value: str) -> str:
    digest = hashlib.sha1(value.encode("utf-8")).hexdigest()[:12]
    safe = "".join(char if char.isalnum() or char in "._-" else "_" for char in value)
    return f"{safe}.{digest}"


def write_raw_json(raw_json_dir: Path | None, prefix: str, key: str, data: dict) -> None:
    if raw_json_dir is None:
        return
    raw_json_dir.mkdir(parents=True, exist_ok=True)
    path = raw_json_dir / f"{prefix}.{stable_id(key)}.json"
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def read_variant_ids(path: Path) -> list[str]:
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
                variants.append(row[normalized["variant_id"]].strip())
            else:
                required = ["chrom", "pos", "ref", "alt"]
                missing = [col for col in required if col not in normalized]
                if missing:
                    raise ValueError(
                        f"{path} needs variant_id or chrom,pos,ref,alt columns; missing {missing}"
                    )
                chrom = row[normalized["chrom"]].strip().removeprefix("chr")
                pos = row[normalized["pos"]].strip()
                ref = row[normalized["ref"]].strip()
                alt = row[normalized["alt"]].strip()
                variants.append(f"{chrom}-{pos}-{ref}-{alt}")
    return [variant for variant in variants if variant]


def read_gene_symbols(path: Path) -> list[str]:
    genes = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            genes.append(line)
    return genes


def popmax_from_populations(block: dict | None) -> tuple[float | None, str | None]:
    if not block:
        return None, None
    best_af = None
    best_pop = None
    for pop in block.get("populations") or []:
        an = pop.get("an") or 0
        ac = pop.get("ac") or 0
        if an <= 0:
            continue
        af = ac / an
        if best_af is None or af > best_af:
            best_af = af
            best_pop = pop.get("id")
    return best_af, best_pop


def pick_transcript(consequences: list[dict] | None) -> dict:
    if not consequences:
        return {}
    for consequence in consequences:
        if consequence.get("is_mane_select"):
            return consequence
    for consequence in consequences:
        if consequence.get("is_canonical"):
            return consequence
    return consequences[0]


def flatten_variant(record: dict | None, requested_variant_id: str) -> dict:
    if not record:
        return {
            "requested_variant_id": requested_variant_id,
            "api_dataset": DATASET,
            "api_reference_genome": REFERENCE_GENOME,
            "query_status": "not_found",
            "found": False,
            "missing_reason": "gnomAD API returned Variant not found",
        }
    exome = record.get("exome") or {}
    genome = record.get("genome") or {}
    joint = record.get("joint") or {}
    consequences = record.get("transcript_consequences") or []
    exome_popmax, exome_popmax_pop = popmax_from_populations(exome)
    genome_popmax, genome_popmax_pop = popmax_from_populations(genome)
    tx = pick_transcript(consequences)
    return {
        "requested_variant_id": requested_variant_id,
        "api_dataset": DATASET,
        "api_reference_genome": REFERENCE_GENOME,
        "query_status": "ok",
        "found": True,
        "missing_reason": "",
        "variant_id": record.get("variant_id"),
        "chrom": record.get("chrom"),
        "pos": record.get("pos"),
        "ref": record.get("ref"),
        "alt": record.get("alt"),
        "variant_flags": ";".join(record.get("flags") or []),
        "has_exome": bool(record.get("exome")),
        "has_genome": bool(record.get("genome")),
        "has_joint": bool(record.get("joint")),
        "has_transcript_consequence": bool(consequences),
        "exome_ac": exome.get("ac"),
        "exome_an": exome.get("an"),
        "exome_af": exome.get("af"),
        "exome_homozygote_count": exome.get("homozygote_count"),
        "exome_hemizygote_count": exome.get("hemizygote_count"),
        "exome_filters": ";".join(exome.get("filters") or []),
        "exome_popmax_af": exome_popmax,
        "exome_popmax_population": exome_popmax_pop,
        "exome_faf95_popmax": (exome.get("faf95") or {}).get("popmax"),
        "exome_faf95_popmax_population": (exome.get("faf95") or {}).get("popmax_population"),
        "exome_faf99_popmax": (exome.get("faf99") or {}).get("popmax"),
        "exome_faf99_popmax_population": (exome.get("faf99") or {}).get("popmax_population"),
        "genome_ac": genome.get("ac"),
        "genome_an": genome.get("an"),
        "genome_af": genome.get("af"),
        "genome_homozygote_count": genome.get("homozygote_count"),
        "genome_hemizygote_count": genome.get("hemizygote_count"),
        "genome_filters": ";".join(genome.get("filters") or []),
        "genome_popmax_af": genome_popmax,
        "genome_popmax_population": genome_popmax_pop,
        "joint_ac": joint.get("ac"),
        "joint_an": joint.get("an"),
        "joint_homozygote_count": joint.get("homozygote_count"),
        "joint_hemizygote_count": joint.get("hemizygote_count"),
        "joint_filters": ";".join(joint.get("filters") or []),
        "joint_faf95_popmax": (joint.get("faf95") or {}).get("popmax"),
        "joint_faf95_popmax_population": (joint.get("faf95") or {}).get("popmax_population"),
        "joint_faf99_popmax": (joint.get("faf99") or {}).get("popmax"),
        "joint_faf99_popmax_population": (joint.get("faf99") or {}).get("popmax_population"),
        "gene_symbol": tx.get("gene_symbol"),
        "gene_id": tx.get("gene_id"),
        "transcript_id": tx.get("transcript_id"),
        "is_mane_select": tx.get("is_mane_select"),
        "is_canonical": tx.get("is_canonical"),
        "major_consequence": tx.get("major_consequence"),
        "consequence_terms": ";".join(tx.get("consequence_terms") or []),
        "hgvsc": tx.get("hgvsc"),
        "hgvsp": tx.get("hgvsp"),
        "lof": tx.get("lof"),
        "lof_filter": tx.get("lof_filter"),
        "lof_flags": tx.get("lof_flags"),
        "polyphen_prediction": tx.get("polyphen_prediction"),
        "sift_prediction": tx.get("sift_prediction"),
    }


def flatten_gene(record: dict | None, gene_symbol: str) -> dict:
    if not record:
        return {
            "requested_gene_symbol": gene_symbol,
            "api_dataset": DATASET,
            "api_reference_genome": REFERENCE_GENOME,
            "query_status": "not_found",
            "found": False,
            "missing_reason": "gnomAD API returned Gene not found",
        }
    constraint = record.get("gnomad_constraint") or {}
    return {
        "requested_gene_symbol": gene_symbol,
        "api_dataset": DATASET,
        "api_reference_genome": REFERENCE_GENOME,
        "query_status": "ok",
        "found": True,
        "missing_reason": "",
        "gene_symbol": record.get("symbol"),
        "gene_id": record.get("gene_id"),
        "chrom": record.get("chrom"),
        "start": record.get("start"),
        "stop": record.get("stop"),
        "canonical_transcript_id": record.get("canonical_transcript_id"),
        "gene_flags": ";".join(record.get("flags") or []),
        "has_constraint": bool(record.get("gnomad_constraint")),
        "pLI": constraint.get("pLI"),
        "pli": constraint.get("pli"),
        "mis_z": constraint.get("mis_z"),
        "lof_z": constraint.get("lof_z"),
        "oe_lof": constraint.get("oe_lof"),
        "oe_lof_lower": constraint.get("oe_lof_lower"),
        "oe_lof_upper": constraint.get("oe_lof_upper"),
        "loeuf": constraint.get("oe_lof_upper"),
        "oe_lof_percentile": constraint.get("oe_lof_percentile"),
        "oe_mis": constraint.get("oe_mis"),
        "oe_mis_lower": constraint.get("oe_mis_lower"),
        "oe_mis_upper": constraint.get("oe_mis_upper"),
        "constraint_flags": ";".join(constraint.get("flags") or []),
    }


def write_tsv(rows: list[dict], output: Path, preferred_fields: list[str]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    extra_fields = sorted({key for row in rows for key in row if key not in preferred_fields})
    fieldnames = preferred_fields + extra_fields
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def print_summary(rows: list[dict], label: str) -> None:
    total = len(rows)
    found = sum(row.get("found") is True for row in rows)
    not_found = sum(row.get("query_status") == "not_found" for row in rows)
    print(f"{label}: total={total} found={found} not_found={not_found}", file=sys.stderr)


def annotate_variants(args: argparse.Namespace) -> None:
    rows = []
    for i, variant_id in enumerate(read_variant_ids(args.input), start=1):
        if args.sleep:
            time.sleep(args.sleep)
        data = post_graphql(
            VARIANT_QUERY,
            {"variantId": variant_id, "dataset": DATASET},
            missing_ok_key="variant",
        )
        write_raw_json(args.raw_json_dir, "variant", variant_id, data)
        rows.append(flatten_variant(data.get("variant"), variant_id))
        if args.progress and i % args.batch_size == 0:
            print(f"annotated {i} variants", file=sys.stderr)
    write_tsv(rows, args.output, VARIANT_FIELDS)
    print_summary(rows, "variants")


def gene_constraints(args: argparse.Namespace) -> None:
    rows = []
    for i, gene_symbol in enumerate(read_gene_symbols(args.input), start=1):
        if args.sleep:
            time.sleep(args.sleep)
        data = post_graphql(
            GENE_QUERY,
            {"gene_symbol": gene_symbol, "reference_genome": REFERENCE_GENOME},
            missing_ok_key="gene",
        )
        write_raw_json(args.raw_json_dir, "gene", gene_symbol, data)
        rows.append(flatten_gene(data.get("gene"), gene_symbol))
        if args.progress and i % args.batch_size == 0:
            print(f"queried {i} genes", file=sys.stderr)
    write_tsv(rows, args.output, GENE_FIELDS)
    print_summary(rows, "genes")


def validate_table(args: argparse.Namespace) -> None:
    with args.input.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{args.input} has no data rows")
    required = ["query_status", "found"]
    missing_columns = [column for column in required if column not in (reader.fieldnames or [])]
    if missing_columns:
        raise ValueError(f"{args.input} missing required columns: {missing_columns}")
    statuses = {}
    for row in rows:
        statuses[row["query_status"]] = statuses.get(row["query_status"], 0) + 1
    found_rows = [row for row in rows if row.get("found") == "True"]
    if args.kind == "variants":
        required_when_found = ["variant_id", "chrom", "pos", "ref", "alt"]
    else:
        required_when_found = ["gene_symbol", "gene_id", "chrom", "start", "stop"]
    bad_found = []
    for i, row in enumerate(found_rows, start=2):
        missing = [column for column in required_when_found if not row.get(column)]
        if missing:
            bad_found.append((i, missing))
    if bad_found:
        preview = "; ".join(f"line {line}: {missing}" for line, missing in bad_found[:10])
        raise ValueError(f"Found rows with missing required fields: {preview}")
    if args.max_not_found_fraction is not None:
        not_found = statuses.get("not_found", 0)
        fraction = not_found / len(rows)
        if fraction > args.max_not_found_fraction:
            raise ValueError(
                f"not_found fraction {fraction:.3f} exceeds limit {args.max_not_found_fraction:.3f}"
            )
    print(
        f"validated {args.kind}: rows={len(rows)} statuses={json.dumps(statuses, sort_keys=True)}",
        file=sys.stderr,
    )


def doctor(args: argparse.Namespace) -> None:
    variant_data = post_graphql(
        VARIANT_QUERY,
        {"variantId": args.variant_id, "dataset": DATASET},
        missing_ok_key="variant",
    )
    gene_data = post_graphql(
        GENE_QUERY,
        {"gene_symbol": args.gene_symbol, "reference_genome": REFERENCE_GENOME},
        missing_ok_key="gene",
    )
    variant_row = flatten_variant(variant_data.get("variant"), args.variant_id)
    gene_row = flatten_gene(gene_data.get("gene"), args.gene_symbol)
    if variant_row["query_status"] != "ok":
        raise RuntimeError(f"Variant doctor query failed status check: {variant_row}")
    if gene_row["query_status"] != "ok":
        raise RuntimeError(f"Gene doctor query failed status check: {gene_row}")
    required_variant_fields = ["variant_id", "chrom", "pos", "ref", "alt"]
    missing_variant = [field for field in required_variant_fields if not variant_row.get(field)]
    required_gene_fields = ["gene_symbol", "gene_id", "chrom", "start", "stop"]
    missing_gene = [field for field in required_gene_fields if not gene_row.get(field)]
    if missing_variant or missing_gene:
        raise RuntimeError(
            f"Doctor query returned incomplete data: variant={missing_variant}, gene={missing_gene}"
        )
    print(API_VERSION_NOTE, file=sys.stderr)
    print(f"doctor variant ok: {variant_row['variant_id']}", file=sys.stderr)
    print(f"doctor gene ok: {gene_row['gene_symbol']} {gene_row['gene_id']}", file=sys.stderr)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(required=True)

    variant_parser = subparsers.add_parser("annotate-variants")
    variant_parser.add_argument("--input", required=True, type=Path)
    variant_parser.add_argument("--output", required=True, type=Path)
    variant_parser.add_argument("--batch-size", default=20, type=int)
    variant_parser.add_argument("--sleep", default=1.0, type=float)
    variant_parser.add_argument("--progress", action="store_true")
    variant_parser.add_argument("--raw-json-dir", type=Path)
    variant_parser.set_defaults(func=annotate_variants)

    gene_parser = subparsers.add_parser("gene-constraints")
    gene_parser.add_argument("--input", required=True, type=Path)
    gene_parser.add_argument("--output", required=True, type=Path)
    gene_parser.add_argument("--batch-size", default=20, type=int)
    gene_parser.add_argument("--sleep", default=1.0, type=float)
    gene_parser.add_argument("--progress", action="store_true")
    gene_parser.add_argument("--raw-json-dir", type=Path)
    gene_parser.set_defaults(func=gene_constraints)

    validate_parser = subparsers.add_parser("validate-table")
    validate_parser.add_argument("--input", required=True, type=Path)
    validate_parser.add_argument("--kind", choices=["variants", "genes"], required=True)
    validate_parser.add_argument("--max-not-found-fraction", type=float)
    validate_parser.set_defaults(func=validate_table)

    doctor_parser = subparsers.add_parser("doctor")
    doctor_parser.add_argument("--variant-id", default="14-23425386-C-T")
    doctor_parser.add_argument("--gene-symbol", default="MYH7")
    doctor_parser.set_defaults(func=doctor)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
