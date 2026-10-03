#!/usr/bin/env python3
"""Join dbNSFP GRCh38 annotations onto a variant list.

Input can be either a master registry-style TSV with `variant_registry_key` or a
smaller table with `variant_id` or `chrom,pos,ref,alt`.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from collections import defaultdict
from pathlib import Path

try:
    import pysam
except ImportError as exc:  # pragma: no cover - environment guard
    raise SystemExit("pysam is required for dbNSFP BGZF/tabix queries") from exc


DEFAULT_DBNSFP = Path(
    "datasets/feature_sources/insilico/raw/dbnsfp/downloads/dbNSFP5.3.1a_grch38.gz"
)

SELECTED_DBNSFP_COLUMNS = [
    "aaref",
    "aaalt",
    "aapos",
    "genename",
    "Ensembl_geneid",
    "Ensembl_transcriptid",
    "Ensembl_proteinid",
    "Uniprot_acc",
    "Uniprot_entry",
    "HGVSc_VEP",
    "HGVSp_VEP",
    "APPRIS",
    "GENCODE_basic",
    "VEP_canonical",
    "MANE",
    "Interpro_domain",
    "SIFT_score",
    "SIFT_pred",
    "SIFT4G_score",
    "SIFT4G_pred",
    "Polyphen2_HDIV_score",
    "Polyphen2_HDIV_pred",
    "Polyphen2_HVAR_score",
    "Polyphen2_HVAR_pred",
    "MutationTaster_score",
    "MutationTaster_pred",
    "MutationAssessor_score",
    "MutationAssessor_pred",
    "PROVEAN_score",
    "PROVEAN_pred",
    "VEST4_score",
    "MetaSVM_score",
    "MetaSVM_pred",
    "MetaLR_score",
    "MetaLR_pred",
    "M-CAP_score",
    "M-CAP_pred",
    "REVEL_score",
    "MutPred2_score",
    "MVP_score",
    "gMVP_score",
    "MPC_score",
    "PrimateAI_score",
    "PrimateAI_pred",
    "DEOGEN2_score",
    "DEOGEN2_pred",
    "BayesDel_addAF_score",
    "BayesDel_addAF_pred",
    "BayesDel_noAF_score",
    "BayesDel_noAF_pred",
    "ClinPred_score",
    "ClinPred_pred",
    "LIST-S2_score",
    "LIST-S2_pred",
    "ESM1b_score",
    "ESM1b_pred",
    "AlphaMissense_score",
    "AlphaMissense_pred",
    "CADD_raw",
    "CADD_phred",
    "DANN_score",
    "fathmm-XF_coding_score",
    "fathmm-XF_coding_pred",
    "Eigen-raw_coding",
    "Eigen-phred_coding",
    "Eigen-PC-raw_coding",
    "Eigen-PC-phred_coding",
    "GERP++_RS",
    "GERP++_RS_rankscore",
    "phyloP100way_vertebrate",
    "phyloP470way_mammalian",
    "phyloP17way_primate",
    "phastCons100way_vertebrate",
    "phastCons470way_mammalian",
    "phastCons17way_primate",
    "gnomAD4.1_joint_AF",
    "gnomAD4.1_joint_nhomalt",
    "gnomAD4.1_joint_POPMAX_AF",
    "gnomAD4.1_joint_POPMAX_nhomalt",
    "dbNSFP_POPMAX_AF",
    "dbNSFP_POPMAX_AC",
    "dbNSFP_POPMAX_POP",
]

OUTPUT_FIELDS = [
    "variant_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "dbnsfp_status",
    "dbnsfp_missing_reason",
    *[f"dbnsfp_{column}" for column in SELECTED_DBNSFP_COLUMNS],
]


def normalize_chrom(chrom: str) -> str:
    chrom = str(chrom or "").strip()
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def normalize_allele(allele: str) -> str:
    return str(allele or "").strip().upper()


def variant_key(chrom: str, pos: str, ref: str, alt: str) -> tuple[str, str, str, str]:
    return (normalize_chrom(chrom), str(pos or "").strip(), normalize_allele(ref), normalize_allele(alt))


def parse_variant_id(variant_id: str) -> tuple[str, str, str, str]:
    parts = str(variant_id or "").strip().split("-")
    if len(parts) != 4:
        raise ValueError(f"Invalid variant id {variant_id!r}; expected chrom-pos-ref-alt")
    return variant_key(*parts)


def read_variants(path: Path, limit: int | None = None) -> list[dict[str, str]]:
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
            if "variant_registry_key" in normalized:
                key = parse_variant_id(row[normalized["variant_registry_key"]])
            elif "variant_id" in normalized:
                key = parse_variant_id(row[normalized["variant_id"]])
            else:
                required = ["chrom", "pos", "ref", "alt"]
                missing = [column for column in required if column not in normalized]
                if missing:
                    raise ValueError(
                        f"{path} needs variant_registry_key, variant_id, or chrom,pos,ref,alt; missing {missing}"
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
            if limit is not None and len(rows) >= limit:
                break
    return rows


def dbnsfp_header(tabix: pysam.TabixFile) -> list[str]:
    header_lines = list(tabix.header)
    if not header_lines:
        raise ValueError("dbNSFP file has no tabix header")
    header = header_lines[-1].lstrip("#").rstrip("\n").split("\t")
    required = ["chr", "pos(1-based)", "ref", "alt"]
    missing = [column for column in required if column not in header]
    if missing:
        raise ValueError(f"dbNSFP header missing required columns: {missing}")
    return header


def as_record(header: list[str], line: str) -> dict[str, str]:
    values = line.rstrip("\n").split("\t")
    return dict(zip(header, values))


def find_match(
    tabix: pysam.TabixFile,
    header: list[str],
    key: tuple[str, str, str, str],
) -> tuple[str, dict[str, str] | None, str]:
    chrom, pos, ref, alt = key
    if not chrom or not pos or not ref or not alt:
        return "not_queryable", None, "Missing chrom, pos, ref, or alt"
    try:
        pos_int = int(pos)
    except ValueError:
        return "not_queryable", None, "Position is not an integer"
    try:
        fetched = tabix.fetch(chrom, pos_int - 1, pos_int)
    except ValueError as exc:
        return "not_found", None, f"dbNSFP contig/region not found: {exc}"

    same_position = 0
    for line in fetched:
        record = as_record(header, line)
        same_position += 1
        if (
            normalize_chrom(record.get("chr", "")) == chrom
            and str(record.get("pos(1-based)", "")).strip() == pos
            and normalize_allele(record.get("ref", "")) == ref
            and normalize_allele(record.get("alt", "")) == alt
        ):
            return "ok", record, ""
    if same_position:
        return "not_found", None, "Coordinate exists in dbNSFP but exact ref/alt was not found"
    return "not_found", None, "No dbNSFP record at coordinate"


def build_intervals(variants: list[dict[str, str]], merge_gap: int) -> dict[str, list[tuple[int, int]]]:
    positions_by_chrom: dict[str, set[int]] = defaultdict(set)
    for row in variants:
        chrom, pos, _ref, _alt = row["_key"]
        try:
            positions_by_chrom[chrom].add(int(pos))
        except ValueError:
            continue

    intervals_by_chrom: dict[str, list[tuple[int, int]]] = {}
    for chrom, positions in positions_by_chrom.items():
        sorted_positions = sorted(positions)
        if not sorted_positions:
            continue
        intervals: list[tuple[int, int]] = []
        start = end = sorted_positions[0]
        for pos in sorted_positions[1:]:
            if pos - end <= merge_gap:
                end = pos
            else:
                intervals.append((start, end))
                start = end = pos
        intervals.append((start, end))
        intervals_by_chrom[chrom] = intervals
    return intervals_by_chrom


def load_batched_matches(
    tabix: pysam.TabixFile,
    header: list[str],
    variants: list[dict[str, str]],
    merge_gap: int,
) -> tuple[dict[tuple[str, str, str, str], dict[str, str]], set[tuple[str, str]]]:
    wanted = {row["_key"] for row in variants}
    wanted_positions = {(chrom, pos) for chrom, pos, _ref, _alt in wanted}
    intervals_by_chrom = build_intervals(variants, merge_gap)
    interval_count = sum(len(intervals) for intervals in intervals_by_chrom.values())
    print(
        f"dbNSFP batched fetch: variants={len(variants)} positions={len(wanted_positions)} intervals={interval_count} merge_gap={merge_gap}",
        file=sys.stderr,
    )

    matches: dict[tuple[str, str, str, str], dict[str, str]] = {}
    seen_positions: set[tuple[str, str]] = set()
    for chrom in sorted(intervals_by_chrom, key=lambda value: (len(value), value)):
        for start, end in intervals_by_chrom[chrom]:
            try:
                fetched = tabix.fetch(chrom, start - 1, end)
            except ValueError:
                continue
            for line in fetched:
                record = as_record(header, line)
                rec_chrom = normalize_chrom(record.get("chr", ""))
                rec_pos = str(record.get("pos(1-based)", "")).strip()
                if (rec_chrom, rec_pos) in wanted_positions:
                    seen_positions.add((rec_chrom, rec_pos))
                key = (
                    rec_chrom,
                    rec_pos,
                    normalize_allele(record.get("ref", "")),
                    normalize_allele(record.get("alt", "")),
                )
                if key in wanted and key not in matches:
                    matches[key] = record
    return matches, seen_positions


def write_join_output(
    variants: list[dict[str, str]],
    output: Path,
    matches: dict[tuple[str, str, str, str], dict[str, str]],
    seen_positions: set[tuple[str, str]] | None = None,
) -> Counter[str]:
    status_counts: Counter[str] = Counter()
    temp_output = output.with_suffix(output.suffix + ".tmp")
    output.parent.mkdir(parents=True, exist_ok=True)
    with temp_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        for row in variants:
            match = matches.get(row["_key"])
            chrom, pos, _ref, _alt = row["_key"]
            if match:
                status = "ok"
                reason = ""
            elif seen_positions is not None and (chrom, pos) in seen_positions:
                status = "not_found"
                reason = "Coordinate exists in dbNSFP but exact ref/alt was not found"
            elif not all(row["_key"]):
                status = "not_queryable"
                reason = "Missing chrom, pos, ref, or alt"
            else:
                status = "not_found"
                reason = "No dbNSFP record at coordinate"
            status_counts[status] += 1
            out = {field: row.get(field, "") for field in ["variant_id", "chrom", "pos", "ref", "alt"]}
            out["dbnsfp_status"] = status
            out["dbnsfp_missing_reason"] = reason
            for column in SELECTED_DBNSFP_COLUMNS:
                out[f"dbnsfp_{column}"] = match.get(column, "") if match else ""
            writer.writerow(out)
    temp_output.replace(output)
    return status_counts


def join_indexed(args: argparse.Namespace, variants: list[dict[str, str]]) -> Counter[str]:
    matches = {}
    seen_positions = set()
    status_counts: Counter[str] = Counter()
    with pysam.TabixFile(str(args.dbnsfp)) as tabix:
        header = dbnsfp_header(tabix)
        missing_selected = [column for column in SELECTED_DBNSFP_COLUMNS if column not in header]
        if missing_selected:
            raise ValueError(f"dbNSFP header missing selected columns: {missing_selected}")
        for row in variants:
            status, match, reason = find_match(tabix, header, row["_key"])
            status_counts[status] += 1
            if match:
                matches[row["_key"]] = match
                chrom, pos, _ref, _alt = row["_key"]
                seen_positions.add((chrom, pos))
            elif reason == "Coordinate exists in dbNSFP but exact ref/alt was not found":
                chrom, pos, _ref, _alt = row["_key"]
                seen_positions.add((chrom, pos))
    write_join_output(variants, args.output, matches, seen_positions)
    return status_counts


def join_batched(args: argparse.Namespace, variants: list[dict[str, str]]) -> Counter[str]:
    with pysam.TabixFile(str(args.dbnsfp)) as tabix:
        header = dbnsfp_header(tabix)
        missing_selected = [column for column in SELECTED_DBNSFP_COLUMNS if column not in header]
        if missing_selected:
            raise ValueError(f"dbNSFP header missing selected columns: {missing_selected}")
        matches, seen_positions = load_batched_matches(tabix, header, variants, args.merge_gap)
    return write_join_output(variants, args.output, matches, seen_positions)


def join(args: argparse.Namespace) -> int:
    variants = read_variants(args.input, args.limit)
    if args.strategy == "indexed":
        status_counts = join_indexed(args, variants)
    else:
        status_counts = join_batched(args, variants)
    print(f"dbNSFP join: total={len(variants)} statuses={dict(status_counts)}", file=sys.stderr)
    return 0


def validate(args: argparse.Namespace) -> int:
    with args.input.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{args.input} has no data rows")
    required = ["variant_id", "dbnsfp_status", "dbnsfp_missing_reason", "dbnsfp_REVEL_score"]
    missing = [field for field in required if field not in (reader.fieldnames or [])]
    if missing:
        raise ValueError(f"{args.input} missing required columns: {missing}")
    status_counts = Counter(row["dbnsfp_status"] for row in rows)
    bad_ok = [
        i
        for i, row in enumerate(rows, start=2)
        if row["dbnsfp_status"] == "ok" and not any(row.get(field, "") for field in ["dbnsfp_genename", "dbnsfp_CADD_phred"])
    ]
    if bad_ok:
        raise ValueError(f"Rows marked ok without core dbNSFP fields: {bad_ok[:10]}")
    print(f"validated dbNSFP rows={len(rows)} statuses={dict(status_counts)}")
    return 0


def doctor(args: argparse.Namespace) -> int:
    with pysam.TabixFile(str(args.dbnsfp)) as tabix:
        header = dbnsfp_header(tabix)
        missing_selected = [column for column in SELECTED_DBNSFP_COLUMNS if column not in header]
        if missing_selected:
            raise ValueError(f"dbNSFP header missing selected columns: {missing_selected}")
        print(f"dbNSFP doctor: contigs={len(tabix.contigs)} header_columns={len(header)}")
        for key_text in args.variant:
            key = parse_variant_id(key_text)
            status, match, reason = find_match(tabix, header, key)
            fields = {
                "variant_id": key_text,
                "status": status,
                "reason": reason,
                "genename": match.get("genename", "") if match else "",
                "REVEL_score": match.get("REVEL_score", "") if match else "",
                "CADD_phred": match.get("CADD_phred", "") if match else "",
                "gnomAD4.1_joint_AF": match.get("gnomAD4.1_joint_AF", "") if match else "",
            }
            print("\t".join(f"{k}={v}" for k, v in fields.items()))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(required=True)

    join_parser = subparsers.add_parser("join")
    join_parser.add_argument("--input", required=True, type=Path)
    join_parser.add_argument("--output", required=True, type=Path)
    join_parser.add_argument("--dbnsfp", default=DEFAULT_DBNSFP, type=Path)
    join_parser.add_argument("--limit", type=int)
    join_parser.add_argument("--strategy", choices=["batched", "indexed"], default="batched")
    join_parser.add_argument("--merge-gap", type=int, default=2000)
    join_parser.set_defaults(func=join)

    validate_parser = subparsers.add_parser("validate")
    validate_parser.add_argument("--input", required=True, type=Path)
    validate_parser.set_defaults(func=validate)

    doctor_parser = subparsers.add_parser("doctor")
    doctor_parser.add_argument("--dbnsfp", default=DEFAULT_DBNSFP, type=Path)
    doctor_parser.add_argument(
        "--variant",
        action="append",
        default=["10-110812298-G-A", "10-110812310-C-T"],
        help="Variant id in chrom-pos-ref-alt format. Can be repeated.",
    )
    doctor_parser.set_defaults(func=doctor)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
