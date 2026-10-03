#!/usr/bin/env python3
"""Run a resumable FoldXPro DDG batch for mapped missense substitutions."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

from run_foldxpro_ddg_pilot import (
    ALPHAFOLD_MANIFEST,
    FOLDXPRO,
    OUTPUT_FIELDS as PILOT_FIELDS,
    PROTEIN_FEATURES,
    ROOT,
    eligible,
    load_manifest,
    run_one,
    write_tsv,
)


OUT_ROOT = ROOT / "datasets/feature_sources/protein_structure/interim/foldxpro_ddg_batch"
OUT_TSV = ROOT / "datasets/feature_sources/protein_structure/interim/foldxpro_ddg_batch.tsv"
OUT_JSON = ROOT / "datasets/feature_sources/protein_structure/interim/foldxpro_ddg_batch.summary.json"
PREFLIGHT_TSV = ROOT / "datasets/feature_sources/protein_structure/interim/foldxpro_ddg_batch_preflight.tsv"

OUTPUT_FIELDS = [
    "variant_id",
    "foldx_batch_key",
    "foldx_batch_mode",
    *PILOT_FIELDS,
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normalize_pos(value: object) -> str:
    text = str(value if value is not None else "").strip()
    return text[:-2] if text.endswith(".0") else text


def variant_id(row: dict[str, str]) -> str:
    chrom = str(row.get("chrom", "")).strip()
    pos = normalize_pos(row.get("pos", ""))
    ref = str(row.get("ref", "")).strip().upper()
    alt = str(row.get("alt", "")).strip().upper()
    return f"{chrom}-{pos}-{ref}-{alt}" if chrom and pos and ref and alt else ""


def batch_key(row: dict[str, str]) -> str:
    return f"{row.get('source_dataset', '')}|{row.get('variant_uid', '')}"


def sort_key(row: dict[str, str]) -> tuple[str, str, str, float, str]:
    return (
        row.get("gene", ""),
        row.get("uniprot_accession", ""),
        normalize_pos(row.get("protein_position", "")),
        -float(row.get("alphafold_residue_plddt") or 0),
        row.get("variant_uid", ""),
    )


def select_rows(rows: list[dict[str, str]], min_plddt: float, limit: int | None) -> list[dict[str, str]]:
    selected = [row for row in rows if eligible(row, min_plddt)]
    selected.sort(key=sort_key)
    if limit is not None:
        selected = selected[:limit]
    return selected


def existing_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    return read_tsv(path)


def summarize(rows: list[dict[str, str]], selected_count: int, output: Path, args: argparse.Namespace) -> dict[str, object]:
    return {
        "input": str(Path(args.protein_features)),
        "output": str(output),
        "work_dir": str(Path(args.work_dir)),
        "mode": args.mode,
        "min_plddt": args.min_plddt,
        "selected_rows": selected_count,
        "output_rows": len(rows),
        "status_counts": dict(Counter(row.get("ddg_status", "") for row in rows)),
        "label_counts": dict(Counter(row.get("target_3class", "") for row in rows)),
        "gene_counts_top20": dict(Counter(row.get("gene", "") for row in rows).most_common(20)),
    }


def write_summary(path: Path, summary: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def preflight(selected: list[dict[str, str]], args: argparse.Namespace) -> int:
    rows = []
    for row in selected:
        out = {
            "variant_id": variant_id(row),
            "foldx_batch_key": batch_key(row),
            "foldx_batch_mode": args.mode,
        }
        for field in PILOT_FIELDS:
            out[field] = row.get(field, "")
        out["foldx_mutation"] = (
            f"{row.get('protein_ref_aa', '')}A{int(float(row.get('protein_position', '0')))}{row.get('protein_alt_aa', '')};"
        )
        out["ddg_status"] = "preflight_selected"
        rows.append(out)
    write_tsv(Path(args.preflight_output), rows, OUTPUT_FIELDS)
    summary = summarize(rows, len(selected), Path(args.preflight_output), args)
    summary["status"] = "preflight_only"
    write_summary(Path(args.summary), summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protein-features", default=str(PROTEIN_FEATURES))
    parser.add_argument("--foldxpro", default=str(FOLDXPRO))
    parser.add_argument("--work-dir", default=str(OUT_ROOT))
    parser.add_argument("--output", default=str(OUT_TSV))
    parser.add_argument("--summary", default=str(OUT_JSON))
    parser.add_argument("--preflight-output", default=str(PREFLIGHT_TSV))
    parser.add_argument("--mode", default="primary_plddt70", choices=["primary_plddt70", "sensitivity_low_plddt", "custom"])
    parser.add_argument("--min-plddt", type=float, default=70.0)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--timeout-seconds", type=int, default=600)
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--force", action="store_true", help="Ignore existing output and rerun selected rows.")
    args = parser.parse_args()

    source_rows = read_tsv(Path(args.protein_features))
    selected = select_rows(source_rows, args.min_plddt, args.limit)
    if args.preflight_only:
        return preflight(selected, args)

    manifest = load_manifest(ALPHAFOLD_MANIFEST)
    output_path = Path(args.output)
    output_rows = [] if args.force else existing_rows(output_path)
    completed = {row.get("foldx_batch_key", "") for row in output_rows if row.get("foldx_batch_key")}

    run_args = SimpleNamespace(
        foldxpro=args.foldxpro,
        work_dir=str(Path(args.work_dir).resolve()),
        timeout_seconds=args.timeout_seconds,
    )

    for row in selected:
        key = batch_key(row)
        if key in completed:
            continue
        names = manifest.get(row["uniprot_accession"], [])
        member_name = next((name for name in names if "-F1-" in name and name.endswith(".pdb.gz")), "")
        if not member_name:
            out = {field: "" for field in PILOT_FIELDS}
            for field in PILOT_FIELDS:
                out[field] = row.get(field, "")
            out["ddg_status"] = "alphafold_missing"
            out["ddg_missing_reason"] = "No AlphaFold F1 PDB member found."
        else:
            out = run_one(row, member_name, run_args)

        out = {
            "variant_id": variant_id(row),
            "foldx_batch_key": key,
            "foldx_batch_mode": args.mode,
            **out,
        }
        output_rows.append(out)
        completed.add(key)
        write_tsv(output_path, output_rows, OUTPUT_FIELDS)
        write_summary(Path(args.summary), summarize(output_rows, len(selected), output_path, args))

    summary = summarize(output_rows, len(selected), output_path, args)
    summary["status"] = "complete"
    write_summary(Path(args.summary), summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
