#!/usr/bin/env python3
"""Calculate FoldX DDG for the frozen primary-model cohort's mapped missense variants."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from run_foldxpro_ddg_pilot import FOLDXPRO, ROOT, run_one


DEFAULT_TABLE = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
DEFAULT_AUDIT = ROOT / "datasets/feature_sources/protein_structure/interim/registry_structure_mapping_audit.tsv"
DEFAULT_STRUCTURE = ROOT / "datasets/feature_sources/protein_structure/interim/registry_structure_features.tsv"
DEFAULT_WORK = ROOT / "datasets/feature_sources/protein_structure/interim/registry_foldxpro_ddg"
DEFAULT_OUTPUT = ROOT / "datasets/feature_sources/protein_structure/interim/registry_foldxpro_ddg.tsv"
DEFAULT_SUMMARY = ROOT / "datasets/feature_sources/protein_structure/interim/registry_foldxpro_ddg.summary.json"

OUTPUT_COLUMNS = [
    "variant_id",
    "split_source_heldout",
    "primary_binary_label",
    "primary_gene",
    "uniprot_accession",
    "protein_variant",
    "protein_global_position",
    "alphafold_fragment_member",
    "alphafold_local_position",
    "alphafold_residue_plddt",
    "foldx_calculation_key",
    "foldx_mutation",
    "ddg_method",
    "ddg_status",
    "ddg_missing_reason",
    "ddg_kcal_mol",
    "ddg_abs",
    "foldxpro_returncode",
    "foldxpro_mutant_pdb",
    "foldxpro_dif_fxout",
    "foldxpro_stdout",
    "foldxpro_stderr",
    "work_dir",
]


def truthy(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def load_candidates(args: argparse.Namespace) -> pd.DataFrame:
    table = pd.read_csv(
        args.table,
        sep="\t",
        usecols=[
            "variant_id",
            "split_source_heldout",
            "training_slice_primary_binary",
            "primary_binary_label",
            "primary_gene",
        ],
        low_memory=False,
    )
    audit = pd.read_csv(
        args.audit,
        sep="\t",
        usecols=[
            "variant_id",
            "audit_uniprot_accession",
            "audit_protein_variant",
            "audit_mapping_status",
            "audit_fragment_member",
            "audit_local_position",
        ],
        low_memory=False,
    )
    structure = pd.read_csv(
        args.structure,
        sep="\t",
        usecols=["variant_id", "structure_plddt"],
        low_memory=False,
    )
    frame = table.merge(audit, on="variant_id", validate="one_to_one").merge(
        structure, on="variant_id", validate="one_to_one"
    )
    frame = frame[
        truthy(frame["training_slice_primary_binary"])
        & frame["audit_mapping_status"].eq("ok")
        & pd.to_numeric(frame["structure_plddt"], errors="coerce").ge(args.min_plddt)
    ].copy()
    parsed = frame["audit_protein_variant"].astype(str).str.extract(r"^([A-Z])(\d+)([A-Z])$")
    frame["protein_ref_aa"] = parsed[0]
    frame["protein_global_position"] = pd.to_numeric(parsed[1], errors="coerce").astype("Int64")
    frame["protein_alt_aa"] = parsed[2]
    frame = frame.dropna(subset=["protein_ref_aa", "protein_global_position", "protein_alt_aa"])
    frame["foldx_calculation_key"] = (
        frame["audit_uniprot_accession"].astype(str)
        + "|"
        + frame["audit_fragment_member"].astype(str)
        + "|"
        + frame["audit_local_position"].astype("Int64").astype(str)
        + "|"
        + frame["protein_ref_aa"]
        + ">"
        + frame["protein_alt_aa"]
    )
    return frame.sort_values(["audit_uniprot_accession", "protein_global_position", "variant_id"])


def run_calculation(row: pd.Series, args: argparse.Namespace) -> dict[str, object]:
    local_position = int(row.audit_local_position)
    pilot_row = {
        "source_dataset": "registry_primary_binary",
        "variant_uid": str(row.foldx_calculation_key),
        "variant_key": str(row.foldx_calculation_key),
        "target_3class": "",
        "gene": str(row.primary_gene),
        "uniprot_accession": str(row.audit_uniprot_accession),
        "protein_position": str(local_position),
        "protein_ref_aa": str(row.protein_ref_aa),
        "protein_alt_aa": str(row.protein_alt_aa),
        "alphafold_residue_plddt": str(row.structure_plddt),
    }
    run_args = SimpleNamespace(
        foldxpro=args.foldxpro,
        work_dir=str(Path(args.work_dir).resolve()),
        timeout_seconds=args.timeout_seconds,
    )
    result = run_one(pilot_row, str(row.audit_fragment_member), run_args)
    result["foldx_calculation_key"] = row.foldx_calculation_key
    return result


def base_output(row: pd.Series) -> dict[str, object]:
    return {
        "variant_id": row.variant_id,
        "split_source_heldout": row.split_source_heldout,
        "primary_binary_label": row.primary_binary_label,
        "primary_gene": row.primary_gene,
        "uniprot_accession": row.audit_uniprot_accession,
        "protein_variant": row.audit_protein_variant,
        "protein_global_position": int(row.protein_global_position),
        "alphafold_fragment_member": row.audit_fragment_member,
        "alphafold_local_position": int(row.audit_local_position),
        "alphafold_residue_plddt": row.structure_plddt,
        "foldx_calculation_key": row.foldx_calculation_key,
    }


def write_outputs(frame: pd.DataFrame, results: dict[str, dict[str, object]], args: argparse.Namespace) -> None:
    rows = []
    for row in frame.itertuples(index=False):
        output = base_output(row)
        output.update(results.get(row.foldx_calculation_key, {}))
        output["protein_global_position"] = int(row.protein_global_position)
        output["alphafold_local_position"] = int(row.audit_local_position)
        rows.append(output)
    output = pd.DataFrame(rows).reindex(columns=OUTPUT_COLUMNS)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)
    unique = output.drop_duplicates("foldx_calculation_key")
    summary = {
        "status": "complete" if len(results) == frame["foldx_calculation_key"].nunique() else "partial",
        "selection_rule": "Frozen primary-binary cohort, registry mapping ok, pLDDT at or above threshold; outcomes not used for selection",
        "min_plddt": args.min_plddt,
        "variant_rows": len(output),
        "unique_calculations": int(frame["foldx_calculation_key"].nunique()),
        "completed_calculations": len(results),
        "status_counts_unique": dict(Counter(unique["ddg_status"].fillna("missing"))),
        "coverage_by_split": {
            str(split): {"rows": len(group), "ok": int(group["ddg_status"].eq("ok").sum())}
            for split, group in output.groupby("split_source_heldout", dropna=False)
        },
        "input_table": str(Path(args.table)),
        "mapping_audit": str(Path(args.audit)),
        "structure_features": str(Path(args.structure)),
        "output": str(Path(args.output)),
        "work_dir": str(Path(args.work_dir)),
        "workers": args.workers,
    }
    Path(args.summary).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", default=str(DEFAULT_TABLE))
    parser.add_argument("--audit", default=str(DEFAULT_AUDIT))
    parser.add_argument("--structure", default=str(DEFAULT_STRUCTURE))
    parser.add_argument("--foldxpro", default=str(FOLDXPRO))
    parser.add_argument("--work-dir", default=str(DEFAULT_WORK))
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--summary", default=str(DEFAULT_SUMMARY))
    parser.add_argument("--min-plddt", type=float, default=70.0)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout-seconds", type=int, default=600)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    frame = load_candidates(args)
    representatives = frame.drop_duplicates("foldx_calculation_key", keep="first")
    if args.limit is not None:
        keys = set(representatives.head(args.limit)["foldx_calculation_key"])
        representatives = representatives[representatives["foldx_calculation_key"].isin(keys)]
        frame = frame[frame["foldx_calculation_key"].isin(keys)]

    if args.preflight_only:
        preflight = frame.apply(base_output, axis=1, result_type="expand").reindex(columns=OUTPUT_COLUMNS)
        preflight["ddg_status"] = "preflight_selected"
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        preflight.to_csv(args.output, sep="\t", index=False)
        print(json.dumps({"variant_rows": len(frame), "unique_calculations": len(representatives)}, indent=2))
        return 0

    results: dict[str, dict[str, object]] = {}
    output_path = Path(args.output)
    if output_path.exists() and not args.force:
        existing = pd.read_csv(output_path, sep="\t", low_memory=False)
        for row in existing.drop_duplicates("foldx_calculation_key").to_dict("records"):
            status = row.get("ddg_status")
            if (
                row.get("foldx_calculation_key")
                and pd.notna(status)
                and str(status) not in {"preflight_selected", ""}
            ):
                results[str(row["foldx_calculation_key"])] = row

    pending = [row for row in representatives.itertuples(index=False) if row.foldx_calculation_key not in results]
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as executor:
        futures = {executor.submit(run_calculation, row, args): row.foldx_calculation_key for row in pending}
        for index, future in enumerate(as_completed(futures), start=1):
            key = futures[future]
            try:
                results[key] = future.result()
            except Exception as exc:
                results[key] = {"foldx_calculation_key": key, "ddg_status": "failed", "ddg_missing_reason": str(exc)[:500]}
            if index % 10 == 0 or index == len(futures):
                write_outputs(frame, results, args)
                print(f"completed {len(results)}/{len(representatives)} unique calculations", flush=True)

    write_outputs(frame, results, args)
    print(Path(args.summary).read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
