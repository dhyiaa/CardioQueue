#!/usr/bin/env python3
"""Combine registry residue-context features with primary-cohort FoldX DDG."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
DEFAULT_STRUCTURE = ROOT / "datasets/feature_sources/protein_structure/interim/registry_structure_features.tsv"
DEFAULT_FOLDX = ROOT / "datasets/feature_sources/protein_structure/interim/registry_foldx51_repaired_ddg.tsv"
DEFAULT_OUTPUT = ROOT / "datasets/feature_sources/protein_structure/interim/registry_structure_foldx_features.tsv"
DEFAULT_FOLDX_OUTPUT = ROOT / "datasets/feature_sources/protein_structure/interim/registry_foldx_only_features.tsv"
DEFAULT_SUMMARY = ROOT / "datasets/feature_sources/protein_structure/interim/registry_structure_foldx_features.summary.json"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", default=str(DEFAULT_STRUCTURE))
    parser.add_argument("--foldx", default=str(DEFAULT_FOLDX))
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--foldx-output", default=str(DEFAULT_FOLDX_OUTPUT))
    parser.add_argument("--summary", default=str(DEFAULT_SUMMARY))
    args = parser.parse_args()

    structure = pd.read_csv(args.structure, sep="\t", low_memory=False)
    foldx = pd.read_csv(
        args.foldx,
        sep="\t",
        usecols=["variant_id", "split_source_heldout", "ddg_status", "ddg_kcal_mol", "ddg_abs"],
        low_memory=False,
    )
    if not structure["variant_id"].is_unique or not foldx["variant_id"].is_unique:
        raise ValueError("Both inputs must contain one row per variant_id")

    foldx["structure_foldx_available"] = foldx["ddg_status"].eq("ok")
    foldx["structure_foldx_ddg"] = pd.to_numeric(foldx["ddg_kcal_mol"], errors="coerce").where(
        foldx["structure_foldx_available"]
    )
    foldx["structure_foldx_abs_ddg"] = pd.to_numeric(foldx["ddg_abs"], errors="coerce").where(
        foldx["structure_foldx_available"]
    )
    foldx_features = foldx[
        ["variant_id", "structure_foldx_available", "structure_foldx_ddg", "structure_foldx_abs_ddg"]
    ]
    output = structure.merge(foldx_features, on="variant_id", how="left", validate="one_to_one")
    output["structure_foldx_available"] = output["structure_foldx_available"].fillna(False).astype(bool)
    output.to_csv(args.output, sep="\t", index=False)
    foldx_features.to_csv(args.foldx_output, sep="\t", index=False)

    split_counts = (
        foldx.assign(ok=foldx["ddg_status"].eq("ok"))
        .groupby("split_source_heldout", dropna=False)
        .agg(rows=("variant_id", "size"), ok=("ok", "sum"))
        .astype(int)
        .to_dict("index")
    )
    summary = {
        "structure_input": str(Path(args.structure)),
        "foldx_input": str(Path(args.foldx)),
        "output": str(Path(args.output)),
        "foldx_only_output": str(Path(args.foldx_output)),
        "rows": len(output),
        "foldx_selected_rows": len(foldx),
        "foldx_numeric_rows": int(output["structure_foldx_available"].sum()),
        "foldx_coverage_by_split": split_counts,
        "added_features": [
            "structure_foldx_available",
            "structure_foldx_ddg",
            "structure_foldx_abs_ddg",
        ],
        "selection_note": "DDG was calculated for all reference-validated, pLDDT>=70 missense variants in the frozen primary-binary cohort without outcome-based selection.",
    }
    Path(args.summary).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
