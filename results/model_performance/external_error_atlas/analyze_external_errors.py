#!/usr/bin/env python3
"""Describe source-heldout errors at the frozen 0.5 operating threshold."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
PREDICTIONS = (
    ROOT
    / "results/models/cardioqueue_v1_source_heldout"
    / "primary_binary_predictions_split_source_heldout.tsv"
)
MODELING_TABLE = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights.tsv"
OUT = ROOT / "results/model_performance/external_error_atlas"
PAPER_TABLES = ROOT / "paper_MS_SI/tables"

DETAIL_COLUMNS = [
    "variant_id",
    "primary_gene",
    "vep_worst_consequence",
    "vep_hgvsp",
    "alphafold_residue_plddt",
    "protein_structure_features_available",
    "foldx_ddg_status",
]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    predictions = pd.read_csv(PREDICTIONS, sep="\t", low_memory=False)
    external = predictions[predictions["split_source_heldout"].str.startswith("external_")].copy()

    details = pd.read_csv(MODELING_TABLE, sep="\t", usecols=DETAIL_COLUMNS, low_memory=False)
    details = details.drop_duplicates("variant_id")
    external = external.merge(details, on=["variant_id", "primary_gene"], how="left", validate="one_to_one")

    external["predicted_label"] = external["pathogenic_probability"].ge(0.5).astype(int)
    external["error_type"] = "correct"
    external.loc[(external["y_true"] == 1) & (external["predicted_label"] == 0), "error_type"] = "false_negative"
    external.loc[(external["y_true"] == 0) & (external["predicted_label"] == 1), "error_type"] = "false_positive"
    external["distance_from_threshold"] = (external["pathogenic_probability"] - 0.5).abs()
    external["source"] = external["split_source_heldout"].str.removeprefix("external_")

    errors = external[external["error_type"] != "correct"].copy()
    errors = errors.sort_values(["error_type", "distance_from_threshold", "source", "primary_gene"])
    errors.to_csv(OUT / "external_error_cases.tsv", sep="\t", index=False)
    errors.to_csv(PAPER_TABLES / "Table_S6b_external_error_cases.tsv", sep="\t", index=False)

    summary = (
        external.groupby(["source", "error_type"], dropna=False)
        .size()
        .rename("rows")
        .reset_index()
    )
    summary.to_csv(OUT / "external_error_counts_by_source.tsv", sep="\t", index=False)
    summary.to_csv(PAPER_TABLES / "Table_S6c_external_error_counts_by_source.tsv", sep="\t", index=False)

    consequence = (
        errors.groupby(["error_type", "vep_worst_consequence"], dropna=False)
        .size()
        .rename("rows")
        .reset_index()
        .sort_values(["error_type", "rows"], ascending=[True, False])
    )
    consequence.to_csv(OUT / "external_error_counts_by_consequence.tsv", sep="\t", index=False)

    gene = (
        errors.groupby(["error_type", "primary_gene"], dropna=False)
        .size()
        .rename("rows")
        .reset_index()
        .sort_values(["error_type", "rows", "primary_gene"], ascending=[True, False, True])
    )
    gene.to_csv(OUT / "external_error_counts_by_gene.tsv", sep="\t", index=False)

    manifest = {
        "threshold": 0.5,
        "external_rows": int(len(external)),
        "false_negatives": int((external["error_type"] == "false_negative").sum()),
        "false_positives": int((external["error_type"] == "false_positive").sum()),
        "errors_with_structure_features": int(
            errors["protein_structure_features_available"].astype(str).str.lower().eq("true").sum()
        ),
        "errors_with_numeric_foldx_in_frozen_matrix": int(errors["foldx_ddg_status"].eq("ok").sum()),
        "note": "These are model/held-out-label discordances for audit, not proven label errors.",
    }
    (OUT / "external_error_atlas_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
