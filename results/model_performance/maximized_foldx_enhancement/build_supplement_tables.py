#!/usr/bin/env python3
"""Build manuscript-ready tables for the completed maximized FoldX analysis."""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
AUDIT = ROOT / "datasets/feature_sources/protein_structure/interim/maximized_foldx_mapping_audit.tsv"
RESULTS = ROOT / "results/model_performance/maximized_foldx_enhancement"
MODELS = ROOT / "results/models/maximized_foldx_2026"
TABLES = ROOT / "paper_MS_SI/tables"
TIERS = ("ge70", "ge50", "any")


def truthy(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def main() -> None:
    audit = pd.read_csv(AUDIT, sep="\t", low_memory=False)
    primary = audit[truthy(audit["training_slice_primary_binary"])].copy()
    coverage_rows = []
    for split, group in primary.groupby("split_source_heldout", sort=False):
        coverage_rows.append({
            "split": split,
            "cohort_rows": len(group),
            "mapped_any_plddt": int(group["foldx_eligible_any_plddt"].sum()),
            "plddt_ge50": int(group["foldx_eligible_plddt_ge50"].sum()),
            "plddt_ge70": int(group["foldx_eligible_plddt_ge70"].sum()),
            "plddt_ge90": int(group["foldx_eligible_plddt_ge90"].sum()),
            "not_exact_missense": int(group["audit_exclusion_category"].eq("not_an_exact_missense_substitution").sum()),
            "missense_mapping_unresolved": int(group["audit_exclusion_category"].eq("missense_without_unambiguous_reviewed_protein_mapping").sum()),
        })
    coverage = pd.DataFrame(coverage_rows)
    coverage.to_csv(TABLES / "Table_S17_foldx_coverage_eligibility.tsv", sep="\t", index=False)

    metrics = []
    paired = []
    for tier in TIERS:
        metric = pd.read_csv(RESULTS / tier / "model_metrics.tsv", sep="\t")
        metric.insert(0, "plddt_tier", tier)
        metrics.append(metric)
        differences = pd.read_csv(RESULTS / tier / "gene_cluster_paired_differences.tsv", sep="\t")
        differences.insert(0, "plddt_tier", tier)
        paired.append(differences)
    pd.concat(metrics, ignore_index=True).to_csv(
        TABLES / "Table_S17b_foldx_tier_metrics.tsv", sep="\t", index=False
    )
    pd.concat(paired, ignore_index=True).to_csv(
        TABLES / "Table_S17c_foldx_paired_bootstrap.tsv", sep="\t", index=False
    )

    importance_rows = []
    for tier in TIERS:
        for model in (f"eligibility_{tier}", f"ddg_{tier}"):
            frame = pd.read_csv(MODELS / model / "feature_importance.tsv", sep="\t")
            frame["overall_rank"] = range(1, len(frame) + 1)
            frame = frame[frame["feature"].str.startswith("maxfoldx_")].copy()
            frame.insert(0, "model", model)
            frame.insert(0, "plddt_tier", tier)
            importance_rows.append(frame)
    pd.concat(importance_rows, ignore_index=True).to_csv(
        TABLES / "Table_S17d_foldx_feature_importance.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
