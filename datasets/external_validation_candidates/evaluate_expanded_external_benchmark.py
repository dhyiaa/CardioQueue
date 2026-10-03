#!/usr/bin/env python3
"""Evaluate expanded external predictions by source with gene-cluster bootstrap CIs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, confusion_matrix, roc_auc_score


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = ROOT / "datasets/external_validation_candidates/expanded_external_2026"
DEFAULT_PREDICTIONS = ROOT / (
    "results/models/cardioqueue_v1_expanded_external_2026_public300/"
    "primary_binary_predictions_split_source_heldout.tsv"
)
DEFAULT_OUT = ROOT / "results/model_performance/expanded_external_2026"


def metrics(df: pd.DataFrame, include_calibration: bool = True) -> dict[str, float]:
    y = df["y_true"].astype(int).to_numpy()
    p = df["pathogenic_probability"].astype(float).to_numpy()
    pred = (p >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    output = {
        "rows": len(df), "pathogenic": int(y.sum()), "benign": int((y == 0).sum()),
        "auroc": roc_auc_score(y, p), "auprc": average_precision_score(y, p),
        "brier": brier_score_loss(y, p),
        "sensitivity": tp / (tp + fn), "specificity": tn / (tn + fp),
        "ppv": tp / (tp + fp), "npv": tn / (tn + fn),
        "tp": int(tp), "fp": int(fp), "tn": int(tn), "fn": int(fn),
    }
    if include_calibration:
        clipped = np.clip(p, 1e-6, 1 - 1e-6)
        logits = np.log(clipped / (1 - clipped)).reshape(-1, 1)
        calibration = LogisticRegression(C=1e6, solver="lbfgs").fit(logits, y)
        output["calibration_intercept"] = float(calibration.intercept_[0])
        output["calibration_slope"] = float(calibration.coef_[0, 0])
    return output


def gene_cluster_bootstrap(df: pd.DataFrame, samples: int, seed: int) -> dict[str, tuple[float, float]]:
    rng = np.random.default_rng(seed)
    genes = df["gene"].fillna("UNKNOWN").astype(str).unique()
    by_gene = {gene: df[df["gene"].fillna("UNKNOWN").astype(str) == gene] for gene in genes}
    values: dict[str, list[float]] = {
        key: [] for key in ["auroc", "auprc", "brier", "sensitivity", "specificity", "ppv", "npv"]
    }
    for _ in range(samples):
        sampled = rng.choice(genes, size=len(genes), replace=True)
        boot = pd.concat([by_gene[gene] for gene in sampled], ignore_index=True)
        if boot["y_true"].nunique() < 2:
            continue
        point = metrics(boot, include_calibration=False)
        for key in values:
            values[key].append(point[key])
    return {
        key: (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))
        for key, vals in values.items() if vals
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--observations", type=Path, default=DEFAULT_DATA / "external_source_observations.tsv")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_DATA / "external_unique_variant_manifest.tsv")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_DATA / "modeling_table_expanded_external_annotated.tsv")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--bootstrap-samples", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260925)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    predictions = pd.read_csv(args.predictions, sep="\t")
    predictions = predictions[predictions["split_source_heldout"] == "external_expanded"].copy()
    manifest = pd.read_csv(args.manifest, sep="\t")
    observations = pd.read_csv(args.observations, sep="\t", low_memory=False)
    matrix = pd.read_csv(args.matrix, sep="\t", low_memory=False)
    benchmark_ids = set(manifest["variant_id"])
    if set(predictions["variant_id"]) != benchmark_ids:
        raise ValueError("Prediction IDs do not exactly equal the frozen unique benchmark manifest")

    combined = predictions.merge(
        manifest[["variant_id", "gene", "combined_label", "external_sources", "already_in_modeling_matrix"]],
        on="variant_id", how="inner", validate="one_to_one",
    )
    combined["y_true"] = combined["combined_label"].map({"Benign": 0, "Pathogenic": 1})
    combined["analysis_subset"] = "combined_unique"

    source = observations[
        observations["variant_id"].isin(benchmark_ids)
        & observations["reference_match"].fillna(False)
        & observations["in_gene_panel"].fillna(False)
    ].merge(
        predictions[["variant_id", "pathogenic_probability"]],
        on="variant_id", how="inner", validate="many_to_one",
    )
    source["y_true"] = source["external_binary_label"].astype(int)
    source["analysis_subset"] = source["source"]

    subsets = [("combined_unique", combined)]
    subsets.extend((name, group.copy()) for name, group in source.groupby("source", sort=True))
    subsets.extend(
        [
            ("existing_matrix_external", combined[combined["already_in_modeling_matrix"]]),
            ("newly_annotated_external", combined[~combined["already_in_modeling_matrix"]]),
        ]
    )
    point_rows = []
    ci_rows = []
    for offset, (name, frame) in enumerate(subsets):
        point = metrics(frame)
        point_rows.append({"subset": name, **point})
        intervals = gene_cluster_bootstrap(frame, args.bootstrap_samples, args.seed + offset)
        for metric, (low, high) in intervals.items():
            ci_rows.append(
                {"subset": name, "metric": metric, "estimate": point[metric], "ci_95_low": low,
                 "ci_95_high": high, "bootstrap_samples": args.bootstrap_samples,
                 "bootstrap_unit": "gene"}
            )

    feature_rows = matrix[matrix["variant_id"].isin(benchmark_ids)].copy()
    coverage_columns = {
        "dbnsfp": (feature_rows["dbnsfp_status"] == "ok"),
        "alphamissense": (feature_rows["alphamissense_direct_status"] == "ok"),
        "vep": (feature_rows["vep_annotation_status"] == "ok"),
        "spliceai": (feature_rows["spliceai_status"] == "ok"),
        "protein": (feature_rows["protein_feature_status"] == "ok"),
        "alphafold_plddt": (feature_rows["alphafold_plddt_status"] == "ok"),
        "foldx_ddg": (feature_rows["foldx_ddg_status"] == "ok"),
        "gnomad_observed": (feature_rows["gnomad_final_status"] == "observed"),
    }
    coverage = []
    for feature, present in coverage_columns.items():
        coverage.append({"feature_group": feature, "present": int(present.sum()), "total": len(feature_rows), "fraction": float(present.mean())})

    point_table = pd.DataFrame(point_rows)
    ci_table = pd.DataFrame(ci_rows)
    point_table.to_csv(args.out_dir / "performance_by_source.tsv", sep="\t", index=False)
    ci_table.to_csv(args.out_dir / "gene_cluster_bootstrap_ci.tsv", sep="\t", index=False)
    pd.DataFrame(coverage).to_csv(args.out_dir / "feature_coverage.tsv", sep="\t", index=False)
    combined.to_csv(args.out_dir / "combined_unique_predictions.tsv", sep="\t", index=False)
    source.to_csv(args.out_dir / "source_observation_predictions.tsv", sep="\t", index=False)
    summary = {
        "unique_external_variants": len(combined),
        "source_observations": len(source),
        "bootstrap": {"unit": "gene", "samples": args.bootstrap_samples, "seed": args.seed},
        "combined": point_rows[0],
        "warnings": [
            "These are retrospective source-held-out cohorts, not prospective clinical validation.",
            "Some external labels draw on public evidence ecosystems that overlap ClinVar conceptually.",
            "Source-specific rows are observations; the combined analysis is deduplicated by GRCh38 allele.",
        ],
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(point_table.to_string(index=False))


if __name__ == "__main__":
    main()
