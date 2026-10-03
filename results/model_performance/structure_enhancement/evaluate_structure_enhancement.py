#!/usr/bin/env python3
"""Compare frozen baseline and registry-wide structure-enhanced models."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, confusion_matrix, roc_auc_score


ROOT = Path(__file__).resolve().parents[3]
BASELINE = ROOT / "results/models/primary_binary_catboost_v0_strict_source_heldout_hiro_emerge_rescued"
ENHANCED = ROOT / "results/models/primary_binary_catboost_structure_enhanced_source_heldout"
STRUCTURE = ROOT / "datasets/feature_sources/protein_structure/interim/registry_structure_features.tsv"
OUT = Path(__file__).resolve().parent
PAPER_TABLES = ROOT / "paper_MS_SI/tables"
SEED = 20260917
N_BOOT = 5000


def load_predictions(directory: Path, name: str) -> pd.DataFrame:
    frame = pd.read_csv(directory / "primary_binary_predictions_split_source_heldout.tsv", sep="\t")
    frame = frame[frame["split_source_heldout"].str.startswith("external_")].copy()
    return frame[["variant_id", "primary_gene", "split_source_heldout", "y_true", "pathogenic_probability"]].rename(
        columns={"pathogenic_probability": name}
    )


def metrics(y: np.ndarray, score: np.ndarray) -> dict[str, float | int]:
    tn, fp, fn, tp = confusion_matrix(y, score >= 0.5, labels=[0, 1]).ravel()
    return {
        "rows": len(y),
        "pathogenic": int(y.sum()),
        "benign": int((y == 0).sum()),
        "auroc": roc_auc_score(y, score),
        "auprc": average_precision_score(y, score),
        "brier": brier_score_loss(y, score),
        "sensitivity_0.5": tp / (tp + fn),
        "specificity_0.5": tn / (tn + fp),
    }


def paired_bootstrap(frame: pd.DataFrame, subset: str) -> pd.DataFrame:
    y = frame["y_true"].to_numpy(dtype=int)
    baseline = frame["baseline_probability"].to_numpy(dtype=float)
    enhanced = frame["structure_probability"].to_numpy(dtype=float)
    groups = [np.flatnonzero(y == value) for value in (0, 1)]
    rng = np.random.default_rng(SEED)
    values = {"auroc": [], "auprc": [], "brier": []}
    for _ in range(N_BOOT):
        idx = np.concatenate([rng.choice(group, len(group), replace=True) for group in groups])
        yb = y[idx]
        values["auroc"].append(roc_auc_score(yb, enhanced[idx]) - roc_auc_score(yb, baseline[idx]))
        values["auprc"].append(
            average_precision_score(yb, enhanced[idx]) - average_precision_score(yb, baseline[idx])
        )
        values["brier"].append(brier_score_loss(yb, enhanced[idx]) - brier_score_loss(yb, baseline[idx]))

    observed = {
        "auroc": roc_auc_score(y, enhanced) - roc_auc_score(y, baseline),
        "auprc": average_precision_score(y, enhanced) - average_precision_score(y, baseline),
        "brier": brier_score_loss(y, enhanced) - brier_score_loss(y, baseline),
    }
    rows = []
    for metric, samples in values.items():
        array = np.asarray(samples)
        rows.append(
            {
                "subset": subset,
                "resampling_unit": "variant_stratified_by_outcome",
                "metric": metric,
                "structure_minus_baseline": observed[metric],
                "ci_95_low": np.quantile(array, 0.025),
                "ci_95_high": np.quantile(array, 0.975),
                "bootstrap_samples": N_BOOT,
                "seed": SEED,
            }
        )
    return pd.DataFrame(rows)


def gene_cluster_bootstrap(frame: pd.DataFrame, subset: str) -> pd.DataFrame:
    genes = frame["primary_gene"].dropna().unique()
    gene_rows = {gene: np.flatnonzero(frame["primary_gene"].eq(gene).to_numpy()) for gene in genes}
    y = frame["y_true"].to_numpy(dtype=int)
    baseline = frame["baseline_probability"].to_numpy(dtype=float)
    enhanced = frame["structure_probability"].to_numpy(dtype=float)
    rng = np.random.default_rng(SEED)
    values = {"auroc": [], "auprc": [], "brier": []}
    attempts = 0
    while len(values["auroc"]) < N_BOOT and attempts < N_BOOT * 2:
        attempts += 1
        sampled_genes = rng.choice(genes, len(genes), replace=True)
        idx = np.concatenate([gene_rows[gene] for gene in sampled_genes])
        yb = y[idx]
        if len(np.unique(yb)) != 2:
            continue
        values["auroc"].append(roc_auc_score(yb, enhanced[idx]) - roc_auc_score(yb, baseline[idx]))
        values["auprc"].append(
            average_precision_score(yb, enhanced[idx]) - average_precision_score(yb, baseline[idx])
        )
        values["brier"].append(brier_score_loss(yb, enhanced[idx]) - brier_score_loss(yb, baseline[idx]))
    if len(values["auroc"]) != N_BOOT:
        raise RuntimeError("Could not obtain the requested number of valid gene-cluster bootstrap samples")

    observed = {
        "auroc": roc_auc_score(y, enhanced) - roc_auc_score(y, baseline),
        "auprc": average_precision_score(y, enhanced) - average_precision_score(y, baseline),
        "brier": brier_score_loss(y, enhanced) - brier_score_loss(y, baseline),
    }
    rows = []
    for metric, samples in values.items():
        array = np.asarray(samples)
        rows.append(
            {
                "subset": subset,
                "resampling_unit": "gene_cluster",
                "metric": metric,
                "structure_minus_baseline": observed[metric],
                "ci_95_low": np.quantile(array, 0.025),
                "ci_95_high": np.quantile(array, 0.975),
                "bootstrap_samples": N_BOOT,
                "seed": SEED,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    PAPER_TABLES.mkdir(parents=True, exist_ok=True)
    baseline = load_predictions(BASELINE, "baseline_probability")
    enhanced = load_predictions(ENHANCED, "structure_probability")
    paired = baseline.merge(
        enhanced,
        on=["variant_id", "primary_gene", "split_source_heldout", "y_true"],
        how="inner",
        validate="one_to_one",
    )
    if len(paired) != 430:
        raise ValueError(f"Expected 430 paired external rows, found {len(paired)}")
    structure = pd.read_csv(STRUCTURE, sep="\t", usecols=["variant_id", "structure_available"])
    paired = paired.merge(structure, on="variant_id", validate="one_to_one")

    metric_rows = []
    subsets = {
        "external_all": paired,
        "external_structure_covered": paired[paired["structure_available"].eq(True)],
        "external_structure_missing": paired[paired["structure_available"].ne(True)],
        "external_hiro": paired[paired["split_source_heldout"].eq("external_hiro")],
        "external_emerge": paired[paired["split_source_heldout"].eq("external_emerge")],
        "external_cardioboost": paired[paired["split_source_heldout"].eq("external_cardioboost")],
    }
    for subset, frame in subsets.items():
        for model, column in [
            ("Baseline", "baseline_probability"),
            ("Structure-enhanced", "structure_probability"),
        ]:
            row = {"subset": subset, "model": model}
            row.update(metrics(frame["y_true"].to_numpy(dtype=int), frame[column].to_numpy(dtype=float)))
            metric_rows.append(row)
    metrics_frame = pd.DataFrame(metric_rows)
    metrics_frame.to_csv(OUT / "structure_enhancement_metrics.tsv", sep="\t", index=False)
    metrics_frame.to_csv(PAPER_TABLES / "Table_S15_structure_enhancement_metrics.tsv", sep="\t", index=False)

    differences = pd.concat(
        [
            paired_bootstrap(subsets["external_all"], "external_all"),
            paired_bootstrap(subsets["external_structure_covered"], "external_structure_covered"),
            gene_cluster_bootstrap(subsets["external_all"], "external_all"),
            gene_cluster_bootstrap(subsets["external_structure_covered"], "external_structure_covered"),
        ],
        ignore_index=True,
    )
    differences.to_csv(OUT / "structure_enhancement_paired_bootstrap.tsv", sep="\t", index=False)
    differences.to_csv(PAPER_TABLES / "Table_S15b_structure_enhancement_paired_bootstrap.tsv", sep="\t", index=False)

    importance = pd.read_csv(ENHANCED / "feature_importance.tsv", sep="\t")
    structure_importance = importance[importance["feature"].str.startswith("structure_")].copy()
    structure_importance.to_csv(OUT / "structure_feature_importance.tsv", sep="\t", index=False)
    structure_importance.to_csv(PAPER_TABLES / "Table_S15c_structure_feature_importance.tsv", sep="\t", index=False)

    all_metrics = metrics_frame[metrics_frame["subset"].eq("external_all")].set_index("model")
    covered_metrics = metrics_frame[metrics_frame["subset"].eq("external_structure_covered")].set_index("model")
    gene_differences = differences[differences["resampling_unit"].eq("gene_cluster")]
    all_delta = gene_differences[gene_differences["subset"].eq("external_all")].set_index("metric")
    covered_delta = gene_differences[
        gene_differences["subset"].eq("external_structure_covered")
    ].set_index("metric")
    report = f"""# Registry-wide protein-structure enhancement

## Coverage

Reference-validated AlphaFold v6, DSSP, and FreeSASA features were available for 2,889 training, 493 validation, and 271 source-held-out variants. Long proteins used overlapping 1,400-residue AlphaFold fragments at 200-residue offsets. The selected fragment maximized distance from its nearest edge. Reference amino acids were checked against both reviewed UniProt sequence and the selected PDB fragment.

## External performance

Across all 430 held-out variants, baseline AUROC was {all_metrics.loc['Baseline', 'auroc']:.3f} and structure-enhanced AUROC was {all_metrics.loc['Structure-enhanced', 'auroc']:.3f}. The gene-clustered paired difference was {all_delta.loc['auroc', 'structure_minus_baseline']:.3f} (95% CI {all_delta.loc['auroc', 'ci_95_low']:.3f} to {all_delta.loc['auroc', 'ci_95_high']:.3f}). Brier score changed from {all_metrics.loc['Baseline', 'brier']:.4f} to {all_metrics.loc['Structure-enhanced', 'brier']:.4f}.

Among 271 structure-covered external variants, AUROC changed from {covered_metrics.loc['Baseline', 'auroc']:.3f} to {covered_metrics.loc['Structure-enhanced', 'auroc']:.3f}. The gene-clustered paired difference was {covered_delta.loc['auroc', 'structure_minus_baseline']:.3f} (95% CI {covered_delta.loc['auroc', 'ci_95_low']:.3f} to {covered_delta.loc['auroc', 'ci_95_high']:.3f}).

## Interpretation

The enhanced model assigned nonzero importance to pLDDT, DSSP, and FreeSASA values. The availability indicator had zero importance. Effect estimates are small and their paired intervals should determine whether the analysis is described as improvement or feasibility. FoldX DDG remains unavailable in development and is not part of the demonstrated incremental signal.
"""
    (OUT / "STRUCTURE_ENHANCEMENT_REPORT.md").write_text(report)
    manifest = {
        "baseline_model": str(BASELINE.relative_to(ROOT)),
        "enhanced_model": str(ENHANCED.relative_to(ROOT)),
        "structure_features": str(STRUCTURE.relative_to(ROOT)),
        "bootstrap_samples": N_BOOT,
        "seed": SEED,
        "external_rows": len(paired),
        "structure_covered_external_rows": int(paired["structure_available"].eq(True).sum()),
    }
    (OUT / "structure_enhancement_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
