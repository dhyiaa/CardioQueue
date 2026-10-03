#!/usr/bin/env python3
"""Compare baseline and structure-enhanced models on gene-stress test rows."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, confusion_matrix, roc_auc_score


ROOT = Path(__file__).resolve().parents[3]
BASELINE = ROOT / "results/models/primary_binary_catboost_v0_strict_gene_stress_hiro_emerge_rescued"
ENHANCED = ROOT / "results/models/primary_binary_catboost_structure_enhanced_gene_stress"
STRUCTURE = ROOT / "datasets/feature_sources/protein_structure/interim/registry_structure_features.tsv"
OUT = Path(__file__).resolve().parent
PAPER_TABLES = ROOT / "paper_MS_SI/tables"
SPLIT = "sparse_gene_stress_test"
SEED = 20260917
N_BOOT = 5000


def load_predictions(directory: Path, name: str) -> pd.DataFrame:
    path = directory / "primary_binary_predictions_split_gene_stress.tsv"
    frame = pd.read_csv(path, sep="\t")
    frame = frame[frame["split_gene_stress"].eq(SPLIT)].copy()
    columns = ["variant_id", "primary_gene", "split_gene_stress", "y_true", "pathogenic_probability"]
    return frame[columns].rename(columns={"pathogenic_probability": name})


def metrics(frame: pd.DataFrame, score_column: str) -> dict[str, float | int]:
    y = frame["y_true"].to_numpy(dtype=int)
    score = frame[score_column].to_numpy(dtype=float)
    tn, fp, fn, tp = confusion_matrix(y, score >= 0.5, labels=[0, 1]).ravel()
    return {
        "rows": len(y),
        "genes": frame["primary_gene"].nunique(),
        "pathogenic": int(y.sum()),
        "benign": int((y == 0).sum()),
        "auroc": roc_auc_score(y, score),
        "auprc": average_precision_score(y, score),
        "brier": brier_score_loss(y, score),
        "sensitivity_0.5": tp / (tp + fn),
        "specificity_0.5": tn / (tn + fp),
        "precision_0.5": tp / (tp + fp),
    }


def metric_differences(frame: pd.DataFrame, indices: np.ndarray) -> tuple[float, float, float]:
    y = frame["y_true"].to_numpy(dtype=int)[indices]
    baseline = frame["baseline_probability"].to_numpy(dtype=float)[indices]
    enhanced = frame["structure_probability"].to_numpy(dtype=float)[indices]
    return (
        roc_auc_score(y, enhanced) - roc_auc_score(y, baseline),
        average_precision_score(y, enhanced) - average_precision_score(y, baseline),
        brier_score_loss(y, enhanced) - brier_score_loss(y, baseline),
    )


def paired_bootstrap(frame: pd.DataFrame, subset: str, unit: str) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    samples: list[tuple[float, float, float]] = []
    if unit == "variant_stratified_by_outcome":
        y = frame["y_true"].to_numpy(dtype=int)
        groups = [np.flatnonzero(y == value) for value in (0, 1)]
        for _ in range(N_BOOT):
            indices = np.concatenate([rng.choice(group, len(group), replace=True) for group in groups])
            samples.append(metric_differences(frame, indices))
    elif unit == "gene_cluster":
        genes = frame["primary_gene"].dropna().unique()
        gene_rows = {gene: np.flatnonzero(frame["primary_gene"].eq(gene).to_numpy()) for gene in genes}
        attempts = 0
        while len(samples) < N_BOOT and attempts < N_BOOT * 3:
            attempts += 1
            selected = rng.choice(genes, len(genes), replace=True)
            indices = np.concatenate([gene_rows[gene] for gene in selected])
            if frame["y_true"].iloc[indices].nunique() == 2:
                samples.append(metric_differences(frame, indices))
        if len(samples) != N_BOOT:
            raise RuntimeError("Could not obtain the requested gene-cluster bootstrap samples")
    else:
        raise ValueError(f"Unknown bootstrap unit: {unit}")

    observed = metric_differences(frame, np.arange(len(frame)))
    array = np.asarray(samples)
    rows = []
    for index, metric in enumerate(["auroc", "auprc", "brier"]):
        rows.append(
            {
                "subset": subset,
                "resampling_unit": unit,
                "metric": metric,
                "structure_minus_baseline": observed[index],
                "ci_95_low": np.quantile(array[:, index], 0.025),
                "ci_95_high": np.quantile(array[:, index], 0.975),
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
        on=["variant_id", "primary_gene", "split_gene_stress", "y_true"],
        how="inner",
        validate="one_to_one",
    )
    if len(paired) != 5890:
        raise ValueError(f"Expected 5,890 paired gene-stress rows, found {len(paired):,}")
    structure = pd.read_csv(STRUCTURE, sep="\t", usecols=["variant_id", "structure_available"])
    paired = paired.merge(structure, on="variant_id", validate="one_to_one")

    subsets = {
        "heldout_genes_all": paired,
        "heldout_genes_structure_covered": paired[paired["structure_available"].eq(True)].copy(),
        "heldout_genes_structure_missing": paired[paired["structure_available"].ne(True)].copy(),
    }
    metric_rows = []
    for subset, frame in subsets.items():
        for model, column in [
            ("Baseline", "baseline_probability"),
            ("Structure-enhanced", "structure_probability"),
        ]:
            row = {"subset": subset, "model": model}
            row.update(metrics(frame, column))
            metric_rows.append(row)
    metrics_frame = pd.DataFrame(metric_rows)
    metrics_frame.to_csv(OUT / "structure_gene_stress_metrics.tsv", sep="\t", index=False)
    metrics_frame.to_csv(PAPER_TABLES / "Table_S15d_structure_gene_stress_metrics.tsv", sep="\t", index=False)

    differences = pd.concat(
        [
            paired_bootstrap(frame, subset, unit)
            for subset, frame in subsets.items()
            if frame["y_true"].nunique() == 2
            for unit in ["variant_stratified_by_outcome", "gene_cluster"]
        ],
        ignore_index=True,
    )
    differences.to_csv(OUT / "structure_gene_stress_paired_bootstrap.tsv", sep="\t", index=False)
    differences.to_csv(PAPER_TABLES / "Table_S15e_structure_gene_stress_paired_bootstrap.tsv", sep="\t", index=False)

    all_metrics = metrics_frame[metrics_frame["subset"].eq("heldout_genes_all")].set_index("model")
    covered_metrics = metrics_frame[
        metrics_frame["subset"].eq("heldout_genes_structure_covered")
    ].set_index("model")
    gene_differences = differences[differences["resampling_unit"].eq("gene_cluster")]
    all_differences = gene_differences[gene_differences["subset"].eq("heldout_genes_all")].set_index(
        "metric"
    )
    covered_differences = gene_differences[
        gene_differences["subset"].eq("heldout_genes_structure_covered")
    ].set_index("metric")
    report = f"""# Structure enhancement under gene holdout

The gene-stress test contained {len(paired):,} variants from {paired['primary_gene'].nunique()} genes absent from training. Structure fields covered {int(paired['structure_available'].eq(True).sum())} rows across {subsets['heldout_genes_structure_covered']['primary_gene'].nunique()} genes.

Across all held-out-gene rows, baseline AUROC was {all_metrics.loc['Baseline', 'auroc']:.5f} and structure-enhanced AUROC was {all_metrics.loc['Structure-enhanced', 'auroc']:.5f}. The gene-clustered difference was {all_differences.loc['auroc', 'structure_minus_baseline']:.5f} (95% CI {all_differences.loc['auroc', 'ci_95_low']:.5f} to {all_differences.loc['auroc', 'ci_95_high']:.5f}). AUPRC changed from {all_metrics.loc['Baseline', 'auprc']:.5f} to {all_metrics.loc['Structure-enhanced', 'auprc']:.5f}.

Among structure-covered rows, AUROC changed from {covered_metrics.loc['Baseline', 'auroc']:.5f} to {covered_metrics.loc['Structure-enhanced', 'auroc']:.5f}. The gene-clustered difference was {covered_differences.loc['auroc', 'structure_minus_baseline']:.5f} (95% CI {covered_differences.loc['auroc', 'ci_95_low']:.5f} to {covered_differences.loc['auroc', 'ci_95_high']:.5f}). The analysis does not establish improved discrimination for unseen genes.

Predictions can change on structure-missing rows because the enhanced model is refitted in full. Those changes cannot be attributed to direct structural measurements.
"""
    (OUT / "STRUCTURE_GENE_STRESS_REPORT.md").write_text(report)

    manifest = {
        "baseline_model": str(BASELINE.relative_to(ROOT)),
        "enhanced_model": str(ENHANCED.relative_to(ROOT)),
        "structure_features": str(STRUCTURE.relative_to(ROOT)),
        "split": SPLIT,
        "rows": len(paired),
        "genes": paired["primary_gene"].nunique(),
        "structure_covered_rows": int(paired["structure_available"].eq(True).sum()),
        "bootstrap_samples": N_BOOT,
        "seed": SEED,
    }
    (OUT / "structure_gene_stress_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
