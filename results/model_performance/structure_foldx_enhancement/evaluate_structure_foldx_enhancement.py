#!/usr/bin/env python3
"""Compare baseline, residue-context, and residue-context plus FoldX models."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


ROOT = Path(__file__).resolve().parents[3]
BASELINE = ROOT / "results/models/primary_binary_catboost_v0_strict_source_heldout_hiro_emerge_rescued"
STRUCTURE = ROOT / "results/models/primary_binary_catboost_structure_enhanced_source_heldout"
FOLDX_MODEL = ROOT / "results/models/primary_binary_catboost_structure_foldx_source_heldout"
DDG_MODEL = ROOT / "results/models/primary_binary_catboost_foldx_only_source_heldout"
AUGMENT = ROOT / "datasets/feature_sources/protein_structure/interim/registry_structure_foldx_features.tsv"
OUT = Path(__file__).resolve().parent
PAPER_TABLES = ROOT / "paper_MS_SI/tables"
SEED = 20260918
N_BOOT = 5000


def load_predictions(directory: Path, score_name: str) -> pd.DataFrame:
    frame = pd.read_csv(directory / "primary_binary_predictions_split_source_heldout.tsv", sep="\t")
    frame = frame[frame["split_source_heldout"].astype(str).str.startswith("external_")].copy()
    return frame[["variant_id", "primary_gene", "split_source_heldout", "y_true", "pathogenic_probability"]].rename(
        columns={"pathogenic_probability": score_name}
    )


def metrics(y: np.ndarray, score: np.ndarray) -> dict[str, float | int]:
    return {
        "rows": len(y),
        "pathogenic": int(y.sum()),
        "benign": int((y == 0).sum()),
        "auroc": float(roc_auc_score(y, score)),
        "auprc": float(average_precision_score(y, score)),
        "brier": float(brier_score_loss(y, score)),
    }


def paired_bootstrap(
    frame: pd.DataFrame,
    subset: str,
    resampling_unit: str,
    reference_column: str,
    comparison_column: str,
    comparison: str,
) -> pd.DataFrame:
    y = frame["y_true"].to_numpy(dtype=int)
    old = frame[reference_column].to_numpy(dtype=float)
    new = frame[comparison_column].to_numpy(dtype=float)
    rng = np.random.default_rng(SEED)
    samples = {"auroc": [], "auprc": [], "brier": []}
    if resampling_unit == "gene_cluster":
        units = frame["primary_gene"].dropna().unique()
        rows = {gene: np.flatnonzero(frame["primary_gene"].eq(gene).to_numpy()) for gene in units}
    else:
        units = None
        class_rows = [np.flatnonzero(y == value) for value in (0, 1)]

    attempts = 0
    while len(samples["auroc"]) < N_BOOT and attempts < N_BOOT * 3:
        attempts += 1
        if units is not None:
            selected = rng.choice(units, len(units), replace=True)
            idx = np.concatenate([rows[unit] for unit in selected])
        else:
            idx = np.concatenate([rng.choice(group, len(group), replace=True) for group in class_rows])
        if len(np.unique(y[idx])) != 2:
            continue
        samples["auroc"].append(roc_auc_score(y[idx], new[idx]) - roc_auc_score(y[idx], old[idx]))
        samples["auprc"].append(
            average_precision_score(y[idx], new[idx]) - average_precision_score(y[idx], old[idx])
        )
        samples["brier"].append(brier_score_loss(y[idx], new[idx]) - brier_score_loss(y[idx], old[idx]))
    if len(samples["auroc"]) != N_BOOT:
        raise RuntimeError(f"Unable to obtain {N_BOOT} valid {resampling_unit} bootstrap samples")

    observed = {
        "auroc": roc_auc_score(y, new) - roc_auc_score(y, old),
        "auprc": average_precision_score(y, new) - average_precision_score(y, old),
        "brier": brier_score_loss(y, new) - brier_score_loss(y, old),
    }
    return pd.DataFrame(
        [
            {
                "subset": subset,
                "resampling_unit": resampling_unit,
                "comparison": comparison,
                "metric": metric,
                "metric_difference": observed[metric],
                "ci_95_low": float(np.quantile(values, 0.025)),
                "ci_95_high": float(np.quantile(values, 0.975)),
                "bootstrap_samples": N_BOOT,
                "seed": SEED,
            }
            for metric, values in samples.items()
        ]
    )


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    PAPER_TABLES.mkdir(parents=True, exist_ok=True)
    frame = load_predictions(BASELINE, "baseline_probability")
    frame = frame.merge(
        load_predictions(STRUCTURE, "structure_probability"),
        on=["variant_id", "primary_gene", "split_source_heldout", "y_true"],
        validate="one_to_one",
    ).merge(
        load_predictions(DDG_MODEL, "ddg_probability"),
        on=["variant_id", "primary_gene", "split_source_heldout", "y_true"],
        validate="one_to_one",
    ).merge(
        load_predictions(FOLDX_MODEL, "foldx_probability"),
        on=["variant_id", "primary_gene", "split_source_heldout", "y_true"],
        validate="one_to_one",
    )
    augment = pd.read_csv(AUGMENT, sep="\t", usecols=["variant_id", "structure_foldx_available"])
    frame = frame.merge(augment, on="variant_id", validate="one_to_one")
    if len(frame) != 430:
        raise ValueError(f"Expected 430 paired held-out variants, found {len(frame)}")

    subsets = {
        "external_all": frame,
        "external_foldx_covered": frame[frame["structure_foldx_available"].eq(True)],
    }
    metric_rows = []
    for subset, data in subsets.items():
        for model, column in [
            ("Baseline", "baseline_probability"),
            ("Structure", "structure_probability"),
            ("FoldX DDG", "ddg_probability"),
            ("Structure+FoldX", "foldx_probability"),
        ]:
            row = {"subset": subset, "model": model}
            row.update(metrics(data["y_true"].to_numpy(dtype=int), data[column].to_numpy(dtype=float)))
            metric_rows.append(row)
    metric_frame = pd.DataFrame(metric_rows)
    metric_frame.to_csv(OUT / "structure_foldx_metrics.tsv", sep="\t", index=False)
    metric_frame.to_csv(PAPER_TABLES / "Table_S17_structure_foldx_metrics.tsv", sep="\t", index=False)

    comparisons = [
        ("baseline_probability", "structure_probability", "Structure minus baseline"),
        ("baseline_probability", "ddg_probability", "FoldX DDG minus baseline"),
        ("baseline_probability", "foldx_probability", "Structure+FoldX minus baseline"),
        ("structure_probability", "foldx_probability", "Structure+FoldX minus structure"),
    ]
    differences = pd.concat(
        [
            paired_bootstrap(data, subset, unit, reference, candidate, label)
            for subset, data in subsets.items()
            for unit in ("variant_stratified_by_outcome", "gene_cluster")
            for reference, candidate, label in comparisons
        ],
        ignore_index=True,
    )
    differences.to_csv(OUT / "structure_foldx_paired_bootstrap.tsv", sep="\t", index=False)
    differences.to_csv(PAPER_TABLES / "Table_S17b_structure_foldx_paired_bootstrap.tsv", sep="\t", index=False)

    importance = pd.read_csv(FOLDX_MODEL / "feature_importance.tsv", sep="\t")
    importance = importance[importance["feature"].str.startswith("structure_")].copy()
    importance.to_csv(OUT / "structure_foldx_feature_importance.tsv", sep="\t", index=False)
    importance.to_csv(PAPER_TABLES / "Table_S17c_structure_foldx_feature_importance.tsv", sep="\t", index=False)

    manifest = {
        "baseline_model": str(BASELINE.relative_to(ROOT)),
        "structure_model": str(STRUCTURE.relative_to(ROOT)),
        "foldx_model": str(FOLDX_MODEL.relative_to(ROOT)),
        "ddg_only_model": str(DDG_MODEL.relative_to(ROOT)),
        "external_rows": len(frame),
        "external_foldx_covered": int(frame["structure_foldx_available"].sum()),
        "bootstrap_samples": N_BOOT,
        "seed": SEED,
    }
    (OUT / "structure_foldx_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(metric_frame.to_string(index=False))
    print(differences.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
