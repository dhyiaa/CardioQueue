#!/usr/bin/env python3
"""Fit logistic recalibration on validation rows and evaluate fixed external predictions."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import brier_score_loss, log_loss


ROOT = Path(__file__).resolve().parents[3]
PREDICTIONS = (
    ROOT
    / "results/models/cardioqueue_v1_source_heldout"
    / "primary_binary_predictions_split_source_heldout.tsv"
)
OUT = Path(__file__).resolve().parent
PAPER_TABLES = ROOT / "paper_MS_SI/tables"
SEED = 20260917
N_BOOT = 5000


def probability_logit(probability: np.ndarray) -> np.ndarray:
    clipped = np.clip(probability, 1e-6, 1 - 1e-6)
    return np.log(clipped / (1 - clipped))


def apply_calibration(probability: np.ndarray, intercept: float, slope: float) -> np.ndarray:
    linear = intercept + slope * probability_logit(probability)
    return 1 / (1 + np.exp(-linear))


def score_rows(frame: pd.DataFrame, subset: str, model: str, column: str) -> dict[str, float | int | str]:
    y = frame["y_true"].to_numpy(dtype=int)
    probability = frame[column].to_numpy(dtype=float)
    calls = probability >= 0.5
    tp = int(((y == 1) & calls).sum())
    fp = int(((y == 0) & calls).sum())
    tn = int(((y == 0) & ~calls).sum())
    fn = int(((y == 1) & ~calls).sum())
    return {
        "subset": subset,
        "model": model,
        "rows": len(frame),
        "pathogenic": int(y.sum()),
        "benign": int((y == 0).sum()),
        "brier": brier_score_loss(y, probability),
        "log_loss": log_loss(y, probability),
        "mean_probability": probability.mean(),
        "sensitivity_0.5": tp / (tp + fn),
        "specificity_0.5": tn / (tn + fp),
    }


def differences(frame: pd.DataFrame, indices: np.ndarray) -> tuple[float, float]:
    y = frame["y_true"].to_numpy(dtype=int)[indices]
    raw = frame["raw_probability"].to_numpy(dtype=float)[indices]
    calibrated = frame["calibrated_probability"].to_numpy(dtype=float)[indices]
    return (
        brier_score_loss(y, calibrated) - brier_score_loss(y, raw),
        log_loss(y, calibrated) - log_loss(y, raw),
    )


def bootstrap(frame: pd.DataFrame, unit: str) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    samples = []
    if unit == "variant_stratified_by_outcome":
        y = frame["y_true"].to_numpy(dtype=int)
        groups = [np.flatnonzero(y == value) for value in (0, 1)]
        for _ in range(N_BOOT):
            indices = np.concatenate([rng.choice(group, len(group), replace=True) for group in groups])
            samples.append(differences(frame, indices))
    elif unit == "gene_cluster":
        genes = frame["primary_gene"].dropna().unique()
        gene_rows = {gene: np.flatnonzero(frame["primary_gene"].eq(gene).to_numpy()) for gene in genes}
        attempts = 0
        while len(samples) < N_BOOT and attempts < N_BOOT * 2:
            attempts += 1
            selected = rng.choice(genes, len(genes), replace=True)
            indices = np.concatenate([gene_rows[gene] for gene in selected])
            if frame["y_true"].iloc[indices].nunique() == 2:
                samples.append(differences(frame, indices))
        if len(samples) != N_BOOT:
            raise RuntimeError("Could not obtain the requested gene-cluster bootstrap samples")
    else:
        raise ValueError(f"Unknown resampling unit: {unit}")

    observed = differences(frame, np.arange(len(frame)))
    array = np.asarray(samples)
    return pd.DataFrame(
        [
            {
                "comparison": "validation_calibrated_minus_raw",
                "resampling_unit": unit,
                "metric": metric,
                "difference": observed[index],
                "ci_95_low": np.quantile(array[:, index], 0.025),
                "ci_95_high": np.quantile(array[:, index], 0.975),
                "bootstrap_samples": N_BOOT,
                "seed": SEED,
            }
            for index, metric in enumerate(["brier", "log_loss"])
        ]
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    PAPER_TABLES.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(PREDICTIONS, sep="\t")
    validation = frame[frame["split_source_heldout"].eq("validation")].copy()
    external = frame[frame["split_source_heldout"].str.startswith("external_")].copy()

    validation_probability = validation["pathogenic_probability"].to_numpy(dtype=float)
    design = sm.add_constant(probability_logit(validation_probability))
    fit = sm.GLM(validation["y_true"].to_numpy(dtype=int), design, family=sm.families.Binomial()).fit(disp=0)
    intercept, slope = (float(value) for value in fit.params)

    external["raw_probability"] = external["pathogenic_probability"].astype(float)
    external["calibrated_probability"] = apply_calibration(
        external["raw_probability"].to_numpy(), intercept, slope
    )
    subsets = {"external_all": external}
    subsets.update(
        {
            name.removeprefix("external_"): subset
            for name, subset in external.groupby("split_source_heldout", sort=True)
        }
    )
    metric_rows = []
    for subset_name, subset in subsets.items():
        metric_rows.append(score_rows(subset, subset_name, "Raw", "raw_probability"))
        metric_rows.append(
            score_rows(subset, subset_name, "Validation-calibrated", "calibrated_probability")
        )
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(OUT / "validation_calibration_metrics.tsv", sep="\t", index=False)
    metrics.to_csv(PAPER_TABLES / "Table_S11b_validation_calibration_metrics.tsv", sep="\t", index=False)

    intervals = pd.concat(
        [bootstrap(external, "variant_stratified_by_outcome"), bootstrap(external, "gene_cluster")],
        ignore_index=True,
    )
    intervals.to_csv(OUT / "validation_calibration_paired_bootstrap.tsv", sep="\t", index=False)
    intervals.to_csv(PAPER_TABLES / "Table_S11c_validation_calibration_paired_bootstrap.tsv", sep="\t", index=False)

    output = external[
        ["variant_id", "primary_gene", "split_source_heldout", "y_true", "raw_probability", "calibrated_probability"]
    ]
    output.to_csv(OUT / "external_validation_calibrated_predictions.tsv", sep="\t", index=False)
    manifest = {
        "prediction_file": str(PREDICTIONS.relative_to(ROOT)),
        "fit_split": "validation",
        "fit_rows": len(validation),
        "evaluation_rows": len(external),
        "intercept": intercept,
        "slope": slope,
        "bootstrap_samples": N_BOOT,
        "seed": SEED,
    }
    (OUT / "validation_calibration_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
