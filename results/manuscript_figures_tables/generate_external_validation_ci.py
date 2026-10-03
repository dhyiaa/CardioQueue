#!/usr/bin/env python3
"""Bootstrap confidence intervals for source-held-out manuscript metrics."""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


ROOT = Path(__file__).resolve().parents[2]
PREDICTIONS = ROOT / (
    "results/models/cardioqueue_v1_source_heldout/"
    "primary_binary_predictions_split_source_heldout.tsv"
)
OUT = ROOT / "results/manuscript_figures_tables/tables/external_validation_bootstrap_ci.tsv"
CLUSTER_OUT = ROOT / "results/manuscript_figures_tables/tables/external_validation_gene_cluster_bootstrap_ci.tsv"
MANIFEST = ROOT / "results/manuscript_figures_tables/external_validation_bootstrap_manifest.json"

SUBSETS = {
    "external_all": ["external_hiro", "external_emerge", "external_cardioboost"],
    "external_hiro": ["external_hiro"],
    "external_emerge": ["external_emerge"],
    "external_cardioboost": ["external_cardioboost"],
}


def safe_div(numerator: int, denominator: int) -> float:
    return float(numerator / denominator) if denominator else float("nan")


def calibration_parameters(y: np.ndarray, p: np.ndarray) -> tuple[float, float]:
    clipped = np.clip(p, 1e-6, 1 - 1e-6)
    logit = np.log(clipped / (1 - clipped))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = sm.GLM(y, sm.add_constant(logit), family=sm.families.Binomial()).fit(disp=0)
    return float(fit.params[0]), float(fit.params[1])


def metrics(y: np.ndarray, p: np.ndarray) -> dict[str, float]:
    binary = p >= 0.5
    tp = int(((y == 1) & binary).sum())
    fp = int(((y == 0) & binary).sum())
    tn = int(((y == 0) & ~binary).sum())
    fn = int(((y == 1) & ~binary).sum())

    benign_call = p <= 0.1
    pathogenic_call = p >= 0.9
    deferred = ~(benign_call | pathogenic_call)
    hc_tp = int(((y == 1) & pathogenic_call).sum())
    hc_fp = int(((y == 0) & pathogenic_call).sum())
    hc_tn = int(((y == 0) & benign_call).sum())
    hc_fn = int(((y == 1) & benign_call).sum())
    called = benign_call | pathogenic_call
    cal_intercept, cal_slope = calibration_parameters(y, p)

    return {
        "auroc": float(roc_auc_score(y, p)),
        "auprc": float(average_precision_score(y, p)),
        "brier": float(brier_score_loss(y, p)),
        "calibration_intercept": cal_intercept,
        "calibration_slope": cal_slope,
        "binary_sensitivity": safe_div(tp, tp + fn),
        "binary_specificity": safe_div(tn, tn + fp),
        "binary_ppv": safe_div(tp, tp + fp),
        "binary_npv": safe_div(tn, tn + fn),
        "defer_sensitivity": safe_div(hc_tp, int((y == 1).sum())),
        "defer_specificity_no_escalation": safe_div(
            int(((y == 0) & ~pathogenic_call).sum()), int((y == 0).sum())
        ),
        "defer_ppv": safe_div(hc_tp, hc_tp + hc_fp),
        "defer_npv": safe_div(hc_tn, hc_tn + hc_fn),
        "deferral_rate": float(deferred.mean()),
        "high_confidence_accuracy": safe_div(
            int((((y == 1) & pathogenic_call) | ((y == 0) & benign_call)).sum()), int(called.sum())
        ),
    }


def stratified_indices(y: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    parts = []
    for label in (0, 1):
        indices = np.flatnonzero(y == label)
        parts.append(rng.choice(indices, size=len(indices), replace=True))
    return np.concatenate(parts)


def gene_cluster_indices(frame: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    genes = frame["primary_gene"].dropna().unique()
    rows = {
        gene: np.flatnonzero(frame["primary_gene"].eq(gene).to_numpy())
        for gene in genes
    }
    selected = rng.choice(genes, size=len(genes), replace=True)
    return np.concatenate([rows[gene] for gene in selected])


def main() -> None:
    n_boot = 5000
    seed = 20260916
    rng = np.random.default_rng(seed)
    predictions = pd.read_csv(PREDICTIONS, sep="\t", low_memory=False)
    rows = []
    cluster_rows = []

    for subset, split_names in SUBSETS.items():
        frame = predictions[predictions["split_source_heldout"].isin(split_names)]
        y = frame["y_true"].astype(int).to_numpy()
        p = frame["pathogenic_probability"].astype(float).to_numpy()
        point = metrics(y, p)
        boot = {name: [] for name in point}
        for _ in range(n_boot):
            idx = stratified_indices(y, rng)
            sampled = metrics(y[idx], p[idx])
            for name, value in sampled.items():
                if np.isfinite(value):
                    boot[name].append(value)

        for name, estimate in point.items():
            values = np.asarray(boot[name], dtype=float)
            rows.append(
                {
                    "subset": subset,
                    "rows": len(frame),
                    "pathogenic": int(y.sum()),
                    "benign": int((y == 0).sum()),
                    "metric": name,
                    "estimate": estimate,
                    "ci_95_lower": float(np.quantile(values, 0.025)),
                    "ci_95_upper": float(np.quantile(values, 0.975)),
                    "successful_replicates": len(values),
                }
            )

        cluster_boot = {name: [] for name in point}
        attempts = 0
        while len(cluster_boot["auroc"]) < n_boot and attempts < n_boot * 3:
            attempts += 1
            idx = gene_cluster_indices(frame, rng)
            if np.unique(y[idx]).size != 2:
                continue
            sampled = metrics(y[idx], p[idx])
            for name, value in sampled.items():
                if np.isfinite(value):
                    cluster_boot[name].append(value)
        if len(cluster_boot["auroc"]) != n_boot:
            raise RuntimeError(f"Could not obtain {n_boot} gene-cluster replicates for {subset}")
        for name, estimate in point.items():
            values = np.asarray(cluster_boot[name], dtype=float)
            cluster_rows.append(
                {
                    "subset": subset,
                    "rows": len(frame),
                    "genes": int(frame["primary_gene"].nunique()),
                    "pathogenic": int(y.sum()),
                    "benign": int((y == 0).sum()),
                    "metric": name,
                    "estimate": estimate,
                    "ci_95_lower": float(np.quantile(values, 0.025)),
                    "ci_95_upper": float(np.quantile(values, 0.975)),
                    "successful_replicates": len(values),
                }
            )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT, sep="\t", index=False)
    pd.DataFrame(cluster_rows).to_csv(CLUSTER_OUT, sep="\t", index=False)
    MANIFEST.write_text(
        json.dumps(
            {
                "predictions": str(PREDICTIONS.relative_to(ROOT)),
                "output": str(OUT.relative_to(ROOT)),
                "subsets": SUBSETS,
                "bootstrap_replicates": n_boot,
                "seed": seed,
                "resampling": "stratified within binary outcome class",
                "gene_cluster_output": str(CLUSTER_OUT.relative_to(ROOT)),
                "gene_cluster_resampling": "genes sampled with replacement; all rows retained per selected gene",
                "interval": "2.5th and 97.5th bootstrap percentiles",
                "thresholds": {"binary": 0.5, "benign_call_max": 0.1, "pathogenic_call_min": 0.9},
            },
            indent=2,
        )
        + "\n"
    )
    print(OUT)


if __name__ == "__main__":
    main()
