#!/usr/bin/env python3
"""Compare the original, public-schema, and weighting-sensitivity CardioQueue fits."""

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, log_loss, roc_auc_score


ROOT = Path(__file__).resolve().parents[2]
MODELS = {
    "Original schema, source and class weighted (338 features)": ROOT
    / "results/models/primary_binary_catboost_v0_strict_source_heldout_hiro_emerge_rescued/primary_binary_predictions_split_source_heldout.tsv",
    "CardioQueue v1.0, source and class weighted (300 features)": ROOT
    / "results/models/cardioqueue_v1_source_heldout/primary_binary_predictions_split_source_heldout.tsv",
    "Source weighted only (300 features)": ROOT
    / "results/models/cardioqueue_v1_source_weight_only/primary_binary_predictions_split_source_heldout.tsv",
}
OUT = ROOT / "results/model_performance/deployment_sensitivities.tsv"
PAPER = ROOT / "paper_MS_SI/tables/Table_S18_deployment_sensitivities.tsv"


def calibration(y: np.ndarray, score: np.ndarray) -> tuple[float, float]:
    clipped = np.clip(score, 1e-6, 1 - 1e-6)
    logit = np.log(clipped / (1 - clipped)).reshape(-1, 1)
    fit = LogisticRegression(C=1e6, solver="lbfgs").fit(logit, y)
    return float(fit.intercept_[0]), float(fit.coef_[0, 0])


def main() -> None:
    rows = []
    for name, path in MODELS.items():
        data = pd.read_csv(path, sep="\t")
        data = data[data["split_source_heldout"].astype(str).str.startswith("external_")]
        y = data["y_true"].astype(int).to_numpy()
        score = data["pathogenic_probability"].astype(float).to_numpy()
        intercept, slope = calibration(y, score)
        rows.append(
            {
                "model": name,
                "rows": len(data),
                "pathogenic": int(y.sum()),
                "mean_predicted_probability": float(score.mean()),
                "observed_prevalence": float(y.mean()),
                "auroc": roc_auc_score(y, score),
                "auprc": average_precision_score(y, score),
                "brier": brier_score_loss(y, score),
                "log_loss": log_loss(y, score),
                "calibration_intercept": intercept,
                "calibration_slope": slope,
            }
        )
    output = pd.DataFrame(rows)
    output.to_csv(OUT, sep="\t", index=False)
    output.to_csv(PAPER, sep="\t", index=False)
    print(output.to_string(index=False))


if __name__ == "__main__":
    main()
