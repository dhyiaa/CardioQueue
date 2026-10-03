#!/usr/bin/env python3
"""Evaluate the released CardioBoost models on their released holdout tables."""

from pathlib import Path

import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


ROOT = Path("results/model_performance/cardioboost_matched_benchmark/tables")
PANELS = {
    "Cardiomyopathy": ROOT / "official_cardioboost_cm_released_holdout_predictions.tsv",
    "Arrhythmia": ROOT / "official_cardioboost_arm_released_holdout_predictions.tsv",
}


def metrics(panel: str, data: pd.DataFrame) -> dict:
    y = pd.to_numeric(data["pathogenic"], errors="raise").astype(int)
    score = pd.to_numeric(data["pathogenicity"], errors="raise")
    pred = (score >= 0.5).astype(int)
    return {
        "panel": panel,
        "rows": len(data),
        "pathogenic": int(y.sum()),
        "benign": int((1 - y).sum()),
        "auroc": roc_auc_score(y, score),
        "auprc": average_precision_score(y, score),
        "brier": brier_score_loss(y, score),
        "sensitivity_at_0_5": float(pred[y.eq(1)].mean()),
        "specificity_at_0_5": float(1 - pred[y.eq(0)].mean()),
    }


def main() -> None:
    frames = {panel: pd.read_csv(path, sep="\t") for panel, path in PANELS.items()}
    rows = [metrics(panel, frame) for panel, frame in frames.items()]
    combined = pd.concat(frames.values(), ignore_index=True)
    rows.append(metrics("Combined", combined))
    output = pd.DataFrame(rows)
    target = ROOT / "cardioboost_released_holdout_positive_control.tsv"
    output.to_csv(target, sep="\t", index=False)
    print(output.to_string(index=False))
    print(f"\nWrote {target}")


if __name__ == "__main__":
    main()
