#!/usr/bin/env python3
"""Primary CardioQueue v1.0 metrics with outcome-stratified bootstrap intervals."""

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


ROOT = Path(__file__).resolve().parents[2]
PRED = ROOT / "results/models/cardioqueue_v1_source_heldout/primary_binary_predictions_split_source_heldout.tsv"
TABLE = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
OUT = ROOT / "results/model_performance/cardioqueue_v1_primary_metrics.tsv"
PAPER = ROOT / "paper_MS_SI/tables/Table_S10_cardioqueue_v1_primary_metrics.tsv"
N_BOOT = 5000
SEED = 20260919


def estimates(y: np.ndarray, score: np.ndarray) -> dict[str, float]:
    pred = score >= 0.5
    return {
        "auroc": roc_auc_score(y, score),
        "auprc": average_precision_score(y, score),
        "brier": brier_score_loss(y, score),
        "sensitivity_at_0_5": pred[y == 1].mean(),
        "specificity_at_0_5": (~pred[y == 0]).mean(),
    }


def summarize(name: str, data: pd.DataFrame, seed_offset: int) -> list[dict]:
    y = data["y_true"].astype(int).to_numpy()
    score = data["pathogenic_probability"].astype(float).to_numpy()
    point = estimates(y, score)
    groups = [np.flatnonzero(y == value) for value in (0, 1)]
    rng = np.random.default_rng(SEED + seed_offset)
    samples = {metric: [] for metric in point}
    for _ in range(N_BOOT):
        index = np.concatenate([rng.choice(group, len(group), replace=True) for group in groups])
        values = estimates(y[index], score[index])
        for metric, value in values.items():
            samples[metric].append(value)
    return [
        {
            "subset": name,
            "rows": len(data),
            "pathogenic": int(y.sum()),
            "benign": int((y == 0).sum()),
            "metric": metric,
            "estimate": value,
            "ci_95_low": float(np.quantile(samples[metric], 0.025)),
            "ci_95_high": float(np.quantile(samples[metric], 0.975)),
            "bootstrap_samples": N_BOOT,
        }
        for metric, value in point.items()
    ]


def main() -> None:
    pred = pd.read_csv(PRED, sep="\t")
    pred = pred[pred["split_source_heldout"].astype(str).str.startswith("external_")].copy()
    annotations = pd.read_csv(
        TABLE,
        sep="\t",
        usecols=["variant_id", "vep_any_missense", "gnomad_final_status"],
        low_memory=False,
    )
    pred = pred.merge(annotations, on="variant_id", validate="one_to_one")
    subsets = {
        "external_all": pred,
        "external_hiro": pred[pred["split_source_heldout"].eq("external_hiro")],
        "external_emerge": pred[pred["split_source_heldout"].eq("external_emerge")],
        "external_cardioboost": pred[pred["split_source_heldout"].eq("external_cardioboost")],
        "missense": pred[pred["vep_any_missense"].astype(str).str.lower().eq("true")],
        "gnomad_confirmed_absent": pred[pred["gnomad_final_status"].eq("confirmed_absent")],
    }
    rows = []
    for offset, (name, data) in enumerate(subsets.items()):
        if data["y_true"].nunique() == 2:
            rows.extend(summarize(name, data, offset))
    output = pd.DataFrame(rows)
    output.to_csv(OUT, sep="\t", index=False)
    output.to_csv(PAPER, sep="\t", index=False)
    print(output.to_string(index=False))


if __name__ == "__main__":
    main()
