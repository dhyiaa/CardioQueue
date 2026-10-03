#!/usr/bin/env python3
"""Paired external evaluation after removing common effect-predictor fields."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


ROOT = Path(__file__).resolve().parents[3]
BASELINE = ROOT / "results/models/cardioqueue_v1_source_heldout"
ABLATION = ROOT / "results/models/cardioqueue_v1_no_precomputed_effect_predictors_source_heldout"
OUT = Path(__file__).resolve().parent
PAPER_TABLES = ROOT / "paper_MS_SI/tables"
N_BOOT = 5000
SEED = 20260924


def load(directory: Path, score_name: str) -> pd.DataFrame:
    path = directory / "primary_binary_predictions_split_source_heldout.tsv"
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    frame = frame[frame["split_source_heldout"].str.startswith("external_")].copy()
    return frame[["variant_id", "primary_gene", "split_source_heldout", "y_true", "pathogenic_probability"]].rename(
        columns={"pathogenic_probability": score_name}
    )


def metrics(y: np.ndarray, score: np.ndarray) -> dict[str, float | int]:
    called = score >= 0.5
    return {
        "rows": int(len(y)),
        "auroc": float(roc_auc_score(y, score)),
        "auprc": float(average_precision_score(y, score)),
        "brier": float(brier_score_loss(y, score)),
        "tp": int(((y == 1) & called).sum()),
        "fp": int(((y == 0) & called).sum()),
        "tn": int(((y == 0) & ~called).sum()),
        "fn": int(((y == 1) & ~called).sum()),
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    frame = load(BASELINE, "baseline_score").merge(
        load(ABLATION, "ablation_score"),
        on=["variant_id", "primary_gene", "split_source_heldout", "y_true"],
        validate="one_to_one",
    )
    y = frame["y_true"].astype(int).to_numpy()
    baseline = frame["baseline_score"].to_numpy()
    ablation = frame["ablation_score"].to_numpy()

    rows = []
    for name, score in [("CardioQueue_v1", baseline), ("No_precomputed_effect_predictors", ablation)]:
        rows.append({"model": name, **metrics(y, score)})
    metric_frame = pd.DataFrame(rows)
    metric_frame.to_csv(OUT / "ablation_external_metrics.tsv", sep="\t", index=False)
    metric_frame.to_csv(PAPER_TABLES / "Table_S10b_no_effect_predictor_ablation.tsv", sep="\t", index=False)

    rng = np.random.default_rng(SEED)
    groups = [np.flatnonzero(y == value) for value in (0, 1)]
    samples = {metric: [] for metric in ("auroc", "auprc", "brier")}
    for _ in range(N_BOOT):
        index = np.concatenate([rng.choice(group, len(group), replace=True) for group in groups])
        base_values = metrics(y[index], baseline[index])
        ablation_values = metrics(y[index], ablation[index])
        for metric in samples:
            samples[metric].append(ablation_values[metric] - base_values[metric])

    differences = []
    point_base = metrics(y, baseline)
    point_ablation = metrics(y, ablation)
    for metric, values in samples.items():
        differences.append(
            {
                "metric": metric,
                "ablation_minus_baseline": point_ablation[metric] - point_base[metric],
                "ci_95_low": float(np.quantile(values, 0.025)),
                "ci_95_high": float(np.quantile(values, 0.975)),
                "bootstrap_samples": N_BOOT,
            }
        )
    difference_frame = pd.DataFrame(differences)
    difference_frame.to_csv(OUT / "ablation_paired_bootstrap.tsv", sep="\t", index=False)
    difference_frame.to_csv(
        PAPER_TABLES / "Table_S10c_no_effect_predictor_paired_bootstrap.tsv", sep="\t", index=False
    )

    manifest = {
        "baseline_model": str(BASELINE.relative_to(ROOT)),
        "ablation_model": str(ABLATION.relative_to(ROOT)),
        "removed_feature_count": 58,
        "removed_families": [
            "SIFT", "PolyPhen-2", "REVEL", "MetaLR", "FATHMM-XF", "CADD", "AlphaMissense", "ESM-1b"
        ],
        "bootstrap_seed": SEED,
        "bootstrap_samples": N_BOOT,
    }
    (OUT / "ablation_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(metric_frame.to_string(index=False))
    print(difference_frame.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
