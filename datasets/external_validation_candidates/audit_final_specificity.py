#!/usr/bin/env python3
"""Audit the specificity change in the final model-held-out benchmark."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, confusion_matrix, log_loss, roc_auc_score


ROOT = Path(__file__).resolve().parents[2]
OLD_PREDICTIONS = ROOT / "results/models/cardioqueue_v1_source_heldout/primary_binary_predictions_split_source_heldout.tsv"
NEW_PREDICTIONS = ROOT / "results/models/cardioqueue_v1_expanded_external_2026_public300/primary_binary_predictions_split_source_heldout.tsv"
FINAL_PREDICTIONS = ROOT / "results/model_performance/final_model_heldout_2026/combined_unique_predictions.tsv"
SOURCE_OBSERVATIONS = ROOT / "results/model_performance/final_model_heldout_2026/source_observation_predictions.tsv"
EXPANDED_MANIFEST = ROOT / "datasets/external_validation_candidates/expanded_external_2026/external_unique_variant_manifest.tsv"
OUTPUT_DIR = ROOT / "results/model_performance/final_model_heldout_2026"


def metrics(y: np.ndarray, probability: np.ndarray, threshold: float = 0.5) -> dict[str, float | int]:
    prediction = (probability >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, prediction, labels=[0, 1]).ravel()
    return {
        "rows": int(len(y)),
        "pathogenic": int(y.sum()),
        "benign": int((y == 0).sum()),
        "auroc": float(roc_auc_score(y, probability)),
        "auprc": float(average_precision_score(y, probability)),
        "brier": float(brier_score_loss(y, probability)),
        "log_loss": float(log_loss(y, probability)),
        "threshold": float(threshold),
        "sensitivity": float(tp / (tp + fn)),
        "specificity": float(tn / (tn + fp)),
        "tp": int(tp),
        "fp": int(fp),
        "tn": int(tn),
        "fn": int(fn),
        "mean_probability": float(probability.mean()),
    }


def read_model_predictions(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    return frame[["variant_id", "split_source_heldout", "y_true", "pathogenic_probability"]]


def best_threshold_for_specificity(y: np.ndarray, probability: np.ndarray, target: float) -> dict[str, float | int]:
    candidates = np.unique(probability)
    eligible: list[dict[str, float | int]] = []
    for threshold in candidates:
        result = metrics(y, probability, float(threshold))
        if result["specificity"] >= target:
            eligible.append(result)
    if not eligible:
        raise RuntimeError(f"No threshold achieved specificity >= {target}")
    return max(eligible, key=lambda item: (item["sensitivity"], -item["threshold"]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    old = read_model_predictions(OLD_PREDICTIONS).rename(columns={"pathogenic_probability": "old_probability"})
    new = read_model_predictions(NEW_PREDICTIONS).rename(columns={"pathogenic_probability": "new_probability"})
    source = pd.read_csv(SOURCE_OBSERVATIONS, sep="\t")
    legacy = source[source["source"].isin(["legacy_hiro", "legacy_emerge", "legacy_cardioboost"])].copy()
    legacy = legacy.drop_duplicates(["source", "variant_id", "external_binary_label"])
    legacy = legacy.merge(old[["variant_id", "old_probability"]], on="variant_id", how="left", validate="many_to_one")
    legacy = legacy.merge(new[["variant_id", "new_probability"]], on="variant_id", how="left", validate="many_to_one")
    if legacy[["old_probability", "new_probability"]].isna().any().any():
        raise RuntimeError("A retained legacy assertion lacks an old or new model prediction")

    same_cohort_rows = []
    for source_name, group in [("legacy_all", legacy), *legacy.groupby("source", sort=True)]:
        y = group["external_binary_label"].to_numpy(dtype=int)
        for model_name, column in [("CardioQueue_v1", "old_probability"), ("expanded_model", "new_probability")]:
            same_cohort_rows.append({"source": source_name, "model": model_name, **metrics(y, group[column].to_numpy())})
    same_cohort = pd.DataFrame(same_cohort_rows)
    same_cohort.to_csv(args.output_dir / "specificity_same_legacy_cohort.tsv", sep="\t", index=False)

    benign = legacy[legacy["external_binary_label"] == 0].copy()
    benign["old_call"] = (benign["old_probability"] >= 0.5).astype(int)
    benign["new_call"] = (benign["new_probability"] >= 0.5).astype(int)
    benign["transition"] = np.select(
        [
            (benign["old_call"] == 0) & (benign["new_call"] == 1),
            (benign["old_call"] == 1) & (benign["new_call"] == 0),
        ],
        ["old_correct_to_new_false_positive", "old_false_positive_to_new_correct"],
        default="unchanged",
    )
    benign.to_csv(args.output_dir / "specificity_benign_score_transitions.tsv", sep="\t", index=False)

    final = pd.read_csv(FINAL_PREDICTIONS, sep="\t")
    expanded = pd.read_csv(EXPANDED_MANIFEST, sep="\t")
    new_lookup = new[["variant_id", "new_probability"]]
    expanded = expanded.merge(new_lookup, on="variant_id", how="left", validate="one_to_one")
    expanded["y_true"] = expanded["combined_label"].map({"Benign": 0, "Pathogenic": 1})
    final_ids = set(final["variant_id"])
    removed = expanded[~expanded["variant_id"].isin(final_ids)].copy()

    composition_rows = []
    for name, frame, column in [
        ("final_775", final, "pathogenic_probability"),
        ("excluded_from_final_eligibility", removed, "new_probability"),
        ("broad_1233", expanded, "new_probability"),
    ]:
        composition_rows.append({"cohort": name, **metrics(frame["y_true"].to_numpy(dtype=int), frame[column].to_numpy())})
    composition = pd.DataFrame(composition_rows)
    composition.to_csv(args.output_dir / "specificity_cohort_composition.tsv", sep="\t", index=False)

    validation = new[new["split_source_heldout"] == "validation"].copy()
    epsilon = 1e-8
    validation_probability = np.clip(validation["new_probability"].to_numpy(), epsilon, 1 - epsilon)
    validation_logit = np.log(validation_probability / (1 - validation_probability)).reshape(-1, 1)
    calibrator = LogisticRegression(C=1e6, solver="lbfgs", max_iter=10_000)
    calibrator.fit(validation_logit, validation["y_true"].to_numpy(dtype=int))
    final_probability = np.clip(final["pathogenic_probability"].to_numpy(), epsilon, 1 - epsilon)
    final_logit = np.log(final_probability / (1 - final_probability)).reshape(-1, 1)
    calibrated_probability = calibrator.predict_proba(final_logit)[:, 1]
    recalibrated = final.copy()
    recalibrated["validation_logistic_probability"] = calibrated_probability
    recalibrated.to_csv(args.output_dir / "validation_only_recalibrated_predictions.tsv", sep="\t", index=False)

    y_final = final["y_true"].to_numpy(dtype=int)
    raw_metrics = metrics(y_final, final_probability)
    calibrated_metrics = metrics(y_final, calibrated_probability)
    post_hoc_90 = best_threshold_for_specificity(y_final, final_probability, 0.90)
    audit = {
        "interpretation": "The specificity change reflects both an upward score/operating-point shift and a harder final benign case mix; ranking discrimination remained stable on identical legacy alleles.",
        "feature_schema_identity": {
            "old_feature_count": 300,
            "new_feature_count": 300,
            "schemas_identical": True,
            "original_rows_compared": 85677,
            "missing_rows": 0,
            "feature_value_differences": 0,
            "note": "Recorded from the full streamed identity audit; see progress log for command provenance.",
        },
        "benign_score_shift_same_legacy_assertions": {
            "rows": int(len(benign)),
            "old_mean": float(benign["old_probability"].mean()),
            "new_mean": float(benign["new_probability"].mean()),
            "old_median": float(benign["old_probability"].median()),
            "new_median": float(benign["new_probability"].median()),
            "old_correct_to_new_false_positive": int((benign["transition"] == "old_correct_to_new_false_positive").sum()),
            "old_false_positive_to_new_correct": int((benign["transition"] == "old_false_positive_to_new_correct").sum()),
        },
        "validation_only_logistic_recalibration": {
            "validation_rows": int(len(validation)),
            "intercept": float(calibrator.intercept_[0]),
            "slope": float(calibrator.coef_[0, 0]),
            "raw_final_775": raw_metrics,
            "recalibrated_final_775": calibrated_metrics,
            "status": "post_hoc sensitivity; not a replacement primary result",
        },
        "post_hoc_threshold_for_at_least_90_percent_specificity": {
            **post_hoc_90,
            "status": "descriptive only; selected on the test set and unsuitable as a validated clinical cutoff",
        },
    }
    with (args.output_dir / "specificity_audit.json").open("w") as handle:
        json.dump(audit, handle, indent=2)

    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
