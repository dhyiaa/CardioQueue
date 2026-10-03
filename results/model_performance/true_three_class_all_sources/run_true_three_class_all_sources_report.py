#!/usr/bin/env python3
"""Score source datasets with the true 3-class CatBoost model."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from catboost import CatBoostClassifier, Pool
from sklearn.metrics import cohen_kappa_score, confusion_matrix, f1_score, matthews_corrcoef


ROOT = Path(__file__).resolve().parents[3]
TABLE = ROOT / "datasets/modeling/ready/modeling_table_three_class_splits_weights_hiro_emerge_rescued.tsv"
MODEL_DIR = ROOT / "results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued"
MODEL_PATH = MODEL_DIR / "three_class_catboost_model.cbm"
FEATURE_SPEC = MODEL_DIR / "feature_columns.json"
OUT_DIR = ROOT / "results/model_performance/true_three_class_all_sources"

SOURCE_RECORD_TABLES = {
    "HiRO source records": ROOT / "datasets/hiro/full_dataset/interim/hiro_linked_source_records.tsv",
    "eMERGE source records": ROOT / "datasets/emerge/full_arrhythmia_gene_dataset/interim/emerge_linked_source_records_grch38_rescued.tsv",
    "CardioBoost source records": ROOT
    / "datasets/cardioboost/public_dataset/processed_for_cardiogenetics/model_inputs/cardioboost_linked_source_records.tsv",
}

LABELS = ["Benign", "VUS", "Pathogenic"]


def prepare_features(df: pd.DataFrame, features: list[str], cat_cols: list[str], numeric_cols: list[str]) -> pd.DataFrame:
    x = df.reindex(columns=features).copy()
    for col in cat_cols:
        if col in x.columns:
            x[col] = x[col].astype("string").fillna("__MISSING__").astype(str)
    for col in numeric_cols:
        if col in x.columns:
            x[col] = pd.to_numeric(x[col], errors="coerce")
    return x[features]


def summarize(name: str, df: pd.DataFrame, true_col: str = "model_label_3class") -> tuple[dict, pd.DataFrame]:
    out = df[df["predicted_label"].notna()].copy()
    labeled = out[out[true_col].isin(LABELS)].copy()
    cm = pd.DataFrame(
        confusion_matrix(labeled[true_col], labeled["predicted_label"], labels=LABELS),
        index=[f"true_{x}" for x in LABELS],
        columns=[f"pred_{x}" for x in LABELS],
    )
    y_true_path = (labeled[true_col] == "Pathogenic").astype(int)
    y_pred_path = (labeled["predicted_label"] == "Pathogenic").astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true_path, y_pred_path, labels=[0, 1]).ravel()

    true_counts = labeled[true_col].value_counts().to_dict()
    pred_counts = out["predicted_label"].value_counts().to_dict()
    total = len(labeled)
    summary = {
        "dataset": name,
        "rows_total": int(len(df)),
        "rows_with_prediction": int(len(out)),
        "rows_without_prediction": int(len(df) - len(out)),
        "labeled_rows_with_prediction": int(total),
        "true_benign": int(true_counts.get("Benign", 0)),
        "true_vus": int(true_counts.get("VUS", 0)),
        "true_pathogenic": int(true_counts.get("Pathogenic", 0)),
        "pred_benign": int(pred_counts.get("Benign", 0)),
        "pred_vus": int(pred_counts.get("VUS", 0)),
        "pred_pathogenic": int(pred_counts.get("Pathogenic", 0)),
        "exact_3class_accuracy": float((labeled[true_col] == labeled["predicted_label"]).mean()) if total else np.nan,
        "weighted_kappa": float(cohen_kappa_score(labeled[true_col], labeled["predicted_label"], labels=LABELS, weights="quadratic"))
        if total
        else np.nan,
        "macro_f1": float(f1_score(labeled[true_col], labeled["predicted_label"], labels=LABELS, average="macro", zero_division=0))
        if total
        else np.nan,
        "plp_sensitivity": float(tp / (tp + fn)) if (tp + fn) else np.nan,
        "specificity_non_plp": float(tn / (tn + fp)) if (tn + fp) else np.nan,
        "ppv_plp": float(tp / (tp + fp)) if (tp + fp) else np.nan,
        "npv_non_plp": float(tn / (tn + fn)) if (tn + fn) else np.nan,
        "binary_mcc_plp_vs_rest": float(matthews_corrcoef(y_true_path, y_pred_path)) if total else np.nan,
        "vus_prediction_rate": float((out["predicted_label"] == "VUS").mean()) if len(out) else np.nan,
        "plp_positive_calls": int((out["predicted_label"] == "Pathogenic").sum()),
        "true_positives": int(tp),
        "false_positives": int(fp),
        "plp_to_vus_deferrals": int(((labeled[true_col] == "Pathogenic") & (labeled["predicted_label"] == "VUS")).sum()),
        "plp_to_benign_errors": int(((labeled[true_col] == "Pathogenic") & (labeled["predicted_label"] == "Benign")).sum()),
        "benign_to_plp_escalations": int(((labeled[true_col] == "Benign") & (labeled["predicted_label"] == "Pathogenic")).sum()),
        "vus_to_plp_escalations": int(((labeled[true_col] == "VUS") & (labeled["predicted_label"] == "Pathogenic")).sum()),
    }
    return summary, cm


def fmt(value: object) -> str:
    if isinstance(value, float):
        if np.isnan(value):
            return "NA"
        return f"{value:.4f}"
    return str(value)


def write_markdown(summaries: list[dict], confusion_paths: dict[str, Path]) -> None:
    metrics = [
        "rows_total",
        "rows_with_prediction",
        "rows_without_prediction",
        "true_benign",
        "true_vus",
        "true_pathogenic",
        "pred_benign",
        "pred_vus",
        "pred_pathogenic",
        "exact_3class_accuracy",
        "weighted_kappa",
        "macro_f1",
        "plp_sensitivity",
        "specificity_non_plp",
        "ppv_plp",
        "npv_non_plp",
        "binary_mcc_plp_vs_rest",
        "vus_prediction_rate",
        "plp_positive_calls",
        "true_positives",
        "false_positives",
        "plp_to_vus_deferrals",
        "plp_to_benign_errors",
        "benign_to_plp_escalations",
        "vus_to_plp_escalations",
    ]
    lines = [
        "# True 3-Class CatBoost All-Sources Report",
        "",
        "Model: `results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued/three_class_catboost_model.cbm`",
        "",
        "Prediction classes are learned directly by the model:",
        "",
        "- `Benign`",
        "- `VUS`",
        "- `Pathogenic`",
        "",
        "Important: this is the exploratory true 3-class model, not the primary CardioBoost-comparable binary + deferral model.",
        "",
        "## Summary Metrics",
        "",
        "| Dataset | " + " | ".join(metrics) + " |",
        "|---|" + "|".join(["---:"] * len(metrics)) + "|",
    ]
    for row in summaries:
        lines.append("| " + row["dataset"] + " | " + " | ".join(fmt(row[m]) for m in metrics) + " |")

    lines += ["", "## Confusion Matrix Files", ""]
    for name, path in confusion_paths.items():
        lines.append(f"- {name}: `{path.relative_to(ROOT)}`")
    lines.append("")
    (OUT_DIR / "TRUE_THREE_CLASS_ALL_SOURCES_REPORT.md").write_text("\n".join(lines))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    spec = json.loads(FEATURE_SPEC.read_text())
    features = spec["features"]
    cat_cols = spec["categorical_features"]
    numeric_cols = spec["numeric_features"]

    df = pd.read_csv(TABLE, sep="\t", low_memory=False)
    score_df = df[df["training_slice_three_class_primary"].astype(str).str.lower().eq("true")].copy()
    x = prepare_features(score_df, features, cat_cols, numeric_cols)

    model = CatBoostClassifier()
    model.load_model(str(MODEL_PATH))
    proba = model.predict_proba(Pool(x, cat_features=cat_cols))
    pred_id = proba.argmax(axis=1)

    scored = score_df[
        [
            "variant_id",
            "primary_gene",
            "model_label_3class",
            "modeling_scope",
            "primary_model_inclusion",
            "in_clinvar",
            "in_hiro",
            "in_emerge",
            "in_cardioboost",
        ]
    ].copy()
    scored["predicted_label_id"] = pred_id
    scored["predicted_label"] = [LABELS[i] for i in pred_id]
    for i, label in enumerate(LABELS):
        scored[f"prob_{label}"] = proba[:, i]
    scored.to_csv(OUT_DIR / "all_variant_rows_true_three_class_predictions.tsv", sep="\t", index=False)

    summaries: list[dict] = []
    confusion_paths: dict[str, Path] = {}
    variant_sources = {
        "All variant rows": scored,
        "ClinVar variant rows": scored[scored["in_clinvar"].astype(bool)],
        "HiRO variant rows": scored[scored["in_hiro"].astype(bool)],
        "eMERGE variant rows": scored[scored["in_emerge"].astype(bool)],
        "CardioBoost variant rows": scored[scored["in_cardioboost"].astype(bool)],
    }
    for name, sub in variant_sources.items():
        stem = name.lower().replace(" ", "_")
        sub.to_csv(OUT_DIR / f"{stem}_predictions.tsv", sep="\t", index=False)
        summary, cm = summarize(name, sub)
        summaries.append(summary)
        cm_path = OUT_DIR / f"{stem}_confusion.tsv"
        cm.to_csv(cm_path, sep="\t")
        confusion_paths[name] = cm_path

    pred_map = scored[["variant_id", "predicted_label", "prob_Benign", "prob_VUS", "prob_Pathogenic"]].drop_duplicates("variant_id")
    for name, path in SOURCE_RECORD_TABLES.items():
        src = pd.read_csv(path, sep="\t", low_memory=False)
        merged = src.merge(pred_map, left_on="resolved_variant_id", right_on="variant_id", how="left")
        out_cols = [
            c
            for c in [
                "resolved_variant_id",
                "variant_id",
                "study",
                "source_row",
                "participant_id",
                "target_3class",
                "target_5class",
                "gene",
                "chrom",
                "pos",
                "ref",
                "alt",
                "predicted_label",
                "prob_Benign",
                "prob_VUS",
                "prob_Pathogenic",
            ]
            if c in merged.columns
        ]
        stem = name.lower().replace(" ", "_")
        merged[out_cols].to_csv(OUT_DIR / f"{stem}_predictions.tsv", sep="\t", index=False)
        summary, cm = summarize(name, merged, true_col="target_3class")
        summaries.append(summary)
        cm_path = OUT_DIR / f"{stem}_confusion.tsv"
        cm.to_csv(cm_path, sep="\t")
        confusion_paths[name] = cm_path

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(OUT_DIR / "true_three_class_summary_metrics.tsv", sep="\t", index=False)
    write_markdown(summaries, confusion_paths)
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
