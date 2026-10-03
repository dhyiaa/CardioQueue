#!/usr/bin/env python3
"""Summarize primary binary model calls on binary-labeled rows only."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix, matthews_corrcoef


ROOT = Path(__file__).resolve().parents[3]
THREE_ZONE_DIR = ROOT / "results/model_performance/primary_binary_three_zone_all_sources"
OUT_DIR = ROOT / "results/model_performance/primary_binary_three_zone_binary_only"

INPUTS = {
    "All variant rows": THREE_ZONE_DIR / "all_variant_rows_predictions.tsv",
    "ClinVar variant rows": THREE_ZONE_DIR / "clinvar_variant_rows_predictions.tsv",
    "HiRO variant rows": THREE_ZONE_DIR / "hiro_variant_rows_predictions.tsv",
    "eMERGE variant rows": THREE_ZONE_DIR / "emerge_variant_rows_predictions.tsv",
    "CardioBoost variant rows": THREE_ZONE_DIR / "cardioboost_variant_rows_predictions.tsv",
    "HiRO source records": THREE_ZONE_DIR / "hiro_source_records_predictions.tsv",
    "eMERGE source records": THREE_ZONE_DIR / "emerge_source_records_predictions.tsv",
    "CardioBoost source records": THREE_ZONE_DIR / "cardioboost_source_records_predictions.tsv",
}

LABEL_MAP = {"Benign": 0, "Pathogenic": 1}
CALLS = ["Benign", "VUS", "Pathogenic"]


def label_column(df: pd.DataFrame) -> str:
    if "model_label_3class" in df.columns:
        return "model_label_3class"
    if "target_3class" in df.columns:
        return "target_3class"
    raise ValueError("No true label column found")


def summarize(name: str, path: Path) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    true_col = label_column(df)
    binary = df[df[true_col].isin(LABEL_MAP)].copy()
    binary = binary[binary["three_zone_call"].isin(CALLS)].copy()

    cm3 = pd.DataFrame(
        confusion_matrix(binary[true_col], binary["three_zone_call"], labels=["Benign", "Pathogenic", "VUS"]),
        index=["true_Benign", "true_Pathogenic", "true_VUS"],
        columns=["pred_Benign", "pred_Pathogenic", "pred_VUS"],
    ).loc[["true_Benign", "true_Pathogenic"], ["pred_Benign", "pred_VUS", "pred_Pathogenic"]]

    callable_binary = binary[binary["three_zone_call"].isin(["Benign", "Pathogenic"])].copy()
    y_true_all = binary[true_col].map(LABEL_MAP).astype(int)
    y_pred_path_all = (binary["three_zone_call"] == "Pathogenic").astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true_all, y_pred_path_all, labels=[0, 1]).ravel()

    if len(callable_binary):
        y_true_called = callable_binary[true_col].map(LABEL_MAP).astype(int)
        y_pred_called = callable_binary["three_zone_call"].map(LABEL_MAP).astype(int)
        ctn, cfp, cfn, ctp = confusion_matrix(y_true_called, y_pred_called, labels=[0, 1]).ravel()
        called_acc = float((y_true_called == y_pred_called).mean())
        called_sens = ctp / (ctp + cfn) if (ctp + cfn) else np.nan
        called_spec = ctn / (ctn + cfp) if (ctn + cfp) else np.nan
        called_ppv = ctp / (ctp + cfp) if (ctp + cfp) else np.nan
        called_npv = ctn / (ctn + cfn) if (ctn + cfn) else np.nan
        called_mcc = matthews_corrcoef(y_true_called, y_pred_called)
    else:
        called_acc = called_sens = called_spec = called_ppv = called_npv = called_mcc = np.nan

    summary = {
        "dataset": name,
        "binary_rows_with_prediction": int(len(binary)),
        "true_benign": int((binary[true_col] == "Benign").sum()),
        "true_pathogenic": int((binary[true_col] == "Pathogenic").sum()),
        "pred_benign": int((binary["three_zone_call"] == "Benign").sum()),
        "pred_vus_defer": int((binary["three_zone_call"] == "VUS").sum()),
        "pred_pathogenic": int((binary["three_zone_call"] == "Pathogenic").sum()),
        "deferral_rate": float((binary["three_zone_call"] == "VUS").mean()) if len(binary) else np.nan,
        "sensitivity_pathogenic_counting_defers_as_negative": float(tp / (tp + fn)) if (tp + fn) else np.nan,
        "specificity_counting_defers_as_non_pathogenic": float(tn / (tn + fp)) if (tn + fp) else np.nan,
        "ppv_pathogenic": float(tp / (tp + fp)) if (tp + fp) else np.nan,
        "npv_counting_defers_as_non_pathogenic": float(tn / (tn + fn)) if (tn + fn) else np.nan,
        "mcc_counting_defers_as_non_pathogenic": float(matthews_corrcoef(y_true_all, y_pred_path_all)) if len(binary) else np.nan,
        "called_rows_excluding_defers": int(len(callable_binary)),
        "called_accuracy_excluding_defers": called_acc,
        "called_sensitivity_excluding_defers": called_sens,
        "called_specificity_excluding_defers": called_spec,
        "called_ppv_excluding_defers": called_ppv,
        "called_npv_excluding_defers": called_npv,
        "called_mcc_excluding_defers": called_mcc,
        "true_positives": int(tp),
        "false_positives": int(fp),
        "true_negatives": int(tn),
        "false_negatives_including_defers": int(fn),
        "pathogenic_to_vus_deferrals": int(((binary[true_col] == "Pathogenic") & (binary["three_zone_call"] == "VUS")).sum()),
        "pathogenic_to_benign_errors": int(((binary[true_col] == "Pathogenic") & (binary["three_zone_call"] == "Benign")).sum()),
        "benign_to_pathogenic_escalations": int(((binary[true_col] == "Benign") & (binary["three_zone_call"] == "Pathogenic")).sum()),
    }
    return summary, cm3, binary


def fmt(x: object) -> str:
    if isinstance(x, float):
        if np.isnan(x):
            return "NA"
        return f"{x:.4f}"
    return str(x)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summaries = []
    confusion_paths = {}
    for name, path in INPUTS.items():
        summary, cm, binary = summarize(name, path)
        summaries.append(summary)
        stem = name.lower().replace(" ", "_")
        binary.to_csv(OUT_DIR / f"{stem}_binary_only_predictions.tsv", sep="\t", index=False)
        cm_path = OUT_DIR / f"{stem}_binary_only_confusion.tsv"
        cm.to_csv(cm_path, sep="\t")
        confusion_paths[name] = cm_path

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(OUT_DIR / "primary_binary_binary_only_summary_metrics.tsv", sep="\t", index=False)

    cols = [
        "binary_rows_with_prediction",
        "true_benign",
        "true_pathogenic",
        "pred_benign",
        "pred_vus_defer",
        "pred_pathogenic",
        "deferral_rate",
        "sensitivity_pathogenic_counting_defers_as_negative",
        "specificity_counting_defers_as_non_pathogenic",
        "ppv_pathogenic",
        "npv_counting_defers_as_non_pathogenic",
        "mcc_counting_defers_as_non_pathogenic",
        "called_rows_excluding_defers",
        "called_accuracy_excluding_defers",
        "called_sensitivity_excluding_defers",
        "called_specificity_excluding_defers",
        "called_ppv_excluding_defers",
        "called_npv_excluding_defers",
        "called_mcc_excluding_defers",
        "pathogenic_to_vus_deferrals",
        "pathogenic_to_benign_errors",
        "benign_to_pathogenic_escalations",
    ]
    lines = [
        "# Primary Binary CatBoost On Binary-Labeled Rows Only",
        "",
        "Truth labels are restricted to `Benign` and `Pathogenic`; true VUS rows are excluded from this report.",
        "",
        "Prediction rule remains:",
        "",
        "- `P(pathogenic) <= 0.1`: Benign-like",
        "- `0.1 < P(pathogenic) < 0.9`: VUS/defer",
        "- `P(pathogenic) >= 0.9`: Pathogenic-like",
        "",
        "Metrics are shown two ways: first counting deferrals as non-pathogenic for binary P/LP-vs-rest metrics, and second on only called rows after excluding deferrals.",
        "",
        "## Summary Metrics",
        "",
        "| Dataset | " + " | ".join(cols) + " |",
        "|---|" + "|".join(["---:"] * len(cols)) + "|",
    ]
    for row in summaries:
        lines.append("| " + row["dataset"] + " | " + " | ".join(fmt(row[c]) for c in cols) + " |")
    lines += ["", "## Confusion Matrix Files", ""]
    for name, path in confusion_paths.items():
        lines.append(f"- {name}: `{path.relative_to(ROOT)}`")
    lines.append("")
    (OUT_DIR / "PRIMARY_BINARY_BINARY_ONLY_REPORT.md").write_text("\n".join(lines))
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
