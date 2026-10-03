#!/usr/bin/env python3
"""Build threshold-adjusted reports for the true 3-class model."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, cohen_kappa_score, confusion_matrix, f1_score, matthews_corrcoef


ROOT = Path(__file__).resolve().parents[3]
BASE_DIR = ROOT / "results/model_performance/true_three_class_all_sources"
OUT_DIR = ROOT / "results/model_performance/three_class_threshold_tuning"
THRESHOLD = 0.20
LABELS = ["Benign", "VUS", "Pathogenic"]

INPUTS = {
    "All variant rows": BASE_DIR / "all_variant_rows_predictions.tsv",
    "ClinVar variant rows": BASE_DIR / "clinvar_variant_rows_predictions.tsv",
    "HiRO variant rows": BASE_DIR / "hiro_variant_rows_predictions.tsv",
    "eMERGE variant rows": BASE_DIR / "emerge_variant_rows_predictions.tsv",
    "CardioBoost variant rows": BASE_DIR / "cardioboost_variant_rows_predictions.tsv",
    "HiRO source records": BASE_DIR / "hiro_source_records_predictions.tsv",
    "eMERGE source records": BASE_DIR / "emerge_source_records_predictions.tsv",
    "CardioBoost source records": BASE_DIR / "cardioboost_source_records_predictions.tsv",
}


def true_col(df: pd.DataFrame) -> str:
    return "model_label_3class" if "model_label_3class" in df.columns else "target_3class"


def threshold_call(df: pd.DataFrame, threshold: float = THRESHOLD) -> pd.Series:
    out = []
    for b, v, p in zip(df["prob_Benign"], df["prob_VUS"], df["prob_Pathogenic"]):
        if pd.isna(p):
            out.append(np.nan)
        elif p >= threshold:
            out.append("Pathogenic")
        elif b >= v:
            out.append("Benign")
        else:
            out.append("VUS")
    return pd.Series(out, index=df.index)


def summarize(name: str, df: pd.DataFrame, pred_col: str = "threshold_adjusted_label") -> tuple[dict, pd.DataFrame]:
    tcol = true_col(df)
    d = df[df[tcol].isin(LABELS) & df[pred_col].isin(LABELS)].copy()
    cm = pd.DataFrame(
        confusion_matrix(d[tcol], d[pred_col], labels=LABELS),
        index=[f"true_{x}" for x in LABELS],
        columns=[f"pred_{x}" for x in LABELS],
    )
    yt = (d[tcol] == "Pathogenic").astype(int)
    yp = (d[pred_col] == "Pathogenic").astype(int)
    tn, fp, fn, tp = confusion_matrix(yt, yp, labels=[0, 1]).ravel()
    pred_counts = d[pred_col].value_counts().to_dict()
    true_counts = d[tcol].value_counts().to_dict()
    return (
        {
            "dataset": name,
            "rows": int(len(d)),
            "true_benign": int(true_counts.get("Benign", 0)),
            "true_vus": int(true_counts.get("VUS", 0)),
            "true_pathogenic": int(true_counts.get("Pathogenic", 0)),
            "pred_benign": int(pred_counts.get("Benign", 0)),
            "pred_vus": int(pred_counts.get("VUS", 0)),
            "pred_pathogenic": int(pred_counts.get("Pathogenic", 0)),
            "accuracy": float(accuracy_score(d[tcol], d[pred_col])),
            "macro_f1": float(f1_score(d[tcol], d[pred_col], labels=LABELS, average="macro", zero_division=0)),
            "weighted_kappa": float(cohen_kappa_score(d[tcol], d[pred_col], labels=LABELS, weights="quadratic")),
            "plp_sensitivity": float(tp / (tp + fn)) if (tp + fn) else np.nan,
            "plp_specificity": float(tn / (tn + fp)) if (tn + fp) else np.nan,
            "plp_ppv": float(tp / (tp + fp)) if (tp + fp) else np.nan,
            "plp_npv": float(tn / (tn + fn)) if (tn + fn) else np.nan,
            "binary_mcc": float(matthews_corrcoef(yt, yp)),
            "plp_to_vus": int(((d[tcol] == "Pathogenic") & (d[pred_col] == "VUS")).sum()),
            "plp_to_benign": int(((d[tcol] == "Pathogenic") & (d[pred_col] == "Benign")).sum()),
            "benign_to_plp": int(((d[tcol] == "Benign") & (d[pred_col] == "Pathogenic")).sum()),
            "vus_to_plp": int(((d[tcol] == "VUS") & (d[pred_col] == "Pathogenic")).sum()),
        },
        cm,
    )


def fmt(x: object) -> str:
    if isinstance(x, float):
        if np.isnan(x):
            return "NA"
        return f"{x:.4f}"
    return str(x)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    confusion_paths = {}
    for name, path in INPUTS.items():
        df = pd.read_csv(path, sep="\t", low_memory=False)
        df["threshold_adjusted_label"] = threshold_call(df)
        stem = name.lower().replace(" ", "_")
        df.to_csv(OUT_DIR / f"{stem}_threshold_{str(THRESHOLD).replace('.', '_')}_predictions.tsv", sep="\t", index=False)
        summary, cm = summarize(name, df)
        rows.append(summary)
        cm_path = OUT_DIR / f"{stem}_threshold_{str(THRESHOLD).replace('.', '_')}_confusion.tsv"
        cm.to_csv(cm_path, sep="\t")
        confusion_paths[name] = cm_path

    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(OUT_DIR / f"threshold_{str(THRESHOLD).replace('.', '_')}_summary_metrics.tsv", sep="\t", index=False)

    cols = [
        "rows",
        "true_benign",
        "true_vus",
        "true_pathogenic",
        "pred_benign",
        "pred_vus",
        "pred_pathogenic",
        "accuracy",
        "macro_f1",
        "weighted_kappa",
        "plp_sensitivity",
        "plp_specificity",
        "plp_ppv",
        "plp_npv",
        "binary_mcc",
        "plp_to_vus",
        "plp_to_benign",
        "benign_to_plp",
        "vus_to_plp",
    ]
    lines = [
        "# 3-Class CatBoost Threshold Adjustment Report",
        "",
        f"Adjustment: call `Pathogenic` whenever `prob_Pathogenic >= {THRESHOLD}`; otherwise choose `Benign` vs `VUS` by the larger of those two probabilities.",
        "",
        "This keeps the true 3-class model but shifts the P/LP operating point toward higher sensitivity.",
        "",
        "## Summary Metrics",
        "",
        "| Dataset | " + " | ".join(cols) + " |",
        "|---|" + "|".join(["---:"] * len(cols)) + "|",
    ]
    for row in rows:
        lines.append("| " + row["dataset"] + " | " + " | ".join(fmt(row[c]) for c in cols) + " |")
    lines += ["", "## Confusion Matrix Files", ""]
    for name, path in confusion_paths.items():
        lines.append(f"- {name}: `{path.relative_to(ROOT)}`")
    lines.append("")
    (OUT_DIR / "THREE_CLASS_THRESHOLD_0_20_REPORT.md").write_text("\n".join(lines))
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
