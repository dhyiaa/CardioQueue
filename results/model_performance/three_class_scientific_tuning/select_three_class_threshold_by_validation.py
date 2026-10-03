#!/usr/bin/env python3
"""Select 3-class operating thresholds using validation only, then test externally."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    cohen_kappa_score,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
)


ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = ROOT / "results/model_performance/three_class_scientific_tuning"
PRED_FILE = ROOT / "results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued/three_class_predictions_split_three_class_source_heldout.tsv"
ALL_SOURCES_DIR = ROOT / "results/model_performance/true_three_class_all_sources"

SOURCE_FILES = {
    "external_hiro_source_records": ALL_SOURCES_DIR / "hiro_source_records_predictions.tsv",
    "external_emerge_source_records": ALL_SOURCES_DIR / "emerge_source_records_predictions.tsv",
    "external_cardioboost_source_records": ALL_SOURCES_DIR / "cardioboost_source_records_predictions.tsv",
    "external_hiro_variant_rows": ALL_SOURCES_DIR / "hiro_variant_rows_predictions.tsv",
    "external_emerge_variant_rows": ALL_SOURCES_DIR / "emerge_variant_rows_predictions.tsv",
    "external_cardioboost_variant_rows": ALL_SOURCES_DIR / "cardioboost_variant_rows_predictions.tsv",
    "clinvar_variant_rows": ALL_SOURCES_DIR / "clinvar_variant_rows_predictions.tsv",
}

LABELS = ["Benign", "VUS", "Pathogenic"]


def make_call(df: pd.DataFrame, path_threshold: float) -> pd.Series:
    calls = []
    for b, v, p in zip(df["prob_Benign"], df["prob_VUS"], df["prob_Pathogenic"]):
        if pd.isna(b) or pd.isna(v) or pd.isna(p):
            calls.append(np.nan)
        elif p >= path_threshold:
            calls.append("Pathogenic")
        elif b >= v:
            calls.append("Benign")
        else:
            calls.append("VUS")
    return pd.Series(calls, index=df.index)


def eval_labels(true: pd.Series, pred: pd.Series) -> dict:
    d = pd.DataFrame({"true": true, "pred": pred})
    d = d[d["true"].isin(LABELS) & d["pred"].isin(LABELS)].copy()
    yt_path = (d["true"] == "Pathogenic").astype(int)
    yp_path = (d["pred"] == "Pathogenic").astype(int)
    tn, fp, fn, tp = confusion_matrix(yt_path, yp_path, labels=[0, 1]).ravel()
    return {
        "rows": int(len(d)),
        "accuracy": float(accuracy_score(d["true"], d["pred"])),
        "balanced_accuracy": float(balanced_accuracy_score(d["true"], d["pred"])),
        "macro_f1": float(f1_score(d["true"], d["pred"], labels=LABELS, average="macro", zero_division=0)),
        "weighted_f1": float(f1_score(d["true"], d["pred"], labels=LABELS, average="weighted", zero_division=0)),
        "weighted_kappa": float(cohen_kappa_score(d["true"], d["pred"], labels=LABELS, weights="quadratic")),
        "binary_mcc_plp": float(matthews_corrcoef(yt_path, yp_path)),
        "plp_sensitivity": float(tp / (tp + fn)) if (tp + fn) else np.nan,
        "plp_specificity": float(tn / (tn + fp)) if (tn + fp) else np.nan,
        "plp_ppv": float(tp / (tp + fp)) if (tp + fp) else np.nan,
        "plp_npv": float(tn / (tn + fn)) if (tn + fn) else np.nan,
        "pred_benign": int((d["pred"] == "Benign").sum()),
        "pred_vus": int((d["pred"] == "VUS").sum()),
        "pred_pathogenic": int((d["pred"] == "Pathogenic").sum()),
        "plp_to_vus": int(((d["true"] == "Pathogenic") & (d["pred"] == "VUS")).sum()),
        "plp_to_benign": int(((d["true"] == "Pathogenic") & (d["pred"] == "Benign")).sum()),
        "benign_to_plp": int(((d["true"] == "Benign") & (d["pred"] == "Pathogenic")).sum()),
        "vus_to_plp": int(((d["true"] == "VUS") & (d["pred"] == "Pathogenic")).sum()),
    }


def true_col(df: pd.DataFrame) -> str:
    return "model_label_3class" if "model_label_3class" in df.columns else "target_3class"


def pick_thresholds(validation_grid: pd.DataFrame) -> list[dict]:
    picks = []
    objectives = [
        ("max_validation_macro_f1", "macro_f1"),
        ("max_validation_weighted_kappa", "weighted_kappa"),
        ("max_validation_balanced_accuracy", "balanced_accuracy"),
        ("max_validation_plp_mcc", "binary_mcc_plp"),
    ]
    for name, col in objectives:
        row = validation_grid.sort_values([col, "macro_f1", "weighted_kappa"], ascending=False).iloc[0]
        picks.append({"selection_rule": name, "path_threshold": float(row["path_threshold"])})

    constrained = validation_grid[
        (validation_grid["plp_sensitivity"] >= 0.90)
        & (validation_grid["plp_specificity"] >= 0.90)
        & (validation_grid["plp_ppv"] >= 0.60)
    ].copy()
    if len(constrained):
        row = constrained.sort_values(["macro_f1", "weighted_kappa"], ascending=False).iloc[0]
        picks.append({"selection_rule": "validation_constrained_sens90_spec90_ppv60_max_macro_f1", "path_threshold": float(row["path_threshold"])})

    # Deduplicate rules that select the same threshold but keep all rule names.
    return picks


def fmt(x: object) -> str:
    if isinstance(x, float):
        if np.isnan(x):
            return "NA"
        return f"{x:.4f}"
    return str(x)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pred = pd.read_csv(PRED_FILE, sep="\t")
    val = pred[pred["split_three_class_source_heldout"].eq("validation")].copy()

    thresholds = np.round(np.arange(0.05, 0.951, 0.025), 3)
    val_rows = []
    for t in thresholds:
        calls = make_call(val, float(t))
        metrics = eval_labels(val["model_label_3class"], calls)
        val_rows.append({"path_threshold": float(t), **metrics})
    val_grid = pd.DataFrame(val_rows)
    val_grid.to_csv(OUT_DIR / "validation_threshold_grid.tsv", sep="\t", index=False)

    picks = pick_thresholds(val_grid)
    rows = []
    confusion_paths = {}
    for pick in picks:
        t = pick["path_threshold"]
        rule = pick["selection_rule"]
        val_calls = make_call(val, t)
        rows.append({"selection_rule": rule, "dataset": "validation_selection_set", "path_threshold": t, **eval_labels(val["model_label_3class"], val_calls)})
        for dataset, path in SOURCE_FILES.items():
            df = pd.read_csv(path, sep="\t", low_memory=False)
            df["validation_selected_call"] = make_call(df, t)
            tc = true_col(df)
            metrics = eval_labels(df[tc], df["validation_selected_call"])
            rows.append({"selection_rule": rule, "dataset": dataset, "path_threshold": t, **metrics})
            cm_df = df[df[tc].isin(LABELS) & df["validation_selected_call"].isin(LABELS)].copy()
            cm = pd.DataFrame(
                confusion_matrix(
                    cm_df[tc],
                    cm_df["validation_selected_call"],
                    labels=LABELS,
                ),
                index=[f"true_{x}" for x in LABELS],
                columns=[f"pred_{x}" for x in LABELS],
            )
            cm_path = OUT_DIR / f"{rule}_{dataset}_threshold_{str(t).replace('.', '_')}_confusion.tsv"
            cm.to_csv(cm_path, sep="\t")
            confusion_paths[(rule, dataset)] = cm_path

    result = pd.DataFrame(rows)
    result.to_csv(OUT_DIR / "validation_selected_threshold_external_results.tsv", sep="\t", index=False)

    cols = [
        "path_threshold",
        "rows",
        "accuracy",
        "macro_f1",
        "weighted_kappa",
        "plp_sensitivity",
        "plp_specificity",
        "plp_ppv",
        "plp_npv",
        "binary_mcc_plp",
        "pred_pathogenic",
        "plp_to_vus",
        "plp_to_benign",
        "benign_to_plp",
        "vus_to_plp",
    ]
    lines = [
        "# Scientifically Controlled 3-Class Threshold Selection",
        "",
        "Thresholds are selected using the internal validation split only. External HiRO, eMERGE, CardioBoost, and ClinVar summaries below are held-out evaluations of those validation-selected rules.",
        "",
        "Decision rule: call `Pathogenic` if `prob_Pathogenic >= threshold`; otherwise choose `Benign` vs `VUS` by the larger probability.",
        "",
    ]
    for rule in result["selection_rule"].drop_duplicates():
        lines += [f"## {rule}", ""]
        sub = result[result["selection_rule"].eq(rule)].copy()
        lines.append("| Dataset | " + " | ".join(cols) + " |")
        lines.append("|---|" + "|".join(["---:"] * len(cols)) + "|")
        for _, row in sub.iterrows():
            lines.append("| " + row["dataset"] + " | " + " | ".join(fmt(row[c]) for c in cols) + " |")
        lines.append("")
    (OUT_DIR / "VALIDATION_SELECTED_THRESHOLD_REPORT.md").write_text("\n".join(lines))

    compact = result[
        result["dataset"].isin(["validation_selection_set", "external_hiro_source_records", "external_emerge_source_records", "external_cardioboost_source_records"])
    ]
    print(compact[["selection_rule", "dataset", *cols]].to_string(index=False))


if __name__ == "__main__":
    main()
