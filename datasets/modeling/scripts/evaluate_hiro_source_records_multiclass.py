#!/usr/bin/env python3
"""Map true 3-class CatBoost predictions back to HiRO source records."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from catboost import CatBoostClassifier, Pool


ROOT = Path(__file__).resolve().parents[3]
MATRIX = ROOT / "datasets/modeling/ready/modeling_table_three_class_splits_weights_hiro_emerge_rescued.tsv"
HIRO = ROOT / "datasets/hiro/full_dataset/interim/hiro_linked_source_records.tsv"
MODEL = ROOT / "results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued/three_class_catboost_model.cbm"
FEATURES = ROOT / "results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued/feature_columns.json"
OUT = ROOT / "results/model_performance/hiro_source_record_true_three_class"

LABELS = ["Benign", "VUS", "Pathogenic"]
PRED_LABELS_ALL = ["Benign", "VUS", "Pathogenic", "No_prediction"]


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def load_feature_spec(path: Path) -> tuple[list[str], list[str]]:
    data = json.loads(path.read_text())
    return data["features"], data["categorical_features"]


def prepare_feature_frame(df: pd.DataFrame, feature_cols: list[str], cat_cols: list[str]) -> pd.DataFrame:
    x = df[feature_cols].copy()
    cat_set = set(cat_cols)
    for col in feature_cols:
        if col in cat_set:
            x[col] = x[col].astype("string").fillna("__MISSING__").astype(str)
        else:
            x[col] = pd.to_numeric(x[col], errors="coerce")
    return x


def confusion(df: pd.DataFrame, pred_col: str, include_no_prediction: bool) -> pd.DataFrame:
    cols = PRED_LABELS_ALL if include_no_prediction else LABELS
    cm = pd.crosstab(df["target_3class"], df[pred_col], dropna=False)
    cm = cm.reindex(index=LABELS, columns=cols, fill_value=0)
    cm.index = [f"true_{x}" for x in cm.index]
    cm.columns = [f"pred_{x}" for x in cm.columns]
    return cm


def plot_cm(cm: pd.DataFrame, title: str, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    arr = cm.to_numpy()
    im = ax.imshow(arr, cmap="Blues")
    ax.set_xticks(range(len(cm.columns)), labels=[c.replace("pred_", "") for c in cm.columns], rotation=30, ha="right")
    ax.set_yticks(range(len(cm.index)), labels=[i.replace("true_", "") for i in cm.index])
    ax.set_xlabel("Predicted call")
    ax.set_ylabel("HiRO source-record label")
    ax.set_title(title)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            ax.text(j, i, str(arr[i, j]), ha="center", va="center", color="black")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)


def metrics(df: pd.DataFrame, pred_col: str) -> dict:
    called = df[df[pred_col] != "No_prediction"].copy()
    out = {
        "rows_total": int(len(df)),
        "rows_with_prediction": int(len(called)),
        "rows_no_prediction": int((df[pred_col] == "No_prediction").sum()),
        "true_label_counts_total": {str(k): int(v) for k, v in df["target_3class"].value_counts(dropna=False).to_dict().items()},
        "predicted_call_counts_total": {str(k): int(v) for k, v in df[pred_col].value_counts(dropna=False).to_dict().items()},
    }
    if len(called):
        out["called_exact_3class_accuracy"] = float((called["target_3class"] == called[pred_col]).mean())
        for label in LABELS:
            denom = int((called["target_3class"] == label).sum())
            out[f"recall_{label}"] = float(((called["target_3class"] == label) & (called[pred_col] == label)).sum() / denom) if denom else None
        out["pathogenic_to_benign_errors"] = int(((called["target_3class"] == "Pathogenic") & (called[pred_col] == "Benign")).sum())
        out["pathogenic_to_vus_errors"] = int(((called["target_3class"] == "Pathogenic") & (called[pred_col] == "VUS")).sum())
        out["benign_to_pathogenic_escalations"] = int(((called["target_3class"] == "Benign") & (called[pred_col] == "Pathogenic")).sum())
        out["vus_to_pathogenic_escalations"] = int(((called["target_3class"] == "VUS") & (called[pred_col] == "Pathogenic")).sum())
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", default=str(MATRIX))
    parser.add_argument("--hiro", default=str(HIRO))
    parser.add_argument("--model", default=str(MODEL))
    parser.add_argument("--features", default=str(FEATURES))
    parser.add_argument("--out-dir", default=str(OUT))
    args = parser.parse_args()

    matrix_path = Path(args.matrix)
    hiro_path = Path(args.hiro)
    model_path = Path(args.model)
    features_path = Path(args.features)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    feature_cols, cat_cols = load_feature_spec(features_path)
    matrix = pd.read_csv(matrix_path, sep="\t", low_memory=False)
    hiro = pd.read_csv(hiro_path, sep="\t", low_memory=False)
    hiro_matrix = matrix[matrix["in_hiro"].astype(str).str.lower().eq("true")].copy()
    x = prepare_feature_frame(hiro_matrix, feature_cols, cat_cols)

    model = CatBoostClassifier()
    model.load_model(str(model_path))
    proba = model.predict_proba(Pool(x, cat_features=cat_cols))
    pred_id = proba.argmax(axis=1)
    hiro_matrix["predicted_3class"] = [LABELS[i] for i in pred_id]
    for i, label in enumerate(LABELS):
        hiro_matrix[f"prob_{label}"] = proba[:, i]

    pred_cols = ["variant_id", "primary_gene", "model_label_3class", "predicted_3class"] + [f"prob_{x}" for x in LABELS]
    variant_predictions = hiro_matrix[pred_cols].drop_duplicates("variant_id")
    variant_predictions.to_csv(out_dir / "hiro_variant_level_true_three_class_predictions.tsv", sep="\t", index=False)

    scored = hiro.merge(
        variant_predictions.add_prefix("matrix_"),
        left_on="resolved_variant_id",
        right_on="matrix_variant_id",
        how="left",
    )
    scored["matrix_predicted_3class"] = scored["matrix_predicted_3class"].fillna("No_prediction")
    scored["prediction_join_status"] = np.where(scored["matrix_variant_id"].notna(), "predicted", "no_prediction_unresolved_or_not_in_matrix")
    scored.to_csv(out_dir / "hiro_source_record_true_three_class_predictions.tsv", sep="\t", index=False)

    cm_all = confusion(scored, "matrix_predicted_3class", include_no_prediction=True)
    cm_all.to_csv(out_dir / "hiro_source_record_true_three_class_confusion_all.tsv", sep="\t")
    plot_cm(cm_all, "HiRO source records: true 3-class model", out_dir / "hiro_source_record_true_three_class_confusion_all.png")
    pred_only = scored[scored["matrix_predicted_3class"] != "No_prediction"].copy()
    cm_pred = confusion(pred_only, "matrix_predicted_3class", include_no_prediction=False)
    cm_pred.to_csv(out_dir / "hiro_source_record_true_three_class_confusion_predicted_only.tsv", sep="\t")
    plot_cm(
        cm_pred,
        "HiRO source records with predictions: true 3-class model",
        out_dir / "hiro_source_record_true_three_class_confusion_predicted_only.png",
    )

    summary = {
        "matrix": rel(matrix_path),
        "hiro_source_records": rel(hiro_path),
        "model": rel(model_path),
        "feature_spec": rel(features_path),
        "out_dir": rel(out_dir),
        "variant_level_hiro_rows_scored": int(len(hiro_matrix)),
        "source_record_metrics": metrics(scored, "matrix_predicted_3class"),
    }
    (out_dir / "hiro_source_record_true_three_class_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
