#!/usr/bin/env python3
"""Train weighted 3-class CatBoost model: Benign / VUS / Pathogenic."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from catboost import CatBoostClassifier, Pool
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    cohen_kappa_score,
    confusion_matrix,
    f1_score,
    log_loss,
    matthews_corrcoef,
    precision_recall_fscore_support,
    roc_auc_score,
)


ROOT = Path(__file__).resolve().parents[3]
TABLE = ROOT / "datasets/modeling/ready/modeling_table_three_class_splits_weights_hiro_emerge_rescued.tsv"
LABELS = ["Benign", "VUS", "Pathogenic"]
LABEL_TO_ID = {label: i for i, label in enumerate(LABELS)}
ID_TO_LABEL = {i: label for label, i in LABEL_TO_ID.items()}


LEAKAGE_PATTERNS = [
    r"label",
    r"raw_label",
    r"review_star",
    r"review_status",
    r"confidence_tier",
    r"clinical_significance",
    r"split_",
    r"training_slice",
    r"training_exclude",
    r"sample_weight",
    r"source_confidence_weight",
    r"class_weight",
    r"model_label",
    r"primary_binary_label",
    r"three_class_label",
    r"target_",
    r"source_leakage_group",
    r"variant_id",
    r"coordinate_key",
    r"modeling_row_id",
    r"clinvar_",
    r"hiro_.*gene",
    r"emerge_",
    r"cardioboost_",
    r"hiro_labels",
    r"hiro_raw_label",
    r"hiro_5class_label",
    r"hiro_label_set",
    r"hiro_pathogenic_record_count",
    r"hiro_benign_record_count",
    r"hiro_vus_record_count",
    r"hiro_missing_label_record_count",
    r"hiro_nonmissing_label_class_count",
    r"hiro_label_discordant_flag",
    r"emerge_labels",
    r"cardioboost_labels",
    r"_caid$",
    r"_rsids$",
    r"existing_variation",
    r"selected_source_dataset",
    r"selected_variant_uid",
    r"source_dataset",
    r"source_record",
    r"source_row",
    r"source_uid",
    r"expert_panels",
    r"clingen_gene_validity_diseases",
]

DROP_EXACT = {
    "chrom",
    "pos",
    "ref",
    "alt",
    "genes",
    "sources",
    "in_clinvar",
    "in_hiro",
    "in_emerge",
    "in_cardioboost",
    "private_or_patient_linked",
    "phenotype_available",
    "source_count",
    "source_row_count",
}


def should_drop(col: str) -> bool:
    if col in DROP_EXACT:
        return True
    return any(re.search(pattern, col, flags=re.IGNORECASE) for pattern in LEAKAGE_PATTERNS)


def prepare_features(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str], list[str], list[str]]:
    feature_cols = [c for c in df.columns if not should_drop(c)]
    feature_df = df[feature_cols].copy()
    cat_cols = []
    numeric_cols = []
    for col in feature_df.columns:
        if feature_df[col].dtype == "object" or str(feature_df[col].dtype).startswith("string") or feature_df[col].dtype == bool:
            cat_cols.append(col)
            feature_df[col] = feature_df[col].astype("string").fillna("__MISSING__").astype(str)
        else:
            numeric_cols.append(col)
            feature_df[col] = pd.to_numeric(feature_df[col], errors="coerce")
    return feature_df, feature_cols, cat_cols, numeric_cols


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def confusion_df(y_true: np.ndarray, y_pred: np.ndarray) -> pd.DataFrame:
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1, 2])
    return pd.DataFrame(
        cm,
        index=[f"true_{x}" for x in LABELS],
        columns=[f"pred_{x}" for x in LABELS],
    )


def evaluate(y_true: np.ndarray, proba: np.ndarray) -> dict:
    pred = proba.argmax(axis=1)
    p, r, f1, support = precision_recall_fscore_support(y_true, pred, labels=[0, 1, 2], zero_division=0)
    out = {
        "rows": int(len(y_true)),
        "label_counts": {LABELS[i]: int((y_true == i).sum()) for i in range(3)},
        "accuracy": float(accuracy_score(y_true, pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, pred)),
        "macro_f1": float(f1_score(y_true, pred, labels=[0, 1, 2], average="macro", zero_division=0)),
        "weighted_f1": float(f1_score(y_true, pred, labels=[0, 1, 2], average="weighted", zero_division=0)),
        "weighted_kappa": float(cohen_kappa_score(y_true, pred, weights="quadratic")),
        "multiclass_mcc": float(matthews_corrcoef(y_true, pred)),
        "log_loss": float(log_loss(y_true, proba, labels=[0, 1, 2])),
        "per_class": {},
    }
    for i, label in enumerate(LABELS):
        out["per_class"][label] = {
            "precision": float(p[i]),
            "recall": float(r[i]),
            "f1": float(f1[i]),
            "support": int(support[i]),
        }
    try:
        out["ovr_macro_auroc"] = float(roc_auc_score(y_true, proba, multi_class="ovr", average="macro", labels=[0, 1, 2]))
        out["ovr_weighted_auroc"] = float(
            roc_auc_score(y_true, proba, multi_class="ovr", average="weighted", labels=[0, 1, 2])
        )
    except ValueError:
        out["ovr_macro_auroc"] = None
        out["ovr_weighted_auroc"] = None
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", default=str(TABLE))
    parser.add_argument("--split-col", default="split_three_class_source_heldout")
    parser.add_argument("--out-dir", default="results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued")
    parser.add_argument("--iterations", type=int, default=1800)
    args = parser.parse_args()

    table = Path(args.table)
    out = ROOT / args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    split_col = args.split_col

    df = pd.read_csv(table, sep="\t", low_memory=False)
    df = df[df["training_slice_three_class_primary"].astype(str).str.lower().eq("true")].copy()
    if split_col not in df.columns:
        raise ValueError(f"Missing split column: {split_col}")

    X, feature_cols, cat_cols, numeric_cols = prepare_features(df)
    y = df["three_class_label"].astype(int).to_numpy()
    weights = df["three_class_sample_weight"].astype(float).to_numpy()
    split = df[split_col].astype(str)
    train_mask = split.eq("train").to_numpy()
    val_mask = split.eq("validation").to_numpy()
    if not train_mask.any() or not val_mask.any():
        raise ValueError(f"{split_col} must contain train and validation rows")

    train_pool = Pool(X.loc[train_mask], y[train_mask], weight=weights[train_mask], cat_features=cat_cols)
    val_pool = Pool(X.loc[val_mask], y[val_mask], weight=weights[val_mask], cat_features=cat_cols)

    model = CatBoostClassifier(
        iterations=args.iterations,
        learning_rate=0.035,
        depth=6,
        loss_function="MultiClass",
        eval_metric="TotalF1:average=Macro",
        random_seed=20260702,
        l2_leaf_reg=8,
        od_type="Iter",
        od_wait=100,
        allow_writing_files=False,
        verbose=100,
    )
    model.fit(train_pool, eval_set=val_pool, use_best_model=True)

    pred_rows = []
    metrics = {
        "input_table": rel(table),
        "output_dir": rel(out),
        "split_col": split_col,
        "labels": LABELS,
        "n_features": len(feature_cols),
        "n_categorical_features": len(cat_cols),
        "n_numeric_features": len(numeric_cols),
        "best_iteration": int(model.get_best_iteration() or 0),
        "splits": {},
        "notes": [
            "Exploratory 3-class model: Benign / VUS / Pathogenic.",
            "VUS is an evidence-status class, so interpret VUS performance cautiously.",
            "Feature selection excludes labels, review stars, source flags, and source-derived label counts.",
        ],
    }

    ordered_splits = ["train", "validation"] + [
        s for s in sorted(split.unique()) if s not in {"train", "validation", "not_eligible"}
    ]
    for split_name in ordered_splits:
        mask = split.eq(split_name).to_numpy()
        if not mask.any():
            continue
        pool = Pool(X.loc[mask], cat_features=cat_cols)
        proba = model.predict_proba(pool)
        y_true = y[mask]
        y_pred = proba.argmax(axis=1)
        metrics["splits"][split_name] = evaluate(y_true, proba)
        cm = confusion_df(y_true, y_pred)
        cm.to_csv(out / f"confusion_matrix_{split_col}_{split_name}.tsv", sep="\t")
        tmp = df.loc[mask, ["variant_id", "primary_gene", "model_label_3class", split_col, "three_class_sample_weight"]].copy()
        tmp["y_true"] = y_true
        tmp["predicted_label_id"] = y_pred
        tmp["predicted_label"] = [ID_TO_LABEL[i] for i in y_pred]
        for i, label in enumerate(LABELS):
            tmp[f"prob_{label}"] = proba[:, i]
        pred_rows.append(tmp)

    predictions = pd.concat(pred_rows, ignore_index=True)
    predictions.to_csv(out / f"three_class_predictions_{split_col}.tsv", sep="\t", index=False)

    importance = pd.DataFrame(
        {
            "feature": feature_cols,
            "importance": model.get_feature_importance(train_pool, type="FeatureImportance"),
        }
    ).sort_values("importance", ascending=False)
    importance.to_csv(out / "feature_importance.tsv", sep="\t", index=False)

    (out / "feature_columns.json").write_text(
        json.dumps({"features": feature_cols, "categorical_features": cat_cols, "numeric_features": numeric_cols}, indent=2)
        + "\n"
    )
    (out / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    model.save_model(str(out / "three_class_catboost_model.cbm"))
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
