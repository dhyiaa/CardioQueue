#!/usr/bin/env python3
"""Train the primary weighted binary CatBoost model."""

from __future__ import annotations

import json
import re
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from catboost import CatBoostClassifier, Pool
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    confusion_matrix,
    log_loss,
    precision_recall_curve,
    roc_auc_score,
)


ROOT = Path(__file__).resolve().parents[3]
TABLE = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights.tsv"


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


def prepare_features(
    df: pd.DataFrame,
    exclude_prefixes: tuple[str, ...] = (),
    exclude_patterns: tuple[str, ...] = (),
) -> tuple[pd.DataFrame, list[str], list[str], list[str]]:
    feature_cols = [
        c for c in df.columns
        if not should_drop(c)
        and not any(c.startswith(prefix) for prefix in exclude_prefixes)
        and not any(re.search(pattern, c, flags=re.IGNORECASE) for pattern in exclude_patterns)
    ]
    feature_df = df[feature_cols].copy()

    # CatBoost accepts NaN for numeric, but categorical/text values must be strings.
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


def metrics_at_threshold(y_true: np.ndarray, proba: np.ndarray, threshold: float) -> dict:
    pred = (proba >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, pred, labels=[0, 1]).ravel()
    sensitivity = tp / (tp + fn) if (tp + fn) else 0.0
    specificity = tn / (tn + fp) if (tn + fp) else 0.0
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    return {
        "threshold": float(threshold),
        "tp": int(tp),
        "fp": int(fp),
        "tn": int(tn),
        "fn": int(fn),
        "sensitivity_recall_pathogenic": float(sensitivity),
        "specificity_benign": float(specificity),
        "precision_ppv_pathogenic": float(precision),
    }


def choose_threshold_for_specificity(y_true: np.ndarray, proba: np.ndarray, target_specificity: float = 0.90) -> dict:
    precision, recall, thresholds = precision_recall_curve(y_true, proba)
    candidates = []
    for t in thresholds:
        m = metrics_at_threshold(y_true, proba, float(t))
        if m["specificity_benign"] >= target_specificity:
            candidates.append(m)
    if not candidates:
        return metrics_at_threshold(y_true, proba, 0.5)
    return max(candidates, key=lambda d: (d["sensitivity_recall_pathogenic"], d["precision_ppv_pathogenic"]))


def evaluate(y_true: np.ndarray, proba: np.ndarray) -> dict:
    out = {
        "rows": int(len(y_true)),
        "positives_pathogenic": int(y_true.sum()),
        "negatives_benign": int((y_true == 0).sum()),
        "auroc": float(roc_auc_score(y_true, proba)),
        "auprc": float(average_precision_score(y_true, proba)),
        "brier": float(brier_score_loss(y_true, proba)),
        "log_loss": float(log_loss(y_true, np.vstack([1 - proba, proba]).T, labels=[0, 1])),
        "threshold_0_5": metrics_at_threshold(y_true, proba, 0.5),
        "threshold_specificity_0_90": choose_threshold_for_specificity(y_true, proba, 0.90),
        "threshold_specificity_0_95": choose_threshold_for_specificity(y_true, proba, 0.95),
    }
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", default=str(TABLE))
    parser.add_argument("--split-col", default="split_internal_grouped")
    parser.add_argument("--slice-col", default="training_slice_primary_binary")
    parser.add_argument(
        "--augment-table",
        help="Optional TSV containing one row per variant_id and additional predictor columns.",
    )
    parser.add_argument("--out-dir", default="results/models/primary_binary_catboost_v0_strict")
    parser.add_argument(
        "--weight-col",
        default="sample_weight",
        help="Modeling-table column used as the CatBoost row weight.",
    )
    parser.add_argument(
        "--exclude-prefix",
        action="append",
        default=[],
        help="Exclude predictor columns beginning with this prefix. May be repeated.",
    )
    parser.add_argument(
        "--exclude-regex",
        action="append",
        default=[],
        help="Exclude predictor columns matching this regular expression. May be repeated.",
    )
    parser.add_argument("--thread-count", type=int, default=-1)
    args = parser.parse_args()
    out = ROOT / args.out_dir
    split_col = args.split_col
    slice_col = args.slice_col
    table = Path(args.table)

    out.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(table, sep="\t", low_memory=False)
    augmentation_features = []
    if args.augment_table:
        augment_path = Path(args.augment_table)
        augment = pd.read_csv(augment_path, sep="\t", low_memory=False)
        if "variant_id" not in augment.columns:
            raise ValueError("Augmentation table must contain variant_id")
        if not augment["variant_id"].is_unique:
            raise ValueError("Augmentation table must contain one row per variant_id")
        augmentation_features = [column for column in augment.columns if column != "variant_id"]
        collisions = sorted(set(augmentation_features) & set(df.columns))
        if collisions:
            raise ValueError(f"Augmentation columns already exist in modeling table: {collisions}")
        before_rows = len(df)
        df = df.merge(augment, on="variant_id", how="left", validate="one_to_one")
        if len(df) != before_rows:
            raise ValueError("Augmentation merge changed modeling row count")
    if slice_col not in df.columns:
        raise ValueError(f"Missing slice column: {slice_col}")
    df = df[df[slice_col].astype(str).str.lower().eq("true")].copy()
    if split_col not in df.columns:
        raise ValueError(f"Missing split column: {split_col}")

    X, feature_cols, cat_cols, numeric_cols = prepare_features(
        df, tuple(args.exclude_prefix), tuple(args.exclude_regex)
    )
    y = df["primary_binary_label"].astype(int).to_numpy()
    if args.weight_col not in df.columns:
        raise ValueError(f"Missing weight column: {args.weight_col}")
    weights = df[args.weight_col].astype(float).to_numpy()

    split = df[split_col].astype(str)
    train_mask = split.eq("train").to_numpy()
    val_mask = split.eq("validation").to_numpy()
    if not train_mask.any() or not val_mask.any():
        raise ValueError(f"{split_col} must contain train and validation rows")

    train_pool = Pool(X.loc[train_mask], y[train_mask], weight=weights[train_mask], cat_features=cat_cols)
    val_pool = Pool(X.loc[val_mask], y[val_mask], weight=weights[val_mask], cat_features=cat_cols)

    model = CatBoostClassifier(
        iterations=1500,
        learning_rate=0.035,
        depth=6,
        loss_function="Logloss",
        eval_metric="PRAUC",
        random_seed=20260702,
        l2_leaf_reg=6,
        od_type="Iter",
        od_wait=80,
        allow_writing_files=False,
        thread_count=args.thread_count,
        verbose=100,
    )
    model.fit(train_pool, eval_set=val_pool, use_best_model=True)

    pred_rows = []
    metrics = {
        "input_table": str(table.relative_to(ROOT)) if table.is_relative_to(ROOT) else str(table),
        "output_dir": str(out.relative_to(ROOT)),
        "split_col": split_col,
        "slice_col": slice_col,
        "weight_col": args.weight_col,
        "augmentation_table": args.augment_table,
        "augmentation_features": augmentation_features,
        "excluded_feature_prefixes": args.exclude_prefix,
        "excluded_feature_regexes": args.exclude_regex,
        "thread_count": args.thread_count,
        "n_features": len(feature_cols),
        "n_categorical_features": len(cat_cols),
        "n_numeric_features": len(numeric_cols),
        "best_iteration": int(model.get_best_iteration() or 0),
        "splits": {},
        "notes": [
            f"Binary model trained from rows selected by {slice_col} on {split_col}.",
            "Feature selection excludes labels, review stars, source flags, and source-derived label counts to reduce leakage.",
            f"Rows weighted with {args.weight_col}.",
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
        proba = model.predict_proba(pool)[:, 1]
        y_true = y[mask]
        metrics["splits"][split_name] = evaluate(y_true, proba)
        tmp = df.loc[mask, ["variant_id", "primary_gene", "model_label_3class", "sample_weight", split_col]].copy()
        tmp["y_true"] = y_true
        tmp["pathogenic_probability"] = proba
        pred_rows.append(tmp)

    predictions = pd.concat(pred_rows, ignore_index=True)
    predictions.to_csv(out / f"primary_binary_predictions_{split_col}.tsv", sep="\t", index=False)

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
    model.save_model(str(out / "primary_binary_catboost_model.cbm"))

    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
