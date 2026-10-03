#!/usr/bin/env python3
"""Train a small grid of 3-class CatBoost class-weight candidates."""

from __future__ import annotations

import json
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from catboost import CatBoostClassifier, Pool
from sklearn.metrics import accuracy_score, cohen_kappa_score, confusion_matrix, f1_score, matthews_corrcoef


ROOT = Path(__file__).resolve().parents[3]
TABLE = ROOT / "datasets/modeling/ready/modeling_table_three_class_splits_weights_hiro_emerge_rescued.tsv"
FEATURE_SPEC = ROOT / "results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued/feature_columns.json"
OUT_DIR = ROOT / "results/model_performance/three_class_weight_tuning"

SOURCE_RECORDS = {
    "HiRO source records": ROOT / "datasets/hiro/full_dataset/interim/hiro_linked_source_records.tsv",
    "eMERGE source records": ROOT / "datasets/emerge/full_arrhythmia_gene_dataset/interim/emerge_linked_source_records_grch38_rescued.tsv",
    "CardioBoost source records": ROOT
    / "datasets/cardioboost/public_dataset/processed_for_cardiogenetics/model_inputs/cardioboost_linked_source_records.tsv",
}

LABELS = ["Benign", "VUS", "Pathogenic"]
LABEL_TO_ID = {label: idx for idx, label in enumerate(LABELS)}

CANDIDATES = [
    ("baseline_like", {"Benign": 1.0, "VUS": 0.9, "Pathogenic": 2.5}),
    ("path_3_5", {"Benign": 1.0, "VUS": 0.9, "Pathogenic": 3.5}),
    ("path_4_vus_0_8", {"Benign": 1.0, "VUS": 0.8, "Pathogenic": 4.0}),
    ("path_5_vus_0_7", {"Benign": 1.0, "VUS": 0.7, "Pathogenic": 5.0}),
    ("path_6_vus_0_6", {"Benign": 1.0, "VUS": 0.6, "Pathogenic": 6.0}),
]


def prepare_features(df: pd.DataFrame, features: list[str], cat_cols: list[str], numeric_cols: list[str]) -> pd.DataFrame:
    x = df.reindex(columns=features).copy()
    for col in cat_cols:
        if col in x.columns:
            x[col] = x[col].astype("string").fillna("__MISSING__").astype(str)
    for col in numeric_cols:
        if col in x.columns:
            x[col] = pd.to_numeric(x[col], errors="coerce")
    return x[features]


def eval_predictions(y_true_labels: pd.Series, pred_labels: pd.Series) -> dict:
    d = pd.DataFrame({"true": y_true_labels, "pred": pred_labels})
    d = d[d["true"].isin(LABELS) & d["pred"].isin(LABELS)]
    yt = (d["true"] == "Pathogenic").astype(int)
    yp = (d["pred"] == "Pathogenic").astype(int)
    tn, fp, fn, tp = confusion_matrix(yt, yp, labels=[0, 1]).ravel()
    return {
        "rows": int(len(d)),
        "accuracy": float(accuracy_score(d["true"], d["pred"])),
        "macro_f1": float(f1_score(d["true"], d["pred"], labels=LABELS, average="macro", zero_division=0)),
        "weighted_kappa": float(cohen_kappa_score(d["true"], d["pred"], labels=LABELS, weights="quadratic")),
        "binary_mcc": float(matthews_corrcoef(yt, yp)),
        "plp_sensitivity": float(tp / (tp + fn)) if (tp + fn) else np.nan,
        "plp_specificity": float(tn / (tn + fp)) if (tn + fp) else np.nan,
        "plp_ppv": float(tp / (tp + fp)) if (tp + fp) else np.nan,
        "plp_npv": float(tn / (tn + fn)) if (tn + fn) else np.nan,
        "pred_benign": int((d["pred"] == "Benign").sum()),
        "pred_vus": int((d["pred"] == "VUS").sum()),
        "pred_pathogenic": int((d["pred"] == "Pathogenic").sum()),
        "pathogenic_to_vus": int(((d["true"] == "Pathogenic") & (d["pred"] == "VUS")).sum()),
        "pathogenic_to_benign": int(((d["true"] == "Pathogenic") & (d["pred"] == "Benign")).sum()),
        "benign_to_pathogenic": int(((d["true"] == "Benign") & (d["pred"] == "Pathogenic")).sum()),
        "vus_to_pathogenic": int(((d["true"] == "VUS") & (d["pred"] == "Pathogenic")).sum()),
    }


def threshold_labels(proba: np.ndarray, threshold: float) -> list[str]:
    out = []
    for b, v, p in proba:
        if p >= threshold:
            out.append("Pathogenic")
        elif b >= v:
            out.append("Benign")
        else:
            out.append("VUS")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", choices=[name for name, _ in CANDIDATES], help="Run one candidate only.")
    parser.add_argument("--iterations", type=int, default=450)
    parser.add_argument("--verbose", type=int, default=100)
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    spec = json.loads(FEATURE_SPEC.read_text())
    features = spec["features"]
    cat_cols = spec["categorical_features"]
    numeric_cols = spec["numeric_features"]

    df = pd.read_csv(TABLE, sep="\t", low_memory=False)
    df = df[df["training_slice_three_class_primary"].astype(str).str.lower().eq("true")].copy()
    x = prepare_features(df, features, cat_cols, numeric_cols)
    y = df["three_class_label"].astype(int).to_numpy()
    split = df["split_three_class_source_heldout"].astype(str)
    train_mask = split.eq("train").to_numpy()
    val_mask = split.eq("validation").to_numpy()

    source_conf = df["three_class_source_confidence_weight"].astype(float)
    true_label = df["model_label_3class"].astype(str)

    source_tables = {name: pd.read_csv(path, sep="\t", low_memory=False) for name, path in SOURCE_RECORDS.items()}

    rows = []
    candidates = CANDIDATES
    if args.candidate:
        candidates = [item for item in CANDIDATES if item[0] == args.candidate]

    for candidate_name, class_weights in candidates:
        print(f"Training candidate {candidate_name}: {class_weights}", flush=True)
        class_weight_series = true_label.map(class_weights).astype(float)
        weights = (source_conf * class_weight_series).to_numpy()

        model = CatBoostClassifier(
            iterations=args.iterations,
            learning_rate=0.035,
            depth=6,
            loss_function="MultiClass",
            eval_metric="TotalF1:average=Macro",
            random_seed=20260702,
            l2_leaf_reg=8,
            od_type="Iter",
            od_wait=80,
            allow_writing_files=False,
            verbose=args.verbose,
        )
        train_pool = Pool(x.loc[train_mask], y[train_mask], weight=weights[train_mask], cat_features=cat_cols)
        val_pool = Pool(x.loc[val_mask], y[val_mask], weight=weights[val_mask], cat_features=cat_cols)
        model.fit(train_pool, eval_set=val_pool, use_best_model=True)

        proba = model.predict_proba(Pool(x, cat_features=cat_cols))
        base_pred = [LABELS[i] for i in proba.argmax(axis=1)]
        pred_map_base = pd.DataFrame(
            {
                "variant_id": df["variant_id"].to_numpy(),
                "predicted_label": base_pred,
                "prob_Benign": proba[:, 0],
                "prob_VUS": proba[:, 1],
                "prob_Pathogenic": proba[:, 2],
            }
        ).drop_duplicates("variant_id")

        for policy_name, pred_labels in [("argmax", base_pred), ("path_threshold_0_20", threshold_labels(proba, 0.20))]:
            pred_series = pd.Series(pred_labels, index=df.index)
            for subset_name, mask in {
                "validation": val_mask,
                "external_hiro_variant": df["in_hiro"].astype(bool).to_numpy(),
                "external_emerge_variant": df["in_emerge"].astype(bool).to_numpy(),
                "external_cardioboost_variant": df["in_cardioboost"].astype(bool).to_numpy(),
            }.items():
                m = eval_predictions(df.loc[mask, "model_label_3class"], pred_series.loc[mask])
                rows.append(
                    {
                        "candidate": candidate_name,
                        "class_weights": json.dumps(class_weights, sort_keys=True),
                        "policy": policy_name,
                        "subset": subset_name,
                        **m,
                    }
                )

            tmp_map = pred_map_base.copy()
            tmp_map["predicted_label"] = pred_labels
            tmp_map = tmp_map[["variant_id", "predicted_label"]].drop_duplicates("variant_id")
            for source_name, source_df in source_tables.items():
                merged = source_df.merge(tmp_map, left_on="resolved_variant_id", right_on="variant_id", how="left")
                merged = merged[merged["predicted_label"].notna()].copy()
                m = eval_predictions(merged["target_3class"], merged["predicted_label"])
                rows.append(
                    {
                        "candidate": candidate_name,
                        "class_weights": json.dumps(class_weights, sort_keys=True),
                        "policy": policy_name,
                        "subset": source_name,
                        **m,
                    }
                )

        model_path = OUT_DIR / f"{candidate_name}.cbm"
        model.save_model(str(model_path))
        print(f"Saved {model_path}", flush=True)

    out = pd.DataFrame(rows)
    if args.candidate:
        result_path = OUT_DIR / f"three_class_weight_tuning_results.{args.candidate}.tsv"
    else:
        result_path = OUT_DIR / "three_class_weight_tuning_results.tsv"
    out.to_csv(result_path, sep="\t", index=False)
    compact = out[
        out["subset"].isin(["validation", "HiRO source records", "eMERGE source records", "CardioBoost source records"])
    ].copy()
    print(
        compact[
            [
                "candidate",
                "policy",
                "subset",
                "rows",
                "accuracy",
                "macro_f1",
                "plp_sensitivity",
                "plp_specificity",
                "plp_ppv",
                "pred_pathogenic",
                "pathogenic_to_vus",
                "pathogenic_to_benign",
                "benign_to_pathogenic",
                "vus_to_pathogenic",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
