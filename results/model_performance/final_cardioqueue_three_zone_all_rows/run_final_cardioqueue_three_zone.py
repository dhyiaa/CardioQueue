#!/usr/bin/env python3
"""Apply the final expanded-quarantine CardioQueue fit to every registry row."""

from pathlib import Path
import json

from catboost import CatBoostClassifier, Pool
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
TABLE = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
MODEL_DIR = ROOT / "results/models/cardioqueue_v1_expanded_external_2026_public300"
OUT = ROOT / "results/model_performance/final_cardioqueue_three_zone_all_rows"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    spec = json.loads((MODEL_DIR / "feature_columns.json").read_text())
    features = spec["features"]
    categorical = spec["categorical_features"]
    numeric = spec["numeric_features"]
    needed = set(features) | {
        "variant_id", "primary_gene", "model_label_3class", "modeling_scope", "primary_model_inclusion",
        "in_clinvar", "in_hiro", "in_emerge", "in_cardioboost"
    }
    frame = pd.read_csv(TABLE, sep="\t", usecols=lambda c: c in needed, low_memory=False)
    if frame["variant_id"].duplicated().any():
        raise ValueError("Registry matrix must contain one row per exact allele")

    x = frame.reindex(columns=features).copy()
    for col in categorical:
        x[col] = x[col].astype("string").fillna("__MISSING__").astype(str)
    for col in numeric:
        x[col] = pd.to_numeric(x[col], errors="coerce")

    model = CatBoostClassifier()
    model.load_model(str(MODEL_DIR / "primary_binary_catboost_model.cbm"))
    score = model.predict_proba(Pool(x, cat_features=categorical))[:, 1]
    if not np.isfinite(score).all() or ((score < 0) | (score > 1)).any():
        raise ValueError("Invalid CardioQueue scores")

    scored = frame[["variant_id", "primary_gene", "model_label_3class", "modeling_scope", "primary_model_inclusion",
                    "in_clinvar", "in_hiro", "in_emerge", "in_cardioboost"]].copy()
    scored["cardioqueue_score"] = score
    scored["review_zone"] = np.where(score <= 0.1, "lower_priority",
                                      np.where(score >= 0.9, "accelerated", "deferred"))
    scored.to_csv(OUT / "all_registry_scores.tsv", sep="\t", index=False)

    vus = scored.loc[scored["model_label_3class"].eq("VUS")].copy()
    vus.to_csv(OUT / "all_registry_vus_scores.tsv", sep="\t", index=False)
    rows = []
    for source, mask in [
        ("All registry", pd.Series(True, index=vus.index)),
        ("ClinVar", vus["in_clinvar"].astype(bool)),
        ("eMERGE", vus["in_emerge"].astype(bool)),
        ("HiRO", vus["in_hiro"].astype(bool)),
    ]:
        subset = vus.loc[mask]
        counts = subset["review_zone"].value_counts()
        rows.append({
            "source": source,
            "vus": len(subset),
            "lower_priority": int(counts.get("lower_priority", 0)),
            "deferred": int(counts.get("deferred", 0)),
            "accelerated": int(counts.get("accelerated", 0)),
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "vus_review_zone_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
