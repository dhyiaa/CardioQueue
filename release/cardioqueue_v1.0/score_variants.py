#!/usr/bin/env python3
"""Score an annotation table with the frozen CardioQueue v1.0 model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from catboost import CatBoostClassifier, Pool


VERSION = "CardioQueue v1.0"
HERE = Path(__file__).resolve().parent
MODEL = HERE / "model/cardioqueue_v1.0.cbm"
SCHEMA = HERE / "model/feature_columns.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Tab-delimited annotation table")
    parser.add_argument("output", help="Output TSV")
    args = parser.parse_args()

    schema = json.loads(SCHEMA.read_text())
    features = schema["features"]
    categorical = schema["categorical_features"]
    numeric = schema["numeric_features"]
    data = pd.read_csv(args.input, sep="\t", low_memory=False)
    if "variant_id" not in data.columns:
        raise ValueError("Input must contain variant_id")
    missing = [column for column in features if column not in data.columns]
    if missing:
        preview = ", ".join(missing[:20])
        raise ValueError(f"Input is missing {len(missing)} required features: {preview}")

    matrix = data[features].copy()
    for column in categorical:
        matrix[column] = matrix[column].astype("string").fillna("__MISSING__").astype(str)
    for column in numeric:
        matrix[column] = pd.to_numeric(matrix[column], errors="coerce")

    model = CatBoostClassifier()
    model.load_model(str(MODEL))
    score = model.predict_proba(Pool(matrix, cat_features=categorical))[:, 1]
    output = pd.DataFrame(
        {
            "variant_id": data["variant_id"].astype(str),
            "cardioqueue_version": VERSION,
            "review_priority_score": score,
        }
    )
    output.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
