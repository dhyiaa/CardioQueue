#!/usr/bin/env python3
"""Build 3-class training flags, weights, and split assignments."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split


ROOT = Path(__file__).resolve().parents[3]
MATRIX = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
OUT_MATRIX = ROOT / "datasets/modeling/ready/modeling_table_three_class_splits_weights_hiro_emerge_rescued.tsv"
TRAIN_TABLE = ROOT / "datasets/modeling/ready/three_class_catboost_table_hiro_emerge_rescued.tsv"
SUMMARY = ROOT / "datasets/modeling/ready/three_class_slices_and_splits_hiro_emerge_rescued.summary.json"

LABEL_TO_ID = {"Benign": 0, "VUS": 1, "Pathogenic": 2}


def as_bool(s: pd.Series) -> pd.Series:
    return s.astype("string").str.lower().isin(["true", "1", "yes", "y", "t"])


def num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def make_strat_label(df: pd.DataFrame) -> pd.Series:
    gene_counts = df["primary_gene"].fillna("UNK").value_counts()
    gene_bin = df["primary_gene"].fillna("UNK").where(df["primary_gene"].fillna("UNK").map(gene_counts) >= 100, "OTHER_LOW_N")
    return df["model_label_3class"].astype(str) + "|" + gene_bin.astype(str)


def safe_strat_label(df: pd.DataFrame) -> pd.Series:
    y = make_strat_label(df)
    if y.value_counts().min() < 2:
        y = df["model_label_3class"].astype(str)
    return y


def assign_internal_grouped_split(df: pd.DataFrame, eligible: pd.Series, seed: int = 20260702) -> pd.Series:
    split = pd.Series("not_eligible", index=df.index, dtype="object")
    idx = df.index[eligible].to_numpy()
    y = safe_strat_label(df.loc[idx])
    train_idx, temp_idx = train_test_split(idx, train_size=0.70, random_state=seed, stratify=y)
    y_temp = safe_strat_label(df.loc[temp_idx])
    val_idx, test_idx = train_test_split(temp_idx, train_size=0.50, random_state=seed + 1, stratify=y_temp)
    split.loc[train_idx] = "train"
    split.loc[val_idx] = "validation"
    split.loc[test_idx] = "test"
    return split


def assign_source_heldout_split(df: pd.DataFrame, eligible: pd.Series, seed: int = 20260702) -> pd.Series:
    split = pd.Series("not_eligible", index=df.index, dtype="object")
    is_hiro = as_bool(df["in_hiro"])
    is_emerge = as_bool(df["in_emerge"])
    is_cardio = as_bool(df["in_cardioboost"])
    external = eligible & (is_hiro | is_emerge | is_cardio)

    split.loc[eligible & is_hiro] = "external_hiro"
    split.loc[eligible & ~is_hiro & is_emerge] = "external_emerge"
    split.loc[eligible & ~is_hiro & ~is_emerge & is_cardio] = "external_cardioboost"

    pool = eligible & ~external
    idx = df.index[pool].to_numpy()
    y = safe_strat_label(df.loc[idx])
    train_idx, val_idx = train_test_split(idx, train_size=0.85, random_state=seed, stratify=y)
    split.loc[train_idx] = "train"
    split.loc[val_idx] = "validation"
    return split


def assign_gene_stress_split(df: pd.DataFrame, eligible: pd.Series, seed: int = 20260702) -> pd.Series:
    split = pd.Series("not_eligible", index=df.index, dtype="object")
    counts = pd.crosstab(df.loc[eligible, "primary_gene"], df.loc[eligible, "model_label_3class"])
    for col in ["Benign", "VUS", "Pathogenic"]:
        if col not in counts:
            counts[col] = 0
    sparse_genes = set(counts.index[counts["Pathogenic"] < 30])
    sparse = eligible & df["primary_gene"].isin(sparse_genes)
    split.loc[sparse] = "sparse_gene_stress_test"

    pool = eligible & ~sparse
    idx = df.index[pool].to_numpy()
    y = safe_strat_label(df.loc[idx])
    train_idx, temp_idx = train_test_split(idx, train_size=0.70, random_state=seed, stratify=y)
    y_temp = safe_strat_label(df.loc[temp_idx])
    val_idx, test_idx = train_test_split(temp_idx, train_size=0.50, random_state=seed + 1, stratify=y_temp)
    split.loc[train_idx] = "train"
    split.loc[val_idx] = "validation"
    split.loc[test_idx] = "well_represented_test"
    return split


def source_confidence_weight(df: pd.DataFrame) -> pd.Series:
    w = pd.Series(1.0, index=df.index)
    stars = num(df["clinvar_review_stars"])
    is_clinvar = as_bool(df["in_clinvar"])
    w.loc[is_clinvar & stars.eq(1)] = 0.70
    w.loc[is_clinvar & stars.eq(2)] = 1.00
    w.loc[is_clinvar & stars.ge(3)] = 1.20
    w = pd.Series(np.maximum(w, np.where(as_bool(df["in_emerge"]), 1.00, 0.0)), index=df.index)
    w = pd.Series(np.maximum(w, np.where(as_bool(df["in_cardioboost"]), 1.00, 0.0)), index=df.index)
    w = pd.Series(np.maximum(w, np.where(as_bool(df["in_hiro"]), 1.40, 0.0)), index=df.index)
    return w.astype(float)


def class_weight(df: pd.DataFrame) -> pd.Series:
    # VUS is abundant and label-noisy; Pathogenic is rare and clinically costly.
    weights = df["model_label_3class"].map({"Benign": 1.0, "VUS": 0.90, "Pathogenic": 2.50})
    return weights.fillna(0.0).astype(float)


def exclusion_reason(df: pd.DataFrame) -> pd.Series:
    reason = pd.Series("", index=df.index, dtype="object")
    reason.loc[~df["primary_model_inclusion"].eq("include")] = "primary_model_inclusion_exclude"
    reason.loc[df["label_conflict_type"].fillna("none") != "none"] = "label_conflict"
    reason.loc[df["ref"].astype(str) == df["alt"].astype(str)] = "ref_equals_alt_noop"
    reason.loc[~df["model_label_3class"].isin(LABEL_TO_ID)] = "missing_or_unknown_label"
    return reason.replace("", "eligible")


def summarize_split(df: pd.DataFrame, split_col: str, eligible_col: str) -> dict:
    tab = pd.crosstab(df[split_col], df["model_label_3class"], dropna=False)
    out = {}
    for idx, row in tab.iterrows():
        out[str(idx)] = {str(k): int(v) for k, v in row.to_dict().items()}
        out[str(idx)]["rows"] = int((df[split_col] == idx).sum())
        out[str(idx)]["eligible_rows"] = int(((df[split_col] == idx) & df[eligible_col]).sum())
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", default=str(MATRIX))
    parser.add_argument("--out-matrix", default=str(OUT_MATRIX))
    parser.add_argument("--train-table", default=str(TRAIN_TABLE))
    parser.add_argument("--summary", default=str(SUMMARY))
    args = parser.parse_args()

    matrix_path = Path(args.matrix)
    out_matrix = Path(args.out_matrix)
    train_table = Path(args.train_table)
    summary_path = Path(args.summary)
    for path in [out_matrix, train_table, summary_path]:
        path.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(matrix_path, sep="\t", low_memory=False)
    clean_base = (
        df["primary_model_inclusion"].eq("include")
        & df["label_conflict_type"].fillna("none").eq("none")
        & (df["ref"].astype(str) != df["alt"].astype(str))
    )
    three_class = clean_base & df["model_label_3class"].isin(LABEL_TO_ID)

    df["training_slice_three_class_primary"] = three_class
    df["three_class_label"] = df["model_label_3class"].map(LABEL_TO_ID)
    df["three_class_training_exclude_reason"] = exclusion_reason(df)
    df["three_class_source_confidence_weight"] = source_confidence_weight(df)
    df["three_class_class_weight"] = class_weight(df)
    df["three_class_sample_weight"] = df["three_class_source_confidence_weight"] * df["three_class_class_weight"]
    df.loc[~three_class, "three_class_sample_weight"] = 0.0

    df["split_three_class_internal_grouped"] = assign_internal_grouped_split(df, three_class)
    df["split_three_class_source_heldout"] = assign_source_heldout_split(df, three_class)
    df["split_three_class_gene_stress"] = assign_gene_stress_split(df, three_class)

    df.to_csv(out_matrix, sep="\t", index=False)
    train_df = df[df["training_slice_three_class_primary"]].copy()
    train_df.to_csv(train_table, sep="\t", index=False)

    summary = {
        "input_matrix": rel(matrix_path),
        "output_matrix": rel(out_matrix),
        "three_class_table": rel(train_table),
        "rows_total": int(len(df)),
        "three_class_rows": int(three_class.sum()),
        "three_class_label_counts": {str(k): int(v) for k, v in train_df["model_label_3class"].value_counts().to_dict().items()},
        "weight_summary": {
            "min": float(train_df["three_class_sample_weight"].min()),
            "median": float(train_df["three_class_sample_weight"].median()),
            "mean": float(train_df["three_class_sample_weight"].mean()),
            "max": float(train_df["three_class_sample_weight"].max()),
        },
        "split_three_class_internal_grouped": summarize_split(df, "split_three_class_internal_grouped", "training_slice_three_class_primary"),
        "split_three_class_source_heldout": summarize_split(df, "split_three_class_source_heldout", "training_slice_three_class_primary"),
        "split_three_class_gene_stress": summarize_split(df, "split_three_class_gene_stress", "training_slice_three_class_primary"),
        "notes": [
            "Three-class model includes Benign, VUS, and Pathogenic rows.",
            "VUS is an evidence-status label, not a stable biological class; this model is exploratory.",
            "Sample weight = source confidence * class weight; Pathogenic is upweighted and VUS is mildly downweighted.",
        ],
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
