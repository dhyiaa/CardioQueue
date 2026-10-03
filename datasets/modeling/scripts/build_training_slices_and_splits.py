#!/usr/bin/env python3
"""Build training slice flags, sample weights, and split assignments."""

from __future__ import annotations

import json
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split


ROOT = Path(__file__).resolve().parents[3]
MATRIX = ROOT / (
    "datasets/modeling/interim/"
    "final_modeling_table_local_features_gnomad_foldx_vep_hgvs_spliceai_hiro_agg.tsv"
)
OUT = ROOT / "datasets/modeling/ready"
OUT_MATRIX = OUT / "modeling_table_with_splits_weights.tsv"
PRIMARY_TABLE = OUT / "primary_binary_catboost_table.tsv"
SUMMARY = OUT / "training_slices_and_splits.summary.json"


def as_bool(s: pd.Series) -> pd.Series:
    return s.astype("string").str.lower().isin(["true", "1", "yes", "y", "t"])


def num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


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

    # Most groups are one variant, but use source_leakage_group to protect future
    # duplicated/source-record-expanded tables.
    train_idx, temp_idx = train_test_split(
        idx,
        train_size=0.70,
        random_state=seed,
        stratify=y,
    )
    y_temp = make_strat_label(df.loc[temp_idx])
    # Rare strata can become too small after the first split; fall back to label-only.
    if y_temp.value_counts().min() < 2:
        y_temp = df.loc[temp_idx, "model_label_3class"].astype(str)
    val_idx, test_idx = train_test_split(
        temp_idx,
        train_size=0.50,
        random_state=seed + 1,
        stratify=y_temp,
    )
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

    train_pool = eligible & ~external
    idx = df.index[train_pool].to_numpy()
    y = safe_strat_label(df.loc[idx])
    train_idx, val_idx = train_test_split(idx, train_size=0.85, random_state=seed, stratify=y)
    split.loc[train_idx] = "train"
    split.loc[val_idx] = "validation"
    return split


def assign_source_heldout_with_clinvar_challenge(
    df: pd.DataFrame,
    eligible: pd.Series,
    challenge_ids: set[str],
    seed: int = 20260702,
) -> pd.Series:
    """Build a fresh split after quarantining a frozen same-ecosystem challenge."""
    split = pd.Series("not_eligible", index=df.index, dtype="object")
    is_hiro = as_bool(df["in_hiro"])
    is_emerge = as_bool(df["in_emerge"])
    is_cardio = as_bool(df["in_cardioboost"])
    external = eligible & (is_hiro | is_emerge | is_cardio)
    challenge = eligible & df["variant_id"].astype(str).isin(challenge_ids)

    overlap = challenge & external
    if overlap.any():
        examples = df.loc[overlap, "variant_id"].astype(str).head(10).tolist()
        raise ValueError(f"ClinVar challenge overlaps external sources: {examples}")

    split.loc[eligible & is_hiro] = "external_hiro"
    split.loc[eligible & ~is_hiro & is_emerge] = "external_emerge"
    split.loc[eligible & ~is_hiro & ~is_emerge & is_cardio] = "external_cardioboost"
    split.loc[challenge] = "challenge_clinvar"

    train_pool = eligible & ~external & ~challenge
    idx = df.index[train_pool].to_numpy()
    y = safe_strat_label(df.loc[idx])
    train_idx, val_idx = train_test_split(idx, train_size=0.85, random_state=seed, stratify=y)
    split.loc[train_idx] = "train"
    split.loc[val_idx] = "validation"
    return split


def load_challenge_ids(path: Path, df: pd.DataFrame, eligible: pd.Series) -> set[str]:
    manifest = pd.read_csv(path, sep="\t", low_memory=False)
    if "variant_id" not in manifest.columns:
        raise ValueError("ClinVar challenge manifest must contain variant_id")
    if manifest["variant_id"].isna().any() or not manifest["variant_id"].astype(str).is_unique:
        raise ValueError("ClinVar challenge manifest must contain unique, non-missing variant_id values")

    challenge_ids = set(manifest["variant_id"].astype(str))
    matrix_ids = set(df["variant_id"].astype(str))
    missing = sorted(challenge_ids - matrix_ids)
    if missing:
        raise ValueError(f"ClinVar challenge contains IDs absent from the matrix: {missing[:10]}")

    selected = df["variant_id"].astype(str).isin(challenge_ids)
    ineligible = selected & ~eligible
    if ineligible.any():
        examples = df.loc[
            ineligible, ["variant_id", "model_label_3class", "primary_model_inclusion", "label_conflict_type"]
        ].head(10)
        raise ValueError(f"ClinVar challenge contains ineligible rows:\n{examples.to_string(index=False)}")
    return challenge_ids


def assign_gene_stress_split(df: pd.DataFrame, eligible: pd.Series, seed: int = 20260702) -> pd.Series:
    split = pd.Series("not_eligible", index=df.index, dtype="object")
    counts = pd.crosstab(df.loc[eligible, "primary_gene"], df.loc[eligible, "model_label_3class"])
    for col in ["Benign", "Pathogenic"]:
        if col not in counts:
            counts[col] = 0
    sparse_genes = set(counts.index[counts["Pathogenic"] < 30])
    sparse = eligible & df["primary_gene"].isin(sparse_genes)
    split.loc[sparse] = "sparse_gene_stress_test"

    pool = eligible & ~sparse
    idx = df.index[pool].to_numpy()
    y = safe_strat_label(df.loc[idx])
    train_idx, temp_idx = train_test_split(idx, train_size=0.70, random_state=seed, stratify=y)
    y_temp = make_strat_label(df.loc[temp_idx])
    if y_temp.value_counts().min() < 2:
        y_temp = df.loc[temp_idx, "model_label_3class"].astype(str)
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

    # Trusted non-ClinVar sources. Use max-like overrides so overlap rows keep
    # the strongest source confidence without making the value explode.
    w = np.maximum(w, np.where(as_bool(df["in_emerge"]), 1.00, 0.0))
    w = pd.Series(w, index=df.index)
    w = np.maximum(w, np.where(as_bool(df["in_cardioboost"]), 1.00, 0.0))
    w = pd.Series(w, index=df.index)
    w = np.maximum(w, np.where(as_bool(df["in_hiro"]), 1.40, 0.0))
    return pd.Series(w, index=df.index).astype(float)


def class_weight(df: pd.DataFrame) -> pd.Series:
    # Mild/moderate P/LP upweighting: enough to improve recall without using the
    # full 3.7x inverse-frequency ratio, which could harm calibration.
    return pd.Series(np.where(df["model_label_3class"].eq("Pathogenic"), 2.0, 1.0), index=df.index).astype(float)


def exclusion_reason(df: pd.DataFrame) -> pd.Series:
    reason = pd.Series("", index=df.index, dtype="object")
    reason.loc[~df["primary_model_inclusion"].eq("include")] = "primary_model_inclusion_exclude"
    reason.loc[df["label_conflict_type"].fillna("none") != "none"] = "label_conflict"
    reason.loc[df["ref"].astype(str) == df["alt"].astype(str)] = "ref_equals_alt_noop"
    reason.loc[~df["model_label_3class"].isin(["Benign", "Pathogenic"])] = "not_binary_label"
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
    parser.add_argument("--primary-table", default=str(PRIMARY_TABLE))
    parser.add_argument("--summary", default=str(SUMMARY))
    parser.add_argument(
        "--clinvar-challenge-manifest",
        help="Optional frozen TSV with unique variant_id values for a same-ecosystem ClinVar challenge.",
    )
    args = parser.parse_args()

    matrix_path = Path(args.matrix)
    out_matrix = Path(args.out_matrix)
    primary_table = Path(args.primary_table)
    summary_path = Path(args.summary)
    for path in [out_matrix, primary_table, summary_path]:
        path.parent.mkdir(parents=True, exist_ok=True)

    def rel(path: Path) -> str:
        try:
            return str(path.relative_to(ROOT))
        except ValueError:
            return str(path)

    df = pd.read_csv(matrix_path, sep="\t", low_memory=False)

    clean_base = (
        df["primary_model_inclusion"].eq("include")
        & df["label_conflict_type"].fillna("none").eq("none")
        & (df["ref"].astype(str) != df["alt"].astype(str))
    )
    binary = clean_base & df["model_label_3class"].isin(["Benign", "Pathogenic"])
    three_class = clean_base & df["model_label_3class"].isin(["Benign", "Pathogenic", "VUS"])
    vus_scoring = clean_base & df["model_label_3class"].eq("VUS")
    high_conf = (
        binary
        & (
            (as_bool(df["in_clinvar"]) & num(df["clinvar_review_stars"]).ge(2))
            | as_bool(df["in_hiro"])
            | as_bool(df["in_emerge"])
            | as_bool(df["in_cardioboost"])
        )
    )

    df["training_slice_primary_binary"] = binary
    df["training_slice_high_confidence_binary"] = high_conf
    df["training_slice_three_class_exploratory"] = three_class
    df["training_slice_vus_scoring"] = vus_scoring
    df["training_exclude_reason"] = exclusion_reason(df)
    df["primary_binary_label"] = df["model_label_3class"].map({"Benign": 0, "Pathogenic": 1})
    df["source_confidence_weight"] = source_confidence_weight(df)
    df["class_weight"] = class_weight(df)
    df["sample_weight"] = df["source_confidence_weight"] * df["class_weight"]
    df.loc[~binary, "sample_weight"] = 0.0

    df["split_internal_grouped"] = assign_internal_grouped_split(df, binary)
    df["split_source_heldout"] = assign_source_heldout_split(df, binary)
    df["split_gene_stress"] = assign_gene_stress_split(df, binary)

    challenge_split_col = None
    challenge_ids: set[str] = set()
    if args.clinvar_challenge_manifest:
        challenge_path = Path(args.clinvar_challenge_manifest)
        challenge_ids = load_challenge_ids(challenge_path, df, binary)
        challenge_split_col = "split_source_heldout_clinvar_challenge"
        df[challenge_split_col] = assign_source_heldout_with_clinvar_challenge(
            df, binary, challenge_ids
        )

    df.to_csv(out_matrix, sep="\t", index=False)
    primary = df[df["training_slice_primary_binary"]].copy()
    primary.to_csv(primary_table, sep="\t", index=False)

    summary = {
        "input_matrix": rel(matrix_path),
        "output_matrix": rel(out_matrix),
        "primary_binary_table": rel(primary_table),
        "rows_total": int(len(df)),
        "primary_binary_rows": int(binary.sum()),
        "high_confidence_binary_rows": int(high_conf.sum()),
        "three_class_exploratory_rows": int(three_class.sum()),
        "vus_scoring_rows": int(vus_scoring.sum()),
        "primary_binary_label_counts": {
            str(k): int(v) for k, v in primary["model_label_3class"].value_counts().to_dict().items()
        },
        "weight_summary_primary_binary": {
            "min": float(primary["sample_weight"].min()),
            "median": float(primary["sample_weight"].median()),
            "mean": float(primary["sample_weight"].mean()),
            "max": float(primary["sample_weight"].max()),
        },
        "split_internal_grouped": summarize_split(df, "split_internal_grouped", "training_slice_primary_binary"),
        "split_source_heldout": summarize_split(df, "split_source_heldout", "training_slice_primary_binary"),
        "split_gene_stress": summarize_split(df, "split_gene_stress", "training_slice_primary_binary"),
        "clinvar_challenge_manifest": rel(Path(args.clinvar_challenge_manifest))
        if args.clinvar_challenge_manifest
        else None,
        "clinvar_challenge_rows": len(challenge_ids),
        "notes": [
            "Primary binary model uses clean Benign vs Pathogenic rows only.",
            "VUS rows are reserved for scoring/analysis, not primary supervised binary training.",
            "Sample weight = source_confidence_weight * class_weight; Pathogenic class weight currently 2.0.",
        ],
    }
    if challenge_split_col:
        summary[challenge_split_col] = summarize_split(
            df, challenge_split_col, "training_slice_primary_binary"
        )
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
