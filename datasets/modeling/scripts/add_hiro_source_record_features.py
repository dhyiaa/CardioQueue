#!/usr/bin/env python3
"""Add HiRO source-record aggregation features to the variant-level matrix.

The final modeling matrix is one row per genomic variant. HiRO, however, is
patient/source-record-level evidence: multiple rows can map to the same variant
because different patients or studies carry the same genotype. This script keeps
the variant-level matrix leakage-safe while adding aggregate columns that
preserve the 482 HiRO source-record events.
"""

from __future__ import annotations

import json
import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
MATRIX = ROOT / "datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx_vep_hgvs_spliceai.tsv"
HIRO = ROOT / "datasets/hiro/full_dataset/interim/hiro_linked_source_records.tsv"

AGG_OUT = ROOT / "datasets/modeling/interim/hiro_source_record_aggregation_features.tsv"
MATRIX_OUT = (
    ROOT
    / "datasets/modeling/interim/"
    "final_modeling_table_local_features_gnomad_foldx_vep_hgvs_spliceai_hiro_agg.tsv"
)
SUMMARY_OUT = (
    ROOT
    / "datasets/modeling/interim/"
    "final_modeling_table_local_features_gnomad_foldx_vep_hgvs_spliceai_hiro_agg.summary.json"
)
SOURCE_RECORD_OUT = ROOT / "datasets/modeling/interim/hiro_source_record_level_table.tsv"


BOOL_TRUE = {"true", "1", "yes", "y", "t"}


def clean_set(series: pd.Series) -> str:
    vals = sorted({str(x) for x in series.dropna() if str(x).strip() and str(x).strip().lower() != "nan"})
    return "|".join(vals)


def truthy(series: pd.Series) -> pd.Series:
    return series.astype("string").str.lower().isin(BOOL_TRUE)


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def count_positive(series: pd.Series) -> int:
    if series.empty:
        return 0
    if series.dtype == bool:
        return int(series.fillna(False).sum())
    text = series.astype("string").str.lower()
    return int(text.isin(["true", "1", "yes", "y", "present", "positive"]).sum())


def any_positive(series: pd.Series) -> bool:
    return count_positive(series) > 0


def build_hiro_aggregation(hiro: pd.DataFrame) -> pd.DataFrame:
    linked = hiro[hiro["resolved_variant_id"].notna()].copy()
    linked["resolved_variant_id"] = linked["resolved_variant_id"].astype(str)

    base = (
        linked.groupby("resolved_variant_id")
        .agg(
            hiro_source_record_count=("variant_uid", "size"),
            hiro_unique_patient_count=("participant_id", pd.Series.nunique),
            hiro_unique_study_count=("study", pd.Series.nunique),
            hiro_studies=("study", clean_set),
            hiro_gene_set=("gene", clean_set),
            hiro_label_set=("target_3class", clean_set),
            hiro_raw_label_set=("target_classification_raw", clean_set),
            hiro_5class_label_set=("target_5class", clean_set),
            hiro_patient_match_count=("patient_match", count_positive),
            hiro_any_patient_match=("patient_match", any_positive),
        )
        .reset_index()
        .rename(columns={"resolved_variant_id": "variant_id"})
    )

    label_counts = (
        linked.assign(target_3class=linked["target_3class"].fillna("Missing").astype(str))
        .pivot_table(
            index="resolved_variant_id",
            columns="target_3class",
            values="variant_uid",
            aggfunc="size",
            fill_value=0,
        )
        .reset_index()
        .rename(columns={"resolved_variant_id": "variant_id"})
    )
    for label in ["Benign", "Pathogenic", "VUS", "Missing"]:
        if label not in label_counts.columns:
            label_counts[label] = 0
    label_counts = label_counts.rename(
        columns={
            "Benign": "hiro_benign_record_count",
            "Pathogenic": "hiro_pathogenic_record_count",
            "VUS": "hiro_vus_record_count",
            "Missing": "hiro_missing_label_record_count",
        }
    )

    agg = base.merge(label_counts, on="variant_id", how="left")
    for col in [
        "hiro_benign_record_count",
        "hiro_pathogenic_record_count",
        "hiro_vus_record_count",
        "hiro_missing_label_record_count",
    ]:
        agg[col] = agg[col].fillna(0).astype(int)
    agg["hiro_nonmissing_label_class_count"] = (
        (agg["hiro_benign_record_count"] > 0).astype(int)
        + (agg["hiro_pathogenic_record_count"] > 0).astype(int)
        + (agg["hiro_vus_record_count"] > 0).astype(int)
    )
    agg["hiro_label_discordant_flag"] = agg["hiro_nonmissing_label_class_count"] > 1
    agg["hiro_duplicate_source_record_flag"] = agg["hiro_source_record_count"] > 1
    agg["hiro_has_casper_wes"] = agg["hiro_studies"].astype(str).str.contains("CASPER_WES", regex=False)
    agg["hiro_has_verdict"] = agg["hiro_studies"].astype(str).str.contains("VERDICT", regex=False)
    agg["hiro_pathogenic_record_fraction"] = (
        agg["hiro_pathogenic_record_count"] / agg["hiro_source_record_count"].replace({0: pd.NA})
    ).astype(float)
    agg["hiro_benign_record_fraction"] = (
        agg["hiro_benign_record_count"] / agg["hiro_source_record_count"].replace({0: pd.NA})
    ).astype(float)
    agg["hiro_vus_record_fraction"] = (
        agg["hiro_vus_record_count"] / agg["hiro_source_record_count"].replace({0: pd.NA})
    ).astype(float)

    symptom_cols = [
        "patient_sym_presyncope",
        "patient_sym_syncope",
        "patient_sym_palpitations",
        "patient_sym_chest_pain",
        "patient_sym_cardiac_arrest",
        "patient_sym_death",
        "patient_sym_other",
        "patient_family_history_flags",
        "patient_followup_presyncope_any",
        "patient_followup_syncope_any",
        "patient_followup_palpitations_any",
        "patient_followup_chest_pain_any",
        "patient_followup_cardiac_arrest_any",
        "patient_followup_death_any",
        "patient_followup_other_symptom_any",
    ]
    symptom_cols = [c for c in symptom_cols if c in linked.columns]
    for col in symptom_cols:
        feature = "hiro_any_" + col.replace("patient_", "").replace("_any", "")
        counts = linked.groupby("resolved_variant_id")[col].agg(any_positive).reset_index()
        counts = counts.rename(columns={"resolved_variant_id": "variant_id", col: feature})
        agg = agg.merge(counts, on="variant_id", how="left")

    count_cols = [c for c in linked.columns if c.startswith("patient_followup_") and c.endswith("_yes_count")]
    for col in count_cols:
        feature_sum = "hiro_" + col.replace("patient_", "") + "_sum"
        values = numeric(linked[col])
        temp = linked[["resolved_variant_id"]].copy()
        temp[col] = values
        counts = temp.groupby("resolved_variant_id")[col].sum(min_count=1).reset_index()
        counts = counts.rename(columns={"resolved_variant_id": "variant_id", col: feature_sum})
        agg = agg.merge(counts, on="variant_id", how="left")

    for col in ["patient_qt_min", "patient_qt_max", "patient_qt_count", "patient_lvef_min", "patient_lvef_max", "patient_lvef_count"]:
        if col not in linked.columns:
            continue
        temp = linked[["resolved_variant_id"]].copy()
        temp[col] = numeric(linked[col])
        grouped = temp.groupby("resolved_variant_id")[col]
        stat = "max" if col.endswith("_max") or col.endswith("_count") else "min"
        values = getattr(grouped, stat)().reset_index()
        values = values.rename(columns={"resolved_variant_id": "variant_id", col: "hiro_" + col.replace("patient_", "")})
        agg = agg.merge(values, on="variant_id", how="left")

    return agg


def build_source_record_level_table(hiro: pd.DataFrame, matrix: pd.DataFrame) -> pd.DataFrame:
    matrix_prefixed = matrix.add_prefix("matrix_")
    out = hiro.merge(
        matrix_prefixed,
        left_on="resolved_variant_id",
        right_on="matrix_variant_id",
        how="left",
    )
    out["matrix_join_status"] = out["matrix_variant_id"].notna().map({True: "joined_to_variant_matrix", False: "not_joined_to_variant_matrix"})
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", default=str(MATRIX))
    parser.add_argument("--hiro", default=str(HIRO))
    parser.add_argument("--agg-out", default=str(AGG_OUT))
    parser.add_argument("--matrix-out", default=str(MATRIX_OUT))
    parser.add_argument("--summary-out", default=str(SUMMARY_OUT))
    parser.add_argument("--source-record-out", default=str(SOURCE_RECORD_OUT))
    args = parser.parse_args()

    matrix_path = Path(args.matrix)
    hiro_path = Path(args.hiro)
    agg_out = Path(args.agg_out)
    matrix_out = Path(args.matrix_out)
    summary_out = Path(args.summary_out)
    source_record_out = Path(args.source_record_out)
    for path in [agg_out, matrix_out, summary_out, source_record_out]:
        path.parent.mkdir(parents=True, exist_ok=True)

    def rel(path: Path) -> str:
        try:
            return str(path.relative_to(ROOT))
        except ValueError:
            return str(path)

    matrix = pd.read_csv(matrix_path, sep="\t", low_memory=False)
    hiro = pd.read_csv(hiro_path, sep="\t", low_memory=False)

    agg = build_hiro_aggregation(hiro)
    agg.to_csv(agg_out, sep="\t", index=False)

    overlap_existing = [c for c in agg.columns if c in matrix.columns and c != "variant_id"]
    if overlap_existing:
        matrix = matrix.drop(columns=overlap_existing)

    merged = matrix.merge(agg, on="variant_id", how="left")
    count_cols = [c for c in merged.columns if c.startswith("hiro_") and c.endswith("_count")]
    for col in count_cols:
        merged[col] = merged[col].fillna(0).astype(int)
    bool_cols = [
        c
        for c in merged.columns
        if c.startswith("hiro_")
        and (c.endswith("_flag") or c.startswith("hiro_any_") or c in ["hiro_has_casper_wes", "hiro_has_verdict"])
    ]
    for col in bool_cols:
        merged[col] = merged[col].fillna(False).astype(bool)

    merged["hiro_source_record_features_available"] = merged["hiro_source_record_count"].fillna(0).astype(int) > 0
    merged.to_csv(matrix_out, sep="\t", index=False)

    source_record_table = build_source_record_level_table(hiro, matrix)
    source_record_table.to_csv(source_record_out, sep="\t", index=False)

    linked = hiro[hiro["resolved_variant_id"].notna()]
    summary = {
        "input_matrix": rel(matrix_path),
        "output_matrix": rel(matrix_out),
        "hiro_source_records": int(len(hiro)),
        "hiro_linked_or_resolved_source_records": int(len(linked)),
        "hiro_unique_resolved_variant_ids": int(linked["resolved_variant_id"].nunique()),
        "hiro_unresolved_no_coordinate_source_records": int((hiro["hiro_registry_link_status"] == "unresolved_no_coordinates").sum()),
        "hiro_aggregation_rows": int(len(agg)),
        "matrix_rows": int(len(merged)),
        "matrix_columns": int(len(merged.columns)),
        "matrix_rows_with_hiro_aggregates": int(merged["hiro_source_record_features_available"].sum()),
        "source_record_level_rows": int(len(source_record_table)),
        "source_record_level_joined_rows": int((source_record_table["matrix_join_status"] == "joined_to_variant_matrix").sum()),
        "source_record_level_not_joined_rows": int((source_record_table["matrix_join_status"] != "joined_to_variant_matrix").sum()),
        "notes": [
            "The variant-level matrix remains one row per variant for leakage-safe modeling.",
            "HiRO repeated patient/source events are preserved as aggregate features and in the source-record-level table.",
            "Unresolved HiRO rows still need coordinate rescue before variant-level feature annotation is possible.",
        ],
    }
    summary_out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
