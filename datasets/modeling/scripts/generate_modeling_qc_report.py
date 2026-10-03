#!/usr/bin/env python3
"""Generate QC reports for the cardiogenetics modeling matrix."""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


BOOL_TRUE = {"true", "1", "yes", "y", "t"}


FEATURE_GROUPS = {
    "dbNSFP": [
        "dbnsfp_status",
        "sift_score",
        "polyphen2_hdiv_score",
        "revel_score",
        "metalr_score",
        "fathmm_xf_coding_score",
        "cadd_phred",
        "alphamissense_dbnsfp_score",
        "esm1b_score",
        "gerp_rs",
        "phylop100way_vertebrate",
        "phastcons100way_vertebrate",
    ],
    "AlphaMissense direct": [
        "alphamissense_direct_status",
        "alphamissense_direct_score",
        "alphamissense_direct_class",
    ],
    "gnomAD": [
        "gnomad_final_status",
        "gnomad_final_af",
        "gnomad_final_popmax_af",
        "gnomad_final_homozygote_count",
        "gnomad_browser_status",
    ],
    "ClinGen": [
        "clingen_gene_validity_status",
        "clingen_dosage_status",
        "clingen_variant_evidence_status",
        "clingen_gene_validity_max_classification",
        "clingen_variant_evidence_classification",
    ],
    "Protein structure": [
        "protein_feature_status",
        "alphafold_plddt_status",
        "alphafold_residue_plddt",
        "dssp_status",
        "dssp_secondary_structure_class",
        "freesasa_status",
        "freesasa_relative",
        "foldx_ddg_status",
        "foldx_ddg_kcal_mol",
    ],
    "VEP": [
        "vep_status",
        "vep_annotation_status",
        "vep_worst_consequence",
        "vep_impact",
        "vep_hgvsc",
        "vep_hgvsp",
    ],
    "SpliceAI": [
        "spliceai_status",
        "vep_SpliceAI_pred_DS_AG",
        "vep_SpliceAI_pred_DS_AL",
        "vep_SpliceAI_pred_DS_DG",
        "vep_SpliceAI_pred_DS_DL",
    ],
}


SOURCE_RECORD_FILES = {
    "HiRO": Path("datasets/hiro/full_dataset/interim/hiro_linked_source_records.tsv"),
    "eMERGE": Path("datasets/emerge/full_arrhythmia_gene_dataset/data/emerge_linked_source_records.tsv"),
    "CardioBoost": Path(
        "datasets/cardioboost/public_dataset/processed_for_cardiogenetics/model_inputs/"
        "cardioboost_linked_source_records.tsv"
    ),
}


def is_missing(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.isna()
    text = series.astype("string")
    return series.isna() | text.str.strip().isin(["", ".", "NA", "NaN", "nan", "None", "null"])


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype("string").str.lower().isin(BOOL_TRUE)


def value_counts_frame(
    df: pd.DataFrame, column: str, name: str | None = None, top_n: int | None = None
) -> pd.DataFrame:
    if column not in df.columns:
        return pd.DataFrame(columns=[column, "rows", "percent"])
    counts = df[column].fillna("<missing>").astype(str).value_counts(dropna=False)
    if top_n is not None:
        counts = counts.head(top_n)
    out = counts.rename_axis(name or column).reset_index(name="rows")
    out["percent"] = (out["rows"] / len(df) * 100).round(3) if len(df) else 0
    return out


def grouped_label_counts(df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    if group_col not in df.columns or "model_label_3class" not in df.columns:
        return pd.DataFrame()
    table = (
        df.groupby([group_col, "model_label_3class"], dropna=False)
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    for col in ["Benign", "Pathogenic", "VUS"]:
        if col not in table.columns:
            table[col] = 0
    table["labeled_3class_rows"] = table[["Benign", "Pathogenic", "VUS"]].sum(axis=1)
    table["supervised_binary_rows"] = table["Benign"] + table["Pathogenic"]
    denom = table["supervised_binary_rows"].astype(float).replace({0.0: np.nan})
    table["pathogenic_fraction_binary"] = (table["Pathogenic"].astype(float) / denom).round(4)
    return table.sort_values("labeled_3class_rows", ascending=False)


def feature_missingness(df: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    rows = []
    for col in columns:
        if col not in df.columns:
            rows.append(
                {
                    "column": col,
                    "present": False,
                    "nonmissing": 0,
                    "missing": len(df),
                    "missing_percent": 100.0,
                    "unique_nonmissing": 0,
                }
            )
            continue
        miss = is_missing(df[col])
        rows.append(
            {
                "column": col,
                "present": True,
                "nonmissing": int((~miss).sum()),
                "missing": int(miss.sum()),
                "missing_percent": round(float(miss.mean() * 100), 3),
                "unique_nonmissing": int(df.loc[~miss, col].nunique(dropna=True)),
            }
        )
    return pd.DataFrame(rows)


def positive_status_mask(series: pd.Series) -> pd.Series:
    values = series.astype("string").str.lower()
    return values.isin(["ok", "observed", "confirmed_absent", "linked_to_registry"])


def source_record_qc(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    summaries = []
    dup_rows = []
    for source_name, rel_path in SOURCE_RECORD_FILES.items():
        path = root / rel_path
        if not path.exists():
            summaries.append({"source": source_name, "path": str(rel_path), "status": "missing"})
            continue
        sdf = pd.read_csv(path, sep="\t", low_memory=False)
        link_col = "hiro_registry_link_status" if "hiro_registry_link_status" in sdf.columns else "source_registry_link_status"
        distinct_col = (
            "is_distinct_hiro_source_record"
            if "is_distinct_hiro_source_record" in sdf.columns
            else "is_distinct_source_record"
        )
        variant_col = "resolved_variant_id"
        linked = int((sdf[link_col].astype("string") == "linked_to_registry").sum()) if link_col in sdf else 0
        unresolved = int(sdf[variant_col].isna().sum()) if variant_col in sdf else 0
        distinct = int(as_bool(sdf[distinct_col]).sum()) if distinct_col in sdf else len(sdf)
        summaries.append(
            {
                "source": source_name,
                "path": str(rel_path),
                "status": "ok",
                "source_rows": len(sdf),
                "distinct_source_records": distinct,
                "linked_source_rows": linked,
                "unresolved_source_rows": unresolved,
                "unique_resolved_variants": int(sdf[variant_col].nunique(dropna=True)) if variant_col in sdf else 0,
            }
        )
        if variant_col in sdf:
            dups = (
                sdf.groupby(variant_col, dropna=True)
                .agg(
                    source_rows=(variant_col, "size"),
                    labels=("target_3class", lambda x: "|".join(sorted(set(x.dropna().astype(str)))) if "target_3class" in sdf else ""),
                    genes=("gene", lambda x: "|".join(sorted(set(x.dropna().astype(str)))) if "gene" in sdf else ""),
                )
                .reset_index()
            )
            dups = dups[dups["source_rows"] > 1].sort_values("source_rows", ascending=False)
            dups.insert(0, "source", source_name)
            dup_rows.append(dups.head(100))
    return pd.DataFrame(summaries), pd.concat(dup_rows, ignore_index=True) if dup_rows else pd.DataFrame()


def write_table(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False)


def md_table(df: pd.DataFrame, max_rows: int = 20) -> str:
    if df.empty:
        return "_No rows._"
    clipped = df.head(max_rows).copy()
    clipped = clipped.fillna("")
    columns = [str(c) for c in clipped.columns]

    def clean_cell(value: object) -> str:
        text = str(value)
        return text.replace("|", "\\|").replace("\n", " ")

    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for _, row in clipped.iterrows():
        lines.append("| " + " | ".join(clean_cell(row[c]) for c in clipped.columns) + " |")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--matrix",
        default="datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx_vep_hgvs_spliceai.tsv",
    )
    parser.add_argument("--out-dir", default=None)
    args = parser.parse_args()

    root = Path.cwd()
    matrix_path = root / args.matrix
    if not matrix_path.exists():
        raise FileNotFoundError(matrix_path)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out_dir) if args.out_dir else root / "datasets/modeling/qc" / f"modeling_matrix_qc_{stamp}"
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(matrix_path, sep="\t", low_memory=False)

    summary = {
        "matrix": str(matrix_path.relative_to(root)),
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "rows": int(len(df)),
        "columns": int(len(df.columns)),
        "duplicate_variant_id_rows": int(df.duplicated("variant_id").sum()) if "variant_id" in df else None,
        "duplicate_coordinate_key_rows": int(df.duplicated("coordinate_key").sum()) if "coordinate_key" in df else None,
        "ref_alt_equal_rows": int((df["ref"].astype(str) == df["alt"].astype(str)).sum())
        if {"ref", "alt"}.issubset(df.columns)
        else None,
    }

    source_cols = [c for c in ["in_clinvar", "in_hiro", "in_emerge", "in_cardioboost"] if c in df.columns]
    for col in source_cols:
        summary[f"{col}_rows"] = int(as_bool(df[col]).sum())

    if "model_label_3class" in df:
        labels = df["model_label_3class"].fillna("<missing>").astype(str).value_counts().to_dict()
        summary["label_counts"] = {k: int(v) for k, v in labels.items()}
        summary["binary_supervised_rows"] = int(df["model_label_3class"].isin(["Benign", "Pathogenic"]).sum())
    if "primary_model_inclusion" in df:
        include = df["primary_model_inclusion"].astype("string") == "include"
        summary["primary_model_include_rows"] = int(include.sum())
        if "model_label_3class" in df:
            summary["primary_include_binary_supervised_rows"] = int(
                (include & df["model_label_3class"].isin(["Benign", "Pathogenic"])).sum()
            )

    tables: dict[str, pd.DataFrame] = {}
    tables["label_counts"] = value_counts_frame(df, "model_label_3class")
    tables["clean_supervised_label_counts"] = value_counts_frame(df, "clean_supervised_label")
    tables["source_overlap_counts"] = value_counts_frame(df, "sources")
    tables["default_split_role_counts"] = value_counts_frame(df, "default_split_role")
    tables["modeling_scope_counts"] = value_counts_frame(df, "modeling_scope")
    tables["gene_panel_decision_counts"] = value_counts_frame(df, "gene_panel_decisions")
    tables["label_conflict_counts"] = value_counts_frame(df, "label_conflict_type")
    tables["gnomad_final_status_counts"] = value_counts_frame(df, "gnomad_final_status")
    tables["vep_status_counts"] = value_counts_frame(df, "vep_status")
    tables["spliceai_status_counts"] = value_counts_frame(df, "spliceai_status")
    tables["foldx_status_counts"] = value_counts_frame(df, "foldx_ddg_status")
    tables["clingen_gene_validity_status_counts"] = value_counts_frame(df, "clingen_gene_validity_status")
    tables["clingen_variant_evidence_status_counts"] = value_counts_frame(df, "clingen_variant_evidence_status")
    tables["protein_feature_status_counts"] = value_counts_frame(df, "protein_feature_status")
    tables["dssp_status_counts"] = value_counts_frame(df, "dssp_status")
    tables["freesasa_status_counts"] = value_counts_frame(df, "freesasa_status")

    for group_col in ["sources", "primary_gene", "modeling_scope", "default_split_role", "gene_panel_decisions"]:
        tables[f"label_counts_by_{group_col}"] = grouped_label_counts(df, group_col)

    if "primary_gene" in df:
        gene_counts = grouped_label_counts(df, "primary_gene")
        if not gene_counts.empty:
            gene_counts["sparse_pathogenic_lt30"] = gene_counts["Pathogenic"] < 30
            tables["gene_label_counts"] = gene_counts
            tables["sparse_pathogenic_genes"] = gene_counts[
                (gene_counts["Pathogenic"] < 30) & (gene_counts["supervised_binary_rows"] > 0)
            ].sort_values(["Pathogenic", "supervised_binary_rows"], ascending=[True, False])

    if {"variant_id", "chrom", "pos", "ref", "alt"}.issubset(df.columns):
        tables["duplicate_variant_id_rows"] = df[df.duplicated("variant_id", keep=False)].sort_values("variant_id")
        tables["ref_alt_equal_rows"] = df[df["ref"].astype(str) == df["alt"].astype(str)][
            [
                c
                for c in [
                    "variant_id",
                    "chrom",
                    "pos",
                    "ref",
                    "alt",
                    "primary_gene",
                    "sources",
                    "model_label_3class",
                    "vep_status",
                    "vep_missing_reason",
                ]
                if c in df.columns
            ]
        ]

    conflict_cols = [
        c
        for c in [
            "variant_id",
            "primary_gene",
            "sources",
            "model_label_3class",
            "label_conflict_type",
            "clinvar_labels",
            "hiro_labels",
            "emerge_labels",
            "cardioboost_labels",
            "source_row_count",
        ]
        if c in df.columns
    ]
    if "label_conflict_type" in df:
        conflicts = df[df["label_conflict_type"].fillna("none") != "none"][conflict_cols]
        tables["label_conflict_rows"] = conflicts.sort_values(["label_conflict_type", "primary_gene", "variant_id"])

    high_source_cols = [
        c
        for c in [
            "variant_id",
            "primary_gene",
            "sources",
            "source_count",
            "source_row_count",
            "model_label_3class",
            "label_conflict_type",
            "clinvar_labels",
            "hiro_labels",
            "emerge_labels",
            "cardioboost_labels",
        ]
        if c in df.columns
    ]
    if "source_row_count" in df:
        high_source = df[df["source_row_count"].fillna(0).astype(float) > 1][high_source_cols]
        tables["multi_source_record_variant_rows"] = high_source.sort_values("source_row_count", ascending=False)

    status_cols = [c for c in df.columns if c.endswith("_status") or c.endswith("_missing_reason")]
    core_feature_cols = sorted(set(sum(FEATURE_GROUPS.values(), [])))
    tables["core_feature_missingness"] = feature_missingness(df, core_feature_cols)
    tables["all_status_missingness"] = feature_missingness(df, status_cols)

    feature_group_rows = []
    for group, cols in FEATURE_GROUPS.items():
        present_cols = [c for c in cols if c in df.columns]
        status_cols_for_group = [c for c in present_cols if c.endswith("_status")]
        value_cols_for_group = [
            c
            for c in present_cols
            if not c.endswith("_status") and not c.endswith("_missing_reason") and not c.endswith("_class")
        ]
        if not present_cols:
            feature_group_rows.append(
                {
                    "feature_group": group,
                    "columns_present": 0,
                    "status_columns": 0,
                    "value_columns": 0,
                    "rows_with_positive_status": 0,
                    "positive_status_percent": 0.0,
                    "rows_with_any_value": 0,
                    "value_coverage_percent": 0.0,
                }
            )
            continue
        positive_status = pd.Series(False, index=df.index)
        for col in status_cols_for_group:
            positive_status |= positive_status_mask(df[col])
        value_present = pd.Series(False, index=df.index)
        for col in value_cols_for_group:
            value_present |= ~is_missing(df[col])
        feature_group_rows.append(
            {
                "feature_group": group,
                "columns_present": len(present_cols),
                "status_columns": len(status_cols_for_group),
                "value_columns": len(value_cols_for_group),
                "rows_with_positive_status": int(positive_status.sum()),
                "positive_status_percent": round(float(positive_status.mean() * 100), 3),
                "rows_with_any_value": int(value_present.sum()),
                "value_coverage_percent": round(float(value_present.mean() * 100), 3),
            }
        )
    tables["feature_group_coverage"] = pd.DataFrame(feature_group_rows)

    source_summary, source_dups = source_record_qc(root)
    tables["linked_source_record_summary"] = source_summary
    tables["linked_source_record_duplicate_variants"] = source_dups

    for name, table in tables.items():
        write_table(table, out_dir / f"{name}.tsv")

    with (out_dir / "summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2)

    report = []
    report.append("# Modeling Matrix QC Report")
    report.append("")
    report.append(f"Generated: `{summary['generated_at']}`")
    report.append(f"Matrix: `{summary['matrix']}`")
    report.append("")
    report.append("## Executive Summary")
    report.append("")
    report.append(f"- Rows: `{summary['rows']:,}`")
    report.append(f"- Columns: `{summary['columns']:,}`")
    report.append(f"- Duplicate `variant_id` rows: `{summary['duplicate_variant_id_rows']}`")
    report.append(f"- Duplicate `coordinate_key` rows: `{summary['duplicate_coordinate_key_rows']}`")
    report.append(f"- `ref == alt` no-op rows: `{summary['ref_alt_equal_rows']}`")
    report.append(f"- Binary supervised rows, before final split decisions: `{summary.get('binary_supervised_rows', 0):,}`")
    report.append(
        f"- Primary-include binary supervised rows: `{summary.get('primary_include_binary_supervised_rows', 0):,}`"
    )
    report.append("")
    report.append("## Label Balance")
    report.append("")
    report.append(md_table(tables["label_counts"]))
    report.append("")
    report.append("## Source Coverage")
    report.append("")
    report.append(md_table(tables["source_overlap_counts"], 20))
    report.append("")
    report.append("## Source Record Preservation")
    report.append("")
    report.append(md_table(tables["linked_source_record_summary"], 20))
    report.append("")
    report.append("## Split Roles And Scope")
    report.append("")
    report.append(md_table(tables["default_split_role_counts"], 20))
    report.append("")
    report.append(md_table(tables["modeling_scope_counts"], 20))
    report.append("")
    report.append("## Label Conflicts")
    report.append("")
    report.append(md_table(tables["label_conflict_counts"], 20))
    report.append("")
    report.append("Top conflict rows are written to `label_conflict_rows.tsv`.")
    report.append("")
    report.append("## Feature Group Coverage")
    report.append("")
    report.append(md_table(tables["feature_group_coverage"], 20))
    report.append("")
    report.append("## Core Feature Missingness")
    report.append("")
    report.append(md_table(tables["core_feature_missingness"], 80))
    report.append("")
    report.append("## Top Genes By Label Count")
    report.append("")
    report.append(md_table(tables.get("gene_label_counts", pd.DataFrame()), 40))
    report.append("")
    report.append("## Sparse Pathogenic Genes")
    report.append("")
    report.append(md_table(tables.get("sparse_pathogenic_genes", pd.DataFrame()), 60))
    report.append("")
    report.append("## Main QC Flags")
    report.append("")
    report.append("- Rows with `label_conflict_type != none` should be excluded from clean supervised training or reviewed.")
    report.append("- `ref == alt` rows explain the VEP-missing rows and should not be treated as annotation failures.")
    report.append("- HiRO/eMERGE/CardioBoost source records are preserved separately from the variant-level matrix.")
    report.append("- `gnomad_final_status == confirmed_absent` is distinct from `not_joined`; only confirmed absent can support an AF=0 indicator.")
    report.append("- FoldX `ref_mismatch` rows should keep DDG missing and use the mismatch flag rather than imputation.")
    report.append("")
    report.append("## Output Files")
    report.append("")
    for name in sorted(tables):
        report.append(f"- `{name}.tsv`")
    report.append("- `summary.json`")

    (out_dir / "QC_REPORT.md").write_text("\n".join(report) + "\n")
    print(out_dir)


if __name__ == "__main__":
    main()
