#!/usr/bin/env python3
"""Create data-characteristics tables and figures for the cardiogenetics project."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "data_characteristics"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"

MATRIX = ROOT / (
    "datasets/modeling/interim/"
    "final_modeling_table_local_features_gnomad_foldx_vep_hgvs_spliceai_hiro_agg.tsv"
)
QC_DIR = ROOT / "datasets/modeling/qc/modeling_matrix_qc_20260702_hiro_agg"
HIRO_SOURCE = ROOT / "datasets/hiro/full_dataset/interim/hiro_linked_source_records.tsv"
EMERGE_SOURCE = ROOT / "datasets/emerge/full_arrhythmia_gene_dataset/data/emerge_linked_source_records.tsv"
CARDIOBOOST_SOURCE = (
    ROOT
    / "datasets/cardioboost/public_dataset/processed_for_cardiogenetics/model_inputs/cardioboost_linked_source_records.tsv"
)


sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update(
    {
        "figure.dpi": 160,
        "savefig.dpi": 240,
        "axes.titlesize": 15,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
    }
)


def ensure_dirs() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)


def read_matrix() -> pd.DataFrame:
    return pd.read_csv(MATRIX, sep="\t", low_memory=False)


def bool_col(df: pd.DataFrame, col: str) -> pd.Series:
    return df[col].astype("string").str.lower().isin(["true", "1", "yes"])


def missing(series: pd.Series) -> pd.Series:
    text = series.astype("string")
    return series.isna() | text.str.strip().isin(["", ".", "NA", "NaN", "nan", "None", "null"])


def pct(n: int | float, d: int | float) -> float:
    return round(float(n) / float(d) * 100, 3) if d else 0.0


def write_table(df: pd.DataFrame, name: str) -> Path:
    path = TABLES / name
    df.to_csv(path, sep="\t", index=False)
    return path


def short_md_table(df: pd.DataFrame, max_rows: int = 30) -> str:
    if df.empty:
        return "_No rows._"
    df = df.head(max_rows).fillna("")
    cols = [str(c) for c in df.columns]
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join(["---"] * len(cols)) + " |",
    ]
    for _, row in df.iterrows():
        vals = [str(row[c]).replace("|", "\\|").replace("\n", " ") for c in df.columns]
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def value_count_table(df: pd.DataFrame, col: str, name: str) -> pd.DataFrame:
    out = df[col].fillna("<missing>").astype(str).value_counts(dropna=False).rename_axis(name).reset_index(name="rows")
    out["percent"] = (out["rows"] / len(df) * 100).round(3)
    return out


def label_counts_by_source(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    source_cols = {
        "ClinVar": "in_clinvar",
        "HiRO": "in_hiro",
        "eMERGE": "in_emerge",
        "CardioBoost": "in_cardioboost",
    }
    for source, col in source_cols.items():
        mask = bool_col(df, col)
        subset = df.loc[mask]
        counts = subset["model_label_3class"].fillna("<missing>").astype(str).value_counts()
        rows.append(
            {
                "source": source,
                "variant_rows_in_matrix": int(mask.sum()),
                "benign": int(counts.get("Benign", 0)),
                "pathogenic": int(counts.get("Pathogenic", 0)),
                "vus": int(counts.get("VUS", 0)),
                "missing_or_conflict": int(counts.get("<missing>", 0)),
                "binary_supervised_rows": int(counts.get("Benign", 0) + counts.get("Pathogenic", 0)),
                "pathogenic_fraction_binary": round(
                    counts.get("Pathogenic", 0) / max(counts.get("Benign", 0) + counts.get("Pathogenic", 0), 1), 4
                ),
            }
        )
    return pd.DataFrame(rows)


def source_record_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    source_files = [
        ("HiRO", HIRO_SOURCE, "hiro_registry_link_status", "is_distinct_hiro_source_record", "in_hiro"),
        ("eMERGE", EMERGE_SOURCE, "source_registry_link_status", "is_distinct_source_record", "in_emerge"),
        ("CardioBoost", CARDIOBOOST_SOURCE, "source_registry_link_status", "is_distinct_source_record", "in_cardioboost"),
    ]
    for source, path, link_col, distinct_col, matrix_col in source_files:
        sdf = pd.read_csv(path, sep="\t", low_memory=False)
        linked = int((sdf[link_col].astype("string") == "linked_to_registry").sum())
        unresolved = int(sdf["resolved_variant_id"].isna().sum())
        unique_resolved = int(sdf["resolved_variant_id"].nunique(dropna=True))
        matrix_mask = bool_col(df, matrix_col)
        overlap_clinvar = int((matrix_mask & bool_col(df, "in_clinvar")).sum())
        rows.append(
            {
                "source": source,
                "source_records": len(sdf),
                "linked_source_records": linked,
                "unresolved_source_records": unresolved,
                "unique_resolved_source_variants": unique_resolved,
                "variant_rows_in_final_matrix": int(matrix_mask.sum()),
                "matrix_rows_also_in_clinvar": overlap_clinvar,
                "matrix_rows_not_in_clinvar": int(matrix_mask.sum() - overlap_clinvar),
                "source_specific_review_stars_nonmissing": int(
                    (~missing(df.loc[matrix_mask, f"{source.lower().replace('boost', 'boost')}_review_stars"]))
                    .sum()
                )
                if source != "CardioBoost"
                else int((~missing(df.loc[matrix_mask, "cardioboost_review_stars"])).sum()),
            }
        )
    rows.insert(
        0,
        {
            "source": "ClinVar",
            "source_records": int(bool_col(df, "in_clinvar").sum()),
            "linked_source_records": int(bool_col(df, "in_clinvar").sum()),
            "unresolved_source_records": 0,
            "unique_resolved_source_variants": int(bool_col(df, "in_clinvar").sum()),
            "variant_rows_in_final_matrix": int(bool_col(df, "in_clinvar").sum()),
            "matrix_rows_also_in_clinvar": int(bool_col(df, "in_clinvar").sum()),
            "matrix_rows_not_in_clinvar": 0,
            "source_specific_review_stars_nonmissing": int((~missing(df.loc[bool_col(df, "in_clinvar"), "clinvar_review_stars"])).sum()),
        },
    )
    return pd.DataFrame(rows)


def feature_coverage_tables(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    groups = [
        ("Genotype / variant identity", ["variant_id", "chrom", "pos", "ref", "alt", "primary_gene"]),
        ("VEP consequence / HGVS", ["vep_worst_consequence", "vep_impact", "vep_hgvsc", "vep_hgvsp"]),
        ("SpliceAI", ["vep_SpliceAI_pred_DS_AG", "vep_SpliceAI_pred_DS_AL", "vep_SpliceAI_pred_DS_DG", "vep_SpliceAI_pred_DS_DL"]),
        ("gnomAD population frequency", ["gnomad_final_af", "gnomad_final_popmax_af", "gnomad_final_homozygote_count"]),
        ("dbNSFP in silico", ["sift_score", "polyphen2_hdiv_score", "revel_score", "metalr_score", "fathmm_xf_coding_score", "cadd_phred"]),
        ("Conservation", ["gerp_rs", "phylop100way_vertebrate", "phastcons100way_vertebrate"]),
        ("AlphaMissense / ESM", ["alphamissense_direct_score", "alphamissense_dbnsfp_score", "esm1b_score"]),
        ("ClinGen", ["clingen_gene_validity_max_classification", "clingen_variant_evidence_classification"]),
        ("Protein domains / AlphaFold", ["uniprot_accession", "protein_position", "alphafold_residue_plddt"]),
        ("DSSP / FreeSASA", ["dssp_secondary_structure_class", "freesasa_relative"]),
        ("FoldX DDG", ["foldx_ddg_kcal_mol"]),
        (
            "HiRO source-record evidence",
            [
                "hiro_source_record_count",
                "hiro_unique_patient_count",
                "hiro_pathogenic_record_count",
                "hiro_benign_record_count",
                "hiro_vus_record_count",
                "hiro_label_discordant_flag",
            ],
        ),
    ]
    group_rows = []
    detail_rows = []
    for group, cols in groups:
        present_cols = [c for c in cols if c in df.columns]
        any_value = pd.Series(False, index=df.index)
        for col in present_cols:
            if group == "HiRO source-record evidence" and "hiro_source_record_features_available" in df.columns:
                if col.endswith("_flag"):
                    col_nonmissing = bool_col(df, "hiro_source_record_features_available") & bool_col(df, col)
                elif col.endswith("_count"):
                    col_nonmissing = pd.to_numeric(df[col], errors="coerce").fillna(0) > 0
                else:
                    col_nonmissing = bool_col(df, "hiro_source_record_features_available") & ~missing(df[col])
            else:
                col_nonmissing = ~missing(df[col])
            any_value |= col_nonmissing
            detail_rows.append(
                {
                    "feature_group": group,
                    "feature": col,
                    "nonmissing_rows": int(col_nonmissing.sum()),
                    "missing_rows": int((~col_nonmissing).sum()),
                    "coverage_percent": pct(col_nonmissing.sum(), len(df)),
                    "meaning": feature_meaning(col),
                }
            )
        group_rows.append(
            {
                "feature_group": group,
                "features_checked": len(present_cols),
                "rows_with_any_value": int(any_value.sum()),
                "coverage_percent": pct(any_value.sum(), len(df)),
            }
        )
    return pd.DataFrame(group_rows), pd.DataFrame(detail_rows)


def feature_meaning(col: str) -> str:
    meanings = {
        "variant_id": "GRCh38 chrom-pos-ref-alt identity; primary join key.",
        "vep_worst_consequence": "Predicted molecular consequence, used for consequence-aware learning.",
        "vep_hgvsc": "Transcript-level HGVS from VEP.",
        "vep_hgvsp": "Protein-level HGVS from VEP.",
        "gnomad_final_af": "Overall population allele frequency; strong benign/prior-probability feature.",
        "gnomad_final_popmax_af": "Maximum population allele frequency; useful for ancestry-aware benign evidence.",
        "gnomad_final_homozygote_count": "Homozygote observation count; helps flag variants inconsistent with severe dominant disease.",
        "sift_score": "Sequence conservation/tolerance predictor; lower values imply deleteriousness.",
        "polyphen2_hdiv_score": "Protein impact predictor trained for Mendelian disease discrimination.",
        "revel_score": "Ensemble missense pathogenicity predictor.",
        "metalr_score": "Ensemble in silico pathogenicity score.",
        "fathmm_xf_coding_score": "Coding variant pathogenicity prediction.",
        "cadd_phred": "Genome-wide deleteriousness score.",
        "gerp_rs": "Evolutionary constraint/conservation.",
        "phylop100way_vertebrate": "PhyloP conservation across vertebrates.",
        "phastcons100way_vertebrate": "PhastCons conservation across vertebrates.",
        "alphamissense_direct_score": "AlphaMissense missense pathogenicity score from direct table join.",
        "alphamissense_dbnsfp_score": "AlphaMissense score via dbNSFP.",
        "esm1b_score": "Protein language model variant-effect score via dbNSFP.",
        "clingen_gene_validity_max_classification": "Strength of ClinGen gene-disease validity evidence.",
        "clingen_variant_evidence_classification": "Variant-level expert/VCEP evidence when matched.",
        "uniprot_accession": "Protein accession used for protein-position features.",
        "protein_position": "Mapped residue position.",
        "alphafold_residue_plddt": "AlphaFold per-residue confidence.",
        "dssp_secondary_structure_class": "Secondary structure class from AlphaFold/DSSP.",
        "freesasa_relative": "Relative solvent accessibility.",
        "foldx_ddg_kcal_mol": "Predicted stability change from FoldX.",
        "hiro_source_record_count": "Number of HiRO patient/source classification records mapped to this variant.",
        "hiro_unique_patient_count": "Number of unique HiRO patients carrying this variant.",
        "hiro_pathogenic_record_count": "Number of HiRO source records classified as pathogenic.",
        "hiro_benign_record_count": "Number of HiRO source records classified as benign.",
        "hiro_vus_record_count": "Number of HiRO source records classified as VUS.",
        "hiro_label_discordant_flag": "Whether HiRO source records mapped to the variant disagree across 3-class labels.",
    }
    return meanings.get(col, "")


def variant_type_table(df: pd.DataFrame) -> pd.DataFrame:
    ref_len = df["ref"].astype(str).str.len()
    alt_len = df["alt"].astype(str).str.len()
    conditions = [
        (ref_len == 1) & (alt_len == 1),
        (ref_len == alt_len) & (ref_len > 1),
        ref_len < alt_len,
        ref_len > alt_len,
    ]
    labels = ["SNV", "MNV", "Insertion/duplication", "Deletion"]
    kind = np.select(conditions, labels, default="Complex/other")
    out = pd.Series(kind).value_counts().rename_axis("variant_type_by_ref_alt").reset_index(name="rows")
    out["percent"] = (out["rows"] / len(df) * 100).round(3)
    return out


def clinvar_confidence_table(df: pd.DataFrame) -> pd.DataFrame:
    clin = df[bool_col(df, "in_clinvar")].copy()
    rows = []
    for label, sub in [("all_clinvar_rows", clin), ("primary_include_clinvar_rows", clin[clin["primary_model_inclusion"] == "include"])]:
        counts = sub["clinvar_review_stars"].fillna("<missing>").astype(str).value_counts()
        for stars, n in counts.items():
            rows.append({"slice": label, "clinvar_review_stars": stars, "rows": int(n), "percent": pct(n, len(sub))})
    return pd.DataFrame(rows)


def hiro_emerge_audit(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "dataset": "HiRO",
                "local_source_identity": "Internal CASPER WES / VERDICT source records with patient phenotype and ACMG-style classification fields.",
                "source_records": 482,
                "variant_rows_in_matrix": int(bool_col(df, "in_hiro").sum()),
                "overlap_with_clinvar_variant_rows": int((bool_col(df, "in_hiro") & bool_col(df, "in_clinvar")).sum()),
                "source_specific_review_stars_nonmissing": int((~missing(df.loc[bool_col(df, "in_hiro"), "hiro_review_stars"])).sum()),
                "interpretation": "Real private HiRO data. Some variants overlap public ClinVar by coordinate, but HiRO labels/phenotypes remain distinct source-record evidence.",
            },
            {
                "dataset": "eMERGE",
                "local_source_identity": "Public eMERGE/Glazer arrhythmia-gene supplement; real cohort variant data from 21,846 participants, genotype-focused in local files.",
                "source_records": 2754,
                "variant_rows_in_matrix": int(bool_col(df, "in_emerge").sum()),
                "overlap_with_clinvar_variant_rows": int((bool_col(df, "in_emerge") & bool_col(df, "in_clinvar")).sum()),
                "source_specific_review_stars_nonmissing": int((~missing(df.loc[bool_col(df, "in_emerge"), "emerge_review_stars"])).sum()),
                "interpretation": "Real eMERGE source variants. Labels may use ClinVar/manual ACMG evidence per manuscript, but local rows are not ClinVar-imported rows.",
            },
            {
                "dataset": "CardioBoost",
                "local_source_identity": "Public CardioBoost coordinate-rich variant dataset reprocessed for this study.",
                "source_records": 355,
                "variant_rows_in_matrix": int(bool_col(df, "in_cardioboost").sum()),
                "overlap_with_clinvar_variant_rows": int((bool_col(df, "in_cardioboost") & bool_col(df, "in_clinvar")).sum()),
                "source_specific_review_stars_nonmissing": int((~missing(df.loc[bool_col(df, "in_cardioboost"), "cardioboost_review_stars"])).sum()),
                "interpretation": "Real public CardioBoost labels/features. Many rows overlap ClinVar coordinates and must be source/leakage controlled.",
            },
        ]
    )


def plot_bar(df: pd.DataFrame, x: str, y: str, title: str, path: Path, hue: str | None = None, rotate: bool = False) -> None:
    plt.figure(figsize=(10, 5.5))
    ax = sns.barplot(data=df, x=x, y=y, hue=hue)
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("Rows")
    if rotate:
        ax.tick_params(axis="x", rotation=35)
        for label in ax.get_xticklabels():
            label.set_horizontalalignment("right")
    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def add_bar_labels(ax) -> None:
    for container in ax.containers:
        ax.bar_label(container, fmt="%.0f", fontsize=8, padding=2)


def make_figures(df: pd.DataFrame, tables: dict[str, pd.DataFrame]) -> None:
    label_long = tables["label_by_source"].melt(
        id_vars=["source"],
        value_vars=["benign", "pathogenic", "vus", "missing_or_conflict"],
        var_name="label",
        value_name="rows",
    )
    plt.figure(figsize=(9, 5.2))
    ax = sns.barplot(data=label_long, x="source", y="rows", hue="label")
    ax.set_title("Label Distribution By Source Flag")
    ax.set_xlabel("")
    ax.set_ylabel("Variant rows")
    plt.tight_layout()
    plt.savefig(FIGURES / "label_distribution_by_source.png")
    plt.close()

    cov = tables["feature_group_coverage"].sort_values("coverage_percent", ascending=True)
    plt.figure(figsize=(9, 6))
    ax = sns.barplot(data=cov, x="coverage_percent", y="feature_group", color="#4C78A8")
    ax.set_title("Feature Coverage By Feature Family")
    ax.set_xlabel("Rows with any value (%)")
    ax.set_ylabel("")
    plt.tight_layout()
    plt.savefig(FIGURES / "feature_group_coverage.png")
    plt.close()

    source_records = tables["source_record_summary"].copy()
    source_records = source_records[source_records["source"].isin(["ClinVar", "HiRO", "eMERGE", "CardioBoost"])]
    source_long = source_records.melt(
        id_vars=["source"],
        value_vars=["source_records", "variant_rows_in_final_matrix"],
        var_name="measure",
        value_name="rows",
    )
    plt.figure(figsize=(10, 5.5))
    ax = sns.barplot(data=source_long, x="source", y="rows", hue="measure")
    ax.set_yscale("log")
    ax.set_title("Source Records Versus Variant-Level Matrix Rows (Log Scale)")
    ax.set_xlabel("")
    ax.set_ylabel("Rows, log scale")
    add_bar_labels(ax)
    plt.tight_layout()
    plt.savefig(FIGURES / "source_records_vs_matrix_rows.png")
    plt.close()

    non_clinvar_source_long = source_long[source_long["source"] != "ClinVar"].copy()
    plt.figure(figsize=(9, 5.2))
    ax = sns.barplot(data=non_clinvar_source_long, x="source", y="rows", hue="measure")
    ax.set_title("HiRO, eMERGE, and CardioBoost Source Rows")
    ax.set_xlabel("")
    ax.set_ylabel("Rows")
    add_bar_labels(ax)
    plt.tight_layout()
    plt.savefig(FIGURES / "source_records_vs_matrix_rows_non_clinvar_zoom.png")
    plt.close()

    top_genes = tables["gene_label_counts"].head(25)
    top_long = top_genes.melt(
        id_vars=["primary_gene"],
        value_vars=["Benign", "Pathogenic", "VUS"],
        var_name="label",
        value_name="rows",
    )
    plot_bar(
        top_long,
        x="primary_gene",
        y="rows",
        hue="label",
        title="Top Genes By Label Distribution",
        path=FIGURES / "top_genes_label_distribution.png",
        rotate=True,
    )

    consequence = tables["vep_consequence_counts"].head(20)
    plt.figure(figsize=(9, 7))
    ax = sns.barplot(data=consequence, y="vep_worst_consequence", x="rows", color="#59A14F")
    ax.set_title("Top VEP Worst Consequences")
    ax.set_xlabel("Rows")
    ax.set_ylabel("")
    plt.tight_layout()
    plt.savefig(FIGURES / "vep_consequence_distribution.png")
    plt.close()

    gnomad = tables["gnomad_status_counts"]
    plot_bar(
        gnomad,
        x="gnomad_final_status",
        y="rows",
        title="gnomAD Final Status",
        path=FIGURES / "gnomad_status_counts.png",
    )

    protein_status = pd.concat(
        [
            value_count_table(df, "protein_feature_status", "status").assign(feature="Protein map"),
            value_count_table(df, "alphafold_plddt_status", "status").assign(feature="AlphaFold pLDDT"),
            value_count_table(df, "dssp_status", "status").assign(feature="DSSP"),
            value_count_table(df, "freesasa_status", "status").assign(feature="FreeSASA"),
            value_count_table(df, "foldx_ddg_status", "status").assign(feature="FoldX"),
        ],
        ignore_index=True,
    )
    protein_ok = protein_status[protein_status["status"].isin(["ok", "not_available", "ref_mismatch", "not_queryable"])]
    plt.figure(figsize=(10, 6))
    ax = sns.barplot(data=protein_ok, x="feature", y="rows", hue="status")
    ax.set_title("Protein / Structure Feature Status")
    ax.set_xlabel("")
    ax.set_ylabel("Rows")
    ax.tick_params(axis="x", rotation=25)
    plt.tight_layout()
    plt.savefig(FIGURES / "protein_structure_status.png")
    plt.close()


def main() -> None:
    ensure_dirs()
    df = read_matrix()

    qc_feature_group = pd.read_csv(QC_DIR / "feature_group_coverage.tsv", sep="\t")
    source_summary = source_record_summary(df)
    label_by_source = label_counts_by_source(df)
    group_cov, detail_cov = feature_coverage_tables(df)
    source_overlap = value_count_table(df, "sources", "source_combination")
    labels = value_count_table(df, "model_label_3class", "model_label_3class")
    variant_types = variant_type_table(df)
    vep_consequence = value_count_table(df, "vep_worst_consequence", "vep_worst_consequence")
    gnomad_status = value_count_table(df, "gnomad_final_status", "gnomad_final_status")
    clinvar_conf = clinvar_confidence_table(df)
    hiro_emerge = hiro_emerge_audit(df)
    gene_labels = pd.read_csv(QC_DIR / "gene_label_counts.tsv", sep="\t")
    sparse_genes = pd.read_csv(QC_DIR / "sparse_pathogenic_genes.tsv", sep="\t")

    inventory = pd.DataFrame(
        [
            {
                "item": "Final modeling matrix",
                "path": str(MATRIX.relative_to(ROOT)),
                "rows": len(df),
                "columns": len(df.columns),
                "size_mb": round(MATRIX.stat().st_size / 1024 / 1024, 2),
                "notes": "One row per unique GRCh38 chrom-pos-ref-alt variant.",
            },
            {
                "item": "QC report",
                "path": str((QC_DIR / "QC_REPORT.md").relative_to(ROOT)),
                "rows": "",
                "columns": "",
                "size_mb": round((QC_DIR / "QC_REPORT.md").stat().st_size / 1024 / 1024, 4),
                "notes": "Full missingness, label, source, duplicate, and leakage QC.",
            },
        ]
    )

    tables = {
        "data_inventory": inventory,
        "source_record_summary": source_summary,
        "source_overlap": source_overlap,
        "label_counts": labels,
        "label_by_source": label_by_source,
        "feature_group_coverage": group_cov,
        "feature_subscore_coverage": detail_cov,
        "qc_feature_group_coverage_status_based": qc_feature_group,
        "variant_type_counts": variant_types,
        "vep_consequence_counts": vep_consequence,
        "gnomad_status_counts": gnomad_status,
        "clinvar_confidence": clinvar_conf,
        "hiro_emerge_cardioboost_source_identity_audit": hiro_emerge,
        "gene_label_counts": gene_labels,
        "sparse_pathogenic_genes": sparse_genes,
    }
    for name, table in tables.items():
        write_table(table, f"{name}.tsv")

    make_figures(df, tables)

    report = []
    report.append("# Data Characteristics Report")
    report.append("")
    report.append("This report summarizes the current final cardiogenetics modeling matrix and its source/feature coverage.")
    report.append("")
    report.append("## Current Matrix")
    report.append("")
    report.append(short_md_table(inventory))
    report.append("")
    report.append("## Source Composition")
    report.append("")
    report.append(short_md_table(source_summary))
    report.append("")
    report.append("## Label Balance")
    report.append("")
    report.append(short_md_table(labels))
    report.append("")
    report.append("## Label Balance By Source Flag")
    report.append("")
    report.append(short_md_table(label_by_source))
    report.append("")
    report.append("## Feature Family Coverage")
    report.append("")
    report.append(short_md_table(group_cov))
    report.append("")
    report.append("## Key Subscore Coverage")
    report.append("")
    report.append(short_md_table(detail_cov, max_rows=80))
    report.append("")
    report.append("## Genotype / Consequence Characteristics")
    report.append("")
    report.append("### Ref/Alt Variant Type")
    report.append(short_md_table(variant_types))
    report.append("")
    report.append("### Top VEP Consequences")
    report.append(short_md_table(vep_consequence, max_rows=25))
    report.append("")
    report.append("## gnomAD Population Coverage")
    report.append("")
    report.append(short_md_table(gnomad_status))
    report.append("")
    report.append("## ClinVar Confidence")
    report.append("")
    report.append(short_md_table(clinvar_conf, max_rows=30))
    report.append("")
    report.append("## HiRO / eMERGE / CardioBoost Source Identity Audit")
    report.append("")
    report.append(short_md_table(hiro_emerge))
    report.append("")
    report.append("Interpretation:")
    report.append("")
    report.append("- HiRO is real internal CASPER WES / VERDICT source data with patient-linked phenotype fields and private source classifications. ClinVar overlap means some exact coordinates are also in ClinVar, not that the HiRO records are ClinVar rows.")
    report.append("- eMERGE is real public cohort/supplement variant data from the Glazer/eMERGE arrhythmia study. The local supplement is genotype-focused; the manuscript states variants were classified using a mix of sequencing-center ACMG/AMP review, ClinVar annotations, manual review, and in vitro evidence for selected VUS. That means its labels can use ClinVar evidence, but the rows are not simply imported ClinVar rows.")
    report.append("- CardioBoost is real public CardioBoost coordinate-rich variant data. It has substantial ClinVar coordinate overlap and therefore needs leakage-aware splitting.")
    report.append("")
    report.append("## Figures")
    report.append("")
    for fig in sorted(FIGURES.glob("*.png")):
        report.append(f"- `figures/{fig.name}`")
    report.append("")
    report.append("## Tables")
    report.append("")
    for table_path in sorted(TABLES.glob("*.tsv")):
        report.append(f"- `tables/{table_path.name}`")
    report.append("")
    report.append("## Short Scientific Takeaway")
    report.append("")
    report.append("The dataset is broad and feature-rich for a cardiogenetics variant model: 86,889 unique variant rows, 43,282 binary supervised rows, 412 columns, and strong coverage for VEP/SpliceAI/gnomAD/ClinGen. In silico missense scores cover about half the full variant matrix because many variants are not dbNSFP-scored missense-like variants. Protein-structure features are intentionally sparse and should be treated as optional high-value features for mapped missense/protein-changing variants rather than universal inputs.")
    (OUT / "DATA_CHARACTERISTICS_REPORT.md").write_text("\n".join(report) + "\n")

    summary = {
        "matrix_rows": int(len(df)),
        "matrix_columns": int(len(df.columns)),
        "binary_supervised_rows": int(df["model_label_3class"].isin(["Benign", "Pathogenic"]).sum()),
        "figures": sorted([p.name for p in FIGURES.glob("*.png")]),
        "tables": sorted([p.name for p in TABLES.glob("*.tsv")]),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(OUT)


if __name__ == "__main__":
    main()
