#!/usr/bin/env python3
"""Generate manuscript-facing figures and tables for the primary binary model."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.calibration import calibration_curve
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    confusion_matrix,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/manuscript_figures_tables"
FIG = OUT / "figures"
TAB = OUT / "tables"
MATRIX = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
PRED = ROOT / "results/models/cardioqueue_v1_source_heldout/primary_binary_predictions_split_source_heldout.tsv"
TRIAGE = ROOT / "results/model_performance/primary_binary_three_zone_all_sources"


SOURCE_LABELS = {
    "validation": "Validation",
    "external_all": "External all",
    "external_hiro": "HiRO",
    "external_emerge": "eMERGE",
    "external_cardioboost": "CardioBoost",
}

GENE_GROUPS = {
    "Arrhythmia": {
        "AKAP9", "ANK2", "CACNA1C", "CACNB2", "CALM1", "CALM2", "CALM3", "CASQ2", "HCN4",
        "KCNE1", "KCNE1B", "KCNE2", "KCNH2", "KCNJ2", "KCNQ1", "RYR2", "SCN1B", "SCN5A",
        "SNTA1", "TRDN", "TRPM4",
    },
    "Cardiomyopathy": {
        "ACTC1", "BAG3", "CSRP3", "DES", "DSC2", "DSG2", "DSP", "EMD", "FHOD3", "FLNC",
        "JUP", "LMNA", "MIB1", "MYBPC3", "MYH7", "PKP2", "PLN", "RBM20", "TCAP", "TMEM43",
        "TNNI3", "TNNI3K", "TNNT2", "TPM1", "TTN",
    },
    "RASopathy_or_phenocopy": {"HRAS", "KRAS", "MAP2K2", "NRAS", "PTPN11", "RAF1", "RIT1", "SOS1"},
    "Metabolic_or_systemic": {"GLA", "PRKAG2", "SLC22A5", "DMD"},
}


def ensure_dirs() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)


def read_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    matrix_cols = [
        "variant_id", "primary_gene", "model_label_3class", "training_slice_primary_binary",
        "split_source_heldout", "source_leakage_group", "sources", "source_count", "source_row_count",
        "in_clinvar", "in_hiro", "in_emerge", "in_cardioboost",
        "dbnsfp_status", "dbnsfp_REVEL_score", "dbnsfp_CADD_phred", "dbnsfp_AlphaMissense_score",
        "alphamissense_direct_status", "alphamissense_direct_score",
        "gnomad_final_status", "gnomad_final_source", "gnomad_final_af", "gnomad_observed_flag",
        "gnomad_confirmed_absent_flag", "gnomad_not_joined_flag",
        "vep_annotation_status", "vep_worst_consequence", "vep_consequence", "vep_impact",
        "vep_any_missense", "vep_any_splice_region", "vep_any_frameshift", "vep_any_stop_gained",
        "foldx_ddg_status", "foldx_ddg_kcal_mol",
        "spliceai_status", "vep_SpliceAI_pred_DS_AG", "vep_SpliceAI_pred_DS_AL",
        "vep_SpliceAI_pred_DS_DG", "vep_SpliceAI_pred_DS_DL",
    ]
    matrix = pd.read_csv(MATRIX, sep="\t", usecols=lambda c: c in matrix_cols, low_memory=False)
    pred = pd.read_csv(PRED, sep="\t", low_memory=False)
    joined = pred.merge(matrix.drop(columns=["primary_gene", "model_label_3class", "split_source_heldout"], errors="ignore"), on="variant_id", how="left")
    return matrix, pred, joined


def subset_frame(df: pd.DataFrame, subset: str) -> pd.DataFrame:
    if subset == "external_all":
        return df[df["split_source_heldout"].isin(["external_hiro", "external_emerge", "external_cardioboost"])].copy()
    return df[df["split_source_heldout"].eq(subset)].copy()


def binary_metrics(y: np.ndarray, p: np.ndarray) -> dict[str, float | int]:
    out: dict[str, float | int] = {
        "rows": int(len(y)),
        "pathogenic": int(y.sum()),
        "benign": int((y == 0).sum()),
    }
    if len(y) == 0 or len(np.unique(y)) < 2:
        out.update({"auroc": np.nan, "auprc": np.nan, "brier": np.nan, "sensitivity": np.nan, "specificity": np.nan, "ppv": np.nan, "npv": np.nan})
        return out
    pred = (p >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    out.update(
        {
            "auroc": float(roc_auc_score(y, p)),
            "auprc": float(average_precision_score(y, p)),
            "brier": float(brier_score_loss(y, p)),
            "tp": int(tp),
            "fp": int(fp),
            "tn": int(tn),
            "fn": int(fn),
            "sensitivity": float(tp / (tp + fn)) if (tp + fn) else np.nan,
            "specificity": float(tn / (tn + fp)) if (tn + fp) else np.nan,
            "ppv": float(tp / (tp + fp)) if (tp + fp) else np.nan,
            "npv": float(tn / (tn + fn)) if (tn + fn) else np.nan,
        }
    )
    zones = np.full(len(p), "defer", dtype=object)
    zones[p <= 0.1] = "benign_like"
    zones[p >= 0.9] = "pathogenic_like"
    out["deferral_rate_0_1_0_9"] = float((zones == "defer").mean())
    called = zones != "defer"
    out["high_conf_accuracy_0_1_0_9"] = float(((zones[called] == "pathogenic_like") == (y[called] == 1)).mean()) if called.any() else np.nan
    return out


def plot_roc_pr(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    subsets = ["validation", "external_all", "external_hiro", "external_emerge", "external_cardioboost"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    for subset in subsets:
        sub = subset_frame(df, subset)
        y = sub["y_true"].astype(int).to_numpy()
        p = sub["pathogenic_probability"].astype(float).to_numpy()
        if len(np.unique(y)) < 2:
            continue
        fpr, tpr, _ = roc_curve(y, p)
        precision, recall, _ = precision_recall_curve(y, p)
        auroc = roc_auc_score(y, p)
        auprc = average_precision_score(y, p)
        axes[0].plot(fpr, tpr, lw=2, label=f"{SOURCE_LABELS[subset]} ({auroc:.3f})")
        axes[1].plot(recall, precision, lw=2, label=f"{SOURCE_LABELS[subset]} ({auprc:.3f})")
        m = binary_metrics(y, p)
        m.update({"subset": subset})
        rows.append(m)
    axes[0].plot([0, 1], [0, 1], color="0.75", lw=1, linestyle="--")
    axes[0].set_title("ROC")
    axes[0].set_xlabel("False positive rate")
    axes[0].set_ylabel("True positive rate")
    axes[1].set_title("Precision-recall")
    axes[1].set_xlabel("Recall")
    axes[1].set_ylabel("Precision")
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "primary_binary_roc_pr_panel.png", dpi=240)
    plt.close(fig)
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "external_validation_performance.tsv", sep="\t", index=False)
    return out


def calibration_analysis(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    curve_rows = []
    subsets = ["validation", "external_all", "external_hiro", "external_emerge", "external_cardioboost"]
    fig, ax = plt.subplots(figsize=(6.8, 5.8))
    ax.plot([0, 1], [0, 1], color="0.65", linestyle="--", lw=1, label="Perfect calibration")
    for subset in subsets:
        sub = subset_frame(df, subset)
        y = sub["y_true"].astype(int).to_numpy()
        p = sub["pathogenic_probability"].astype(float).to_numpy()
        if len(np.unique(y)) < 2:
            continue
        n_bins = 10 if len(sub) >= 500 else 5
        frac_pos, mean_pred = calibration_curve(y, p, n_bins=n_bins, strategy="quantile")
        ax.plot(mean_pred, frac_pos, marker="o", lw=2, label=SOURCE_LABELS[subset])
        rows.append(
            {
                "subset": subset,
                "rows": int(len(sub)),
                "brier": float(brier_score_loss(y, p)),
                "mean_predicted_probability": float(np.mean(p)),
                "observed_pathogenic_fraction": float(np.mean(y)),
                "calibration_intercept_like_mean_error": float(np.mean(p) - np.mean(y)),
            }
        )
        for i, (mp, fp) in enumerate(zip(mean_pred, frac_pos), start=1):
            curve_rows.append({"subset": subset, "bin": i, "mean_predicted_probability": mp, "observed_pathogenic_fraction": fp})
    ax.set_title("Reliability plot")
    ax.set_xlabel("Mean predicted pathogenic probability")
    ax.set_ylabel("Observed pathogenic fraction")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "primary_binary_calibration_reliability.png", dpi=240)
    plt.close(fig)
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "calibration_brier_summary.tsv", sep="\t", index=False)
    pd.DataFrame(curve_rows).to_csv(TAB / "calibration_curve_points.tsv", sep="\t", index=False)
    return out


def source_flow_and_coverage(matrix: pd.DataFrame, pred: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    all_rows = len(matrix)
    binary = matrix[matrix["training_slice_primary_binary"].astype(str).str.lower().eq("true")]
    vus = matrix[matrix["model_label_3class"].eq("VUS")]
    flow_rows = [
        {"stage": "Final matrix rows", "rows": all_rows},
        {"stage": "Primary binary trainable rows", "rows": len(binary)},
        {"stage": "VUS scoring pool rows", "rows": len(vus)},
        {"stage": "Source-held-out train rows", "rows": int(pred["split_source_heldout"].eq("train").sum())},
        {"stage": "Source-held-out validation rows", "rows": int(pred["split_source_heldout"].eq("validation").sum())},
        {"stage": "External HiRO binary rows", "rows": int(pred["split_source_heldout"].eq("external_hiro").sum())},
        {"stage": "External eMERGE binary rows", "rows": int(pred["split_source_heldout"].eq("external_emerge").sum())},
        {"stage": "External CardioBoost binary rows", "rows": int(pred["split_source_heldout"].eq("external_cardioboost").sum())},
    ]
    flow = pd.DataFrame(flow_rows)
    flow.to_csv(TAB / "source_flow_modeling_rows.tsv", sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(10, 5.2))
    ax.barh(flow["stage"], flow["rows"], color="#3576A8")
    ax.invert_yaxis()
    ax.set_xlabel("Rows")
    ax.set_title("Modeling row flow")
    for i, v in enumerate(flow["rows"]):
        ax.text(v, i, f" {v:,}", va="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "source_flow_modeling_rows.png", dpi=240)
    plt.close(fig)

    coverage_defs = {
        "dbNSFP": matrix["dbnsfp_status"].astype(str).str.lower().eq("ok"),
        "gnomAD final AF": pd.to_numeric(matrix["gnomad_final_af"], errors="coerce").notna(),
        "gnomAD observed": matrix["gnomad_observed_flag"].astype(str).str.lower().eq("true"),
        "VEP": matrix["vep_annotation_status"].astype(str).str.lower().eq("ok"),
        "SpliceAI": matrix[["vep_SpliceAI_pred_DS_AG", "vep_SpliceAI_pred_DS_AL", "vep_SpliceAI_pred_DS_DG", "vep_SpliceAI_pred_DS_DL"]].apply(pd.to_numeric, errors="coerce").notna().any(axis=1),
        "AlphaMissense direct": matrix["alphamissense_direct_status"].astype(str).str.lower().eq("ok"),
        "FoldX DDG": pd.to_numeric(matrix["foldx_ddg_kcal_mol"], errors="coerce").notna(),
    }
    coverage = pd.DataFrame(
        [
            {
                "feature_group": k,
                "covered_rows": int(v.sum()),
                "total_rows": int(len(matrix)),
                "coverage_fraction": float(v.mean()),
            }
            for k, v in coverage_defs.items()
        ]
    ).sort_values("coverage_fraction", ascending=True)
    coverage.to_csv(TAB / "feature_group_coverage_manuscript.tsv", sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    ax.barh(coverage["feature_group"], coverage["coverage_fraction"], color="#2F8F83")
    ax.set_xlim(0, 1)
    ax.set_xlabel("Coverage across final matrix")
    ax.set_title("Feature coverage")
    for i, r in coverage.reset_index(drop=True).iterrows():
        ax.text(r["coverage_fraction"], i, f" {r['covered_rows']:,}", va="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "feature_group_coverage_manuscript.png", dpi=240)
    plt.close(fig)
    return flow, coverage


def simplify_consequence(value: object) -> str:
    text = str(value).lower()
    if text in {"", "nan", "none"}:
        return "missing"
    if "missense" in text:
        return "missense"
    if "frameshift" in text:
        return "frameshift"
    if "stop_gained" in text or "nonsense" in text:
        return "stop_gained"
    if "splice" in text:
        return "splice"
    if "synonymous" in text:
        return "synonymous"
    if "inframe" in text:
        return "inframe_indel"
    if "utr" in text or "intron" in text:
        return "noncoding_or_intronic"
    return text.split("&")[0].split(",")[0][:40]


def assign_gene_group(gene: object) -> str:
    g = str(gene)
    for group, genes in GENE_GROUPS.items():
        if g in genes:
            return group
    return "Other_or_panel_artifact"


def stratified_performance(df: pd.DataFrame) -> pd.DataFrame:
    work = df.copy()
    path_counts = work[work["model_label_3class"].eq("Pathogenic")].groupby("primary_gene").size()
    work["gene_group"] = work["primary_gene"].map(assign_gene_group)
    work["sparse_gene_group"] = work["primary_gene"].map(lambda g: "sparse_PLP_lt30" if path_counts.get(g, 0) < 30 else "well_represented_PLP_ge30")
    work["consequence_group"] = work["vep_worst_consequence"].map(simplify_consequence)
    work["dbnsfp_coverage_group"] = np.where(work["dbnsfp_status"].astype(str).str.lower().eq("ok"), "dbNSFP_matched", "dbNSFP_missing")
    work["gnomad_status_group"] = np.select(
        [
            work["gnomad_observed_flag"].astype(str).str.lower().eq("true"),
            work["gnomad_confirmed_absent_flag"].astype(str).str.lower().eq("true"),
            work["gnomad_not_joined_flag"].astype(str).str.lower().eq("true"),
        ],
        ["gnomAD_observed", "gnomAD_confirmed_absent", "gnomAD_not_joined"],
        default=work["gnomad_final_status"].fillna("gnomAD_unknown").astype(str),
    )
    strata = {
        "gene_group": "gene_group",
        "sparse_gene_group": "sparse_gene_group",
        "consequence_group": "consequence_group",
        "dbnsfp_coverage_group": "dbnsfp_coverage_group",
        "gnomad_status_group": "gnomad_status_group",
    }
    rows = []
    subsets = ["validation", "external_all", "external_hiro", "external_emerge", "external_cardioboost"]
    for subset in subsets:
        base = subset_frame(work, subset)
        for stratum_type, col in strata.items():
            for stratum, sub in base.groupby(col, dropna=False):
                y = sub["y_true"].astype(int).to_numpy()
                p = sub["pathogenic_probability"].astype(float).to_numpy()
                m = binary_metrics(y, p)
                m.update(
                    {
                        "subset": subset,
                        "stratum_type": stratum_type,
                        "stratum": str(stratum),
                        "metric_status": "ok" if len(sub) >= 20 and len(np.unique(y)) == 2 else "small_or_single_class",
                    }
                )
                rows.append(m)
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "stratified_performance.tsv", sep="\t", index=False)

    heat = out[(out["subset"].eq("external_all")) & (out["metric_status"].eq("ok"))].copy()
    heat["label"] = heat["stratum_type"].str.replace("_group", "", regex=False) + ": " + heat["stratum"]
    heat = heat.sort_values(["stratum_type", "rows"], ascending=[True, False]).head(30)
    fig, ax = plt.subplots(figsize=(8.5, max(5.5, len(heat) * 0.28)))
    vals = heat[["auroc", "auprc", "sensitivity", "specificity", "ppv"]].astype(float).to_numpy()
    im = ax.imshow(vals, aspect="auto", vmin=0.5, vmax=1.0, cmap="viridis")
    ax.set_yticks(np.arange(len(heat)))
    ax.set_yticklabels(heat["label"], fontsize=8)
    ax.set_xticks(np.arange(5))
    ax.set_xticklabels(["AUROC", "AUPRC", "Sens", "Spec", "PPV"])
    ax.set_title("External-all stratified performance")
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            ax.text(j, i, f"{vals[i, j]:.2f}", ha="center", va="center", color="white" if vals[i, j] < 0.75 else "black", fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    fig.tight_layout()
    fig.savefig(FIG / "external_all_stratified_performance_heatmap.png", dpi=240)
    plt.close(fig)
    return out


def all_label_triage_analysis() -> pd.DataFrame:
    summary = pd.read_csv(TRIAGE / "primary_binary_three_zone_summary_metrics.tsv", sep="\t")
    wanted = [
        "dataset", "rows_total", "rows_with_prediction", "rows_without_prediction",
        "true_benign", "true_vus", "true_pathogenic", "pred_benign", "pred_vus_defer", "pred_pathogenic",
        "exact_3class_accuracy", "weighted_kappa", "macro_f1", "plp_sensitivity", "specificity_non_plp",
        "ppv_plp", "npv_non_plp", "binary_mcc_plp_vs_rest", "vus_deferral_rate",
        "true_positives", "false_positives", "plp_to_vus_deferrals", "plp_to_benign_errors",
        "benign_to_plp_escalations", "vus_to_plp_escalations",
    ]
    summary = summary[wanted].copy()
    summary.to_csv(TAB / "primary_binary_three_zone_all_labels_summary.tsv", sep="\t", index=False)

    confusion_files = sorted(TRIAGE.glob("*_confusion.tsv"))
    combined = []
    for path in confusion_files:
        name = path.name.replace("_confusion.tsv", "")
        mat = pd.read_csv(path, sep="\t", index_col=0)
        mat.to_csv(TAB / f"triage_confusion_{name}.tsv", sep="\t")
        long = mat.reset_index(names="true_label").melt(id_vars="true_label", var_name="predicted_zone", value_name="rows")
        long["dataset_key"] = name
        combined.append(long)
    if combined:
        pd.concat(combined, ignore_index=True).to_csv(TAB / "primary_binary_three_zone_all_labels_confusions_long.tsv", sep="\t", index=False)

    vus = summary[summary["true_vus"].fillna(0).gt(0)].copy()
    if not vus.empty:
        vus["vus_to_benign_like"] = np.nan
        vus["vus_to_defer"] = np.nan
        vus["vus_to_pathogenic_like"] = np.nan
        for i, row in vus.iterrows():
            key = row["dataset"].lower().replace(" ", "_")
            candidates = list(TRIAGE.glob(f"{key}_confusion.tsv"))
            if not candidates:
                continue
            mat = pd.read_csv(candidates[0], sep="\t", index_col=0)
            if "true_VUS" in mat.index:
                total = mat.loc["true_VUS"].sum()
                if total:
                    vus.at[i, "vus_to_benign_like"] = mat.loc["true_VUS", "pred_Benign"] / total
                    vus.at[i, "vus_to_defer"] = mat.loc["true_VUS", "pred_VUS"] / total
                    vus.at[i, "vus_to_pathogenic_like"] = mat.loc["true_VUS", "pred_Pathogenic"] / total
        plot = vus.dropna(subset=["vus_to_benign_like", "vus_to_defer", "vus_to_pathogenic_like"])
        if not plot.empty:
            fig, ax = plt.subplots(figsize=(10, 5.2))
            labels = plot["dataset"].tolist()
            benign = plot["vus_to_benign_like"].to_numpy()
            defer = plot["vus_to_defer"].to_numpy()
            path = plot["vus_to_pathogenic_like"].to_numpy()
            y = np.arange(len(labels))
            ax.barh(y, benign, label="Benign-like", color="#4C78A8")
            ax.barh(y, defer, left=benign, label="Deferred/VUS zone", color="#F2CF5B")
            ax.barh(y, path, left=benign + defer, label="Pathogenic-like", color="#C44E52")
            ax.set_yticks(y)
            ax.set_yticklabels(labels)
            ax.set_xlim(0, 1)
            ax.set_xlabel("Fraction of true VUS rows")
            ax.set_title("How the binary model triages VUS rows")
            ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3)
            fig.tight_layout()
            fig.savefig(FIG / "legacy_vus_three_zone_triage_distribution.png", dpi=240)
            plt.close(fig)

    return summary


def write_report(perf: pd.DataFrame, cal: pd.DataFrame, flow: pd.DataFrame, coverage: pd.DataFrame, strat: pd.DataFrame, triage: pd.DataFrame) -> None:
    def md_table(df: pd.DataFrame, cols: list[str], n: int | None = None) -> list[str]:
        d = df[cols].copy()
        if n is not None:
            d = d.head(n)
        for c in d.columns:
            if pd.api.types.is_float_dtype(d[c]):
                d[c] = d[c].map(lambda x: "" if pd.isna(x) else f"{x:.3f}")
        lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
        for _, r in d.iterrows():
            lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
        return lines

    lines = [
        "# Manuscript Figures, Calibration, and Stratified Performance",
        "",
        "This bundle summarizes the primary binary CatBoost P/LP vs B/LB model with 0.1/0.9 deferral. The figures and tables are generated from the source-held-out model predictions and the rescue-aware ready matrix.",
        "",
        "## External Validation",
        "",
    ]
    lines += md_table(perf, ["subset", "rows", "pathogenic", "benign", "auroc", "auprc", "brier", "sensitivity", "specificity", "ppv", "deferral_rate_0_1_0_9", "high_conf_accuracy_0_1_0_9"])
    lines += [
        "",
        "## All-Label Clinical Triage",
        "",
        "This table applies the same binary model to Benign, VUS, and Pathogenic rows using three probability zones: `<=0.1` benign-like, `0.1-0.9` deferred/VUS-zone, and `>=0.9` pathogenic-like. This is not a true three-class classifier. It is a clinical triage analysis.",
        "",
    ]
    triage_cols = [
        "dataset", "rows_with_prediction", "true_benign", "true_vus", "true_pathogenic",
        "pred_benign", "pred_vus_defer", "pred_pathogenic", "exact_3class_accuracy",
        "plp_sensitivity", "specificity_non_plp", "ppv_plp", "vus_deferral_rate",
        "plp_to_benign_errors", "benign_to_plp_escalations", "vus_to_plp_escalations",
    ]
    lines += md_table(triage, triage_cols)
    lines += [
        "",
        "## Calibration",
        "",
    ]
    lines += md_table(cal, ["subset", "rows", "brier", "mean_predicted_probability", "observed_pathogenic_fraction", "calibration_intercept_like_mean_error"])
    lines += [
        "",
        "## Feature Coverage",
        "",
    ]
    lines += md_table(coverage.sort_values("coverage_fraction", ascending=False), ["feature_group", "covered_rows", "total_rows", "coverage_fraction"])
    lines += [
        "",
        "## Key Stratified External-All Results",
        "",
    ]
    key = strat[(strat["subset"].eq("external_all")) & (strat["metric_status"].eq("ok"))].copy()
    key = key.sort_values(["stratum_type", "rows"], ascending=[True, False])
    lines += md_table(key, ["stratum_type", "stratum", "rows", "pathogenic", "benign", "auroc", "auprc", "sensitivity", "specificity", "ppv"], n=30)
    lines += [
        "",
        "## Generated Files",
        "",
        "- `figures/primary_binary_roc_pr_panel.png`",
        "- `figures/primary_binary_calibration_reliability.png`",
        "- `figures/source_flow_modeling_rows.png`",
        "- `figures/feature_group_coverage_manuscript.png`",
        "- `figures/external_all_stratified_performance_heatmap.png`",
        "- `figures/legacy_vus_three_zone_triage_distribution.png`",
        "- `tables/external_validation_performance.tsv`",
        "- `tables/primary_binary_three_zone_all_labels_summary.tsv`",
        "- `tables/primary_binary_three_zone_all_labels_confusions_long.tsv`",
        "- `tables/calibration_brier_summary.tsv`",
        "- `tables/calibration_curve_points.tsv`",
        "- `tables/source_flow_modeling_rows.tsv`",
        "- `tables/feature_group_coverage_manuscript.tsv`",
        "- `tables/stratified_performance.tsv`",
    ]
    (OUT / "MANUSCRIPT_RESULTS_BUNDLE_REPORT.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    ensure_dirs()
    matrix, pred, joined = read_inputs()
    perf = plot_roc_pr(joined)
    cal = calibration_analysis(joined)
    flow, coverage = source_flow_and_coverage(matrix, pred)
    strat = stratified_performance(joined)
    triage = all_label_triage_analysis()
    write_report(perf, cal, flow, coverage, strat, triage)
    manifest = {
        "matrix": str(MATRIX.relative_to(ROOT)),
        "predictions": str(PRED.relative_to(ROOT)),
        "output": str(OUT.relative_to(ROOT)),
        "figures": sorted(p.name for p in FIG.glob("*.png")),
        "tables": sorted(p.name for p in TAB.glob("*.tsv")),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
