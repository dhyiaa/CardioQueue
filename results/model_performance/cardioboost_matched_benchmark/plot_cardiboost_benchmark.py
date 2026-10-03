#!/usr/bin/env python3
"""Create CardioBoost matched benchmark ROC and PR figures."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import auc, precision_recall_curve, roc_curve


BASE = Path(__file__).resolve().parent
TABLE_DIR = BASE / "tables"
FIG_DIR = BASE / "figures"

BENCHMARKS = {
    "hgvs_cdot_primary_external_all": (
        "External matched CardioBoost-eligible rows",
        lambda df: df["split_source_heldout"].astype(str).str.startswith("external_"),
    ),
    "hgvs_cdot_primary_external_cardioboost": (
        "External CardioBoost-source matched rows",
        lambda df: df["split_source_heldout"].astype(str).eq("external_cardioboost"),
    ),
}

SCORE_COLUMNS = {
    "our_catboost_score": "CardioQueue",
    "cardioboost_score": "Official CardioBoost",
    "revel_score_numeric": "REVEL",
    "alphamissense_best_score": "AlphaMissense",
    "cadd_phred": "CADD PHRED",
}

FIG_DPI = 600
FIG_SIZE = (8.0, 6.0)
LINE_WIDTH = 2.4
FONT_SIZE = 11

plt.rcParams.update(
    {
        "font.size": FONT_SIZE,
        "axes.titlesize": 13,
        "axes.labelsize": 12,
        "legend.fontsize": 8.5,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "savefig.dpi": FIG_DPI,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def _load_metrics() -> pd.DataFrame:
    return pd.read_csv(TABLE_DIR / "hgvs_cardioBoost_matched_metrics.tsv", sep="\t")


def _load_benchmark() -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "hgvs_cdot_cardioBoost_matched_benchmark.tsv", sep="\t")
    if "y_true" not in df.columns:
        raise ValueError("Expected y_true column in benchmark table.")
    df = df[df["y_true"].isin([0, 1])].copy()
    df["y_true"] = df["y_true"].astype(int)
    return df


def _plot_roc(df: pd.DataFrame, metrics: pd.DataFrame, benchmark: str, title: str, mask_fn) -> None:
    subset = df[mask_fn(df)].copy()
    metric_subset = metrics[metrics["benchmark"] == benchmark].copy()

    fig, ax = plt.subplots(figsize=FIG_SIZE, dpi=FIG_DPI)
    for score_col, label in SCORE_COLUMNS.items():
        if score_col not in subset.columns:
            continue
        curve_df = subset[["y_true", score_col]].dropna()
        if curve_df["y_true"].nunique() < 2 or curve_df.empty:
            continue
        fpr, tpr, _ = roc_curve(curve_df["y_true"], curve_df[score_col])
        auroc_row = metric_subset.loc[metric_subset["model"] == label, "AUROC"]
        suffix = f" (AUROC {auroc_row.iloc[0]:.3f})" if not auroc_row.empty else ""
        ax.plot(fpr, tpr, linewidth=LINE_WIDTH, label=f"{label}{suffix}")

    ax.plot([0, 1], [0, 1], linestyle=":", color="#777777", linewidth=1.2)
    ax.set_title(title)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_xticks([i / 10 for i in range(0, 11)])
    ax.set_yticks([i / 10 for i in range(0, 11)])
    ax.minorticks_on()
    ax.grid(which="major", alpha=0.28, linewidth=0.8)
    ax.grid(which="minor", alpha=0.12, linewidth=0.5)
    ax.legend(loc="lower right", frameon=True, framealpha=0.94)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(FIG_DIR / f"{benchmark}_roc.{ext}", dpi=FIG_DPI)
    plt.close(fig)


def _plot_pr(df: pd.DataFrame, metrics: pd.DataFrame, benchmark: str, title: str, mask_fn) -> None:
    subset = df[mask_fn(df)].copy()
    metric_subset = metrics[metrics["benchmark"] == benchmark].copy()

    fig, ax = plt.subplots(figsize=FIG_SIZE, dpi=FIG_DPI)
    prevalence = subset["y_true"].mean()
    ax.axhline(prevalence, linestyle=":", color="#777777", linewidth=1.2, label=f"Prevalence {prevalence:.3f}")

    for score_col, label in SCORE_COLUMNS.items():
        if score_col not in subset.columns:
            continue
        curve_df = subset[["y_true", score_col]].dropna()
        if curve_df["y_true"].nunique() < 2 or curve_df.empty:
            continue
        precision, recall, _ = precision_recall_curve(curve_df["y_true"], curve_df[score_col])
        auprc_value = auc(recall, precision)
        auprc_row = metric_subset.loc[metric_subset["model"] == label, "AUPRC"]
        if not auprc_row.empty:
            auprc_value = auprc_row.iloc[0]
        ax.plot(recall, precision, linewidth=LINE_WIDTH, label=f"{label} (AUPRC {auprc_value:.3f})")

    ax.set_title(title)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_xticks([i / 10 for i in range(0, 11)])
    ax.set_yticks([i / 10 for i in range(0, 11)])
    ax.minorticks_on()
    ax.grid(which="major", alpha=0.28, linewidth=0.8)
    ax.grid(which="minor", alpha=0.12, linewidth=0.5)
    ax.legend(loc="lower left", frameon=True, framealpha=0.94)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(FIG_DIR / f"{benchmark}_precision_recall.{ext}", dpi=FIG_DPI)
    plt.close(fig)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    metrics = _load_metrics()
    benchmark_df = _load_benchmark()
    for benchmark, (title, mask_fn) in BENCHMARKS.items():
        _plot_roc(benchmark_df, metrics, benchmark, title, mask_fn)
        _plot_pr(benchmark_df, metrics, benchmark, title, mask_fn)


if __name__ == "__main__":
    main()
