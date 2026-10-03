#!/usr/bin/env python3
"""Plot ROC curves for corrected CatBoost model outputs."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import auc, roc_curve


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/model_performance/source_heldout_hiro_emerge_rescued/roc_curves"

PREDICTION_FILES = {
    "internal_grouped": ROOT
    / "results/models/primary_binary_catboost_v0_strict_internal_hiro_emerge_rescued/primary_binary_predictions_split_internal_grouped.tsv",
    "source_heldout": ROOT
    / "results/models/primary_binary_catboost_v0_strict_source_heldout_hiro_emerge_rescued/primary_binary_predictions_split_source_heldout.tsv",
    "gene_stress": ROOT
    / "results/models/primary_binary_catboost_v0_strict_gene_stress_hiro_emerge_rescued/primary_binary_predictions_split_gene_stress.tsv",
}

SPLIT_COLUMNS = {
    "internal_grouped": "split_internal_grouped",
    "source_heldout": "split_source_heldout",
    "gene_stress": "split_gene_stress",
}

PLOT_ORDER = {
    "internal_grouped": ["validation", "test"],
    "source_heldout": ["validation", "external_hiro", "external_emerge", "external_cardioboost", "external_all"],
    "gene_stress": ["validation", "well_represented_test", "sparse_gene_stress_test"],
}


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def roc_for_subset(df: pd.DataFrame) -> dict:
    fpr, tpr, thresholds = roc_curve(df["y_true"].astype(int), df["pathogenic_probability"].astype(float))
    return {
        "fpr": fpr,
        "tpr": tpr,
        "thresholds": thresholds,
        "auroc": float(auc(fpr, tpr)),
        "rows": int(len(df)),
        "positives_pathogenic": int(df["y_true"].astype(int).sum()),
        "negatives_benign": int((df["y_true"].astype(int) == 0).sum()),
    }


def plot_group(model_name: str, curves: dict[str, dict], out_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    for split_name, curve in curves.items():
        ax.plot(
            curve["fpr"],
            curve["tpr"],
            linewidth=2.2,
            label=f"{split_name} (AUC={curve['auroc']:.3f}, n={curve['rows']})",
        )
    ax.plot([0, 1], [0, 1], linestyle="--", linewidth=1.2, color="0.45", label="random")
    ax.set_title(f"ROC Curves: {model_name.replace('_', ' ').title()}")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate / P/LP sensitivity")
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)
    ax.grid(True, linestyle=":", linewidth=0.8, alpha=0.7)
    ax.legend(loc="lower right", frameon=True)
    fig.tight_layout()
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    manifest = {"outputs": {}, "inputs": {k: rel(v) for k, v in PREDICTION_FILES.items()}}

    all_external_curves = {}
    for model_name, pred_path in PREDICTION_FILES.items():
        df = pd.read_csv(pred_path, sep="\t", low_memory=False)
        split_col = SPLIT_COLUMNS[model_name]
        curves = {}
        for split_name in PLOT_ORDER[model_name]:
            if split_name == "external_all":
                sub = df[df[split_col].astype(str).str.startswith("external_")].copy()
            else:
                sub = df[df[split_col].eq(split_name)].copy()
            if sub.empty:
                continue
            curve = roc_for_subset(sub)
            curves[split_name] = curve
            summary_rows.append(
                {
                    "model": model_name,
                    "split": split_name,
                    "rows": curve["rows"],
                    "positives_pathogenic": curve["positives_pathogenic"],
                    "negatives_benign": curve["negatives_benign"],
                    "auroc": curve["auroc"],
                }
            )
            if model_name == "source_heldout" and split_name.startswith("external"):
                all_external_curves[split_name] = curve

        out_png = OUT / f"roc_{model_name}.png"
        plot_group(model_name, curves, out_png)
        manifest["outputs"][f"roc_{model_name}"] = rel(out_png)

    if all_external_curves:
        out_png = OUT / "roc_external_validation_sets.png"
        plot_group("external_validation_sets", all_external_curves, out_png)
        manifest["outputs"]["roc_external_validation_sets"] = rel(out_png)

    summary = pd.DataFrame(summary_rows)
    summary_path = OUT / "roc_auc_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    manifest["outputs"]["summary"] = rel(summary_path)
    (OUT / "roc_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
