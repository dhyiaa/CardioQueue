#!/usr/bin/env python3
"""Evaluate corrected source-held-out model predictions.

This produces both strict binary metrics and a clinical-style three-call view
where middle probabilities are deferred as VUS.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    cohen_kappa_score,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
)


ROOT = Path(__file__).resolve().parents[3]
PREDICTIONS = ROOT / (
    "results/models/primary_binary_catboost_v0_strict_source_heldout_hiro_emerge_rescued/"
    "primary_binary_predictions_split_source_heldout.tsv"
)
OUT_DIR = ROOT / "results/model_performance/source_heldout_hiro_emerge_rescued"


LABEL_ORDER = ["Benign", "VUS", "Pathogenic"]
LABEL_TO_NUM = {"Benign": 0, "VUS": 1, "Pathogenic": 2}


def assign_calls(proba: pd.Series, lower: float, upper: float) -> pd.Series:
    if lower == upper:
        return pd.Series(np.where(proba >= upper, "Pathogenic", "Benign"), index=proba.index)
    return pd.Series(
        np.select([proba < lower, proba >= upper], ["Benign", "Pathogenic"], default="VUS"),
        index=proba.index,
    )


def safe_div(num: float, den: float) -> float:
    return float(num / den) if den else float("nan")


def evaluate_subset(df: pd.DataFrame, subset: str, policy: str, lower: float, upper: float) -> tuple[dict, pd.DataFrame]:
    work = df.copy()
    work["predicted_3class"] = assign_calls(work["pathogenic_probability"], lower, upper)
    y_true = work["model_label_3class"].astype(str)
    y_pred = work["predicted_3class"].astype(str)

    cm = pd.DataFrame(
        confusion_matrix(y_true, y_pred, labels=LABEL_ORDER),
        index=[f"true_{x}" for x in LABEL_ORDER],
        columns=[f"pred_{x}" for x in LABEL_ORDER],
    )

    true_p = y_true.eq("Pathogenic")
    true_b = y_true.eq("Benign")
    pred_p = y_pred.eq("Pathogenic")
    pred_b = y_pred.eq("Benign")
    pred_vus = y_pred.eq("VUS")

    tp = int((true_p & pred_p).sum())
    fp = int((true_b & pred_p).sum())
    tn_strict = int((true_b & pred_b).sum())
    fn_strict = int((true_p & pred_b).sum())
    p_to_vus = int((true_p & pred_vus).sum())
    b_to_vus = int((true_b & pred_vus).sum())
    p_total = int(true_p.sum())
    b_total = int(true_b.sum())
    p_calls = int(pred_p.sum())
    b_calls = int(pred_b.sum())
    vus_calls = int(pred_vus.sum())

    # Binary MCC here treats "called P/LP" as positive and "not called P/LP"
    # as negative, so deferrals count as negative calls for MCC/specificity.
    y_true_binary = true_p.astype(int)
    y_pred_binary = pred_p.astype(int)

    metrics = {
        "subset": subset,
        "policy": policy,
        "lower_benign_threshold": lower,
        "upper_pathogenic_threshold": upper,
        "rows": int(len(work)),
        "true_benign": b_total,
        "true_pathogenic": p_total,
        "exact_3class_accuracy": float(accuracy_score(y_true, y_pred)),
        "weighted_kappa": float(
            cohen_kappa_score(
                y_true.map(LABEL_TO_NUM),
                y_pred.map(LABEL_TO_NUM),
                weights="quadratic",
            )
        ),
        "macro_f1": float(f1_score(y_true, y_pred, labels=LABEL_ORDER, average="macro", zero_division=0)),
        "p_lp_sensitivity": safe_div(tp, p_total),
        "specificity_no_p_lp_escalation": safe_div(int((true_b & ~pred_p).sum()), b_total),
        "ppv_p_lp_calls": safe_div(tp, p_calls),
        "npv_b_lb_calls": safe_div(tn_strict, b_calls),
        "binary_mcc_p_call_vs_not_p": float(matthews_corrcoef(y_true_binary, y_pred_binary)),
        "vus_deferral_rate": safe_div(vus_calls, len(work)),
        "p_lp_positive_calls": p_calls,
        "true_positives": tp,
        "false_positives": fp,
        "p_lp_to_vus_deferrals": p_to_vus,
        "p_lp_to_b_lb_errors": fn_strict,
        "b_lb_to_p_lp_escalations": fp,
        "b_lb_to_vus_deferrals": b_to_vus,
        "true_negatives_strict_b_call": tn_strict,
        "false_negatives_strict_b_call": fn_strict,
    }
    return metrics, cm


def plot_confusion(cm: pd.DataFrame, title: str, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))
    arr = cm.to_numpy()
    im = ax.imshow(arr, cmap="Blues")
    ax.set_xticks(range(len(cm.columns)), labels=[c.replace("pred_", "") for c in cm.columns], rotation=35, ha="right")
    ax.set_yticks(range(len(cm.index)), labels=[i.replace("true_", "") for i in cm.index])
    ax.set_xlabel("Predicted call")
    ax.set_ylabel("True label")
    ax.set_title(title)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            ax.text(j, i, str(arr[i, j]), ha="center", va="center", color="black")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", default=str(PREDICTIONS))
    parser.add_argument("--out-dir", default=str(OUT_DIR))
    parser.add_argument("--defer-lower", type=float, default=0.10)
    parser.add_argument("--defer-upper", type=float, default=0.90)
    args = parser.parse_args()

    pred_path = Path(args.predictions)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(pred_path, sep="\t", low_memory=False)
    subsets = {
        "validation": ["validation"],
        "external_all": ["external_hiro", "external_emerge", "external_cardioboost"],
        "external_hiro": ["external_hiro"],
        "external_emerge": ["external_emerge"],
        "external_cardioboost": ["external_cardioboost"],
    }
    policies = {
        "binary_threshold_0_5": (0.5, 0.5),
        f"defer_{args.defer_lower:g}_{args.defer_upper:g}": (args.defer_lower, args.defer_upper),
    }

    metrics_rows = []
    manifest = {
        "predictions": str(pred_path.relative_to(ROOT)) if pred_path.is_relative_to(ROOT) else str(pred_path),
        "out_dir": str(out_dir.relative_to(ROOT)) if out_dir.is_relative_to(ROOT) else str(out_dir),
        "subsets": subsets,
        "policies": policies,
        "notes": [
            "Primary model is binary P/LP vs B/LB.",
            "The deferral policy is a post-hoc three-call view: B/LB below lower threshold, P/LP at or above upper threshold, VUS/defer otherwise.",
            "Exact 3-class accuracy is computed on binary-labeled rows only here, so VUS calls are counted as deferrals/errors rather than true VUS hits.",
        ],
    }

    for subset, split_names in subsets.items():
        sub = df[df["split_source_heldout"].isin(split_names)].copy()
        for policy, (lower, upper) in policies.items():
            metrics, cm = evaluate_subset(sub, subset, policy, lower, upper)
            metrics_rows.append(metrics)
            cm_path = out_dir / f"confusion_matrix_{subset}_{policy}.tsv"
            cm.to_csv(cm_path, sep="\t")
            plot_confusion(cm, f"{subset} / {policy}", out_dir / f"confusion_matrix_{subset}_{policy}.png")

    metrics_df = pd.DataFrame(metrics_rows)
    metrics_df.to_csv(out_dir / "performance_metrics_by_subset_policy.tsv", sep="\t", index=False)
    (out_dir / "evaluation_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"out_dir": manifest["out_dir"], "metrics_rows": len(metrics_df)}, indent=2))


if __name__ == "__main__":
    main()
