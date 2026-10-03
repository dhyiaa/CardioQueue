#!/usr/bin/env python3
"""Paired bootstrap confidence intervals for the CardioBoost matched benchmark."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parent
TABLES = BASE / "tables"
INPUT = TABLES / "hgvs_cdot_cardioBoost_matched_benchmark.tsv"
N_BOOT = 5000
SEED = 20260703

BENCHMARKS = {
    "hgvs_cdot_pdot_primary_external_all": lambda df: (
        df["split_source_heldout"].astype(str).str.startswith("external_")
        & df["protein_concordant"]
    ),
    "hgvs_cdot_pdot_primary_external_cardioboost": lambda df: (
        df["split_source_heldout"].astype(str).eq("external_cardioboost")
        & df["protein_concordant"]
    ),
    "hgvs_cdot_pdot_primary_external_cm": lambda df: (
        df["split_source_heldout"].astype(str).str.startswith("external_")
        & df["protein_concordant"]
        & df["cardioboost_panel"].astype(str).eq("CM")
    ),
    "hgvs_cdot_pdot_primary_external_arm": lambda df: (
        df["split_source_heldout"].astype(str).str.startswith("external_")
        & df["protein_concordant"]
        & df["cardioboost_panel"].astype(str).eq("ARM")
    ),
    "hgvs_cdot_primary_external_all": lambda df: df["split_source_heldout"].astype(str).str.startswith("external_"),
    "hgvs_cdot_primary_external_cardioboost": lambda df: df["split_source_heldout"].astype(str).eq("external_cardioboost"),
}

MODELS = {
    "CardioQueue": "our_catboost_score",
    "Official CardioBoost": "cardioboost_score",
}

METRICS = [
    "AUROC",
    "AUPRC",
    "sensitivity",
    "specificity",
    "PPV",
    "NPV",
    "deferral_rate",
    "high_conf_accuracy",
    "high_conf_classified_rate",
    "MCC_at_0_5",
]


def three_zone(score: np.ndarray) -> np.ndarray:
    calls = np.full(score.shape, "VUS", dtype=object)
    calls[score <= 0.1] = "Benign"
    calls[score >= 0.9] = "Pathogenic"
    return calls


def fast_rankdata(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    sorted_values = values[order]
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and sorted_values[end] == sorted_values[start]:
            end += 1
        avg_rank = (start + 1 + end) / 2.0
        ranks[order[start:end]] = avg_rank
        start = end
    return ranks


def fast_auroc(y: np.ndarray, s: np.ndarray) -> float:
    n_pos = int(y.sum())
    n_neg = int(len(y) - n_pos)
    if n_pos == 0 or n_neg == 0:
        return np.nan
    ranks = fast_rankdata(s)
    sum_pos_ranks = ranks[y == 1].sum()
    return float((sum_pos_ranks - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def fast_average_precision(y: np.ndarray, s: np.ndarray) -> float:
    n_pos = int(y.sum())
    if n_pos == 0:
        return np.nan
    order = np.argsort(-s, kind="mergesort")
    y_sorted = y[order]
    tp = np.cumsum(y_sorted == 1)
    ranks = np.arange(1, len(y_sorted) + 1)
    precision = tp / ranks
    return float(precision[y_sorted == 1].sum() / n_pos)


def fast_mcc(y: np.ndarray, pred: np.ndarray) -> float:
    tp = int(((y == 1) & (pred == 1)).sum())
    tn = int(((y == 0) & (pred == 0)).sum())
    fp = int(((y == 0) & (pred == 1)).sum())
    fn = int(((y == 1) & (pred == 0)).sum())
    denom = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    return float(((tp * tn) - (fp * fn)) / np.sqrt(denom)) if denom else np.nan


def metric_values(df: pd.DataFrame, score_col: str) -> dict[str, float]:
    d = df[["y_true", score_col]].copy()
    d[score_col] = pd.to_numeric(d[score_col], errors="coerce")
    d = d.dropna()
    y = d["y_true"].astype(int).to_numpy()
    s = d[score_col].to_numpy(dtype=float)
    out = {m: np.nan for m in METRICS}
    if len(y) == 0:
        return out
    if len(np.unique(y)) == 2:
        out["AUROC"] = fast_auroc(y, s)
        out["AUPRC"] = fast_average_precision(y, s)
        out["MCC_at_0_5"] = fast_mcc(y, (s >= 0.5).astype(int))

    calls = three_zone(s)
    pred_p = calls == "Pathogenic"
    pred_b = calls == "Benign"
    pred_v = calls == "VUS"
    tp = int(((y == 1) & pred_p).sum())
    fp = int(((y == 0) & pred_p).sum())
    tn = int(((y == 0) & pred_b).sum())
    fn = int(((y == 1) & pred_b).sum())
    p_defer = int(((y == 1) & pred_v).sum())
    b_defer = int(((y == 0) & pred_v).sum())

    out["sensitivity"] = tp / (tp + fn + p_defer) if (tp + fn + p_defer) else np.nan
    out["specificity"] = tn / (tn + fp + b_defer) if (tn + fp + b_defer) else np.nan
    out["PPV"] = tp / (tp + fp) if (tp + fp) else np.nan
    out["NPV"] = tn / (tn + fn) if (tn + fn) else np.nan
    out["deferral_rate"] = float(pred_v.mean())
    out["high_conf_classified_rate"] = float((pred_p | pred_b).mean())
    out["high_conf_accuracy"] = (tp + tn) / int((pred_p | pred_b).sum()) if (pred_p | pred_b).sum() else np.nan
    return out


def summarize(values: np.ndarray) -> tuple[float, float, float]:
    values = values[~np.isnan(values)]
    if len(values) == 0:
        return np.nan, np.nan, np.nan
    return (
        float(np.nanmean(values)),
        float(np.nanpercentile(values, 2.5)),
        float(np.nanpercentile(values, 97.5)),
    )


def stratified_row_indices(df: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    pieces = []
    for label in [0, 1]:
        idx = df.index[df["y_true"].eq(label)].to_numpy()
        pieces.append(rng.choice(idx, size=len(idx), replace=True))
    return np.concatenate(pieces)


def make_variant_cluster_sampler(df: pd.DataFrame):
    group_label = df.groupby("variant_id")["y_true"].nunique()
    if int(group_label.max()) != 1:
        raise ValueError("Variant cluster bootstrap requires one binary label per variant_id.")
    cluster_to_indices = {
        variant_id: rows.index.to_numpy(dtype=int)
        for variant_id, rows in df.groupby("variant_id", sort=False)
    }
    clusters_by_label = {}
    for label in [0, 1]:
        clusters_by_label[label] = (
            df.loc[df["y_true"].eq(label), ["variant_id"]]
            .drop_duplicates()["variant_id"]
            .to_numpy()
        )

    def sampler(_df: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
        pieces = []
        for label in [0, 1]:
            clusters = clusters_by_label[label]
            sampled = rng.choice(clusters, size=len(clusters), replace=True)
            pieces.append(np.concatenate([cluster_to_indices[variant_id] for variant_id in sampled]))
        return np.concatenate(pieces)

    return sampler


def run_bootstrap(df: pd.DataFrame, benchmark: str, policy: str, index_fn) -> tuple[list[dict], list[dict]]:
    rng = np.random.default_rng(SEED + abs(hash((benchmark, policy))) % 1_000_000)
    model_boot = {
        model: {metric: [] for metric in METRICS}
        for model in MODELS
        if MODELS[model] in df.columns
    }
    diff_boot = {metric: [] for metric in METRICS}

    for _ in range(N_BOOT):
        sample_idx = index_fn(df, rng)
        sample = df.loc[sample_idx].copy()
        model_metrics = {}
        for model, col in MODELS.items():
            if col not in sample.columns:
                continue
            model_metrics[model] = metric_values(sample, col)
            for metric in METRICS:
                model_boot[model][metric].append(model_metrics[model][metric])
        if "CardioQueue" in model_metrics and "Official CardioBoost" in model_metrics:
            for metric in METRICS:
                diff_boot[metric].append(
                    model_metrics["CardioQueue"][metric] - model_metrics["Official CardioBoost"][metric]
                )

    model_rows = []
    point_metrics = {
        model: metric_values(df, col)
        for model, col in MODELS.items()
        if col in df.columns
    }
    for model, metric_map in model_boot.items():
        for metric, vals in metric_map.items():
            arr = np.asarray(vals, dtype=float)
            boot_mean, ci_low, ci_high = summarize(arr)
            model_rows.append(
                {
                    "benchmark": benchmark,
                    "bootstrap_policy": policy,
                    "model": model,
                    "metric": metric,
                    "point_estimate": point_metrics[model][metric],
                    "bootstrap_mean": boot_mean,
                    "ci_95_low": ci_low,
                    "ci_95_high": ci_high,
                    "n_boot": N_BOOT,
                    "rows": len(df),
                    "unique_variants": df["variant_id"].nunique(),
                    "pathogenic": int(df["y_true"].sum()),
                    "benign": int((df["y_true"] == 0).sum()),
                }
            )

    diff_rows = []
    point_our = point_metrics["CardioQueue"]
    point_cb = point_metrics["Official CardioBoost"]
    for metric, vals in diff_boot.items():
        arr = np.asarray(vals, dtype=float)
        boot_mean, ci_low, ci_high = summarize(arr)
        diff_rows.append(
            {
                "benchmark": benchmark,
                "bootstrap_policy": policy,
                "comparison": "CardioQueue minus Official CardioBoost",
                "metric": metric,
                "point_difference": point_our[metric] - point_cb[metric],
                "bootstrap_mean_difference": boot_mean,
                "ci_95_low": ci_low,
                "ci_95_high": ci_high,
                "n_boot": N_BOOT,
                "rows": len(df),
                "unique_variants": df["variant_id"].nunique(),
                "pathogenic": int(df["y_true"].sum()),
                "benign": int((df["y_true"] == 0).sum()),
            }
        )
    return model_rows, diff_rows


def main() -> None:
    df = pd.read_csv(INPUT, sep="\t", low_memory=False)
    df = df[df["y_true"].isin([0, 1])].copy()
    df["y_true"] = df["y_true"].astype(int)
    df["protein_concordant"] = (
        df["our_pdot"].fillna("").astype(str).ne("")
        & df["our_pdot"].fillna("").astype(str).eq(df["cardioboost_pdot"].fillna("").astype(str))
    )
    df.loc[df["protein_concordant"]].to_csv(
        TABLES / "hgvs_cdot_pdot_cardioBoost_matched_benchmark.tsv", sep="\t", index=False
    )
    all_model_rows = []
    all_diff_rows = []

    for benchmark, mask_fn in BENCHMARKS.items():
        subset = df.loc[mask_fn(df)].copy().reset_index(drop=True)
        if subset.empty:
            continue
        for policy, index_fn in [
            ("row_stratified_paired", stratified_row_indices),
            ("variant_cluster_stratified_paired", make_variant_cluster_sampler(subset)),
        ]:
            model_rows, diff_rows = run_bootstrap(subset, benchmark, policy, index_fn)
            all_model_rows.extend(model_rows)
            all_diff_rows.extend(diff_rows)

    model_ci = pd.DataFrame(all_model_rows)
    diff_ci = pd.DataFrame(all_diff_rows)
    model_ci.to_csv(TABLES / "hgvs_cardioBoost_matched_bootstrap_model_ci.tsv", sep="\t", index=False)
    diff_ci.to_csv(TABLES / "hgvs_cardioBoost_matched_bootstrap_paired_difference_ci.tsv", sep="\t", index=False)

    summary = {
        "input": str(INPUT),
        "n_boot": N_BOOT,
        "seed": SEED,
        "policies": ["row_stratified_paired", "variant_cluster_stratified_paired"],
        "benchmarks": list(BENCHMARKS),
        "outputs": {
            "model_ci": str(TABLES / "hgvs_cardioBoost_matched_bootstrap_model_ci.tsv"),
            "paired_difference_ci": str(TABLES / "hgvs_cardioBoost_matched_bootstrap_paired_difference_ci.tsv"),
        },
    }
    (TABLES / "hgvs_cardioBoost_matched_bootstrap_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
