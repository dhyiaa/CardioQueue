#!/usr/bin/env python3
"""Evaluate score-guided review queues on source-held-out variants."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
PREDICTIONS = (
    ROOT
    / "results/models/cardioqueue_v1_source_heldout"
    / "primary_binary_predictions_split_source_heldout.tsv"
)
TEMPORAL_PREDICTIONS = (
    ROOT
    / "results/model_performance/temporal_clinvar_2024"
    / "all_historical_vus_predictions.tsv"
)
MATCHED = (
    ROOT
    / "results/model_performance/cardioboost_matched_benchmark/tables"
    / "hgvs_cdot_cardioBoost_matched_benchmark.tsv"
)
OUT = Path(__file__).resolve().parent
PAPER_TABLES = ROOT / "paper_MS_SI/tables"
PAPER_FIGURES = ROOT / "paper_MS_SI/figures"
SEED = 20260917
N_BOOT = 5000
N_RANDOM = 10000
BUDGETS = (0.05, 0.10, 0.20, 0.30, 0.40, 0.50)
TARGETS = (0.50, 0.80, 0.90, 0.95)
CURVE_GRID = np.linspace(0.0, 1.0, 101)


def ranked(frame: pd.DataFrame, score_column: str) -> pd.DataFrame:
    return frame.sort_values([score_column, "variant_id"], ascending=[False, True], kind="mergesort")


def budget_values(frame: pd.DataFrame, score_column: str, fraction: float) -> tuple[int, int, float, float, float]:
    ordered = ranked(frame, score_column)
    rows = max(1, int(np.ceil(fraction * len(ordered))))
    found = int(ordered.iloc[:rows]["y_true"].sum())
    positives = int(ordered["y_true"].sum())
    prevalence = positives / len(ordered)
    recall = found / positives
    precision = found / rows
    enrichment = precision / prevalence
    return rows, found, recall, precision, enrichment


def rows_to_target(frame: pd.DataFrame, score_column: str, target: float) -> int:
    ordered = ranked(frame, score_column)
    required = int(np.ceil(target * ordered["y_true"].sum()))
    cumulative = ordered["y_true"].cumsum().to_numpy()
    return int(np.flatnonzero(cumulative >= required)[0] + 1)


def stratified_bootstrap(frame: pd.DataFrame, score_column: str) -> dict[float, dict[str, np.ndarray]]:
    rng = np.random.default_rng(SEED)
    y = frame["y_true"].to_numpy(dtype=int)
    groups = [np.flatnonzero(y == value) for value in (0, 1)]
    values = {
        budget: {"recall": [], "precision": [], "enrichment": []}
        for budget in BUDGETS
    }
    for _ in range(N_BOOT):
        indices = np.concatenate([rng.choice(group, len(group), replace=True) for group in groups])
        sample = frame.iloc[indices].copy()
        for budget in BUDGETS:
            _, _, recall, precision, enrichment = budget_values(sample, score_column, budget)
            values[budget]["recall"].append(recall)
            values[budget]["precision"].append(precision)
            values[budget]["enrichment"].append(enrichment)
    return {
        budget: {metric: np.asarray(samples) for metric, samples in metrics.items()}
        for budget, metrics in values.items()
    }


def random_order_simulation(frame: pd.DataFrame) -> tuple[dict[float, dict[str, np.ndarray]], dict[float, np.ndarray], np.ndarray]:
    rng = np.random.default_rng(SEED)
    y = frame["y_true"].to_numpy(dtype=int)
    positives = int(y.sum())
    prevalence = positives / len(y)
    budget_samples = {
        budget: {"recall": [], "precision": [], "enrichment": []}
        for budget in BUDGETS
    }
    target_samples = {target: [] for target in TARGETS}
    curves = []
    grid_rows = np.ceil(CURVE_GRID * len(y)).astype(int)
    for _ in range(N_RANDOM):
        permuted = y[rng.permutation(len(y))]
        cumulative = np.concatenate([[0], np.cumsum(permuted)])
        curves.append(cumulative[grid_rows] / positives)
        for budget in BUDGETS:
            rows = max(1, int(np.ceil(budget * len(y))))
            found = int(cumulative[rows])
            precision = found / rows
            budget_samples[budget]["recall"].append(found / positives)
            budget_samples[budget]["precision"].append(precision)
            budget_samples[budget]["enrichment"].append(precision / prevalence)
        for target in TARGETS:
            required = int(np.ceil(target * positives))
            target_samples[target].append(int(np.flatnonzero(cumulative >= required)[0]))
    budgets = {
        budget: {metric: np.asarray(samples) for metric, samples in metrics.items()}
        for budget, metrics in budget_samples.items()
    }
    targets = {target: np.asarray(samples) for target, samples in target_samples.items()}
    return budgets, targets, np.asarray(curves)


def score_curve_bootstrap(frame: pd.DataFrame, score_column: str) -> np.ndarray:
    rng = np.random.default_rng(SEED)
    y = frame["y_true"].to_numpy(dtype=int)
    groups = [np.flatnonzero(y == value) for value in (0, 1)]
    curves = []
    for _ in range(N_BOOT):
        indices = np.concatenate([rng.choice(group, len(group), replace=True) for group in groups])
        ordered = ranked(frame.iloc[indices].copy(), score_column)
        cumulative = np.concatenate([[0], np.cumsum(ordered["y_true"].to_numpy(dtype=int))])
        grid_rows = np.ceil(CURVE_GRID * len(ordered)).astype(int)
        curves.append(cumulative[grid_rows] / ordered["y_true"].sum())
    return np.asarray(curves)


def queue_tables(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray]:
    score_boot = stratified_bootstrap(frame, "pathogenic_probability")
    random_budgets, random_targets, random_curves = random_order_simulation(frame)
    budget_rows = []
    positives = int(frame["y_true"].sum())
    prevalence = positives / len(frame)
    for budget in BUDGETS:
        reviews, found, recall, precision, enrichment = budget_values(
            frame, "pathogenic_probability", budget
        )
        row = {
            "queue": "CardioQueue score",
            "budget_fraction": budget,
            "review_rows": reviews,
            "pathogenic_found": found,
            "pathogenic_recall": recall,
            "precision": precision,
            "enrichment_over_prevalence": enrichment,
            "pathogenic_recall_ci_95_low": np.quantile(score_boot[budget]["recall"], 0.025),
            "pathogenic_recall_ci_95_high": np.quantile(score_boot[budget]["recall"], 0.975),
            "precision_ci_95_low": np.quantile(score_boot[budget]["precision"], 0.025),
            "precision_ci_95_high": np.quantile(score_boot[budget]["precision"], 0.975),
            "resampling_interpretation": "outcome-stratified bootstrap",
        }
        budget_rows.append(row)
        random = random_budgets[budget]
        budget_rows.append(
            {
                "queue": "Unprioritized random order",
                "budget_fraction": budget,
                "review_rows": reviews,
                "pathogenic_found": reviews * prevalence,
                "pathogenic_recall": budget,
                "precision": prevalence,
                "enrichment_over_prevalence": 1.0,
                "pathogenic_recall_ci_95_low": np.quantile(random["recall"], 0.025),
                "pathogenic_recall_ci_95_high": np.quantile(random["recall"], 0.975),
                "precision_ci_95_low": np.quantile(random["precision"], 0.025),
                "precision_ci_95_high": np.quantile(random["precision"], 0.975),
                "resampling_interpretation": "random-order simulation interval",
            }
        )

    target_rows = []
    for target in TARGETS:
        score_rows = rows_to_target(frame, "pathogenic_probability", target)
        random = random_targets[target]
        random_median = float(np.median(random))
        target_rows.append(
            {
                "target_pathogenic_recall": target,
                "required_pathogenic_variants": int(np.ceil(target * positives)),
                "score_ranked_review_rows": score_rows,
                "score_ranked_review_fraction": score_rows / len(frame),
                "random_review_rows_median": random_median,
                "random_review_rows_ci_95_low": np.quantile(random, 0.025),
                "random_review_rows_ci_95_high": np.quantile(random, 0.975),
                "reviews_avoided_vs_random_median": random_median - score_rows,
                "relative_workload_reduction": 1 - score_rows / random_median,
            }
        )
    return pd.DataFrame(budget_rows), pd.DataFrame(target_rows), score_curve_bootstrap(
        frame, "pathogenic_probability"
    ), random_curves


def source_budget_table(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for source, subset in frame.groupby("split_source_heldout", sort=True):
        for budget in (0.10, 0.20, 0.30, 0.50):
            reviews, found, recall, precision, enrichment = budget_values(
                subset, "pathogenic_probability", budget
            )
            rows.append(
                {
                    "source": source.removeprefix("external_"),
                    "rows": len(subset),
                    "pathogenic": int(subset["y_true"].sum()),
                    "budget_fraction": budget,
                    "review_rows": reviews,
                    "pathogenic_found": found,
                    "pathogenic_recall": recall,
                    "precision": precision,
                    "enrichment_over_prevalence": enrichment,
                }
            )
    return pd.DataFrame(rows)


def matched_comparator_table() -> pd.DataFrame:
    frame = pd.read_csv(MATCHED, sep="\t", low_memory=False)
    frame = frame[frame["split_source_heldout"].astype(str).str.startswith("external_")].copy()
    models = {
        "CardioQueue": "our_catboost_score",
        "Official CardioBoost": "cardioboost_score",
    }
    rows = []
    for model, column in models.items():
        for budget in (0.10, 0.20, 0.30, 0.50):
            reviews, found, recall, precision, enrichment = budget_values(frame, column, budget)
            rows.append(
                {
                    "model": model,
                    "metric_type": "fixed_review_budget",
                    "target": budget,
                    "review_rows": reviews,
                    "pathogenic_found": found,
                    "pathogenic_recall": recall,
                    "precision": precision,
                    "enrichment_over_prevalence": enrichment,
                    "rows": len(frame),
                    "unique_variants": frame["variant_id"].nunique(),
                }
            )
        for target in TARGETS:
            reviews = rows_to_target(frame, column, target)
            rows.append(
                {
                    "model": model,
                    "metric_type": "target_pathogenic_recall",
                    "target": target,
                    "review_rows": reviews,
                    "pathogenic_found": int(np.ceil(target * frame["y_true"].sum())),
                    "pathogenic_recall": target,
                    "precision": np.nan,
                    "enrichment_over_prevalence": np.nan,
                    "rows": len(frame),
                    "unique_variants": frame["variant_id"].nunique(),
                }
            )
    return pd.DataFrame(rows)


def plot_gain_curve(frame: pd.DataFrame, score_curves: np.ndarray, random_curves: np.ndarray) -> None:
    ordered = ranked(frame, "pathogenic_probability")
    cumulative = np.concatenate([[0], np.cumsum(ordered["y_true"].to_numpy(dtype=int))])
    grid_rows = np.ceil(CURVE_GRID * len(ordered)).astype(int)
    score_curve = cumulative[grid_rows] / ordered["y_true"].sum()
    score_low, score_high = np.quantile(score_curves, [0.025, 0.975], axis=0)
    random_low, random_high = np.quantile(random_curves, [0.025, 0.975], axis=0)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    temporal = pd.read_csv(TEMPORAL_PREDICTIONS, sep="\t", low_memory=False)
    temporal = temporal[temporal["primary_model_inclusion"].eq("include")].copy()
    temporal["y_true"] = temporal["transition_outcome"].eq("later_pathogenic").astype(int)
    temporal["pathogenic_probability"] = temporal["score_full_300"]
    temporal = ranked(temporal, "pathogenic_probability")
    temporal_cumulative = np.concatenate([[0], np.cumsum(temporal["y_true"].to_numpy(dtype=int))])
    temporal_grid_rows = np.ceil(CURVE_GRID * len(temporal)).astype(int)
    temporal_curve = temporal_cumulative[temporal_grid_rows] / temporal["y_true"].sum()

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8), sharex=True, sharey=True)
    ax = axes[0]
    ax.fill_between(CURVE_GRID, score_low, score_high, color="#0B6E75", alpha=0.16, linewidth=0)
    ax.plot(CURVE_GRID, score_curve, color="#0B6E75", linewidth=2.4, label="CardioQueue score")
    ax.fill_between(CURVE_GRID, random_low, random_high, color="#707070", alpha=0.14, linewidth=0)
    ax.plot(CURVE_GRID, CURVE_GRID, color="#707070", linewidth=1.8, linestyle="--", label="Unprioritized order")
    external_reviews = rows_to_target(frame, "pathogenic_probability", 0.80)
    ax.scatter(external_reviews / len(frame), 0.80, color="#B2472D", s=38, zorder=5)
    ax.annotate(
        f"80% recovered\n{external_reviews}/{len(frame)} reviewed",
        xy=(external_reviews / len(frame), 0.80), xytext=(0.58, 0.68),
        arrowprops={"arrowstyle": "->", "color": "#5A5A5A", "lw": 1}, fontsize=9,
    )
    ax.set_title("A  Source-heldout variants", loc="left", fontweight="bold")
    ax.set(xlim=(0, 1), ylim=(0, 1), xlabel="Fraction of queue reviewed", ylabel="Fraction of P/LP variants recovered")
    ax.grid(color="#D8D8D8", linewidth=0.7, alpha=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="lower right")

    ax = axes[1]
    ax.plot(CURVE_GRID, temporal_curve, color="#3C6E47", linewidth=2.4, label="CardioQueue temporal score")
    ax.plot(CURVE_GRID, CURVE_GRID, color="#707070", linewidth=1.8, linestyle="--", label="Unprioritized order")
    temporal_5_rows = int(np.ceil(0.05 * len(temporal)))
    temporal_5_found = int(temporal.iloc[:temporal_5_rows]["y_true"].sum())
    temporal_5_recall = temporal_5_found / temporal["y_true"].sum()
    temporal_80_reviews = rows_to_target(temporal, "pathogenic_probability", 0.80)
    ax.scatter(temporal_5_rows / len(temporal), temporal_5_recall, color="#B2472D", s=38, zorder=5)
    ax.annotate(
        f"First 5%: {temporal_5_found}/{int(temporal['y_true'].sum())}\nlater-P/LP recovered",
        xy=(temporal_5_rows / len(temporal), temporal_5_recall), xytext=(0.18, 0.45),
        arrowprops={"arrowstyle": "->", "color": "#5A5A5A", "lw": 1}, fontsize=9,
    )
    ax.scatter(temporal_80_reviews / len(temporal), 0.80, color="#B2472D", s=38, zorder=5)
    ax.annotate(
        f"80% recovered\n{temporal_80_reviews}/{len(temporal)} reviewed",
        xy=(temporal_80_reviews / len(temporal), 0.80), xytext=(0.42, 0.78),
        arrowprops={"arrowstyle": "->", "color": "#5A5A5A", "lw": 1}, fontsize=9,
    )
    ax.set_title(f"B  January 2024 VUS cohort (n={len(temporal):,})", loc="left", fontweight="bold")
    ax.set(xlim=(0, 1), ylim=(0, 1), xlabel="Fraction of queue reviewed")
    ax.grid(color="#D8D8D8", linewidth=0.7, alpha=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    for path in [OUT / "review_queue_cumulative_gain.png", PAPER_FIGURES / "Figure_5_review_queue_cumulative_gain.png"]:
        fig.savefig(path, dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "review_queue_cumulative_gain.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    PAPER_TABLES.mkdir(parents=True, exist_ok=True)
    PAPER_FIGURES.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(PREDICTIONS, sep="\t")
    frame = frame[frame["split_source_heldout"].str.startswith("external_")].copy()
    if len(frame) != 430 or int(frame["y_true"].sum()) != 255:
        raise ValueError("Expected 430 held-out variants, including 255 P/LP")

    budget, target, score_curves, random_curves = queue_tables(frame)
    source = source_budget_table(frame)
    matched = matched_comparator_table()
    outputs = [
        (budget, "review_queue_budget_metrics.tsv", "Table_S16_review_queue_budget_metrics.tsv"),
        (target, "review_queue_target_workload.tsv", "Table_S16b_review_queue_target_workload.tsv"),
        (source, "review_queue_source_metrics.tsv", "Table_S16c_review_queue_source_metrics.tsv"),
        (matched, "review_queue_matched_cardioboost.tsv", "Table_S16d_review_queue_matched_cardioboost.tsv"),
    ]
    for table, local_name, paper_name in outputs:
        table.to_csv(OUT / local_name, sep="\t", index=False)
        table.to_csv(PAPER_TABLES / paper_name, sep="\t", index=False)
    plot_gain_curve(frame, score_curves, random_curves)

    budget_index = budget[budget["queue"].eq("CardioQueue score")].set_index("budget_fraction")
    target_index = target.set_index("target_pathogenic_recall")
    matched_budget = matched[matched["metric_type"].eq("fixed_review_budget")]
    matched_target = matched[matched["metric_type"].eq("target_pathogenic_recall")]
    report = f"""# Held-out variant review-queue simulation

This retrospective simulation ordered 430 source-held-out variants by the frozen CardioQueue score. The comparison queue used 10,000 random permutations as a proxy for arrival-order review when arrival order is unrelated to pathogenicity. No observed timestamps, staffing data, turnaround times, or patient outcomes were available.

At a 20% review budget, the score-ranked queue reviewed {int(budget_index.loc[0.20, 'review_rows'])} variants and recovered {int(budget_index.loc[0.20, 'pathogenic_found'])} P/LP variants, corresponding to {budget_index.loc[0.20, 'pathogenic_recall']:.1%} recall and {budget_index.loc[0.20, 'precision']:.1%} precision. At a 50% budget, it recovered {int(budget_index.loc[0.50, 'pathogenic_found'])} of 255 P/LP variants ({budget_index.loc[0.50, 'pathogenic_recall']:.1%}) with {budget_index.loc[0.50, 'precision']:.1%} precision.

Recovering 80% of P/LP variants required {int(target_index.loc[0.80, 'score_ranked_review_rows'])} score-ranked reviews, compared with a median {int(target_index.loc[0.80, 'random_review_rows_median'])} reviews under random ordering. The observed workload reduction was {target_index.loc[0.80, 'relative_workload_reduction']:.1%}. At 95% recovery, the corresponding counts were {int(target_index.loc[0.95, 'score_ranked_review_rows'])} and {int(target_index.loc[0.95, 'random_review_rows_median'])}, a {target_index.loc[0.95, 'relative_workload_reduction']:.1%} reduction.

On matched observations, CardioQueue recovered {matched_budget[(matched_budget['model'].eq('CardioQueue')) & (matched_budget['target'].eq(0.50))]['pathogenic_recall'].iloc[0]:.1%} of P/LP observations in the first half of the queue, compared with {matched_budget[(matched_budget['model'].eq('Official CardioBoost')) & (matched_budget['target'].eq(0.50))]['pathogenic_recall'].iloc[0]:.1%} for official CardioBoost. Reaching 90% recovery required {int(matched_target[(matched_target['model'].eq('CardioQueue')) & (matched_target['target'].eq(0.90))]['review_rows'].iloc[0])} versus {int(matched_target[(matched_target['model'].eq('Official CardioBoost')) & (matched_target['target'].eq(0.90))]['review_rows'].iloc[0])} observations.

These are worklist-efficiency estimates on a consensus-filtered retrospective cohort. They do not measure emergency triage, diagnostic yield, clinical action, or turnaround time.
"""
    (OUT / "REVIEW_QUEUE_REPORT.md").write_text(report)

    manifest = {
        "prediction_file": str(PREDICTIONS.relative_to(ROOT)),
        "external_rows": len(frame),
        "pathogenic_rows": int(frame["y_true"].sum()),
        "benign_rows": int((frame["y_true"] == 0).sum()),
        "budgets": BUDGETS,
        "targets": TARGETS,
        "bootstrap_samples": N_BOOT,
        "random_order_simulations": N_RANDOM,
        "seed": SEED,
        "random_order_interpretation": "Proxy for an arrival-order queue only when arrival order is unrelated to pathogenicity.",
    }
    (OUT / "review_queue_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
