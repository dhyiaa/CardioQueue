#!/usr/bin/env python3
"""Score every January-2024 VUS and evaluate fixed-horizon review queues."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from catboost import CatBoostClassifier, Pool
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


ROOT = Path(__file__).resolve().parents[3]
TABLE = ROOT / "datasets/modeling/ready/modeling_table_temporal_clinvar_2024.tsv"
OUT = Path(__file__).resolve().parent
PAPER_TABLES = ROOT / "paper_MS_SI/tables"
PAPER_FIGURES = ROOT / "paper_MS_SI/figures"
MODELS = {
    "full_300": ROOT / "results/models/cardioqueue_temporal_clinvar_2024_in_scope_public300",
    "no_precomputed_effect_242": ROOT
    / "results/models/cardioqueue_temporal_clinvar_2024_in_scope_no_effect_predictors",
    "low_circularity": ROOT
    / "results/models/cardioqueue_temporal_clinvar_2024_in_scope_low_circularity",
    "static_molecular": ROOT
    / "results/models/cardioqueue_temporal_clinvar_2024_in_scope_static19",
}
BUDGET_ROWS = (100, 250, 500, 1000, 2500, 5000)
TARGET_RECALLS = (0.50, 0.80, 0.90, 0.95)
BOOTSTRAP_ROWS = BUDGET_ROWS
N_BOOT = 5000
N_RANDOM = 10000
SEED = 20260927


def prepare_model_matrix(frame: pd.DataFrame, model_dir: Path) -> tuple[Pool, list[str]]:
    schema = json.loads((model_dir / "feature_columns.json").read_text())
    features = schema["features"]
    categorical = schema["categorical_features"]
    missing = sorted(set(features) - set(frame.columns))
    if missing:
        raise ValueError(f"Missing model features for {model_dir.name}: {missing}")
    matrix = frame[features].copy()
    categorical_set = set(categorical)
    for column in features:
        if column in categorical_set:
            matrix[column] = matrix[column].astype("string").fillna("__MISSING__").astype(str)
        else:
            matrix[column] = pd.to_numeric(matrix[column], errors="coerce")
    return Pool(matrix, cat_features=categorical), features


def transition_outcome(current: object) -> str:
    if current == "Pathogenic":
        return "later_pathogenic"
    if current == "Benign":
        return "later_benign"
    if current == "VUS":
        return "still_vus"
    return "no_unambiguous_current_class"


def ranked(frame: pd.DataFrame, score_column: str, direction: str) -> pd.DataFrame:
    ascending = direction == "later_benign"
    return frame.sort_values(
        [score_column, "variant_id"], ascending=[ascending, True], kind="mergesort"
    )


def queue_budget_rows(
    frame: pd.DataFrame, model: str, direction: str, cohort: str
) -> list[dict]:
    score_column = f"score_{model}"
    ordered = ranked(frame, score_column, direction)
    target = ordered["transition_outcome"].eq(direction).to_numpy(dtype=int)
    total_target = int(target.sum())
    prevalence = total_target / len(ordered)
    rows = []
    for requested_rows in BUDGET_ROWS:
        reviewed = min(requested_rows, len(ordered))
        fraction = reviewed / len(ordered)
        found = int(target[:reviewed].sum())
        reviewed_outcomes = ordered.iloc[:reviewed]["transition_outcome"].value_counts()
        observed_yield = found / reviewed
        rows.append(
            {
                "model": model,
                "cohort": cohort,
                "queue_direction": direction,
                "historical_vus_rows": len(ordered),
                "target_transitions": total_target,
                "budget_fraction": fraction,
                "review_rows": reviewed,
                "target_found": found,
                "target_recall": found / total_target,
                "later_pathogenic_in_reviewed": int(reviewed_outcomes.get("later_pathogenic", 0)),
                "later_benign_in_reviewed": int(reviewed_outcomes.get("later_benign", 0)),
                "still_vus_in_reviewed": int(reviewed_outcomes.get("still_vus", 0)),
                "no_unambiguous_current_class_in_reviewed": int(
                    reviewed_outcomes.get("no_unambiguous_current_class", 0)
                ),
                "observed_transition_yield_lower_bound": observed_yield,
                "enrichment_over_transition_prevalence": observed_yield / prevalence,
                "random_expected_recall": fraction,
            }
        )
    return rows


@lru_cache(maxsize=None)
def random_target_workload(
    population_rows: int, total_target: int, required: int
) -> tuple[float, float, float]:
    rng = np.random.default_rng(SEED + total_target + required)
    values = np.empty(N_RANDOM, dtype=int)
    for index in range(N_RANDOM):
        target_positions = rng.choice(population_rows, total_target, replace=False)
        values[index] = np.partition(target_positions, required - 1)[required - 1] + 1
    return tuple(float(value) for value in np.quantile(values, [0.50, 0.025, 0.975]))


def queue_target_rows(
    frame: pd.DataFrame, model: str, direction: str, cohort: str
) -> list[dict]:
    score_column = f"score_{model}"
    ordered = ranked(frame, score_column, direction)
    target = ordered["transition_outcome"].eq(direction).to_numpy(dtype=int)
    total_target = int(target.sum())
    prevalence = total_target / len(ordered)
    cumulative = np.cumsum(target)
    rows = []
    for target_recall in TARGET_RECALLS:
        required = int(np.ceil(target_recall * total_target))
        reviews = int(np.flatnonzero(cumulative >= required)[0] + 1)
        random_median, random_low, random_high = random_target_workload(
            len(ordered), total_target, required
        )
        rows.append(
            {
                "model": model,
                "cohort": cohort,
                "queue_direction": direction,
                "target_recall": target_recall,
                "required_transitions": required,
                "review_rows": reviews,
                "review_fraction": reviews / len(ordered),
                "oracle_minimum_reviews": required,
                "random_expected_reviews": required / prevalence,
                "random_median_reviews": random_median,
                "random_95_low_reviews": random_low,
                "random_95_high_reviews": random_high,
                "relative_reduction_vs_random_expectation": 1 - reviews / (required / prevalence),
            }
        )
    return rows


def two_ended_budget_rows(frame: pd.DataFrame, model: str, cohort: str) -> list[dict]:
    score_column = f"score_{model}"
    ordered = frame.sort_values(
        [score_column, "variant_id"], ascending=[False, True], kind="mergesort"
    )
    total_resolved = int(
        ordered["transition_outcome"].isin(["later_pathogenic", "later_benign"]).sum()
    )
    rows = []
    for requested_rows in BUDGET_ROWS:
        reviewed = min(requested_rows, len(ordered))
        high_rows = (reviewed + 1) // 2
        low_rows = reviewed // 2
        selected = pd.concat([ordered.iloc[:high_rows], ordered.iloc[len(ordered) - low_rows :]])
        outcomes = selected["transition_outcome"].value_counts()
        pathogenic = int(outcomes.get("later_pathogenic", 0))
        benign = int(outcomes.get("later_benign", 0))
        resolved = pathogenic + benign
        rows.append(
            {
                "model": model,
                "cohort": cohort,
                "review_rows": reviewed,
                "high_score_rows": high_rows,
                "low_score_rows": low_rows,
                "later_pathogenic_found": pathogenic,
                "later_benign_found": benign,
                "any_binary_transition_found": resolved,
                "any_binary_transition_recall": resolved / total_resolved,
                "observed_resolution_yield_lower_bound": resolved / reviewed,
                "still_vus_in_reviewed": int(outcomes.get("still_vus", 0)),
                "no_unambiguous_current_class_in_reviewed": int(
                    outcomes.get("no_unambiguous_current_class", 0)
                ),
            }
        )
    return rows


def score_decile_audit(frame: pd.DataFrame) -> pd.DataFrame:
    ordered = frame.sort_values(
        ["score_full_300", "variant_id"], ascending=[False, True], kind="mergesort"
    ).copy()
    ordered["score_decile"] = np.minimum(
        10, np.floor(np.arange(len(ordered)) * 10 / len(ordered)).astype(int) + 1
    )
    rows = []
    for decile, group in ordered.groupby("score_decile", sort=True):
        outcomes = group["transition_outcome"].value_counts()
        rows.append(
            {
                "cohort": "in_scope_primary",
                "score_decile": int(decile),
                "queue_position": "highest scores" if decile == 1 else ("lowest scores" if decile == 10 else "middle"),
                "rows": len(group),
                "score_min": group["score_full_300"].min(),
                "score_max": group["score_full_300"].max(),
                "later_pathogenic": int(outcomes.get("later_pathogenic", 0)),
                "later_benign": int(outcomes.get("later_benign", 0)),
                "still_vus": int(outcomes.get("still_vus", 0)),
                "no_unambiguous_current_class": int(
                    outcomes.get("no_unambiguous_current_class", 0)
                ),
            }
        )
    return pd.DataFrame(rows)


def resolved_metrics(frame: pd.DataFrame, model: str, cohort: str) -> dict:
    resolved = frame[frame["transition_outcome"].isin(["later_pathogenic", "later_benign"])]
    y = resolved["transition_outcome"].eq("later_pathogenic").to_numpy(dtype=int)
    score = resolved[f"score_{model}"].to_numpy(dtype=float)
    return {
        "model": model,
        "cohort": cohort,
        "rows": len(resolved),
        "later_pathogenic": int(y.sum()),
        "later_benign": int((1 - y).sum()),
        "auroc": roc_auc_score(y, score),
        "auprc": average_precision_score(y, score),
        "brier": brier_score_loss(y, score),
        "sensitivity_at_0_5": float(((score >= 0.5) & (y == 1)).sum() / y.sum()),
        "specificity_at_0_5": float(((score < 0.5) & (y == 0)).sum() / (y == 0).sum()),
    }


def gene_cluster_bootstrap(frame: pd.DataFrame, model: str) -> pd.DataFrame:
    score_column = f"score_{model}"
    genes = frame["primary_gene"].fillna("__MISSING_GENE__").astype(str)
    gene_names = genes.unique()
    positions = {gene: np.flatnonzero(genes.to_numpy() == gene) for gene in gene_names}
    rng = np.random.default_rng(SEED)
    scores = frame[score_column].to_numpy(dtype=float)
    outcomes = frame["transition_outcome"].eq("later_pathogenic").to_numpy(dtype=int)
    samples = {rows: [] for rows in BOOTSTRAP_ROWS}
    for _ in range(N_BOOT):
        sampled_genes = rng.choice(gene_names, len(gene_names), replace=True)
        index = np.concatenate([positions[gene] for gene in sampled_genes])
        order = np.argsort(-scores[index], kind="stable")
        target = outcomes[index][order]
        positives = int(target.sum())
        if positives == 0:
            continue
        cumulative = np.cumsum(target)
        for requested_rows in BOOTSTRAP_ROWS:
            reviewed = min(requested_rows, len(index))
            samples[requested_rows].append(cumulative[reviewed - 1] / positives)
    point_rows = queue_budget_rows(frame, model, "later_pathogenic", "in_scope_primary")
    point = {row["review_rows"]: row["target_recall"] for row in point_rows}
    return pd.DataFrame(
        [
            {
                "model": model,
                "cohort": "in_scope_primary",
                "queue_direction": "later_pathogenic",
                "review_rows": requested_rows,
                "budget_fraction": requested_rows / len(frame),
                "target_recall": point[requested_rows],
                "gene_cluster_ci_low": float(np.quantile(samples[requested_rows], 0.025)),
                "gene_cluster_ci_high": float(np.quantile(samples[requested_rows], 0.975)),
                "bootstrap_replicates": len(samples[requested_rows]),
                "seed": SEED,
            }
            for requested_rows in BOOTSTRAP_ROWS
        ]
    )


def cumulative_curve(frame: pd.DataFrame, model: str, direction: str, grid: np.ndarray) -> np.ndarray:
    ordered = ranked(frame, f"score_{model}", direction)
    target = ordered["transition_outcome"].eq(direction).to_numpy(dtype=int)
    cumulative = np.concatenate([[0], np.cumsum(target)])
    rows = np.ceil(grid * len(ordered)).astype(int)
    return cumulative[rows] / target.sum()


def make_figure(frame: pd.DataFrame) -> None:
    grid = np.linspace(0, 0.20, 201)
    colors = {
        "full_300": "#0B6E75",
        "no_precomputed_effect_242": "#B2472D",
        "low_circularity": "#3C6E47",
        "static_molecular": "#6B5B95",
    }
    labels = {
        "full_300": "Full temporal model (300 features)",
        "no_precomputed_effect_242": "Without 8 predictor families (242 features)",
        "low_circularity": "Low-circularity sensitivity",
        "static_molecular": "Static molecular sensitivity",
    }
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    fig, axes = plt.subplots(1, 3, figsize=(14.6, 4.8))
    outcome_order = [
        "later_pathogenic",
        "later_benign",
        "still_vus",
        "no_unambiguous_current_class",
    ]
    outcome_labels = ["Later P/LP", "Later B/LB", "Still VUS", "No current class"]
    outcome_colors = ["#B2472D", "#3C6E47", "#707070", "#C7A24B"]
    counts = frame["transition_outcome"].value_counts()
    values = [int(counts.get(outcome, 0)) for outcome in outcome_order]
    axes[0].bar(outcome_labels, values, color=outcome_colors, width=0.72)
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Variants (log scale)")
    axes[0].set_title("A  Fixed-horizon cohort accounting", loc="left", fontweight="bold")
    axes[0].tick_params(axis="x", rotation=28)
    axes[0].grid(axis="y", color="#D8D8D8", linewidth=0.7, alpha=0.8)
    axes[0].spines[["top", "right"]].set_visible(False)
    for index, value in enumerate(values):
        axes[0].text(index, value * 1.12, f"{value:,}", ha="center", va="bottom", fontsize=8.5)

    for ax, direction, title in [
        (axes[1], "later_pathogenic", "B  Future P/LP from top of queue"),
        (axes[2], "later_benign", "C  Future B/LB from bottom of queue"),
    ]:
        for model in MODELS:
            ax.plot(
                grid,
                cumulative_curve(frame, model, direction, grid),
                color=colors[model],
                linewidth=2.3,
                label=labels[model],
            )
        ax.plot(grid, grid, color="#707070", linewidth=1.6, linestyle="--", label="Unprioritized order")
        ax.set_title(title, loc="left", fontweight="bold")
        ax.set_xlim(0, 0.20)
        ax.set_ylim(0, 1.0)
        ax.set_xlabel("Fraction of all January-2024 VUS reviewed")
        ax.grid(color="#D8D8D8", linewidth=0.7, alpha=0.8)
        ax.spines[["top", "right"]].set_visible(False)
    axes[1].set_ylabel("Fraction of later transitions recovered")
    axes[2].legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=7.8)
    fig.tight_layout()
    for path in [
        OUT / "all_historical_vus_bidirectional_gain.png",
        PAPER_FIGURES / "Figure_S6_all_historical_vus_bidirectional_gain.png",
    ]:
        fig.savefig(path, dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "all_historical_vus_bidirectional_gain.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    PAPER_TABLES.mkdir(parents=True, exist_ok=True)
    PAPER_FIGURES.mkdir(parents=True, exist_ok=True)
    table = pd.read_csv(TABLE, sep="\t", low_memory=False)
    frame = table[table["clinvar_2024_class"].eq("VUS")].copy()
    if len(frame) != 25_737 or not frame["variant_id"].is_unique:
        raise ValueError("Expected 25,737 unique January-2024 VUS")
    frame["transition_outcome"] = frame["clinvar_current_class"].map(transition_outcome)

    feature_counts = {}
    for name, model_dir in MODELS.items():
        model = CatBoostClassifier()
        model.load_model(str(model_dir / "primary_binary_catboost_model.cbm"))
        pool, features = prepare_model_matrix(frame, model_dir)
        frame[f"score_{name}"] = model.predict_proba(pool)[:, 1]
        feature_counts[name] = len(features)

    prediction_columns = [
        "variant_id",
        "primary_gene",
        "primary_model_inclusion",
        "clinvar_2024_class",
        "clinvar_current_class",
        "transition_outcome",
        *[f"score_{name}" for name in MODELS],
    ]
    frame[prediction_columns].to_csv(OUT / "all_historical_vus_predictions.tsv", sep="\t", index=False)

    cohorts = {
        "in_scope_primary": frame[frame["primary_model_inclusion"].eq("include")].copy(),
        "all_historical_vus_sensitivity": frame,
    }
    count_rows = []
    for cohort, cohort_frame in cohorts.items():
        for outcome, variants in cohort_frame["transition_outcome"].value_counts().items():
            count_rows.append(
                {
                    "cohort": cohort,
                    "transition_outcome": outcome,
                    "variants": int(variants),
                    "fraction_of_historical_vus": variants / len(cohort_frame),
                }
            )
    counts = pd.DataFrame(count_rows)
    counts["included_in_resolved_binary_metrics"] = counts["transition_outcome"].isin(
        ["later_pathogenic", "later_benign"]
    )
    budget = pd.DataFrame(
        [
            row
            for cohort, cohort_frame in cohorts.items()
            for model in MODELS
            for direction in ("later_pathogenic", "later_benign")
            for row in queue_budget_rows(cohort_frame, model, direction, cohort)
        ]
    )
    targets = pd.DataFrame(
        [
            row
            for cohort, cohort_frame in cohorts.items()
            for model in MODELS
            for direction in ("later_pathogenic", "later_benign")
            for row in queue_target_rows(cohort_frame, model, direction, cohort)
        ]
    )
    resolved = pd.DataFrame(
        [
            resolved_metrics(cohort_frame, model, cohort)
            for cohort, cohort_frame in cohorts.items()
            for model in MODELS
        ]
    )
    primary_frame = cohorts["in_scope_primary"]
    two_ended = pd.DataFrame(
        [
            row
            for cohort, cohort_frame in cohorts.items()
            for model in MODELS
            for row in two_ended_budget_rows(cohort_frame, model, cohort)
        ]
    )
    deciles = score_decile_audit(primary_frame)
    bootstrap = pd.concat(
        [
            gene_cluster_bootstrap(primary_frame, model)
            for model in ("full_300", "low_circularity")
        ],
        ignore_index=True,
    )

    outputs = [
        (counts, "all_historical_vus_transition_counts.tsv", "Table_S19_all_vus_transition_counts.tsv"),
        (budget, "all_historical_vus_queue_budgets.tsv", "Table_S19b_all_vus_queue_budgets.tsv"),
        (targets, "all_historical_vus_queue_targets.tsv", "Table_S19c_all_vus_queue_targets.tsv"),
        (resolved, "resolved_vus_model_metrics.tsv", "Table_S19d_resolved_vus_model_metrics.tsv"),
        (bootstrap, "all_historical_vus_gene_bootstrap.tsv", "Table_S19e_all_vus_gene_bootstrap.tsv"),
        (deciles, "all_historical_vus_score_deciles.tsv", "Table_S19f_all_vus_score_deciles.tsv"),
        (two_ended, "all_historical_vus_two_ended_queue.tsv", "Table_S19g_all_vus_two_ended_queue.tsv"),
    ]
    for data, local_name, paper_name in outputs:
        data.to_csv(OUT / local_name, sep="\t", index=False)
        data.to_csv(PAPER_TABLES / paper_name, sep="\t", index=False)

    make_figure(primary_frame)
    summary = {
        "historical_vus_rows": len(frame),
        "in_scope_historical_vus_rows": len(primary_frame),
        "transition_counts_all": frame["transition_outcome"].value_counts().to_dict(),
        "transition_counts_in_scope": primary_frame["transition_outcome"].value_counts().to_dict(),
        "models": feature_counts,
        "bootstrap_replicates": N_BOOT,
        "bootstrap_unit": "primary_gene",
        "random_orderings": N_RANDOM,
        "primary_endpoint": "Recovery of variants observed as P/LP by the current release from the complete January-2024 VUS queue.",
        "resolved_only_endpoint": "Secondary discrimination of later P/LP versus later B/LB among variants with a current binary class.",
        "caution": "Still-VUS and currently unresolved records are not treated as benign; reclassification is influenced by scrutiny and evidence availability.",
    }
    (OUT / "all_historical_vus_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print("\nBudget results\n", budget.to_string(index=False))
    print("\nResolved-only metrics\n", resolved.to_string(index=False))


if __name__ == "__main__":
    main()
