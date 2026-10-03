#!/usr/bin/env python3
"""Evaluate review ordering among January 2024 ClinVar VUS resolved by the current release."""

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
PRED = ROOT / "results/models/cardioqueue_temporal_clinvar_2024_public300/primary_binary_predictions_split_temporal_2024.tsv"
OUT = Path(__file__).resolve().parent
PAPER = ROOT / "paper_MS_SI/tables"
TARGETS = (0.50, 0.80, 0.90, 0.95)
BUDGETS = (0.05, 0.10, 0.20, 0.30, 0.50)


def main() -> None:
    data = pd.read_csv(PRED, sep="\t")
    data = data[data["split_temporal_2024"].eq("external_temporal_vus")].copy()
    data = data.sort_values(
        ["pathogenic_probability", "variant_id"], ascending=[False, True], kind="mergesort"
    )
    positives = int(data["y_true"].sum())
    prevalence = positives / len(data)
    cumulative = data["y_true"].cumsum().to_numpy()
    budget_rows = []
    for budget in BUDGETS:
        reviews = int(np.ceil(budget * len(data)))
        found = int(cumulative[reviews - 1])
        budget_rows.append(
            {
                "budget_fraction": budget,
                "review_rows": reviews,
                "pathogenic_found": found,
                "pathogenic_recall": found / positives,
                "precision": found / reviews,
                "enrichment_over_resolved_vus_prevalence": (found / reviews) / prevalence,
            }
        )
    target_rows = []
    for target in TARGETS:
        required = int(np.ceil(target * positives))
        reviews = int(np.flatnonzero(cumulative >= required)[0] + 1)
        random_expected = required / prevalence
        target_rows.append(
            {
                "target_pathogenic_recall": target,
                "required_pathogenic": required,
                "cardioqueue_reviews": reviews,
                "review_fraction": reviews / len(data),
                "oracle_minimum_reviews": required,
                "random_expected_reviews": random_expected,
                "relative_reduction_vs_random_expectation": 1 - reviews / random_expected,
            }
        )
    budget_table = pd.DataFrame(budget_rows)
    target_table = pd.DataFrame(target_rows)
    budget_table.to_csv(OUT / "temporal_vus_queue_budgets.tsv", sep="\t", index=False)
    target_table.to_csv(OUT / "temporal_vus_queue_targets.tsv", sep="\t", index=False)
    budget_table.to_csv(PAPER / "Table_S19_temporal_vus_queue_budgets.tsv", sep="\t", index=False)
    target_table.to_csv(PAPER / "Table_S19b_temporal_vus_queue_targets.tsv", sep="\t", index=False)
    summary = {
        "rows": len(data),
        "current_pathogenic": positives,
        "current_benign": len(data) - positives,
        "interpretation": "Retrospective label-temporal evaluation; predictor annotations were not uniformly frozen to January 2024.",
    }
    (OUT / "temporal_queue_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print("\n", budget_table.to_string(index=False))
    print("\n", target_table.to_string(index=False))


if __name__ == "__main__":
    main()
