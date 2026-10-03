#!/usr/bin/env python3
"""Paired gene-cluster bootstrap comparison of two expanded-external models."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


def score(y: np.ndarray, p: np.ndarray) -> dict[str, float]:
    return {
        "auroc": roc_auc_score(y, p),
        "auprc": average_precision_score(y, p),
        "brier": brier_score_loss(y, p),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--alternative", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--samples", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260925)
    args = parser.parse_args()

    use = ["variant_id", "y_true", "pathogenic_probability", "split_source_heldout"]
    base = pd.read_csv(args.baseline, sep="\t", usecols=use)
    alt = pd.read_csv(args.alternative, sep="\t", usecols=use)
    base = base[base["split_source_heldout"] == "external_expanded"]
    alt = alt[alt["split_source_heldout"] == "external_expanded"]
    merged = base.merge(alt, on=["variant_id", "y_true"], suffixes=("_baseline", "_alternative"), validate="one_to_one")
    manifest = pd.read_csv(args.manifest, sep="\t", usecols=["variant_id", "gene"])
    merged = merged.merge(manifest, on="variant_id", validate="one_to_one")
    y = merged["y_true"].to_numpy()
    baseline_point = score(y, merged["pathogenic_probability_baseline"].to_numpy())
    alternative_point = score(y, merged["pathogenic_probability_alternative"].to_numpy())

    rng = np.random.default_rng(args.seed)
    genes = merged["gene"].fillna("UNKNOWN").unique()
    by_gene = {gene: merged[merged["gene"].fillna("UNKNOWN") == gene] for gene in genes}
    deltas = {metric: [] for metric in baseline_point}
    for _ in range(args.samples):
        sampled = rng.choice(genes, size=len(genes), replace=True)
        boot = pd.concat([by_gene[gene] for gene in sampled], ignore_index=True)
        if boot["y_true"].nunique() < 2:
            continue
        yb = boot["y_true"].to_numpy()
        b = score(yb, boot["pathogenic_probability_baseline"].to_numpy())
        a = score(yb, boot["pathogenic_probability_alternative"].to_numpy())
        for metric in deltas:
            deltas[metric].append(a[metric] - b[metric])

    rows = []
    for metric, values in deltas.items():
        rows.append({
            "metric": metric,
            "baseline": baseline_point[metric],
            "alternative": alternative_point[metric],
            "alternative_minus_baseline": alternative_point[metric] - baseline_point[metric],
            "ci_95_low": np.percentile(values, 2.5),
            "ci_95_high": np.percentile(values, 97.5),
            "bootstrap_samples": args.samples,
            "bootstrap_unit": "gene",
        })
    output = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)
    print(output.to_string(index=False))


if __name__ == "__main__":
    main()
