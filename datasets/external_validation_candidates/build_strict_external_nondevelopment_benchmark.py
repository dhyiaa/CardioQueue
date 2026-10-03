#!/usr/bin/env python3
"""Build and evaluate the final multi-source model-held-out benchmark."""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[2]
EXPANDED = ROOT / "datasets/external_validation_candidates/expanded_external_2026"
ORIGINAL_MATRIX = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
PREDICTIONS = ROOT / "results/model_performance/expanded_external_2026/combined_unique_predictions.tsv"
OUT_DATA = ROOT / "datasets/external_validation_candidates/final_model_heldout_2026"
OUT_RESULTS = ROOT / "results/model_performance/final_model_heldout_2026"

PRIOR_MODEL_HELDOUT = {"legacy_hiro", "legacy_emerge", "legacy_cardioboost"}
INSTITUTIONALLY_EXTERNAL = {"legacy_emerge", "legacy_cardioboost"}
ADDED_EXTERNAL = {"Fernandez_Falgueras_PLOS_2024", "SHaRe_HCM_2026Q1"}
DEVELOPMENT_SPLITS = {"train", "validation"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    output = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + 1 + end) / 2
        for index in order[start:end]:
            output[index] = rank
        start = end
    return output


def auroc(labels: list[int], scores: list[float]) -> float:
    positive = sum(labels)
    negative = len(labels) - positive
    if not positive or not negative:
        return math.nan
    score_ranks = ranks(scores)
    positive_rank_sum = sum(rank for rank, label in zip(score_ranks, labels) if label)
    return (positive_rank_sum - positive * (positive + 1) / 2) / (positive * negative)


def average_precision(labels: list[int], scores: list[float]) -> float:
    positive = sum(labels)
    if not positive:
        return math.nan
    ordered = sorted(zip(scores, labels), reverse=True)
    found = reviewed = 0
    output = 0.0
    start = 0
    while start < len(ordered):
        end = start + 1
        while end < len(ordered) and ordered[end][0] == ordered[start][0]:
            end += 1
        group_positive = sum(label for _, label in ordered[start:end])
        found += group_positive
        reviewed += end - start
        output += (group_positive / positive) * (found / reviewed)
        start = end
    return output


def metrics(rows: list[dict]) -> dict[str, float | int]:
    labels = [int(row["y_true"]) for row in rows]
    scores = [float(row["pathogenic_probability"]) for row in rows]
    tp = fp = tn = fn = 0
    for label, score in zip(labels, scores):
        predicted = int(score >= 0.5)
        if label and predicted:
            tp += 1
        elif label:
            fn += 1
        elif predicted:
            fp += 1
        else:
            tn += 1
    return {
        "rows": len(rows),
        "pathogenic": sum(labels),
        "benign": len(labels) - sum(labels),
        "auroc": auroc(labels, scores),
        "auprc": average_precision(labels, scores),
        "brier": sum((score - label) ** 2 for score, label in zip(scores, labels)) / len(rows),
        "sensitivity": tp / (tp + fn),
        "specificity": tn / (tn + fp),
        "ppv": tp / (tp + fp),
        "npv": tn / (tn + fn),
        "tp": tp,
        "fp": fp,
        "tn": tn,
        "fn": fn,
    }


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] * (upper - position) + ordered[upper] * (position - lower)


def gene_bootstrap(rows: list[dict], samples: int, seed: int) -> list[dict]:
    generator = random.Random(seed)
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[row.get("gene") or "UNKNOWN"].append(row)
    genes = sorted(grouped)
    values: dict[str, list[float]] = defaultdict(list)
    for _ in range(samples):
        sample = []
        for gene in generator.choices(genes, k=len(genes)):
            sample.extend(grouped[gene])
        if len({row["y_true"] for row in sample}) < 2:
            continue
        point = metrics(sample)
        for name in ("auroc", "auprc", "brier", "sensitivity", "specificity", "ppv", "npv"):
            value = float(point[name])
            if not math.isnan(value):
                values[name].append(value)
    point = metrics(rows)
    return [
        {
            "metric": name,
            "estimate": point[name],
            "ci_95_low": percentile(metric_values, 0.025),
            "ci_95_high": percentile(metric_values, 0.975),
            "bootstrap_samples": samples,
            "bootstrap_unit": "gene",
        }
        for name, metric_values in values.items()
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original-matrix", type=Path, default=ORIGINAL_MATRIX)
    parser.add_argument("--expanded-manifest", type=Path, default=EXPANDED / "external_unique_variant_manifest.tsv")
    parser.add_argument("--expanded-observations", type=Path, default=EXPANDED / "external_source_observations.tsv")
    parser.add_argument("--expanded-conflicts", type=Path, default=EXPANDED / "external_label_conflicts.tsv")
    parser.add_argument("--predictions", type=Path, default=PREDICTIONS)
    parser.add_argument("--out-data", type=Path, default=OUT_DATA)
    parser.add_argument("--out-results", type=Path, default=OUT_RESULTS)
    parser.add_argument("--bootstrap-samples", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260925)
    args = parser.parse_args()

    original = {
        row["variant_id"]: row.get("split_source_heldout", "")
        for row in read_tsv(args.original_matrix)
    }
    expanded_manifest = {row["variant_id"]: row for row in read_tsv(args.expanded_manifest)}
    observations = read_tsv(args.expanded_observations)
    conflict_rows = read_tsv(args.expanded_conflicts)
    predictions = {row["variant_id"]: row for row in read_tsv(args.predictions)}

    eligible_observations = []
    audit_rows = []
    for row in observations:
        variant_id = row["variant_id"]
        source = row["source"]
        if variant_id not in expanded_manifest:
            continue
        prior_split = original.get(variant_id, "not_in_original_matrix")
        include = source in PRIOR_MODEL_HELDOUT or (
            source in ADDED_EXTERNAL and prior_split not in DEVELOPMENT_SPLITS
        )
        if source in PRIOR_MODEL_HELDOUT:
            reason = "included_prior_model_heldout_source"
        elif source in ADDED_EXTERNAL and include:
            reason = "included_absent_from_prior_development"
        elif source in ADDED_EXTERNAL:
            reason = "excluded_prior_train_or_validation"
        else:
            reason = "excluded_internal_source"
        audit_rows.append({
            "source": source,
            "variant_id": variant_id,
            "prior_split": prior_split,
            "included": str(include),
            "reason": reason,
        })
        if include:
            eligible_observations.append({**row, "prior_split": prior_split})

    benchmark_ids = {row["variant_id"] for row in eligible_observations}
    manifest_rows = []
    for variant_id in sorted(benchmark_ids):
        base = expanded_manifest[variant_id]
        sources = sorted({row["source"] for row in eligible_observations if row["variant_id"] == variant_id})
        manifest_rows.append({
            **base,
            "eligible_external_sources": ";".join(sources),
            "prior_split": original.get(variant_id, "not_in_original_matrix"),
        })

    if benchmark_ids - predictions.keys():
        raise ValueError("Eligible benchmark variants are missing predictions")
    if any(original.get(variant_id) in DEVELOPMENT_SPLITS for variant_id in benchmark_ids if not any(
        row["source"] in PRIOR_MODEL_HELDOUT for row in eligible_observations if row["variant_id"] == variant_id
    )):
        raise ValueError("A PLOS/SHaRE-only benchmark allele overlaps prior development")

    combined = []
    for row in manifest_rows:
        prediction = predictions[row["variant_id"]]
        combined.append({
            **row,
            "y_true": 1 if row["combined_label"] == "Pathogenic" else 0,
            "pathogenic_probability": prediction["pathogenic_probability"],
            "analysis_subset": "combined_unique",
        })
    source_predictions = []
    for row in eligible_observations:
        source_predictions.append({
            **row,
            "y_true": row["external_binary_label"],
            "pathogenic_probability": predictions[row["variant_id"]]["pathogenic_probability"],
            "analysis_subset": row["source"],
        })

    subsets = [("combined_unique", combined)]
    for source in sorted(PRIOR_MODEL_HELDOUT | ADDED_EXTERNAL):
        subset = [row for row in source_predictions if row["source"] == source]
        if subset:
            subsets.append((source, subset))
    consensus = [
        {**row, "analysis_subset": "multi_source_consensus"}
        for row in combined
        if len(row["eligible_external_sources"].split(";")) > 1
    ]
    subsets.append(("multi_source_consensus", consensus))
    performance = []
    intervals = []
    for offset, (name, rows) in enumerate(subsets):
        performance.append({"subset": name, **metrics(rows)})
        for interval in gene_bootstrap(rows, args.bootstrap_samples, args.seed + offset):
            intervals.append({"subset": name, **interval})

    args.out_data.mkdir(parents=True, exist_ok=True)
    args.out_results.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out_data / "external_unique_variant_manifest.tsv", manifest_rows, manifest_rows[0].keys())
    write_tsv(args.out_data / "external_source_observations.tsv", eligible_observations, eligible_observations[0].keys())
    write_tsv(args.out_data / "external_label_conflicts.tsv", conflict_rows, conflict_rows[0].keys())
    write_tsv(args.out_data / "eligibility_audit.tsv", audit_rows, audit_rows[0].keys())
    write_tsv(args.out_results / "combined_unique_predictions.tsv", combined, combined[0].keys())
    write_tsv(args.out_results / "source_observation_predictions.tsv", source_predictions, source_predictions[0].keys())
    write_tsv(args.out_results / "multi_source_consensus_predictions.tsv", consensus, consensus[0].keys())
    write_tsv(args.out_results / "performance_by_source.tsv", performance, performance[0].keys())
    write_tsv(args.out_results / "gene_cluster_bootstrap_ci.tsv", intervals, intervals[0].keys())

    source_sets: dict[str, set[str]] = defaultdict(set)
    for row in eligible_observations:
        source_sets[row["variant_id"]].add(row["source"])
    overlap_combinations = []
    for combination in sorted({tuple(sorted(sources)) for sources in source_sets.values()}):
        variant_ids = {
            variant_id
            for variant_id, sources in source_sets.items()
            if tuple(sorted(sources)) == combination
        }
        labels = [row["combined_label"] for row in manifest_rows if row["variant_id"] in variant_ids]
        overlap_combinations.append({
            "source_combination": ";".join(combination),
            "source_count": len(combination),
            "unique_alleles": len(variant_ids),
            "pathogenic": labels.count("Pathogenic"),
            "benign": labels.count("Benign"),
        })
    write_tsv(
        args.out_data / "source_overlap_combinations.tsv",
        overlap_combinations,
        overlap_combinations[0].keys(),
    )

    source_display = {
        "legacy_hiro": "HiRO",
        "legacy_emerge": "eMERGE",
        "legacy_cardioboost": "CardioBoost_strict",
        "SHaRe_HCM_2026Q1": "SHaRe_eligible",
        "Fernandez_Falgueras_PLOS_2024": "PLOS_eligible",
    }
    source_status = {
        "legacy_hiro": "model-external; institutionally internal",
        "legacy_emerge": "model-external; institutionally external",
        "legacy_cardioboost": "model-external; institutionally external",
        "SHaRe_HCM_2026Q1": "model-external; institutionally external; retained alleles absent from original train/validation",
        "Fernandez_Falgueras_PLOS_2024": "model-external; institutionally external; retained alleles absent from original train/validation",
    }
    accounting = []
    ordered_sources = [
        "legacy_hiro",
        "legacy_emerge",
        "legacy_cardioboost",
        "SHaRe_HCM_2026Q1",
        "Fernandez_Falgueras_PLOS_2024",
    ]
    for source in ordered_sources:
        candidates = [row for row in audit_rows if row["source"] == source]
        retained = [row for row in eligible_observations if row["source"] == source]
        retained_ids = {row["variant_id"] for row in retained}
        accounting.append({
            "source": source_display[source],
            "strict_candidates_before_prior_development_filter": len(candidates),
            "excluded_prior_train_or_validation": sum(
                row["reason"] == "excluded_prior_train_or_validation" for row in candidates
            ),
            "retained_source_assertions": len(retained),
            "pathogenic": sum(int(row["external_binary_label"]) for row in retained),
            "benign": sum(not int(row["external_binary_label"]) for row in retained),
            "alleles_also_in_another_retained_source": sum(
                len(source_sets[variant_id]) > 1 for variant_id in retained_ids
            ),
            "source_exclusive_alleles": sum(
                len(source_sets[variant_id]) == 1 for variant_id in retained_ids
            ),
            "external_status": source_status[source],
        })
    accounting.append({
        "source": "combined_source_assertions",
        "strict_candidates_before_prior_development_filter": sum(
            int(row["strict_candidates_before_prior_development_filter"]) for row in accounting
        ),
        "excluded_prior_train_or_validation": sum(
            int(row["excluded_prior_train_or_validation"]) for row in accounting
        ),
        "retained_source_assertions": len(eligible_observations),
        "pathogenic": sum(int(row["external_binary_label"]) for row in eligible_observations),
        "benign": sum(not int(row["external_binary_label"]) for row in eligible_observations),
        "alleles_also_in_another_retained_source": sum(
            len(sources) for sources in source_sets.values() if len(sources) > 1
        ),
        "source_exclusive_alleles": sum(len(sources) == 1 for sources in source_sets.values()),
        "external_status": "source counts overlap",
    })
    accounting.append({
        "source": "combined_unique_alleles",
        "strict_candidates_before_prior_development_filter": "NA",
        "excluded_prior_train_or_validation": "NA",
        "retained_source_assertions": len(combined),
        "pathogenic": sum(int(row["y_true"]) for row in combined),
        "benign": sum(not int(row["y_true"]) for row in combined),
        "alleles_also_in_another_retained_source": len(consensus),
        "source_exclusive_alleles": sum(len(sources) == 1 for sources in source_sets.values()),
        "external_status": "one row per exact allele for pooled performance",
    })
    write_tsv(args.out_data / "source_accounting.tsv", accounting, accounting[0].keys())

    source_counts = {}
    for source in sorted(PRIOR_MODEL_HELDOUT | ADDED_EXTERNAL):
        rows = [row for row in eligible_observations if row["source"] == source]
        if rows:
            source_counts[source] = {
                "observations": len(rows),
                "pathogenic": sum(int(row["external_binary_label"]) for row in rows),
                "benign": sum(not int(row["external_binary_label"]) for row in rows),
            }
    summary = {
        "design": "HiRO, eMERGE, and CardioBoost plus PLOS/SHaRE alleles absent from original train and validation",
        "unique_variants": len(combined),
        "pathogenic": sum(int(row["y_true"]) for row in combined),
        "benign": sum(not int(row["y_true"]) for row in combined),
        "source_counts": source_counts,
        "combined_performance": performance[0],
        "multi_source_consensus_performance": next(
            row for row in performance if row["subset"] == "multi_source_consensus"
        ),
        "cross_source_label_conflicts": len(conflict_rows),
        "bootstrap": {"unit": "gene", "samples": args.bootstrap_samples, "seed": args.seed},
        "source_interpretation": {
            "legacy_hiro": "model-external/source-held-out but institutionally internal",
            "legacy_emerge": "model-external and institutionally external",
            "legacy_cardioboost": "model-external and institutionally external",
            "Fernandez_Falgueras_PLOS_2024": "model-external and institutionally external; study-team ACMG/AMP reinterpretation distinct from the current ClinVar aggregate",
            "SHaRe_HCM_2026Q1": "model-external and institutionally external; SHaRE classification retained separately from ClinVar",
        },
        "important_limitation": (
            "Prior-nondevelopment status is not equivalent to complete absence from ClinVar or public evidence."
        ),
    }
    (args.out_data / "build_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (args.out_results / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
