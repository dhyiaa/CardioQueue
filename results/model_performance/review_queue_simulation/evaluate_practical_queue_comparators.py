#!/usr/bin/env python3
"""Compare CardioQueue review queues with routinely available annotation rankings."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[3]
PREDICTIONS = ROOT / "results/models/cardioqueue_v1_source_heldout/primary_binary_predictions_split_source_heldout.tsv"
TABLE = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
OUT = Path(__file__).resolve().parent
PAPER_TABLES = ROOT / "paper_MS_SI/tables"
PAPER_FIGURES = ROOT / "paper_MS_SI/figures"
SEED = 20260919
N_BOOT = 5000
N_TIE = 10000
BUDGETS = (0.10, 0.20, 0.30, 0.50)

CONSEQUENCE_RANK = {
    "transcript_ablation": 18, "splice_acceptor_variant": 17, "splice_donor_variant": 17,
    "stop_gained": 16, "frameshift_variant": 15, "stop_lost": 14, "start_lost": 14,
    "transcript_amplification": 13, "inframe_insertion": 12, "inframe_deletion": 12,
    "missense_variant": 11, "protein_altering_variant": 10, "splice_region_variant": 9,
    "incomplete_terminal_codon_variant": 8, "start_retained_variant": 7, "stop_retained_variant": 7,
    "synonymous_variant": 6, "coding_sequence_variant": 5, "mature_miRNA_variant": 5,
    "5_prime_UTR_variant": 4, "3_prime_UTR_variant": 4, "non_coding_transcript_exon_variant": 3,
    "intron_variant": 2, "NMD_transcript_variant": 2, "upstream_gene_variant": 1,
    "downstream_gene_variant": 1, "intergenic_variant": 0,
}


def parse_max(value: object) -> float:
    values = pd.to_numeric(pd.Series(str(value).split(";")), errors="coerce").dropna()
    return float(values.max()) if len(values) else np.nan


def ordered_indices(scores: np.ndarray, rng: np.random.Generator | None = None) -> np.ndarray:
    if rng is None:
        secondary = np.arange(len(scores))
    else:
        secondary = rng.random(len(scores))
    return np.lexsort((secondary, -scores))


def queue_metrics(y: np.ndarray, scores: np.ndarray, budget: float, rng: np.random.Generator | None = None) -> tuple[float, float]:
    order = ordered_indices(scores, rng)
    rows = max(1, int(np.ceil(budget * len(y))))
    found = int(y[order[:rows]].sum())
    return found / int(y.sum()), found / rows


def reviews_to_recall(y: np.ndarray, scores: np.ndarray, target: float, rng: np.random.Generator | None = None) -> int:
    order = ordered_indices(scores, rng)
    required = int(np.ceil(target * y.sum()))
    return int(np.flatnonzero(np.cumsum(y[order]) >= required)[0] + 1)


def load_frame() -> pd.DataFrame:
    pred = pd.read_csv(PREDICTIONS, sep="\t", low_memory=False)
    pred = pred[pred["split_source_heldout"].astype(str).str.startswith("external_")].copy()
    columns = [
        "variant_id", "primary_gene", "vep_worst_consequence", "revel_score_max",
        "alphamissense_direct_score", "alphamissense_dbnsfp_score_max", "dbnsfp_CADD_phred",
        "gnomad_final_status", "gnomad_final_popmax_af",
    ]
    table = pd.read_csv(TABLE, sep="\t", usecols=columns, low_memory=False)
    frame = pred.merge(table, on=["variant_id", "primary_gene"], validate="one_to_one")
    frame["cardioqueue"] = pd.to_numeric(frame["pathogenic_probability"], errors="coerce")
    frame["vep_consequence"] = frame["vep_worst_consequence"].map(CONSEQUENCE_RANK).fillna(-1).astype(float)
    frame["revel"] = pd.to_numeric(frame["revel_score_max"], errors="coerce")
    frame["alphamissense_direct"] = pd.to_numeric(frame["alphamissense_direct_score"], errors="coerce")
    fallback = pd.to_numeric(frame["alphamissense_dbnsfp_score_max"], errors="coerce")
    frame["alphamissense_composite"] = frame["alphamissense_direct"].fillna(fallback)
    frame["cadd"] = pd.to_numeric(frame["dbnsfp_CADD_phred"], errors="coerce")
    af = pd.to_numeric(frame["gnomad_final_popmax_af"], errors="coerce")
    valid_gnomad = frame["gnomad_final_status"].isin(["observed", "confirmed_absent"]) & af.notna()
    frame["gnomad_rarity"] = np.nan
    frame.loc[valid_gnomad, "gnomad_rarity"] = -np.log10(np.maximum(af[valid_gnomad], 1e-8))
    return frame


def comparator_specs(frame: pd.DataFrame) -> list[tuple[str, str, pd.Series]]:
    return [
        ("VEP consequence severity", "vep_consequence", pd.Series(True, index=frame.index)),
        ("CADD PHRED", "cadd", frame["cadd"].notna()),
        ("gnomAD rarity", "gnomad_rarity", frame["gnomad_rarity"].notna()),
        ("REVEL", "revel", frame["revel"].notna()),
        ("AlphaMissense direct", "alphamissense_direct", frame["alphamissense_direct"].notna()),
        ("AlphaMissense direct+fallback", "alphamissense_composite", frame["alphamissense_composite"].notna()),
    ]


def paired_bootstrap(data: pd.DataFrame, comparator_col: str) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    strata = [group.index.to_numpy() for _, group in data.groupby(["split_source_heldout", "y_true"], sort=True)]
    values = {budget: [] for budget in BUDGETS}
    y_full = data["y_true"].to_numpy(dtype=int)
    model_full = data["cardioqueue"].to_numpy(dtype=float)
    comparator_full = data[comparator_col].to_numpy(dtype=float)
    index_to_position = {index: position for position, index in enumerate(data.index)}
    for _ in range(N_BOOT):
        sampled_labels = np.concatenate([rng.choice(group, len(group), replace=True) for group in strata])
        positions = np.asarray([index_to_position[index] for index in sampled_labels])
        y = y_full[positions]
        model = model_full[positions]
        comparator = comparator_full[positions]
        for budget in BUDGETS:
            model_recall, _ = queue_metrics(y, model, budget)
            comparator_recall, _ = queue_metrics(y, comparator, budget, rng)
            values[budget].append(model_recall - comparator_recall)
    rows = []
    for budget, samples in values.items():
        array = np.asarray(samples)
        model_recall, _ = queue_metrics(y_full, model_full, budget)
        tie_recall = [queue_metrics(y_full, comparator_full, budget, rng)[0] for _ in range(N_TIE)]
        rows.append(
            {
                "budget_fraction": budget,
                "cardioqueue_recall": model_recall,
                "comparator_recall_tie_mean": float(np.mean(tie_recall)),
                "cardioqueue_minus_comparator_recall": model_recall - float(np.mean(tie_recall)),
                "paired_ci_95_low": float(np.quantile(array, 0.025)),
                "paired_ci_95_high": float(np.quantile(array, 0.975)),
                "bootstrap_samples": N_BOOT,
            }
        )
    return pd.DataFrame(rows)


def evaluate_comparator(frame: pd.DataFrame, name: str, column: str, mask: pd.Series) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    data = frame.loc[mask].copy()
    y = data["y_true"].to_numpy(dtype=int)
    model = data["cardioqueue"].to_numpy(dtype=float)
    comparator = data[column].to_numpy(dtype=float)
    rng = np.random.default_rng(SEED)
    coverage = {
        "comparator": name, "score_column": column, "rows": len(data),
        "coverage_fraction": len(data) / len(frame), "pathogenic": int(y.sum()),
        "benign": int((y == 0).sum()), "genes": data["primary_gene"].nunique(),
        "sources": ";".join(f"{key}:{value}" for key, value in data["split_source_heldout"].value_counts().sort_index().items()),
        "cardioqueue_auroc": roc_auc_score(y, model), "comparator_auroc": roc_auc_score(y, comparator),
        "cardioqueue_auprc": average_precision_score(y, model), "comparator_auprc": average_precision_score(y, comparator),
    }
    paired = paired_bootstrap(data, column)
    paired.insert(0, "comparator", name)
    target_model = reviews_to_recall(y, model, 0.80)
    target_comparator = np.asarray([reviews_to_recall(y, comparator, 0.80, rng) for _ in range(N_TIE)])
    target = pd.DataFrame([{
        "comparator": name, "rows": len(data), "pathogenic": int(y.sum()), "target_recall": 0.80,
        "cardioqueue_review_rows": target_model,
        "comparator_review_rows_tie_median": float(np.median(target_comparator)),
        "comparator_review_rows_tie_ci_95_low": float(np.quantile(target_comparator, 0.025)),
        "comparator_review_rows_tie_ci_95_high": float(np.quantile(target_comparator, 0.975)),
        "reviews_avoided_vs_comparator_median": float(np.median(target_comparator) - target_model),
    }])
    return coverage, paired, target


def plot_primary(differences: pd.DataFrame) -> None:
    data = differences[differences["budget_fraction"].eq(0.20)].copy()
    data = data.sort_values("cardioqueue_minus_comparator_recall")
    lower = data["cardioqueue_minus_comparator_recall"] - data["paired_ci_95_low"]
    upper = data["paired_ci_95_high"] - data["cardioqueue_minus_comparator_recall"]
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    colors = ["#0B6E75" if value >= 0 else "#B54736" for value in data["cardioqueue_minus_comparator_recall"]]
    ax.barh(data["comparator"], data["cardioqueue_minus_comparator_recall"], color=colors)
    ax.errorbar(data["cardioqueue_minus_comparator_recall"], data["comparator"], xerr=[lower, upper], fmt="none", color="#222222", capsize=3)
    ax.axvline(0, color="#555555", linewidth=1)
    ax.set_xlabel("Difference in P/LP recall at 20% review budget")
    ax.set_ylabel("")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#DDDDDD", linewidth=0.7)
    fig.tight_layout()
    for path in [OUT / "practical_queue_comparator_differences.png", PAPER_FIGURES / "Figure_6_practical_queue_comparators.png"]:
        fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    PAPER_TABLES.mkdir(parents=True, exist_ok=True)
    PAPER_FIGURES.mkdir(parents=True, exist_ok=True)
    frame = load_frame()
    if len(frame) != 430 or int(frame["y_true"].sum()) != 255:
        raise ValueError("Expected 430 external rows with 255 P/LP")
    coverage_rows, difference_frames, target_frames = [], [], []
    for name, column, mask in comparator_specs(frame):
        coverage, differences, target = evaluate_comparator(frame, name, column, mask)
        coverage_rows.append(coverage)
        difference_frames.append(differences)
        target_frames.append(target)
    coverage = pd.DataFrame(coverage_rows)
    differences = pd.concat(difference_frames, ignore_index=True)
    targets = pd.concat(target_frames, ignore_index=True)
    for data, local, paper in [
        (coverage, "practical_queue_comparator_coverage.tsv", "Table_S16e_practical_queue_comparator_coverage.tsv"),
        (differences, "practical_queue_comparator_differences.tsv", "Table_S16f_practical_queue_comparator_differences.tsv"),
        (targets, "practical_queue_comparator_targets.tsv", "Table_S16g_practical_queue_comparator_targets.tsv"),
    ]:
        data.to_csv(OUT / local, sep="\t", index=False)
        data.to_csv(PAPER_TABLES / paper, sep="\t", index=False)
    plot_primary(differences)
    manifest = {
        "prediction_file": str(PREDICTIONS.relative_to(ROOT)), "modeling_table": str(TABLE.relative_to(ROOT)),
        "external_rows": len(frame), "bootstrap_samples": N_BOOT, "tie_permutations": N_TIE, "seed": SEED,
        "primary_endpoint": "Paired difference in P/LP recall at a 20% review budget on each comparator-specific shared cohort.",
        "missing_score_policy": "Comparator-specific complete cases; missing scores were never imputed as benign.",
    }
    (OUT / "practical_queue_comparator_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(coverage.to_string(index=False))
    print(differences[differences["budget_fraction"].eq(0.20)].to_string(index=False))
    print(targets.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
