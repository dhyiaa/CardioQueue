#!/usr/bin/env python3
"""Queue ceiling, prevalence transport, and rare-missense sensitivity analyses."""

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
PAPER = ROOT / "paper_MS_SI/tables"
PRED = (
    ROOT
    / "results/models/cardioqueue_v1_source_heldout"
    / "primary_binary_predictions_split_source_heldout.tsv"
)
MODEL_TABLE = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
TARGETS = (0.50, 0.80, 0.90, 0.95)
PREVALENCES = (0.02, 0.05, 0.10, 0.20)


def order(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.sort_values(
        ["pathogenic_probability", "variant_id"], ascending=[False, True], kind="mergesort"
    )


def operating_point(frame: pd.DataFrame, target: float) -> dict:
    ranked = order(frame)
    required = int(np.ceil(target * ranked["y_true"].sum()))
    reviewed = int(np.flatnonzero(ranked["y_true"].cumsum().to_numpy() >= required)[0] + 1)
    selected = ranked.iloc[:reviewed]
    tp = int(selected["y_true"].sum())
    fp = reviewed - tp
    positives = int(ranked["y_true"].sum())
    negatives = len(ranked) - positives
    return {
        "target_recall": target,
        "required_pathogenic": required,
        "review_rows": reviewed,
        "review_fraction": reviewed / len(ranked),
        "sensitivity": tp / positives,
        "false_positive_rate": fp / negatives,
        "precision": tp / reviewed,
    }


def ceiling_table(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    prevalence = frame["y_true"].mean()
    for target in TARGETS:
        point = operating_point(frame, target)
        oracle = point["required_pathogenic"]
        random_expected = point["required_pathogenic"] / prevalence
        denominator = random_expected - oracle
        attained = (random_expected - point["review_rows"]) / denominator if denominator > 0 else np.nan
        rows.append(
            {
                **point,
                "oracle_minimum_reviews": oracle,
                "excess_reviews_above_oracle": point["review_rows"] - oracle,
                "random_expected_reviews": random_expected,
                "fraction_of_attainable_reduction": attained,
            }
        )
    return pd.DataFrame(rows)


def transport_table(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for target in (0.80, 0.90, 0.95):
        point = operating_point(frame, target)
        sensitivity = point["sensitivity"]
        fpr = point["false_positive_rate"]
        for prevalence in PREVALENCES:
            selected_fraction = sensitivity * prevalence + fpr * (1 - prevalence)
            ppv = sensitivity * prevalence / selected_fraction if selected_fraction else np.nan
            rows.append(
                {
                    "target_recall": target,
                    "heldout_threshold_sensitivity": sensitivity,
                    "heldout_threshold_fpr": fpr,
                    "assumed_pathogenic_prevalence": prevalence,
                    "expected_review_fraction": selected_fraction,
                    "expected_reviews_per_1000": 1000 * selected_fraction,
                    "expected_pathogenic_found_per_1000": 1000 * sensitivity * prevalence,
                    "transported_ppv": ppv,
                    "interpretation": "Algebraic prevalence transport with held-out sensitivity and FPR fixed; case mix and calibration are assumed unchanged.",
                }
            )
    return pd.DataFrame(rows)


def rare_missense_table(frame: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "variant_id",
        "vep_any_missense",
        "gnomad_final_status",
        "gnomad_final_af",
    ]
    annotations = pd.read_csv(MODEL_TABLE, sep="\t", usecols=columns, low_memory=False)
    merged = frame.merge(annotations, on="variant_id", how="left", validate="one_to_one")
    missense = merged["vep_any_missense"].astype(str).str.lower().eq("true")
    af = pd.to_numeric(merged["gnomad_final_af"], errors="coerce")
    absent = merged["gnomad_final_status"].astype(str).str.contains("absent", case=False, na=False)
    rare = merged.loc[missense & (absent | af.lt(1e-4))].copy()
    rows = []
    for name, subset in [("all_external", merged), ("rare_missense_af_lt_1e-4_or_absent", rare)]:
        for target in TARGETS:
            rows.append({"stratum": name, "rows": len(subset), "pathogenic": int(subset.y_true.sum()), **operating_point(subset, target)})
    return pd.DataFrame(rows)


def main() -> None:
    frame = pd.read_csv(PRED, sep="\t")
    frame = frame[frame["split_source_heldout"].astype(str).str.startswith("external_")].copy()
    outputs = {
        "queue_oracle_ceiling.tsv": ceiling_table(frame),
        "queue_prevalence_transport.tsv": transport_table(frame),
        "queue_rare_missense_sensitivity.tsv": rare_missense_table(frame),
    }
    paper_names = {
        "queue_oracle_ceiling.tsv": "Table_S16h_queue_oracle_ceiling.tsv",
        "queue_prevalence_transport.tsv": "Table_S16i_queue_prevalence_transport.tsv",
        "queue_rare_missense_sensitivity.tsv": "Table_S16j_queue_rare_missense_sensitivity.tsv",
    }
    for name, table in outputs.items():
        table.to_csv(OUT / name, sep="\t", index=False)
        table.to_csv(PAPER / paper_names[name], sep="\t", index=False)
        print(f"\n{name}\n{table.to_string(index=False)}")


if __name__ == "__main__":
    main()
