#!/usr/bin/env python3

from pathlib import Path
import json

import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    matthews_corrcoef,
    precision_recall_curve,
    roc_auc_score,
)

OUTDIR = Path("results/model_performance/cardioboost_matched_benchmark")
TABLES = OUTDIR / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

OUR_PRED = Path(
    "results/models/cardioqueue_v1_source_heldout/"
    "primary_binary_predictions_split_source_heldout.tsv"
)
MATRIX = Path("datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv")
CB_CM = TABLES / "official_cardioboost_cm_all_rare_predictions.tsv"
CB_ARM = TABLES / "official_cardioboost_arm_all_rare_predictions.tsv"


def three_zone(score: pd.Series) -> pd.Series:
    return np.select(
        [score <= 0.1, score >= 0.9],
        ["Benign", "Pathogenic"],
        default="VUS",
    )


def parse_score_series(series: pd.Series) -> pd.Series:
    def parse_one(value):
        if pd.isna(value):
            return np.nan
        vals = []
        for token in str(value).replace(",", ";").split(";"):
            token = token.strip()
            if token in {"", ".", "nan", "NA", "None"}:
                continue
            try:
                vals.append(float(token))
            except ValueError:
                continue
        return max(vals) if vals else np.nan

    return series.apply(parse_one)


def metrics_for(df: pd.DataFrame, score_col: str, label_col: str = "y_true") -> dict:
    d = df[[label_col, score_col]].dropna()
    y = d[label_col].astype(int).to_numpy()
    s = pd.to_numeric(d[score_col], errors="coerce").to_numpy()
    ok = ~np.isnan(s)
    y = y[ok]
    s = s[ok]
    out = {
        "rows_scored": int(len(y)),
        "pathogenic": int(y.sum()),
        "benign": int((y == 0).sum()),
    }
    if len(np.unique(y)) == 2:
        out["AUROC"] = float(roc_auc_score(y, s))
        out["AUPRC"] = float(average_precision_score(y, s))
    else:
        out["AUROC"] = np.nan
        out["AUPRC"] = np.nan

    calls = np.asarray(three_zone(pd.Series(s)))
    pred_p = calls == "Pathogenic"
    pred_b = calls == "Benign"
    pred_v = calls == "VUS"
    tp = int(((y == 1) & pred_p).sum())
    fp = int(((y == 0) & pred_p).sum())
    tn = int(((y == 0) & pred_b).sum())
    fn = int(((y == 1) & pred_b).sum())
    p_defer = int(((y == 1) & pred_v).sum())
    b_defer = int(((y == 0) & pred_v).sum())
    out.update(
        {
            "TP_high_conf": tp,
            "FP_high_conf": fp,
            "TN_high_conf": tn,
            "FN_hard": fn,
            "P_to_defer": p_defer,
            "B_to_defer": b_defer,
            "pred_pathogenic": int(pred_p.sum()),
            "pred_benign": int(pred_b.sum()),
            "deferred": int(pred_v.sum()),
            "deferral_rate": float(pred_v.mean()) if len(y) else np.nan,
            "high_conf_classified_rate": float((pred_p | pred_b).mean()) if len(y) else np.nan,
            "sensitivity": float(tp / (tp + fn + p_defer)) if (tp + fn + p_defer) else np.nan,
            "specificity": float(tn / (tn + fp + b_defer)) if (tn + fp + b_defer) else np.nan,
            "PPV": float(tp / (tp + fp)) if (tp + fp) else np.nan,
            "NPV": float(tn / (tn + fn)) if (tn + fn) else np.nan,
            "high_conf_accuracy": float((tp + tn) / (pred_p | pred_b).sum())
            if (pred_p | pred_b).sum()
            else np.nan,
            "overall_threshold_accuracy": float((tp + tn) / len(y)) if len(y) else np.nan,
        }
    )
    binary_pred = (s >= 0.5).astype(int)
    out["MCC_at_0_5"] = float(matthews_corrcoef(y, binary_pred)) if len(np.unique(y)) == 2 else np.nan
    return out


def read_cardio(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df["cardioboost_score"] = pd.to_numeric(df["pathogenicity"], errors="coerce")
    df["cardioboost_cdot"] = df["HGVSc"].astype(str).str.split(":", n=1).str[-1]
    df["cardioboost_pdot"] = df["HGVSp"].astype(str).str.split(":", n=1).str[-1]
    return df[
        [
            "variant_id",
            "cardioboost_panel",
            "gene",
            "HGVSc",
            "HGVSp",
            "gnomAD_AF",
            "MCAP",
            "REVEL",
            "cardioboost_score",
            "cardioboost_cdot",
            "cardioboost_pdot",
        ]
    ].rename(
        columns={
            "gene": "cardioboost_gene",
            "HGVSc": "cardioboost_hgvsc",
            "HGVSp": "cardioboost_hgvsp",
            "gnomAD_AF": "cardioboost_gnomad_af",
            "MCAP": "mcap_score",
            "REVEL": "cardioboost_revel_score",
        }
    )


def main() -> None:
    our = pd.read_csv(OUR_PRED, sep="\t", low_memory=False)
    our = our.rename(columns={"pathogenic_probability": "our_catboost_score"})
    cols = [
        "variant_id",
        "primary_gene",
        "gnomad_final_af",
        "gnomad_final_popmax_af",
        "gnomad_final_status",
        "vep_worst_consequence",
        "vep_hgvsc",
        "vep_hgvsp",
        "vep_any_missense",
        "revel_score",
        "cadd_phred",
        "alphamissense_direct_score",
        "alphamissense_dbnsfp_score",
        "in_cardioboost",
        "in_clinvar",
        "in_hiro",
        "in_emerge",
        "split_source_heldout",
    ]
    matrix = pd.read_csv(MATRIX, sep="\t", usecols=lambda c: c in cols, low_memory=False)
    our = our.merge(matrix, on="variant_id", how="left")
    if "primary_gene_x" in our.columns:
        our["primary_gene"] = our["primary_gene_x"].fillna(our.get("primary_gene_y"))
        our = our.drop(columns=[c for c in ["primary_gene_x", "primary_gene_y"] if c in our.columns])
    for flag in ["in_cardioboost", "in_clinvar", "in_hiro", "in_emerge"]:
        x, y = f"{flag}_x", f"{flag}_y"
        if x in our.columns:
            our[flag] = our[x].fillna(our.get(y))
            our = our.drop(columns=[c for c in [x, y] if c in our.columns])
    if "split_source_heldout_x" in our.columns:
        our["split_source_heldout"] = our["split_source_heldout_x"].fillna(
            our.get("split_source_heldout_y")
        )
        our = our.drop(
            columns=[
                column
                for column in ["split_source_heldout_x", "split_source_heldout_y"]
                if column in our.columns
            ]
        )

    cb = pd.concat([read_cardio(CB_CM), read_cardio(CB_ARM)], ignore_index=True)
    cb = cb.sort_values(["variant_id", "cardioboost_panel"]).drop_duplicates(
        ["variant_id", "cardioboost_panel"]
    )
    matched = our.merge(cb, on="variant_id", how="inner")
    matched = matched[matched["model_label_3class"].isin(["Benign", "Pathogenic"])].copy()
    matched["y_true"] = (matched["model_label_3class"] == "Pathogenic").astype(int)
    matched["our_three_zone"] = three_zone(pd.to_numeric(matched["our_catboost_score"], errors="coerce"))
    matched["cardioboost_three_zone"] = three_zone(
        pd.to_numeric(matched["cardioboost_score"], errors="coerce")
    )

    # Prefer direct AlphaMissense, then dbNSFP AlphaMissense fallback.
    matched["alphamissense_best_score"] = pd.to_numeric(
        matched["alphamissense_direct_score"], errors="coerce"
    )
    matched["alphamissense_best_score"] = matched["alphamissense_best_score"].fillna(
        parse_score_series(matched["alphamissense_dbnsfp_score"])
    )
    matched["revel_score_numeric"] = parse_score_series(matched["revel_score"])

    matched.to_csv(TABLES / "exact_coordinate_cardioBoost_matched_benchmark.tsv", sep="\t", index=False)

    metric_rows = []
    score_cols = {
        "CardioQueue": "our_catboost_score",
        "Official CardioBoost": "cardioboost_score",
        "REVEL local": "revel_score_numeric",
        "REVEL CardioBoost-file": "cardioboost_revel_score",
        "CADD PHRED": "cadd_phred",
        "AlphaMissense best": "alphamissense_best_score",
        "M-CAP CardioBoost-file": "mcap_score",
    }
    for name, col in score_cols.items():
        if col not in matched.columns:
            continue
        row = {"model": name}
        row.update(metrics_for(matched, col))
        metric_rows.append(row)
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(TABLES / "exact_coordinate_cardioBoost_matched_metrics.tsv", sep="\t", index=False)

    summary = {
        "input_prediction_file": str(OUR_PRED),
        "input_matrix": str(MATRIX),
        "official_cardioboost_cm_predictions": str(CB_CM),
        "official_cardioboost_arm_predictions": str(CB_ARM),
        "exact_matched_binary_rows": int(len(matched)),
        "exact_matched_unique_variants": int(matched["variant_id"].nunique()),
        "exact_matched_panels": matched["cardioboost_panel"].value_counts(dropna=False).to_dict(),
    }
    with open(TABLES / "exact_coordinate_cardioBoost_matched_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    print(metrics.to_string(index=False))
    print(json.dumps(summary, indent=2))

    build_hgvs_benchmarks(our, cb)


def build_hgvs_benchmarks(our: pd.DataFrame, cb: pd.DataFrame) -> None:
    d = our.copy()
    d = d[d["model_label_3class"].isin(["Benign", "Pathogenic"])].copy()
    d["y_true"] = (d["model_label_3class"] == "Pathogenic").astype(int)
    d["our_cdot"] = d["vep_hgvsc"].astype(str).str.split(":", n=1).str[-1]
    d["our_pdot"] = d["vep_hgvsp"].astype(str).str.split(":", n=1).str[-1]
    d.loc[d["our_cdot"].isin(["nan", "", "-"]), "our_cdot"] = np.nan
    d.loc[d["our_pdot"].isin(["nan", "", "-"]), "our_pdot"] = np.nan

    ckey = ["primary_gene", "our_cdot"]
    cb_c = cb.rename(columns={"cardioboost_gene": "primary_gene", "cardioboost_cdot": "our_cdot"})
    cdna = d.merge(cb_c, on=ckey, how="inner", suffixes=("", "_cb"))
    cdna["match_method"] = "gene_cdot"

    # Protein-only fallback, excluding rows already matched by gene+cDot+panel.
    cb_p = cb.rename(columns={"cardioboost_gene": "primary_gene", "cardioboost_pdot": "our_pdot"})
    prot = d.merge(cb_p, on=["primary_gene", "our_pdot"], how="inner", suffixes=("", "_cb"))
    prot["match_method"] = "gene_pdot_fallback"
    existing = set(zip(cdna["variant_id"], cdna["cardioboost_panel"]))
    prot = prot[
        ~pd.Series(list(zip(prot["variant_id"], prot["cardioboost_panel"])), index=prot.index).isin(existing)
    ].copy()

    both = pd.concat([cdna, prot], ignore_index=True)
    both = both.sort_values(["variant_id", "cardioboost_panel", "match_method"])
    both = both.drop_duplicates(["variant_id", "cardioboost_panel", "match_method"])
    both["our_three_zone"] = three_zone(pd.to_numeric(both["our_catboost_score"], errors="coerce"))
    both["cardioboost_three_zone"] = three_zone(pd.to_numeric(both["cardioboost_score"], errors="coerce"))
    both["alphamissense_best_score"] = pd.to_numeric(
        both["alphamissense_direct_score"], errors="coerce"
    ).fillna(parse_score_series(both["alphamissense_dbnsfp_score"]))
    both["revel_score_numeric"] = parse_score_series(both["revel_score"])
    cdna["alphamissense_best_score"] = pd.to_numeric(
        cdna["alphamissense_direct_score"], errors="coerce"
    ).fillna(parse_score_series(cdna["alphamissense_dbnsfp_score"]))
    cdna["revel_score_numeric"] = parse_score_series(cdna["revel_score"])

    both.to_csv(TABLES / "hgvs_cardioBoost_matched_benchmark.tsv", sep="\t", index=False)
    cdna.to_csv(TABLES / "hgvs_cdot_cardioBoost_matched_benchmark.tsv", sep="\t", index=False)
    prot.to_csv(TABLES / "hgvs_pdot_fallback_cardioBoost_matched_benchmark.tsv", sep="\t", index=False)

    rows = []
    subset_defs = [
        ("hgvs_cdot_primary_all", cdna),
        ("hgvs_cdot_primary_validation", cdna[cdna["split_source_heldout"].eq("validation")]),
        (
            "hgvs_cdot_primary_external_all",
            cdna[cdna["split_source_heldout"].astype(str).str.startswith("external_")],
        ),
        (
            "hgvs_cdot_primary_external_cardioboost",
            cdna[cdna["split_source_heldout"].eq("external_cardioboost")],
        ),
        ("hgvs_cdot_plus_pdot_fallback_all", both),
        (
            "hgvs_cdot_plus_pdot_fallback_external_all",
            both[both["split_source_heldout"].astype(str).str.startswith("external_")],
        ),
        (
            "hgvs_cdot_plus_pdot_fallback_external_cardioboost",
            both[both["split_source_heldout"].eq("external_cardioboost")],
        ),
    ]
    for subset_name, subset in subset_defs:
        for name, col in {
            "CardioQueue": "our_catboost_score",
            "Official CardioBoost": "cardioboost_score",
            "REVEL local": "revel_score_numeric",
            "REVEL CardioBoost-file": "cardioboost_revel_score",
            "CADD PHRED": "cadd_phred",
            "AlphaMissense best": "alphamissense_best_score",
            "M-CAP CardioBoost-file": "mcap_score",
        }.items():
            if col not in subset.columns:
                continue
            row = {"benchmark": subset_name, "model": name}
            row.update(metrics_for(subset, col))
            rows.append(row)
    hgvs_metrics = pd.DataFrame(rows)
    hgvs_metrics.to_csv(TABLES / "hgvs_cardioBoost_matched_metrics.tsv", sep="\t", index=False)

    summary = {
        "cdna_panel_rows": int(len(cdna)),
        "cdna_unique_variants": int(cdna["variant_id"].nunique()) if len(cdna) else 0,
        "protein_fallback_panel_rows": int(len(prot)),
        "protein_fallback_unique_variants": int(prot["variant_id"].nunique()) if len(prot) else 0,
        "combined_panel_rows": int(len(both)),
        "combined_unique_variants": int(both["variant_id"].nunique()) if len(both) else 0,
        "cdna_panels": cdna["cardioboost_panel"].value_counts(dropna=False).to_dict()
        if len(cdna)
        else {},
        "combined_panels": both["cardioboost_panel"].value_counts(dropna=False).to_dict()
        if len(both)
        else {},
        "cdna_split_counts": cdna["split_source_heldout"].value_counts(dropna=False).to_dict()
        if len(cdna)
        else {},
        "combined_split_counts": both["split_source_heldout"].value_counts(dropna=False).to_dict()
        if len(both)
        else {},
    }
    with open(TABLES / "hgvs_cardioBoost_matched_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    print("\nHGVS matched metrics")
    print(hgvs_metrics.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
