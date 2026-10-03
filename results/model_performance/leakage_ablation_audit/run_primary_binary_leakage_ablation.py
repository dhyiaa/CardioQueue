#!/usr/bin/env python3
"""Run leakage-oriented feature-family ablations for the primary binary model."""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from catboost import CatBoostClassifier, Pool
from sklearn.metrics import average_precision_score, confusion_matrix, roc_auc_score


ROOT = Path(__file__).resolve().parents[3]
TABLE = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
TRAIN_SCRIPT = ROOT / "datasets/modeling/scripts/train_primary_binary_catboost.py"
OUT = ROOT / "results/model_performance/leakage_ablation_audit"
SPLIT_COL = "split_source_heldout"
SEED = 20260703


def load_primary_trainer():
    spec = importlib.util.spec_from_file_location("primary_binary_trainer", TRAIN_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {TRAIN_SCRIPT}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


TRAINER = load_primary_trainer()


@dataclass(frozen=True)
class Ablation:
    name: str
    description: str
    patterns: tuple[str, ...] = ()
    exact: tuple[str, ...] = ()


GENE_PATTERNS = (
    r"(^|_)gene($|_)",
    r"genename",
    r"gene_name",
    r"gene_symbol",
    r"primary_gene",
    r"symbol$",
    r"gene_panel",
    r"clingen_gene_validity",
)

CONSEQUENCE_PATTERNS = (
    r"^vep_(annotation_status|transcript_count|symbol|gene_id|feature|feature_type)$",
    r"^vep_(consequence|worst_consequence|impact|canonical|mane_select|mane_plus_clinical)$",
    r"^vep_(hgvsc|hgvsp|protein_position|amino_acids|codons|all_consequences|all_impacts)$",
    r"^vep_any_",
    r"gnomad_browser_.*(major_consequence|hgvsc|hgvsp|transcript_id|gene_symbol|consequence_terms)",
)

CLINVAR_PATTERNS = (
    r"clinvar",
    r"review",
    r"star",
    r"confidence",
    r"clinical_significance",
    r"expert",
)

SOURCE_AGG_PATTERNS = (
    r"^hiro_",
    r"^emerge_",
    r"^cardioboost_",
    r"source_record",
    r"source_row",
    r"selected_source",
    r"selected_variant",
    r"source_dataset",
    r"source_uid",
    r"source_count",
    r"source_label",
    r"source_confidence",
    r"phenotype",
    r"patient",
)

HGVS_ID_PATTERNS = (
    r"hgv",
    r"transcript",
    r"feature$",
    r"uniprot_accession",
    r"protein_variant",
    r"existing_variation",
    r"rsid",
    r"caid",
    r"coordinate",
    r"variant_id",
)

PRECOMPUTED_EFFECT_PREDICTOR_PATTERNS = (
    r"sift",
    r"polyphen",
    r"revel",
    r"metalr",
    r"fathmm",
    r"cadd",
    r"alphamissense",
    r"esm1b",
)

ABLATIONS = (
    Ablation(
        "baseline_current",
        "Current primary feature-selection rules.",
    ),
    Ablation(
        "no_gene_identity",
        "Remove direct and gene-level identity features, including gene symbols and gene validity flags.",
        patterns=GENE_PATTERNS,
        exact=("primary_gene", "gene_panel_decisions"),
    ),
    Ablation(
        "no_vep_consequence_hgvs",
        "Remove VEP/gnomAD consequence, impact, transcript, and HGVS-like annotations while keeping SpliceAI scores.",
        patterns=CONSEQUENCE_PATTERNS,
    ),
    Ablation(
        "no_clinvar_metadata",
        "Remove any residual ClinVar/review/confidence metadata beyond the default leakage filter.",
        patterns=CLINVAR_PATTERNS,
    ),
    Ablation(
        "no_source_record_aggregates",
        "Remove source-record aggregates and HiRO phenotype/source evidence fields.",
        patterns=SOURCE_AGG_PATTERNS,
    ),
    Ablation(
        "no_hgvs_transcript_identifiers",
        "Remove HGVS, transcript, protein-accession, rsID/CAID, and direct identifier-like fields.",
        patterns=HGVS_ID_PATTERNS,
    ),
    Ablation(
        "strict_low_leakage",
        "Remove gene identity, consequence/HGVS/transcript identifiers, ClinVar metadata, and source-record aggregates together.",
        patterns=GENE_PATTERNS + CONSEQUENCE_PATTERNS + CLINVAR_PATTERNS + SOURCE_AGG_PATTERNS + HGVS_ID_PATTERNS,
        exact=("primary_gene", "gene_panel_decisions"),
    ),
    Ablation(
        "no_precomputed_effect_predictors",
        "Remove SIFT, PolyPhen-2, REVEL, MetaLR, FATHMM-XF, CADD, AlphaMissense, and ESM-1b fields.",
        patterns=PRECOMPUTED_EFFECT_PREDICTOR_PATTERNS,
    ),
)


def ablation_drops(col: str, ablation: Ablation) -> bool:
    if col in ablation.exact:
        return True
    return any(re.search(pattern, col, flags=re.IGNORECASE) for pattern in ablation.patterns)


def feature_matrix(df: pd.DataFrame, ablation: Ablation) -> tuple[pd.DataFrame, list[str], list[str], list[str], list[str]]:
    base_features = [c for c in df.columns if not TRAINER.should_drop(c)]
    feature_cols = [c for c in base_features if not ablation_drops(c, ablation)]
    dropped = [c for c in base_features if c not in feature_cols]
    x = df[feature_cols].copy()

    cat_cols: list[str] = []
    numeric_cols: list[str] = []
    for col in x.columns:
        if x[col].dtype == "object" or str(x[col].dtype).startswith("string") or x[col].dtype == bool:
            cat_cols.append(col)
            x[col] = x[col].astype("string").fillna("__MISSING__").astype(str)
        else:
            numeric_cols.append(col)
            x[col] = pd.to_numeric(x[col], errors="coerce")
    return x, feature_cols, cat_cols, numeric_cols, dropped


def metrics_at_threshold(y_true: np.ndarray, proba: np.ndarray, threshold: float = 0.5) -> dict[str, float | int]:
    pred = (proba >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, pred, labels=[0, 1]).ravel()
    return {
        "tp": int(tp),
        "fp": int(fp),
        "tn": int(tn),
        "fn": int(fn),
        "sensitivity": float(tp / (tp + fn)) if (tp + fn) else np.nan,
        "specificity": float(tn / (tn + fp)) if (tn + fp) else np.nan,
        "ppv": float(tp / (tp + fp)) if (tp + fp) else np.nan,
        "npv": float(tn / (tn + fn)) if (tn + fn) else np.nan,
    }


def metrics_three_zone(y_true: np.ndarray, proba: np.ndarray) -> dict[str, float | int]:
    pred_zone = np.full(len(proba), "defer", dtype=object)
    pred_zone[proba <= 0.1] = "benign"
    pred_zone[proba >= 0.9] = "pathogenic"
    called = pred_zone != "defer"
    tp = int(((y_true == 1) & (pred_zone == "pathogenic")).sum())
    fp = int(((y_true == 0) & (pred_zone == "pathogenic")).sum())
    tn = int(((y_true == 0) & (pred_zone == "benign")).sum())
    fn = int(((y_true == 1) & (pred_zone == "benign")).sum())
    return {
        "defer_count": int((~called).sum()),
        "deferral_rate": float((~called).mean()),
        "three_zone_tp": tp,
        "three_zone_fp": fp,
        "three_zone_tn": tn,
        "three_zone_fn": fn,
        "three_zone_sensitivity": float(tp / (y_true == 1).sum()) if (y_true == 1).sum() else np.nan,
        "three_zone_specificity": float(tn / (y_true == 0).sum()) if (y_true == 0).sum() else np.nan,
        "three_zone_ppv": float(tp / (tp + fp)) if (tp + fp) else np.nan,
        "three_zone_high_conf_accuracy": float((tp + tn) / called.sum()) if called.sum() else np.nan,
    }


def evaluate_subset(y_true: np.ndarray, proba: np.ndarray) -> dict[str, float | int]:
    out = {
        "rows": int(len(y_true)),
        "positives_pathogenic": int(y_true.sum()),
        "negatives_benign": int((y_true == 0).sum()),
        "auroc": float(roc_auc_score(y_true, proba)) if len(np.unique(y_true)) == 2 else np.nan,
        "auprc": float(average_precision_score(y_true, proba)) if len(np.unique(y_true)) == 2 else np.nan,
    }
    out.update({f"threshold_0_5_{k}": v for k, v in metrics_at_threshold(y_true, proba).items()})
    out.update(metrics_three_zone(y_true, proba))
    return out


def train_one(df: pd.DataFrame, ablation: Ablation) -> tuple[pd.DataFrame, dict]:
    x, feature_cols, cat_cols, numeric_cols, dropped_cols = feature_matrix(df, ablation)
    y = df["primary_binary_label"].astype(int).to_numpy()
    weights = df["sample_weight"].astype(float).to_numpy()
    split = df[SPLIT_COL].astype(str)

    train_mask = split.eq("train").to_numpy()
    val_mask = split.eq("validation").to_numpy()
    train_pool = Pool(x.loc[train_mask], y[train_mask], weight=weights[train_mask], cat_features=cat_cols)
    val_pool = Pool(x.loc[val_mask], y[val_mask], weight=weights[val_mask], cat_features=cat_cols)

    model = CatBoostClassifier(
        iterations=1500,
        learning_rate=0.035,
        depth=6,
        loss_function="Logloss",
        eval_metric="PRAUC",
        random_seed=SEED,
        l2_leaf_reg=6,
        od_type="Iter",
        od_wait=80,
        allow_writing_files=False,
        verbose=False,
    )
    model.fit(train_pool, eval_set=val_pool, use_best_model=True)

    rows = []
    metric_rows = []
    split_names = ["train", "validation", "external_hiro", "external_emerge", "external_cardioboost"]
    for split_name in split_names:
        mask = split.eq(split_name).to_numpy()
        if not mask.any():
            continue
        proba = model.predict_proba(Pool(x.loc[mask], cat_features=cat_cols))[:, 1]
        y_true = y[mask]
        m = evaluate_subset(y_true, proba)
        m.update({"ablation": ablation.name, "subset": split_name})
        metric_rows.append(m)
        tmp = df.loc[mask, ["variant_id", "primary_gene", "model_label_3class", SPLIT_COL]].copy()
        tmp["ablation"] = ablation.name
        tmp["y_true"] = y_true
        tmp["pathogenic_probability"] = proba
        rows.append(tmp)

    external_mask = split.isin(["external_hiro", "external_emerge", "external_cardioboost"]).to_numpy()
    if external_mask.any():
        proba = model.predict_proba(Pool(x.loc[external_mask], cat_features=cat_cols))[:, 1]
        y_true = y[external_mask]
        m = evaluate_subset(y_true, proba)
        m.update({"ablation": ablation.name, "subset": "external_all"})
        metric_rows.append(m)

    spec = {
        "ablation": ablation.name,
        "description": ablation.description,
        "n_features": len(feature_cols),
        "n_categorical_features": len(cat_cols),
        "n_numeric_features": len(numeric_cols),
        "n_dropped_from_baseline": len(dropped_cols),
        "best_iteration": int(model.get_best_iteration() or 0),
        "features": feature_cols,
        "dropped_from_baseline": dropped_cols,
        "metrics": metric_rows,
    }
    pred = pd.concat(rows, ignore_index=True)
    return pred, spec


def format_metric(value: object, digits: int = 3) -> str:
    if value is None or pd.isna(value):
        return ""
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return f"{float(value):.{digits}f}"


def write_report(metrics: pd.DataFrame, specs: list[dict]) -> None:
    external = metrics[metrics["subset"].isin(["validation", "external_all", "external_hiro", "external_emerge", "external_cardioboost"])].copy()
    lines = [
        "# Leakage and Ablation Audit",
        "",
        "This audit retrained the primary binary CatBoost model after removing reviewer-sensitive feature families. All runs used the same source-held-out split, sample weights, CatBoost recipe, and binary P/LP vs B/LB training slice.",
        "",
        "The purpose is to test whether performance depends on obvious shortcuts such as gene identity, ClinVar/review metadata, VEP consequence/HGVS annotations, or source-record aggregate fields.",
        "",
        "## Ablation Definitions",
        "",
        "| Ablation | Features | Dropped From Baseline | Description |",
        "|---|---:|---:|---|",
    ]
    for spec in specs:
        lines.append(
            f"| {spec['ablation']} | {spec['n_features']} | {spec['n_dropped_from_baseline']} | {spec['description']} |"
        )

    lines.extend(
        [
            "",
            "## Main Performance Summary",
            "",
            "| Ablation | Subset | Rows | AUROC | AUPRC | Sens @0.5 | Spec @0.5 | PPV @0.5 | Deferral | High-conf Accuracy |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for _, r in external.iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    str(r["ablation"]),
                    str(r["subset"]),
                    format_metric(r["rows"], 0),
                    format_metric(r["auroc"]),
                    format_metric(r["auprc"]),
                    format_metric(r["threshold_0_5_sensitivity"]),
                    format_metric(r["threshold_0_5_specificity"]),
                    format_metric(r["threshold_0_5_ppv"]),
                    format_metric(r["deferral_rate"]),
                    format_metric(r["three_zone_high_conf_accuracy"]),
                ]
            )
            + " |"
        )

    baseline = metrics[(metrics["ablation"] == "baseline_current") & (metrics["subset"] == "external_all")]
    lines.extend(["", "## Interpretation", ""])
    if not baseline.empty:
        b = baseline.iloc[0]
        lines.append(
            f"The baseline external-all AUROC was {format_metric(b['auroc'])} and AUPRC was {format_metric(b['auprc'])}. "
            "Compare each ablation against that row. Large drops after removing one feature family would indicate possible dependence on that family."
        )
    lines.append(
        "The strict_low_leakage run is intentionally harsh. It should not be treated as the final clinical model, but it is useful as a reviewer-facing stress test."
    )
    lines.append("")
    lines.append("Saved machine-readable outputs:")
    lines.append("")
    lines.append("- `ablation_metrics.tsv`")
    lines.append("- `ablation_feature_set_summary.tsv`")
    lines.append("- `ablation_specs.json`")
    lines.append("- `predictions_<ablation>.tsv`")
    (OUT / "LEAKAGE_ABLATION_AUDIT_REPORT.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--only",
        action="append",
        default=[],
        help="Run only the named ablation; may be repeated. Existing combined outputs are not rewritten.",
    )
    args = parser.parse_args()
    selected = [ablation for ablation in ABLATIONS if not args.only or ablation.name in args.only]
    missing = sorted(set(args.only) - {ablation.name for ablation in ABLATIONS})
    if missing:
        raise ValueError(f"Unknown ablation(s): {missing}")

    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(TABLE, sep="\t", low_memory=False)
    df = df[df["training_slice_primary_binary"].astype(str).str.lower().eq("true")].copy()
    if SPLIT_COL not in df.columns:
        raise ValueError(f"Missing split column: {SPLIT_COL}")

    all_metrics = []
    specs = []
    feature_summaries = []
    for ablation in selected:
        print(f"[leakage-audit] training {ablation.name}", flush=True)
        pred, spec = train_one(df, ablation)
        pred.to_csv(OUT / f"predictions_{ablation.name}.tsv", sep="\t", index=False)
        (OUT / f"features_{ablation.name}.json").write_text(
            json.dumps(
                {
                    "ablation": spec["ablation"],
                    "description": spec["description"],
                    "features": spec["features"],
                    "dropped_from_baseline": spec["dropped_from_baseline"],
                },
                indent=2,
            )
            + "\n"
        )
        all_metrics.extend(spec["metrics"])
        specs.append(spec)
        feature_summaries.append(
            {
                "ablation": spec["ablation"],
                "description": spec["description"],
                "n_features": spec["n_features"],
                "n_categorical_features": spec["n_categorical_features"],
                "n_numeric_features": spec["n_numeric_features"],
                "n_dropped_from_baseline": spec["n_dropped_from_baseline"],
                "best_iteration": spec["best_iteration"],
            }
        )

    if args.only:
        focused_name = "_".join(args.only)
        pd.DataFrame(all_metrics).to_csv(OUT / f"focused_metrics_{focused_name}.tsv", sep="\t", index=False)
        (OUT / f"focused_specs_{focused_name}.json").write_text(json.dumps(specs, indent=2) + "\n")
        print(f"[leakage-audit] completed focused run(s): {', '.join(args.only)}", flush=True)
        return

    metrics = pd.DataFrame(all_metrics)
    metrics.to_csv(OUT / "ablation_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(feature_summaries).to_csv(OUT / "ablation_feature_set_summary.tsv", sep="\t", index=False)
    (OUT / "ablation_specs.json").write_text(json.dumps(specs, indent=2) + "\n")
    write_report(metrics, specs)
    print(f"[leakage-audit] wrote {OUT / 'LEAKAGE_ABLATION_AUDIT_REPORT.md'}", flush=True)


if __name__ == "__main__":
    main()
