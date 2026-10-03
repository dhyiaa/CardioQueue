#!/usr/bin/env python3
"""Build AF3-draft features and run a primary binary CatBoost smoke test.

This is intentionally a prototype layer. It does not replace the established
modeling table; it writes a new table with af3Draft-prefixed columns.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from catboost import CatBoostClassifier, Pool
from sklearn.metrics import average_precision_score, confusion_matrix, roc_auc_score


ROOT = Path(__file__).resolve().parents[3]
BASE_TABLE = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
TRAIN_SCRIPT = ROOT / "datasets/modeling/scripts/train_primary_binary_catboost.py"
OUT = ROOT / "results/model_performance/af3Draft_smoke_test"
AF3DRAFT_RESULTS = ROOT / "results/af3Draft"
AF3DRAFT_FEATURE_DIR = ROOT / "datasets/feature_sources/af3Draft/processed"
AF3_LONG_FILES = [
    ROOT / "datasets/feature_sources/alphafold3/processed/pilot_v1/af3_pilot_v1_batch_variant_features_long.tsv",
    ROOT / "datasets/feature_sources/alphafold3/processed/folds_2026_07_04_22_36/af3_folds_2026_07_04_22_36_variant_features_long.tsv",
]
BASELINE_PREDICTIONS = (
    ROOT
    / "results/models/primary_binary_catboost_v0_strict_source_heldout_hiro_emerge_rescued/primary_binary_predictions_split_source_heldout.tsv"
)
BASELINE_METRICS = ROOT / "results/models/primary_binary_catboost_v0_strict_source_heldout_hiro_emerge_rescued/metrics.json"
SPLIT_COL = "split_source_heldout"
SEED = 20260702
AF3DRAFT_PROVENANCE_FEATURES = {
    "af3Draft_status",
    "af3Draft_n_long_rows",
    "af3Draft_n_jobs",
    "af3Draft_n_complexes",
    "af3Draft_n_chain_gene_contexts",
    "af3Draft_best_pair_support_class",
    "af3Draft_nearest_partner_genes",
}


def load_primary_trainer():
    spec = importlib.util.spec_from_file_location("primary_binary_trainer", TRAIN_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {TRAIN_SCRIPT}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


TRAINER = load_primary_trainer()


def support_rank(value: object) -> int:
    text = str(value)
    if text.startswith("strong"):
        return 2
    if text.startswith("moderate"):
        return 1
    return 0


def support_label(rank: int) -> str:
    if rank >= 2:
        return "strong_supported_interface"
    if rank == 1:
        return "moderate_supported_interface"
    return "unsupported_or_missing"


def load_af3_long() -> pd.DataFrame:
    frames = []
    for path in AF3_LONG_FILES:
        if not path.exists():
            continue
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        frame["af3Draft_source_file"] = str(path.relative_to(ROOT))
        frames.append(frame)
    if not frames:
        raise FileNotFoundError("No AF3 long feature files were found.")
    df = pd.concat(frames, ignore_index=True)
    df = df[df["variant_id"].notna()].copy()
    df["af3Draft_pair_support_rank_row"] = df["af3_best_pair_support_class"].map(support_rank)
    df["af3Draft_has_residue_feature_row"] = df["af3_residue_plddt"].notna().astype(int)
    return df


def aggregate_af3_features(long_df: pd.DataFrame) -> pd.DataFrame:
    numeric_cols = [
        "af3_residue_plddt",
        "af3_min_distance_to_partner_angstrom",
        "af3_best_pair_iptm",
        "af3_best_pair_min_pae",
        "af3_best_pair_max_contact_probability",
    ]
    for col in numeric_cols:
        long_df[col] = pd.to_numeric(long_df[col], errors="coerce")

    rows = []
    for variant_id, group in long_df.groupby("variant_id", dropna=False):
        ranks = group["af3Draft_pair_support_rank_row"].astype(int)
        trusted = group[ranks >= 1]
        best_source = trusted if not trusted.empty else group
        best_rank = int(ranks.max()) if len(ranks) else 0
        rows.append(
            {
                "variant_id": variant_id,
                "af3Draft_status": "has_af3Draft_output",
                "af3Draft_n_long_rows": int(len(group)),
                "af3Draft_n_jobs": int(group["af3_job_id"].nunique(dropna=True)),
                "af3Draft_n_complexes": int(group["af3_complex_group_id"].nunique(dropna=True)),
                "af3Draft_n_chain_gene_contexts": int(group["af3_chain_gene"].nunique(dropna=True)),
                "af3Draft_has_residue_feature": int(group["af3Draft_has_residue_feature_row"].max()),
                "af3Draft_best_pair_support_rank": best_rank,
                "af3Draft_best_pair_support_class": support_label(best_rank),
                "af3Draft_has_strong_or_moderate_pair": int(best_rank >= 1),
                "af3Draft_has_strong_pair": int(best_rank >= 2),
                "af3Draft_best_residue_plddt": float(best_source["af3_residue_plddt"].max())
                if best_source["af3_residue_plddt"].notna().any()
                else np.nan,
                "af3Draft_mean_residue_plddt": float(best_source["af3_residue_plddt"].mean())
                if best_source["af3_residue_plddt"].notna().any()
                else np.nan,
                "af3Draft_min_distance_to_partner_angstrom": float(
                    best_source["af3_min_distance_to_partner_angstrom"].min()
                )
                if best_source["af3_min_distance_to_partner_angstrom"].notna().any()
                else np.nan,
                "af3Draft_within_partner_interface_5A": int(
                    pd.to_numeric(group["af3_within_partner_interface_5A"], errors="coerce").fillna(0).max()
                ),
                "af3Draft_within_partner_interface_8A": int(
                    pd.to_numeric(group["af3_within_partner_interface_8A"], errors="coerce").fillna(0).max()
                ),
                "af3Draft_best_pair_iptm": float(best_source["af3_best_pair_iptm"].max())
                if best_source["af3_best_pair_iptm"].notna().any()
                else np.nan,
                "af3Draft_best_pair_min_pae": float(best_source["af3_best_pair_min_pae"].min())
                if best_source["af3_best_pair_min_pae"].notna().any()
                else np.nan,
                "af3Draft_best_pair_max_contact_probability": float(
                    best_source["af3_best_pair_max_contact_probability"].max()
                )
                if best_source["af3_best_pair_max_contact_probability"].notna().any()
                else np.nan,
                "af3Draft_jobs": ";".join(sorted(map(str, group["af3_job_id"].dropna().unique()))),
                "af3Draft_complexes": ";".join(sorted(map(str, group["af3_complex_name"].dropna().unique()))),
                "af3Draft_nearest_partner_genes": ";".join(
                    sorted(map(str, group["af3_nearest_partner_gene"].dropna().unique()))
                ),
                "af3Draft_layer_note": "draft_AF3_layer_not_publication_final",
            }
        )
    out = pd.DataFrame(rows)
    flag_cols = [
        "af3Draft_n_long_rows",
        "af3Draft_n_jobs",
        "af3Draft_n_complexes",
        "af3Draft_n_chain_gene_contexts",
        "af3Draft_has_residue_feature",
        "af3Draft_best_pair_support_rank",
        "af3Draft_has_strong_or_moderate_pair",
        "af3Draft_has_strong_pair",
        "af3Draft_within_partner_interface_5A",
        "af3Draft_within_partner_interface_8A",
    ]
    out[flag_cols] = out[flag_cols].fillna(0).astype(int)
    return out


def merge_af3(base: pd.DataFrame, af3: pd.DataFrame) -> pd.DataFrame:
    merged = base.merge(af3, on="variant_id", how="left")
    for col in [
        "af3Draft_n_long_rows",
        "af3Draft_n_jobs",
        "af3Draft_n_complexes",
        "af3Draft_n_chain_gene_contexts",
        "af3Draft_has_residue_feature",
        "af3Draft_best_pair_support_rank",
        "af3Draft_has_strong_or_moderate_pair",
        "af3Draft_has_strong_pair",
        "af3Draft_within_partner_interface_5A",
        "af3Draft_within_partner_interface_8A",
    ]:
        merged[col] = merged[col].fillna(0).astype(int)
    text_defaults = {
        "af3Draft_status": "no_af3Draft_output",
        "af3Draft_best_pair_support_class": "unsupported_or_missing",
        "af3Draft_jobs": "",
        "af3Draft_complexes": "",
        "af3Draft_nearest_partner_genes": "",
        "af3Draft_layer_note": "draft_AF3_layer_not_publication_final",
    }
    for col, value in text_defaults.items():
        merged[col] = merged[col].fillna(value)
    return merged


def prepare_features(df: pd.DataFrame, strict_structural_only: bool) -> tuple[pd.DataFrame, list[str], list[str], list[str]]:
    feature_cols = [c for c in df.columns if not TRAINER.should_drop(c)]
    # Keep direct job names out of the model. They are useful for QC, but too close
    # to batch/gene provenance for a feature smoke test.
    feature_cols = [c for c in feature_cols if c not in {"af3Draft_jobs", "af3Draft_complexes", "af3Draft_layer_note"}]
    if strict_structural_only:
        feature_cols = [c for c in feature_cols if c not in AF3DRAFT_PROVENANCE_FEATURES]
    x = df[feature_cols].copy()
    cat_cols = []
    numeric_cols = []
    for col in x.columns:
        if x[col].dtype == "object" or str(x[col].dtype).startswith("string") or x[col].dtype == bool:
            cat_cols.append(col)
            x[col] = x[col].astype("string").fillna("__MISSING__").astype(str)
        else:
            numeric_cols.append(col)
            x[col] = pd.to_numeric(x[col], errors="coerce")
    return x, feature_cols, cat_cols, numeric_cols


def metric_block(y_true: np.ndarray, proba: np.ndarray, threshold: float = 0.5) -> dict[str, float | int]:
    pred = (proba >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, pred, labels=[0, 1]).ravel()
    return {
        "rows": int(len(y_true)),
        "positives_pathogenic": int(y_true.sum()),
        "negatives_benign": int((y_true == 0).sum()),
        "auroc": float(roc_auc_score(y_true, proba)) if len(np.unique(y_true)) == 2 else np.nan,
        "auprc": float(average_precision_score(y_true, proba)) if len(np.unique(y_true)) == 2 else np.nan,
        "tp": int(tp),
        "fp": int(fp),
        "tn": int(tn),
        "fn": int(fn),
        "sensitivity": float(tp / (tp + fn)) if (tp + fn) else np.nan,
        "specificity": float(tn / (tn + fp)) if (tn + fp) else np.nan,
        "ppv": float(tp / (tp + fp)) if (tp + fp) else np.nan,
        "npv": float(tn / (tn + fn)) if (tn + fn) else np.nan,
    }


def train_af3_model(
    df: pd.DataFrame, model_name: str, strict_structural_only: bool
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    train_df = df[df["training_slice_primary_binary"].astype(str).str.lower().eq("true")].copy()
    x, feature_cols, cat_cols, numeric_cols = prepare_features(train_df, strict_structural_only=strict_structural_only)
    y = train_df["primary_binary_label"].astype(int).to_numpy()
    weights = train_df["sample_weight"].astype(float).to_numpy()
    split = train_df[SPLIT_COL].astype(str)
    train_mask = split.eq("train").to_numpy()
    val_mask = split.eq("validation").to_numpy()

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
        verbose=100,
    )
    model.fit(
        Pool(x.loc[train_mask], y[train_mask], weight=weights[train_mask], cat_features=cat_cols),
        eval_set=Pool(x.loc[val_mask], y[val_mask], weight=weights[val_mask], cat_features=cat_cols),
        use_best_model=True,
    )

    pred_frames = []
    metric_rows = []
    for split_name in ["train", "validation", "external_cardioboost", "external_emerge", "external_hiro"]:
        mask = split.eq(split_name).to_numpy()
        if not mask.any():
            continue
        proba = model.predict_proba(Pool(x.loc[mask], cat_features=cat_cols))[:, 1]
        y_true = y[mask]
        block = metric_block(y_true, proba)
        block.update({"model": model_name, "subset": split_name, "scope": "all_rows"})
        metric_rows.append(block)
        covered = train_df.loc[mask, "af3Draft_has_strong_or_moderate_pair"].to_numpy() == 1
        if covered.any():
            covered_block = metric_block(y_true[covered], proba[covered])
            covered_block.update(
                {"model": model_name, "subset": split_name, "scope": "af3Draft_trusted_pair_rows"}
            )
            metric_rows.append(covered_block)
        tmp = train_df.loc[
            mask,
            [
                "variant_id",
                "primary_gene",
                "model_label_3class",
                SPLIT_COL,
                "af3Draft_status",
                "af3Draft_best_pair_support_class",
                "af3Draft_has_strong_or_moderate_pair",
                "af3Draft_within_partner_interface_5A",
                "af3Draft_within_partner_interface_8A",
            ],
        ].copy()
        tmp["y_true"] = y_true
        tmp["pathogenic_probability"] = proba
        pred_frames.append(tmp)

    predictions = pd.concat(pred_frames, ignore_index=True)
    metrics = pd.DataFrame(metric_rows)
    meta = {
        "model": model_name,
        "input_table": str(BASE_TABLE.relative_to(ROOT)),
        "split_col": SPLIT_COL,
        "n_features": len(feature_cols),
        "n_categorical_features": len(cat_cols),
        "n_numeric_features": len(numeric_cols),
        "n_af3Draft_features_used": sum(c.startswith("af3Draft_") for c in feature_cols),
        "af3Draft_features_used": [c for c in feature_cols if c.startswith("af3Draft_")],
        "strict_structural_only": strict_structural_only,
        "excluded_af3Draft_provenance_features": sorted(AF3DRAFT_PROVENANCE_FEATURES)
        if strict_structural_only
        else [],
        "best_iteration": int(model.get_best_iteration() or 0),
        "notes": [
            "Prototype AF3-draft smoke test.",
            "Same CatBoost hyperparameters and source-held-out split as the primary binary model.",
            "AF3 job IDs and complex names are not used as model features.",
            "Strict structural-only mode also removes AF3 status/count/partner-gene provenance fields.",
        ],
    }
    model.save_model(str(OUT / f"{model_name}_catboost_model.cbm"))
    pd.DataFrame(
        {"feature": feature_cols, "importance": model.get_feature_importance(type="FeatureImportance")}
    ).sort_values("importance", ascending=False).to_csv(OUT / f"{model_name}_feature_importance.tsv", sep="\t", index=False)
    return predictions, metrics, meta


def baseline_covered_metrics(af3_features: pd.DataFrame) -> pd.DataFrame:
    if not BASELINE_PREDICTIONS.exists():
        return pd.DataFrame()
    pred = pd.read_csv(BASELINE_PREDICTIONS, sep="\t")
    cols = ["variant_id", "af3Draft_has_strong_or_moderate_pair", "af3Draft_best_pair_support_class"]
    pred = pred.merge(af3_features[cols].drop_duplicates("variant_id"), on="variant_id", how="left")
    pred["af3Draft_has_strong_or_moderate_pair"] = pred["af3Draft_has_strong_or_moderate_pair"].fillna(0).astype(int)
    rows = []
    split_col = SPLIT_COL if SPLIT_COL in pred.columns else [c for c in pred.columns if c.startswith("split_")][0]
    for split_name, group in pred.groupby(split_col):
        y_true = group["y_true"].astype(int).to_numpy()
        proba = group["pathogenic_probability"].astype(float).to_numpy()
        block = metric_block(y_true, proba)
        block.update({"model": "baseline_primary_binary_existing", "subset": split_name, "scope": "all_rows"})
        rows.append(block)
        covered = group["af3Draft_has_strong_or_moderate_pair"].to_numpy() == 1
        if covered.any():
            covered_block = metric_block(y_true[covered], proba[covered])
            covered_block.update(
                {"model": "baseline_primary_binary_existing", "subset": split_name, "scope": "af3Draft_trusted_pair_rows"}
            )
            rows.append(covered_block)
    return pd.DataFrame(rows)


def write_report(coverage: pd.DataFrame, metrics: pd.DataFrame, metas: list[dict]) -> None:
    def markdown_table(df: pd.DataFrame) -> str:
        if df.empty:
            return "_No rows._"
        cols = list(df.columns)
        lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
        for _, row in df.iterrows():
            values = [str(row[col]) for col in cols]
            lines.append("| " + " | ".join(values) + " |")
        return "\n".join(lines)

    all_metric = metrics.copy()
    show_cols = ["model", "subset", "scope", "rows", "auroc", "auprc", "sensitivity", "specificity", "ppv", "npv"]
    display = all_metric[show_cols].copy()
    for col in ["auroc", "auprc", "sensitivity", "specificity", "ppv", "npv"]:
        display[col] = display[col].map(lambda x: "" if pd.isna(x) else f"{x:.4f}")

    cov_display = coverage.copy()
    for col in ["pct_all_rows", "pct_binary_rows", "pct_vus_rows"]:
        cov_display[col] = cov_display[col].map(lambda x: f"{x:.2f}%" if pd.notna(x) else "")

    feature_sections = []
    for meta in metas:
        feature_sections.append(
            f"""### {meta["model"]}

```text
{chr(10).join(meta["af3Draft_features_used"])}
```"""
        )

    report = f"""# AF3Draft Smoke Test

This is a prototype analysis. The current AlphaFold 3 layer is named **AF3Draft** because the partner list and AF3 job coverage are not final.

## What Was Tested

- Base model: current primary weighted binary CatBoost model.
- Smoke-test models: same CatBoost setup and same `{SPLIT_COL}` split, with only added `af3Draft_` features.
- AF3 job IDs and complex names were kept for QC but excluded from training features.
- A stricter AF3Draft model also excludes AF3 status/count/partner-gene provenance fields and keeps only structure-like columns.
- This test is not a final claim. It asks whether the draft AF3 pathway has enough signal to justify scaling tomorrow.

## AF3Draft Coverage

{markdown_table(cov_display)}

## Performance Comparison

{markdown_table(display)}

## AF3Draft Features Used

{chr(10).join(feature_sections)}

## Interpretation

The most important comparison is not only global performance, because current AF3Draft coverage is small. The key smoke-test readout is whether AF3Draft-covered rows improve or stay stable when these features are added. If global metrics are unchanged, that is expected at this stage because most modelable variants have no AF3Draft output yet.

"""
    (OUT / "AF3DRAFT_SMOKE_TEST_REPORT.md").write_text(report)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    AF3DRAFT_RESULTS.mkdir(parents=True, exist_ok=True)
    AF3DRAFT_FEATURE_DIR.mkdir(parents=True, exist_ok=True)

    base = pd.read_csv(BASE_TABLE, sep="\t", low_memory=False)
    af3_long = load_af3_long()
    af3_features = aggregate_af3_features(af3_long)
    af3_features.to_csv(AF3DRAFT_FEATURE_DIR / "af3Draft_variant_features_aggregated.tsv", sep="\t", index=False)
    af3_long.to_csv(AF3DRAFT_FEATURE_DIR / "af3Draft_variant_features_long_combined.tsv", sep="\t", index=False)

    merged = merge_af3(base, af3_features)
    merged_path = ROOT / "datasets/modeling/interim/modeling_table_with_splits_weights_plus_af3Draft.tsv"
    merged.to_csv(merged_path, sep="\t", index=False)

    coverage_rows = []
    scopes = {
        "all_modeling_rows": merged.index == merged.index,
        "primary_binary_rows": merged["training_slice_primary_binary"].astype(str).str.lower().eq("true"),
        "vus_rows": merged["model_label_3class"].astype(str).eq("VUS"),
    }
    for name, mask in scopes.items():
        sub = merged[mask].copy()
        coverage_rows.append(
            {
                "scope": name,
                "rows": len(sub),
                "genes": sub["primary_gene"].nunique(dropna=True),
                "rows_with_any_af3Draft_output": int(sub["af3Draft_status"].eq("has_af3Draft_output").sum()),
                "rows_with_strong_or_moderate_af3Draft_pair": int(sub["af3Draft_has_strong_or_moderate_pair"].sum()),
                "rows_with_strong_af3Draft_pair": int(sub["af3Draft_has_strong_pair"].sum()),
                "pct_all_rows": 100 * sub["af3Draft_status"].eq("has_af3Draft_output").mean() if len(sub) else np.nan,
                "pct_binary_rows": np.nan,
                "pct_vus_rows": np.nan,
            }
        )
    coverage = pd.DataFrame(coverage_rows)
    binary_den = coverage.loc[coverage["scope"].eq("primary_binary_rows"), "rows"].iloc[0]
    vus_den = coverage.loc[coverage["scope"].eq("vus_rows"), "rows"].iloc[0]
    coverage.loc[:, "pct_binary_rows"] = coverage["rows_with_strong_or_moderate_af3Draft_pair"] / binary_den * 100
    coverage.loc[:, "pct_vus_rows"] = coverage["rows_with_strong_or_moderate_af3Draft_pair"] / vus_den * 100

    by_gene = (
        merged.groupby("primary_gene", dropna=False)
        .agg(
            rows=("variant_id", "count"),
            binary_rows=("training_slice_primary_binary", lambda s: s.astype(str).str.lower().eq("true").sum()),
            vus_rows=("model_label_3class", lambda s: s.astype(str).eq("VUS").sum()),
            rows_with_any_af3Draft_output=("af3Draft_status", lambda s: s.eq("has_af3Draft_output").sum()),
            rows_with_strong_or_moderate_af3Draft_pair=("af3Draft_has_strong_or_moderate_pair", "sum"),
            rows_with_strong_af3Draft_pair=("af3Draft_has_strong_pair", "sum"),
        )
        .reset_index()
        .sort_values(["rows_with_strong_or_moderate_af3Draft_pair", "rows_with_any_af3Draft_output"], ascending=False)
    )
    coverage.to_csv(AF3DRAFT_RESULTS / "af3Draft_variant_coverage_summary.tsv", sep="\t", index=False)
    by_gene.to_csv(AF3DRAFT_RESULTS / "af3Draft_variant_coverage_by_gene.tsv", sep="\t", index=False)

    full_predictions, full_metrics, full_meta = train_af3_model(
        merged, "af3Draft_plus_primary_binary_fullDraft", strict_structural_only=False
    )
    full_predictions.to_csv(OUT / "af3Draft_plus_primary_binary_fullDraft_predictions.tsv", sep="\t", index=False)
    strict_predictions, strict_metrics, strict_meta = train_af3_model(
        merged, "af3Draft_plus_primary_binary_strictStructural", strict_structural_only=True
    )
    strict_predictions.to_csv(
        OUT / "af3Draft_plus_primary_binary_strictStructural_predictions.tsv", sep="\t", index=False
    )
    baseline_metrics = baseline_covered_metrics(af3_features)
    combined_metrics = pd.concat([baseline_metrics, full_metrics, strict_metrics], ignore_index=True)
    combined_metrics.to_csv(OUT / "af3Draft_smoke_test_metrics.tsv", sep="\t", index=False)

    for meta in [full_meta, strict_meta]:
        meta["merged_table"] = str(merged_path.relative_to(ROOT))
        meta["af3Draft_aggregated_features"] = str(
            (AF3DRAFT_FEATURE_DIR / "af3Draft_variant_features_aggregated.tsv").relative_to(ROOT)
        )
        meta["coverage"] = coverage.to_dict(orient="records")
        if BASELINE_METRICS.exists():
            meta["baseline_metrics_file"] = str(BASELINE_METRICS.relative_to(ROOT))
    manifest = {"models": [full_meta, strict_meta]}
    (OUT / "af3Draft_smoke_test_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    write_report(coverage, combined_metrics, [full_meta, strict_meta])

    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
