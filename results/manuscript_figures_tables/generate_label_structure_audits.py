#!/usr/bin/env python3
"""Generate ClinVar review-quality and protein-structure audit tables."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, confusion_matrix, roc_auc_score


ROOT = Path(__file__).resolve().parents[2]
TABLE = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
PRIMARY_DIR = ROOT / "results/models/cardioqueue_v1_source_heldout"
HIGH_STAR_DIR = ROOT / "results/models/cardioqueue_v1_clinvar_2plus_source_heldout"
OUT = ROOT / "paper_MS_SI/tables"
REPORT = ROOT / "results/model_performance/clinvar_structure_audit/CLINVAR_STRUCTURE_AUDIT.md"
SEED = 20260916
N_BOOT = 5000


def binary_metrics(frame: pd.DataFrame) -> dict[str, float | int]:
    y = frame["y_true"].to_numpy(dtype=int)
    p = frame["pathogenic_probability"].to_numpy(dtype=float)
    tn, fp, fn, tp = confusion_matrix(y, p >= 0.5, labels=[0, 1]).ravel()
    return {
        "rows": len(frame),
        "pathogenic": int(y.sum()),
        "benign": int((y == 0).sum()),
        "auroc": roc_auc_score(y, p),
        "auprc": average_precision_score(y, p),
        "brier": brier_score_loss(y, p),
        "sensitivity_0.5": tp / (tp + fn),
        "specificity_0.5": tn / (tn + fp),
    }


def external_predictions(directory: Path, score_name: str) -> pd.DataFrame:
    path = directory / "primary_binary_predictions_split_source_heldout.tsv"
    frame = pd.read_csv(path, sep="\t")
    frame = frame[frame["split_source_heldout"].str.startswith("external_")].copy()
    return frame[["variant_id", "split_source_heldout", "y_true", "pathogenic_probability"]].rename(
        columns={"pathogenic_probability": score_name}
    )


def paired_bootstrap(primary: pd.DataFrame, high_star: pd.DataFrame) -> pd.DataFrame:
    paired = primary.merge(
        high_star,
        on=["variant_id", "split_source_heldout", "y_true"],
        how="inner",
        validate="one_to_one",
    )
    if len(paired) != len(primary) or len(paired) != len(high_star):
        raise ValueError("Sensitivity predictions do not form an exact paired external set")

    y = paired["y_true"].to_numpy(dtype=int)
    p0 = paired["primary_probability"].to_numpy(dtype=float)
    p1 = paired["high_star_probability"].to_numpy(dtype=float)
    class_indices = [np.flatnonzero(y == value) for value in (0, 1)]
    rng = np.random.default_rng(SEED)
    samples = {"auroc": [], "auprc": [], "brier": []}
    for _ in range(N_BOOT):
        idx = np.concatenate([rng.choice(group, size=len(group), replace=True) for group in class_indices])
        yb = y[idx]
        samples["auroc"].append(roc_auc_score(yb, p1[idx]) - roc_auc_score(yb, p0[idx]))
        samples["auprc"].append(average_precision_score(yb, p1[idx]) - average_precision_score(yb, p0[idx]))
        samples["brier"].append(brier_score_loss(yb, p1[idx]) - brier_score_loss(yb, p0[idx]))

    rows = []
    for metric, values in samples.items():
        arr = np.asarray(values)
        if metric == "auroc":
            observed = roc_auc_score(y, p1) - roc_auc_score(y, p0)
        elif metric == "auprc":
            observed = average_precision_score(y, p1) - average_precision_score(y, p0)
        else:
            observed = brier_score_loss(y, p1) - brier_score_loss(y, p0)
        rows.append(
            {
                "metric": metric,
                "high_star_minus_primary": observed,
                "ci_95_low": np.quantile(arr, 0.025),
                "ci_95_high": np.quantile(arr, 0.975),
                "bootstrap_samples": N_BOOT,
                "seed": SEED,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "training_slice_primary_binary",
        "split_source_heldout",
        "primary_binary_label",
        "clinvar_review_stars",
        "uniprot_accession",
        "alphafold_residue_plddt",
        "dssp_status",
        "freesasa_status",
        "foldx_ddg_kcal_mol",
    ]
    matrix = pd.read_csv(TABLE, sep="\t", usecols=columns, low_memory=False)
    primary_slice = matrix[matrix["training_slice_primary_binary"].eq(True)].copy()

    development = primary_slice[primary_slice["split_source_heldout"].isin(["train", "validation"])].copy()
    stars = (
        development.groupby(["split_source_heldout", "clinvar_review_stars", "primary_binary_label"])
        .size()
        .rename("rows")
        .reset_index()
    )
    stars["label"] = stars["primary_binary_label"].map({0.0: "B/LB", 1.0: "P/LP"})
    stars = stars.drop(columns="primary_binary_label")
    stars["clinvar_review_stars"] = stars["clinvar_review_stars"].astype(int)
    stars.to_csv(OUT / "Table_S13_clinvar_star_distribution.tsv", sep="\t", index=False)

    primary = external_predictions(PRIMARY_DIR, "primary_probability")
    high_star = external_predictions(HIGH_STAR_DIR, "high_star_probability")
    model_rows = []
    for name, frame, train_rows, validation_rows, minimum_stars in [
        ("Primary weighted model", primary.rename(columns={"primary_probability": "pathogenic_probability"}), 36176, 6384, 1),
        ("ClinVar >=2-star sensitivity model", high_star.rename(columns={"high_star_probability": "pathogenic_probability"}), 12635, 2308, 2),
    ]:
        row = {
            "model": name,
            "minimum_clinvar_stars_for_development": minimum_stars,
            "training_rows": train_rows,
            "validation_rows": validation_rows,
        }
        row.update(binary_metrics(frame))
        model_rows.append(row)
    sensitivity = pd.DataFrame(model_rows)
    sensitivity.to_csv(OUT / "Table_S13b_clinvar_star_sensitivity.tsv", sep="\t", index=False)
    paired = paired_bootstrap(primary, high_star)
    paired.to_csv(OUT / "Table_S13c_clinvar_star_paired_differences.tsv", sep="\t", index=False)

    final_coverage = []
    coverage_columns = {
        "Reviewed UniProt accession": "uniprot_accession",
        "AlphaFold residue pLDDT": "alphafold_residue_plddt",
        "FoldX delta-delta-G": "foldx_ddg_kcal_mol",
    }
    split_order = ["train", "validation", "external_hiro", "external_emerge", "external_cardioboost", "not_eligible"]
    for feature, column in coverage_columns.items():
        covered = matrix[column].notna() & matrix[column].astype(str).ne("")
        row = {"feature": feature, "all_matrix_rows": int(covered.sum())}
        for split in split_order:
            row[split] = int((covered & matrix["split_source_heldout"].eq(split)).sum())
        final_coverage.append(row)
    final_coverage.extend(
        [
            {
                "feature": "DSSP numeric residue features",
                "all_matrix_rows": 0,
                **{split: 0 for split in split_order},
            },
            {
                "feature": "FreeSASA numeric residue features",
                "all_matrix_rows": 0,
                **{split: 0 for split in split_order},
            },
        ]
    )
    pd.DataFrame(final_coverage).to_csv(OUT / "Table_S14_structure_final_matrix_coverage.tsv", sep="\t", index=False)

    upstream = pd.DataFrame(
        [
            ["Source protein-feature input", 3591, "Rows from HiRO, eMERGE, and CardioBoost"],
            ["Reviewed UniProt mapping successful", 3284, "Source-level mapping"],
            ["AlphaFold residue pLDDT successful", 2336, "Source-level residue lookup"],
            ["DSSP and FreeSASA calculated", 2336, "Upstream source rows; numeric fields were not retained in the final ready matrix"],
            ["FoldX selected at pLDDT >=70", 1230, "Source rows submitted for stability calculation"],
            ["FoldX successful", 1169, "Source-level calculations"],
            ["FoldX reference mismatch", 61, "Excluded from numeric delta-delta-G; no imputation"],
            ["FoldX numeric value in final matrix", 313, "Variant rows after collapse, joins, and final filtering"],
        ],
        columns=["stage", "rows", "interpretation"],
    )
    upstream.to_csv(OUT / "Table_S14b_structure_pipeline_audit.tsv", sep="\t", index=False)

    importance = pd.read_csv(PRIMARY_DIR / "feature_importance.tsv", sep="\t")
    structure_mask = importance["feature"].str.contains(
        r"^(?:uniprot_|alphafold_|dssp_|freesasa_|foldx_)", case=False, regex=True
    )
    structure_importance = importance[structure_mask].copy()
    structure_importance.to_csv(OUT / "Table_S14c_structure_feature_importance.tsv", sep="\t", index=False)

    p = sensitivity.set_index("model")
    delta = paired.set_index("metric")
    report = f"""# ClinVar review-quality and structure-feature audit

## ClinVar review quality

The primary training set contained 36,176 ClinVar rows: 23,541 with one review star, 12,436 with two stars, and 199 with three stars. One-star rows were retained at weight 0.70; two-star and three-star rows received weights 1.00 and 1.20. Review status was never a predictor.

A separate model trained only on ClinVar records with at least two stars used 12,635 training rows and 2,308 validation rows. On the identical 430 source-held-out variants, its AUROC was {p.loc['ClinVar >=2-star sensitivity model', 'auroc']:.3f}, AUPRC {p.loc['ClinVar >=2-star sensitivity model', 'auprc']:.3f}, and Brier score {p.loc['ClinVar >=2-star sensitivity model', 'brier']:.3f}. The primary model values were {p.loc['Primary weighted model', 'auroc']:.3f}, {p.loc['Primary weighted model', 'auprc']:.3f}, and {p.loc['Primary weighted model', 'brier']:.3f}. The paired AUROC difference was {delta.loc['auroc', 'high_star_minus_primary']:.3f} (95% CI {delta.loc['auroc', 'ci_95_low']:.3f} to {delta.loc['auroc', 'ci_95_high']:.3f}).

## Protein structure

The upstream source pipeline mapped 3,284 of 3,591 rows to reviewed UniProt entries. AlphaFold residue pLDDT was recovered for 2,336 rows. DSSP and FreeSASA were calculated upstream for those rows. FoldX processed 1,230 high-confidence residue mappings, returned 1,169 values, and flagged 61 reference mismatches that were excluded without imputation.

After source rows were collapsed and joined to the final 85,677-row matrix, 524 rows had AlphaFold residue pLDDT and 313 had numeric FoldX values. Neither field was present in training or validation. DSSP and FreeSASA numeric outputs were absent from the final ready matrix. Every retained structure-related feature had zero frozen-model importance. These annotations did not provide learned evidence for the frozen 338-feature baseline. The separate registry-wide enhancement is evaluated in `results/model_performance/structure_enhancement/`.
"""
    REPORT.write_text(report)

    manifest = {
        "matrix": str(TABLE.relative_to(ROOT)),
        "primary_model": str(PRIMARY_DIR.relative_to(ROOT)),
        "high_star_model": str(HIGH_STAR_DIR.relative_to(ROOT)),
        "bootstrap_samples": N_BOOT,
        "seed": SEED,
        "outputs": [
            str(path.relative_to(ROOT))
            for path in sorted(OUT.glob("Table_S1[34]*"))
        ] + [str(REPORT.relative_to(ROOT))],
    }
    (REPORT.parent / "audit_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
