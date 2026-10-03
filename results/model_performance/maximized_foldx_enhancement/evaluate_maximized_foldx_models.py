#!/usr/bin/env python3
"""Compare matched baseline, FoldX-eligibility, and DDG models."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MANIFEST = ROOT / "datasets/external_validation_candidates/final_model_heldout_2026/external_unique_variant_manifest.tsv"
DEFAULT_AUDIT = ROOT / "datasets/feature_sources/protein_structure/interim/maximized_foldx_mapping_audit.tsv"
DEFAULT_OUT = ROOT / "results/model_performance/maximized_foldx_enhancement"
MODEL_NAMES = ("baseline", "eligibility", "ddg")


def load_predictions(path: str, name: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    required = {"variant_id", "primary_gene", "y_true", "pathogenic_probability", "split_source_heldout"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{name} predictions missing columns: {sorted(missing)}")
    return frame[list(required)].rename(columns={
        "primary_gene": "gene",
        "pathogenic_probability": f"probability_{name}",
        "y_true": f"y_{name}",
        "split_source_heldout": f"split_{name}",
    })


def score(y: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    return {
        "auroc": float(roc_auc_score(y, probability)),
        "auprc": float(average_precision_score(y, probability)),
        "brier": float(brier_score_loss(y, probability)),
    }


def bootstrap_differences(frame: pd.DataFrame, replicates: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    genes = frame["gene"].dropna().astype(str).unique()
    rows = []
    for replicate in range(replicates):
        sampled = rng.choice(genes, size=len(genes), replace=True)
        parts = [frame[frame["gene"].astype(str).eq(gene)] for gene in sampled]
        sample = pd.concat(parts, ignore_index=True)
        y = sample["y_true"].astype(int).to_numpy()
        if np.unique(y).size < 2:
            continue
        values = {name: score(y, sample[f"probability_{name}"].to_numpy()) for name in MODEL_NAMES}
        for comparison, left, right in [
            ("ddg_minus_baseline", "ddg", "baseline"),
            ("ddg_minus_eligibility", "ddg", "eligibility"),
            ("eligibility_minus_baseline", "eligibility", "baseline"),
        ]:
            for metric in ("auroc", "auprc", "brier"):
                rows.append({
                    "replicate": replicate,
                    "comparison": comparison,
                    "metric": metric,
                    "difference": values[left][metric] - values[right][metric],
                })
    return pd.DataFrame(rows)


def summarize_bootstrap(samples: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (comparison, metric), group in samples.groupby(["comparison", "metric"]):
        values = group["difference"].to_numpy()
        rows.append({
            "comparison": comparison,
            "metric": metric,
            "mean_difference": float(values.mean()),
            "ci_low": float(np.quantile(values, 0.025)),
            "ci_high": float(np.quantile(values, 0.975)),
            "bootstrap_replicates": len(values),
        })
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-predictions", required=True)
    parser.add_argument("--eligibility-predictions", required=True)
    parser.add_argument("--ddg-predictions", required=True)
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    parser.add_argument("--audit", default=str(DEFAULT_AUDIT))
    parser.add_argument(
        "--eligibility-column",
        default="foldx_eligible_plddt_ge70",
        help="Mapping-audit boolean defining the tier-specific covered subset.",
    )
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT))
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260926)
    args = parser.parse_args()

    paths = {
        "baseline": args.baseline_predictions,
        "eligibility": args.eligibility_predictions,
        "ddg": args.ddg_predictions,
    }
    merged = None
    for name in MODEL_NAMES:
        predictions = load_predictions(paths[name], name)
        merged = predictions if merged is None else merged.merge(
            predictions, on=["variant_id", "gene"], how="inner", validate="one_to_one"
        )
    for name in MODEL_NAMES:
        if not merged[f"y_{name}"].equals(merged["y_baseline"]):
            raise ValueError(f"Outcome mismatch in {name} predictions")
    merged["y_true"] = merged["y_baseline"].astype(int)

    audit = pd.read_csv(
        args.audit,
        sep="\t",
        usecols=["variant_id", args.eligibility_column],
        low_memory=False,
    )
    merged = merged.merge(audit, on="variant_id", how="left", validate="one_to_one")
    merged[args.eligibility_column] = merged[args.eligibility_column].fillna(False).astype(bool)
    manifest = pd.read_csv(args.manifest, sep="\t", low_memory=False)
    manifest["manifest_y"] = manifest["combined_label"].map({"Benign": 0, "Pathogenic": 1})

    validation = merged[merged["split_baseline"].eq("validation")].copy()
    heldout = manifest[["variant_id", "manifest_y"]].merge(
        merged, on="variant_id", how="inner", validate="one_to_one"
    )
    if len(heldout) != len(manifest):
        raise ValueError(f"Only {len(heldout)} of {len(manifest)} held-out variants have predictions")
    if not heldout["manifest_y"].astype(int).equals(heldout["y_true"].astype(int)):
        raise ValueError("Held-out manifest labels disagree with model outcome labels")

    subsets = {
        "validation_all": validation,
        "validation_foldx_covered": validation[validation[args.eligibility_column]],
        "heldout_775_all": heldout,
        "heldout_foldx_covered": heldout[heldout[args.eligibility_column]],
    }
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    metric_rows = []
    bootstrap_summaries = []
    for subset_name, subset in subsets.items():
        y = subset["y_true"].astype(int).to_numpy()
        for name in MODEL_NAMES:
            values = score(y, subset[f"probability_{name}"].to_numpy())
            metric_rows.append({
                "subset": subset_name,
                "model": name,
                "rows": len(subset),
                "pathogenic": int(y.sum()),
                "benign": int((y == 0).sum()),
                **values,
            })
        samples = bootstrap_differences(subset, args.bootstrap_replicates, args.seed)
        samples.insert(0, "subset", subset_name)
        samples.to_csv(out_dir / f"gene_bootstrap_samples_{subset_name}.tsv.gz", sep="\t", index=False)
        summary = summarize_bootstrap(samples)
        summary.insert(0, "subset", subset_name)
        bootstrap_summaries.append(summary)

    metrics = pd.DataFrame(metric_rows)
    bootstrap = pd.concat(bootstrap_summaries, ignore_index=True)
    metrics.to_csv(out_dir / "model_metrics.tsv", sep="\t", index=False)
    bootstrap.to_csv(out_dir / "gene_cluster_paired_differences.tsv", sep="\t", index=False)
    heldout.to_csv(out_dir / "heldout_775_predictions.tsv", sep="\t", index=False)
    summary = {
        "models": paths,
        "heldout_manifest": args.manifest,
        "mapping_audit": args.audit,
        "eligibility_column": args.eligibility_column,
        "primary_test": "DDG versus eligibility-only control on heldout_foldx_ge70",
        "bootstrap_unit": "gene",
        "bootstrap_replicates": args.bootstrap_replicates,
        "metrics_output": str(out_dir / "model_metrics.tsv"),
        "paired_output": str(out_dir / "gene_cluster_paired_differences.tsv"),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(metrics.to_string(index=False))
    print("\nPaired gene-cluster bootstrap\n", bootstrap.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
