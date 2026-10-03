#!/usr/bin/env python3
"""Audit SHaRe and PLOS candidates before building an external manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
PLOS = ROOT / (
    "datasets/external_validation_candidates/plos_cardiovascular_reinterpretation/"
    "pone.0297914.s001.xlsx"
)
PANEL = ROOT / "datasets/gene_panels/cardiogenetics_classifier_genes.keep.txt"


def cdna_key(gene: pd.Series, hgvs: pd.Series) -> pd.Series:
    cdna = hgvs.astype("string").str.extract(r"(c\.[^; )]+)", expand=False)
    cdna = cdna.str.replace(r"[()]", "", regex=True)
    return gene.astype("string") + ":" + cdna


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--share-csv", required=True, help="SHaRe HCM annotation CSV")
    parser.add_argument(
        "--out",
        default=str(ROOT / "results/model_performance/public_external_dataset_audit.summary.json"),
    )
    args = parser.parse_args()

    share_path = Path(args.share_csv)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    panel = {line.strip() for line in PANEL.read_text().splitlines() if line.strip()}
    model = pd.read_csv(
        MODEL,
        sep="\t",
        usecols=["variant_id", "split_source_heldout"],
        low_memory=False,
    )
    share = pd.read_csv(share_path, low_memory=False)
    share_binary = share[share["VarClass"].isin(["B/LB", "P/LP"])].copy()
    share_binary["binary_label"] = share_binary["VarClass"]
    share_binary["cdna_key"] = cdna_key(
        share_binary["external_gene_name"], share_binary["HGVSc"]
    )
    share_local = share_binary.merge(
        model, left_on="VariantID", right_on="variant_id", how="left", validate="one_to_one"
    )

    plos = pd.read_excel(PLOS, header=2)
    plos.columns = ["variant", "original_class", "final_class", "acmg_criteria"]
    plos["gene"] = plos["variant"].astype("string").str.extract(r"^([^:]+):", expand=False)
    plos["cdna_key"] = cdna_key(plos["gene"], plos["variant"])
    plos_binary = plos[
        plos["gene"].isin(panel) & plos["final_class"].isin(["B", "LB", "P", "LP"])
    ].copy()
    plos_binary["binary_label"] = plos_binary["final_class"].map(
        {"B": "B/LB", "LB": "B/LB", "P": "P/LP", "LP": "P/LP"}
    )

    share_key_labels = share_binary.groupby("cdna_key")["binary_label"].nunique()
    plos_key_labels = plos_binary.groupby("cdna_key")["binary_label"].nunique()
    if (share_key_labels > 1).any() or (plos_key_labels > 1).any():
        raise ValueError("A gene/cDNA key has conflicting binary labels within a source")
    share_keys = share_binary.loc[
        share_binary["cdna_key"].notna(), ["cdna_key", "binary_label"]
    ].drop_duplicates()
    plos_keys = plos_binary.loc[
        plos_binary["cdna_key"].notna(), ["cdna_key", "binary_label"]
    ].drop_duplicates()
    cross = share_keys.merge(
        plos_keys,
        on="cdna_key",
        suffixes=("_share", "_plos"),
        validate="one_to_one",
    )
    discordant = cross[cross["binary_label_share"] != cross["binary_label_plos"]]
    existing_external_overlap = share_local["split_source_heldout"].astype("string").str.startswith(
        "external_"
    )
    development_overlap = share_local["split_source_heldout"].isin(["train", "validation"])

    summary = {
        "inputs": {
            "share_csv": str(share_path),
            "share_sha256": sha256(share_path),
            "plos_xlsx": str(PLOS.relative_to(ROOT)),
            "plos_sha256": sha256(PLOS),
        },
        "share": {
            "rows": int(len(share)),
            "binary_rows": int(len(share_binary)),
            "binary_labels": share_binary["binary_label"].value_counts().to_dict(),
            "genes": int(share_binary["external_gene_name"].nunique()),
            "exact_modeling_matrix_overlap": int(share_local["variant_id"].notna().sum()),
            "current_development_overlap": int(development_overlap.sum()),
            "current_external_overlap": int(existing_external_overlap.sum()),
            "potential_new_rows_after_existing_external_precedence": int(
                len(share_binary) - existing_external_overlap.sum()
            ),
        },
        "plos": {
            "listed_rows": int(len(plos)),
            "all_final_classes": plos["final_class"].value_counts(dropna=False).to_dict(),
            "panel_binary_rows": int(len(plos_binary)),
            "panel_binary_labels": plos_binary["binary_label"].value_counts().to_dict(),
            "panel_binary_genes": int(plos_binary["gene"].nunique()),
        },
        "cross_source": {
            "share_plos_gene_cdna_overlap": int(len(cross)),
            "binary_label_discordances": int(len(discordant)),
            "preliminary_nonduplicate_share_plus_plos": int(
                len(share_binary) + len(plos_binary) - cross["cdna_key"].nunique()
            ),
        },
        "warnings": [
            "Gene/cDNA overlap is preliminary; final deduplication requires normalized GRCh38 alleles.",
            "The PLOS source lacks transcript accessions and genomic coordinates.",
            "SHaRe publication reuse requires confirmation because the repository states all rights reserved.",
            "Every evaluation allele must be quarantined before rebuilding train/validation and retraining.",
        ],
    }
    out_path.write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
