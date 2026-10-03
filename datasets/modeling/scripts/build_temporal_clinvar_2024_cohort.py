#!/usr/bin/env python3
"""Build a label-temporal ClinVar cohort from January 2024 and the current release."""

from __future__ import annotations

import json
import hashlib
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split


ROOT = Path(__file__).resolve().parents[3]
MODEL_TABLE = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
SNAPSHOT = ROOT / "datasets/clinvar/archive/variant_summary_2024-01.txt.gz"
CURRENT = ROOT / "datasets/clinvar/variant_summary.txt.gz"
OUT = ROOT / "datasets/modeling/ready/modeling_table_temporal_clinvar_2024.tsv"
AUDIT = ROOT / "results/model_performance/temporal_clinvar_2024"
USECOLS = [
    "ClinicalSignificance",
    "Assembly",
    "Chromosome",
    "PositionVCF",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF",
    "ReviewStatus",
]
SEED = 20260702


def broad_class(value: object) -> str | None:
    text = str(value).strip().lower().replace("_", " ")
    if not text or text == "nan" or "conflict" in text:
        return None
    has_uncertain = "uncertain significance" in text
    has_pathogenic = "pathogenic" in text
    has_benign = "benign" in text
    if has_uncertain and not has_pathogenic and not has_benign:
        return "VUS"
    if has_pathogenic and not has_benign and not has_uncertain:
        return "Pathogenic"
    if has_benign and not has_pathogenic and not has_uncertain:
        return "Benign"
    return None


def review_stars(value: object) -> int:
    text = str(value).lower().replace("_", " ")
    if "practice guideline" in text:
        return 4
    if "expert panel" in text:
        return 3
    if "multiple submitters" in text and "no conflicts" in text:
        return 2
    if "single submitter" in text or "criteria provided" in text:
        return 1
    return 0


def read_release(path: Path, allowed_variants: set[str]) -> pd.DataFrame:
    pieces = []
    for chunk in pd.read_csv(path, sep="\t", usecols=USECOLS, chunksize=250_000, low_memory=False):
        chunk = chunk[chunk["Assembly"].eq("GRCh38")].copy()
        chunk["class"] = chunk["ClinicalSignificance"].map(broad_class)
        chunk = chunk[chunk["class"].notna()]
        pos = pd.to_numeric(chunk["PositionVCF"], errors="coerce")
        valid = pos.notna() & chunk["ReferenceAlleleVCF"].notna() & chunk["AlternateAlleleVCF"].notna()
        chunk = chunk.loc[valid].copy()
        pos = pos.loc[valid].astype(int).astype(str)
        chunk["variant_id"] = (
            chunk["Chromosome"].astype(str)
            + "-"
            + pos
            + "-"
            + chunk["ReferenceAlleleVCF"].astype(str)
            + "-"
            + chunk["AlternateAlleleVCF"].astype(str)
        )
        chunk = chunk[chunk["variant_id"].isin(allowed_variants)]
        if len(chunk):
            chunk["review_stars"] = chunk["ReviewStatus"].map(review_stars)
            pieces.append(chunk[["variant_id", "class", "review_stars"]])
    if not pieces:
        return pd.DataFrame(columns=["variant_id", "release_class", "review_stars"])
    data = pd.concat(pieces, ignore_index=True)
    rows = []
    for variant_id, group in data.groupby("variant_id", sort=False):
        classes = sorted(group["class"].unique())
        rows.append(
            {
                "variant_id": variant_id,
                "release_class": classes[0] if len(classes) == 1 else None,
                "release_class_count": len(classes),
                "review_stars": int(group["review_stars"].max()),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    AUDIT.mkdir(parents=True, exist_ok=True)
    table = pd.read_csv(MODEL_TABLE, sep="\t", low_memory=False)
    allowed = set(table["variant_id"].astype(str))
    old = read_release(SNAPSHOT, allowed).rename(
        columns={
            "release_class": "clinvar_2024_class",
            "release_class_count": "clinvar_2024_class_count",
            "review_stars": "clinvar_2024_review_stars",
        }
    )
    current = read_release(CURRENT, allowed).rename(
        columns={
            "release_class": "clinvar_current_class",
            "release_class_count": "clinvar_current_class_count",
            "review_stars": "clinvar_current_review_stars",
        }
    )
    labels = old.merge(current, on="variant_id", how="outer", validate="one_to_one")
    labels.to_csv(AUDIT / "clinvar_2024_to_current_exact_allele_labels.tsv", sep="\t", index=False)
    table = table.merge(labels, on="variant_id", how="left", validate="one_to_one")

    old_binary = table["clinvar_2024_class"].isin(["Benign", "Pathogenic"])
    resolved_vus = table["clinvar_2024_class"].eq("VUS") & table["clinvar_current_class"].isin(
        ["Benign", "Pathogenic"]
    )
    in_scope = table["primary_model_inclusion"].eq("include")
    table["training_slice_temporal_binary"] = old_binary | resolved_vus
    table["split_temporal_2024"] = "not_eligible"

    binary_index = table.index[old_binary].to_numpy()
    stratify = table.loc[binary_index, "clinvar_2024_class"].astype(str)
    train_index, validation_index = train_test_split(
        binary_index, test_size=0.15, random_state=SEED, stratify=stratify
    )
    table.loc[train_index, "split_temporal_2024"] = "train"
    table.loc[validation_index, "split_temporal_2024"] = "validation"
    table.loc[resolved_vus, "split_temporal_2024"] = "external_temporal_vus"

    old_binary_in_scope = old_binary & in_scope
    resolved_vus_in_scope = resolved_vus & in_scope
    table["training_slice_temporal_binary_in_scope"] = old_binary_in_scope | resolved_vus_in_scope
    table["split_temporal_2024_in_scope"] = "not_eligible"
    in_scope_binary_index = table.index[old_binary_in_scope].to_numpy()
    in_scope_stratify = table.loc[in_scope_binary_index, "clinvar_2024_class"].astype(str)
    in_scope_train_index, in_scope_validation_index = train_test_split(
        in_scope_binary_index, test_size=0.15, random_state=SEED, stratify=in_scope_stratify
    )
    table.loc[in_scope_train_index, "split_temporal_2024_in_scope"] = "train"
    table.loc[in_scope_validation_index, "split_temporal_2024_in_scope"] = "validation"
    table.loc[resolved_vus_in_scope, "split_temporal_2024_in_scope"] = "external_temporal_vus"

    old_label = table["clinvar_2024_class"].map({"Benign": 0, "Pathogenic": 1})
    current_label = table["clinvar_current_class"].map({"Benign": 0, "Pathogenic": 1})
    table.loc[old_binary, "primary_binary_label"] = old_label.loc[old_binary]
    table.loc[resolved_vus, "primary_binary_label"] = current_label.loc[resolved_vus]
    star_weight = table["clinvar_2024_review_stars"].map({0: 0.5, 1: 0.7, 2: 1.0, 3: 1.2, 4: 1.2})
    table["temporal_source_weight"] = star_weight.fillna(1.0)
    table.to_csv(OUT, sep="\t", index=False)

    transitions = (
        table.loc[table["clinvar_2024_class"].notna()]
        .groupby(["clinvar_2024_class", "clinvar_current_class"], dropna=False)
        .size()
        .rename("rows")
        .reset_index()
    )
    transitions.to_csv(AUDIT / "clinvar_2024_to_current_transitions.tsv", sep="\t", index=False)
    summary = {
        "snapshot": str(SNAPSHOT.relative_to(ROOT)),
        "current": str(CURRENT.relative_to(ROOT)),
        "current_file_sha256": hashlib.sha256(CURRENT.read_bytes()).hexdigest(),
        "current_file_mtime_utc": datetime.fromtimestamp(
            CURRENT.stat().st_mtime, tz=timezone.utc
        ).isoformat(),
        "model_table": str(MODEL_TABLE.relative_to(ROOT)),
        "output": str(OUT.relative_to(ROOT)),
        "exact_allele_matches_2024": int(table["clinvar_2024_class"].notna().sum()),
        "training_rows": int(len(train_index)),
        "validation_rows": int(len(validation_index)),
        "heldout_2024_vus_resolved_currently": int(resolved_vus.sum()),
        "heldout_pathogenic": int(current_label.loc[resolved_vus].sum()),
        "heldout_benign": int((1 - current_label.loc[resolved_vus]).sum()),
        "in_scope_training_rows": int(len(in_scope_train_index)),
        "in_scope_validation_rows": int(len(in_scope_validation_index)),
        "in_scope_heldout_2024_vus_resolved_currently": int(resolved_vus_in_scope.sum()),
        "in_scope_heldout_pathogenic": int(current_label.loc[resolved_vus_in_scope].sum()),
        "in_scope_heldout_benign": int((1 - current_label.loc[resolved_vus_in_scope]).sum()),
        "label_time_boundary": "Only January 2024 B/LB and P/LP labels enter model fitting; every January 2024 VUS is excluded.",
        "limitation": "Predictor annotations were not archived uniformly and may postdate January 2024.",
    }
    (AUDIT / "temporal_cohort_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print("\nTransitions\n", transitions.to_string(index=False))


if __name__ == "__main__":
    main()
