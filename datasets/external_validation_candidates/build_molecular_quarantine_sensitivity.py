#!/usr/bin/env python3
"""Remove development rows sharing a genomic site or protein change with the external set."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "datasets/external_validation_candidates/expanded_external_2026"


def usable_hgvsp(series: pd.Series) -> pd.Series:
    text = series.fillna("").astype(str).str.strip()
    return text.where(~text.isin(["", "-", "nan", "NA"]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DATA / "modeling_table_expanded_external_annotated.tsv")
    parser.add_argument("--output", type=Path, default=DATA / "modeling_table_expanded_external_molecular_quarantine.tsv")
    parser.add_argument("--audit", type=Path, default=DATA / "molecular_quarantine_rows.tsv")
    parser.add_argument("--summary", type=Path, default=DATA / "molecular_quarantine_summary.json")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", low_memory=False)
    external = df["split_source_heldout"].eq("external_expanded")
    development = df["split_source_heldout"].isin(["train", "validation"])
    site = df["chrom"].astype(str) + ":" + df["pos"].astype(str)
    hgvsp = usable_hgvsp(df["vep_hgvsp"])
    external_sites = set(site[external])
    external_hgvsp = set(hgvsp[external].dropna())
    same_site = development & site.isin(external_sites)
    same_protein = development & hgvsp.isin(external_hgvsp)
    quarantine = same_site | same_protein

    audit = df.loc[quarantine, [
        "variant_id", "primary_gene", "model_label_3class", "split_source_heldout",
        "chrom", "pos", "ref", "alt", "vep_hgvsc", "vep_hgvsp",
    ]].copy()
    audit["same_external_genomic_site"] = same_site[quarantine].to_numpy()
    audit["same_external_hgvsp"] = same_protein[quarantine].to_numpy()
    audit.to_csv(args.audit, sep="\t", index=False)

    df.loc[quarantine, "split_source_heldout"] = "not_eligible"
    df.loc[quarantine, "training_exclude_reason"] = "external_molecular_neighborhood_quarantine"
    df.to_csv(args.output, sep="\t", index=False)
    summary = {
        "external_variants_unchanged": int(external.sum()),
        "development_rows_quarantined": int(quarantine.sum()),
        "same_genomic_site": int(same_site.sum()),
        "same_hgvsp": int(same_protein.sum()),
        "both": int((same_site & same_protein).sum()),
        "remaining_train": int((df["split_source_heldout"] == "train").sum()),
        "remaining_validation": int((df["split_source_heldout"] == "validation").sum()),
        "external_in_development_after_quarantine": int(
            (external & df["split_source_heldout"].isin(["train", "validation"])).sum()
        ),
        "interpretation": "Sensitivity split prevents exact-allele, same-site, and exact-protein-change familiarity.",
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
