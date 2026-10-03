#!/usr/bin/env python3
"""Create matched eligibility-control and DDG augmentation tables by pLDDT tier."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
DEFAULT_AUDIT = ROOT / "datasets/feature_sources/protein_structure/interim/maximized_foldx_mapping_audit.tsv"
DEFAULT_DDG = ROOT / "datasets/feature_sources/protein_structure/interim/maximized_foldx51_ddg.tsv"
DEFAULT_OUT = ROOT / "datasets/feature_sources/protein_structure/interim/maximized_foldx_augmentations"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", default=str(DEFAULT_AUDIT))
    parser.add_argument("--ddg", default=str(DEFAULT_DDG))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    audit = pd.read_csv(
        args.audit,
        sep="\t",
        usecols=["variant_id", "audit_mapping_status", "audit_residue_plddt"],
        low_memory=False,
    )
    ddg = pd.read_csv(
        args.ddg,
        sep="\t",
        usecols=["variant_id", "ddg_status", "ddg_kcal_mol", "ddg_abs"],
        low_memory=False,
    )
    if not audit["variant_id"].is_unique or not ddg["variant_id"].is_unique:
        raise ValueError("Audit and DDG inputs must each contain one row per variant_id")
    frame = audit.merge(ddg, on="variant_id", how="left", validate="one_to_one")
    frame["audit_residue_plddt"] = pd.to_numeric(frame["audit_residue_plddt"], errors="coerce")
    frame["ddg_kcal_mol"] = pd.to_numeric(frame["ddg_kcal_mol"], errors="coerce")
    frame["ddg_abs"] = pd.to_numeric(frame["ddg_abs"], errors="coerce")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = {"audit": str(Path(args.audit)), "ddg": str(Path(args.ddg)), "tiers": {}}
    for name, minimum in [("ge70", 70.0), ("ge50", 50.0), ("any", 0.0)]:
        mapped = frame["audit_mapping_status"].eq("ok") & frame["audit_residue_plddt"].ge(minimum)
        numeric = mapped & frame["ddg_status"].eq("ok") & frame["ddg_kcal_mol"].notna()
        prefix = f"maxfoldx_{name}"
        eligibility = pd.DataFrame({
            "variant_id": frame["variant_id"],
            f"{prefix}_eligible": mapped,
        })
        features = eligibility.copy()
        features[f"{prefix}_numeric_available"] = numeric
        features[f"{prefix}_ddg_kcal_mol"] = frame["ddg_kcal_mol"].where(numeric)
        features[f"{prefix}_abs_ddg"] = frame["ddg_abs"].where(numeric)
        eligibility_path = out_dir / f"eligibility_only_{name}.tsv"
        features_path = out_dir / f"ddg_features_{name}.tsv"
        eligibility.to_csv(eligibility_path, sep="\t", index=False)
        features.to_csv(features_path, sep="\t", index=False)
        summary["tiers"][name] = {
            "minimum_plddt": minimum,
            "eligible_rows": int(mapped.sum()),
            "numeric_ddg_rows": int(numeric.sum()),
            "eligibility_output": str(eligibility_path),
            "ddg_output": str(features_path),
        }
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
