#!/usr/bin/env python3
"""Merge parsed VEP features into the current modeling matrix."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--vep-features", required=True, type=Path)
    parser.add_argument("--out-matrix", required=True, type=Path)
    parser.add_argument("--summary-json", required=True, type=Path)
    args = parser.parse_args()

    matrix = pd.read_csv(args.matrix, sep="\t", dtype=str)
    vep = pd.read_csv(args.vep_features, sep="\t", dtype=str)
    if matrix["variant_id"].duplicated().any():
        raise SystemExit("Input matrix has duplicate variant_id rows")
    if vep["variant_id"].duplicated().any():
        raise SystemExit("VEP feature table has duplicate variant_id rows")

    drop_cols = [c for c in vep.columns if c in matrix.columns and c != "variant_id"]
    matrix = matrix.drop(columns=drop_cols)
    out = matrix.merge(vep, on="variant_id", how="left")

    status = out.get("vep_annotation_status")
    if status is not None:
        out["vep_status"] = out["vep_annotation_status"].where(out["vep_annotation_status"].eq("ok"), "missing")
        out["vep_missing_reason"] = out["vep_annotation_missing_reason"].fillna("")
        out.loc[out["vep_status"].eq("ok"), "vep_missing_reason"] = ""

    splice_cols = [c for c in out.columns if c.startswith("vep_SpliceAI")]
    if splice_cols:
        nonmissing = out[splice_cols].notna().any(axis=1)
        out["spliceai_status"] = nonmissing.map({True: "ok", False: "missing"})
        out["spliceai_missing_reason"] = ""
        out.loc[~nonmissing, "spliceai_missing_reason"] = "no_spliceai_score_returned"

    args.out_matrix.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_matrix, sep="\t", index=False)

    summary = {
        "input_rows": int(len(matrix)),
        "output_rows": int(len(out)),
        "output_columns": int(out.shape[1]),
        "vep_feature_columns": int(len([c for c in vep.columns if c != "variant_id"])),
        "vep_status_counts": out.get("vep_status", pd.Series(dtype=str)).fillna("missing").value_counts().to_dict(),
        "spliceai_status_counts": out.get("spliceai_status", pd.Series(dtype=str)).fillna("missing").value_counts().to_dict(),
    }
    args.summary_json.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
