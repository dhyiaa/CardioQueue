#!/usr/bin/env python3
"""Attach corrected VEP consequence QC to lifted eMERGE coordinates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
RESCUE = ROOT / "datasets/emerge/full_arrhythmia_gene_dataset/interim/emerge_coordinate_rescue_liftover.tsv"
MATRIX = ROOT / (
    "datasets/modeling/interim/"
    "final_modeling_table_local_features_hiro_emerge_rescued_gnomad_foldx_vep_hgvs_spliceai.tsv"
)
OUT = ROOT / "datasets/emerge/full_arrhythmia_gene_dataset/interim/emerge_coordinate_rescue_liftover_with_vep_qc.tsv"
SUMMARY = ROOT / "datasets/emerge/full_arrhythmia_gene_dataset/interim/emerge_coordinate_rescue_liftover_with_vep_qc.summary.json"


MISSENSE = {"missense_variant"}
SYNONYMOUS = {"synonymous_variant"}
LOF = {
    "frameshift_variant",
    "stop_gained",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "stop_lost",
    "transcript_ablation",
}
INFRAME = {"inframe_insertion", "inframe_deletion"}
SPLICE_REGION = {"splice_region_variant"}
INTRONIC = {"intron_variant", "non_coding_transcript_intron_variant"}


def consequence_terms(value: object) -> set[str]:
    if pd.isna(value):
        return set()
    out: set[str] = set()
    for chunk in str(value).replace("&", ",").replace("|", ",").split(","):
        term = chunk.strip()
        if term:
            out.add(term)
    return out


def broad_class(value: object) -> str:
    terms = consequence_terms(value)
    if not terms:
        return "missing"
    if terms & LOF:
        return "lof"
    if terms & MISSENSE:
        return "missense"
    if terms & INFRAME:
        return "inframe_indel"
    if terms & SYNONYMOUS:
        return "synonymous"
    if terms & SPLICE_REGION:
        return "splice_region"
    if terms & INTRONIC:
        return "intronic"
    return "other"


def source_class(value: object) -> str:
    text = "" if pd.isna(value) else str(value).lower()
    if not text:
        return "missing"
    if any(k in text for k in ["frameshift", "stop", "nonsense", "splice_acceptor", "splice_donor", "canonical_splice", "loss_of_function", "lof"]):
        return "lof"
    if "missense" in text:
        return "missense"
    if "inframe" in text:
        return "inframe_indel"
    if "synonymous" in text or "silent" in text:
        return "synonymous"
    if "splice" in text:
        return "splice_region"
    if "intron" in text or "intronic" in text:
        return "intronic"
    return "other"


def agreement(source: object, vep: object) -> str:
    src = source_class(source)
    ann = broad_class(vep)
    if src == "missing" or ann == "missing":
        return "missing_source_or_vep"
    if src == ann:
        return "broad_class_match"
    if src == "lof" and ann in {"splice_region", "inframe_indel"}:
        return "partial_lof_related"
    if src == "splice_region" and ann == "lof":
        return "partial_lof_related"
    return "broad_class_mismatch"


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rescue", default=str(RESCUE))
    parser.add_argument("--matrix", default=str(MATRIX))
    parser.add_argument("--out", default=str(OUT))
    parser.add_argument("--summary", default=str(SUMMARY))
    args = parser.parse_args()

    rescue_path = Path(args.rescue)
    matrix_path = Path(args.matrix)
    out_path = Path(args.out)
    summary_path = Path(args.summary)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    rescue = pd.read_csv(rescue_path, sep="\t", low_memory=False)
    matrix_cols = ["variant_id", "vep_status", "vep_consequence", "vep_impact"]
    matrix = pd.read_csv(matrix_path, sep="\t", usecols=lambda c: c in matrix_cols, low_memory=False)
    matrix = matrix.drop_duplicates("variant_id")

    out = rescue.merge(
        matrix.rename(
            columns={
                "variant_id": "grch38_variant_id",
                "vep_status": "vep_status_after_rescue",
                "vep_consequence": "vep_consequence_after_rescue",
                "vep_impact": "vep_impact_after_rescue",
            }
        ),
        on="grch38_variant_id",
        how="left",
        suffixes=("", "_matrix"),
    )
    for col in ["vep_consequence_after_rescue", "vep_impact_after_rescue"]:
        matrix_col = f"{col}_matrix"
        if matrix_col in out.columns:
            out[col] = out[matrix_col].combine_first(out.get(col))
            out = out.drop(columns=[matrix_col])

    out["source_consequence_broad_class"] = out["source_inferred_consequence"].map(source_class)
    out["vep_consequence_broad_class_after_rescue"] = out["vep_consequence_after_rescue"].map(broad_class)
    out["source_vs_vep_consequence_agreement"] = [
        agreement(src, ann) for src, ann in zip(out["source_inferred_consequence"], out["vep_consequence_after_rescue"])
    ]
    out["vep_consequence_qc_flag"] = out["source_vs_vep_consequence_agreement"].map(
        {
            "broad_class_match": "pass",
            "partial_lof_related": "review",
            "broad_class_mismatch": "review",
            "missing_source_or_vep": "missing",
        }
    )
    out.to_csv(out_path, sep="\t", index=False)

    summary = {
        "rescue_input": rel(rescue_path),
        "matrix_input": rel(matrix_path),
        "output": rel(out_path),
        "rows": int(len(out)),
        "liftover_status_counts": {str(k): int(v) for k, v in out["liftover_status"].value_counts(dropna=False).to_dict().items()},
        "coordinate_qc_flag_counts": {
            str(k): int(v) for k, v in out["coordinate_qc_flag"].value_counts(dropna=False).to_dict().items()
        },
        "allele_ref_validation_counts": {
            str(k): int(v) for k, v in out["allele_ref_validation"].value_counts(dropna=False).to_dict().items()
        },
        "vep_status_after_rescue_counts": {
            str(k): int(v) for k, v in out["vep_status_after_rescue"].value_counts(dropna=False).to_dict().items()
        },
        "source_vs_vep_consequence_agreement_counts": {
            str(k): int(v) for k, v in out["source_vs_vep_consequence_agreement"].value_counts(dropna=False).to_dict().items()
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
