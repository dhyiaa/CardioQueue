#!/usr/bin/env python3
"""Create a model-ready selected dbNSFP feature table.

This script is intentionally non-destructive. It reads the current final
modeling table skeleton and writes a new derived table with clean, stable
feature names. Raw dbNSFP files and earlier joined/interim tables are not
modified.
"""

from __future__ import annotations

import csv
import json
import argparse
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]

INPUT = ROOT / "datasets/modeling/interim/final_modeling_table_skeleton.tsv"
OUT_DIR = ROOT / "datasets/modeling/interim"
OUTPUT = OUT_DIR / "dbnsfp_selected_features.tsv"
SUMMARY = OUT_DIR / "dbnsfp_selected_features.summary.json"

MISSING_VALUES = {"", ".", "NA", "N/A", "na", "nan", "None", "Missing"}

# output_name, input_name, kind
FEATURE_MAP = [
    ("dbnsfp_status", "dbnsfp_status", "metadata"),
    ("dbnsfp_missing_reason", "dbnsfp_missing_reason", "metadata"),
    ("dbnsfp_ref_aa", "dbnsfp_aaref", "annotation"),
    ("dbnsfp_alt_aa", "dbnsfp_aaalt", "annotation"),
    ("dbnsfp_protein_position_raw", "dbnsfp_aapos", "annotation"),
    ("dbnsfp_gene_name_raw", "dbnsfp_genename", "annotation"),
    ("dbnsfp_ensembl_transcript_raw", "dbnsfp_Ensembl_transcriptid", "annotation"),
    ("dbnsfp_uniprot_accession_raw", "dbnsfp_Uniprot_acc", "annotation"),
    ("dbnsfp_hgvsc_vep_raw", "dbnsfp_HGVSc_VEP", "annotation"),
    ("dbnsfp_hgvsp_vep_raw", "dbnsfp_HGVSp_VEP", "annotation"),
    ("dbnsfp_interpro_domain_raw", "dbnsfp_Interpro_domain", "annotation"),
    ("sift_score", "dbnsfp_SIFT_score", "numeric"),
    ("sift_pred", "dbnsfp_SIFT_pred", "categorical"),
    ("polyphen2_hdiv_score", "dbnsfp_Polyphen2_HDIV_score", "numeric"),
    ("polyphen2_hdiv_pred", "dbnsfp_Polyphen2_HDIV_pred", "categorical"),
    ("revel_score", "dbnsfp_REVEL_score", "numeric"),
    ("metalr_score", "dbnsfp_MetaLR_score", "numeric"),
    ("metalr_pred", "dbnsfp_MetaLR_pred", "categorical"),
    ("fathmm_xf_coding_score", "dbnsfp_fathmm-XF_coding_score", "numeric"),
    ("fathmm_xf_coding_pred", "dbnsfp_fathmm-XF_coding_pred", "categorical"),
    ("cadd_raw", "dbnsfp_CADD_raw", "numeric"),
    ("cadd_phred", "dbnsfp_CADD_phred", "numeric"),
    ("alphamissense_dbnsfp_score", "dbnsfp_AlphaMissense_score", "numeric"),
    ("alphamissense_dbnsfp_pred", "dbnsfp_AlphaMissense_pred", "categorical"),
    ("esm1b_score", "dbnsfp_ESM1b_score", "numeric"),
    ("esm1b_pred", "dbnsfp_ESM1b_pred", "categorical"),
    ("gerp_rs", "dbnsfp_GERP++_RS", "numeric"),
    ("phylop100way_vertebrate", "dbnsfp_phyloP100way_vertebrate", "numeric"),
    ("phastcons100way_vertebrate", "dbnsfp_phastCons100way_vertebrate", "numeric"),
    ("gnomad41_joint_af_dbnsfp", "dbnsfp_gnomAD4.1_joint_AF", "numeric"),
    ("gnomad41_joint_nhomalt_dbnsfp", "dbnsfp_gnomAD4.1_joint_nhomalt", "numeric"),
    ("gnomad41_joint_popmax_af_dbnsfp", "dbnsfp_gnomAD4.1_joint_POPMAX_AF", "numeric"),
    ("gnomad41_joint_popmax_nhomalt_dbnsfp", "dbnsfp_gnomAD4.1_joint_POPMAX_nhomalt", "numeric"),
    ("dbnsfp_popmax_af", "dbnsfp_dbNSFP_POPMAX_AF", "numeric"),
    ("dbnsfp_popmax_pop", "dbnsfp_dbNSFP_POPMAX_POP", "categorical"),
]

KEY_COLUMNS = [
    "variant_id",
    "genome_build",
    "chrom",
    "pos",
    "ref",
    "alt",
    "primary_gene",
    "genes",
]

MULTIVALUE_NUMERIC_COLUMNS = {
    "sift_score",
    "revel_score",
    "alphamissense_dbnsfp_score",
    "esm1b_score",
}


def rel(path: Path) -> str:
    path = path.resolve()
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def is_missing(value: object) -> bool:
    return str(value if value is not None else "").strip() in MISSING_VALUES


def clean(value: object) -> str:
    text = str(value if value is not None else "").strip()
    return "" if text in MISSING_VALUES else text


def numeric_summary(value: str, mode: str) -> str:
    """Summarize semicolon-delimited numeric fields conservatively.

    dbNSFP often stores transcript-level values separated by semicolons. For
    damagingness scores where larger usually means more deleterious, max is a
    useful transcript-collapsed feature. For SIFT, smaller is more damaging, so
    min is used.
    """

    values: list[float] = []
    for piece in str(value or "").replace("|", ";").split(";"):
        piece = piece.strip()
        if piece in MISSING_VALUES:
            continue
        try:
            values.append(float(piece))
        except ValueError:
            continue
    if not values:
        return ""
    result = min(values) if mode == "min" else max(values)
    return f"{result:.8g}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--summary", type=Path, default=SUMMARY)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    counts: Counter = Counter()
    feature_kinds = {out: kind for out, _, kind in FEATURE_MAP}

    fieldnames = KEY_COLUMNS + [out for out, _, _ in FEATURE_MAP]
    fieldnames += [f"{out}_is_missing" for out, _, kind in FEATURE_MAP if kind in {"numeric", "categorical"}]
    fieldnames += [
        "sift_score_min",
        "revel_score_max",
        "alphamissense_dbnsfp_score_max",
        "esm1b_score_max",
    ]

    with args.input.open(newline="") as in_handle, args.output.open("w", newline="") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row in reader:
            counts["rows"] += 1
            out_row = {col: clean(row.get(col, "")) for col in KEY_COLUMNS}

            for out_name, in_name, kind in FEATURE_MAP:
                value = clean(row.get(in_name, ""))
                out_row[out_name] = value
                if kind in {"numeric", "categorical"}:
                    missing = is_missing(value)
                    out_row[f"{out_name}_is_missing"] = "true" if missing else "false"
                    counts[f"{out_name}_{'missing' if missing else 'present'}"] += 1

            out_row["sift_score_min"] = numeric_summary(out_row["sift_score"], "min")
            out_row["revel_score_max"] = numeric_summary(out_row["revel_score"], "max")
            out_row["alphamissense_dbnsfp_score_max"] = numeric_summary(
                out_row["alphamissense_dbnsfp_score"], "max"
            )
            out_row["esm1b_score_max"] = numeric_summary(out_row["esm1b_score"], "max")

            counts[f"dbnsfp_status_{out_row['dbnsfp_status'] or 'blank'}"] += 1
            writer.writerow(out_row)

    summary = {
        "inputs": {"modeling_skeleton": rel(args.input)},
        "outputs": {
            "selected_features": rel(args.output),
            "summary": rel(args.summary),
        },
        "counts": dict(sorted(counts.items())),
        "selected_columns": [
            {"output": out, "input": inp, "kind": kind}
            for out, inp, kind in FEATURE_MAP
        ],
        "derived_columns": {
            "sift_score_min": "minimum transcript-level SIFT score parsed from dbNSFP semicolon field",
            "revel_score_max": "maximum transcript-level REVEL score parsed from dbNSFP semicolon field",
            "alphamissense_dbnsfp_score_max": "maximum transcript-level AlphaMissense score parsed from dbNSFP semicolon field",
            "esm1b_score_max": "maximum transcript-level ESM1b score parsed from dbNSFP semicolon field",
        },
        "notes": [
            "Non-destructive derived table; no source dbNSFP or modeling skeleton files are edited.",
            "Original raw dbNSFP columns remain in the modeling skeleton and local feature registry.",
            "Missing flags are included for numeric and categorical model features.",
        ],
    }
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
