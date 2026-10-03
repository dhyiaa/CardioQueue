#!/usr/bin/env python3
"""Collapse VEP tab output to one feature row per input variant."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


CONSEQUENCE_SEVERITY = {
    "transcript_ablation": 1,
    "splice_acceptor_variant": 2,
    "splice_donor_variant": 3,
    "stop_gained": 4,
    "frameshift_variant": 5,
    "stop_lost": 6,
    "start_lost": 7,
    "transcript_amplification": 8,
    "inframe_insertion": 9,
    "inframe_deletion": 10,
    "missense_variant": 11,
    "protein_altering_variant": 12,
    "splice_region_variant": 13,
    "incomplete_terminal_codon_variant": 14,
    "start_retained_variant": 15,
    "stop_retained_variant": 16,
    "synonymous_variant": 17,
    "coding_sequence_variant": 18,
    "mature_miRNA_variant": 19,
    "5_prime_UTR_variant": 20,
    "3_prime_UTR_variant": 21,
    "non_coding_transcript_exon_variant": 22,
    "intron_variant": 23,
    "NMD_transcript_variant": 24,
    "non_coding_transcript_variant": 25,
    "upstream_gene_variant": 26,
    "downstream_gene_variant": 27,
    "TFBS_ablation": 28,
    "TFBS_amplification": 29,
    "TF_binding_site_variant": 30,
    "regulatory_region_ablation": 31,
    "regulatory_region_amplification": 32,
    "feature_elongation": 33,
    "regulatory_region_variant": 34,
    "feature_truncation": 35,
    "intergenic_variant": 36,
}

IMPACT_SEVERITY = {"HIGH": 1, "MODERATE": 2, "LOW": 3, "MODIFIER": 4}


def split_terms(value: str) -> list[str]:
    if not isinstance(value, str) or value in {"", "-"}:
        return []
    terms: list[str] = []
    for part in value.replace("&", ",").split(","):
        part = part.strip()
        if part:
            terms.append(part)
    return terms


def worst_consequence(value: str) -> tuple[str, int]:
    terms = split_terms(value)
    if not terms:
        return "", 999
    scored = sorted((CONSEQUENCE_SEVERITY.get(term, 998), term) for term in terms)
    return scored[0][1], scored[0][0]


def collapse_join(series: pd.Series, limit: int = 8) -> str:
    vals = []
    seen = set()
    for val in series.dropna().astype(str):
        if val in {"", "-"}:
            continue
        for part in val.split(","):
            part = part.strip()
            if part and part not in seen:
                vals.append(part)
                seen.add(part)
            if len(vals) >= limit:
                return "|".join(vals)
    return "|".join(vals)


def first_nonmissing(series: pd.Series) -> str:
    for val in series.dropna().astype(str):
        if val not in {"", "-"}:
            return val
    return ""


def numeric_max(series: pd.Series) -> float | None:
    vals = pd.to_numeric(series.replace("-", pd.NA), errors="coerce").dropna()
    if vals.empty:
        return None
    return float(vals.max())


def choose_representative(group: pd.DataFrame) -> pd.Series:
    tmp = group.copy()
    worst = tmp["Consequence"].apply(worst_consequence)
    tmp["_consequence_rank"] = [x[1] for x in worst]
    tmp["_canonical_rank"] = (tmp.get("CANONICAL", "") != "YES").astype(int)
    tmp["_mane_rank"] = (~tmp.get("MANE_SELECT", pd.Series("", index=tmp.index)).fillna("").astype(str).ne("")).astype(int)
    tmp["_impact_rank"] = tmp.get("IMPACT", pd.Series("", index=tmp.index)).map(IMPACT_SEVERITY).fillna(999).astype(int)
    tmp = tmp.sort_values(["_canonical_rank", "_mane_rank", "_impact_rank", "_consequence_rank"], kind="mergesort")
    return tmp.iloc[0]


def read_vep_tab(path: Path) -> pd.DataFrame:
    """Read VEP tab output while preserving the #Uploaded_variation header."""
    header: list[str] | None = None
    header_idx: int | None = None
    with path.open() as fh:
        for idx, line in enumerate(fh):
            line = line.rstrip("\n")
            if line.startswith("#Uploaded_variation"):
                header = line.lstrip("#").split("\t")
                header_idx = idx
                break
    if header is None or header_idx is None:
        raise SystemExit("VEP output is missing #Uploaded_variation header")
    return pd.read_csv(path, sep="\t", dtype=str, names=header, skiprows=header_idx + 1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vep-tab", required=True, type=Path)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--out-features", required=True, type=Path)
    parser.add_argument("--summary-json", required=True, type=Path)
    args = parser.parse_args()

    matrix_ids = pd.read_csv(args.matrix, sep="\t", dtype=str, usecols=["variant_id"])
    matrix_ids = matrix_ids.drop_duplicates("variant_id")

    vep = read_vep_tab(args.vep_tab)
    if "Uploaded_variation" not in vep.columns:
        raise SystemExit("VEP output is missing Uploaded_variation column")

    rows = []
    grouped = vep.groupby("Uploaded_variation", dropna=False, sort=False)
    splice_cols = [c for c in vep.columns if c.startswith("SpliceAI")]
    for variant_id, group in grouped:
        rep = choose_representative(group)
        worst_terms = [worst_consequence(x)[0] for x in group["Consequence"].dropna()]
        worst_terms = [x for x in worst_terms if x]
        worst_term = sorted(worst_terms, key=lambda x: CONSEQUENCE_SEVERITY.get(x, 998))[0] if worst_terms else ""
        row = {
            "variant_id": variant_id,
            "vep_annotation_status": "ok",
            "vep_transcript_count": len(group),
            "vep_symbol": first_nonmissing(group.get("SYMBOL", pd.Series(dtype=str))),
            "vep_gene_id": first_nonmissing(group.get("Gene", pd.Series(dtype=str))),
            "vep_feature": rep.get("Feature", ""),
            "vep_feature_type": rep.get("Feature_type", ""),
            "vep_consequence": rep.get("Consequence", ""),
            "vep_worst_consequence": worst_term,
            "vep_impact": rep.get("IMPACT", ""),
            "vep_canonical": rep.get("CANONICAL", ""),
            "vep_mane_select": rep.get("MANE_SELECT", ""),
            "vep_mane_plus_clinical": rep.get("MANE_PLUS_CLINICAL", ""),
            "vep_hgvsc": rep.get("HGVSc", ""),
            "vep_hgvsp": rep.get("HGVSp", ""),
            "vep_protein_position": rep.get("Protein_position", ""),
            "vep_amino_acids": rep.get("Amino_acids", ""),
            "vep_codons": rep.get("Codons", ""),
            "vep_existing_variation": collapse_join(group.get("Existing_variation", pd.Series(dtype=str))),
            "vep_all_consequences": collapse_join(group["Consequence"], limit=16),
            "vep_all_impacts": collapse_join(group.get("IMPACT", pd.Series(dtype=str))),
            "vep_any_high_impact": bool((group.get("IMPACT", pd.Series(dtype=str)) == "HIGH").any()),
            "vep_any_moderate_impact": bool((group.get("IMPACT", pd.Series(dtype=str)) == "MODERATE").any()),
            "vep_any_missense": bool(group["Consequence"].fillna("").str.contains("missense_variant").any()),
            "vep_any_splice_region": bool(group["Consequence"].fillna("").str.contains("splice").any()),
            "vep_any_frameshift": bool(group["Consequence"].fillna("").str.contains("frameshift_variant").any()),
            "vep_any_stop_gained": bool(group["Consequence"].fillna("").str.contains("stop_gained").any()),
        }
        for col in splice_cols:
            row[f"vep_{col}"] = numeric_max(group[col])
        rows.append(row)

    features = pd.DataFrame(rows)
    merged = matrix_ids.merge(features, on="variant_id", how="left")
    merged["vep_annotation_status"] = merged["vep_annotation_status"].fillna("not_returned")
    merged["vep_annotation_missing_reason"] = ""
    merged.loc[merged["vep_annotation_status"] == "not_returned", "vep_annotation_missing_reason"] = "variant_not_returned_by_vep"
    bool_cols = [c for c in merged.columns if c.startswith("vep_any_")]
    for col in bool_cols:
        merged[col] = merged[col].fillna(False).astype(bool)

    args.out_features.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.out_features, sep="\t", index=False)

    summary = {
        "matrix_variants": int(len(matrix_ids)),
        "vep_raw_rows": int(len(vep)),
        "vep_unique_variants": int(vep["Uploaded_variation"].nunique()),
        "feature_rows": int(len(merged)),
        "not_returned": int((merged["vep_annotation_status"] == "not_returned").sum()),
        "spliceai_columns": splice_cols,
    }
    args.summary_json.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
