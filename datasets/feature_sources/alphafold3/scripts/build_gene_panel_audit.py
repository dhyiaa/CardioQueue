#!/usr/bin/env python3
"""Build the AlphaFold 3 gene-panel audit table.

The audit is intentionally gene-level, not variant-level. It defines the AF3
review universe before partner selection, sequence retrieval, or server jobs.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


DEFAULT_READY = Path(
    "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
)
DEFAULT_OUT = Path("results/af3/gene_panel_audit.tsv")


def truthy(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes", "y"})


def nonmissing(series: pd.Series) -> pd.Series:
    if series is None:
        return pd.Series(dtype=bool)
    s = series.astype("string")
    return s.notna() & ~s.str.strip().str.lower().isin({"", "nan", "none", "missing", "."})


def collapse_values(series: pd.Series) -> str:
    vals = sorted(
        {
            str(v)
            for v in series.dropna().unique()
            if str(v).strip() and str(v).strip().lower() not in {"nan", "none"}
        }
    )
    return ";".join(vals)


def infer_panel_bucket(scopes: set[str], inclusions: set[str], gene: str) -> str:
    scopes_l = {s.lower() for s in scopes}
    inclusions_l = {s.lower() for s in inclusions}

    if gene == "TTN" or "ttn_separate" in scopes_l:
        return "ttn_separate"
    if "include" in inclusions_l:
        if "primary" in scopes_l:
            return "primary"
        if "primary_monitor" in scopes_l:
            return "primary_monitor"
        if "expanded_scope" in scopes_l:
            return "expanded_scope"
        return "included_other"
    if "sensitivity_only" in scopes_l:
        return "sensitivity_only"
    if "exclude_or_quarantine" in scopes_l:
        return "exclude_or_quarantine"
    return "other_or_unassigned"


def initial_af3_status(panel_bucket: str) -> tuple[str, bool]:
    if panel_bucket == "exclude_or_quarantine":
        return "excluded_or_quarantine", True
    if panel_bucket == "ttn_separate":
        return "ttn_domain_only_or_defer", True
    if panel_bucket == "sensitivity_only":
        return "manual_review_needed", True
    return "manual_review_needed", False


def af3_inclusion_decision(gene: str, panel_bucket: str) -> str:
    if gene == "KCNE1":
        return "include_for_af3_user_decision"
    if panel_bucket == "exclude_or_quarantine":
        return "exclude_for_af3_initial"
    if panel_bucket == "ttn_separate":
        return "domain_only_or_defer"
    return "include_for_af3_review"


def build_audit(input_path: Path) -> pd.DataFrame:
    df = pd.read_csv(input_path, sep="\t", low_memory=False)

    required = {"primary_gene", "modeling_scope", "primary_model_inclusion", "model_label_3class"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise SystemExit(f"Missing required columns in {input_path}: {missing}")

    df["gene_raw"] = df["primary_gene"].astype("string")
    df = df[nonmissing(df["gene_raw"])].copy()
    df["gene_normalized"] = df["gene_raw"].str.upper()

    if "vep_any_missense" in df.columns:
        missense_mask = truthy(df["vep_any_missense"])
    elif "vep_consequence" in df.columns:
        missense_mask = df["vep_consequence"].astype(str).str.contains("missense_variant", case=False, na=False)
    elif "protein_variant_type" in df.columns:
        missense_mask = df["protein_variant_type"].astype(str).str.contains("substitution|missense", case=False, na=False)
    else:
        missense_mask = pd.Series(False, index=df.index)
    df["_is_missense"] = missense_mask

    hgvsp_cols = [c for c in ["vep_hgvsp", "dbnsfp_HGVSp_VEP", "gnomad_browser_canonical_hgvsp"] if c in df.columns]
    if hgvsp_cols:
        df["_hgvsp_nonmissing"] = False
        for c in hgvsp_cols:
            df["_hgvsp_nonmissing"] |= nonmissing(df[c])
    else:
        df["_hgvsp_nonmissing"] = False

    protein_pos_cols = [c for c in ["protein_position", "vep_protein_position", "foldx_ddg_protein_position"] if c in df.columns]
    if protein_pos_cols:
        df["_protein_position_nonmissing"] = False
        for c in protein_pos_cols:
            df["_protein_position_nonmissing"] |= nonmissing(df[c])
    else:
        df["_protein_position_nonmissing"] = False

    uniprot_cols = [c for c in ["uniprot_accession", "dbnsfp_Uniprot_acc", "foldx_ddg_uniprot_accession"] if c in df.columns]
    if uniprot_cols:
        df["_uniprot_mapped"] = False
        for c in uniprot_cols:
            df["_uniprot_mapped"] |= nonmissing(df[c])
    else:
        df["_uniprot_mapped"] = False

    label = df["model_label_3class"].astype("string")
    df["_is_plp"] = label.eq("Pathogenic")
    df["_is_blb"] = label.eq("Benign")
    df["_is_vus"] = label.eq("VUS")

    rows = []
    for gene, sub in df.groupby("gene_normalized", dropna=True):
        raw_genes = sorted(sub["gene_raw"].dropna().unique())
        scopes = {str(v) for v in sub["modeling_scope"].dropna().unique()}
        inclusions = {str(v) for v in sub["primary_model_inclusion"].dropna().unique()}
        panel_bucket = infer_panel_bucket(scopes, inclusions, gene)
        status, manual = initial_af3_status(panel_bucket)

        case_alias_warnings = []
        if len(raw_genes) > 1:
            case_alias_warnings.append("multiple_raw_gene_strings")
        if gene == "RYR2" and any(r != "RYR2" for r in raw_genes):
            case_alias_warnings.append("case_alias_ryr2")
        if gene == "KCNH2":
            case_alias_warnings.append("kcnH2_present_check_knch2_absent")
        if gene == "KNCH2":
            case_alias_warnings.append("possible_typo_for_KCNH2")
        if gene == "KCNE1" and "exclude_or_quarantine" in {s.lower() for s in scopes} and "include" in {i.lower() for i in inclusions}:
            case_alias_warnings.append("mixed_include_and_quarantine_scope")
            case_alias_warnings.append("include_for_af3_by_user_decision")

        rows.append(
            {
                "gene_normalized": gene,
                "gene_raw_values": ";".join(raw_genes),
                "n_total_variants": int(len(sub)),
                "n_missense": int(sub["_is_missense"].sum()),
                "n_plp": int(sub["_is_plp"].sum()),
                "n_blb": int(sub["_is_blb"].sum()),
                "n_vus": int(sub["_is_vus"].sum()),
                "n_hgvsp_nonmissing": int(sub["_hgvsp_nonmissing"].sum()),
                "n_protein_position_nonmissing": int(sub["_protein_position_nonmissing"].sum()),
                "n_uniprot_mapped": int(sub["_uniprot_mapped"].sum()),
                "modeling_scope_values": collapse_values(sub["modeling_scope"]),
                "primary_model_inclusion_values": collapse_values(sub["primary_model_inclusion"]),
                "panel_bucket": panel_bucket,
                "ttn_separate_flag": gene == "TTN" or panel_bucket == "ttn_separate",
                "excluded_or_quarantine_flag": "exclude_or_quarantine" in {s.lower() for s in scopes},
                "case_alias_warning": ";".join(case_alias_warnings) if case_alias_warnings else "",
                "af3_inclusion_decision": af3_inclusion_decision(gene, panel_bucket),
                "af3_coverage_status": status,
                "af3_complex_group": "",
                "af3_manual_review_needed": bool(manual),
            }
        )

    out = pd.DataFrame(rows).sort_values(
        ["panel_bucket", "gene_normalized"], kind="stable"
    )
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_READY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    audit = build_audit(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(args.output, sep="\t", index=False)

    summary = {
        "input": str(args.input),
        "output": str(args.output),
        "n_genes_normalized": int(len(audit)),
        "n_raw_gene_strings_estimated": int(
            sum(len(str(v).split(";")) for v in audit["gene_raw_values"])
        ),
        "included_normalized_genes": int(
            audit["primary_model_inclusion_values"].str.contains("include", na=False).sum()
        ),
        "panel_bucket_counts": audit["panel_bucket"].value_counts().to_dict(),
        "af3_coverage_status_counts": audit["af3_coverage_status"].value_counts().to_dict(),
        "warnings": {
            "KCNH2_present": bool((audit["gene_normalized"] == "KCNH2").any()),
            "KNCH2_present": bool((audit["gene_normalized"] == "KNCH2").any()),
            "RYR2_raw_values": audit.loc[
                audit["gene_normalized"] == "RYR2", "gene_raw_values"
            ].tolist(),
            "KCNE1_rows": audit.loc[
                audit["gene_normalized"] == "KCNE1",
                [
                    "modeling_scope_values",
                    "primary_model_inclusion_values",
                    "case_alias_warning",
                ],
            ].to_dict("records"),
        },
    }
    summary_path = args.output.with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
