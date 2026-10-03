#!/usr/bin/env python3
"""Audit HiRO source records that do not link to the variant registry."""

from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
LINKED = ROOT / "datasets/hiro/full_dataset/interim/hiro_linked_source_records.tsv"
RESCUE = ROOT / "datasets/hiro/full_dataset/interim/hiro_coordinate_rescue_candidates.tsv"
MATRIX = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
OUTDIR = ROOT / "results/data_characteristics/hiro_unresolved_record_audit"
TABLES = OUTDIR / "tables"

MISSING = {"", "missing", "nan", "none", "na", "n/a", "-", "."}


def clean(value: object) -> str:
    text = str(value if value is not None else "").strip()
    return "" if text.lower() in MISSING else text


def split_genes(value: object) -> set[str]:
    return {part.strip().upper() for part in re.split(r"[;,|]", clean(value)) if part.strip()}


def hgvs_terms(*values: object) -> set[str]:
    terms: set[str] = set()
    for value in values:
        text = clean(value)
        if not text:
            continue
        for match in re.finditer(r"([cp]\.[A-Za-z0-9_+\-*>?=\[\]/]+)", text, flags=re.IGNORECASE):
            term = re.sub(r"\s+", "", match.group(1).rstrip(",;|").lower())
            if len(term) >= 4:
                terms.add(term)
        for part in re.split(r"[:;,|]", text):
            term = re.sub(r"\s+", "", part.strip().lower())
            if term.startswith(("c.", "p.")) and len(term) >= 4:
                terms.add(term)
    return terms


def build_matrix_hgvs_index() -> dict[tuple[str, str], list[dict[str, str]]]:
    needed = {
        "variant_id",
        "primary_gene",
        "sources",
        "model_label_3class",
        "vep_hgvsc",
        "vep_hgvsp",
        "vep_symbol",
        "dbnsfp_HGVSc_VEP",
        "dbnsfp_HGVSp_VEP",
        "dbnsfp_genename",
        "alphamissense_direct_protein_variant",
    }
    matrix = pd.read_csv(MATRIX, sep="\t", usecols=lambda c: c in needed, low_memory=False)
    index: dict[tuple[str, str], list[dict[str, str]]] = {}
    for _, row in matrix.iterrows():
        genes = (
            split_genes(row.get("primary_gene", ""))
            | split_genes(row.get("vep_symbol", ""))
            | split_genes(row.get("dbnsfp_genename", ""))
        )
        terms = hgvs_terms(
            row.get("vep_hgvsc", ""),
            row.get("vep_hgvsp", ""),
            row.get("dbnsfp_HGVSc_VEP", ""),
            row.get("dbnsfp_HGVSp_VEP", ""),
            row.get("alphamissense_direct_protein_variant", ""),
        )
        payload = {
            "variant_id": clean(row.get("variant_id", "")),
            "primary_gene": clean(row.get("primary_gene", "")),
            "sources": clean(row.get("sources", "")),
            "model_label_3class": clean(row.get("model_label_3class", "")),
        }
        for gene in genes:
            for term in terms:
                index.setdefault((gene, term), []).append(payload)
    return index


def matrix_hgvs_matches(unresolved: pd.DataFrame) -> pd.DataFrame:
    index = build_matrix_hgvs_index()
    rows = []
    for _, row in unresolved.iterrows():
        hits = []
        for gene in split_genes(row.get("gene", "")):
            for term in hgvs_terms(row.get("hgvs_c", ""), row.get("hgvs_p", ""), row.get("gene_refgene", ""), row.get("variant_key", "")):
                for payload in index.get((gene, term), []):
                    hit = dict(payload)
                    hit["matched_gene"] = gene
                    hit["matched_term"] = term
                    hits.append(hit)

        unique_ids = sorted({hit["variant_id"] for hit in hits if hit["variant_id"]})
        cdna_ids = sorted({hit["variant_id"] for hit in hits if hit["matched_term"].startswith("c.") and hit["variant_id"]})
        protein_ids = sorted({hit["variant_id"] for hit in hits if hit["matched_term"].startswith("p.") and hit["variant_id"]})
        rows.append(
            {
                "variant_uid": row["variant_uid"],
                "matrix_candidate_variant_id_count": len(unique_ids),
                "matrix_candidate_variant_ids": ";".join(unique_ids[:20]),
                "matrix_cdna_candidate_variant_id_count": len(cdna_ids),
                "matrix_cdna_candidate_variant_ids": ";".join(cdna_ids[:20]),
                "matrix_protein_candidate_variant_id_count": len(protein_ids),
                "matrix_protein_candidate_variant_ids": ";".join(protein_ids[:20]),
                "matrix_matched_terms": ";".join(sorted({hit["matched_term"] for hit in hits})[:50]),
                "matrix_matched_sources": ";".join(sorted({hit["sources"] for hit in hits if hit["sources"]})[:20]),
            }
        )
    return pd.DataFrame(rows)


def classify_blocker(row: pd.Series) -> tuple[str, str]:
    status = clean(row.get("hiro_registry_link_status"))
    rescue_status = clean(row.get("coordinate_rescue_status"))
    consequence = clean(row.get("inferred_consequence")).lower()
    hgvs_c = clean(row.get("hgvs_c"))
    hgvs_p = clean(row.get("hgvs_p"))
    gene_refgene = clean(row.get("gene_refgene"))
    chrom = clean(row.get("chrom"))
    pos = clean(row.get("pos"))
    ref = clean(row.get("ref"))
    alt = clean(row.get("alt"))
    gene = clean(row.get("gene"))

    if status == "rescued_but_not_in_promoted_registry":
        return "rescued_coordinate_not_promoted", "Candidate coordinate exists but was not in the promoted registry; review promotion logic."
    if rescue_status == "review_candidate_clinvar_hgvs":
        return "clinvar_review_candidate_not_promoted", "ClinVar candidate exists but evidence was not strong enough for automatic promotion."
    if any([chrom, pos, ref, alt]) and not all([chrom, pos, ref, alt]):
        return "partial_coordinate_or_allele", "Some coordinate fields are present but the VCF allele key is incomplete."
    if not any([hgvs_c, hgvs_p, gene_refgene]):
        return "no_variant_text_or_hgvs", "No usable HGVS/variant text is present in the processed HiRO row."
    if gene.upper() in {"FPGT", "FPGT-TNNI3K", "FPGT;FPGT-TNNI3K", "KNCH2"} or "FPGT" in gene.upper() or gene == "RyR2":
        return "gene_symbol_or_readthrough_cleanup", "Gene symbol/readthrough/alias issue should be resolved before coordinate promotion."
    if consequence in {"frameshift", "indel", "splice_region", "nonsense"}:
        return "needs_hgvs_normalization_indel_splice_lof", "HGVS-only indel/splice/LoF record needs transcript-aware genomic normalization."
    if consequence == "missense":
        return "needs_hgvs_normalization_missense", "HGVS-only missense record needs transcript-aware genomic normalization."
    return "needs_manual_review", "Record has variant text but does not fit a simple rescue category."


def recommend(row: pd.Series) -> str:
    if row.get("matrix_cdna_candidate_variant_id_count", 0) == 1:
        return "Review for promotion: unique internal matrix gene+cDNA candidate."
    if row.get("matrix_protein_candidate_variant_id_count", 0) == 1:
        return "Manual review only: unique protein-level candidate may be isoform-dependent."
    blocker = row["blocker_category"]
    if blocker == "clinvar_review_candidate_not_promoted":
        return "Manual ClinVar/HGVS review; do not auto-promote."
    if blocker == "partial_coordinate_or_allele":
        return "Recover missing allele fields from original VCF/report or normalize HGVS."
    if blocker == "no_variant_text_or_hgvs":
        return "Cannot rescue from current processed table; return to raw source/report."
    if blocker.startswith("needs_hgvs_normalization"):
        return "Run online VEP HGVS, VariantValidator, Mutalyzer, or a local transcript-aware HGVS normalizer."
    if blocker == "gene_symbol_or_readthrough_cleanup":
        return "Clean gene alias/readthrough first, then run HGVS normalization."
    return "Manual review."


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)

    linked = pd.read_csv(LINKED, sep="\t", low_memory=False)
    rescue = pd.read_csv(RESCUE, sep="\t", low_memory=False)
    unresolved = linked[linked["hiro_registry_link_status"].ne("linked_to_registry")].copy()

    match = matrix_hgvs_matches(unresolved)
    audit = unresolved.merge(match, on="variant_uid", how="left")
    rescue_keep = rescue[
        [
            "variant_uid",
            "candidate_count",
            "rescue_status",
            "rescue_confidence",
            "rescued_chrom",
            "rescued_pos",
            "rescued_ref",
            "rescued_alt",
            "rescue_clinvar_variation_id",
            "rescue_clinvar_name",
            "rescue_clinvar_classification",
            "rescue_clinvar_review_status",
            "rescue_matched_terms",
            "rescue_matched_term_types",
        ]
    ]
    audit = audit.merge(rescue_keep, on="variant_uid", how="left", suffixes=("", "_rescue_table"))

    classified = audit.apply(classify_blocker, axis=1, result_type="expand")
    audit["blocker_category"] = classified[0]
    audit["blocker_explanation"] = classified[1]
    audit["recommended_next_action"] = audit.apply(recommend, axis=1)

    cols = [
        "variant_uid",
        "study",
        "source_row",
        "participant_id",
        "target_3class",
        "target_5class",
        "gene",
        "variant_key",
        "hgvs_c",
        "hgvs_p",
        "gene_refgene",
        "inferred_consequence",
        "hiro_registry_link_status",
        "coordinate_rescue_status",
        "coordinate_rescue_confidence",
        "coordinate_rescue_note",
        "chrom",
        "pos",
        "ref",
        "alt",
        "candidate_count",
        "rescue_status",
        "rescue_confidence",
        "rescued_chrom",
        "rescued_pos",
        "rescued_ref",
        "rescued_alt",
        "rescue_clinvar_variation_id",
        "rescue_clinvar_name",
        "rescue_clinvar_classification",
        "rescue_clinvar_review_status",
        "rescue_matched_terms",
        "rescue_matched_term_types",
        "matrix_candidate_variant_id_count",
        "matrix_candidate_variant_ids",
        "matrix_cdna_candidate_variant_id_count",
        "matrix_cdna_candidate_variant_ids",
        "matrix_protein_candidate_variant_id_count",
        "matrix_protein_candidate_variant_ids",
        "matrix_matched_terms",
        "matrix_matched_sources",
        "blocker_category",
        "blocker_explanation",
        "recommended_next_action",
    ]
    audit[cols].to_csv(TABLES / "hiro_unresolved_74_record_audit.tsv", sep="\t", index=False)

    high_priority = audit[
        (audit["matrix_cdna_candidate_variant_id_count"].fillna(0).astype(int).eq(1))
        | audit["coordinate_rescue_status"].eq("review_candidate_clinvar_hgvs")
        | audit["matrix_protein_candidate_variant_id_count"].fillna(0).astype(int).eq(1)
    ].copy()
    high_priority[cols].to_csv(TABLES / "hiro_unresolved_review_candidates.tsv", sep="\t", index=False)

    summary_tables = {
        "status_by_label": pd.crosstab(audit["hiro_registry_link_status"], audit["target_3class"], dropna=False),
        "blocker_by_label": pd.crosstab(audit["blocker_category"], audit["target_3class"], dropna=False),
        "blocker_by_consequence": pd.crosstab(audit["blocker_category"], audit["inferred_consequence"], dropna=False),
        "gene_counts": audit["gene"].value_counts(dropna=False).rename_axis("gene").reset_index(name="rows"),
        "recommended_actions": audit["recommended_next_action"].value_counts(dropna=False).rename_axis("recommended_next_action").reset_index(name="rows"),
    }
    for name, table in summary_tables.items():
        table.to_csv(TABLES / f"{name}.tsv", sep="\t")

    summary = {
        "input_linked_source_records": str(LINKED.relative_to(ROOT)),
        "input_rescue_candidates": str(RESCUE.relative_to(ROOT)),
        "input_modeling_matrix": str(MATRIX.relative_to(ROOT)),
        "output_dir": str(OUTDIR.relative_to(ROOT)),
        "total_hiro_source_records": int(len(linked)),
        "linked_to_registry": int(linked["hiro_registry_link_status"].eq("linked_to_registry").sum()),
        "unresolved_or_not_promoted": int(len(audit)),
        "unresolved_status_counts": Counter(audit["hiro_registry_link_status"]).most_common(),
        "label_counts": Counter(audit["target_3class"]).most_common(),
        "blocker_counts": Counter(audit["blocker_category"]).most_common(),
        "recommendation_counts": Counter(audit["recommended_next_action"]).most_common(),
        "unique_internal_cdna_candidates": int(audit["matrix_cdna_candidate_variant_id_count"].fillna(0).astype(int).eq(1).sum()),
        "unique_internal_protein_candidates": int(audit["matrix_protein_candidate_variant_id_count"].fillna(0).astype(int).eq(1).sum()),
        "clinvar_review_candidates": int(audit["coordinate_rescue_status"].eq("review_candidate_clinvar_hgvs").sum()),
        "notes": [
            "No registry/modeling table is modified by this audit.",
            "The two unique internal cDNA candidates are review/promote candidates, not automatically promoted here.",
            "VEP offline cannot parse HGVS input; remaining records need online VEP, VariantValidator, Mutalyzer, or a local transcript-aware HGVS normalizer.",
        ],
    }
    (OUTDIR / "hiro_unresolved_74_audit.summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_report(audit, summary)
    print(json.dumps(summary, indent=2))


def md_table(df: pd.DataFrame) -> str:
    work = df.copy()
    lines = ["| " + " | ".join(map(str, work.columns)) + " |", "| " + " | ".join(["---"] * len(work.columns)) + " |"]
    for _, row in work.iterrows():
        lines.append("| " + " | ".join(str(row[col]) for col in work.columns) + " |")
    return "\n".join(lines)


def write_report(audit: pd.DataFrame, summary: dict[str, object]) -> None:
    status = pd.crosstab(audit["hiro_registry_link_status"], audit["target_3class"], dropna=False).reset_index()
    blockers = audit["blocker_category"].value_counts().rename_axis("blocker_category").reset_index(name="rows")
    actions = audit["recommended_next_action"].value_counts().rename_axis("recommended_next_action").reset_index(name="rows")
    review = audit[
        (audit["matrix_cdna_candidate_variant_id_count"].fillna(0).astype(int).eq(1))
        | audit["coordinate_rescue_status"].eq("review_candidate_clinvar_hgvs")
        | audit["matrix_protein_candidate_variant_id_count"].fillna(0).astype(int).eq(1)
    ][
        [
            "variant_uid",
            "target_3class",
            "gene",
            "hgvs_c",
            "hgvs_p",
            "inferred_consequence",
            "matrix_cdna_candidate_variant_ids",
            "matrix_protein_candidate_variant_ids",
            "rescue_clinvar_name",
            "recommended_next_action",
        ]
    ]

    report = [
        "# HiRO Unresolved Source-Record Audit",
        "",
        "This audit explains the HiRO records that still do not link to a promoted variant-registry row. It preserves every HiRO source record and does not modify the registry or modeling matrix.",
        "",
        "## Summary",
        "",
        f"- HiRO source records: {summary['total_hiro_source_records']}",
        f"- Linked to registry: {summary['linked_to_registry']}",
        f"- Unresolved or not promoted: {summary['unresolved_or_not_promoted']}",
        f"- Unique internal gene+cDNA candidates: {summary['unique_internal_cdna_candidates']}",
        f"- Unique internal protein-only candidates: {summary['unique_internal_protein_candidates']}",
        f"- ClinVar review candidates from the previous rescue pass: {summary['clinvar_review_candidates']}",
        "",
        "## Status By Label",
        "",
        md_table(status),
        "",
        "## Main Blockers",
        "",
        md_table(blockers),
        "",
        "## Recommended Actions",
        "",
        md_table(actions),
        "",
        "## Review Candidates",
        "",
        md_table(review),
        "",
        "## Interpretation",
        "",
        "Most unresolved rows have HGVS-like variant text but no safe genomic VCF key. The largest blocker is transcript-aware HGVS normalization, especially for missense, indel, frameshift, and splice records. Offline VEP cannot parse HGVS input, so the remaining rescue work needs online VEP, VariantValidator, Mutalyzer, or a local transcript-aware HGVS normalizer.",
        "",
        "The current audit finds two unique internal gene+cDNA candidates and one protein-only internal candidate. These should be reviewed before promotion because transcript choice changes cDNA numbering for some genes. The three older ClinVar review candidates remain manual-review candidates and should not be automatically promoted.",
        "",
        "## Files",
        "",
        "- `tables/hiro_unresolved_74_record_audit.tsv`",
        "- `tables/hiro_unresolved_review_candidates.tsv`",
        "- `tables/blocker_by_label.tsv`",
        "- `tables/blocker_by_consequence.tsv`",
        "- `hiro_unresolved_74_audit.summary.json`",
    ]
    (OUTDIR / "HIRO_UNRESOLVED_74_AUDIT_REPORT.md").write_text("\n".join(report) + "\n")


if __name__ == "__main__":
    main()
