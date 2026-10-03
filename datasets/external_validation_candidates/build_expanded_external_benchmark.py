#!/usr/bin/env python3
"""Build a leakage-controlled, multi-source external benchmark and modeling table."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd
import pysam
from sklearn.model_selection import train_test_split


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MATRIX = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
DEFAULT_PLOS = ROOT / (
    "datasets/external_validation_candidates/plos_cardiovascular_reinterpretation/"
    "pone.0297914.s001.xlsx"
)
DEFAULT_CLINVAR = ROOT / "datasets/clinvar/variant_summary.txt.gz"
DEFAULT_PANEL = ROOT / "datasets/gene_panels/cardiogenetics_classifier_genes.keep.txt"
DEFAULT_FASTA = ROOT / (
    "datasets/feature_sources/insilico/raw/vep/fasta/"
    "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
)
DEFAULT_OUT = ROOT / "datasets/external_validation_candidates/expanded_external_2026"

LABEL_NUMERIC = {"Benign": 0, "Pathogenic": 1}
SHARE_LABELS = {"B/LB": "Benign", "P/LP": "Pathogenic"}
PLOS_LABELS = {"B": "Benign", "LB": "Benign", "P": "Pathogenic", "LP": "Pathogenic"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_variant_id(value: str) -> tuple[str, int, str, str]:
    parts = str(value).strip().split("-", 3)
    if len(parts) != 4:
        raise ValueError(f"Invalid chrom-pos-ref-alt identifier: {value!r}")
    chrom, pos, ref, alt = parts
    return chrom.removeprefix("chr"), int(pos), ref.upper(), alt.upper()


def reference_matches(fasta: pysam.FastaFile, variant_id: str) -> tuple[bool, str]:
    try:
        chrom, pos, ref, _alt = parse_variant_id(variant_id)
    except (TypeError, ValueError):
        return False, "invalid_chrom_pos_ref_alt_identifier"
    contig = chrom if chrom in fasta.references else f"chr{chrom}"
    if contig not in fasta.references:
        return False, "contig_not_in_reference"
    observed = fasta.fetch(contig, pos - 1, pos - 1 + len(ref)).upper()
    if observed != ref:
        return False, f"reference_mismatch_expected_{ref}_observed_{observed}"
    return True, ""


def load_share(path: Path, panel: set[str]) -> pd.DataFrame:
    raw = pd.read_csv(path, low_memory=False)
    selected = raw[raw["VarClass"].isin(SHARE_LABELS)].copy()
    selected["variant_id"] = selected["VariantID"].astype(str).str.strip()
    selected["source"] = "SHaRe_HCM_2026Q1"
    selected["source_record_id"] = selected["variant_id"]
    selected["external_label"] = selected["VarClass"].map(SHARE_LABELS)
    selected["source_class"] = selected["VarClass"]
    selected["gene"] = selected["external_gene_name"].astype(str).str.strip()
    selected["reported_hgvsc"] = selected["HGVSc"].fillna("").astype(str)
    selected["reported_hgvsp"] = selected["HGVSp"].fillna("").astype(str)
    selected["mapping_method"] = "source_GRCh38_VariantID"
    selected["coordinate_mapping_status"] = "mapped"
    selected["in_gene_panel"] = selected["gene"].isin(panel)
    if selected["variant_id"].duplicated().any():
        raise ValueError("SHaRe binary table contains duplicate VariantID values")
    return selected


def load_plos(path: Path, panel: set[str]) -> pd.DataFrame:
    raw = pd.read_excel(path, header=2)
    raw.columns = ["variant", "original_class", "final_class", "acmg_criteria"]
    raw = raw[raw["variant"].notna()].copy()
    raw["gene"] = raw["variant"].astype(str).str.extract(r"^([^:]+):", expand=False)
    raw["reported_cdna"] = raw["variant"].astype(str).str.extract(r":(c\.[^; ]+)", expand=False)
    raw["reported_protein"] = raw["variant"].astype(str).str.extract(r"p\.\(([^)]+)\)", expand=False)
    raw = raw[raw["gene"].isin(panel) & raw["final_class"].isin(PLOS_LABELS)].copy()
    raw["external_label"] = raw["final_class"].map(PLOS_LABELS)
    raw["source"] = "Fernandez_Falgueras_PLOS_2024"
    raw["source_record_id"] = raw["variant"].astype(str)
    raw["source_class"] = raw["final_class"]
    raw["reported_hgvsc"] = raw["reported_cdna"].fillna("")
    raw["reported_hgvsp"] = raw["reported_protein"].fillna("")
    raw["mapping_method"] = "unmapped"
    raw["coordinate_mapping_status"] = "unmapped"
    raw["in_gene_panel"] = True
    return raw


def resolve_plos_with_clinvar(plos: pd.DataFrame, clinvar_path: Path) -> pd.DataFrame:
    wanted = set((plos["gene"] + ":" + plos["reported_cdna"]).dropna())
    hits: dict[str, set[tuple[str, str, str, str, str]]] = defaultdict(set)
    with gzip.open(clinvar_path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["Assembly"] != "GRCh38":
                continue
            match = re.search(r":(c\.[^ ]+)", row["Name"])
            if not match:
                continue
            key = f"{row['GeneSymbol']}:{match.group(1)}"
            if key not in wanted:
                continue
            coordinate = (
                row["Chromosome"],
                row["PositionVCF"],
                row["ReferenceAlleleVCF"],
                row["AlternateAlleleVCF"],
                row["VariationID"],
            )
            if all(value not in {"", "na", "-"} for value in coordinate[:4]):
                hits[key].add(coordinate)

    out = plos.copy()
    out["variant_id"] = pd.NA
    out["clinvar_variation_id_for_coordinate"] = pd.NA
    for idx, row in out.iterrows():
        key = f"{row['gene']}:{row['reported_cdna']}"
        coordinates = hits.get(key, set())
        unique_alleles = {(c, p, r, a) for c, p, r, a, _ in coordinates}
        if len(unique_alleles) == 1:
            chrom, pos, ref, alt = next(iter(unique_alleles))
            out.at[idx, "variant_id"] = f"{chrom}-{pos}-{ref.upper()}-{alt.upper()}"
            variation_ids = sorted({variation_id for *_, variation_id in coordinates})
            out.at[idx, "clinvar_variation_id_for_coordinate"] = "|".join(variation_ids)
            out.at[idx, "mapping_method"] = "ClinVar_GRCh38_exact_gene_cdna"
            out.at[idx, "coordinate_mapping_status"] = "mapped_unique"
        elif len(unique_alleles) > 1:
            out.at[idx, "coordinate_mapping_status"] = "ambiguous_multiple_GRCh38_alleles"
    return out


def legacy_observations(matrix: pd.DataFrame) -> pd.DataFrame:
    selected = matrix[matrix["split_source_heldout"].astype(str).str.startswith("external_")].copy()
    source = selected["split_source_heldout"].str.replace("external_", "legacy_", regex=False)
    return pd.DataFrame(
        {
            "source": source,
            "source_record_id": selected["variant_id"],
            "variant_id": selected["variant_id"],
            "gene": selected["primary_gene"],
            "external_label": selected["model_label_3class"],
            "source_class": selected["model_label_3class"],
            "reported_hgvsc": selected["vep_hgvsc"].fillna(""),
            "reported_hgvsp": selected["vep_hgvsp"].fillna(""),
            "mapping_method": "existing_exact_GRCh38_matrix_allele",
            "coordinate_mapping_status": "mapped",
            "in_gene_panel": True,
        }
    )


def safe_strata(df: pd.DataFrame) -> pd.Series:
    gene = df["primary_gene"].fillna("UNK")
    counts = gene.value_counts()
    gene_bin = gene.where(gene.map(counts) >= 100, "OTHER_LOW_N")
    strata = df["model_label_3class"].astype(str) + "|" + gene_bin.astype(str)
    if strata.value_counts().min() < 2:
        strata = df["model_label_3class"].astype(str)
    return strata


def populate_new_row(columns: list[str], row: pd.Series) -> dict[str, object]:
    chrom, pos, ref, alt = parse_variant_id(row["variant_id"])
    out: dict[str, object] = {column: pd.NA for column in columns}
    out.update(
        {
            "variant_id": row["variant_id"],
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "modeling_row_id": row["variant_id"],
            "genome_build": "GRCh38",
            "coordinate_key": row["variant_id"],
            "primary_gene": row["gene"],
            "genes": row["gene"],
            "modeling_scope": "cardiogenetics_classifier",
            "primary_model_inclusion": "include",
            "model_label_3class": row["combined_label"],
            "model_label_numeric": LABEL_NUMERIC[row["combined_label"]],
            "clean_supervised_label": row["combined_label"],
            "label_conflict_type": "none",
            "source_leakage_group": row["variant_id"],
            "sources": row["external_sources"],
            "source_count": row["source_count"],
            "source_row_count": row["source_observation_count"],
            "in_clinvar": False,
            "in_hiro": False,
            "in_emerge": False,
            "in_cardioboost": False,
            "phenotype_available": False,
            "private_or_patient_linked": False,
            "labels_3class": row["combined_label"],
            "raw_labels": row["combined_label"],
            "gene_panel_decisions": "keep",
            "training_slice_primary_binary": True,
            "training_slice_high_confidence_binary": True,
            "training_slice_three_class_exploratory": True,
            "training_slice_vus_scoring": False,
            "training_exclude_reason": "eligible_external_only",
            "primary_binary_label": LABEL_NUMERIC[row["combined_label"]],
            "source_confidence_weight": 1.0,
            "class_weight": 2.0 if row["combined_label"] == "Pathogenic" else 1.0,
            "sample_weight": 2.0 if row["combined_label"] == "Pathogenic" else 1.0,
            "split_internal_grouped": "not_eligible",
            "split_source_heldout": "external_expanded",
            "split_gene_stress": "not_eligible",
        }
    )
    status_defaults = {
        "gnomad_browser_status": "not_annotated",
        "vep_status": "pending_external_annotation",
        "spliceai_status": "pending_external_annotation",
        "clingen_gene_validity_status": "not_annotated",
        "clingen_variant_evidence_status": "not_annotated",
        "dssp_status": "not_annotated",
        "freesasa_status": "not_annotated",
        "foldx_ddg_status": "not_available",
        "alphamissense_direct_status": "pending_external_annotation",
        "dbnsfp_status": "pending_external_annotation",
        "protein_feature_status": "pending_external_annotation",
        "alphafold_plddt_status": "pending_external_annotation",
    }
    for column, value in status_defaults.items():
        if column in out:
            out[column] = value
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--share-csv", type=Path, required=True)
    parser.add_argument("--plos-xlsx", type=Path, default=DEFAULT_PLOS)
    parser.add_argument("--clinvar", type=Path, default=DEFAULT_CLINVAR)
    parser.add_argument("--panel", type=Path, default=DEFAULT_PANEL)
    parser.add_argument("--reference-fasta", type=Path, default=DEFAULT_FASTA)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--seed", type=int, default=20260925)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    panel = {line.strip() for line in args.panel.read_text().splitlines() if line.strip()}
    matrix = pd.read_csv(args.matrix, sep="\t", low_memory=False)
    if not matrix["variant_id"].is_unique:
        raise ValueError("Modeling matrix must have one row per variant_id")

    share = load_share(args.share_csv, panel)
    plos = resolve_plos_with_clinvar(load_plos(args.plos_xlsx, panel), args.clinvar)
    observations = pd.concat(
        [legacy_observations(matrix), share, plos], ignore_index=True, sort=False
    )
    observations["reference_match"] = False
    observations["reference_validation_reason"] = "coordinate_unmapped"
    fasta = pysam.FastaFile(str(args.reference_fasta))
    for idx, value in observations["variant_id"].dropna().items():
        ok, reason = reference_matches(fasta, str(value))
        observations.at[idx, "reference_match"] = ok
        observations.at[idx, "reference_validation_reason"] = reason
    fasta.close()
    observations["external_binary_label"] = observations["external_label"].map(LABEL_NUMERIC)

    mapped = observations[
        observations["variant_id"].notna()
        & observations["reference_match"]
        & observations["in_gene_panel"]
    ].copy()
    if mapped["external_binary_label"].isna().any():
        raise ValueError("Mapped external observations contain a non-binary label")

    unique_rows = []
    for variant_id, group in mapped.groupby("variant_id", sort=True):
        labels = sorted(set(group["external_label"]))
        genes = sorted(set(group["gene"].dropna().astype(str)))
        sources = sorted(set(group["source"].astype(str)))
        unique_rows.append(
            {
                "variant_id": variant_id,
                "gene": genes[0] if len(genes) == 1 else "|".join(genes),
                "external_sources": "|".join(sources),
                "source_count": len(sources),
                "source_observation_count": len(group),
                "external_labels": "|".join(labels),
                "label_concordant": len(labels) == 1,
                "combined_label": labels[0] if len(labels) == 1 else pd.NA,
                "already_in_modeling_matrix": variant_id in set(matrix["variant_id"]),
            }
        )
    unique = pd.DataFrame(unique_rows)
    conflicts = unique[~unique["label_concordant"]].copy()
    benchmark = unique[unique["label_concordant"]].copy()
    benchmark_ids = set(benchmark["variant_id"])
    quarantine_ids = set(unique["variant_id"])

    expanded = matrix.copy()
    selected_existing = expanded["variant_id"].isin(benchmark_ids)
    label_lookup = benchmark.set_index("variant_id")["combined_label"]
    expanded.loc[selected_existing, "model_label_3class"] = expanded.loc[
        selected_existing, "variant_id"
    ].map(label_lookup)
    expanded.loc[selected_existing, "primary_binary_label"] = expanded.loc[
        selected_existing, "model_label_3class"
    ].map(LABEL_NUMERIC)
    expanded.loc[selected_existing, "training_slice_primary_binary"] = True
    expanded.loc[selected_existing, "training_exclude_reason"] = "eligible_external_only"
    expanded.loc[selected_existing, "split_source_heldout"] = "external_expanded"

    conflict_existing = expanded["variant_id"].isin(set(conflicts["variant_id"]))
    expanded.loc[conflict_existing, "training_slice_primary_binary"] = False
    expanded.loc[conflict_existing, "training_exclude_reason"] = "external_source_label_conflict"
    expanded.loc[conflict_existing, "split_source_heldout"] = "external_label_conflict"

    missing = benchmark[~benchmark["already_in_modeling_matrix"]].copy()
    if len(missing):
        appended = pd.DataFrame(
            [populate_new_row(list(expanded.columns), row) for _, row in missing.iterrows()],
            columns=expanded.columns,
        )
        expanded = pd.concat([expanded, appended], ignore_index=True)

    binary = (
        expanded["training_slice_primary_binary"].astype(str).str.lower().eq("true")
        & expanded["model_label_3class"].isin(LABEL_NUMERIC)
    )
    development = binary & ~expanded["variant_id"].isin(quarantine_ids)
    dev_idx = expanded.index[development].to_numpy()
    train_idx, validation_idx = train_test_split(
        dev_idx,
        train_size=0.85,
        random_state=args.seed,
        stratify=safe_strata(expanded.loc[dev_idx]),
    )
    expanded.loc[development, "split_source_heldout"] = "not_assigned"
    expanded.loc[train_idx, "split_source_heldout"] = "train"
    expanded.loc[validation_idx, "split_source_heldout"] = "validation"
    if expanded.loc[expanded["variant_id"].isin(quarantine_ids), "split_source_heldout"].isin(
        ["train", "validation"]
    ).any():
        raise AssertionError("An external benchmark allele remained in development")

    observation_columns = [
        "source", "source_record_id", "variant_id", "gene", "external_label",
        "external_binary_label", "source_class", "reported_hgvsc", "reported_hgvsp",
        "mapping_method", "coordinate_mapping_status", "reference_match",
        "reference_validation_reason", "in_gene_panel",
    ]
    observations[observation_columns].to_csv(
        args.out_dir / "external_source_observations.tsv", sep="\t", index=False
    )
    benchmark.to_csv(args.out_dir / "external_unique_variant_manifest.tsv", sep="\t", index=False)
    conflicts.to_csv(args.out_dir / "external_label_conflicts.tsv", sep="\t", index=False)
    missing[["variant_id", "gene", "combined_label", "external_sources"]].to_csv(
        args.out_dir / "new_variants_requiring_annotation.tsv", sep="\t", index=False
    )
    expanded.to_csv(
        args.out_dir / "modeling_table_expanded_external_unannotated.tsv", sep="\t", index=False
    )

    mapped_counts = mapped.groupby("source").agg(
        observations=("variant_id", "size"),
        unique_variants=("variant_id", "nunique"),
        pathogenic=("external_binary_label", "sum"),
    )
    mapped_counts["benign"] = mapped_counts["observations"] - mapped_counts["pathogenic"]
    summary = {
        "inputs": {
            "matrix": str(args.matrix),
            "matrix_sha256": sha256(args.matrix),
            "share_csv": str(args.share_csv),
            "share_sha256": sha256(args.share_csv),
            "plos_xlsx": str(args.plos_xlsx),
            "plos_sha256": sha256(args.plos_xlsx),
            "clinvar_coordinate_resolver": str(args.clinvar),
            "reference_fasta": str(args.reference_fasta),
        },
        "source_counts": mapped_counts.astype(int).to_dict(orient="index"),
        "unmapped_or_reference_invalid_by_source": (
            observations.loc[
                observations["variant_id"].isna() | ~observations["reference_match"]
            ]["source"].value_counts().astype(int).to_dict()
        ),
        "unique_benchmark": {
            "variants": int(len(benchmark)),
            "pathogenic": int((benchmark["combined_label"] == "Pathogenic").sum()),
            "benign": int((benchmark["combined_label"] == "Benign").sum()),
            "source_label_conflicts_excluded": int(len(conflicts)),
            "already_in_matrix": int(benchmark["already_in_modeling_matrix"].sum()),
            "new_requiring_annotation": int((~benchmark["already_in_modeling_matrix"]).sum()),
        },
        "expanded_modeling_table": {
            "rows": int(len(expanded)),
            "train": int((expanded["split_source_heldout"] == "train").sum()),
            "validation": int((expanded["split_source_heldout"] == "validation").sum()),
            "external_expanded": int((expanded["split_source_heldout"] == "external_expanded").sum()),
            "quarantine_overlap_with_development": 0,
            "split_seed": args.seed,
        },
        "interpretation": [
            "SHaRe and PLOS are retrospective source-held-out validation cohorts, not prospective clinical validation.",
            "ClinVar supplies coordinates for PLOS records but does not supply their outcome label.",
            "Source-specific observations are retained separately; the combined manifest excludes label-discordant alleles.",
        ],
    }
    (args.out_dir / "build_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
