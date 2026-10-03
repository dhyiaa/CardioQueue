#!/usr/bin/env python3
"""Merge annotations for new external alleles into the frozen modeling schema."""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
from pathlib import Path

import pandas as pd
import pysam


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DIR = ROOT / "datasets/external_validation_candidates/expanded_external_2026"
MISSING = {"", ".", "NA", "N/A", "nan", "None", "Missing", pd.NA}

DBNSFP_MAP = {
    "dbnsfp_ref_aa": "dbnsfp_aaref", "dbnsfp_alt_aa": "dbnsfp_aaalt",
    "dbnsfp_protein_position_raw": "dbnsfp_aapos", "dbnsfp_gene_name_raw": "dbnsfp_genename",
    "dbnsfp_ensembl_transcript_raw": "dbnsfp_Ensembl_transcriptid",
    "dbnsfp_uniprot_accession_raw": "dbnsfp_Uniprot_acc",
    "dbnsfp_hgvsc_vep_raw": "dbnsfp_HGVSc_VEP", "dbnsfp_hgvsp_vep_raw": "dbnsfp_HGVSp_VEP",
    "dbnsfp_interpro_domain_raw": "dbnsfp_Interpro_domain", "sift_score": "dbnsfp_SIFT_score",
    "sift_pred": "dbnsfp_SIFT_pred", "polyphen2_hdiv_score": "dbnsfp_Polyphen2_HDIV_score",
    "polyphen2_hdiv_pred": "dbnsfp_Polyphen2_HDIV_pred", "revel_score": "dbnsfp_REVEL_score",
    "metalr_score": "dbnsfp_MetaLR_score", "metalr_pred": "dbnsfp_MetaLR_pred",
    "fathmm_xf_coding_score": "dbnsfp_fathmm-XF_coding_score",
    "fathmm_xf_coding_pred": "dbnsfp_fathmm-XF_coding_pred", "cadd_raw": "dbnsfp_CADD_raw",
    "cadd_phred": "dbnsfp_CADD_phred", "alphamissense_dbnsfp_score": "dbnsfp_AlphaMissense_score",
    "alphamissense_dbnsfp_pred": "dbnsfp_AlphaMissense_pred", "esm1b_score": "dbnsfp_ESM1b_score",
    "esm1b_pred": "dbnsfp_ESM1b_pred", "gerp_rs": "dbnsfp_GERP++_RS",
    "phylop100way_vertebrate": "dbnsfp_phyloP100way_vertebrate",
    "phastcons100way_vertebrate": "dbnsfp_phastCons100way_vertebrate",
    "gnomad41_joint_af_dbnsfp": "dbnsfp_gnomAD4.1_joint_AF",
    "gnomad41_joint_nhomalt_dbnsfp": "dbnsfp_gnomAD4.1_joint_nhomalt",
    "gnomad41_joint_popmax_af_dbnsfp": "dbnsfp_gnomAD4.1_joint_POPMAX_AF",
    "gnomad41_joint_popmax_nhomalt_dbnsfp": "dbnsfp_gnomAD4.1_joint_POPMAX_nhomalt",
    "dbnsfp_popmax_af": "dbnsfp_dbNSFP_POPMAX_AF", "dbnsfp_popmax_pop": "dbnsfp_dbNSFP_POPMAX_POP",
}
AM_MAP = {
    "alphamissense_direct_status": "alphamissense_status",
    "alphamissense_direct_missing_reason": "alphamissense_missing_reason",
    "alphamissense_direct_score": "am_pathogenicity", "alphamissense_direct_class": "am_class",
    "alphamissense_direct_uniprot_id": "am_uniprot_id",
    "alphamissense_direct_transcript_id": "am_transcript_id",
    "alphamissense_direct_protein_variant": "am_protein_variant",
}
SPLICE_FIELDS = ["SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL"]


def clean(value: object) -> object:
    if pd.isna(value):
        return pd.NA
    text = str(value).strip()
    return pd.NA if text in {"", ".", "NA", "N/A", "nan", "None", "Missing"} else text


def numeric_summary(value: object, mode: str) -> float | object:
    values = []
    for piece in str(value if pd.notna(value) else "").replace("|", ";").split(";"):
        try:
            values.append(float(piece))
        except ValueError:
            continue
    if not values:
        return pd.NA
    return min(values) if mode == "min" else max(values)


def first_numeric(value: object) -> float | object:
    return numeric_summary(value, "max")


def load_protein_helpers():
    path = ROOT / "datasets/feature_sources/protein_structure/scripts/build_protein_features.py"
    spec = importlib.util.spec_from_file_location("cardioqueue_protein_helpers", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def spliceai_for_variants(path: Path, variants: pd.DataFrame) -> dict[str, dict[str, object]]:
    vcf = pysam.VariantFile(str(path))
    output: dict[str, dict[str, object]] = {}
    for row in variants.itertuples(index=False):
        if len(str(row.ref)) != 1 or len(str(row.alt)) != 1:
            output[row.variant_id] = {"status": "not_applicable_non_snv"}
            continue
        try:
            records = vcf.fetch(str(row.chrom).removeprefix("chr"), int(row.pos) - 1, int(row.pos))
        except ValueError:
            output[row.variant_id] = {"status": "contig_not_found"}
            continue
        matches = []
        for record in records:
            if record.pos != int(row.pos) or record.ref != row.ref or row.alt not in record.alts:
                continue
            for entry in record.info.get("SpliceAI", []):
                pieces = str(entry).split("|")
                if len(pieces) == 10 and pieces[0] == row.alt:
                    matches.append(dict(zip(["ALLELE", *SPLICE_FIELDS], pieces)))
        if matches:
            best = max(matches, key=lambda item: max(float(item[field]) for field in SPLICE_FIELDS[1:5]))
            output[row.variant_id] = {"status": "ok", **best}
        else:
            output[row.variant_id] = {"status": "not_found"}
    vcf.close()
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_DIR / "modeling_table_expanded_external_unannotated.tsv")
    parser.add_argument("--new-variants", type=Path, default=DEFAULT_DIR / "new_variants_requiring_annotation.tsv")
    parser.add_argument("--dbnsfp", type=Path, default=DEFAULT_DIR / "annotations/dbnsfp.tsv")
    parser.add_argument("--alphamissense", type=Path, default=DEFAULT_DIR / "annotations/alphamissense.tsv")
    parser.add_argument("--vep", type=Path, default=DEFAULT_DIR / "annotations/vep_features.tsv")
    parser.add_argument("--spliceai", type=Path, default=ROOT / "datasets/feature_sources/insilico/raw/spliceai/downloads/spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz")
    parser.add_argument("--uniprot", type=Path, default=ROOT / "datasets/feature_sources/protein_structure/raw/uniprot_human_reviewed_features_2026_02.tsv")
    parser.add_argument("--output", type=Path, default=DEFAULT_DIR / "modeling_table_expanded_external_annotated.tsv")
    parser.add_argument("--summary", type=Path, default=DEFAULT_DIR / "feature_merge_summary.json")
    parser.add_argument("--structure", type=Path, default=DEFAULT_DIR / "annotations/structure_features.tsv")
    parser.add_argument("--structure-input", type=Path, default=DEFAULT_DIR / "annotations/structure_input.tsv")
    args = parser.parse_args()

    matrix = pd.read_csv(args.matrix, sep="\t", low_memory=False)
    numeric_columns = list(matrix.select_dtypes(include="number").columns)
    new_manifest = pd.read_csv(args.new_variants, sep="\t", dtype=str)
    dbnsfp = pd.read_csv(args.dbnsfp, sep="\t", dtype=str).set_index("variant_id", drop=False)
    am = pd.read_csv(args.alphamissense, sep="\t", dtype=str).set_index("variant_id", drop=False)
    vep = pd.read_csv(args.vep, sep="\t", low_memory=False).set_index("variant_id", drop=False)
    new_ids = set(new_manifest["variant_id"])
    target = matrix["variant_id"].isin(new_ids)
    if int(target.sum()) != len(new_ids):
        raise ValueError("Not every new external variant occurs once in the expanded matrix")

    coordinate_rows = matrix.loc[target, ["variant_id", "chrom", "pos", "ref", "alt"]].copy()
    splice = spliceai_for_variants(args.spliceai, coordinate_rows)
    protein_helpers = load_protein_helpers()
    uniprot_by_gene = protein_helpers.load_uniprot_features(args.uniprot)

    for idx in matrix.index[target]:
        variant_id = matrix.at[idx, "variant_id"]
        db = dbnsfp.loc[variant_id]
        alpha = am.loc[variant_id]
        vep_row = vep.loc[variant_id]
        for column in [column for column in dbnsfp.columns if column.startswith("dbnsfp_")]:
            if column in matrix.columns:
                matrix.at[idx, column] = clean(db[column])
        for output_column, input_column in DBNSFP_MAP.items():
            matrix.at[idx, output_column] = clean(db[input_column])
            missing_column = f"{output_column}_is_missing"
            if missing_column in matrix.columns:
                matrix.at[idx, missing_column] = pd.isna(clean(db[input_column]))
        matrix.at[idx, "sift_score_min"] = numeric_summary(db.get("dbnsfp_SIFT_score"), "min")
        matrix.at[idx, "revel_score_max"] = numeric_summary(db.get("dbnsfp_REVEL_score"), "max")
        matrix.at[idx, "alphamissense_dbnsfp_score_max"] = numeric_summary(db.get("dbnsfp_AlphaMissense_score"), "max")
        matrix.at[idx, "esm1b_score_max"] = numeric_summary(db.get("dbnsfp_ESM1b_score"), "max")
        for output_column, input_column in AM_MAP.items():
            matrix.at[idx, output_column] = clean(alpha[input_column])
        matrix.at[idx, "alphamissense_direct_score_is_missing"] = pd.isna(clean(alpha["am_pathogenicity"]))
        matrix.at[idx, "alphamissense_direct_class_is_missing"] = pd.isna(clean(alpha["am_class"]))
        for column in [column for column in vep.columns if column.startswith("vep_")]:
            if column in matrix.columns:
                matrix.at[idx, column] = vep_row[column]
        matrix.at[idx, "vep_status"] = vep_row["vep_annotation_status"]
        matrix.at[idx, "vep_missing_reason"] = vep_row["vep_annotation_missing_reason"]

        sp = splice[variant_id]
        matrix.at[idx, "spliceai_status"] = sp["status"]
        matrix.at[idx, "spliceai_missing_reason"] = "" if sp["status"] == "ok" else sp["status"]
        if sp["status"] == "ok":
            for field in SPLICE_FIELDS:
                column = f"vep_SpliceAI_pred_{field}"
                matrix.at[idx, column] = sp[field]

        af = first_numeric(db.get("dbnsfp_gnomAD4.1_joint_AF"))
        popmax = first_numeric(db.get("dbnsfp_gnomAD4.1_joint_POPMAX_AF"))
        nhom = first_numeric(db.get("dbnsfp_gnomAD4.1_joint_nhomalt"))
        if pd.notna(af):
            matrix.at[idx, "gnomad_final_status"] = "observed"
            matrix.at[idx, "gnomad_final_source"] = "dbNSFP_gnomAD4.1_joint"
            matrix.at[idx, "gnomad_final_missing_reason"] = ""
            matrix.at[idx, "gnomad_final_af"] = af
            matrix.at[idx, "gnomad_final_popmax_af"] = popmax
            matrix.at[idx, "gnomad_final_homozygote_count"] = nhom
            matrix.at[idx, "gnomad_final_af_is_missing"] = False
            matrix.at[idx, "gnomad_observed_flag"] = True
            matrix.at[idx, "gnomad_not_joined_flag"] = False
        else:
            matrix.at[idx, "gnomad_final_status"] = "not_joined"
            matrix.at[idx, "gnomad_final_source"] = "none"
            matrix.at[idx, "gnomad_final_missing_reason"] = "No exact dbNSFP gnomAD4.1 AF; absence not inferred"
            matrix.at[idx, "gnomad_final_af_is_missing"] = True
            matrix.at[idx, "gnomad_observed_flag"] = False
            matrix.at[idx, "gnomad_not_joined_flag"] = True
        matrix.at[idx, "gnomad_confirmed_absent_flag"] = False

        gene = str(matrix.at[idx, "primary_gene"])
        protein_variant = clean(alpha["am_protein_variant"])
        parsed = protein_helpers.parse_hgvs_p(
            f"p.{protein_variant}" if pd.notna(protein_variant) else "",
            f"p.{protein_variant}" if pd.notna(protein_variant) else "",
        )
        record = uniprot_by_gene.get(gene)
        if parsed["status"] == "ok" and record:
            position = int(parsed["position"])
            matrix.at[idx, "protein_feature_status"] = "ok"
            matrix.at[idx, "protein_feature_missing_reason"] = ""
            matrix.at[idx, "uniprot_accession"] = clean(alpha["am_uniprot_id"]) if pd.notna(clean(alpha["am_uniprot_id"])) else record["accession"]
            matrix.at[idx, "uniprot_length"] = record["length"]
            matrix.at[idx, "protein_position"] = position
            matrix.at[idx, "protein_ref_aa"] = parsed["ref_aa"]
            matrix.at[idx, "protein_alt_aa"] = parsed["alt_aa"]
            matrix.at[idx, "protein_variant_type"] = parsed["variant_type"]
            matrix.at[idx, "protein_position_normalized"] = position / record["length"] if record["length"] else pd.NA
            any_hit = False
            for prefix in ["uniprot_domain", "uniprot_region", "uniprot_motif", "uniprot_binding_site", "uniprot_active_site"]:
                hit, names = protein_helpers.overlap_features(record, position, prefix)
                matrix.at[idx, f"{prefix}_hit"] = hit
                matrix.at[idx, f"{prefix}_names"] = names
                any_hit = any_hit or hit
            matrix.at[idx, "uniprot_any_functional_feature_hit"] = any_hit
            matrix.at[idx, "alphafold_plddt_status"] = "pending_structure_lookup"
        else:
            matrix.at[idx, "protein_feature_status"] = "not_mapped"
            matrix.at[idx, "protein_feature_missing_reason"] = "No validated missense protein mapping"
            matrix.at[idx, "alphafold_plddt_status"] = "not_queryable"
            matrix.at[idx, "alphafold_plddt_missing_reason"] = "Protein position unavailable"
        matrix.at[idx, "foldx_ddg_status"] = "not_available"
        matrix.at[idx, "foldx_ddg_missing_reason"] = "No precomputed external-variant FoldX result"
        matrix.at[idx, "foldx_ddg_is_missing"] = True
        matrix.at[idx, "foldx_ddg_ref_mismatch_flag"] = False

    structure_joined = 0
    if args.structure.exists():
        structure = pd.read_csv(args.structure, sep="\t", low_memory=False).set_index("variant_id")
        for idx in matrix.index[target]:
            variant_id = matrix.at[idx, "variant_id"]
            if variant_id not in structure.index:
                continue
            row = structure.loc[variant_id]
            if bool(row.get("structure_available", False)):
                structure_joined += 1
                matrix.at[idx, "alphafold_plddt_status"] = "ok"
                matrix.at[idx, "alphafold_plddt_missing_reason"] = ""
                matrix.at[idx, "alphafold_residue_plddt"] = row.get("structure_plddt")
                matrix.at[idx, "alphafold_plddt_bin"] = row.get("structure_plddt_bin")
                matrix.at[idx, "dssp_status"] = "ok" if pd.notna(row.get("structure_dssp_code")) else "not_returned"
                matrix.at[idx, "dssp_missing_reason"] = "" if matrix.at[idx, "dssp_status"] == "ok" else "DSSP did not return residue"
                matrix.at[idx, "freesasa_status"] = "ok" if pd.notna(row.get("structure_freesasa_total")) else "not_returned"
                matrix.at[idx, "freesasa_missing_reason"] = "" if matrix.at[idx, "freesasa_status"] == "ok" else "FreeSASA did not return residue"
            elif matrix.at[idx, "alphafold_plddt_status"] == "pending_structure_lookup":
                matrix.at[idx, "alphafold_plddt_status"] = "not_mapped"
                matrix.at[idx, "alphafold_plddt_missing_reason"] = "No validated AlphaFold residue mapping"

    for column in numeric_columns:
        matrix[column] = pd.to_numeric(matrix[column], errors="coerce")

    matrix.loc[target, [
        "variant_id", "split_source_heldout", "training_slice_primary_binary",
        "alphamissense_direct_uniprot_id", "alphamissense_direct_protein_variant",
    ]].to_csv(args.structure_input, sep="\t", index=False)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(args.output, sep="\t", index=False)
    new_rows = matrix[target]
    summary = {
        "new_external_variants": int(len(new_rows)),
        "dbnsfp_ok": int((new_rows["dbnsfp_status"] == "ok").sum()),
        "alphamissense_ok": int((new_rows["alphamissense_direct_status"] == "ok").sum()),
        "vep_ok": int((new_rows["vep_annotation_status"] == "ok").sum()),
        "spliceai_ok": int((new_rows["spliceai_status"] == "ok").sum()),
        "protein_mapped": int((new_rows["protein_feature_status"] == "ok").sum()),
        "alphafold_structure_joined": int(structure_joined),
        "gnomad_observed_from_dbnsfp": int((new_rows["gnomad_final_status"] == "observed").sum()),
        "foldx_new_external_variants": 0,
        "notes": [
            "Missing predictors remain explicit; confirmed population absence is never inferred from a failed join.",
            "FoldX remains unavailable for newly introduced external alleles and is represented as missing.",
        ],
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
