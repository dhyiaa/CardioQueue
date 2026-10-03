#!/usr/bin/env python3
"""Join local feature layers onto the clean combined registry."""

from __future__ import annotations

import csv
import json
import argparse
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]

CLEAN_REGISTRY = ROOT / "datasets/variant_registry/interim/clean_combined_variant_registry.tsv"
DBNSFP_FEATURES = ROOT / "datasets/feature_sources/insilico/interim/dbnsfp_features_master_variant_registry.tsv"
PROTEIN_FEATURES = ROOT / "datasets/feature_sources/protein_structure/interim/protein_features.tsv"

OUT_DIR = ROOT / "datasets/variant_registry/interim"
FEATURE_REGISTRY_OUT = OUT_DIR / "clean_combined_local_feature_registry.tsv"
SUMMARY_OUT = OUT_DIR / "clean_combined_local_feature_registry.summary.json"

DBNSFP_COLUMNS = [
    "dbnsfp_status",
    "dbnsfp_missing_reason",
    "dbnsfp_aaref",
    "dbnsfp_aaalt",
    "dbnsfp_aapos",
    "dbnsfp_genename",
    "dbnsfp_Ensembl_transcriptid",
    "dbnsfp_Uniprot_acc",
    "dbnsfp_HGVSc_VEP",
    "dbnsfp_HGVSp_VEP",
    "dbnsfp_Interpro_domain",
    "dbnsfp_SIFT_score",
    "dbnsfp_SIFT_pred",
    "dbnsfp_Polyphen2_HDIV_score",
    "dbnsfp_Polyphen2_HDIV_pred",
    "dbnsfp_REVEL_score",
    "dbnsfp_MetaLR_score",
    "dbnsfp_MetaLR_pred",
    "dbnsfp_fathmm-XF_coding_score",
    "dbnsfp_fathmm-XF_coding_pred",
    "dbnsfp_CADD_raw",
    "dbnsfp_CADD_phred",
    "dbnsfp_AlphaMissense_score",
    "dbnsfp_AlphaMissense_pred",
    "dbnsfp_ESM1b_score",
    "dbnsfp_ESM1b_pred",
    "dbnsfp_GERP++_RS",
    "dbnsfp_phyloP100way_vertebrate",
    "dbnsfp_phastCons100way_vertebrate",
    "dbnsfp_gnomAD4.1_joint_AF",
    "dbnsfp_gnomAD4.1_joint_nhomalt",
    "dbnsfp_gnomAD4.1_joint_POPMAX_AF",
    "dbnsfp_gnomAD4.1_joint_POPMAX_nhomalt",
    "dbnsfp_dbNSFP_POPMAX_AF",
    "dbnsfp_dbNSFP_POPMAX_POP",
]

PROTEIN_COLUMNS = [
    "protein_feature_status",
    "protein_feature_missing_reason",
    "uniprot_accession",
    "uniprot_length",
    "protein_position",
    "protein_ref_aa",
    "protein_alt_aa",
    "protein_variant_type",
    "protein_position_normalized",
    "uniprot_domain_hit",
    "uniprot_domain_names",
    "uniprot_region_hit",
    "uniprot_region_names",
    "uniprot_motif_hit",
    "uniprot_motif_names",
    "uniprot_binding_site_hit",
    "uniprot_binding_site_names",
    "uniprot_active_site_hit",
    "uniprot_active_site_names",
    "uniprot_any_functional_feature_hit",
    "alphafold_plddt_status",
    "alphafold_plddt_missing_reason",
    "alphafold_residue_plddt",
    "alphafold_plddt_bin",
    "alphafold_fragment_count",
]


def rel(path: Path) -> str:
    path = path.resolve()
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def variant_key(row: dict[str, str]) -> str:
    return f"{row.get('chrom', '').strip()}-{row.get('pos', '').strip()}-{row.get('ref', '').strip().upper()}-{row.get('alt', '').strip().upper()}"


def load_dbnsfp(path: Path) -> dict[str, dict[str, str]]:
    features = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            features[row["variant_id"]] = {col: row.get(col, "") for col in DBNSFP_COLUMNS}
    return features


def protein_rank(row: dict[str, str]) -> tuple[int, int]:
    status = row.get("protein_feature_status", "")
    af_status = row.get("alphafold_plddt_status", "")
    return (0 if status == "ok" else 1, 0 if af_status == "ok" else 1)


def load_protein_features(path: Path) -> dict[str, dict[str, str]]:
    best: dict[str, dict[str, str]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = variant_key(row)
            if "----" in key or not row.get("chrom") or not row.get("pos"):
                continue
            current = best.get(key)
            if current is None or protein_rank(row) < protein_rank(current):
                best[key] = {col: row.get(col, "") for col in PROTEIN_COLUMNS}
    return best


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clean-registry", type=Path, default=CLEAN_REGISTRY)
    parser.add_argument("--dbnsfp-features", type=Path, default=DBNSFP_FEATURES)
    parser.add_argument("--protein-features", type=Path, default=PROTEIN_FEATURES)
    parser.add_argument("--out", type=Path, default=FEATURE_REGISTRY_OUT)
    parser.add_argument("--summary", type=Path, default=SUMMARY_OUT)
    args = parser.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    dbnsfp = load_dbnsfp(args.dbnsfp_features)
    protein = load_protein_features(args.protein_features)
    counts: Counter = Counter()

    with args.clean_registry.open(newline="") as in_handle, args.out.open("w", newline="") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or []) + DBNSFP_COLUMNS + PROTEIN_COLUMNS
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in reader:
            counts["rows"] += 1
            key = row["variant_id"]
            db = dbnsfp.get(key)
            if db:
                counts[f"dbnsfp_{db.get('dbnsfp_status') or 'unknown'}"] += 1
                row.update(db)
            else:
                counts["dbnsfp_not_in_feature_table"] += 1
                row.update({col: "" for col in DBNSFP_COLUMNS})
                row["dbnsfp_status"] = "not_joined"
                row["dbnsfp_missing_reason"] = "variant not present in dbNSFP master-registry feature table"

            prot = protein.get(key)
            if prot:
                counts[f"protein_{prot.get('protein_feature_status') or 'unknown'}"] += 1
                counts[f"alphafold_{prot.get('alphafold_plddt_status') or 'unknown'}"] += 1
                row.update(prot)
            else:
                counts["protein_no_source_row_feature"] += 1
                row.update({col: "" for col in PROTEIN_COLUMNS})
                row["protein_feature_status"] = "not_available"
                row["protein_feature_missing_reason"] = (
                    "no local protein feature row; current protein table covers HiRO/eMERGE/CardioBoost only"
                )
                row["alphafold_plddt_status"] = "not_available"
                row["alphafold_plddt_missing_reason"] = (
                    "no local protein feature row; current protein table covers HiRO/eMERGE/CardioBoost only"
                )

            writer.writerow(row)

    summary = {
        "inputs": {
            "clean_registry": rel(args.clean_registry),
            "dbnsfp_features": rel(args.dbnsfp_features),
            "protein_features": rel(args.protein_features),
        },
        "output": rel(args.out),
        "counts": dict(sorted(counts.items())),
        "notes": [
            "This is a local-only feature registry: dbNSFP plus precomputed protein/AlphaFold/UniProt features.",
            "gnomAD API, ClinGen API, CADD API, VEP/SpliceAI, DSSP/FreeSASA batch, and FoldXPro batch are not included yet.",
            "Protein features are currently available only for source rows already parsed from HiRO/eMERGE/CardioBoost.",
        ],
    }
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
