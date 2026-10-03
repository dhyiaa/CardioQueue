#!/usr/bin/env python3
"""Build basic protein features from local UniProt and AlphaFold sources."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import tarfile
from pathlib import Path
from typing import Any


AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Sec": "U", "Pyl": "O",
}


INPUTS = {
    "hiro": Path("datasets/hiro/full_dataset/model_inputs/ml_baseline_features.csv"),
    "emerge": Path("datasets/emerge/full_arrhythmia_gene_dataset/data/emerge_arrhythmia_ml_baseline_features.csv"),
    "cardioboost": Path(
        "datasets/cardioboost/public_dataset/processed_for_cardiogenetics/model_inputs/"
        "reannotated_ml_baseline_features.csv"
    ),
}

OUTPUT_FIELDS = [
    "source_dataset",
    "variant_uid",
    "variant_key",
    "target_3class",
    "gene",
    "gene_input",
    "hgvs_p",
    "hgvs_p_one_letter",
    "inferred_consequence",
    "chrom",
    "pos",
    "ref",
    "alt",
    "protein_feature_status",
    "protein_feature_missing_reason",
    "uniprot_accession",
    "uniprot_gene_names",
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


FEATURE_SPECS = [
    ("Domain [FT]", "DOMAIN", "uniprot_domain"),
    ("Region", "REGION", "uniprot_region"),
    ("Motif", "MOTIF", "uniprot_motif"),
    ("Binding site", "BINDING", "uniprot_binding_site"),
    ("Active site", "ACT_SITE", "uniprot_active_site"),
]


def clean_missing(value: Any) -> str:
    text = "" if value is None else str(value).strip()
    if text in {"", "Missing", "NA", "NaN", "nan", "None"}:
        return ""
    return text


def normalize_gene(gene: str) -> str:
    gene = clean_missing(gene).upper()
    if ";" in gene:
        parts = [part for part in gene.split(";") if part and not part.startswith("FPGT-")]
        if parts:
            return parts[-1]
    return gene


def parse_hgvs_p(hgvs_p: str, hgvs_p_one_letter: str) -> dict[str, Any]:
    one = clean_missing(hgvs_p_one_letter)
    three = clean_missing(hgvs_p)
    for text in [one, three]:
        if not text:
            continue
        text = text.strip()
        match = re.search(r"p\.([A-Z*])(\d+)([A-Z*])$", text)
        if match:
            return {
                "status": "ok",
                "position": int(match.group(2)),
                "ref_aa": match.group(1),
                "alt_aa": match.group(3),
                "variant_type": "substitution",
            }
        match = re.search(r"p\.([A-Z])(\d+)fs", text, flags=re.IGNORECASE)
        if match:
            return {
                "status": "ok",
                "position": int(match.group(2)),
                "ref_aa": match.group(1).upper(),
                "alt_aa": "fs",
                "variant_type": "frameshift",
            }
        match = re.search(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter)$", text)
        if match:
            return {
                "status": "ok",
                "position": int(match.group(2)),
                "ref_aa": AA3_TO_1.get(match.group(1), ""),
                "alt_aa": AA3_TO_1.get(match.group(3), ""),
                "variant_type": "substitution",
            }
    return {
        "status": "not_mapped",
        "position": "",
        "ref_aa": "",
        "alt_aa": "",
        "variant_type": "",
    }


def parse_note(entry: str) -> str:
    note = re.search(r'/note="([^"]+)"', entry)
    if note:
        return note.group(1)
    ligand = re.search(r'/ligand="([^"]+)"', entry)
    if ligand:
        return f"ligand:{ligand.group(1)}"
    return ""


def parse_feature_text(text: str, feature_code: str) -> list[dict[str, Any]]:
    if not text:
        return []
    pattern = re.compile(rf"({feature_code})\s+(\d+)(?:\.\.(\d+))?;([^;]*(?:;(?!\s+[A-Z_]+\s+\d)[^;]*)*)")
    features = []
    for match in pattern.finditer(text):
        start = int(match.group(2))
        end = int(match.group(3) or match.group(2))
        features.append({
            "start": start,
            "end": end,
            "name": parse_note(match.group(0)),
        })
    return features


def load_uniprot_features(path: Path) -> dict[str, dict[str, Any]]:
    by_gene: dict[str, dict[str, Any]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            genes = [gene.upper() for gene in row.get("Gene Names", "").split()]
            if not genes:
                continue
            record = {
                "accession": row.get("Entry", ""),
                "gene_names": row.get("Gene Names", ""),
                "length": int(row.get("Length") or 0),
                "features": {},
            }
            for col, code, prefix in FEATURE_SPECS:
                record["features"][prefix] = parse_feature_text(row.get(col, ""), code)
            for gene in genes:
                by_gene.setdefault(gene, record)
    return by_gene


def load_alphafold_manifest(path: Path) -> dict[str, list[str]]:
    by_accession: dict[str, list[str]] = {}
    pattern = re.compile(r"^AF-(?P<acc>.+)-F(?P<frag>\d+)-model_v\d+\.pdb\.gz$")
    with path.open() as handle:
        for raw_line in handle:
            name = raw_line.strip()
            match = pattern.match(name)
            if not match:
                continue
            by_accession.setdefault(match.group("acc"), []).append(name)
    for names in by_accession.values():
        names.sort(key=lambda name: int(re.search(r"-F(\d+)-", name).group(1)))
    return by_accession


def load_plddt_for_member(tar_path: Path, member_name: str) -> dict[int, float]:
    scores: dict[int, float] = {}
    with tarfile.open(tar_path, "r") as archive:
        member = archive.getmember(member_name)
        payload = archive.extractfile(member).read()
    text = gzip.decompress(payload).decode("utf-8", "ignore")
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        try:
            residue = int(line[22:26])
            b_factor = float(line[60:66])
        except ValueError:
            continue
        scores[residue] = b_factor
    return scores


def plddt_bin(score: Any) -> str:
    if score == "":
        return ""
    score = float(score)
    if score < 50:
        return "very_low"
    if score < 70:
        return "low"
    if score < 90:
        return "confident"
    return "very_high"


def overlap_features(record: dict[str, Any], position: int, prefix: str) -> tuple[bool, str]:
    hits = [
        feature.get("name", "")
        for feature in record["features"].get(prefix, [])
        if feature["start"] <= position <= feature["end"]
    ]
    names = "|".join(hit for hit in hits if hit)
    return bool(hits), names


def read_source_rows() -> list[dict[str, str]]:
    rows = []
    for source, path in INPUTS.items():
        with path.open(newline="", errors="replace") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                rows.append({
                    "source_dataset": source,
                    "variant_uid": row.get("variant_uid", ""),
                    "variant_key": row.get("variant_key", ""),
                    "target_3class": row.get("target_3class", ""),
                    "gene": row.get("gene", ""),
                    "hgvs_p": row.get("hgvs_p", ""),
                    "hgvs_p_one_letter": row.get("hgvs_p_one_letter", ""),
                    "inferred_consequence": row.get("inferred_consequence", ""),
                    "chrom": row.get("chrom", ""),
                    "pos": row.get("pos", ""),
                    "ref": row.get("ref", ""),
                    "alt": row.get("alt", ""),
                })
    return rows


def build_features(args: argparse.Namespace) -> dict[str, Any]:
    uniprot = load_uniprot_features(Path(args.uniprot_features))
    af_manifest = load_alphafold_manifest(Path(args.alphafold_manifest))
    af_tar = Path(args.alphafold_tar)
    plddt_cache: dict[str, dict[int, float]] = {}

    output_rows = []
    status_counts: dict[str, int] = {}
    plddt_status_counts: dict[str, int] = {}

    for source_row in read_source_rows():
        row = {field: "" for field in OUTPUT_FIELDS}
        row.update(source_row)
        gene = normalize_gene(source_row["gene"])
        row["gene_input"] = source_row["gene"]
        row["gene"] = gene

        parsed = parse_hgvs_p(source_row["hgvs_p"], source_row["hgvs_p_one_letter"])
        if parsed["status"] != "ok":
            row["protein_feature_status"] = "not_mapped"
            row["protein_feature_missing_reason"] = "No parseable protein position in hgvs_p/hgvs_p_one_letter."
            row["alphafold_plddt_status"] = "not_queryable"
            row["alphafold_plddt_missing_reason"] = "Protein position unavailable."
            output_rows.append(row)
            status_counts[row["protein_feature_status"]] = status_counts.get(row["protein_feature_status"], 0) + 1
            plddt_status_counts[row["alphafold_plddt_status"]] = plddt_status_counts.get(row["alphafold_plddt_status"], 0) + 1
            continue

        position = int(parsed["position"])
        row["protein_position"] = position
        row["protein_ref_aa"] = parsed["ref_aa"]
        row["protein_alt_aa"] = parsed["alt_aa"]
        row["protein_variant_type"] = parsed["variant_type"]

        record = uniprot.get(gene)
        if not record:
            row["protein_feature_status"] = "no_uniprot_mapping"
            row["protein_feature_missing_reason"] = "Gene not mapped to reviewed human UniProt entry."
            row["alphafold_plddt_status"] = "not_queryable"
            row["alphafold_plddt_missing_reason"] = "UniProt accession unavailable."
            output_rows.append(row)
            status_counts[row["protein_feature_status"]] = status_counts.get(row["protein_feature_status"], 0) + 1
            plddt_status_counts[row["alphafold_plddt_status"]] = plddt_status_counts.get(row["alphafold_plddt_status"], 0) + 1
            continue

        row["protein_feature_status"] = "ok"
        row["uniprot_accession"] = record["accession"]
        row["uniprot_gene_names"] = record["gene_names"]
        row["uniprot_length"] = record["length"]
        if record["length"]:
            row["protein_position_normalized"] = round(position / record["length"], 6)

        any_hit = False
        for _, _, prefix in FEATURE_SPECS:
            hit, names = overlap_features(record, position, prefix)
            row[f"{prefix}_hit"] = hit
            row[f"{prefix}_names"] = names
            any_hit = any_hit or hit
        row["uniprot_any_functional_feature_hit"] = any_hit

        names = af_manifest.get(record["accession"], [])
        row["alphafold_fragment_count"] = len(names)
        if not names:
            row["alphafold_plddt_status"] = "not_found"
            row["alphafold_plddt_missing_reason"] = "Accession not present in AlphaFold human bulk manifest."
        elif len(names) > 1 and position > 1400:
            row["alphafold_plddt_status"] = "fragment_mapping_unresolved"
            row["alphafold_plddt_missing_reason"] = "AlphaFold model has multiple fragments; non-F1 global residue mapping not implemented."
        else:
            member_name = names[0]
            if member_name not in plddt_cache:
                plddt_cache[member_name] = load_plddt_for_member(af_tar, member_name)
            score = plddt_cache[member_name].get(position, "")
            if score == "":
                row["alphafold_plddt_status"] = "position_not_found"
                row["alphafold_plddt_missing_reason"] = "Residue position absent from selected AlphaFold fragment."
            else:
                row["alphafold_plddt_status"] = "ok"
                row["alphafold_residue_plddt"] = score
                row["alphafold_plddt_bin"] = plddt_bin(score)

        output_rows.append(row)
        status_counts[row["protein_feature_status"]] = status_counts.get(row["protein_feature_status"], 0) + 1
        plddt_status_counts[row["alphafold_plddt_status"]] = plddt_status_counts.get(row["alphafold_plddt_status"], 0) + 1

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    summary = {
        "input_rows": len(output_rows),
        "output": str(output_path),
        "protein_feature_status_counts": status_counts,
        "alphafold_plddt_status_counts": plddt_status_counts,
        "sources": {source: str(path) for source, path in INPUTS.items()},
    }
    Path(args.summary).write_text(json.dumps(summary, indent=2))
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--uniprot-features",
        default="datasets/feature_sources/protein_structure/raw/uniprot_human_reviewed_features_2026_02.tsv",
    )
    parser.add_argument(
        "--alphafold-manifest",
        default="datasets/feature_sources/protein_structure/raw/alphafold/UP000005640_9606_HUMAN_v6.manifest.txt",
    )
    parser.add_argument(
        "--alphafold-tar",
        default="datasets/feature_sources/protein_structure/raw/alphafold/UP000005640_9606_HUMAN_v6.tar",
    )
    parser.add_argument(
        "--output",
        default="datasets/feature_sources/protein_structure/interim/protein_features.tsv",
    )
    parser.add_argument(
        "--summary",
        default="datasets/feature_sources/protein_structure/interim/protein_features.summary.json",
    )
    args = parser.parse_args()
    print(json.dumps(build_features(args), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

