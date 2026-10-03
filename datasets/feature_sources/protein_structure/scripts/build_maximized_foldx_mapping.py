#!/usr/bin/env python3
"""Build an outcome-blind, maximal defensible FoldX mapping for the accepted cohort.

Mappings are accepted in descending order of specificity: direct AlphaMissense
protein mapping, aligned dbNSFP transcript fields, then an exact VEP amino-acid
substitution on the gene's reviewed canonical UniProt protein. Every accepted
mapping must agree with both the reviewed UniProt sequence and the selected
AlphaFold v6 fragment.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import tarfile
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
DEFAULT_TABLE = ROOT / "datasets/external_validation_candidates/expanded_external_2026/modeling_table_expanded_external_annotated.tsv"
DEFAULT_FASTA = ROOT / "datasets/feature_sources/protein_structure/raw/uniprot_human_reviewed_2026_02.fasta"
DEFAULT_MANIFEST = ROOT / "datasets/feature_sources/protein_structure/raw/alphafold/UP000005640_9606_HUMAN_v6.manifest.txt"
DEFAULT_TAR = ROOT / "datasets/feature_sources/protein_structure/raw/alphafold/UP000005640_9606_HUMAN_v6.tar"
DEFAULT_VEP_TAB = ROOT / "datasets/feature_sources/insilico/interim/vep/final_modeling_table_local_features_gnomad_foldx.hiro_emerge_rescued_vep_hgvs_spliceai.tsv"
DEFAULT_VEP_REST_JSON = ROOT / "datasets/external_validation_candidates/expanded_external_2026/annotations/vep_rest_raw.json"
DEFAULT_WORK = ROOT / "datasets/feature_sources/protein_structure/interim/registry_structure"
DEFAULT_AUDIT = ROOT / "datasets/feature_sources/protein_structure/interim/maximized_foldx_mapping_audit.tsv"
DEFAULT_FEATURES = ROOT / "datasets/feature_sources/protein_structure/interim/maximized_foldx_structure_features.tsv"
DEFAULT_SUMMARY = ROOT / "datasets/feature_sources/protein_structure/interim/maximized_foldx_mapping.summary.json"

AA = "ACDEFGHIKLMNPQRSTVWY"
VARIANT_RE = re.compile(rf"^([{AA}])(\d+)([{AA}])$")
VEP_AA_RE = re.compile(rf"^([{AA}])/([{AA}])$")
MEMBER_RE = re.compile(r"^AF-(?P<accession>.+)-F(?P<fragment>\d+)-model_v(?P<version>\d+)\.pdb\.gz$")
FRAGMENT_OFFSET = 200
AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def truthy(value: object) -> bool:
    return str(value).lower() in {"true", "1", "yes"}


def clean(value: object) -> str:
    return "" if pd.isna(value) else str(value).strip()


def canonical_accession(value: object) -> str:
    accession = clean(value)
    if not accession or accession == ".":
        return ""
    return accession.split("-")[0]


def split_field(value: object) -> list[str]:
    text = clean(value)
    return text.split(";") if text else []


def load_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    accession = ""
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if accession:
                sequences[accession] = "".join(chunks)
            parts = line.split("|")
            accession = parts[1] if len(parts) > 1 else line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if accession:
        sequences[accession] = "".join(chunks)
    return sequences


def load_manifest(path: Path) -> dict[str, list[tuple[int, str]]]:
    members: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for name in path.read_text().splitlines():
        match = MEMBER_RE.match(name.strip())
        if match:
            members[match.group("accession")].append((int(match.group("fragment")), name.strip()))
    return {accession: sorted(values) for accession, values in members.items()}


def parse_pdb(text: str) -> tuple[dict[int, float], dict[int, str]]:
    plddt: dict[int, float] = {}
    residues: dict[int, str] = {}
    for line in text.splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        try:
            position = int(line[22:26])
            plddt[position] = float(line[60:66])
        except ValueError:
            continue
        residues[position] = AA3_TO_1.get(line[17:20].strip(), "")
    return plddt, residues


def valid_sequence_mapping(accession: str, ref: str, position: int, sequences: dict[str, str]) -> bool:
    sequence = sequences.get(accession, "")
    return bool(sequence) and 1 <= position <= len(sequence) and sequence[position - 1] == ref


def direct_candidate(row: object, sequences: dict[str, str]) -> tuple[str, str, int, str] | None:
    accession = canonical_accession(row.alphamissense_direct_uniprot_id)
    match = VARIANT_RE.match(clean(row.alphamissense_direct_protein_variant))
    if not accession or not match:
        return None
    ref, position, alt = match.group(1), int(match.group(2)), match.group(3)
    if ref == alt or not valid_sequence_mapping(accession, ref, position, sequences):
        return None
    return accession, ref, position, alt


def dbnsfp_candidate(row: object, sequences: dict[str, str]) -> tuple[str, str, int, str] | None:
    ref, alt = clean(row.dbnsfp_aaref), clean(row.dbnsfp_aaalt)
    if len(ref) != 1 or ref not in AA or len(alt) != 1 or alt not in AA or ref == alt:
        return None
    accessions = split_field(row.dbnsfp_Uniprot_acc)
    positions = split_field(row.dbnsfp_aapos)
    if not accessions or len(accessions) != len(positions):
        return None
    candidates = set()
    for raw_accession, raw_position in zip(accessions, positions):
        accession = canonical_accession(raw_accession)
        try:
            position = int(raw_position)
        except ValueError:
            continue
        if accession and valid_sequence_mapping(accession, ref, position, sequences):
            candidates.add((accession, ref, position, alt))
    return next(iter(candidates)) if len(candidates) == 1 else None


def build_gene_accessions(frame: pd.DataFrame, sequences: dict[str, str]) -> dict[str, str]:
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in frame.itertuples(index=False):
        candidate = direct_candidate(row, sequences)
        gene = clean(row.primary_gene)
        if candidate and gene:
            counts[gene][candidate[0]] += 1
    result = {}
    for gene, accession_counts in counts.items():
        ordered = accession_counts.most_common()
        if ordered and (len(ordered) == 1 or ordered[0][1] > ordered[1][1]):
            result[gene] = ordered[0][0]
    return result


def vep_candidate(
    row: object, sequences: dict[str, str], gene_accessions: dict[str, str]
) -> tuple[str, str, int, str] | None:
    aa_match = VEP_AA_RE.match(clean(row.vep_amino_acids))
    position_text = clean(row.vep_protein_position)
    if not aa_match or not re.fullmatch(r"\d+", position_text):
        return None
    ref, alt = aa_match.group(1), aa_match.group(2)
    position = int(position_text)
    accession = gene_accessions.get(clean(row.primary_gene), "")
    if ref == alt or not accession or not valid_sequence_mapping(accession, ref, position, sequences):
        return None
    return accession, ref, position, alt


def load_transcript_vep_rescues(
    path: Path,
    targets: dict[str, str],
    sequences: dict[str, str],
    gene_accessions: dict[str, str],
) -> dict[str, tuple[str, str, int, str]]:
    """Recover exact missense consequences hidden by one-row VEP collapsing."""
    if not path.exists() or not targets:
        return {}
    header: list[str] | None = None
    candidates: dict[str, set[tuple[str, str, int, str]]] = defaultdict(set)
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("#Uploaded_variation"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                break
        if header is None:
            raise ValueError(f"VEP tab file lacks #Uploaded_variation header: {path}")
        index = {name: header.index(name) for name in [
            "Uploaded_variation", "SYMBOL", "Consequence", "Protein_position", "Amino_acids"
        ]}
        reader = csv.reader(handle, delimiter="\t")
        for fields in reader:
            if len(fields) < len(header):
                continue
            variant_id = fields[index["Uploaded_variation"]]
            gene = targets.get(variant_id)
            if not gene or fields[index["SYMBOL"]] != gene:
                continue
            if "missense_variant" not in fields[index["Consequence"]]:
                continue
            aa_match = VEP_AA_RE.match(fields[index["Amino_acids"]])
            position_text = fields[index["Protein_position"]]
            if not aa_match or not re.fullmatch(r"\d+", position_text):
                continue
            ref, alt = aa_match.group(1), aa_match.group(2)
            position = int(position_text)
            accession = gene_accessions.get(gene, "")
            if ref != alt and accession and valid_sequence_mapping(accession, ref, position, sequences):
                candidates[variant_id].add((accession, ref, position, alt))
    return {variant_id: next(iter(values)) for variant_id, values in candidates.items() if len(values) == 1}


def load_rest_vep_rescues(
    path: Path,
    targets: dict[str, str],
    sequences: dict[str, str],
    gene_accessions: dict[str, str],
) -> dict[str, tuple[str, str, int, str]]:
    if not path.exists() or not targets:
        return {}
    raw = json.loads(path.read_text())
    candidates: dict[str, set[tuple[str, str, int, str]]] = defaultdict(set)
    for variant_id, result in raw.items():
        gene = targets.get(variant_id)
        if not gene:
            continue
        for transcript in result.get("transcript_consequences", []) or []:
            if transcript.get("gene_symbol") != gene:
                continue
            if "missense_variant" not in (transcript.get("consequence_terms") or []):
                continue
            aa_match = VEP_AA_RE.match(clean(transcript.get("amino_acids")))
            position_value = transcript.get("protein_start")
            if not aa_match or not isinstance(position_value, int):
                continue
            ref, alt = aa_match.group(1), aa_match.group(2)
            accession = gene_accessions.get(gene, "")
            if ref != alt and accession and valid_sequence_mapping(accession, ref, position_value, sequences):
                candidates[variant_id].add((accession, ref, position_value, alt))
    return {variant_id: next(iter(values)) for variant_id, values in candidates.items() if len(values) == 1}


def select_fragment(
    candidates: list[tuple[int, str]], global_position: int, fragment_lengths: dict[str, int]
) -> tuple[int, str, int, int] | None:
    valid = []
    for fragment, member in candidates:
        local_position = global_position - FRAGMENT_OFFSET * (fragment - 1)
        length = fragment_lengths.get(member, 0)
        if 1 <= local_position <= length:
            edge_distance = min(local_position - 1, length - local_position)
            valid.append((edge_distance, -fragment, fragment, member, local_position))
    if not valid:
        return None
    edge_distance, _, fragment, member, local_position = max(valid)
    return fragment, member, local_position, edge_distance


def plddt_tier(value: float) -> str:
    if value < 50:
        return "very_low_lt50"
    if value < 70:
        return "low_50_69"
    if value < 90:
        return "confident_70_89"
    return "very_high_ge90"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", default=str(DEFAULT_TABLE))
    parser.add_argument("--fasta", default=str(DEFAULT_FASTA))
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    parser.add_argument("--alphafold-tar", default=str(DEFAULT_TAR))
    parser.add_argument("--vep-tab", default=str(DEFAULT_VEP_TAB))
    parser.add_argument("--vep-rest-json", default=str(DEFAULT_VEP_REST_JSON))
    parser.add_argument("--work-dir", default=str(DEFAULT_WORK))
    parser.add_argument("--audit-output", default=str(DEFAULT_AUDIT))
    parser.add_argument("--features-output", default=str(DEFAULT_FEATURES))
    parser.add_argument("--summary", default=str(DEFAULT_SUMMARY))
    args = parser.parse_args()

    columns = [
        "variant_id", "primary_gene", "split_source_heldout", "training_slice_primary_binary",
        "primary_binary_label", "alphamissense_direct_uniprot_id",
        "alphamissense_direct_protein_variant", "dbnsfp_aaref", "dbnsfp_aaalt",
        "dbnsfp_aapos", "dbnsfp_Uniprot_acc", "vep_protein_position", "vep_amino_acids",
        "vep_any_missense",
    ]
    frame = pd.read_csv(args.table, sep="\t", usecols=columns, low_memory=False)
    if not frame["variant_id"].is_unique:
        raise ValueError("Input must contain one row per variant_id")
    sequences = load_fasta(Path(args.fasta))
    manifest = load_manifest(Path(args.manifest))
    gene_accessions = build_gene_accessions(frame, sequences)

    initial = []
    transcript_targets = {}
    for row in frame.itertuples(index=False):
        direct = direct_candidate(row, sequences)
        dbnsfp = dbnsfp_candidate(row, sequences) if direct is None else None
        vep = vep_candidate(row, sequences, gene_accessions) if direct is None and dbnsfp is None else None
        initial.append((direct, dbnsfp, vep))
        if direct is None and dbnsfp is None and vep is None and truthy(row.vep_any_missense):
            transcript_targets[row.variant_id] = clean(row.primary_gene)
    transcript_rescues = load_transcript_vep_rescues(
        Path(args.vep_tab), transcript_targets, sequences, gene_accessions
    )
    rest_targets = {
        variant_id: gene for variant_id, gene in transcript_targets.items()
        if variant_id not in transcript_rescues
    }
    rest_rescues = load_rest_vep_rescues(
        Path(args.vep_rest_json), rest_targets, sequences, gene_accessions
    )

    proposed = []
    relevant_accessions: set[str] = set()
    for row, (direct, dbnsfp, vep) in zip(frame.itertuples(index=False), initial):
        transcript_vep = transcript_rescues.get(row.variant_id)
        rest_vep = rest_rescues.get(row.variant_id)
        candidate = direct or dbnsfp or vep or transcript_vep or rest_vep
        method = (
            "alphamissense_direct" if direct else
            "dbnsfp_aligned" if dbnsfp else
            "vep_gene_canonical" if vep else
            "vep_transcript_canonical_sequence" if transcript_vep else
            "vep_rest_transcript_canonical_sequence" if rest_vep else ""
        )
        if candidate:
            relevant_accessions.add(candidate[0])
        proposed.append((candidate, method))

    payloads: dict[str, tuple[str, dict[int, float], dict[int, str]]] = {}
    fragment_lengths: dict[str, int] = {}
    with tarfile.open(args.alphafold_tar, "r") as archive:
        for accession in sorted(relevant_accessions):
            for _, member in manifest.get(accession, []):
                handle = archive.extractfile(member)
                if handle is None:
                    continue
                text = gzip.decompress(handle.read()).decode("utf-8", "ignore")
                plddt, residues = parse_pdb(text)
                payloads[member] = (text, plddt, residues)
                fragment_lengths[member] = max(residues, default=0)

    audit_rows = []
    feature_rows = []
    requested_members: set[str] = set()
    for row, (candidate, method) in zip(frame.itertuples(index=False), proposed):
        audit = {
            "variant_id": row.variant_id,
            "primary_gene": row.primary_gene,
            "split_source_heldout": row.split_source_heldout,
            "training_slice_primary_binary": row.training_slice_primary_binary,
            "primary_binary_label": row.primary_binary_label,
            "audit_mapping_method": method,
            "audit_uniprot_accession": "",
            "audit_protein_variant": "",
            "audit_mapping_status": "not_mapped",
            "audit_missing_reason": "No exact, reference-validated amino-acid substitution mapping",
            "audit_exclusion_category": "not_mapped",
            "audit_global_position": pd.NA,
            "audit_fragment_index": pd.NA,
            "audit_fragment_member": "",
            "audit_local_position": pd.NA,
            "audit_fragment_edge_distance": pd.NA,
            "audit_residue_plddt": pd.NA,
            "audit_plddt_tier": "not_mapped",
            "foldx_eligible_any_plddt": False,
            "foldx_eligible_plddt_ge50": False,
            "foldx_eligible_plddt_ge70": False,
            "foldx_eligible_plddt_ge90": False,
        }
        if candidate is not None:
            accession, ref, position, alt = candidate
            audit.update({
                "audit_uniprot_accession": accession,
                "audit_protein_variant": f"{ref}{position}{alt}",
                "audit_global_position": position,
            })
            selected = select_fragment(manifest.get(accession, []), position, fragment_lengths)
            if selected is None:
                audit["audit_missing_reason"] = "No AlphaFold v6 fragment contains the protein position"
                audit["audit_exclusion_category"] = "no_alphafold_fragment_at_position"
            else:
                fragment, member, local_position, edge_distance = selected
                plddt, residues = payloads[member][1], payloads[member][2]
                if residues.get(local_position) != ref:
                    audit["audit_missing_reason"] = "Reference amino acid mismatch against selected AlphaFold fragment"
                    audit["audit_exclusion_category"] = "alphafold_reference_mismatch"
                else:
                    confidence = float(plddt[local_position])
                    audit.update({
                        "audit_mapping_status": "ok",
                        "audit_missing_reason": "",
                        "audit_exclusion_category": "eligible",
                        "audit_fragment_index": fragment,
                        "audit_fragment_member": member,
                        "audit_local_position": local_position,
                        "audit_fragment_edge_distance": edge_distance,
                        "audit_residue_plddt": confidence,
                        "audit_plddt_tier": plddt_tier(confidence),
                        "foldx_eligible_any_plddt": True,
                        "foldx_eligible_plddt_ge50": confidence >= 50,
                        "foldx_eligible_plddt_ge70": confidence >= 70,
                        "foldx_eligible_plddt_ge90": confidence >= 90,
                    })
                    requested_members.add(member)
        elif truthy(row.vep_any_missense):
            audit["audit_exclusion_category"] = "missense_without_unambiguous_reviewed_protein_mapping"
        else:
            audit["audit_exclusion_category"] = "not_an_exact_missense_substitution"
        audit_rows.append(audit)
        feature_rows.append({
            "variant_id": row.variant_id,
            "structure_available": audit["audit_mapping_status"] == "ok",
            "structure_plddt": audit["audit_residue_plddt"],
            "structure_plddt_bin": audit["audit_plddt_tier"],
            "structure_mapping_method": audit["audit_mapping_method"] or "not_mapped",
        })

    pdb_dir = Path(args.work_dir) / "pdb"
    pdb_dir.mkdir(parents=True, exist_ok=True)
    for member in sorted(requested_members):
        path = pdb_dir / member.removesuffix(".gz")
        if not path.exists() or not path.stat().st_size:
            path.write_text(payloads[member][0], encoding="utf-8")

    audit_frame = pd.DataFrame(audit_rows)
    features = pd.DataFrame(feature_rows)
    Path(args.audit_output).parent.mkdir(parents=True, exist_ok=True)
    audit_frame.to_csv(args.audit_output, sep="\t", index=False)
    features.to_csv(args.features_output, sep="\t", index=False)

    eligible = audit_frame[audit_frame["training_slice_primary_binary"].map(truthy)].copy()
    coverage = {}
    for split, group in eligible.groupby("split_source_heldout", dropna=False):
        coverage[str(split)] = {
            "cohort_rows": len(group),
            "mapped_any_plddt": int(group["foldx_eligible_any_plddt"].sum()),
            "plddt_ge50": int(group["foldx_eligible_plddt_ge50"].sum()),
            "plddt_ge70": int(group["foldx_eligible_plddt_ge70"].sum()),
            "plddt_ge90": int(group["foldx_eligible_plddt_ge90"].sum()),
            "mapping_methods": dict(Counter(group.loc[group["audit_mapping_status"].eq("ok"), "audit_mapping_method"])),
        }
    summary = {
        "input_table": str(Path(args.table)),
        "rows": len(frame),
        "primary_binary_rows": len(eligible),
        "selection_used_outcome_labels": False,
        "mapping_precedence": [
            "alphamissense_direct", "dbnsfp_aligned", "vep_gene_canonical",
            "vep_transcript_canonical_sequence", "vep_rest_transcript_canonical_sequence",
        ],
        "transcript_vep_rescues": len(transcript_rescues),
        "rest_vep_rescues": len(rest_rescues),
        "mapping_requirements": [
            "single amino-acid substitution", "reviewed UniProt reference agreement",
            "AlphaFold v6 fragment coverage", "AlphaFold fragment reference agreement",
        ],
        "confidence_analysis_plan": {
            "primary": "pLDDT>=70", "sensitivity": "pLDDT>=50", "exploratory_ceiling": "all mapped residues",
        },
        "coverage_by_split": coverage,
        "requested_fragments": len(requested_members),
        "audit_output": str(Path(args.audit_output)),
        "features_output": str(Path(args.features_output)),
        "pdb_dir": str(pdb_dir),
    }
    Path(args.summary).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
