#!/usr/bin/env python3
"""Build leakage-free AlphaFold, DSSP, and FreeSASA features for registry variants."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import subprocess
import tarfile
from collections import Counter
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
DEFAULT_TABLE = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
DEFAULT_FASTA = ROOT / "datasets/feature_sources/protein_structure/raw/uniprot_human_reviewed_2026_02.fasta"
DEFAULT_MANIFEST = ROOT / "datasets/feature_sources/protein_structure/raw/alphafold/UP000005640_9606_HUMAN_v6.manifest.txt"
DEFAULT_TAR = ROOT / "datasets/feature_sources/protein_structure/raw/alphafold/UP000005640_9606_HUMAN_v6.tar"
DEFAULT_ENV = ROOT / ".tools/envs/protein-structure"
VARIANT_RE = re.compile(r"^([A-Z*])(\d+)([A-Z*])$")
MEMBER_RE = re.compile(r"^AF-(?P<accession>.+)-F(?P<fragment>\d+)-model_v(?P<version>\d+)\.pdb\.gz$")
FRAGMENT_OFFSET = 200

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
}

MAX_ASA_TIEN = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "Q": 223.0, "E": 225.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}

SS_CLASS = {
    "H": "helix", "G": "helix", "I": "helix", "E": "strand", "B": "strand",
    "T": "turn", "S": "bend", "": "coil", " ": "coil",
}

FEATURE_COLUMNS = [
    "structure_available",
    "structure_plddt",
    "structure_plddt_bin",
    "structure_dssp_code",
    "structure_dssp_class",
    "structure_dssp_acc",
    "structure_dssp_phi",
    "structure_dssp_psi",
    "structure_freesasa_total",
    "structure_freesasa_polar",
    "structure_freesasa_apolar",
    "structure_freesasa_relative",
    "structure_freesasa_exposure_bin",
]


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
    members: dict[str, list[tuple[int, str]]] = {}
    for name in path.read_text().splitlines():
        match = MEMBER_RE.match(name.strip())
        if not match:
            continue
        members.setdefault(match.group("accession"), []).append((int(match.group("fragment")), name.strip()))
    for values in members.values():
        values.sort()
    return members


def parse_variant(value: object) -> tuple[str, int, str] | None:
    match = VARIANT_RE.match(str(value))
    if not match:
        return None
    return match.group(1), int(match.group(2)), match.group(3)


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


def select_fragment(
    candidates: list[tuple[int, str]], global_position: int, fragment_lengths: dict[str, int]
) -> tuple[int, str, int, int] | None:
    valid = []
    for fragment, member in candidates:
        local_position = global_position - FRAGMENT_OFFSET * (fragment - 1)
        length = fragment_lengths[member]
        if 1 <= local_position <= length:
            edge_distance = min(local_position - 1, length - local_position)
            valid.append((edge_distance, -fragment, fragment, member, local_position))
    if not valid:
        return None
    edge_distance, _, fragment, member, local_position = max(valid)
    return fragment, member, local_position, edge_distance


def run_dssp(mkdssp: Path, pdb_path: Path, output: Path, libcifpp: Path) -> tuple[str, str]:
    if output.exists() and output.stat().st_size:
        return "ok", ""
    env = os.environ.copy()
    if libcifpp.exists():
        env["LIBCIFPP_DATA_DIR"] = str(libcifpp)
    result = subprocess.run([str(mkdssp), str(pdb_path), str(output)], capture_output=True, text=True, env=env)
    if result.returncode or not output.exists():
        return "failed", (result.stderr or result.stdout or "mkdssp failed").strip()[:500]
    return "ok", ""


def parse_dssp(path: Path) -> dict[int, dict[str, object]]:
    rows: dict[int, dict[str, object]] = {}
    started = False
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("  #  RESIDUE"):
            started = True
            continue
        if not started or len(line) < 115:
            continue
        try:
            position = int(line[5:10])
        except ValueError:
            continue
        if line[11].strip() not in {"", "A"}:
            continue
        code = line[16].strip()
        rows[position] = {
            "structure_dssp_code": code or "coil",
            "structure_dssp_class": SS_CLASS.get(code, "other"),
            "structure_dssp_acc": pd.to_numeric(line[34:38].strip(), errors="coerce"),
            "structure_dssp_phi": pd.to_numeric(line[103:109].strip(), errors="coerce"),
            "structure_dssp_psi": pd.to_numeric(line[109:115].strip(), errors="coerce"),
        }
    return rows


def run_freesasa(python_bin: Path, pdb_path: Path) -> tuple[str, str, dict[int, dict[str, float]]]:
    code = """
import json, sys, freesasa
s = freesasa.Structure(sys.argv[1]); r = freesasa.calc(s); out = {}
for chain, residues in r.residueAreas().items():
    if chain != 'A': continue
    for pos, area in residues.items():
        out[str(int(pos))] = {'total': area.total, 'polar': area.polar, 'apolar': area.apolar}
print(json.dumps(out))
"""
    result = subprocess.run([str(python_bin), "-c", code, str(pdb_path)], capture_output=True, text=True)
    if result.returncode:
        return "failed", (result.stderr or result.stdout or "FreeSASA failed").strip()[:500], {}
    raw = json.loads(result.stdout or "{}")
    return "ok", "", {int(key): value for key, value in raw.items()}


def plddt_bin(value: float) -> str:
    if value < 50:
        return "very_low"
    if value < 70:
        return "low"
    if value < 90:
        return "confident"
    return "very_high"


def exposure_bin(value: float) -> str:
    if value < 0.05:
        return "buried"
    if value < 0.25:
        return "partly_buried"
    if value < 0.50:
        return "intermediate"
    return "exposed"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", default=str(DEFAULT_TABLE))
    parser.add_argument("--fasta", default=str(DEFAULT_FASTA))
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    parser.add_argument("--alphafold-tar", default=str(DEFAULT_TAR))
    parser.add_argument("--mkdssp", default=str(DEFAULT_ENV / "bin/mkdssp"))
    parser.add_argument("--python-bin", default=str(DEFAULT_ENV / "bin/python"))
    parser.add_argument("--libcifpp", default=str(DEFAULT_ENV / "share/libcifpp"))
    parser.add_argument("--work-dir", default=str(ROOT / "datasets/feature_sources/protein_structure/interim/registry_structure"))
    parser.add_argument("--output", default=str(ROOT / "datasets/feature_sources/protein_structure/interim/registry_structure_features.tsv"))
    parser.add_argument("--audit-output", default=str(ROOT / "datasets/feature_sources/protein_structure/interim/registry_structure_mapping_audit.tsv"))
    parser.add_argument("--summary", default=str(ROOT / "datasets/feature_sources/protein_structure/interim/registry_structure_features.summary.json"))
    args = parser.parse_args()

    work_dir = Path(args.work_dir)
    pdb_dir = work_dir / "pdb"
    dssp_dir = work_dir / "dssp"
    pdb_dir.mkdir(parents=True, exist_ok=True)
    dssp_dir.mkdir(parents=True, exist_ok=True)

    columns = [
        "variant_id", "split_source_heldout", "training_slice_primary_binary",
        "alphamissense_direct_uniprot_id", "alphamissense_direct_protein_variant",
    ]
    registry = pd.read_csv(args.table, sep="\t", usecols=columns, low_memory=False)
    sequences = load_fasta(Path(args.fasta))
    manifest = load_manifest(Path(args.manifest))

    audit_rows = []
    requested_members: set[str] = set()
    with tarfile.open(args.alphafold_tar, "r") as archive:
        member_payload: dict[str, tuple[str, dict[int, float], dict[int, str]]] = {}
        fragment_lengths: dict[str, int] = {}
        relevant_accessions = set(registry["alphamissense_direct_uniprot_id"].dropna().astype(str))
        for accession in sorted(relevant_accessions):
            for _, member in manifest.get(accession, []):
                handle = archive.extractfile(member)
                if handle is None:
                    continue
                text = gzip.decompress(handle.read()).decode("utf-8", "ignore")
                plddt, residues = parse_pdb(text)
                fragment_lengths[member] = max(residues, default=0)
                member_payload[member] = (text, plddt, residues)

        for row in registry.itertuples(index=False):
            accession = str(row.alphamissense_direct_uniprot_id) if pd.notna(row.alphamissense_direct_uniprot_id) else ""
            parsed = parse_variant(row.alphamissense_direct_protein_variant)
            audit = {
                "variant_id": row.variant_id,
                "split_source_heldout": row.split_source_heldout,
                "training_slice_primary_binary": row.training_slice_primary_binary,
                "audit_uniprot_accession": accession,
                "audit_protein_variant": row.alphamissense_direct_protein_variant,
                "audit_mapping_status": "not_mapped",
                "audit_missing_reason": "AlphaMissense protein mapping unavailable",
                "audit_global_position": pd.NA,
                "audit_fragment_index": pd.NA,
                "audit_fragment_member": "",
                "audit_local_position": pd.NA,
                "audit_fragment_edge_distance": pd.NA,
            }
            if not accession or parsed is None:
                audit_rows.append(audit)
                continue
            ref, global_position, _ = parsed
            audit["audit_global_position"] = global_position
            sequence = sequences.get(accession, "")
            if not sequence or global_position > len(sequence) or sequence[global_position - 1] != ref:
                audit["audit_missing_reason"] = "Reference amino acid mismatch against reviewed UniProt sequence"
                audit_rows.append(audit)
                continue
            selected = select_fragment(manifest.get(accession, []), global_position, fragment_lengths)
            if selected is None:
                audit["audit_missing_reason"] = "No AlphaFold fragment contains the global position"
                audit_rows.append(audit)
                continue
            fragment, member, local_position, edge_distance = selected
            _, _, pdb_residues = member_payload[member]
            if pdb_residues.get(local_position) != ref:
                audit["audit_missing_reason"] = "Reference amino acid mismatch against selected AlphaFold fragment"
                audit_rows.append(audit)
                continue
            audit.update(
                {
                    "audit_mapping_status": "ok",
                    "audit_missing_reason": "",
                    "audit_fragment_index": fragment,
                    "audit_fragment_member": member,
                    "audit_local_position": local_position,
                    "audit_fragment_edge_distance": edge_distance,
                }
            )
            requested_members.add(member)
            audit_rows.append(audit)

        for member in sorted(requested_members):
            text, _, _ = member_payload[member]
            pdb_path = pdb_dir / member.replace(".gz", "")
            if not pdb_path.exists() or not pdb_path.stat().st_size:
                pdb_path.write_text(text)

    fragment_features: dict[str, dict[str, object]] = {}
    tool_counts: Counter[str] = Counter()
    for member in sorted(requested_members):
        pdb_path = pdb_dir / member.replace(".gz", "")
        dssp_path = dssp_dir / f"{pdb_path.stem}.dssp"
        dssp_status, dssp_reason = run_dssp(Path(args.mkdssp), pdb_path, dssp_path, Path(args.libcifpp))
        dssp_rows = parse_dssp(dssp_path) if dssp_status == "ok" else {}
        sasa_status, sasa_reason, sasa_rows = run_freesasa(Path(args.python_bin), pdb_path)
        fragment_features[member] = {
            "dssp_status": dssp_status,
            "dssp_reason": dssp_reason,
            "dssp_rows": dssp_rows,
            "sasa_status": sasa_status,
            "sasa_reason": sasa_reason,
            "sasa_rows": sasa_rows,
        }
        tool_counts[f"dssp_{dssp_status}"] += 1
        tool_counts[f"freesasa_{sasa_status}"] += 1

    audit_frame = pd.DataFrame(audit_rows)
    feature_rows = []
    status_counts: Counter[str] = Counter()
    member_payload_cache: dict[str, tuple[dict[int, float], dict[int, str]]] = {}
    with tarfile.open(args.alphafold_tar, "r") as archive:
        for member in sorted(requested_members):
            handle = archive.extractfile(member)
            text = gzip.decompress(handle.read()).decode("utf-8", "ignore") if handle else ""
            member_payload_cache[member] = parse_pdb(text)

    for row in audit_frame.itertuples(index=False):
        out = {"variant_id": row.variant_id, **{column: pd.NA for column in FEATURE_COLUMNS}}
        out["structure_available"] = False
        if row.audit_mapping_status != "ok":
            status_counts[row.audit_mapping_status] += 1
            feature_rows.append(out)
            continue
        member = row.audit_fragment_member
        local_position = int(row.audit_local_position)
        plddt = member_payload_cache[member][0].get(local_position)
        calculated = fragment_features[member]
        dssp = calculated["dssp_rows"].get(local_position)
        sasa = calculated["sasa_rows"].get(local_position)
        if plddt is None or dssp is None or sasa is None:
            status_counts["calculation_missing"] += 1
            feature_rows.append(out)
            continue
        ref = str(row.audit_protein_variant)[0]
        relative = float(sasa["total"]) / MAX_ASA_TIEN[ref]
        out.update(
            {
                "structure_available": True,
                "structure_plddt": plddt,
                "structure_plddt_bin": plddt_bin(plddt),
                **dssp,
                "structure_freesasa_total": sasa["total"],
                "structure_freesasa_polar": sasa["polar"],
                "structure_freesasa_apolar": sasa["apolar"],
                "structure_freesasa_relative": relative,
                "structure_freesasa_exposure_bin": exposure_bin(relative),
            }
        )
        status_counts["ok"] += 1
        feature_rows.append(out)

    features = pd.DataFrame(feature_rows)
    if not features["variant_id"].is_unique or not audit_frame["variant_id"].is_unique:
        raise ValueError("Registry structure outputs must contain one row per variant_id")
    features.to_csv(args.output, sep="\t", index=False)
    audit_frame.to_csv(args.audit_output, sep="\t", index=False)

    joined = registry[["variant_id", "split_source_heldout", "training_slice_primary_binary"]].merge(
        features, on="variant_id", validate="one_to_one"
    )
    coverage = (
        joined.groupby("split_source_heldout")["structure_available"].agg(["sum", "count"]).astype(int).to_dict("index")
    )
    binary_coverage = (
        joined[joined["training_slice_primary_binary"].eq(True)]
        .groupby("split_source_heldout")["structure_available"]
        .agg(["sum", "count"])
        .astype(int)
        .to_dict("index")
    )
    summary = {
        "input_table": str(Path(args.table)),
        "rows": len(registry),
        "candidate_substitutions": int(registry["alphamissense_direct_protein_variant"].notna().sum()),
        "mapped_accessions": int(
            audit_frame.loc[audit_frame["audit_mapping_status"].eq("ok"), "audit_uniprot_accession"].nunique()
        ),
        "reference_check_failures": int(
            audit_frame["audit_missing_reason"].str.contains("Reference amino acid mismatch", na=False).sum()
        ),
        "feature_columns": FEATURE_COLUMNS,
        "fragment_offset": FRAGMENT_OFFSET,
        "fragment_selection": "maximal distance from fragment edge; lower fragment index breaks ties",
        "requested_fragments": len(requested_members),
        "mapping_status_counts": dict(Counter(audit_frame["audit_mapping_status"])),
        "calculation_status_counts": dict(status_counts),
        "tool_status_counts": dict(tool_counts),
        "coverage_by_split": coverage,
        "primary_binary_coverage_by_split": binary_coverage,
        "output": str(Path(args.output)),
        "audit_output": str(Path(args.audit_output)),
    }
    Path(args.summary).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
