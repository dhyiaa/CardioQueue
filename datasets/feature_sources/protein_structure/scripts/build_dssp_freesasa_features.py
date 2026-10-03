#!/usr/bin/env python3
"""Batch DSSP secondary-structure and FreeSASA exposure features.

The script uses the local AlphaFold human bulk tar and the existing
protein_features.tsv mapping table. It runs tools once per UniProt accession and
joins residue-level features back to variant rows.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import subprocess
import tarfile
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[4]
DEFAULT_ENV = ROOT / ".tools/envs/protein-structure"
DEFAULT_MKDSSP = DEFAULT_ENV / "bin/mkdssp"
DEFAULT_PYTHON = DEFAULT_ENV / "bin/python"
DEFAULT_LIBCIFPP = DEFAULT_ENV / "share/libcifpp"

MAX_ASA_TIEN = {
    "A": 129.0,
    "R": 274.0,
    "N": 195.0,
    "D": 193.0,
    "C": 167.0,
    "Q": 223.0,
    "E": 225.0,
    "G": 104.0,
    "H": 224.0,
    "I": 197.0,
    "L": 201.0,
    "K": 236.0,
    "M": 224.0,
    "F": 240.0,
    "P": 159.0,
    "S": 155.0,
    "T": 172.0,
    "W": 285.0,
    "Y": 263.0,
    "V": 174.0,
}

SS_CLASS = {
    "H": "helix",
    "G": "helix",
    "I": "helix",
    "E": "strand",
    "B": "strand",
    "T": "turn",
    "S": "bend",
    "": "coil",
    " ": "coil",
}

OUTPUT_FIELDS = [
    "source_dataset",
    "variant_uid",
    "variant_key",
    "target_3class",
    "gene",
    "uniprot_accession",
    "protein_position",
    "protein_ref_aa",
    "protein_alt_aa",
    "alphafold_residue_plddt",
    "dssp_status",
    "dssp_missing_reason",
    "dssp_secondary_structure",
    "dssp_secondary_structure_class",
    "dssp_acc",
    "dssp_phi",
    "dssp_psi",
    "freesasa_status",
    "freesasa_missing_reason",
    "freesasa_total",
    "freesasa_polar",
    "freesasa_apolar",
    "freesasa_relative",
    "freesasa_exposure_bin",
    "alphafold_pdb_path",
    "dssp_output_path",
]


def clean(value: Any) -> str:
    text = "" if value is None else str(value).strip()
    return "" if text in {"", "Missing", "nan", "NaN", "None"} else text


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_manifest(path: Path) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    with path.open() as handle:
        for line in handle:
            name = line.strip()
            if not name.endswith(".pdb.gz") or not name.startswith("AF-"):
                continue
            accession = name.split("-")[1]
            out.setdefault(accession, []).append(name)
    for names in out.values():
        names.sort()
    return out


def extract_pdb(tar_path: Path, member_name: str, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = out_dir / member_name.replace(".gz", "")
    if pdb_path.exists() and pdb_path.stat().st_size > 0:
        return pdb_path
    with tarfile.open(tar_path, "r") as archive:
        member = archive.extractfile(member_name)
        if member is None:
            raise FileNotFoundError(member_name)
        pdb_path.write_bytes(gzip.decompress(member.read()))
    return pdb_path


def run_dssp(mkdssp: Path, pdb_path: Path, dssp_path: Path, libcifpp_dir: Path) -> tuple[str, str]:
    dssp_path.parent.mkdir(parents=True, exist_ok=True)
    if dssp_path.exists() and dssp_path.stat().st_size > 0:
        return "ok", ""
    env = os.environ.copy()
    if libcifpp_dir.exists():
        env["LIBCIFPP_DATA_DIR"] = str(libcifpp_dir)
    result = subprocess.run(
        [str(mkdssp), str(pdb_path), str(dssp_path)],
        text=True,
        capture_output=True,
        check=False,
        env=env,
    )
    if result.returncode != 0 or not dssp_path.exists():
        reason = (result.stderr or result.stdout or "mkdssp failed").strip().replace("\n", " ")[:500]
        return "failed", reason
    return "ok", ""


def parse_dssp(path: Path) -> dict[int, dict[str, str]]:
    residues: dict[int, dict[str, str]] = {}
    if not path.exists():
        return residues
    started = False
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("  #  RESIDUE"):
            started = True
            continue
        if not started or len(line) < 115:
            continue
        try:
            pos = int(line[5:10])
        except ValueError:
            continue
        chain = line[11].strip()
        if chain and chain != "A":
            continue
        ss = line[16].strip()
        acc = line[34:38].strip()
        phi = line[103:109].strip()
        psi = line[109:115].strip()
        residues[pos] = {
            "dssp_secondary_structure": ss,
            "dssp_secondary_structure_class": SS_CLASS.get(ss, "other"),
            "dssp_acc": acc,
            "dssp_phi": phi,
            "dssp_psi": psi,
        }
    return residues


def run_freesasa(python_bin: Path, pdb_path: Path) -> tuple[str, str, dict[int, dict[str, str]]]:
    code = r"""
import json
import sys
import freesasa

structure = freesasa.Structure(sys.argv[1])
result = freesasa.calc(structure)
areas = result.residueAreas()
payload = {}
for chain, residues in areas.items():
    if chain != "A":
        continue
    for residue_number, area in residues.items():
        payload[str(int(residue_number))] = {
            "total": area.total,
            "polar": area.polar,
            "apolar": area.apolar,
        }
print(json.dumps(payload))
"""
    result = subprocess.run(
        [str(python_bin), "-c", code, str(pdb_path)],
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        reason = (result.stderr or result.stdout or "FreeSASA failed").strip().replace("\n", " ")[:500]
        return "failed", reason, {}
    raw = json.loads(result.stdout or "{}")
    rows = {
        int(pos): {
            "freesasa_total": f"{values['total']:.4f}",
            "freesasa_polar": f"{values['polar']:.4f}",
            "freesasa_apolar": f"{values['apolar']:.4f}",
        }
        for pos, values in raw.items()
    }
    return "ok", "", rows


def exposure_bin(relative: str) -> str:
    if not relative:
        return ""
    value = float(relative)
    if value < 0.05:
        return "buried"
    if value < 0.25:
        return "partly_buried"
    if value < 0.50:
        return "intermediate"
    return "exposed"


def eligible_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [
        row
        for row in rows
        if row.get("protein_feature_status") == "ok"
        and row.get("alphafold_plddt_status") == "ok"
        and clean(row.get("uniprot_accession"))
        and clean(row.get("protein_position"))
    ]


def build(args: argparse.Namespace) -> dict[str, Any]:
    protein_rows = read_tsv(Path(args.protein_features))
    rows = eligible_rows(protein_rows)
    manifest = load_manifest(Path(args.alphafold_manifest))
    tar_path = Path(args.alphafold_tar)
    work_dir = Path(args.work_dir)
    pdb_dir = work_dir / "pdb"
    dssp_dir = work_dir / "dssp"

    accessions = sorted({row["uniprot_accession"] for row in rows})
    protein_cache: dict[str, dict[str, Any]] = {}
    protein_status_counts: Counter[str] = Counter()

    for accession in accessions:
        names = manifest.get(accession, [])
        if not names:
            protein_cache[accession] = {"status": "alphafold_missing", "reason": "accession absent from manifest"}
            protein_status_counts["alphafold_missing"] += 1
            continue
        member_name = next((name for name in names if "-F1-" in name), names[0])
        try:
            pdb_path = extract_pdb(tar_path, member_name, pdb_dir)
        except Exception as exc:
            protein_cache[accession] = {"status": "extract_failed", "reason": str(exc)}
            protein_status_counts["extract_failed"] += 1
            continue

        dssp_path = dssp_dir / f"{pdb_path.stem}.dssp"
        dssp_status, dssp_reason = run_dssp(Path(args.mkdssp), pdb_path, dssp_path, Path(args.libcifpp_dir))
        dssp_rows = parse_dssp(dssp_path) if dssp_status == "ok" else {}

        freesasa_status, freesasa_reason, sasa_rows = run_freesasa(Path(args.python_bin), pdb_path)

        protein_cache[accession] = {
            "status": "ok",
            "pdb_path": str(pdb_path),
            "dssp_path": str(dssp_path),
            "dssp_status": dssp_status,
            "dssp_reason": dssp_reason,
            "dssp_rows": dssp_rows,
            "freesasa_status": freesasa_status,
            "freesasa_reason": freesasa_reason,
            "sasa_rows": sasa_rows,
        }
        protein_status_counts["ok"] += 1

    output_rows: list[dict[str, Any]] = []
    dssp_status_counts: Counter[str] = Counter()
    freesasa_status_counts: Counter[str] = Counter()

    for row in rows:
        out = {field: "" for field in OUTPUT_FIELDS}
        for field in [
            "source_dataset",
            "variant_uid",
            "variant_key",
            "target_3class",
            "gene",
            "uniprot_accession",
            "protein_position",
            "protein_ref_aa",
            "protein_alt_aa",
            "alphafold_residue_plddt",
        ]:
            out[field] = row.get(field, "")

        accession = row["uniprot_accession"]
        position = int(float(row["protein_position"]))
        cached = protein_cache.get(accession, {})
        out["alphafold_pdb_path"] = cached.get("pdb_path", "")
        out["dssp_output_path"] = cached.get("dssp_path", "")

        if cached.get("status") != "ok":
            out["dssp_status"] = cached.get("status", "not_run")
            out["dssp_missing_reason"] = cached.get("reason", "")
            out["freesasa_status"] = cached.get("status", "not_run")
            out["freesasa_missing_reason"] = cached.get("reason", "")
        else:
            dssp_rows = cached["dssp_rows"]
            dssp_hit = dssp_rows.get(position)
            if dssp_hit:
                out["dssp_status"] = "ok"
                out.update(dssp_hit)
            else:
                out["dssp_status"] = cached.get("dssp_status", "not_found")
                out["dssp_missing_reason"] = cached.get("dssp_reason") or "residue absent from DSSP output"

            sasa_rows = cached["sasa_rows"]
            sasa_hit = sasa_rows.get(position)
            if sasa_hit:
                out["freesasa_status"] = "ok"
                out.update(sasa_hit)
                max_asa = MAX_ASA_TIEN.get(row.get("protein_ref_aa", ""))
                if max_asa:
                    relative = float(out["freesasa_total"]) / max_asa
                    out["freesasa_relative"] = f"{relative:.4f}"
                    out["freesasa_exposure_bin"] = exposure_bin(out["freesasa_relative"])
            else:
                out["freesasa_status"] = cached.get("freesasa_status", "not_found")
                out["freesasa_missing_reason"] = cached.get("freesasa_reason") or "residue absent from FreeSASA output"

        dssp_status_counts[out["dssp_status"]] += 1
        freesasa_status_counts[out["freesasa_status"]] += 1
        output_rows.append(out)

    write_tsv(Path(args.output), output_rows, OUTPUT_FIELDS)
    summary = {
        "input": args.protein_features,
        "output": args.output,
        "work_dir": args.work_dir,
        "eligible_variant_rows": len(rows),
        "unique_uniprot_accessions": len(accessions),
        "protein_status_counts": dict(protein_status_counts),
        "dssp_status_counts": dict(dssp_status_counts),
        "freesasa_status_counts": dict(freesasa_status_counts),
    }
    Path(args.summary).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protein-features", default="datasets/feature_sources/protein_structure/interim/protein_features.tsv")
    parser.add_argument("--alphafold-manifest", default="datasets/feature_sources/protein_structure/raw/alphafold/UP000005640_9606_HUMAN_v6.manifest.txt")
    parser.add_argument("--alphafold-tar", default="datasets/feature_sources/protein_structure/raw/alphafold/UP000005640_9606_HUMAN_v6.tar")
    parser.add_argument("--mkdssp", default=str(DEFAULT_MKDSSP))
    parser.add_argument("--python-bin", default=str(DEFAULT_PYTHON))
    parser.add_argument("--libcifpp-dir", default=str(DEFAULT_LIBCIFPP))
    parser.add_argument("--work-dir", default="datasets/feature_sources/protein_structure/interim/dssp_freesasa_batch")
    parser.add_argument("--output", default="datasets/feature_sources/protein_structure/interim/dssp_freesasa_features.tsv")
    parser.add_argument("--summary", default="datasets/feature_sources/protein_structure/interim/dssp_freesasa_features.summary.json")
    args = parser.parse_args()
    print(json.dumps(build(args), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
