#!/usr/bin/env python3
"""Run a small FoldXPro DDG pilot on clean mapped missense variants."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import shutil
import subprocess
import tarfile
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[4]
FOLDXPRO = ROOT / ".tools/foldx/bin/foldxpro"
ALPHAFOLD_TAR = ROOT / "datasets/feature_sources/protein_structure/raw/alphafold/UP000005640_9606_HUMAN_v6.tar"
ALPHAFOLD_MANIFEST = ROOT / "datasets/feature_sources/protein_structure/raw/alphafold/UP000005640_9606_HUMAN_v6.manifest.txt"
PROTEIN_FEATURES = ROOT / "datasets/feature_sources/protein_structure/interim/protein_features.tsv"
OUT_ROOT = ROOT / "datasets/feature_sources/protein_structure/interim/foldxpro_ddg_pilot"
OUT_TSV = ROOT / "datasets/feature_sources/protein_structure/interim/foldxpro_ddg_pilot.tsv"
OUT_JSON = ROOT / "datasets/feature_sources/protein_structure/interim/foldxpro_ddg_pilot.summary.json"


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
    "foldx_mutation",
    "ddg_method",
    "ddg_status",
    "ddg_missing_reason",
    "ddg_kcal_mol",
    "ddg_abs",
    "foldxpro_returncode",
    "foldxpro_mutant_pdb",
    "foldxpro_dif_fxout",
    "foldxpro_stdout",
    "foldxpro_stderr",
    "work_dir",
]

AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "SEC": "U",
    "PYL": "O",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)[:120].strip("_")


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


def extract_pdb(member_name: str, work_dir: Path) -> Path:
    input_dir = work_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = input_dir / member_name.replace(".gz", "")
    if pdb_path.exists() and pdb_path.stat().st_size > 0:
        return pdb_path
    with tarfile.open(ALPHAFOLD_TAR) as archive:
        member = archive.extractfile(member_name)
        if member is None:
            raise FileNotFoundError(member_name)
        pdb_path.write_bytes(gzip.decompress(member.read()))
    return pdb_path


def pdb_residue_aa(pdb_path: Path, position: int, chain: str = "A") -> str:
    with pdb_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            if line[21].strip() != chain:
                continue
            try:
                residue_position = int(line[22:26])
            except ValueError:
                continue
            if residue_position == position:
                return AA3_TO_1.get(line[17:20].strip(), "")
    return ""


def parse_dif_fxout(work_dir: Path, pdb_name: str) -> dict[str, str]:
    dif_path = work_dir / "output" / f"Dif_{Path(pdb_name).stem}.fxout"
    if not dif_path.exists():
        raise FileNotFoundError(dif_path)
    rows = [
        line.rstrip("\n")
        for line in dif_path.read_text(errors="replace").splitlines()
        if line.strip() and not line.startswith("FoldX") and not line.startswith("by ")
    ]
    header_index = next(i for i, line in enumerate(rows) if line.startswith("Pdb\t"))
    header = rows[header_index].split("\t")
    values = rows[header_index + 1].split("\t")
    parsed = dict(zip(header, values))
    parsed["dif_fxout"] = str(dif_path.relative_to(ROOT))
    return parsed


def eligible(row: dict[str, str], min_plddt: float) -> bool:
    if row.get("protein_feature_status") != "ok":
        return False
    if row.get("protein_variant_type") != "substitution":
        return False
    if row.get("alphafold_plddt_status") != "ok":
        return False
    ref = row.get("protein_ref_aa", "")
    alt = row.get("protein_alt_aa", "")
    if len(ref) != 1 or len(alt) != 1 or "*" in {ref, alt}:
        return False
    try:
        if float(row.get("alphafold_residue_plddt", "0")) < min_plddt:
            return False
        int(float(row.get("protein_position", "")))
    except ValueError:
        return False
    return bool(row.get("uniprot_accession"))


def select_pilot_rows(rows: list[dict[str, str]], limit: int, min_plddt: float) -> list[dict[str, str]]:
    candidates = [row for row in rows if eligible(row, min_plddt)]
    preferred_sources = {"cardioboost": 0, "hiro": 1, "emerge": 2}
    class_order = {"Pathogenic": 0, "Benign": 1, "VUS": 2}
    candidates.sort(
        key=lambda row: (
            class_order.get(row.get("target_3class", ""), 9),
            preferred_sources.get(row.get("source_dataset", ""), 9),
            row.get("gene", ""),
            row.get("variant_uid", ""),
        )
    )

    selected: list[dict[str, str]] = []
    seen = set()
    quotas = {"Pathogenic": max(1, limit // 3), "Benign": max(1, limit // 3), "VUS": max(1, limit - 2 * (limit // 3))}
    for label in ["Pathogenic", "Benign", "VUS"]:
        for row in candidates:
            key = (row["uniprot_accession"], row["protein_position"], row["protein_ref_aa"], row["protein_alt_aa"])
            if row.get("target_3class") != label or key in seen:
                continue
            selected.append(row)
            seen.add(key)
            if sum(1 for item in selected if item.get("target_3class") == label) >= quotas[label]:
                break
    for row in candidates:
        if len(selected) >= limit:
            break
        key = (row["uniprot_accession"], row["protein_position"], row["protein_ref_aa"], row["protein_alt_aa"])
        if key not in seen:
            selected.append(row)
            seen.add(key)
    return selected[:limit]


def run_one(row: dict[str, str], member_name: str, args: argparse.Namespace) -> dict[str, Any]:
    position = int(float(row["protein_position"]))
    mutation = f"{row['protein_ref_aa']}A{position}{row['protein_alt_aa']};"
    work_dir = Path(args.work_dir) / f"{safe_name(row['gene'])}_{safe_name(row['variant_uid'])}"
    output_dir = work_dir / "output"
    work_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(exist_ok=True)

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
    out["foldx_mutation"] = mutation
    out["ddg_method"] = "FoldXPro_2025_BuildModel_AlphaFold_monomer"
    out["work_dir"] = str(work_dir.relative_to(ROOT))

    try:
        pdb_path = extract_pdb(member_name, work_dir)
        pdb_ref = pdb_residue_aa(pdb_path, position)
        if not pdb_ref:
            out["ddg_status"] = "residue_not_found"
            out["ddg_missing_reason"] = f"Residue A{position} not found in AlphaFold PDB."
            return out
        if pdb_ref != row["protein_ref_aa"]:
            out["ddg_status"] = "ref_mismatch"
            out["ddg_missing_reason"] = (
                f"Variant ref AA {row['protein_ref_aa']} does not match AlphaFold PDB "
                f"chain A residue {position}={pdb_ref}; likely transcript/isoform mismatch."
            )
            return out
        (work_dir / "individual_list.txt").write_text(mutation + "\n", encoding="utf-8")
        (work_dir / "foldxpro_buildmodel.cfg").write_text(
            "\n".join(
                [
                    "command=BuildModel",
                    f"pdb={pdb_path.name}",
                    "pdb-dir=input",
                    "output-dir=output",
                    "mutant-file=individual_list.txt",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        result = subprocess.run(
            [str(args.foldxpro), "--config", "foldxpro_buildmodel.cfg"],
            cwd=work_dir,
            text=True,
            capture_output=True,
            check=False,
            timeout=args.timeout_seconds,
        )
        stdout = work_dir / "foldxpro.stdout.log"
        stderr = work_dir / "foldxpro.stderr.log"
        stdout.write_text(result.stdout, encoding="utf-8")
        stderr.write_text(result.stderr, encoding="utf-8")
        out["foldxpro_returncode"] = str(result.returncode)
        out["foldxpro_stdout"] = str(stdout.relative_to(ROOT))
        out["foldxpro_stderr"] = str(stderr.relative_to(ROOT))
        if result.returncode != 0:
            out["ddg_status"] = "failed"
            out["ddg_missing_reason"] = (result.stderr or result.stdout or "FoldXPro failed").strip().replace("\n", " ")[:500]
            return out
        parsed = parse_dif_fxout(work_dir, pdb_path.name)
        ddg = parsed.get("total energy", "")
        out["ddg_status"] = "ok"
        out["ddg_kcal_mol"] = ddg
        out["ddg_abs"] = str(abs(float(ddg))) if ddg else ""
        out["foldxpro_mutant_pdb"] = parsed.get("Pdb", "")
        out["foldxpro_dif_fxout"] = parsed.get("dif_fxout", "")
        return out
    except subprocess.TimeoutExpired:
        out["ddg_status"] = "timeout"
        out["ddg_missing_reason"] = f"FoldXPro exceeded timeout_seconds={args.timeout_seconds}"
        return out
    except Exception as exc:
        out["ddg_status"] = "failed"
        out["ddg_missing_reason"] = str(exc)[:500]
        return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protein-features", default=str(PROTEIN_FEATURES))
    parser.add_argument("--foldxpro", default=str(FOLDXPRO))
    parser.add_argument("--work-dir", default=str(OUT_ROOT))
    parser.add_argument("--output", default=str(OUT_TSV))
    parser.add_argument("--summary", default=str(OUT_JSON))
    parser.add_argument("--limit", type=int, default=9)
    parser.add_argument("--min-plddt", type=float, default=70.0)
    parser.add_argument("--timeout-seconds", type=int, default=600)
    args = parser.parse_args()

    rows = read_tsv(Path(args.protein_features))
    selected = select_pilot_rows(rows, args.limit, args.min_plddt)
    manifest = load_manifest(ALPHAFOLD_MANIFEST)

    output_rows = []
    for row in selected:
        names = manifest.get(row["uniprot_accession"], [])
        member_name = next((name for name in names if "-F1-" in name and name.endswith(".pdb.gz")), "")
        if not member_name:
            out = {field: "" for field in OUTPUT_FIELDS}
            out.update({field: row.get(field, "") for field in out})
            out["ddg_status"] = "alphafold_missing"
            out["ddg_missing_reason"] = "No AlphaFold F1 PDB member found."
        else:
            out = run_one(row, member_name, args)
        output_rows.append(out)
        write_tsv(Path(args.output), output_rows, OUTPUT_FIELDS)

    summary = {
        "input": args.protein_features,
        "output": args.output,
        "work_dir": args.work_dir,
        "limit": args.limit,
        "min_plddt": args.min_plddt,
        "selected_rows": len(selected),
        "status_counts": dict(Counter(row.get("ddg_status", "") for row in output_rows)),
        "label_counts": dict(Counter(row.get("target_3class", "") for row in output_rows)),
        "gene_counts": dict(Counter(row.get("gene", "") for row in output_rows)),
    }
    Path(args.summary).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
