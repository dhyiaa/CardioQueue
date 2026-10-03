#!/usr/bin/env python3
"""Run a one-variant FoldXPro DDG smoke test.

This intentionally uses a small, high-confidence ACTC1/CardioBoost variant so
we can validate command syntax and output parsing before scoring a cohort.
"""

from __future__ import annotations

import csv
import gzip
import json
import subprocess
import tarfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
FOLDXPRO = ROOT / ".tools/foldx/bin/foldxpro"
ALPHAFOLD_TAR = (
    ROOT
    / "datasets/feature_sources/protein_structure/raw/alphafold/UP000005640_9606_HUMAN_v6.tar"
)
OUT_ROOT = ROOT / "datasets/feature_sources/protein_structure/interim/foldxpro_smoke"
SMOKE_TSV = ROOT / "datasets/feature_sources/protein_structure/interim/foldxpro_smoke_test.tsv"
SMOKE_JSON = ROOT / "datasets/feature_sources/protein_structure/interim/foldxpro_smoke_test.summary.json"


SMOKE_VARIANT = {
    "source_dataset": "cardioboost",
    "variant_uid": "CardioBoost_reannotated_cardiomyopathy|ACTC1|c.401T>C|1",
    "gene": "ACTC1",
    "uniprot_accession": "P68032",
    "alphafold_member": "AF-P68032-F1-model_v6.pdb.gz",
    "pdb_name": "AF-P68032-F1-model_v6.pdb",
    "protein_change": "p.M134T",
    "foldx_mutation": "MA134T;",
    "protein_position": "134",
    "protein_ref_aa": "M",
    "protein_alt_aa": "T",
    "alphafold_residue_plddt": "98.12",
}


def extract_alphafold_pdb(work_dir: Path) -> Path:
    input_dir = work_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = input_dir / SMOKE_VARIANT["pdb_name"]
    with tarfile.open(ALPHAFOLD_TAR) as tar:
        member = tar.extractfile(SMOKE_VARIANT["alphafold_member"])
        if member is None:
            raise FileNotFoundError(SMOKE_VARIANT["alphafold_member"])
        with gzip.GzipFile(fileobj=member) as gz, pdb_path.open("wb") as out:
            out.write(gz.read())
    return pdb_path


def write_foldx_inputs(work_dir: Path) -> None:
    (work_dir / "output").mkdir(parents=True, exist_ok=True)
    (work_dir / "individual_list.txt").write_text(SMOKE_VARIANT["foldx_mutation"] + "\n")
    (work_dir / "foldxpro_buildmodel.cfg").write_text(
        "\n".join(
            [
                "command=BuildModel",
                f"pdb={SMOKE_VARIANT['pdb_name']}",
                "pdb-dir=input",
                "output-dir=output",
                "mutant-file=individual_list.txt",
                "",
            ]
        )
    )


def run_foldxpro(work_dir: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(FOLDXPRO), "--config", "foldxpro_buildmodel.cfg"],
        cwd=work_dir,
        text=True,
        capture_output=True,
        check=False,
    )


def parse_dif_fxout(work_dir: Path) -> dict[str, str]:
    dif_path = work_dir / "output" / f"Dif_{Path(SMOKE_VARIANT['pdb_name']).stem}.fxout"
    if not dif_path.exists():
        raise FileNotFoundError(dif_path)
    rows = [
        line.rstrip("\n")
        for line in dif_path.read_text().splitlines()
        if line.strip() and not line.startswith("FoldX") and not line.startswith("by ")
    ]
    header_index = next(i for i, line in enumerate(rows) if line.startswith("Pdb\t"))
    header = rows[header_index].split("\t")
    values = rows[header_index + 1].split("\t")
    parsed = dict(zip(header, values))
    parsed["dif_fxout"] = str(dif_path.relative_to(ROOT))
    return parsed


def main() -> None:
    work_dir = OUT_ROOT / "ACTC1_M134T"
    work_dir.mkdir(parents=True, exist_ok=True)
    extract_alphafold_pdb(work_dir)
    write_foldx_inputs(work_dir)
    result = run_foldxpro(work_dir)

    row = dict(SMOKE_VARIANT)
    row["ddg_method"] = "FoldXPro_2025_BuildModel"
    row["foldxpro_returncode"] = str(result.returncode)
    row["foldxpro_stdout"] = str((work_dir / "foldxpro.stdout.log").relative_to(ROOT))
    row["foldxpro_stderr"] = str((work_dir / "foldxpro.stderr.log").relative_to(ROOT))
    (work_dir / "foldxpro.stdout.log").write_text(result.stdout)
    (work_dir / "foldxpro.stderr.log").write_text(result.stderr)

    if result.returncode == 0:
        parsed = parse_dif_fxout(work_dir)
        row["ddg_status"] = "ok"
        row["ddg_missing_reason"] = ""
        row["ddg_kcal_mol"] = parsed.get("total energy", "")
        row["ddg_abs"] = str(abs(float(row["ddg_kcal_mol"]))) if row["ddg_kcal_mol"] else ""
        row["foldxpro_mutant_pdb"] = parsed.get("Pdb", "")
        row["foldxpro_dif_fxout"] = parsed.get("dif_fxout", "")
    else:
        row["ddg_status"] = "failed"
        row["ddg_missing_reason"] = "FoldXPro returned non-zero exit status."
        row["ddg_kcal_mol"] = ""
        row["ddg_abs"] = ""
        row["foldxpro_mutant_pdb"] = ""
        row["foldxpro_dif_fxout"] = ""

    fields = [
        "source_dataset",
        "variant_uid",
        "gene",
        "uniprot_accession",
        "protein_change",
        "foldx_mutation",
        "protein_position",
        "protein_ref_aa",
        "protein_alt_aa",
        "alphafold_residue_plddt",
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
    ]
    with SMOKE_TSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow({field: row.get(field, "") for field in fields})

    summary = {
        "status": row["ddg_status"],
        "output": str(SMOKE_TSV.relative_to(ROOT)),
        "work_dir": str(work_dir.relative_to(ROOT)),
        "ddg_kcal_mol": row["ddg_kcal_mol"],
        "foldxpro_returncode": result.returncode,
    }
    SMOKE_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
