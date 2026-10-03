#!/usr/bin/env python3
"""Run repaired, batched FoldX 5.1 DDG for mapped primary-cohort missense variants."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
DEFAULT_TABLE = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
DEFAULT_AUDIT = ROOT / "datasets/feature_sources/protein_structure/interim/registry_structure_mapping_audit.tsv"
DEFAULT_STRUCTURE = ROOT / "datasets/feature_sources/protein_structure/interim/registry_structure_features.tsv"
DEFAULT_PDB = ROOT / "datasets/feature_sources/protein_structure/interim/registry_structure/pdb"
DEFAULT_FOLDX = ROOT / ".tools/foldx/bin/foldx"
DEFAULT_ROTABASE = ROOT / ".tools/foldx/foldx5.1_MacIntel/rotabase.txt"
DEFAULT_WORK = ROOT / "datasets/feature_sources/protein_structure/interim/registry_foldx51_repaired"
DEFAULT_OUTPUT = ROOT / "datasets/feature_sources/protein_structure/interim/registry_foldx51_repaired_ddg.tsv"
DEFAULT_SUMMARY = ROOT / "datasets/feature_sources/protein_structure/interim/registry_foldx51_repaired_ddg.summary.json"
SAFE_RE = re.compile(r"[^A-Za-z0-9_.-]+")
TERMINAL_STATUSES = {"ok", "repair_failed", "build_failed", "parse_failed", "timeout"}

OUTPUT_COLUMNS = [
    "variant_id", "split_source_heldout", "primary_binary_label", "primary_gene",
    "uniprot_accession", "protein_variant", "protein_global_position",
    "alphafold_fragment_member", "alphafold_local_position", "alphafold_residue_plddt",
    "fragment_edge_distance", "foldx_calculation_key", "foldx_mutation", "ddg_method",
    "ddg_status", "ddg_missing_reason", "ddg_kcal_mol", "ddg_abs", "repair_pdb",
    "repair_log", "build_log", "dif_fxout", "chunk_id",
]


def safe(value: object) -> str:
    return SAFE_RE.sub("_", str(value)).strip("_")


def root_relative(path: Path) -> str:
    resolved = path.resolve()
    return str(resolved.relative_to(ROOT)) if resolved.is_relative_to(ROOT) else str(resolved)


def truthy(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def load_candidates(args: argparse.Namespace) -> pd.DataFrame:
    table = pd.read_csv(
        args.table,
        sep="\t",
        usecols=["variant_id", "split_source_heldout", "training_slice_primary_binary", "primary_binary_label", "primary_gene"],
        low_memory=False,
    )
    audit = pd.read_csv(
        args.audit,
        sep="\t",
        usecols=[
            "variant_id", "audit_uniprot_accession", "audit_protein_variant", "audit_mapping_status",
            "audit_fragment_member", "audit_local_position", "audit_fragment_edge_distance",
        ],
        low_memory=False,
    )
    structure = pd.read_csv(args.structure, sep="\t", usecols=["variant_id", "structure_plddt"], low_memory=False)
    frame = table.merge(audit, on="variant_id", validate="one_to_one").merge(
        structure, on="variant_id", validate="one_to_one"
    )
    frame = frame[
        truthy(frame["training_slice_primary_binary"])
        & frame["audit_mapping_status"].eq("ok")
        & pd.to_numeric(frame["structure_plddt"], errors="coerce").ge(args.min_plddt)
    ].copy()
    parsed = frame["audit_protein_variant"].astype(str).str.extract(r"^([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])$")
    frame["protein_ref_aa"] = parsed[0]
    frame["protein_global_position"] = pd.to_numeric(parsed[1], errors="coerce").astype("Int64")
    frame["protein_alt_aa"] = parsed[2]
    frame = frame.dropna(subset=["protein_ref_aa", "protein_global_position", "protein_alt_aa"])
    frame = frame[frame["protein_ref_aa"].ne(frame["protein_alt_aa"])].copy()
    frame["foldx_mutation"] = (
        frame["protein_ref_aa"] + "A" + frame["audit_local_position"].astype("Int64").astype(str) + frame["protein_alt_aa"] + ";"
    )
    frame["foldx_calculation_key"] = (
        frame["audit_fragment_member"].astype(str) + "|" + frame["foldx_mutation"]
    )
    return frame.sort_values(["audit_fragment_member", "audit_local_position", "protein_alt_aa", "variant_id"])


def run_command(command: list[str], log_path: Path, timeout: int) -> tuple[int, str]:
    try:
        result = subprocess.run(command, text=True, capture_output=True, timeout=timeout, check=False)
        text = result.stdout + ("\n" + result.stderr if result.stderr else "")
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(text, encoding="utf-8")
        return result.returncode, text
    except subprocess.TimeoutExpired as exc:
        text = f"Timeout after {timeout} seconds\n{exc.stdout or ''}\n{exc.stderr or ''}"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(text, encoding="utf-8")
        return 124, text


def repair_fragment(member: str, args: argparse.Namespace) -> tuple[str, str, str]:
    pdb_name = member.removesuffix(".gz")
    stem = Path(pdb_name).stem
    fragment_dir = Path(args.work_dir) / "fragments" / safe(stem)
    repaired_dir = fragment_dir / "repaired"
    repaired_dir.mkdir(parents=True, exist_ok=True)
    repaired_pdb = repaired_dir / f"{stem}_Repair.pdb"
    repair_log = fragment_dir / "repair.log"
    if repaired_pdb.exists() and repaired_pdb.stat().st_size:
        return member, "ok", str(repaired_pdb)
    command = [
        str(Path(args.foldx).resolve()), "--command=RepairPDB", f"--pdb={pdb_name}",
        f"--pdb-dir={Path(args.pdb_dir).resolve()}", f"--output-dir={repaired_dir.resolve()}",
        f"--rotabaseLocation={Path(args.rotabase).resolve()}", "--screen=false",
    ]
    code, _ = run_command(command, repair_log, args.repair_timeout_seconds)
    status = "ok" if code == 0 and repaired_pdb.exists() and repaired_pdb.stat().st_size else "repair_failed"
    return member, status, str(repaired_pdb)


def parse_dif(path: Path) -> list[float]:
    lines = [line for line in path.read_text(errors="replace").splitlines() if line.strip()]
    header_index = next(i for i, line in enumerate(lines) if line.startswith("Pdb\t"))
    header = lines[header_index].split("\t")
    energy_index = header.index("total energy")
    return [float(line.split("\t")[energy_index]) for line in lines[header_index + 1:]]


def build_chunk(chunk: pd.DataFrame, repaired_pdb: Path, chunk_id: str, args: argparse.Namespace) -> pd.DataFrame:
    chunk_dir = Path(args.work_dir) / "chunks" / chunk_id
    output_dir = chunk_dir / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    result_path = chunk_dir / "chunk_results.tsv"
    if result_path.exists() and result_path.stat().st_size:
        cached = pd.read_csv(result_path, sep="\t", low_memory=False)
        cache_matches = (
            len(cached) == len(chunk)
            and cached["ddg_status"].eq("ok").all()
            and cached["foldx_calculation_key"].astype(str).tolist()
            == chunk["foldx_calculation_key"].astype(str).tolist()
            and cached["foldx_mutation"].astype(str).tolist()
            == chunk["foldx_mutation"].astype(str).tolist()
        )
        if cache_matches:
            return cached
    mutation_file = chunk_dir / "individual_list.txt"
    mutation_file.write_text("\n".join(chunk["foldx_mutation"]) + "\n", encoding="utf-8")
    build_log = chunk_dir / "build.log"
    for stale in output_dir.glob(f"*_{repaired_pdb.stem}.fxout"):
        stale.unlink()
    command = [
        str(Path(args.foldx).resolve()), "--command=BuildModel", f"--pdb={repaired_pdb.name}",
        f"--pdb-dir={repaired_pdb.parent.resolve()}", f"--output-dir={output_dir.resolve()}",
        f"--mutant-file={mutation_file.resolve()}", f"--rotabaseLocation={Path(args.rotabase).resolve()}",
        "--out-pdb=false", "--screen=false",
    ]
    code, command_text = run_command(command, build_log, args.build_timeout_seconds)
    result = chunk[["foldx_calculation_key", "foldx_mutation"]].copy()
    result["chunk_id"] = chunk_id
    result["repair_pdb"] = root_relative(repaired_pdb)
    result["repair_log"] = root_relative(repaired_pdb.parent.parent / "repair.log")
    result["build_log"] = root_relative(build_log)
    result["ddg_method"] = "FoldX_5.1_RepairPDB_BuildModel_AlphaFold_v6"
    if code != 0:
        result["ddg_status"] = "timeout" if code == 124 else "build_failed"
        result["ddg_missing_reason"] = command_text.replace("\n", " ")[-500:]
        result["ddg_kcal_mol"] = pd.NA
        result["ddg_abs"] = pd.NA
        result["dif_fxout"] = ""
    else:
        dif_path = output_dir / f"Dif_{repaired_pdb.stem}.fxout"
        try:
            values = parse_dif(dif_path)
            if len(values) != len(chunk):
                raise ValueError(f"Expected {len(chunk)} DDG rows, parsed {len(values)}")
            result["ddg_status"] = "ok"
            result["ddg_missing_reason"] = ""
            result["ddg_kcal_mol"] = values
            result["ddg_abs"] = [abs(value) for value in values]
            result["dif_fxout"] = root_relative(dif_path)
        except Exception as exc:
            result["ddg_status"] = "parse_failed"
            result["ddg_missing_reason"] = str(exc)[:500]
            result["ddg_kcal_mol"] = pd.NA
            result["ddg_abs"] = pd.NA
            result["dif_fxout"] = root_relative(dif_path) if dif_path.exists() else ""
    result.to_csv(result_path, sep="\t", index=False)
    return result


def load_cached_calculations(work_dir: Path) -> pd.DataFrame:
    cached_frames = []
    for path in work_dir.glob("chunks/*/chunk_results.tsv"):
        try:
            frame = pd.read_csv(path, sep="\t", low_memory=False)
        except Exception:
            continue
        required = {"foldx_calculation_key", "foldx_mutation", "ddg_status"}
        if required.issubset(frame.columns):
            cached_frames.append(frame[frame["ddg_status"].eq("ok")])
    if not cached_frames:
        return pd.DataFrame()
    return (
        pd.concat(cached_frames, ignore_index=True)
        .drop_duplicates("foldx_calculation_key", keep="last")
    )


def make_output(frame: pd.DataFrame, calculations: pd.DataFrame) -> pd.DataFrame:
    columns = {
        "audit_uniprot_accession": "uniprot_accession",
        "audit_protein_variant": "protein_variant",
        "audit_fragment_member": "alphafold_fragment_member",
        "audit_local_position": "alphafold_local_position",
        "structure_plddt": "alphafold_residue_plddt",
        "audit_fragment_edge_distance": "fragment_edge_distance",
    }
    output = frame.rename(columns=columns).merge(
        calculations.drop(columns=["foldx_mutation"], errors="ignore"),
        on="foldx_calculation_key",
        how="left",
        validate="many_to_one",
    )
    return output.reindex(columns=OUTPUT_COLUMNS)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", default=str(DEFAULT_TABLE))
    parser.add_argument("--audit", default=str(DEFAULT_AUDIT))
    parser.add_argument("--structure", default=str(DEFAULT_STRUCTURE))
    parser.add_argument("--pdb-dir", default=str(DEFAULT_PDB))
    parser.add_argument("--foldx", default=str(DEFAULT_FOLDX))
    parser.add_argument("--rotabase", default=str(DEFAULT_ROTABASE))
    parser.add_argument("--work-dir", default=str(DEFAULT_WORK))
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--summary", default=str(DEFAULT_SUMMARY))
    parser.add_argument("--min-plddt", type=float, default=70.0)
    parser.add_argument("--batch-size", type=int, default=50)
    parser.add_argument("--repair-workers", type=int, default=8)
    parser.add_argument("--build-workers", type=int, default=12)
    parser.add_argument("--repair-timeout-seconds", type=int, default=3600)
    parser.add_argument("--build-timeout-seconds", type=int, default=7200)
    parser.add_argument("--limit-fragments", type=int)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()

    frame = load_candidates(args)
    calculations = frame.drop_duplicates("foldx_calculation_key", keep="first").copy()
    members = sorted(calculations["audit_fragment_member"].unique())
    if args.limit_fragments is not None:
        members = members[: args.limit_fragments]
        calculations = calculations[calculations["audit_fragment_member"].isin(members)].copy()
        frame = frame[frame["foldx_calculation_key"].isin(calculations["foldx_calculation_key"])].copy()
    preflight = {
        "variant_rows": len(frame), "unique_calculations": len(calculations), "fragments": len(members),
        "batches": int(sum((len(group) + args.batch_size - 1) // args.batch_size for _, group in calculations.groupby("audit_fragment_member"))),
    }
    if args.preflight_only:
        print(json.dumps(preflight, indent=2))
        return 0

    repairs: dict[str, tuple[str, str]] = {}
    with ThreadPoolExecutor(max_workers=max(1, args.repair_workers)) as executor:
        futures = {executor.submit(repair_fragment, member, args): member for member in members}
        for future in as_completed(futures):
            member, status, repaired = future.result()
            repairs[member] = (status, repaired)
            print(f"repair {len(repairs)}/{len(members)} {status} {member}", flush=True)

    cached = load_cached_calculations(Path(args.work_dir))
    cached_keys = set(cached["foldx_calculation_key"].astype(str)) if not cached.empty else set()
    calculations_pending = calculations[
        ~calculations["foldx_calculation_key"].astype(str).isin(cached_keys)
    ].copy()
    build_jobs = []
    failed_rows = []
    for member, group in calculations_pending.groupby("audit_fragment_member", sort=True):
        repair_status, repaired = repairs[member]
        if repair_status != "ok":
            failed = group[["foldx_calculation_key", "foldx_mutation"]].copy()
            failed["ddg_status"] = "repair_failed"
            failed["ddg_missing_reason"] = f"RepairPDB failed for {member}"
            failed_rows.append(failed)
            continue
        group = group.reset_index(drop=True)
        for start in range(0, len(group), args.batch_size):
            chunk = group.iloc[start:start + args.batch_size].copy()
            digest = hashlib.sha1(
                "\n".join(chunk["foldx_calculation_key"].astype(str)).encode("utf-8")
            ).hexdigest()[:12]
            chunk_id = f"{safe(Path(member).stem)}_{digest}"
            build_jobs.append((chunk, Path(repaired), chunk_id))

    built = [cached] if not cached.empty else []
    completed_jobs = 0
    with ThreadPoolExecutor(max_workers=max(1, args.build_workers)) as executor:
        futures = {
            executor.submit(build_chunk, chunk, repaired, chunk_id, args): chunk_id
            for chunk, repaired, chunk_id in build_jobs
        }
        for future in as_completed(futures):
            built.append(future.result())
            completed_jobs += 1
            print(f"build {completed_jobs}/{len(build_jobs)} {futures[future]}", flush=True)
    all_results = built + failed_rows
    calculations_out = pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()
    output = make_output(frame, calculations_out)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)

    unique = output.drop_duplicates("foldx_calculation_key")
    summary = {
        **preflight,
        "status": "complete" if len(unique) == len(calculations) and unique["ddg_status"].isin(TERMINAL_STATUSES).all() else "partial",
        "status_counts_unique": dict(Counter(unique["ddg_status"].fillna("missing"))),
        "coverage_by_split": {
            str(split): {"rows": len(group), "ok": int(group["ddg_status"].eq("ok").sum())}
            for split, group in output.groupby("split_source_heldout", dropna=False)
        },
        "selection_rule": (
            "All canonical, non-synonymous missense substitutions in the primary-binary cohort "
            f"with validated protein/structure mapping and pLDDT>={args.min_plddt:g}; class values "
            "were not used to select within that cohort."
        ),
        "method": "FoldX 5.1 RepairPDB once per AlphaFold v6 fragment, then batched BuildModel; mutant PDB output disabled.",
        "foldx_binary": str(Path(args.foldx)), "rotabase": str(Path(args.rotabase)),
        "input_table": str(Path(args.table)), "mapping_audit": str(Path(args.audit)),
        "structure_features": str(Path(args.structure)), "output": str(Path(args.output)),
        "work_dir": str(Path(args.work_dir)), "min_plddt": args.min_plddt, "batch_size": args.batch_size,
        "cached_successful_calculations_reused": int(len(cached_keys & set(calculations["foldx_calculation_key"].astype(str)))),
        "new_calculations_submitted": int(len(calculations_pending)),
    }
    Path(args.summary).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
