#!/usr/bin/env python3
"""Analyze AlphaFold Server batch outputs and extract feature scaffolds."""

from __future__ import annotations

import csv
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[4]
DEFAULT_BATCH_DIR = (
    ROOT
    / "datasets/feature_sources/alphafold3/raw/pilot_v1_server_outputs/folds_2026_07_04_21_37"
)
PILOT_V1_MANIFEST = ROOT / "datasets/feature_sources/alphafold3/jobs/pilot_v1/af3_pilot_v1_submission_manifest.tsv"
PILOT_V2_MANIFEST = ROOT / "results/af3/pilot_v2_second20/af3_pilot_v2_second20_job_summary.tsv"
VARIANT_SCAFFOLD = (
    ROOT
    / "datasets/feature_sources/alphafold3/processed/pilot_v1/af3_pilot_v1_variant_feature_scaffold.tsv"
)


def resolve_batch_dir() -> Path:
    if len(sys.argv) > 1:
        path = Path(sys.argv[1]).expanduser()
        if not path.is_absolute():
            path = ROOT / path
        return path
    return DEFAULT_BATCH_DIR


BATCH_DIR = resolve_batch_dir()
BATCH_NAME = BATCH_DIR.name
OUT_DIR = ROOT / "results/af3/batch_output_qc" / BATCH_NAME
FEATURE_OUT = (
    ROOT
    / "datasets/feature_sources/alphafold3/processed"
    / BATCH_NAME
    / f"af3_{BATCH_NAME}_variant_features_long.tsv"
)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_json(path: Path) -> object:
    with path.open() as handle:
        return json.load(handle)


def read_manifest() -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    if PILOT_V1_MANIFEST.exists():
        df = pd.read_csv(PILOT_V1_MANIFEST, sep="\t")
        for _, row in df.iterrows():
            payload = row.to_dict()
            payload.setdefault("complex_group_id", payload.get("job_id", ""))
            rows[payload["job_id"]] = payload
    if PILOT_V2_MANIFEST.exists():
        df = pd.read_csv(PILOT_V2_MANIFEST, sep="\t")
        for _, row in df.iterrows():
            payload = row.to_dict()
            payload["complex_group_id"] = payload.get("job_id", "")
            payload["complex_name"] = payload.get("complex_name", payload.get("job_id", ""))
            rows[payload["job_id"]] = payload
    return rows


def ordered_chain_ids(full_data: dict) -> list[str]:
    seen = []
    for chain_id in full_data.get("token_chain_ids", []):
        if chain_id not in seen:
            seen.append(chain_id)
    return seen


def expand_protein_entities(job_request: dict, manifest_row: dict[str, str]) -> list[str]:
    manifest_entities = [x.strip() for x in str(manifest_row["protein_entities"]).split(",") if x.strip()]
    expanded = []
    manifest_index = 0
    for item in job_request.get("sequences", []):
        if "proteinChain" not in item:
            continue
        count = int(item["proteinChain"].get("count") or 1)
        for _ in range(count):
            if manifest_index < len(manifest_entities):
                expanded.append(manifest_entities[manifest_index])
            else:
                expanded.append(f"UNKNOWN_PROTEIN_{manifest_index + 1}")
            manifest_index += 1
    return expanded


def chain_map_for_job(job_request: dict, full_data: dict, manifest_row: dict[str, str]) -> dict[str, str]:
    chain_ids = ordered_chain_ids(full_data)
    protein_entities = expand_protein_entities(job_request, manifest_row)
    mapping = {}
    for chain_id, gene in zip(chain_ids, protein_entities):
        mapping[chain_id] = gene
    return mapping


def best_model_index(job_dir: Path) -> tuple[int, dict]:
    best_idx = -1
    best_payload: dict | None = None
    best_score = -math.inf
    for path in sorted(job_dir.glob("*_summary_confidences_*.json")):
        payload = load_json(path)
        idx = int(path.stem.rsplit("_", 1)[-1])
        score = payload.get("ranking_score")
        score_float = float(score) if score is not None else -math.inf
        if score_float > best_score:
            best_idx = idx
            best_payload = payload
            best_score = score_float
    if best_payload is None:
        raise FileNotFoundError(f"No summary confidence JSON found in {job_dir}")
    return best_idx, best_payload


def parse_cif_atoms(path: Path) -> list[dict[str, object]]:
    cols = []
    rows: list[dict[str, object]] = []
    in_atom_loop = False
    for line in path.read_text().splitlines():
        if line.startswith("_atom_site."):
            cols.append(line.strip().replace("_atom_site.", ""))
            in_atom_loop = True
            continue
        if in_atom_loop:
            if line.startswith("#"):
                if rows:
                    break
                continue
            if line.startswith("ATOM") or line.startswith("HETATM"):
                parts = line.split()
                if len(parts) < len(cols):
                    continue
                raw = dict(zip(cols, parts))
                chain = raw.get("label_asym_id", "")
                seq_id = raw.get("label_seq_id", "")
                if not chain or seq_id in {"", ".", "?"}:
                    continue
                try:
                    residue = int(float(str(seq_id)))
                    rows.append(
                        {
                            "chain_id": chain,
                            "residue": residue,
                            "atom_name": raw.get("label_atom_id", ""),
                            "x": float(raw["Cartn_x"]),
                            "y": float(raw["Cartn_y"]),
                            "z": float(raw["Cartn_z"]),
                            "plddt": float(raw["B_iso_or_equiv"]),
                            "group": raw.get("group_PDB", ""),
                        }
                    )
                except (KeyError, ValueError):
                    continue
    return rows


def residue_plddt_and_coords(atoms: list[dict[str, object]]) -> tuple[dict[tuple[str, int], float], dict[str, dict[int, list[tuple[float, float, float]]]]]:
    plddt_values: dict[tuple[str, int], list[float]] = defaultdict(list)
    coords: dict[str, dict[int, list[tuple[float, float, float]]]] = defaultdict(lambda: defaultdict(list))
    for atom in atoms:
        key = (str(atom["chain_id"]), int(atom["residue"]))
        plddt_values[key].append(float(atom["plddt"]))
        coords[str(atom["chain_id"])][int(atom["residue"])].append(
            (float(atom["x"]), float(atom["y"]), float(atom["z"]))
        )
    plddt = {key: statistics.mean(values) for key, values in plddt_values.items()}
    return plddt, coords


def p10(values: list[float]) -> float:
    if not values:
        return float("nan")
    values = sorted(values)
    return values[int(0.10 * (len(values) - 1))]


def pair_matrix_value(matrix: list, i: int, j: int) -> object:
    try:
        return matrix[i][j]
    except Exception:
        return None


def pair_metrics(
    full_data: dict, chain_a: str, chain_b: str, pair_iptm_matrix: list, pair_pae_min_matrix: list, chain_order: list[str]
) -> dict[str, object]:
    token_chains = full_data["token_chain_ids"]
    pae = full_data["pae"]
    contact = full_data["contact_probs"]
    ixs = [i for i, c in enumerate(token_chains) if c == chain_a]
    jxs = [i for i, c in enumerate(token_chains) if c == chain_b]
    pae_values = []
    contact_values = []
    for i in ixs:
        pae_row = pae[i]
        contact_row = contact[i]
        for j in jxs:
            pae_values.append(float(pae_row[j]))
            contact_values.append(float(contact_row[j]))
    idx_a = chain_order.index(chain_a) if chain_a in chain_order else -1
    idx_b = chain_order.index(chain_b) if chain_b in chain_order else -1
    pair_iptm = pair_matrix_value(pair_iptm_matrix, idx_a, idx_b) if idx_a >= 0 and idx_b >= 0 else None
    pair_pae_min = pair_matrix_value(pair_pae_min_matrix, idx_a, idx_b) if idx_a >= 0 and idx_b >= 0 else None
    min_pae = min(pae_values) if pae_values else None
    median_pae = statistics.median(pae_values) if pae_values else None
    max_contact = max(contact_values) if contact_values else None
    mean_contact = statistics.mean(contact_values) if contact_values else None
    support = classify_pair(pair_iptm, pair_pae_min, min_pae, median_pae, max_contact)
    return {
        "pair_iptm": pair_iptm,
        "pair_pae_min_summary": pair_pae_min,
        "min_pae": min_pae,
        "median_pae": median_pae,
        "p10_pae": p10(pae_values),
        "max_contact_probability": max_contact,
        "mean_contact_probability": mean_contact,
        "pair_support_class": support,
    }


def classify_pair(pair_iptm: object, pair_pae_min: object, min_pae: object, median_pae: object, max_contact: object) -> str:
    try:
        iptm = float(pair_iptm)
        pae_min = float(pair_pae_min if pair_pae_min is not None else min_pae)
        med = float(median_pae)
        contact = float(max_contact)
    except (TypeError, ValueError):
        return "not_applicable"
    if iptm >= 0.70 and pae_min <= 5 and contact >= 0.60:
        return "strong_supported_interface"
    if iptm >= 0.45 and pae_min <= 10 and contact >= 0.35:
        return "moderate_supported_interface"
    if iptm >= 0.30 and pae_min <= 15 and contact >= 0.20 and med <= 25:
        return "weak_possible_interface"
    return "unsupported_or_low_confidence"


def min_distances_to_partner(
    coords: dict[str, dict[int, list[tuple[float, float, float]]]], chain_a: str, chain_b: str
) -> dict[int, float]:
    partner_points = []
    for points in coords.get(chain_b, {}).values():
        partner_points.extend(points)
    if not partner_points:
        return {}
    tree = cKDTree(partner_points)
    distances = {}
    for residue, points in coords.get(chain_a, {}).items():
        if not points:
            continue
        dist, _idx = tree.query(points, k=1)
        if hasattr(dist, "tolist"):
            dist_values = dist.tolist()
            distances[residue] = float(min(dist_values))
        else:
            distances[residue] = float(dist)
    return distances


def clean_num(value: object, digits: int = 4) -> object:
    if value is None:
        return ""
    try:
        if math.isnan(float(value)):
            return ""
        return round(float(value), digits)
    except (TypeError, ValueError):
        return value


def protein_position_to_int(value: object) -> int | None:
    if value is None or str(value).strip() in {"", "-", "nan", "NaN"}:
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def analyze_job(job_dir: Path, manifest: dict[str, dict[str, str]]) -> tuple[list[dict], list[dict], list[dict], dict, dict]:
    job_req_path = next(job_dir.glob("*_job_request.json"))
    job_request = load_json(job_req_path)
    if isinstance(job_request, list):
        job_request = job_request[0]
    job_id = job_request["name"]
    manifest_row = manifest[job_id]
    best_idx, best_summary = best_model_index(job_dir)
    full_data = load_json(next(job_dir.glob(f"*_full_data_{best_idx}.json")))
    model_cif = next(job_dir.glob(f"*_model_{best_idx}.cif"))
    atoms = parse_cif_atoms(model_cif)
    residue_plddt, coords = residue_plddt_and_coords(atoms)
    chain_order = ordered_chain_ids(full_data)
    chain_to_gene = chain_map_for_job(job_request, full_data, manifest_row)

    chain_rows = []
    for chain_id, gene in chain_to_gene.items():
        residues = sorted(coords.get(chain_id, {}))
        values = [residue_plddt[(chain_id, r)] for r in residues if (chain_id, r) in residue_plddt]
        if not values:
            continue
        chain_rows.append(
            {
                "job_id": job_id,
                "job_folder": job_dir.name,
                "complex_group_id": manifest_row["complex_group_id"],
                "complex_name": manifest_row["complex_name"],
                "best_model_index": best_idx,
                "chain_id": chain_id,
                "gene": gene,
                "n_residues": len(values),
                "mean_plddt": clean_num(statistics.mean(values), 2),
                "median_plddt": clean_num(statistics.median(values), 2),
                "p10_plddt": clean_num(p10(values), 2),
                "n_residues_plddt_lt50": sum(v < 50 for v in values),
                "n_residues_plddt_gte70": sum(v >= 70 for v in values),
            }
        )

    pair_rows = []
    pair_distance_cache: dict[tuple[str, str], dict[int, float]] = {}
    chains = list(chain_to_gene)
    for i, chain_a in enumerate(chains):
        for chain_b in chains[i + 1 :]:
            metrics = pair_metrics(
                full_data,
                chain_a,
                chain_b,
                best_summary.get("chain_pair_iptm") or [],
                best_summary.get("chain_pair_pae_min") or [],
                chain_order,
            )
            d_ab = min_distances_to_partner(coords, chain_a, chain_b)
            d_ba = min_distances_to_partner(coords, chain_b, chain_a)
            pair_distance_cache[(chain_a, chain_b)] = d_ab
            pair_distance_cache[(chain_b, chain_a)] = d_ba
            min_distance = min([*d_ab.values(), *d_ba.values()]) if d_ab or d_ba else None
            n_res_5 = sum(d <= 5 for d in d_ab.values()) + sum(d <= 5 for d in d_ba.values())
            n_res_8 = sum(d <= 8 for d in d_ab.values()) + sum(d <= 8 for d in d_ba.values())
            pair_rows.append(
                {
                    "job_id": job_id,
                    "job_folder": job_dir.name,
                    "complex_group_id": manifest_row["complex_group_id"],
                    "complex_name": manifest_row["complex_name"],
                    "best_model_index": best_idx,
                    "chain_pair": f"{chain_a}-{chain_b}",
                    "gene_pair": f"{chain_to_gene[chain_a]}-{chain_to_gene[chain_b]}",
                    "pair_iptm": clean_num(metrics["pair_iptm"], 3),
                    "pair_pae_min_summary": clean_num(metrics["pair_pae_min_summary"], 2),
                    "min_pae": clean_num(metrics["min_pae"], 2),
                    "median_pae": clean_num(metrics["median_pae"], 2),
                    "p10_pae": clean_num(metrics["p10_pae"], 2),
                    "max_contact_probability": clean_num(metrics["max_contact_probability"], 4),
                    "mean_contact_probability": clean_num(metrics["mean_contact_probability"], 6),
                    "min_interchain_atom_distance_angstrom": clean_num(min_distance, 2),
                    "n_residues_at_interface_5A": n_res_5,
                    "n_residues_near_interface_8A": n_res_8,
                    "pair_support_class": metrics["pair_support_class"],
                }
            )

    usable_pairs = [r for r in pair_rows if r["pair_support_class"] in {"strong_supported_interface", "moderate_supported_interface"}]
    possible_pairs = [r for r in pair_rows if r["pair_support_class"] == "weak_possible_interface"]
    job_call = "usable_complex_features" if usable_pairs else ("local_or_exploratory_only" if possible_pairs else "not_interface_usable")
    job_summary = {
        "job_id": job_id,
        "job_folder": job_dir.name,
        "complex_group_id": manifest_row["complex_group_id"],
        "complex_name": manifest_row["complex_name"],
        "target_genes": manifest_row["target_genes"],
        "protein_entities": manifest_row["protein_entities"],
        "best_model_index": best_idx,
        "ranking_score": clean_num(best_summary.get("ranking_score"), 3),
        "iptm": clean_num(best_summary.get("iptm"), 3),
        "ptm": clean_num(best_summary.get("ptm"), 3),
        "fraction_disordered": clean_num(best_summary.get("fraction_disordered"), 3),
        "has_clash": clean_num(best_summary.get("has_clash"), 3),
        "n_chain_pairs": len(pair_rows),
        "n_strong_pairs": sum(r["pair_support_class"] == "strong_supported_interface" for r in pair_rows),
        "n_moderate_pairs": sum(r["pair_support_class"] == "moderate_supported_interface" for r in pair_rows),
        "n_weak_pairs": sum(r["pair_support_class"] == "weak_possible_interface" for r in pair_rows),
        "job_feature_call": job_call,
        "chain_map": ";".join(f"{c}:{g}" for c, g in chain_to_gene.items()),
    }
    aux = {
        "chain_to_gene": chain_to_gene,
        "residue_plddt": residue_plddt,
        "coords": coords,
        "pair_rows": pair_rows,
        "pair_distance_cache": pair_distance_cache,
        "manifest_row": manifest_row,
        "best_idx": best_idx,
    }
    return [job_summary], chain_rows, pair_rows, job_summary, aux


def build_variant_features(job_aux_list: list[dict]) -> list[dict]:
    scaffold = pd.read_csv(VARIANT_SCAFFOLD, sep="\t", low_memory=False)
    rows = []
    for aux in job_aux_list:
        chain_to_gene = aux["chain_to_gene"]
        gene_to_chains: dict[str, list[tuple[str, str]]] = defaultdict(list)
        target_genes = {
            x.strip()
            for x in str(aux["manifest_row"].get("target_genes", "")).split(",")
            if x.strip()
        }
        for chain_id, gene in chain_to_gene.items():
            gene_to_chains[gene].append((chain_id, gene))
            # CALM1/2/3 encode highly similar/identical calmodulin proteins in
            # this pilot design. Keep the physical AF3 chain labeled CALM1, but
            # allow CALM2/CALM3 variant rows to use the same residue context.
            if gene == "CALM1":
                for alias in ["CALM2", "CALM3"]:
                    if alias in target_genes:
                        gene_to_chains[alias].append((chain_id, gene))
        for gene, chains in gene_to_chains.items():
            sub = scaffold[scaffold["primary_gene"] == gene].copy()
            if sub.empty:
                continue
            for _, var in sub.iterrows():
                pos = protein_position_to_int(var.get("protein_position"))
                for chain_id, chain_gene in chains:
                    residue_plddt = aux["residue_plddt"].get((chain_id, pos)) if pos is not None else None
                    partner_summaries = []
                    nearest_partner_gene = ""
                    nearest_partner_chain = ""
                    nearest_partner_dist = None
                    best_pair_support = "not_applicable"
                    best_pair_iptm = ""
                    best_pair_min_pae = ""
                    best_pair_max_contact = ""
                    for partner_chain, partner_gene in chain_to_gene.items():
                        if partner_chain == chain_id:
                            continue
                        dist = aux["pair_distance_cache"].get((chain_id, partner_chain), {}).get(pos) if pos is not None else None
                        pair_key_1 = f"{chain_id}-{partner_chain}"
                        pair_key_2 = f"{partner_chain}-{chain_id}"
                        pair_row = next(
                            (r for r in aux["pair_rows"] if r["chain_pair"] in {pair_key_1, pair_key_2}),
                            None,
                        )
                        if pair_row:
                            support = pair_row["pair_support_class"]
                            partner_summaries.append(f"{partner_gene}:{support}:dist={clean_num(dist,2)}")
                            support_rank = {
                                "strong_supported_interface": 4,
                                "moderate_supported_interface": 3,
                                "weak_possible_interface": 2,
                                "unsupported_or_low_confidence": 1,
                                "not_applicable": 0,
                            }
                            if support_rank.get(support, 0) > support_rank.get(best_pair_support, 0):
                                best_pair_support = support
                                best_pair_iptm = pair_row["pair_iptm"]
                                best_pair_min_pae = pair_row["min_pae"]
                                best_pair_max_contact = pair_row["max_contact_probability"]
                        if dist is not None and (nearest_partner_dist is None or dist < nearest_partner_dist):
                            nearest_partner_dist = dist
                            nearest_partner_gene = partner_gene
                            nearest_partner_chain = partner_chain
                    missing_reason = ""
                    if pos is None:
                        missing_reason = "missing_protein_position"
                    elif residue_plddt is None:
                        missing_reason = "protein_position_not_in_af3_chain"
                    rows.append(
                        {
                            "variant_id": var.get("variant_id", ""),
                            "primary_gene": gene,
                            "model_label_3class": var.get("model_label_3class", ""),
                            "primary_model_inclusion": var.get("primary_model_inclusion", ""),
                            "sources": var.get("sources", ""),
                            "protein_position": "" if pos is None else pos,
                            "protein_ref_aa": var.get("protein_ref_aa", ""),
                            "protein_alt_aa": var.get("protein_alt_aa", ""),
                            "af3_job_id": aux["manifest_row"]["job_id"],
                            "af3_complex_group_id": aux["manifest_row"]["complex_group_id"],
                            "af3_complex_name": aux["manifest_row"]["complex_name"],
                            "af3_best_model_index": aux["best_idx"],
                            "af3_chain_id": chain_id,
                            "af3_chain_gene": chain_gene,
                            "af3_residue_plddt": clean_num(residue_plddt, 2),
                            "af3_nearest_partner_gene": nearest_partner_gene,
                            "af3_nearest_partner_chain": nearest_partner_chain,
                            "af3_min_distance_to_partner_angstrom": clean_num(nearest_partner_dist, 2),
                            "af3_within_partner_interface_5A": int(nearest_partner_dist is not None and nearest_partner_dist <= 5),
                            "af3_within_partner_interface_8A": int(nearest_partner_dist is not None and nearest_partner_dist <= 8),
                            "af3_best_pair_support_class": best_pair_support,
                            "af3_best_pair_iptm": best_pair_iptm,
                            "af3_best_pair_min_pae": best_pair_min_pae,
                            "af3_best_pair_max_contact_probability": best_pair_max_contact,
                            "af3_partner_distance_summary": ";".join(partner_summaries),
                            "af3_feature_missing_reason": missing_reason,
                        }
                    )
    return rows


def build_gene_feature_coverage(feature_rows: list[dict]) -> list[dict]:
    if not feature_rows:
        return []
    df = pd.DataFrame(feature_rows)
    covered = []
    for (job_id, gene), sub in df.groupby(["af3_job_id", "primary_gene"], dropna=False):
        usable = sub[sub["af3_feature_missing_reason"].fillna("") == ""]
        covered.append(
            {
                "af3_job_id": job_id,
                "primary_gene": gene,
                "feature_rows_long": len(sub),
                "rows_with_residue_feature": len(usable),
                "rows_with_strong_or_moderate_pair": int(
                    sub["af3_best_pair_support_class"]
                    .isin(["strong_supported_interface", "moderate_supported_interface"])
                    .sum()
                ),
                "rows_within_partner_interface_5A": int(sub["af3_within_partner_interface_5A"].sum()),
                "rows_within_partner_interface_8A": int(sub["af3_within_partner_interface_8A"].sum()),
                "missing_protein_position": int((sub["af3_feature_missing_reason"] == "missing_protein_position").sum()),
                "protein_position_not_in_af3_chain": int(
                    (sub["af3_feature_missing_reason"] == "protein_position_not_in_af3_chain").sum()
                ),
            }
        )
    return covered


def write_report(job_rows: list[dict], pair_rows: list[dict], chain_rows: list[dict], feature_rows: list[dict]) -> None:
    report = [
        "# AF3 Pilot V1 Batch Output QC",
        "",
        f"Batch folder: `{BATCH_DIR}`",
        "",
        "This report treats AF3 as an additive feature layer. It does not modify the existing variant registry or CatBoost feature matrices.",
        "",
        "## Job-Level Calls",
        "",
        "| Job | Proteins | Ranking | ipTM | pTM | Strong | Moderate | Weak | Feature call |",
        "|---|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in sorted(job_rows, key=lambda r: r["job_id"]):
        report.append(
            f"| {row['job_id']} | {row['protein_entities']} | {row['ranking_score']} | {row['iptm']} | {row['ptm']} | "
            f"{row['n_strong_pairs']} | {row['n_moderate_pairs']} | {row['n_weak_pairs']} | {row['job_feature_call']} |"
        )
    report.extend(
        [
            "",
            "## Supported Or Possible Pair Interfaces",
            "",
            "| Job | Pair | Pair ipTM | Min PAE | Max contact | Min distance A | Interface residues 5A | Support |",
            "|---|---|---:|---:|---:|---:|---:|---|",
        ]
    )
    show_pairs = [
        r
        for r in pair_rows
        if r["pair_support_class"]
        in {"strong_supported_interface", "moderate_supported_interface", "weak_possible_interface"}
    ]
    for row in sorted(show_pairs, key=lambda r: (r["job_id"], r["gene_pair"])):
        report.append(
            f"| {row['job_id']} | {row['gene_pair']} | {row['pair_iptm']} | {row['min_pae']} | "
            f"{row['max_contact_probability']} | {row['min_interchain_atom_distance_angstrom']} | "
            f"{row['n_residues_at_interface_5A']} | {row['pair_support_class']} |"
        )
    report.extend(
        [
            "",
            "## Variant Feature Extraction",
            "",
            f"Long AF3 feature rows extracted: `{len(feature_rows):,}`.",
            "",
            "A row is one variant-position in one AF3 chain/job. Genes present in multiple AF3 complexes deliberately have multiple rows so aggregation can be reviewed later.",
        ]
    )
    (OUT_DIR / "AF3_BATCH_OUTPUT_QC_REPORT.md").write_text("\n".join(report) + "\n")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = read_manifest()
    job_rows: list[dict] = []
    chain_rows: list[dict] = []
    pair_rows: list[dict] = []
    aux_list: list[dict] = []
    for job_dir in sorted(p for p in BATCH_DIR.iterdir() if p.is_dir()):
        jr, cr, pr, _summary, aux = analyze_job(job_dir, manifest)
        job_rows.extend(jr)
        chain_rows.extend(cr)
        pair_rows.extend(pr)
        aux_list.append(aux)
    feature_rows = build_variant_features(aux_list)
    gene_feature_rows = build_gene_feature_coverage(feature_rows)

    write_tsv(
        OUT_DIR / "af3_batch_job_qc_summary.tsv",
        job_rows,
        [
            "job_id",
            "job_folder",
            "complex_group_id",
            "complex_name",
            "target_genes",
            "protein_entities",
            "best_model_index",
            "ranking_score",
            "iptm",
            "ptm",
            "fraction_disordered",
            "has_clash",
            "n_chain_pairs",
            "n_strong_pairs",
            "n_moderate_pairs",
            "n_weak_pairs",
            "job_feature_call",
            "chain_map",
        ],
    )
    write_tsv(
        OUT_DIR / "af3_batch_chain_confidence.tsv",
        chain_rows,
        [
            "job_id",
            "job_folder",
            "complex_group_id",
            "complex_name",
            "best_model_index",
            "chain_id",
            "gene",
            "n_residues",
            "mean_plddt",
            "median_plddt",
            "p10_plddt",
            "n_residues_plddt_lt50",
            "n_residues_plddt_gte70",
        ],
    )
    write_tsv(
        OUT_DIR / "af3_batch_pair_interaction_qc.tsv",
        pair_rows,
        [
            "job_id",
            "job_folder",
            "complex_group_id",
            "complex_name",
            "best_model_index",
            "chain_pair",
            "gene_pair",
            "pair_iptm",
            "pair_pae_min_summary",
            "min_pae",
            "median_pae",
            "p10_pae",
            "max_contact_probability",
            "mean_contact_probability",
            "min_interchain_atom_distance_angstrom",
            "n_residues_at_interface_5A",
            "n_residues_near_interface_8A",
            "pair_support_class",
        ],
    )
    feature_fields = [
        "variant_id",
        "primary_gene",
        "model_label_3class",
        "primary_model_inclusion",
        "sources",
        "protein_position",
        "protein_ref_aa",
        "protein_alt_aa",
        "af3_job_id",
        "af3_complex_group_id",
        "af3_complex_name",
        "af3_best_model_index",
        "af3_chain_id",
        "af3_chain_gene",
        "af3_residue_plddt",
        "af3_nearest_partner_gene",
        "af3_nearest_partner_chain",
        "af3_min_distance_to_partner_angstrom",
        "af3_within_partner_interface_5A",
        "af3_within_partner_interface_8A",
        "af3_best_pair_support_class",
        "af3_best_pair_iptm",
        "af3_best_pair_min_pae",
        "af3_best_pair_max_contact_probability",
        "af3_partner_distance_summary",
        "af3_feature_missing_reason",
    ]
    write_tsv(FEATURE_OUT, feature_rows, feature_fields)
    write_tsv(
        OUT_DIR / "af3_batch_gene_feature_coverage.tsv",
        gene_feature_rows,
        [
            "af3_job_id",
            "primary_gene",
            "feature_rows_long",
            "rows_with_residue_feature",
            "rows_with_strong_or_moderate_pair",
            "rows_within_partner_interface_5A",
            "rows_within_partner_interface_8A",
            "missing_protein_position",
            "protein_position_not_in_af3_chain",
        ],
    )
    summary = {
        "batch_dir": str(BATCH_DIR),
        "n_jobs": len(job_rows),
        "n_chain_rows": len(chain_rows),
        "n_pair_rows": len(pair_rows),
        "n_variant_feature_rows_long": len(feature_rows),
        "n_usable_jobs": sum(r["job_feature_call"] == "usable_complex_features" for r in job_rows),
        "n_local_or_exploratory_jobs": sum(r["job_feature_call"] == "local_or_exploratory_only" for r in job_rows),
        "n_not_interface_usable_jobs": sum(r["job_feature_call"] == "not_interface_usable" for r in job_rows),
        "outputs": {
            "job_qc": str(OUT_DIR / "af3_batch_job_qc_summary.tsv"),
            "chain_confidence": str(OUT_DIR / "af3_batch_chain_confidence.tsv"),
            "pair_qc": str(OUT_DIR / "af3_batch_pair_interaction_qc.tsv"),
            "gene_feature_coverage": str(OUT_DIR / "af3_batch_gene_feature_coverage.tsv"),
            "variant_features_long": str(FEATURE_OUT),
            "report": str(OUT_DIR / "AF3_BATCH_OUTPUT_QC_REPORT.md"),
        },
    }
    (OUT_DIR / "af3_batch_output_qc_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_report(job_rows, pair_rows, chain_rows, feature_rows)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
