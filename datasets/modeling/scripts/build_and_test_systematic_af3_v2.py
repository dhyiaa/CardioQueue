#!/usr/bin/env python3
"""Build systematic AF3 V2 features from gene_001-gene_015 and train CatBoost.

This is an additive experiment. It does not replace the baseline modeling table.
The primary V2 model uses structural AF3 features only and keeps AF3 job/gene
provenance out of the training features.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from catboost import CatBoostClassifier, Pool
from scipy.spatial import cKDTree
from sklearn.metrics import average_precision_score, confusion_matrix, roc_auc_score


ROOT = Path(__file__).resolve().parents[3]
SYSTEMATIC_DIR = ROOT / "AF3-systematic protein/partner list"
BASE_TABLE = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
TRAIN_SCRIPT = ROOT / "datasets/modeling/scripts/train_primary_binary_catboost.py"
PROCESSED_DIR = ROOT / "datasets/feature_sources/alphafold3/processed/systematic_v2_gene001_015"
OUT = ROOT / "results/model_performance/systematic_AF3_v2_gene001_015"
RESULTS_META = ROOT / "results/af3/systematic_v2_gene001_015"
MERGED_TABLE = ROOT / "datasets/modeling/interim/modeling_table_with_splits_weights_plus_systematicAF3_v2_gene001_015.tsv"
BASELINE_PRED = (
    ROOT
    / "results/models/primary_binary_catboost_v0_strict_source_heldout_hiro_emerge_rescued"
    / "primary_binary_predictions_split_source_heldout.tsv"
)
SPLIT_COL = "split_source_heldout"
SEED = 20260702
GENE_DIR_RE = re.compile(r"gene_(\d{3})_([A-Z0-9-]+)$")
STRUCTURAL_FEATURE_PREFIX = "systematicAF3_"


def load_primary_trainer():
    spec = importlib.util.spec_from_file_location("primary_binary_trainer", TRAIN_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {TRAIN_SCRIPT}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


TRAINER = load_primary_trainer()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_json(path: Path) -> object:
    with path.open() as handle:
        return json.load(handle)


def best_model_index(job_dir: Path) -> tuple[int, dict]:
    best_idx = -1
    best_payload: dict | None = None
    best_score = -math.inf
    for path in sorted(job_dir.glob("*_summary_confidences_*.json")):
        payload = load_json(path)
        idx = int(path.stem.rsplit("_", 1)[-1])
        score = float(payload.get("ranking_score") or -math.inf)
        if score > best_score:
            best_score = score
            best_idx = idx
            best_payload = payload
    if best_payload is None:
        raise FileNotFoundError(f"No summary_confidences JSON in {job_dir}")
    return best_idx, best_payload


def ordered_chain_ids(full_data: dict) -> list[str]:
    seen = []
    for chain_id in full_data.get("token_chain_ids", []):
        if chain_id not in seen:
            seen.append(chain_id)
    return seen


def parse_cif_atoms(path: Path) -> list[dict[str, object]]:
    cols: list[str] = []
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
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            parts = line.split()
            if len(parts) < len(cols):
                continue
            raw = dict(zip(cols, parts))
            chain = raw.get("label_asym_id", "")
            seq_id = raw.get("label_seq_id", "")
            if not chain or seq_id in {"", ".", "?"}:
                continue
            try:
                rows.append(
                    {
                        "chain_id": chain,
                        "residue": int(float(str(seq_id))),
                        "atom_name": raw.get("label_atom_id", ""),
                        "x": float(raw["Cartn_x"]),
                        "y": float(raw["Cartn_y"]),
                        "z": float(raw["Cartn_z"]),
                        "plddt": float(raw["B_iso_or_equiv"]),
                    }
                )
            except (KeyError, ValueError):
                continue
    return rows


def residue_plddt_and_coords(
    atoms: list[dict[str, object]],
) -> tuple[dict[tuple[str, int], float], dict[str, dict[int, list[tuple[float, float, float]]]]]:
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


def min_distances_to_partner(
    coords: dict[str, dict[int, list[tuple[float, float, float]]]], chain_a: str, partner_chains: list[str]
) -> dict[int, float]:
    partner_points = []
    for chain_b in partner_chains:
        for points in coords.get(chain_b, {}).values():
            partner_points.extend(points)
    if not partner_points:
        return {}
    tree = cKDTree(partner_points)
    out = {}
    for residue, points in coords.get(chain_a, {}).items():
        if not points:
            continue
        distances, _ = tree.query(points, k=1)
        out[residue] = float(np.min(distances))
    return out


def matrix_value(matrix: list, i: int, j: int) -> object:
    try:
        return matrix[i][j]
    except Exception:
        return None


def pair_metrics(full_data: dict, summary: dict, chain_a: str, chain_b: str, chain_order: list[str]) -> dict[str, object]:
    token_chains = full_data["token_chain_ids"]
    pae = full_data["pae"]
    contact = full_data["contact_probs"]
    ixs = [i for i, c in enumerate(token_chains) if c == chain_a]
    jxs = [i for i, c in enumerate(token_chains) if c == chain_b]
    pae_values = []
    contact_values = []
    for i in ixs:
        for j in jxs:
            pae_values.append(float(pae[i][j]))
            contact_values.append(float(contact[i][j]))
    idx_a = chain_order.index(chain_a) if chain_a in chain_order else -1
    idx_b = chain_order.index(chain_b) if chain_b in chain_order else -1
    pair_iptm = matrix_value(summary.get("chain_pair_iptm", []), idx_a, idx_b) if idx_a >= 0 and idx_b >= 0 else None
    pair_pae_min = (
        matrix_value(summary.get("chain_pair_pae_min", []), idx_a, idx_b) if idx_a >= 0 and idx_b >= 0 else None
    )
    min_pae = min(pae_values) if pae_values else np.nan
    median_pae = statistics.median(pae_values) if pae_values else np.nan
    max_contact = max(contact_values) if contact_values else np.nan
    mean_contact = statistics.mean(contact_values) if contact_values else np.nan
    return {
        "pair_iptm": float(pair_iptm) if pair_iptm is not None else np.nan,
        "pair_pae_min_summary": float(pair_pae_min) if pair_pae_min is not None else np.nan,
        "min_pae": float(min_pae),
        "median_pae": float(median_pae),
        "max_contact_probability": float(max_contact),
        "mean_contact_probability": float(mean_contact),
    }


def classify_pair(row: dict[str, object]) -> str:
    iptm = float(row.get("pair_iptm", np.nan))
    pae_min = float(row.get("pair_pae_min_summary", np.nan))
    if np.isnan(pae_min):
        pae_min = float(row.get("min_pae", np.nan))
    max_contact = float(row.get("max_contact_probability", np.nan))
    median_pae = float(row.get("median_pae", np.nan))
    if iptm >= 0.70 and pae_min <= 5 and max_contact >= 0.60:
        return "strong_supported_interface"
    if iptm >= 0.60 and pae_min <= 8 and max_contact >= 0.45:
        return "exploratory_supported_interface"
    if iptm >= 0.45 and pae_min <= 10 and max_contact >= 0.35:
        return "moderate_possible_interface"
    if iptm >= 0.30 and pae_min <= 15 and max_contact >= 0.20 and median_pae <= 25:
        return "weak_possible_interface"
    return "unsupported_or_low_confidence"


def support_rank(value: object) -> int:
    text = str(value)
    if text.startswith("strong"):
        return 3
    if text.startswith("exploratory"):
        return 2
    if text.startswith("moderate"):
        return 1
    return 0


def is_trusted_support(value: object) -> bool:
    return support_rank(value) >= 2


def first_positive_position(value: object) -> float:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return np.nan
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "<na>", ".", "-"}:
        return np.nan
    candidates = []
    for token in re.findall(r"-?\d+(?:\.\d+)?", text):
        try:
            val = float(token)
        except ValueError:
            continue
        if val > 0:
            candidates.append(val)
    if not candidates:
        return np.nan
    return float(candidates[0])


def resolve_variant_protein_position(row: pd.Series) -> tuple[float, str]:
    # Prefer the locally curated protein parser when present, but rescue ClinVar
    # and other rows with VEP/dbNSFP/FoldX positions. This is critical for AF3:
    # a variant can only be mapped to an AF3 structure if its residue position is
    # known.
    for col in [
        "protein_position",
        "vep_protein_position",
        "dbnsfp_protein_position_raw",
        "dbnsfp_aapos",
        "foldx_ddg_protein_position",
    ]:
        if col not in row:
            continue
        pos = first_positive_position(row[col])
        if pd.notna(pos):
            return pos, col
    return np.nan, "missing"


def parse_segment_text(text: str, fallback_gene: str | None = None) -> list[dict[str, object]]:
    segments: list[dict[str, object]] = []
    for gene, start, end, xcount in re.findall(r"\b([A-Z0-9-]+)\s+(\d+)\s*-\s*(\d+)(?:\s*x\s*(\d+))?", text):
        count = int(xcount) if xcount else 1
        for _ in range(count):
            segments.append({"gene": gene, "start": int(start), "end": int(end)})
    if not segments and fallback_gene:
        for start, end, xcount in re.findall(r"(\d+)\s*-\s*(\d+)(?:\s*x\s*(\d+))?", text):
            count = int(xcount) if xcount else 1
            for _ in range(count):
                segments.append({"gene": fallback_gene, "start": int(start), "end": int(end)})
    if not segments and fallback_gene and re.search(r"full[- ]length", text, flags=re.IGNORECASE):
        segments.append({"gene": fallback_gene, "start": 1, "end": None})
    return segments


def parse_partner_segments(text: str) -> list[dict[str, object]]:
    segments = parse_segment_text(text)
    if segments:
        return segments
    blocked = {
        "CA",
        "DNA",
        "RNA",
        "IONS",
        "ION",
        "FULL",
        "LENGTH",
        "SELF",
        "DRAFT",
        "OR",
        "AND",
        "N",
        "C",
    }
    for token in re.findall(r"\b[A-Z][A-Z0-9]{2,9}\b", text):
        if token in blocked:
            continue
        segments.append({"gene": token, "start": 1, "end": None})
    return segments


def normalized_tokens(text: object) -> set[str]:
    tokens = set(re.findall(r"[A-Z0-9]+", str(text).upper()))
    return {
        t
        for t in tokens
        if t
        and t
        not in {
            "AF3",
            "SYS",
            "DRAFT",
            "TOP3",
            "D",
            "G",
            "FULL",
            "LENGTH",
            "CA",
            "CA2",
            "CA4",
            "CA16",
            "PLUS",
        }
    }


def design_row_for_job(design_id: str, designs: list[dict[str, str]]) -> dict[str, str]:
    if not designs:
        return {}
    lower_lookup = {}
    for row in designs:
        for key in ["design_id", "short_design", "design_label"]:
            if row.get(key):
                lower_lookup[str(row[key]).lower()] = row
    if design_id.lower() in lower_lookup:
        return lower_lookup[design_id.lower()]
    job_tokens = normalized_tokens(design_id)
    best_row: dict[str, str] = {}
    best_score = 0
    for row in designs:
        row_text = " ".join(str(row.get(k, "")) for k in ["design_id", "short_design", "design_label", "entities"])
        row_tokens = normalized_tokens(row_text)
        if not row_tokens:
            continue
        score = len(job_tokens & row_tokens)
        if score > best_score:
            best_score = score
            best_row = row
    return best_row if best_score >= 2 else {}


def fallback_target_segment_from_job_id(design_id: str, target_gene: str, sequence_length: int) -> list[dict[str, object]]:
    text = design_id.upper()
    escaped = re.escape(target_gene.upper())
    match = re.search(rf"{escaped}.*?(\d+)_(\d+)", text)
    if match:
        return [{"gene": target_gene, "start": int(match.group(1)), "end": int(match.group(2))}]
    return [{"gene": target_gene, "start": 1, "end": sequence_length}]


def infer_chain_segments(
    job_request: dict, design_row: dict[str, str], target_gene: str, design_id: str
) -> dict[str, dict[str, object]]:
    chain_ids = []
    chain_lengths = []
    for item in job_request.get("sequences", []):
        if "proteinChain" not in item:
            continue
        count = int(item["proteinChain"].get("count") or 1)
        seq_len = len(item["proteinChain"].get("sequence", ""))
        for _ in range(count):
            chain_lengths.append(seq_len)
    chain_ids = [chr(ord("A") + i) for i in range(len(chain_lengths))]

    target_text = design_row.get("target_residues") or design_row.get("residues") or ""
    target_segments = parse_segment_text(target_text, fallback_gene=target_gene)
    if not target_segments and chain_lengths:
        target_segments = fallback_target_segment_from_job_id(design_id, target_gene, chain_lengths[0])
    partner_text = design_row.get("partner") or design_row.get("entities") or ""
    partner_segments = parse_partner_segments(partner_text)
    segments = target_segments + partner_segments

    out: dict[str, dict[str, object]] = {}
    for i, chain_id in enumerate(chain_ids):
        seq_len = chain_lengths[i]
        seg = segments[i].copy() if i < len(segments) else {"gene": target_gene if i == 0 else "UNKNOWN_PROTEIN", "start": 1, "end": None}
        seg["start"] = int(seg.get("start") or 1)
        seg["end"] = int(seg.get("end") or (seg["start"] + seq_len - 1))
        seg["sequence_length"] = seq_len
        out[chain_id] = seg
    return out


def find_gene_dirs() -> list[Path]:
    dirs = []
    for path in SYSTEMATIC_DIR.iterdir():
        if not path.is_dir():
            continue
        m = GENE_DIR_RE.match(path.name)
        if not m:
            continue
        if 1 <= int(m.group(1)) <= 15:
            dirs.append(path)
    return sorted(dirs)


def collect_job_qc() -> tuple[pd.DataFrame, pd.DataFrame]:
    job_rows = []
    pair_rows = []
    for gene_dir in find_gene_dirs():
        m = GENE_DIR_RE.match(gene_dir.name)
        assert m
        gene_index = int(m.group(1))
        target_gene = m.group(2)
        design_path = gene_dir / "04_af3_design_decisions.tsv"
        designs = read_tsv(design_path) if design_path.exists() else []
        for job_dir in sorted((gene_dir / "af3_job_results").glob("*/*")):
            if not job_dir.is_dir():
                continue
            try:
                best_idx, summary = best_model_index(job_dir)
                full_data = load_json(next(job_dir.glob(f"*_full_data_{best_idx}.json")))
                job_request_raw = load_json(next(job_dir.glob("*_job_request.json")))
                job_request = job_request_raw[0] if isinstance(job_request_raw, list) else job_request_raw
            except Exception as exc:
                job_rows.append(
                    {
                        "gene_index": gene_index,
                        "target_gene": target_gene,
                        "job_id": job_dir.name,
                        "job_dir": str(job_dir.relative_to(ROOT)),
                        "parse_status": "failed",
                        "parse_error": str(exc),
                    }
                )
                continue
            design_id = str(job_request.get("name") or job_dir.name)
            design_row = design_row_for_job(design_id, designs)
            chain_segments = infer_chain_segments(job_request, design_row, target_gene, design_id)
            chain_order = ordered_chain_ids(full_data)
            protein_chain_order = [c for c in chain_order if c in chain_segments]
            job_rows.append(
                {
                    "gene_index": gene_index,
                    "target_gene": target_gene,
                    "job_id": design_id,
                    "job_dir": str(job_dir.relative_to(ROOT)),
                    "best_model_index": best_idx,
                    "ranking_score": summary.get("ranking_score"),
                    "iptm": summary.get("iptm"),
                    "ptm": summary.get("ptm"),
                    "fraction_disordered": summary.get("fraction_disordered"),
                    "has_clash": summary.get("has_clash"),
                    "chain_segments_json": json.dumps(chain_segments, sort_keys=True),
                    "parse_status": "ok",
                    "parse_error": "",
                }
            )
            for i, chain_a in enumerate(protein_chain_order):
                for chain_b in protein_chain_order[i + 1 :]:
                    metrics = pair_metrics(full_data, summary, chain_a, chain_b, chain_order)
                    support = classify_pair(metrics)
                    seg_a = chain_segments.get(chain_a, {})
                    seg_b = chain_segments.get(chain_b, {})
                    pair_rows.append(
                        {
                            "gene_index": gene_index,
                            "target_gene": target_gene,
                            "job_id": design_id,
                            "chain_a": chain_a,
                            "chain_b": chain_b,
                            "gene_a": seg_a.get("gene", ""),
                            "gene_b": seg_b.get("gene", ""),
                            "start_a": seg_a.get("start", ""),
                            "end_a": seg_a.get("end", ""),
                            "start_b": seg_b.get("start", ""),
                            "end_b": seg_b.get("end", ""),
                            "pair_support_class": support,
                            **metrics,
                        }
                    )
    return pd.DataFrame(job_rows), pd.DataFrame(pair_rows)


def build_variant_long_features(base: pd.DataFrame, job_qc: pd.DataFrame, pair_qc: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keep_cols = [
        "variant_id",
        "primary_gene",
        "model_label_3class",
        "primary_model_inclusion",
        "sources",
        "protein_position",
        "vep_protein_position",
        "dbnsfp_protein_position_raw",
        "dbnsfp_aapos",
        "foldx_ddg_protein_position",
        "protein_ref_aa",
        "protein_alt_aa",
    ]
    base_positions = base[[c for c in keep_cols if c in base.columns]].copy()
    resolved = base_positions.apply(resolve_variant_protein_position, axis=1, result_type="expand")
    base_positions["systematicAF3_variant_protein_position"] = resolved[0]
    base_positions["systematicAF3_variant_protein_position_source"] = resolved[1]
    base_positions = base_positions[base_positions["systematicAF3_variant_protein_position"].notna()].copy()
    pair_lookup = pair_qc.copy()
    pair_lookup["support_rank"] = pair_lookup["pair_support_class"].map(support_rank)

    for _, job in job_qc[job_qc["parse_status"].eq("ok")].iterrows():
        job_dir = ROOT / str(job["job_dir"])
        best_idx = int(job["best_model_index"])
        full_data = load_json(next(job_dir.glob(f"*_full_data_{best_idx}.json")))
        cif = next(job_dir.glob(f"*_model_{best_idx}.cif"))
        plddt, coords = residue_plddt_and_coords(parse_cif_atoms(cif))
        chain_segments = json.loads(job["chain_segments_json"])
        chain_order = ordered_chain_ids(full_data)
        for chain_id, seg in chain_segments.items():
            gene = str(seg.get("gene", ""))
            if not gene or gene == "UNKNOWN_PROTEIN":
                continue
            start = int(seg.get("start") or 1)
            end = int(seg.get("end") or start)
            gene_variants = base_positions[
                base_positions["primary_gene"].astype(str).eq(gene)
                & base_positions["systematicAF3_variant_protein_position"].between(start, end, inclusive="both")
            ].copy()
            if gene_variants.empty:
                continue
            partner_chains = [c for c in chain_segments if c != chain_id]
            dists = min_distances_to_partner(coords, chain_id, partner_chains)
            pair_sub = pair_lookup[
                (pair_lookup["job_id"].astype(str).eq(str(job["job_id"])))
                & ((pair_lookup["chain_a"].eq(chain_id)) | (pair_lookup["chain_b"].eq(chain_id)))
            ].copy()
            trusted_pairs = pair_sub[pair_sub["support_rank"] >= 2].copy()
            best_pairs = trusted_pairs if not trusted_pairs.empty else pair_sub
            if not best_pairs.empty:
                best_pair = best_pairs.sort_values(
                    ["support_rank", "pair_iptm", "max_contact_probability"],
                    ascending=[False, False, False],
                ).iloc[0]
                best_support = best_pair["pair_support_class"]
                best_iptm = best_pair["pair_iptm"]
                best_pae = best_pair["pair_pae_min_summary"]
                best_contact = best_pair["max_contact_probability"]
                nearest_partner_gene = best_pair["gene_b"] if best_pair["chain_a"] == chain_id else best_pair["gene_a"]
            else:
                best_support = "not_applicable"
                best_iptm = np.nan
                best_pae = np.nan
                best_contact = np.nan
                nearest_partner_gene = ""
            for _, var in gene_variants.iterrows():
                local_res = int(var["systematicAF3_variant_protein_position"] - start + 1)
                dist = dists.get(local_res, np.nan)
                rows.append(
                    {
                        "variant_id": var["variant_id"],
                        "primary_gene": var["primary_gene"],
                        "model_label_3class": var["model_label_3class"],
                        "primary_model_inclusion": var["primary_model_inclusion"],
                        "sources": var["sources"],
                        "protein_position": var["systematicAF3_variant_protein_position"],
                        "systematicAF3_variant_protein_position_source": var[
                            "systematicAF3_variant_protein_position_source"
                        ],
                        "protein_ref_aa": var["protein_ref_aa"],
                        "protein_alt_aa": var["protein_alt_aa"],
                        "systematicAF3_job_id": job["job_id"],
                        "systematicAF3_gene_index": job["gene_index"],
                        "systematicAF3_chain_id": chain_id,
                        "systematicAF3_chain_gene": gene,
                        "systematicAF3_chain_start": start,
                        "systematicAF3_chain_end": end,
                        "systematicAF3_local_residue": local_res,
                        "systematicAF3_residue_plddt": plddt.get((chain_id, local_res), np.nan),
                        "systematicAF3_min_distance_to_partner_angstrom": dist,
                        "systematicAF3_within_partner_interface_5A": int(pd.notna(dist) and dist <= 5),
                        "systematicAF3_within_partner_interface_8A": int(pd.notna(dist) and dist <= 8),
                        "systematicAF3_best_pair_support_class": best_support,
                        "systematicAF3_best_pair_support_rank": support_rank(best_support),
                        "systematicAF3_best_pair_iptm": best_iptm,
                        "systematicAF3_best_pair_min_pae": best_pae,
                        "systematicAF3_best_pair_max_contact_probability": best_contact,
                        "systematicAF3_nearest_partner_gene": nearest_partner_gene,
                        "systematicAF3_feature_missing_reason": "",
                    }
                )
    return pd.DataFrame(rows)


def aggregate_variant_features(long_df: pd.DataFrame) -> pd.DataFrame:
    if long_df.empty:
        return pd.DataFrame({"variant_id": []})
    numeric_cols = [
        "systematicAF3_residue_plddt",
        "systematicAF3_min_distance_to_partner_angstrom",
        "systematicAF3_best_pair_iptm",
        "systematicAF3_best_pair_min_pae",
        "systematicAF3_best_pair_max_contact_probability",
    ]
    for col in numeric_cols:
        long_df[col] = pd.to_numeric(long_df[col], errors="coerce")
    rows = []
    for variant_id, group in long_df.groupby("variant_id", dropna=False):
        ranks = group["systematicAF3_best_pair_support_rank"].astype(int)
        trusted = group[ranks >= 2]
        source = trusted if not trusted.empty else group
        rows.append(
            {
                "variant_id": variant_id,
                "systematicAF3_status": "has_systematicAF3_output",
                "systematicAF3_n_long_rows": int(len(group)),
                "systematicAF3_n_jobs": int(group["systematicAF3_job_id"].nunique()),
                "systematicAF3_best_pair_support_rank": int(ranks.max()),
                "systematicAF3_has_exploratory_or_strong_pair": int(ranks.max() >= 2),
                "systematicAF3_has_strong_pair": int(ranks.max() >= 3),
                "systematicAF3_best_residue_plddt": float(source["systematicAF3_residue_plddt"].max())
                if source["systematicAF3_residue_plddt"].notna().any()
                else np.nan,
                "systematicAF3_mean_residue_plddt": float(source["systematicAF3_residue_plddt"].mean())
                if source["systematicAF3_residue_plddt"].notna().any()
                else np.nan,
                "systematicAF3_min_distance_to_partner_angstrom": float(
                    source["systematicAF3_min_distance_to_partner_angstrom"].min()
                )
                if source["systematicAF3_min_distance_to_partner_angstrom"].notna().any()
                else np.nan,
                "systematicAF3_within_partner_interface_5A": int(
                    pd.to_numeric(source["systematicAF3_within_partner_interface_5A"], errors="coerce").fillna(0).max()
                ),
                "systematicAF3_within_partner_interface_8A": int(
                    pd.to_numeric(source["systematicAF3_within_partner_interface_8A"], errors="coerce").fillna(0).max()
                ),
                "systematicAF3_best_pair_iptm": float(source["systematicAF3_best_pair_iptm"].max())
                if source["systematicAF3_best_pair_iptm"].notna().any()
                else np.nan,
                "systematicAF3_best_pair_min_pae": float(source["systematicAF3_best_pair_min_pae"].min())
                if source["systematicAF3_best_pair_min_pae"].notna().any()
                else np.nan,
                "systematicAF3_best_pair_max_contact_probability": float(
                    source["systematicAF3_best_pair_max_contact_probability"].max()
                )
                if source["systematicAF3_best_pair_max_contact_probability"].notna().any()
                else np.nan,
                "systematicAF3_jobs": ";".join(sorted(map(str, group["systematicAF3_job_id"].dropna().unique()))),
                "systematicAF3_nearest_partner_genes": ";".join(
                    sorted(map(str, source["systematicAF3_nearest_partner_gene"].dropna().unique()))
                ),
                "systematicAF3_layer_note": "systematic_AF3_gene001_015_draft_layer",
            }
        )
    out = pd.DataFrame(rows)
    int_cols = [
        "systematicAF3_n_long_rows",
        "systematicAF3_n_jobs",
        "systematicAF3_best_pair_support_rank",
        "systematicAF3_has_exploratory_or_strong_pair",
        "systematicAF3_has_strong_pair",
        "systematicAF3_within_partner_interface_5A",
        "systematicAF3_within_partner_interface_8A",
    ]
    out[int_cols] = out[int_cols].fillna(0).astype(int)
    return out


def merge_features(base: pd.DataFrame, af3: pd.DataFrame) -> pd.DataFrame:
    merged = base.merge(af3, on="variant_id", how="left")
    for col in [
        "systematicAF3_n_long_rows",
        "systematicAF3_n_jobs",
        "systematicAF3_best_pair_support_rank",
        "systematicAF3_has_exploratory_or_strong_pair",
        "systematicAF3_has_strong_pair",
        "systematicAF3_within_partner_interface_5A",
        "systematicAF3_within_partner_interface_8A",
    ]:
        merged[col] = merged[col].fillna(0).astype(int)
    text_defaults = {
        "systematicAF3_status": "no_systematicAF3_output",
        "systematicAF3_jobs": "",
        "systematicAF3_nearest_partner_genes": "",
        "systematicAF3_layer_note": "systematic_AF3_gene001_015_draft_layer",
    }
    for col, value in text_defaults.items():
        merged[col] = merged[col].fillna(value)
    return merged


def prepare_features(df: pd.DataFrame, include_coverage_flags: bool = False) -> tuple[pd.DataFrame, list[str], list[str], list[str]]:
    feature_cols = [c for c in df.columns if not TRAINER.should_drop(c)]
    # Keep job names, partner gene names, and AF3 coverage/provenance out of the main V2 model.
    blocked = {
        "systematicAF3_status",
        "systematicAF3_jobs",
        "systematicAF3_nearest_partner_genes",
        "systematicAF3_layer_note",
        "systematicAF3_n_long_rows",
        "systematicAF3_n_jobs",
    }
    if not include_coverage_flags:
        blocked.update(
            {
                "systematicAF3_has_exploratory_or_strong_pair",
                "systematicAF3_has_strong_pair",
                "systematicAF3_best_pair_support_rank",
            }
        )
    feature_cols = [c for c in feature_cols if c not in blocked]
    x = df[feature_cols].copy()
    cat_cols = []
    numeric_cols = []
    for col in x.columns:
        if x[col].dtype == "object" or str(x[col].dtype).startswith("string") or x[col].dtype == bool:
            cat_cols.append(col)
            x[col] = x[col].astype("string").fillna("__MISSING__").astype(str)
        else:
            numeric_cols.append(col)
            x[col] = pd.to_numeric(x[col], errors="coerce")
    return x, feature_cols, cat_cols, numeric_cols


def metric_block(y_true: np.ndarray, proba: np.ndarray, threshold: float = 0.5) -> dict[str, float | int]:
    pred = (proba >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, pred, labels=[0, 1]).ravel()
    return {
        "rows": int(len(y_true)),
        "positives_pathogenic": int(y_true.sum()),
        "negatives_benign": int((y_true == 0).sum()),
        "auroc": float(roc_auc_score(y_true, proba)) if len(np.unique(y_true)) == 2 else np.nan,
        "auprc": float(average_precision_score(y_true, proba)) if len(np.unique(y_true)) == 2 else np.nan,
        "tp": int(tp),
        "fp": int(fp),
        "tn": int(tn),
        "fn": int(fn),
        "sensitivity": float(tp / (tp + fn)) if (tp + fn) else np.nan,
        "specificity": float(tn / (tn + fp)) if (tn + fp) else np.nan,
        "ppv": float(tp / (tp + fp)) if (tp + fp) else np.nan,
        "npv": float(tn / (tn + fn)) if (tn + fn) else np.nan,
    }


def train_model(df: pd.DataFrame, model_name: str, include_coverage_flags: bool = False) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    train_df = df[df["training_slice_primary_binary"].astype(str).str.lower().eq("true")].copy()
    x, feature_cols, cat_cols, numeric_cols = prepare_features(train_df, include_coverage_flags=include_coverage_flags)
    y = train_df["primary_binary_label"].astype(int).to_numpy()
    weights = train_df["sample_weight"].astype(float).to_numpy()
    split = train_df[SPLIT_COL].astype(str)
    train_mask = split.eq("train").to_numpy()
    val_mask = split.eq("validation").to_numpy()
    model = CatBoostClassifier(
        iterations=1500,
        learning_rate=0.035,
        depth=6,
        loss_function="Logloss",
        eval_metric="PRAUC",
        random_seed=SEED,
        l2_leaf_reg=6,
        od_type="Iter",
        od_wait=80,
        allow_writing_files=False,
        verbose=100,
    )
    model.fit(
        Pool(x.loc[train_mask], y[train_mask], weight=weights[train_mask], cat_features=cat_cols),
        eval_set=Pool(x.loc[val_mask], y[val_mask], weight=weights[val_mask], cat_features=cat_cols),
        use_best_model=True,
    )
    pred_frames = []
    metric_rows = []
    for split_name in ["train", "validation", "external_cardioboost", "external_emerge", "external_hiro"]:
        mask = split.eq(split_name).to_numpy()
        if not mask.any():
            continue
        proba = model.predict_proba(Pool(x.loc[mask], cat_features=cat_cols))[:, 1]
        y_true = y[mask]
        for scope_name, scope_mask in {
            "all_rows": np.ones(len(y_true), dtype=bool),
            "systematicAF3_any_output_rows": train_df.loc[mask, "systematicAF3_status"].eq("has_systematicAF3_output").to_numpy(),
            "systematicAF3_trusted_pair_rows": train_df.loc[mask, "systematicAF3_has_exploratory_or_strong_pair"].eq(1).to_numpy(),
            "systematicAF3_strong_pair_rows": train_df.loc[mask, "systematicAF3_has_strong_pair"].eq(1).to_numpy(),
        }.items():
            if not scope_mask.any():
                continue
            block = metric_block(y_true[scope_mask], proba[scope_mask])
            block.update({"model": model_name, "subset": split_name, "scope": scope_name})
            metric_rows.append(block)
        tmp = train_df.loc[
            mask,
            [
                "variant_id",
                "primary_gene",
                "model_label_3class",
                SPLIT_COL,
                "systematicAF3_status",
                "systematicAF3_has_exploratory_or_strong_pair",
                "systematicAF3_has_strong_pair",
                "systematicAF3_best_pair_iptm",
                "systematicAF3_min_distance_to_partner_angstrom",
                "systematicAF3_within_partner_interface_5A",
                "systematicAF3_within_partner_interface_8A",
            ],
        ].copy()
        tmp["y_true"] = y_true
        tmp["pathogenic_probability"] = proba
        pred_frames.append(tmp)
    predictions = pd.concat(pred_frames, ignore_index=True)
    metrics = pd.DataFrame(metric_rows)
    model.save_model(str(OUT / f"{model_name}_catboost_model.cbm"))
    pd.DataFrame({"feature": feature_cols, "importance": model.get_feature_importance(type="FeatureImportance")}).sort_values(
        "importance", ascending=False
    ).to_csv(OUT / f"{model_name}_feature_importance.tsv", sep="\t", index=False)
    meta = {
        "model": model_name,
        "input_table": str(MERGED_TABLE.relative_to(ROOT)),
        "split_col": SPLIT_COL,
        "n_features": len(feature_cols),
        "n_categorical_features": len(cat_cols),
        "n_numeric_features": len(numeric_cols),
        "n_systematicAF3_features_used": sum(c.startswith(STRUCTURAL_FEATURE_PREFIX) for c in feature_cols),
        "systematicAF3_features_used": [c for c in feature_cols if c.startswith(STRUCTURAL_FEATURE_PREFIX)],
        "include_coverage_flags": include_coverage_flags,
        "best_iteration": int(model.get_best_iteration() or 0),
    }
    return predictions, metrics, meta


def baseline_metrics(af3_features: pd.DataFrame) -> pd.DataFrame:
    if not BASELINE_PRED.exists():
        return pd.DataFrame()
    pred = pd.read_csv(BASELINE_PRED, sep="\t")
    cols = ["variant_id", "systematicAF3_status", "systematicAF3_has_exploratory_or_strong_pair", "systematicAF3_has_strong_pair"]
    pred = pred.merge(af3_features[cols].drop_duplicates("variant_id"), on="variant_id", how="left")
    pred["systematicAF3_status"] = pred["systematicAF3_status"].fillna("no_systematicAF3_output")
    for col in ["systematicAF3_has_exploratory_or_strong_pair", "systematicAF3_has_strong_pair"]:
        pred[col] = pred[col].fillna(0).astype(int)
    rows = []
    split_col = SPLIT_COL if SPLIT_COL in pred.columns else [c for c in pred.columns if c.startswith("split_")][0]
    for split_name, group in pred.groupby(split_col):
        y_true = group["y_true"].astype(int).to_numpy()
        proba = group["pathogenic_probability"].astype(float).to_numpy()
        for scope_name, scope_mask in {
            "all_rows": np.ones(len(group), dtype=bool),
            "systematicAF3_any_output_rows": group["systematicAF3_status"].eq("has_systematicAF3_output").to_numpy(),
            "systematicAF3_trusted_pair_rows": group["systematicAF3_has_exploratory_or_strong_pair"].eq(1).to_numpy(),
            "systematicAF3_strong_pair_rows": group["systematicAF3_has_strong_pair"].eq(1).to_numpy(),
        }.items():
            if not scope_mask.any():
                continue
            block = metric_block(y_true[scope_mask], proba[scope_mask])
            block.update({"model": "baseline_primary_binary_existing", "subset": split_name, "scope": scope_name})
            rows.append(block)
    return pd.DataFrame(rows)


def write_report(coverage: pd.DataFrame, metrics: pd.DataFrame, job_qc: pd.DataFrame, pair_qc: pd.DataFrame, metas: list[dict]) -> None:
    def md_table(df: pd.DataFrame) -> str:
        if df.empty:
            return "_No rows._"
        lines = ["| " + " | ".join(df.columns) + " |", "| " + " | ".join(["---"] * len(df.columns)) + " |"]
        for _, row in df.iterrows():
            lines.append("| " + " | ".join(str(row[c]) for c in df.columns) + " |")
        return "\n".join(lines)

    display = metrics[["model", "subset", "scope", "rows", "auroc", "auprc", "sensitivity", "specificity", "ppv", "npv"]].copy()
    for col in ["auroc", "auprc", "sensitivity", "specificity", "ppv", "npv"]:
        display[col] = display[col].map(lambda x: "" if pd.isna(x) else f"{x:.4f}")

    cov = coverage.copy()
    for col in [c for c in cov.columns if c.startswith("pct_")]:
        cov[col] = cov[col].map(lambda x: "" if pd.isna(x) else f"{x:.2f}%")

    top_pairs = pair_qc.sort_values(["pair_support_class", "pair_iptm"], ascending=[True, False]).copy()
    top_pairs = top_pairs[
        [
            "target_gene",
            "job_id",
            "gene_a",
            "gene_b",
            "pair_support_class",
            "pair_iptm",
            "pair_pae_min_summary",
            "max_contact_probability",
        ]
    ].head(30)
    for col in ["pair_iptm", "pair_pae_min_summary", "max_contact_probability"]:
        top_pairs[col] = pd.to_numeric(top_pairs[col], errors="coerce").map(lambda x: "" if pd.isna(x) else f"{x:.3f}")

    feature_text = "\n\n".join(
        [
            f"### {m['model']}\n\n```text\n" + "\n".join(m["systematicAF3_features_used"]) + "\n```"
            for m in metas
        ]
    )
    report = f"""# Systematic AF3 V2 Gene 001-015 CatBoost Test

This is an additive V2 experiment using systematic AlphaFold 3 outputs from genes 1-15 only. The established baseline tables and models were not modified.

## AF3 Interaction Threshold Policy

- Strong interface: pair ipTM >= 0.70, pair/min PAE <= 5, max contact probability >= 0.60.
- Exploratory supported interface: pair ipTM >= 0.60, pair/min PAE <= 8, max contact probability >= 0.45.
- The main model uses structural columns only. AF3 job IDs, partner names, status text, and count/provenance columns are excluded from training.
- Unsupported jobs are retained in QC files but are not treated as positive interaction evidence.

## Coverage

{md_table(cov)}

## Performance

{md_table(display)}

## Pair QC Preview

{md_table(top_pairs)}

## AF3 Features Used

{feature_text}

## Outputs

- Aggregated AF3 variant features: `datasets/feature_sources/alphafold3/processed/systematic_v2_gene001_015/systematicAF3_v2_variant_features_aggregated.tsv`
- Long AF3 variant features: `datasets/feature_sources/alphafold3/processed/systematic_v2_gene001_015/systematicAF3_v2_variant_features_long.tsv`
- Merged modeling table: `datasets/modeling/interim/modeling_table_with_splits_weights_plus_systematicAF3_v2_gene001_015.tsv`
- Metrics: `results/model_performance/systematic_AF3_v2_gene001_015/systematic_AF3_v2_metrics.tsv`
"""
    (OUT / "SYSTEMATIC_AF3_V2_GENE001_015_REPORT.md").write_text(report)


def main() -> None:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    RESULTS_META.mkdir(parents=True, exist_ok=True)
    base = pd.read_csv(BASE_TABLE, sep="\t", low_memory=False)

    job_qc, pair_qc = collect_job_qc()
    job_qc.to_csv(RESULTS_META / "systematicAF3_v2_job_qc.tsv", sep="\t", index=False)
    pair_qc.to_csv(RESULTS_META / "systematicAF3_v2_pair_qc.tsv", sep="\t", index=False)

    long_df = build_variant_long_features(base, job_qc, pair_qc)
    long_df.to_csv(PROCESSED_DIR / "systematicAF3_v2_variant_features_long.tsv", sep="\t", index=False)
    agg = aggregate_variant_features(long_df)
    agg.to_csv(PROCESSED_DIR / "systematicAF3_v2_variant_features_aggregated.tsv", sep="\t", index=False)

    merged = merge_features(base, agg)
    merged.to_csv(MERGED_TABLE, sep="\t", index=False)

    scopes = {
        "all_modeling_rows": merged.index == merged.index,
        "primary_binary_rows": merged["training_slice_primary_binary"].astype(str).str.lower().eq("true"),
        "vus_rows": merged["model_label_3class"].astype(str).eq("VUS"),
    }
    cov_rows = []
    for scope, mask in scopes.items():
        sub = merged[mask]
        cov_rows.append(
            {
                "scope": scope,
                "rows": len(sub),
                "genes": sub["primary_gene"].nunique(dropna=True),
                "rows_with_any_systematicAF3_output": int(sub["systematicAF3_status"].eq("has_systematicAF3_output").sum()),
                "rows_with_exploratory_or_strong_pair": int(sub["systematicAF3_has_exploratory_or_strong_pair"].sum()),
                "rows_with_strong_pair": int(sub["systematicAF3_has_strong_pair"].sum()),
                "pct_any_output": 100 * sub["systematicAF3_status"].eq("has_systematicAF3_output").mean() if len(sub) else np.nan,
                "pct_exploratory_or_strong_pair": 100 * sub["systematicAF3_has_exploratory_or_strong_pair"].mean()
                if len(sub)
                else np.nan,
                "pct_strong_pair": 100 * sub["systematicAF3_has_strong_pair"].mean() if len(sub) else np.nan,
            }
        )
    coverage = pd.DataFrame(cov_rows)
    coverage.to_csv(RESULTS_META / "systematicAF3_v2_coverage_summary.tsv", sep="\t", index=False)
    by_gene = (
        merged.groupby("primary_gene", dropna=False)
        .agg(
            rows=("variant_id", "count"),
            binary_rows=("training_slice_primary_binary", lambda s: s.astype(str).str.lower().eq("true").sum()),
            vus_rows=("model_label_3class", lambda s: s.astype(str).eq("VUS").sum()),
            rows_with_any_systematicAF3_output=("systematicAF3_status", lambda s: s.eq("has_systematicAF3_output").sum()),
            rows_with_exploratory_or_strong_pair=("systematicAF3_has_exploratory_or_strong_pair", "sum"),
            rows_with_strong_pair=("systematicAF3_has_strong_pair", "sum"),
        )
        .reset_index()
        .sort_values(["rows_with_exploratory_or_strong_pair", "rows_with_any_systematicAF3_output"], ascending=False)
    )
    by_gene.to_csv(RESULTS_META / "systematicAF3_v2_coverage_by_gene.tsv", sep="\t", index=False)

    main_pred, main_metrics, main_meta = train_model(
        merged, "systematicAF3_v2_primary_binary_structuralOnly", include_coverage_flags=False
    )
    main_pred.to_csv(OUT / "systematicAF3_v2_primary_binary_structuralOnly_predictions.tsv", sep="\t", index=False)
    flag_pred, flag_metrics, flag_meta = train_model(
        merged, "systematicAF3_v2_primary_binary_structuralPlusFlags", include_coverage_flags=True
    )
    flag_pred.to_csv(OUT / "systematicAF3_v2_primary_binary_structuralPlusFlags_predictions.tsv", sep="\t", index=False)
    base_metrics = baseline_metrics(agg)
    metrics = pd.concat([base_metrics, main_metrics, flag_metrics], ignore_index=True)
    metrics.to_csv(OUT / "systematic_AF3_v2_metrics.tsv", sep="\t", index=False)
    manifest = {
        "base_table": str(BASE_TABLE.relative_to(ROOT)),
        "merged_table": str(MERGED_TABLE.relative_to(ROOT)),
        "gene_scope": "gene_001 through gene_015 systematic AF3 folders",
        "trusted_threshold": "pair ipTM >= 0.60, pair/min PAE <= 8, max contact probability >= 0.45",
        "strong_threshold": "pair ipTM >= 0.70, pair/min PAE <= 5, max contact probability >= 0.60",
        "models": [main_meta, flag_meta],
        "coverage": coverage.to_dict(orient="records"),
    }
    (OUT / "systematic_AF3_v2_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    write_report(coverage, metrics, job_qc, pair_qc, [main_meta, flag_meta])
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
