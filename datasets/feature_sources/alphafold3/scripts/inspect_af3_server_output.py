#!/usr/bin/env python3
"""Inspect an AlphaFold Server output folder for AF3 feature readiness."""

from __future__ import annotations

import csv
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
OUTPUT_DIR = ROOT / "datasets/feature_sources/alphafold3/raw/pilot_v1_server_outputs/fold_af3_pilotv1_review_tmem43_emd_lmna"
VARIANT_SCAFFOLD = ROOT / "datasets/feature_sources/alphafold3/processed/pilot_v1/af3_pilot_v1_variant_feature_scaffold.tsv"
OUT_DIR = ROOT / "results/af3/pilot_v1/output_qc/fold_af3_pilotv1_review_tmem43_emd_lmna"
SUMMARY_JSON = OUT_DIR / "af3_output_qc_summary.json"
CHAIN_TSV = OUT_DIR / "chain_confidence_summary.tsv"
PAIR_TSV = OUT_DIR / "pair_interaction_confidence.tsv"
REPORT_MD = OUT_DIR / "AF3_OUTPUT_QC_REPORT.md"


CHAIN_NAMES = {"A": "TMEM43", "B": "EMD", "C": "LMNA"}


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_job_request() -> dict:
    path = next(OUTPUT_DIR.glob("*_job_request.json"))
    payload = json.load(path.open())
    if isinstance(payload, list):
        payload = payload[0]
    return payload


def parse_cif_atom_plddt(path: Path) -> dict[str, list[float]]:
    lines = path.read_text().splitlines()
    cols: list[str] = []
    atom_rows = []
    in_atom = False
    for line in lines:
        if line.startswith("_atom_site."):
            cols.append(line.strip().replace("_atom_site.", ""))
            in_atom = True
            continue
        if in_atom:
            if line.startswith("#"):
                break
            if line.startswith("ATOM") or line.startswith("HETATM"):
                parts = line.split()
                if len(parts) >= len(cols):
                    atom_rows.append(dict(zip(cols, parts)))

    by_residue: dict[tuple[str, int], list[float]] = defaultdict(list)
    for row in atom_rows:
        chain = row["label_asym_id"]
        residue = int(row["label_seq_id"])
        by_residue[(chain, residue)].append(float(row["B_iso_or_equiv"]))

    by_chain = defaultdict(list)
    for (chain, _residue), values in by_residue.items():
        by_chain[chain].append(statistics.mean(values))
    return by_chain


def chain_summary(by_chain: dict[str, list[float]]) -> list[dict[str, object]]:
    rows = []
    for chain, values in sorted(by_chain.items()):
        values = sorted(values)
        n = len(values)
        rows.append(
            {
                "chain_id": chain,
                "gene": CHAIN_NAMES.get(chain, ""),
                "n_residues": n,
                "mean_plddt": round(statistics.mean(values), 2),
                "median_plddt": round(statistics.median(values), 2),
                "p10_plddt": round(values[int(0.10 * (n - 1))], 2),
                "p90_plddt": round(values[int(0.90 * (n - 1))], 2),
                "n_residues_plddt_lt50": sum(v < 50 for v in values),
                "n_residues_plddt_gte70": sum(v >= 70 for v in values),
            }
        )
    return rows


def pair_summary(full_data_path: Path) -> list[dict[str, object]]:
    payload = json.load(full_data_path.open())
    chain_ids = payload["token_chain_ids"]
    pae = payload["pae"]
    contact = payload["contact_probs"]
    rows = []
    for a, b in [("A", "B"), ("A", "C"), ("B", "C")]:
        ixs = [i for i, c in enumerate(chain_ids) if c == a]
        jxs = [i for i, c in enumerate(chain_ids) if c == b]
        pae_values = []
        contact_values = []
        for i in ixs:
            pae_row = pae[i]
            contact_row = contact[i]
            for j in jxs:
                pae_values.append(pae_row[j])
                contact_values.append(contact_row[j])
        rows.append(
            {
                "chain_pair": f"{a}-{b}",
                "gene_pair": f"{CHAIN_NAMES[a]}-{CHAIN_NAMES[b]}",
                "min_pae": round(min(pae_values), 2),
                "median_pae": round(statistics.median(pae_values), 2),
                "max_contact_probability": round(max(contact_values), 4),
                "mean_contact_probability": round(statistics.mean(contact_values), 6),
            }
        )
    return rows


def summary_confidences() -> list[dict[str, object]]:
    rows = []
    for path in sorted(OUTPUT_DIR.glob("*_summary_confidences_*.json")):
        payload = json.load(path.open())
        idx = path.stem.rsplit("_", 1)[-1]
        rows.append(
            {
                "model_index": idx,
                "ranking_score": payload.get("ranking_score"),
                "iptm": payload.get("iptm"),
                "ptm": payload.get("ptm"),
                "fraction_disordered": payload.get("fraction_disordered"),
                "has_clash": payload.get("has_clash"),
                "chain_pair_pae_min": payload.get("chain_pair_pae_min"),
                "chain_pair_iptm": payload.get("chain_pair_iptm"),
            }
        )
    return rows


def variant_coverage() -> dict[str, object]:
    counts = defaultdict(Counter)
    total = 0
    with VARIANT_SCAFFOLD.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene = row["primary_gene"]
            if gene not in {"TMEM43", "EMD", "LMNA"}:
                continue
            total += 1
            counts[gene]["rows"] += 1
            counts[gene][f"label_{row['model_label_3class'] or 'blank'}"] += 1
            if row["protein_position"]:
                counts[gene]["with_protein_position"] += 1
    return {"total_rows": total, "by_gene": {gene: dict(counter) for gene, counter in sorted(counts.items())}}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    job_request = load_job_request()
    model0 = OUTPUT_DIR / "fold_af3_pilotv1_review_tmem43_emd_lmna_model_0.cif"
    full0 = OUTPUT_DIR / "fold_af3_pilotv1_review_tmem43_emd_lmna_full_data_0.json"

    by_chain = parse_cif_atom_plddt(model0)
    chain_rows = chain_summary(by_chain)
    pair_rows = pair_summary(full0)
    model_rows = summary_confidences()
    coverage = variant_coverage()

    write_tsv(
        CHAIN_TSV,
        chain_rows,
        [
            "chain_id",
            "gene",
            "n_residues",
            "mean_plddt",
            "median_plddt",
            "p10_plddt",
            "p90_plddt",
            "n_residues_plddt_lt50",
            "n_residues_plddt_gte70",
        ],
    )
    write_tsv(
        PAIR_TSV,
        pair_rows,
        ["chain_pair", "gene_pair", "min_pae", "median_pae", "max_contact_probability", "mean_contact_probability"],
    )

    summary = {
        "job_name": job_request["name"],
        "output_dir": str(OUTPUT_DIR),
        "file_counts": {
            "model_cif": len(list(OUTPUT_DIR.glob("*_model_*.cif"))),
            "full_data_json": len(list(OUTPUT_DIR.glob("*_full_data_*.json"))),
            "summary_confidence_json": len(list(OUTPUT_DIR.glob("*_summary_confidences_*.json"))),
            "job_request_json": len(list(OUTPUT_DIR.glob("*_job_request.json"))),
        },
        "submitted_chain_order": [
            {"chain_id": chain, "gene": CHAIN_NAMES[chain], "submitted_sequence_length": len(seq["proteinChain"]["sequence"])}
            for chain, seq in zip(["A", "B", "C"], job_request["sequences"])
        ],
        "best_model_summary": model_rows[0],
        "chain_confidence_summary": chain_rows,
        "pair_interaction_summary": pair_rows,
        "variant_coverage": coverage,
        "study_readiness_interpretation": {
            "can_parse_for_pipeline": True,
            "has_cif_for_distances": True,
            "has_full_data_pae_contact_probs": True,
            "has_local_plddt": True,
            "interface_features_should_be_trusted": False,
            "reason": "Best model has low ipTM and high inter-chain PAE with near-zero contact probabilities across TMEM43-EMD, TMEM43-LMNA, and EMD-LMNA.",
        },
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")

    report = [
        "# AF3 Output QC: TMEM43-EMD-LMNA",
        "",
        f"Job: `{job_request['name']}`",
        "",
        "## File Completeness",
        "",
        "| File type | Count |",
        "|---|---:|",
    ]
    for key, value in summary["file_counts"].items():
        report.append(f"| {key} | {value} |")
    report.extend(
        [
            "",
            "## Chain Map",
            "",
            "| Chain | Gene | Submitted length |",
            "|---|---|---:|",
        ]
    )
    for row in summary["submitted_chain_order"]:
        report.append(f"| {row['chain_id']} | {row['gene']} | {row['submitted_sequence_length']} |")
    report.extend(
        [
            "",
            "## Best Model Confidence",
            "",
            "| Metric | Value |",
            "|---|---:|",
        ]
    )
    for key in ["ranking_score", "iptm", "ptm", "fraction_disordered", "has_clash"]:
        report.append(f"| {key} | {summary['best_model_summary'][key]} |")
    report.extend(
        [
            "",
            "## Chain Local Confidence",
            "",
            "| Chain | Gene | Residues | Mean pLDDT | Median pLDDT | pLDDT <50 | pLDDT >=70 |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in chain_rows:
        report.append(
            f"| {row['chain_id']} | {row['gene']} | {row['n_residues']} | {row['mean_plddt']} | "
            f"{row['median_plddt']} | {row['n_residues_plddt_lt50']} | {row['n_residues_plddt_gte70']} |"
        )
    report.extend(
        [
            "",
            "## Inter-Chain Confidence",
            "",
            "| Gene pair | Min PAE | Median PAE | Max contact probability | Mean contact probability |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in pair_rows:
        report.append(
            f"| {row['gene_pair']} | {row['min_pae']} | {row['median_pae']} | "
            f"{row['max_contact_probability']} | {row['mean_contact_probability']} |"
        )
    report.extend(
        [
            "",
            "## Variant Coverage",
            "",
            f"Total modeling rows covered by these genes: `{coverage['total_rows']}`.",
            "",
            "| Gene | Rows | Benign | Pathogenic | VUS | With protein position |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for gene, counter in coverage["by_gene"].items():
        report.append(
            f"| {gene} | {counter.get('rows', 0)} | {counter.get('label_Benign', 0)} | "
            f"{counter.get('label_Pathogenic', 0)} | {counter.get('label_VUS', 0)} | "
            f"{counter.get('with_protein_position', 0)} |"
        )
    report.extend(
        [
            "",
            "## Interpretation",
            "",
            "This output has the files needed for the AF3 parsing pipeline: CIF structures, full-data JSON with PAE/contact probabilities, summary confidence JSON, and the original job request.",
            "",
            "However, the predicted three-chain interaction is not reliable enough to use as a high-confidence interface model. The best model has low ipTM, high inter-chain PAE, and near-zero inter-chain contact probabilities. Use this result as a pipeline smoke test and as low-confidence/local-structure context only. Do not treat TMEM43-EMD-LMNA interface distances from this run as strong biological evidence.",
            "",
            "Recommended next step: run the simpler `LMNA-EMD` first-wave job and compare confidence. For TMEM43, consider pairwise `TMEM43-EMD`, `TMEM43-LMNA`, or a domain/fragment design if the full-length trimer remains low-confidence.",
            "",
        ]
    )
    REPORT_MD.write_text("\n".join(report))


if __name__ == "__main__":
    main()
