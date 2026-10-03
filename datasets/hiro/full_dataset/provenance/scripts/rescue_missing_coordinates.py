#!/usr/bin/env python3
"""Create coordinate rescue candidates for HiRO rows missing chrom-pos-ref-alt.

This script is intentionally non-destructive. It reads the frozen HiRO model
input table, searches local evidence sources, and writes a review table of
candidate GRCh38 coordinates for rows whose coordinate join key is incomplete.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


MISSING_TOKENS = {"", "missing", "nan", "none", "na", "n/a", "-", "."}
COORD_COLUMNS = ("chrom", "pos", "ref", "alt")


def is_present(value: object) -> bool:
    return str(value or "").strip().lower() not in MISSING_TOKENS


def has_complete_coords(row: dict[str, object]) -> bool:
    return all(is_present(row.get(col, "")) for col in COORD_COLUMNS)


def coord_status(row: dict[str, object]) -> str:
    present = [is_present(row.get(col, "")) for col in COORD_COLUMNS]
    if all(present):
        return "complete"
    if any(present):
        return "partial"
    return "missing"


def clean(value: object) -> str:
    text = str(value or "").strip()
    return "" if text.lower() in MISSING_TOKENS else text


def normalize_hgvs(value: str) -> str:
    return re.sub(r"\s+", "", value.strip().lower())


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def split_gene_symbols(value: str) -> set[str]:
    return {part.strip().upper() for part in re.split(r"[;,]", value or "") if part.strip()}


def variant_text_terms(row: dict[str, str]) -> set[str]:
    terms: set[str] = set()

    for col in ("hgvs_c", "hgvs_p", "hgvs_p_one_letter", "variant_key"):
        value = clean(row.get(col, ""))
        if not value:
            continue
        for part in re.split(r"[:;,|]", value):
            part = part.strip()
            if part.lower().startswith(("c.", "p.")):
                terms.add(normalize_hgvs(part))

    gene_refgene = clean(row.get("gene_refgene", ""))
    if gene_refgene:
        for match in re.finditer(r"([cp]\.[^,;|]+)", gene_refgene, flags=re.IGNORECASE):
            term = match.group(1).strip()
            term = re.sub(r":.*$", "", term)
            if term.lower().startswith(("c.", "p.")):
                terms.add(normalize_hgvs(term))

    return {term for term in terms if len(term) >= 4}


def term_has_specific_allele(term: str) -> bool:
    lowered = term.lower()
    return any(token in lowered for token in (">", "del", "dup", "ins"))


def term_type(term: str) -> str:
    if term.startswith("c."):
        return "cdna"
    if term.startswith("p."):
        return "protein"
    return "other"


def load_existing_annotation_candidates(paths: Iterable[Path], target_uids: set[str]) -> dict[str, list[dict[str, str]]]:
    candidates: dict[str, list[dict[str, str]]] = defaultdict(list)
    for path in paths:
        if not path.exists():
            continue
        for row in read_csv_rows(path):
            uid = str(row.get("variant_uid", ""))
            if uid in target_uids and has_complete_coords(row):
                candidate = dict(row)
                candidate["_evidence_source"] = path.name
                candidates[uid].append(candidate)
    return candidates


def load_cache_candidates(path: Path, target_uids: set[str]) -> dict[str, list[dict[str, str]]]:
    candidates: dict[str, list[dict[str, str]]] = defaultdict(list)
    if not path.exists():
        return candidates
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            uid = str(row.get("variant_uid", ""))
            if uid in target_uids and has_complete_coords(row):
                candidate = {key: "" if value is None else str(value) for key, value in row.items()}
                candidate["_evidence_source"] = path.name
                candidates[uid].append(candidate)
    return candidates


def build_clinvar_candidates(clinvar_path: Path, queries: dict[str, dict[str, object]]) -> dict[str, list[dict[str, str]]]:
    by_gene: dict[str, list[tuple[str, set[str]]]] = defaultdict(list)
    for uid, payload in queries.items():
        for gene in payload["genes"]:
            by_gene[gene].append((uid, payload["terms"]))

    matches: dict[str, list[dict[str, str]]] = defaultdict(list)
    if not clinvar_path.exists() or not by_gene:
        return matches

    with gzip.open(clinvar_path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("Assembly") != "GRCh38":
                continue
            gene = str(row.get("GeneSymbol", "")).upper()
            if gene not in by_gene:
                continue

            name_norm = normalize_hgvs(row.get("Name", ""))
            for uid, terms in by_gene[gene]:
                hit_terms = sorted(term for term in terms if term and term in name_norm)
                if not hit_terms:
                    continue
                candidate = {
                    "_evidence_source": "clinvar_variant_summary",
                    "variant_uid": uid,
                    "gene": gene,
                    "chrom": row.get("Chromosome", ""),
                    "pos": row.get("PositionVCF", ""),
                    "ref": row.get("ReferenceAlleleVCF", ""),
                    "alt": row.get("AlternateAlleleVCF", ""),
                    "clinvar_variation_id": row.get("VariationID", ""),
                    "clinvar_name": row.get("Name", ""),
                    "clinvar_classification": row.get("ClinicalSignificance", ""),
                    "clinvar_review_status": row.get("ReviewStatus", ""),
                    "clinvar_matched_terms": ";".join(hit_terms),
                    "clinvar_matched_term_types": ";".join(sorted({term_type(term) for term in hit_terms})),
                    "clinvar_specific_allele_match": str(any(term_has_specific_allele(term) for term in hit_terms)).lower(),
                }
                if has_complete_coords(candidate):
                    matches[uid].append(candidate)
    return matches


def choose_candidate(candidates: list[dict[str, str]]) -> tuple[str, dict[str, str] | None, str]:
    if not candidates:
        return "not_rescued", None, "no local exact HGVS coordinate candidate"

    coords = {(clean(c.get("chrom")), clean(c.get("pos")), clean(c.get("ref")), clean(c.get("alt"))) for c in candidates}
    if len(coords) != 1:
        return "conflicting_candidates", None, f"{len(coords)} coordinate candidates"

    candidate = candidates[0]
    source_names = sorted({c.get("_evidence_source", "") for c in candidates})
    if any(source != "clinvar_variant_summary" for source in source_names):
        return "rescued_local_annotation", candidate, ";".join(source_names)

    matched_types = set(";".join(c.get("clinvar_matched_term_types", "") for c in candidates).split(";"))
    specific = any(c.get("clinvar_specific_allele_match") == "true" for c in candidates)
    if "cdna" in matched_types and specific:
        return "rescued_clinvar_hgvs_cdna", candidate, f"{len(candidates)} ClinVar exact HGVS candidate(s)"
    return "review_candidate_clinvar_hgvs", candidate, "ClinVar match lacks a specific cDNA allele-change term"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hiro-input", default="datasets/hiro/full_dataset/model_inputs/ml_baseline_features.csv")
    parser.add_argument("--clinvar", default="datasets/clinvar/variant_summary.txt.gz")
    parser.add_argument("--annotation-dir", default="datasets/hiro/full_dataset/annotations")
    parser.add_argument("--output", default="datasets/hiro/full_dataset/interim/hiro_coordinate_rescue_candidates.tsv")
    parser.add_argument("--summary", default="datasets/hiro/full_dataset/interim/hiro_coordinate_rescue_summary.json")
    args = parser.parse_args()

    hiro_path = Path(args.hiro_input)
    annotation_dir = Path(args.annotation_dir)
    output_path = Path(args.output)
    summary_path = Path(args.summary)

    rows = read_csv_rows(hiro_path)
    target_rows = [row for row in rows if not has_complete_coords(row)]
    target_uids = {row["variant_uid"] for row in target_rows}

    queries = {
        row["variant_uid"]: {
            "genes": split_gene_symbols(row.get("gene", "")),
            "terms": variant_text_terms(row),
        }
        for row in target_rows
    }

    local_candidates = load_existing_annotation_candidates(
        [
            annotation_dir / "variant_annotations.csv",
            annotation_dir / "annotated_cohort.csv",
        ],
        target_uids,
    )
    cache_candidates = load_cache_candidates(annotation_dir / "annotation_cache.jsonl", target_uids)
    clinvar_candidates = build_clinvar_candidates(Path(args.clinvar), queries)

    output_rows: list[dict[str, object]] = []
    status_counts: Counter[str] = Counter()
    confidence_counts: Counter[str] = Counter()

    for row in target_rows:
        uid = row["variant_uid"]
        candidates = []
        candidates.extend(local_candidates.get(uid, []))
        candidates.extend(cache_candidates.get(uid, []))
        candidates.extend(clinvar_candidates.get(uid, []))

        status, chosen, note = choose_candidate(candidates)
        status_counts[status] += 1
        confidence = "none"
        if status in {"rescued_local_annotation", "rescued_clinvar_hgvs_cdna"}:
            confidence = "high"
        elif status == "review_candidate_clinvar_hgvs":
            confidence = "review"
        elif status == "conflicting_candidates":
            confidence = "conflict"
        confidence_counts[confidence] += 1

        original = {f"original_{col}": clean(row.get(col, "")) for col in COORD_COLUMNS}
        rescued = {f"rescued_{col}": "" for col in COORD_COLUMNS}
        evidence = {
            "rescue_source": "",
            "rescue_clinvar_variation_id": "",
            "rescue_clinvar_name": "",
            "rescue_clinvar_classification": "",
            "rescue_clinvar_review_status": "",
            "rescue_matched_terms": "",
            "rescue_matched_term_types": "",
        }
        if chosen:
            rescued = {f"rescued_{col}": clean(chosen.get(col, "")) for col in COORD_COLUMNS}
            evidence.update(
                {
                    "rescue_source": chosen.get("_evidence_source", ""),
                    "rescue_clinvar_variation_id": chosen.get("clinvar_variation_id", chosen.get("clinvar_id", "")),
                    "rescue_clinvar_name": chosen.get("clinvar_name", ""),
                    "rescue_clinvar_classification": chosen.get("clinvar_classification", ""),
                    "rescue_clinvar_review_status": chosen.get("clinvar_review_status", ""),
                    "rescue_matched_terms": chosen.get("clinvar_matched_terms", ""),
                    "rescue_matched_term_types": chosen.get("clinvar_matched_term_types", ""),
                }
            )

        output_rows.append(
            {
                "variant_uid": uid,
                "study": row.get("study", ""),
                "source_row": row.get("source_row", ""),
                "participant_id": row.get("participant_id", ""),
                "gene": row.get("gene", ""),
                "variant_key": row.get("variant_key", ""),
                "gene_refgene": row.get("gene_refgene", ""),
                "hgvs_c": row.get("hgvs_c", ""),
                "hgvs_p": row.get("hgvs_p", ""),
                "hgvs_p_one_letter": row.get("hgvs_p_one_letter", ""),
                "inferred_consequence": row.get("inferred_consequence", ""),
                "coordinate_status": coord_status(row),
                "rescue_status": status,
                "rescue_confidence": confidence,
                "rescue_note": note,
                "query_genes": ";".join(sorted(queries[uid]["genes"])),
                "query_hgvs_terms": ";".join(sorted(queries[uid]["terms"])),
                "candidate_count": len(candidates),
                **original,
                **rescued,
                **evidence,
            }
        )

    columns = [
        "variant_uid",
        "study",
        "source_row",
        "participant_id",
        "gene",
        "variant_key",
        "gene_refgene",
        "hgvs_c",
        "hgvs_p",
        "hgvs_p_one_letter",
        "inferred_consequence",
        "coordinate_status",
        "rescue_status",
        "rescue_confidence",
        "rescue_note",
        "query_genes",
        "query_hgvs_terms",
        "candidate_count",
        "original_chrom",
        "original_pos",
        "original_ref",
        "original_alt",
        "rescued_chrom",
        "rescued_pos",
        "rescued_ref",
        "rescued_alt",
        "rescue_source",
        "rescue_clinvar_variation_id",
        "rescue_clinvar_name",
        "rescue_clinvar_classification",
        "rescue_clinvar_review_status",
        "rescue_matched_terms",
        "rescue_matched_term_types",
    ]
    write_tsv(output_path, output_rows, columns)

    summary = {
        "input": str(hiro_path),
        "clinvar": str(args.clinvar),
        "output": str(output_path),
        "total_hiro_rows": len(rows),
        "complete_coordinate_rows": sum(1 for row in rows if has_complete_coords(row)),
        "incomplete_coordinate_rows": len(target_rows),
        "coordinate_status_counts": dict(Counter(coord_status(row) for row in target_rows)),
        "rescue_status_counts": dict(status_counts),
        "rescue_confidence_counts": dict(confidence_counts),
        "rows_with_hgvs_query_terms": sum(1 for payload in queries.values() if payload["terms"]),
        "local_annotation_complete_candidate_rows": len(local_candidates),
        "cache_complete_candidate_rows": len(cache_candidates),
        "clinvar_candidate_rows": len(clinvar_candidates),
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("HiRO coordinate rescue candidate table written")
    print(f"- output: {output_path}")
    print(f"- incomplete coordinate rows: {len(target_rows)}")
    print(f"- rescue status counts: {dict(status_counts)}")
    print(f"- rescue confidence counts: {dict(confidence_counts)}")


if __name__ == "__main__":
    main()
