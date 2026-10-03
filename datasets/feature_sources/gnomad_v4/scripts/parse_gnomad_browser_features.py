#!/usr/bin/env python3
"""Parse gnomAD browser/Hail JSON rows into flat model-ready features."""

from __future__ import annotations

import csv
import gzip
import json
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[4]

INPUT = ROOT / "datasets/feature_sources/gnomad_v4/interim/gnomad_browser_v4_1_missing_combined_from_mixed_chunks.tsv.gz"
OUT_DIR = ROOT / "datasets/feature_sources/gnomad_v4/interim"
OUTPUT = OUT_DIR / "gnomad_browser_v4_1_parsed_features.tsv"
SUMMARY = OUT_DIR / "gnomad_browser_v4_1_parsed_features.summary.json"

DATASETS = ("joint", "exome", "genome")


BASE_FIELDS = [
    "variant_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "genes",
    "sources",
    "labels_3class",
    "gnomad_browser_status",
    "gnomad_browser_missing_reason",
    "gnomad_browser_variant_id",
    "gnomad_browser_caid",
    "gnomad_browser_rsids",
    "gnomad_browser_has_joint",
    "gnomad_browser_has_exome",
    "gnomad_browser_has_genome",
    "gnomad_browser_any_flags",
    "gnomad_browser_any_filters",
    "gnomad_browser_max_af",
    "gnomad_browser_max_af_source",
    "gnomad_browser_max_popmax_af",
    "gnomad_browser_max_popmax_source",
    "gnomad_browser_max_homozygote_count",
    "gnomad_browser_max_hemizygote_count",
    "gnomad_browser_canonical_gene_symbol",
    "gnomad_browser_canonical_transcript_id",
    "gnomad_browser_canonical_major_consequence",
    "gnomad_browser_canonical_hgvsc",
    "gnomad_browser_canonical_hgvsp",
    "gnomad_browser_canonical_lof",
    "gnomad_browser_mane_gene_symbol",
    "gnomad_browser_mane_transcript_id",
    "gnomad_browser_mane_major_consequence",
    "gnomad_browser_mane_hgvsc",
    "gnomad_browser_mane_hgvsp",
    "gnomad_browser_mane_lof",
    "gnomad_browser_consequence_terms",
]


def dataset_fields(dataset: str) -> list[str]:
    prefix = f"gnomad_browser_{dataset}"
    return [
        f"{prefix}_status",
        f"{prefix}_ac",
        f"{prefix}_ac_raw",
        f"{prefix}_an",
        f"{prefix}_af",
        f"{prefix}_homozygote_count",
        f"{prefix}_hemizygote_count",
        f"{prefix}_popmax_af",
        f"{prefix}_popmax_pop",
        f"{prefix}_popmax_ac",
        f"{prefix}_popmax_an",
        f"{prefix}_popmax_homozygote_count",
        f"{prefix}_popmax_hemizygote_count",
        f"{prefix}_grpmax_af",
        f"{prefix}_grpmax_pop",
        f"{prefix}_faf95_max",
        f"{prefix}_faf95_max_pop",
        f"{prefix}_faf99_max",
        f"{prefix}_faf99_max_pop",
        f"{prefix}_flags",
        f"{prefix}_filters",
        f"{prefix}_subsets",
    ]


FIELDS = BASE_FIELDS + [field for dataset in DATASETS for field in dataset_fields(dataset)]


def clean(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def num(value: Any) -> str:
    if value is None:
        return ""
    return str(value)


def af(ac: Any, an: Any) -> str:
    try:
        ac_f = float(ac)
        an_f = float(an)
    except (TypeError, ValueError):
        return ""
    if an_f <= 0:
        return ""
    return repr(ac_f / an_f)


def pipe(values: list[Any] | tuple[Any, ...] | None) -> str:
    if not values:
        return ""
    return "|".join(str(v) for v in values if v not in (None, ""))


def all_freq(dataset_obj: dict[str, Any] | None) -> dict[str, Any]:
    if not dataset_obj:
        return {}
    freq = dataset_obj.get("freq") or {}
    return freq.get("all") or {}


def is_population_group(group_id: str) -> bool:
    if not group_id:
        return False
    if group_id.endswith("_XX") or group_id.endswith("_XY"):
        return False
    if group_id in {"XX", "XY", "remaining"}:
        return False
    return True


def computed_popmax(freq_all: dict[str, Any]) -> dict[str, Any]:
    best: dict[str, Any] = {}
    best_af = -1.0
    for group in freq_all.get("ancestry_groups") or []:
        group_id = clean(group.get("id"))
        if not is_population_group(group_id):
            continue
        group_af = af(group.get("ac"), group.get("an"))
        if not group_af:
            continue
        group_af_f = float(group_af)
        if group_af_f > best_af:
            best_af = group_af_f
            best = {
                "af": group_af,
                "pop": group_id,
                "ac": num(group.get("ac")),
                "an": num(group.get("an")),
                "homozygote_count": num(group.get("homozygote_count")),
                "hemizygote_count": num(group.get("hemizygote_count")),
            }
    return best


def fafmax_values(dataset: str, dataset_obj: dict[str, Any] | None) -> dict[str, str]:
    if not dataset_obj:
        return {}
    raw = dataset_obj.get("fafmax")
    if not raw:
        return {}
    if dataset == "joint":
        src = raw
    elif isinstance(raw, dict):
        src = raw.get("gnomad") or raw.get("all") or raw.get("non_ukb") or raw.get("tgp") or raw.get("hgdp") or {}
    else:
        src = {}
    return {
        "faf95_max": num(src.get("faf95_max")),
        "faf95_max_pop": clean(src.get("faf95_max_gen_anc")),
        "faf99_max": num(src.get("faf99_max")),
        "faf99_max_pop": clean(src.get("faf99_max_gen_anc")),
    }


def dataset_row(dataset: str, dataset_obj: dict[str, Any] | None) -> dict[str, str]:
    prefix = f"gnomad_browser_{dataset}"
    out = {field: "" for field in dataset_fields(dataset)}
    out[f"{prefix}_status"] = "ok" if dataset_obj else "not_available"
    if not dataset_obj:
        return out

    freq_all = all_freq(dataset_obj)
    popmax = computed_popmax(freq_all)
    faf = fafmax_values(dataset, dataset_obj)
    grpmax = dataset_obj.get("grpmax") or {}

    out.update(
        {
            f"{prefix}_ac": num(freq_all.get("ac")),
            f"{prefix}_ac_raw": num(freq_all.get("ac_raw")),
            f"{prefix}_an": num(freq_all.get("an")),
            f"{prefix}_af": af(freq_all.get("ac"), freq_all.get("an")),
            f"{prefix}_homozygote_count": num(freq_all.get("homozygote_count")),
            f"{prefix}_hemizygote_count": num(freq_all.get("hemizygote_count")),
            f"{prefix}_popmax_af": popmax.get("af", ""),
            f"{prefix}_popmax_pop": popmax.get("pop", ""),
            f"{prefix}_popmax_ac": popmax.get("ac", ""),
            f"{prefix}_popmax_an": popmax.get("an", ""),
            f"{prefix}_popmax_homozygote_count": popmax.get("homozygote_count", ""),
            f"{prefix}_popmax_hemizygote_count": popmax.get("hemizygote_count", ""),
            f"{prefix}_grpmax_af": num(grpmax.get("AF")),
            f"{prefix}_grpmax_pop": clean(grpmax.get("gen_anc")),
            f"{prefix}_faf95_max": faf.get("faf95_max", ""),
            f"{prefix}_faf95_max_pop": faf.get("faf95_max_pop", ""),
            f"{prefix}_faf99_max": faf.get("faf99_max", ""),
            f"{prefix}_faf99_max_pop": faf.get("faf99_max_pop", ""),
            f"{prefix}_flags": pipe(dataset_obj.get("flags")),
            f"{prefix}_filters": pipe(dataset_obj.get("filters")),
            f"{prefix}_subsets": pipe(dataset_obj.get("subsets")),
        }
    )
    return out


def max_numeric(candidates: list[tuple[str, str]]) -> tuple[str, str]:
    best_value = ""
    best_source = ""
    best_float = -1.0
    for value, source in candidates:
        if value == "":
            continue
        try:
            value_float = float(value)
        except ValueError:
            continue
        if value_float > best_float:
            best_float = value_float
            best_value = value
            best_source = source
    return best_value, best_source


def selected_consequence(obj: dict[str, Any], selector: str) -> dict[str, str]:
    consequences = obj.get("transcript_consequences") or []
    if selector == "canonical":
        matches = [tc for tc in consequences if tc.get("is_canonical")]
    else:
        matches = [tc for tc in consequences if tc.get("is_mane_select")]
    if not matches and consequences:
        matches = [consequences[0]]
    tc = matches[0] if matches else {}
    return {
        "gene_symbol": clean(tc.get("gene_symbol")),
        "transcript_id": clean(tc.get("transcript_id")),
        "major_consequence": clean(tc.get("major_consequence")),
        "hgvsc": clean(tc.get("hgvsc")),
        "hgvsp": clean(tc.get("hgvsp")),
        "lof": clean(tc.get("lof")),
    }


def parse_row(row: dict[str, str], counts: Counter[str]) -> dict[str, str]:
    out = {field: "" for field in FIELDS}
    for field in ["variant_id", "chrom", "pos", "ref", "alt", "genes", "sources", "labels_3class"]:
        out[field] = row.get(field, "")
    out["gnomad_browser_status"] = row.get("gnomad_browser_status", "")
    out["gnomad_browser_missing_reason"] = row.get("gnomad_browser_missing_reason", "")

    if row.get("gnomad_browser_status") != "ok" or not row.get("gnomad_browser_row_json"):
        for dataset in DATASETS:
            out.update(dataset_row(dataset, None))
        return out

    try:
        obj = json.loads(row["gnomad_browser_row_json"])
    except json.JSONDecodeError as exc:
        out["gnomad_browser_status"] = "parse_error"
        out["gnomad_browser_missing_reason"] = f"invalid JSON: {exc}"
        counts["parse_error"] += 1
        for dataset in DATASETS:
            out.update(dataset_row(dataset, None))
        return out

    out["gnomad_browser_variant_id"] = clean(obj.get("variant_id"))
    out["gnomad_browser_caid"] = clean(obj.get("caid"))
    out["gnomad_browser_rsids"] = pipe(obj.get("rsids"))

    all_flags = []
    all_filters = []
    max_af_candidates = []
    max_popmax_candidates = []
    hom_candidates = []
    hemi_candidates = []
    for dataset in DATASETS:
        dataset_obj = obj.get(dataset)
        out[f"gnomad_browser_has_{dataset}"] = bool_text(bool(dataset_obj))
        ds_out = dataset_row(dataset, dataset_obj)
        out.update(ds_out)
        if dataset_obj:
            all_flags.extend(dataset_obj.get("flags") or [])
            all_filters.extend(dataset_obj.get("filters") or [])
            max_af_candidates.append((ds_out.get(f"gnomad_browser_{dataset}_af", ""), dataset))
            max_popmax_candidates.append((ds_out.get(f"gnomad_browser_{dataset}_popmax_af", ""), dataset))
            hom_candidates.append((ds_out.get(f"gnomad_browser_{dataset}_homozygote_count", ""), dataset))
            hemi_candidates.append((ds_out.get(f"gnomad_browser_{dataset}_hemizygote_count", ""), dataset))
    out["gnomad_browser_any_flags"] = pipe(sorted(set(all_flags)))
    out["gnomad_browser_any_filters"] = pipe(sorted(set(all_filters)))
    out["gnomad_browser_max_af"], out["gnomad_browser_max_af_source"] = max_numeric(max_af_candidates)
    out["gnomad_browser_max_popmax_af"], out["gnomad_browser_max_popmax_source"] = max_numeric(max_popmax_candidates)
    out["gnomad_browser_max_homozygote_count"], _ = max_numeric(hom_candidates)
    out["gnomad_browser_max_hemizygote_count"], _ = max_numeric(hemi_candidates)

    canonical = selected_consequence(obj, "canonical")
    mane = selected_consequence(obj, "mane")
    for key, value in canonical.items():
        out[f"gnomad_browser_canonical_{key}"] = value
    for key, value in mane.items():
        out[f"gnomad_browser_mane_{key}"] = value

    terms = set()
    for tc in obj.get("transcript_consequences") or []:
        if tc.get("major_consequence"):
            terms.add(clean(tc.get("major_consequence")))
        for term in tc.get("consequence_terms") or []:
            terms.add(clean(term))
    out["gnomad_browser_consequence_terms"] = pipe(sorted(terms))
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    counts: Counter[str] = Counter()
    status_counts: Counter[str] = Counter()
    with gzip.open(INPUT, "rt", newline="") as in_handle, OUTPUT.open("w", newline="") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        writer = csv.DictWriter(out_handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            counts["input_rows"] += 1
            parsed = parse_row(row, counts)
            status_counts[parsed.get("gnomad_browser_status", "") or "blank"] += 1
            for dataset in DATASETS:
                if parsed.get(f"gnomad_browser_has_{dataset}") == "true":
                    counts[f"has_{dataset}"] += 1
                if parsed.get(f"gnomad_browser_{dataset}_af"):
                    counts[f"nonmissing_{dataset}_af"] += 1
                if parsed.get(f"gnomad_browser_{dataset}_popmax_af"):
                    counts[f"nonmissing_{dataset}_popmax_af"] += 1
            if parsed.get("gnomad_browser_max_af"):
                counts["nonmissing_max_af"] += 1
            if parsed.get("gnomad_browser_any_filters"):
                counts["has_any_filter"] += 1
            if parsed.get("gnomad_browser_any_flags"):
                counts["has_any_flag"] += 1
            writer.writerow(parsed)

    summary = {
        "input": str(INPUT.relative_to(ROOT)),
        "output": str(OUTPUT.relative_to(ROOT)),
        "output_rows": counts["input_rows"],
        "output_columns": len(FIELDS),
        "status_counts": dict(sorted(status_counts.items())),
        "counts": dict(sorted(counts.items())),
        "notes": [
            "This parser flattens gnomAD browser JSON into population-frequency, popmax, filter/flag, and compact transcript-consequence features.",
            "Popmax AF is computed from non-sex-stratified ancestry groups in each dataset-specific all-frequency block.",
            "The original combined browser output remains unchanged and retains full JSON rows.",
        ],
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
