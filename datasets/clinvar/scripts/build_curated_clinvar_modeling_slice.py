#!/usr/bin/env python3
"""Build a curated ClinVar cardiogenetics modeling slice.

This is a label-cleaning layer, not the final train/test split. It keeps
GRCh38 coordinate variants in the current registry gene panel, restricts to
germline ClinVar assertions, removes somatic/cancer-context rows, and separates
high-confidence rows from lower-confidence sensitivity rows.
"""

from __future__ import annotations

import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
CLINVAR = ROOT / "datasets/clinvar/variant_summary.txt.gz"
GENE_PANEL = ROOT / "datasets/gene_panels/cardiogenetics_classifier_genes.keep.txt"
OUT_DIR = ROOT / "datasets/clinvar/interim"
COMBINED_OUT = OUT_DIR / "clinvar_cardiogenetics_modeling_slice.tsv"
HIGH_CONF_OUT = OUT_DIR / "clinvar_cardiogenetics_high_confidence.tsv"
LOWER_CONF_OUT = OUT_DIR / "clinvar_cardiogenetics_lower_confidence_star1.tsv"
SUMMARY_OUT = OUT_DIR / "clinvar_cardiogenetics_modeling_slice.summary.json"

MISSING = {"", "-", ".", "na", "nan", "none", "not provided"}
CANCER_RE = re.compile(
    r"cancer|carcinoma|tumou?r|neoplasm|leukemia|lymphoma|melanoma|oncolog|glioma|blastoma|sarcoma",
    re.IGNORECASE,
)

OUTPUT_FIELDS = [
    "variant_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "gene",
    "label_3class",
    "clinical_significance",
    "confidence_tier",
    "review_stars",
    "review_status",
    "number_submitters",
    "last_evaluated",
    "origin",
    "origin_simple",
    "variation_id",
    "allele_id",
    "rcv_accession",
    "name",
    "variant_type",
    "phenotype_ids",
    "phenotype_list",
    "guidelines",
    "tested_in_gtr",
    "source",
]


def clean(value: object) -> str:
    return str(value if value is not None else "").strip()


def is_present(value: object) -> bool:
    return clean(value).lower() not in MISSING


def norm_chrom(chrom: object) -> str:
    value = clean(chrom)
    if value.lower().startswith("chr"):
        value = value[3:]
    return value.upper() if value.upper() in {"X", "Y", "MT", "M"} else value


def norm_allele(allele: object) -> str:
    return clean(allele).upper()


def variant_id(row: dict[str, str]) -> str | None:
    chrom = norm_chrom(row.get("Chromosome", ""))
    pos = clean(row.get("PositionVCF", ""))
    ref = norm_allele(row.get("ReferenceAlleleVCF", ""))
    alt = norm_allele(row.get("AlternateAlleleVCF", ""))
    if not all(is_present(x) for x in [chrom, pos, ref, alt]):
        return None
    try:
        if int(pos) <= 0:
            return None
    except ValueError:
        return None
    return f"{chrom}-{pos}-{ref}-{alt}"


def split_genes(raw: object) -> list[str]:
    text = clean(raw).replace(";", ",").replace("|", ",")
    return [piece.strip().upper() for piece in text.split(",") if piece.strip()]


def label_3class(raw: object) -> str:
    text = clean(raw).lower()
    if not text:
        return ""
    if "conflict" in text or "drug response" in text or "association" in text:
        return "Conflict_or_other"
    has_path = "pathogenic" in text
    has_benign = "benign" in text
    has_vus = "uncertain" in text or "vus" in text
    if has_path and has_benign:
        return "Conflict_or_other"
    if has_path:
        return "Pathogenic"
    if has_benign:
        return "Benign"
    if has_vus:
        return "VUS"
    return "Other"


def clinvar_stars(review_status: object) -> str:
    text = clean(review_status).lower()
    if "practice guideline" in text:
        return "4"
    if "reviewed by expert panel" in text:
        return "3"
    if "multiple submitters" in text and "no conflicts" in text:
        return "2"
    if "criteria provided" in text:
        return "1"
    if "no assertion" in text or "no classification" in text:
        return "0"
    return ""


def has_somatic_or_oncogenic_context(row: dict[str, str]) -> bool:
    somatic_fields = [
        "SomaticClinicalImpact",
        "SomaticClinicalImpactLastEvaluated",
        "ReviewStatusClinicalImpact",
        "Oncogenicity",
        "OncogenicityLastEvaluated",
        "ReviewStatusOncogenicity",
        "SCVsForAggregateSomaticClinicalImpact",
        "SCVsForAggregateOncogenicityClassification",
    ]
    if "somatic" in clean(row.get("Origin", "")).lower():
        return True
    return any(is_present(row.get(field, "")) for field in somatic_fields)


def has_cancer_context(row: dict[str, str]) -> bool:
    text = " | ".join(
        [
            clean(row.get("PhenotypeList", "")),
            clean(row.get("PhenotypeIDS", "")),
            clean(row.get("Name", "")),
        ]
    )
    return bool(CANCER_RE.search(text))


def load_gene_panel() -> set[str]:
    with GENE_PANEL.open() as handle:
        return {line.strip().upper() for line in handle if line.strip()}


def row_to_output(row: dict[str, str], gene: str, label: str, stars: str, tier: str, vid: str) -> dict[str, str]:
    chrom, pos, ref, alt = vid.split("-", 3)
    return {
        "variant_id": vid,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "gene": gene,
        "label_3class": label,
        "clinical_significance": clean(row.get("ClinicalSignificance", "")),
        "confidence_tier": tier,
        "review_stars": stars,
        "review_status": clean(row.get("ReviewStatus", "")),
        "number_submitters": clean(row.get("NumberSubmitters", "")),
        "last_evaluated": clean(row.get("LastEvaluated", "")),
        "origin": clean(row.get("Origin", "")),
        "origin_simple": clean(row.get("OriginSimple", "")),
        "variation_id": clean(row.get("VariationID", "")),
        "allele_id": clean(row.get("#AlleleID", "")),
        "rcv_accession": clean(row.get("RCVaccession", "")),
        "name": clean(row.get("Name", "")),
        "variant_type": clean(row.get("Type", "")),
        "phenotype_ids": clean(row.get("PhenotypeIDS", "")),
        "phenotype_list": clean(row.get("PhenotypeList", "")),
        "guidelines": clean(row.get("Guidelines", "")),
        "tested_in_gtr": clean(row.get("TestedInGTR", "")),
        "source": "ClinVar",
    }


def choose_best(rows: list[dict[str, str]]) -> dict[str, str]:
    label_priority = {"Pathogenic": 3, "Benign": 2, "VUS": 1}

    def key(row: dict[str, str]) -> tuple[int, int, int]:
        try:
            stars = int(row["review_stars"])
        except ValueError:
            stars = -1
        try:
            submitters = int(row["number_submitters"])
        except ValueError:
            submitters = 0
        return (stars, submitters, label_priority.get(row["label_3class"], 0))

    return sorted(rows, key=key, reverse=True)[0]


def main() -> int:
    panel = load_gene_panel()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary: Counter[str] = Counter()
    excluded: Counter[str] = Counter()
    label_counts: Counter[str] = Counter()
    tier_counts: Counter[str] = Counter()
    star_counts: Counter[str] = Counter()
    gene_counts: Counter[str] = Counter()
    by_variant: dict[str, list[dict[str, str]]] = defaultdict(list)

    with gzip.open(CLINVAR, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            summary["clinvar_rows_seen"] += 1
            if row.get("Assembly") != "GRCh38":
                excluded["not_grch38"] += 1
                continue
            genes = [gene for gene in split_genes(row.get("GeneSymbol", "")) if gene in panel]
            if not genes:
                excluded["outside_current_gene_panel"] += 1
                continue
            vid = variant_id(row)
            if not vid:
                excluded["missing_coordinate_key"] += 1
                continue
            if clean(row.get("OriginSimple", "")).lower() != "germline":
                excluded["origin_simple_not_germline"] += 1
                continue
            if has_somatic_or_oncogenic_context(row):
                excluded["somatic_or_oncogenic_context"] += 1
                continue
            if has_cancer_context(row):
                excluded["cancer_context"] += 1
                continue
            label = label_3class(row.get("ClinicalSignificance", ""))
            if label in {"Conflict_or_other", "Other", ""}:
                excluded[f"label_{label or 'blank'}"] += 1
                continue
            stars = clinvar_stars(row.get("ReviewStatus", ""))
            if stars in {"2", "3", "4"}:
                tier = "high_confidence_stars_ge2"
            elif stars == "1":
                tier = "lower_confidence_star1"
            else:
                excluded[f"review_stars_{stars or 'blank'}"] += 1
                continue
            for gene in genes:
                out = row_to_output(row, gene, label, stars, tier, vid)
                by_variant[vid].append(out)
                label_counts[label] += 1
                tier_counts[tier] += 1
                star_counts[stars] += 1
                gene_counts[gene] += 1

    selected = [choose_best(rows) for rows in by_variant.values()]
    selected.sort(key=lambda row: (row["chrom"], int(row["pos"]) if row["pos"].isdigit() else 0, row["ref"], row["alt"]))

    high_rows = [row for row in selected if row["confidence_tier"] == "high_confidence_stars_ge2"]
    lower_rows = [row for row in selected if row["confidence_tier"] == "lower_confidence_star1"]

    for path, rows in [(COMBINED_OUT, selected), (HIGH_CONF_OUT, high_rows), (LOWER_CONF_OUT, lower_rows)]:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    selected_label_counts = Counter(row["label_3class"] for row in selected)
    selected_tier_counts = Counter(row["confidence_tier"] for row in selected)
    selected_star_counts = Counter(row["review_stars"] for row in selected)
    selected_gene_counts = Counter(row["gene"] for row in selected)
    summary_payload = {
        "inputs": {
            "clinvar": str(CLINVAR.relative_to(ROOT)),
            "gene_panel": str(GENE_PANEL.relative_to(ROOT)),
            "gene_panel_size": len(panel),
        },
        "outputs": {
            "combined": str(COMBINED_OUT.relative_to(ROOT)),
            "high_confidence": str(HIGH_CONF_OUT.relative_to(ROOT)),
            "lower_confidence_star1": str(LOWER_CONF_OUT.relative_to(ROOT)),
        },
        "filters": [
            "Assembly == GRCh38",
            "GeneSymbol intersects audited cardiogenetics classifier keep panel",
            "complete chrom-pos-ref-alt key from ClinVar VCF columns",
            "OriginSimple == germline",
            "no somatic/oncogenic ClinVar fields",
            "no cancer-context phenotype/name keyword",
            "label in Benign, Pathogenic, VUS",
            "review stars >= 1",
        ],
        "source_rows_passing_before_variant_dedup": sum(tier_counts.values()),
        "deduplicated_variant_rows": len(selected),
        "high_confidence_rows": len(high_rows),
        "lower_confidence_star1_rows": len(lower_rows),
        "selected_label_counts": dict(sorted(selected_label_counts.items())),
        "selected_tier_counts": dict(sorted(selected_tier_counts.items())),
        "selected_star_counts": dict(sorted(selected_star_counts.items())),
        "selected_top_genes": dict(selected_gene_counts.most_common(25)),
        "source_row_label_counts_before_dedup": dict(sorted(label_counts.items())),
        "source_row_tier_counts_before_dedup": dict(sorted(tier_counts.items())),
        "source_row_star_counts_before_dedup": dict(sorted(star_counts.items())),
        "excluded_counts": dict(sorted(excluded.items())),
        "notes": [
            "This is a curated ClinVar label pool, not a train/test split.",
            "VUS rows are retained as a clinically meaningful third class but should be treated as noisier than B/P labels.",
            "Current gene panel is inherited from the master registry and still needs final cardiogenetics-panel freeze.",
        ],
    }
    SUMMARY_OUT.write_text(json.dumps(summary_payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary_payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
