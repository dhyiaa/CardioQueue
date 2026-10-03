#!/usr/bin/env python3
"""Apply current gene-panel QC decisions to curated ClinVar rows."""

from __future__ import annotations

import csv
import json
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
CLINVAR_CURATED = ROOT / "datasets/clinvar/interim/clinvar_cardiogenetics_modeling_slice.tsv"
SOURCE_COUNTS = ROOT / "datasets/gene_panels/interim/source_gene_counts.tsv"
OUT_DIR = ROOT / "datasets/gene_panels/interim"

DECISIONS_OUT = OUT_DIR / "gene_panel_decisions.tsv"
PRIMARY_PANEL_OUT = ROOT / "datasets/gene_panels/cardiogenetics_classifier_genes.keep.txt"
REMOVE_PANEL_OUT = ROOT / "datasets/gene_panels/cardiogenetics_classifier_genes.remove.txt"
SENSITIVITY_PANEL_OUT = ROOT / "datasets/gene_panels/cardiogenetics_classifier_genes.sensitivity_or_separate.txt"

PRIMARY_CLINVAR_OUT = ROOT / "datasets/clinvar/interim/clinvar_primary_after_gene_qc.tsv"
HOLDOUT_CLINVAR_OUT = ROOT / "datasets/clinvar/interim/clinvar_gene_qc_holdout_or_sensitivity.tsv"
REMOVED_CLINVAR_OUT = ROOT / "datasets/clinvar/interim/clinvar_removed_gene_contamination.tsv"
TTN_CLINVAR_OUT = ROOT / "datasets/clinvar/interim/clinvar_ttn_separate_review.tsv"
SUMMARY_OUT = OUT_DIR / "gene_panel_decisions.summary.json"


PRIMARY_EXCLUDE = {
    "AKAP9": (
        "remove_contamination",
        "ClinGen LQTS disputed and only 3 pathogenic rows in curated ClinVar; exclude from primary model.",
    ),
    "DMD": (
        "remove_contamination",
        "Dystrophinopathy/neuromuscular primary label context; cardiac involvement is phenotype-specific.",
    ),
    "FPGT": (
        "remove_contamination",
        "Non-cardiac FPGT rows entered through nearby/readthrough FPGT-TNNI3K/TNNI3K locus.",
    ),
}

ALIAS_OR_QUARANTINE = {
    "FPGT-TNNI3K": (
        "alias_or_quarantine",
        "HGNC readthrough locus; map transcript-aware to TNNI3K when supported, otherwise quarantine.",
    ),
    "KCNE1B": (
        "alias_or_quarantine",
        "Withdrawn/false-duplication symbol; collapse evidence-supported rows to KCNE1, otherwise quarantine.",
    ),
    "KNCH2": (
        "alias_or_quarantine",
        "Likely KCNH2 typo; correct only when variant/transcript evidence supports it.",
    ),
}

SEPARATE_OR_SENSITIVITY = {
    "TTN": (
        "separate_ttn",
        "Do not include in primary missense-like classifier; analyze TTN separately, especially truncating variants for DCM.",
    ),
    "ANK2": (
        "sensitivity_only",
        "ClinGen disputes ANK2 for LQTS/Brugada/CPVT; large VUS burden and low pathogenic yield.",
    ),
    "SOS1": (
        "sensitivity_only",
        "RASopathy/Noonan gene; HCM can occur syndromically but this is not primary cardiac-specific context.",
    ),
    "TNNI3K": (
        "sensitivity_or_pool",
        "Real cardiac gene with ClinGen moderate DCM evidence but very sparse pathogenic counts; pool with DCM/sensitivity.",
    ),
    "TRPM4": (
        "sensitivity_only",
        "ClinGen refutes TRPM4-Brugada; use only in a documented conduction-disease sensitivity scope.",
    ),
}

CONDITIONAL_PHENOCOPY = {
    "EMD": "Muscular dystrophy/cardiac conduction overlap; keep only if expanded cardiomyopathy-overlap scope is intended.",
    "GLA": "Fabry cardiac phenocopy; CardioBoost included it, retain with phenocopy flag.",
    "HRAS": "RASopathy gene; keep only in expanded phenocopy/syndromic sensitivity analyses.",
    "KRAS": "RASopathy gene; keep only in expanded phenocopy/syndromic sensitivity analyses.",
    "LAMP2": "Danon cardiac phenocopy; CardioBoost included it, retain with phenocopy flag.",
    "MAP2K2": "RASopathy gene; keep only in expanded phenocopy/syndromic sensitivity analyses.",
    "NRAS": "RASopathy gene; keep only in expanded phenocopy/syndromic sensitivity analyses.",
    "PTPN11": "Noonan/RASopathy HCM-overlap gene; CardioBoost included it, retain with phenocopy flag.",
    "RAF1": "RASopathy gene with HCM overlap; keep only in expanded phenocopy/syndromic sensitivity analyses.",
    "RIT1": "RASopathy gene; keep only in expanded phenocopy/syndromic sensitivity analyses.",
    "SLC22A5": "Primary carnitine deficiency/metabolic cardiomyopathy phenocopy; retain with phenocopy flag if scope includes phenocopies.",
}

MONITOR_EVIDENCE = {
    "CACNB2": "Channelopathy panel gene with low pathogenic yield; retain but monitor.",
    "FHOD3": "Cardiomyopathy gene with limited/smaller evidence base; retain but monitor.",
    "MIB1": "Cardiomyopathy/LVNC-associated panel gene; retain but monitor.",
    "SCN1B": "Channelopathy panel gene with modest counts; retain but monitor.",
    "SNTA1": "Channelopathy panel gene with low pathogenic yield; retain but monitor.",
}

TRUNCATING_RE = re.compile(
    r"frameshift|\bfs\b|\bfs\*|\*|ter|stop|nonsense|splice|donor|acceptor",
    re.IGNORECASE,
)


def load_source_genes() -> set[str]:
    genes = set()
    with SOURCE_COUNTS.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            genes.add(row["gene"])
    return genes


def decision_for_gene(gene: str) -> tuple[str, str]:
    if gene in PRIMARY_EXCLUDE:
        return PRIMARY_EXCLUDE[gene]
    if gene in ALIAS_OR_QUARANTINE:
        return ALIAS_OR_QUARANTINE[gene]
    if gene in SEPARATE_OR_SENSITIVITY:
        return SEPARATE_OR_SENSITIVITY[gene]
    if gene in CONDITIONAL_PHENOCOPY:
        return "conditional_phenocopy", CONDITIONAL_PHENOCOPY[gene]
    if gene in MONITOR_EVIDENCE:
        return "primary_keep_monitor", MONITOR_EVIDENCE[gene]
    return "primary_keep", "Retain in current primary cardiogenetics feature set."


def write_decisions(genes: set[str]) -> dict[str, tuple[str, str]]:
    decisions = {gene: decision_for_gene(gene) for gene in sorted(genes)}
    fields = ["gene", "decision", "primary_model_inclusion", "reason"]
    with DECISIONS_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for gene, (decision, reason) in decisions.items():
            primary = "include" if decision in {"primary_keep", "primary_keep_monitor", "conditional_phenocopy"} else "exclude"
            writer.writerow(
                {
                    "gene": gene,
                    "decision": decision,
                    "primary_model_inclusion": primary,
                    "reason": reason,
                }
            )
    return decisions


def write_panels(decisions: dict[str, tuple[str, str]]) -> None:
    primary = sorted(
        gene
        for gene, (decision, _reason) in decisions.items()
        if decision in {"primary_keep", "primary_keep_monitor", "conditional_phenocopy"}
    )
    remove = sorted(
        gene
        for gene, (decision, _reason) in decisions.items()
        if decision in {"remove_contamination", "alias_or_quarantine"}
    )
    sensitivity = sorted(
        gene
        for gene, (decision, _reason) in decisions.items()
        if decision in {"separate_ttn", "sensitivity_only", "sensitivity_or_pool"}
    )
    PRIMARY_PANEL_OUT.write_text("\n".join(primary) + "\n")
    REMOVE_PANEL_OUT.write_text("\n".join(remove) + "\n")
    SENSITIVITY_PANEL_OUT.write_text("\n".join(sensitivity) + "\n")


def is_likely_ttn_truncating(row: dict[str, str]) -> bool:
    text = " ".join([row.get("name", ""), row.get("variant_type", ""), row.get("clinical_significance", "")])
    return bool(TRUNCATING_RE.search(text))


def write_clinvar_slices(decisions: dict[str, tuple[str, str]]) -> dict[str, object]:
    with CLINVAR_CURATED.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        extra_fields = ["gene_panel_decision", "gene_panel_reason", "ttn_likely_truncating"]
        out_fields = fields + extra_fields
        buckets = {
            "primary": [],
            "holdout": [],
            "removed": [],
            "ttn": [],
        }
        for row in reader:
            gene = row["gene"]
            decision, reason = decisions.get(gene, decision_for_gene(gene))
            out = dict(row)
            out["gene_panel_decision"] = decision
            out["gene_panel_reason"] = reason
            out["ttn_likely_truncating"] = "true" if gene == "TTN" and is_likely_ttn_truncating(row) else "false"
            if decision in {"primary_keep", "primary_keep_monitor", "conditional_phenocopy"}:
                buckets["primary"].append(out)
            elif decision == "remove_contamination":
                buckets["removed"].append(out)
                buckets["holdout"].append(out)
            elif decision == "separate_ttn":
                buckets["ttn"].append(out)
                buckets["holdout"].append(out)
            else:
                buckets["holdout"].append(out)

    for path, rows in [
        (PRIMARY_CLINVAR_OUT, buckets["primary"]),
        (HOLDOUT_CLINVAR_OUT, buckets["holdout"]),
        (REMOVED_CLINVAR_OUT, buckets["removed"]),
        (TTN_CLINVAR_OUT, buckets["ttn"]),
    ]:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=out_fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def summarize(rows: list[dict[str, str]]) -> dict[str, object]:
        return {
            "rows": len(rows),
            "label_counts": dict(sorted(Counter(row["label_3class"] for row in rows).items())),
            "gene_counts": dict(Counter(row["gene"] for row in rows).most_common()),
            "decision_counts": dict(sorted(Counter(row["gene_panel_decision"] for row in rows).items())),
        }

    return {name: summarize(rows) for name, rows in buckets.items()}


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    genes = load_source_genes()
    decisions = write_decisions(genes)
    write_panels(decisions)
    slice_summary = write_clinvar_slices(decisions)
    summary = {
        "inputs": {
            "source_counts": str(SOURCE_COUNTS.relative_to(ROOT)),
            "clinvar_curated": str(CLINVAR_CURATED.relative_to(ROOT)),
        },
        "outputs": {
            "decisions": str(DECISIONS_OUT.relative_to(ROOT)),
            "primary_panel": str(PRIMARY_PANEL_OUT.relative_to(ROOT)),
            "remove_or_quarantine_panel": str(REMOVE_PANEL_OUT.relative_to(ROOT)),
            "sensitivity_or_separate_panel": str(SENSITIVITY_PANEL_OUT.relative_to(ROOT)),
            "primary_clinvar": str(PRIMARY_CLINVAR_OUT.relative_to(ROOT)),
            "holdout_clinvar": str(HOLDOUT_CLINVAR_OUT.relative_to(ROOT)),
            "removed_clinvar": str(REMOVED_CLINVAR_OUT.relative_to(ROOT)),
            "ttn_clinvar": str(TTN_CLINVAR_OUT.relative_to(ROOT)),
        },
        "decision_counts": dict(sorted(Counter(decision for decision, _reason in decisions.values()).items())),
        "clinvar_slices": slice_summary,
        "notes": [
            "No raw source rows are deleted.",
            "Primary ClinVar slice excludes requested Tier B/C genes from primary training.",
            "TTN is separated for consequence-aware review, not discarded.",
            "Alias/quarantine symbols are absent from the curated ClinVar slice but present in source gene counts.",
        ],
    }
    SUMMARY_OUT.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
