#!/usr/bin/env python3
"""Rescue eMERGE source hg19 coordinates to GRCh38 with explicit QC.

The current modeling matrix treated eMERGE source coordinates as GRCh38. The
eMERGE provenance code shows those coordinates are hg19/GRCh37-like source
coordinates, so this script creates a side-by-side rescue table:

- original/source hg19 coordinate and variant id
- lifted GRCh38 coordinate and variant id
- liftover status and allele/ref validation against the local GRCh38 FASTA
- whether GRCh38 VEP consequence later agrees with eMERGE's source consequence

This script intentionally does not overwrite the existing eMERGE tables.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pandas as pd
import pysam
from pyliftover import LiftOver


ROOT = Path(__file__).resolve().parents[4]
EMERGE = ROOT / "datasets/emerge/full_arrhythmia_gene_dataset"
SOURCE = EMERGE / "data/emerge_arrhythmia_ml_enriched_features.csv"
LINKED = EMERGE / "data/emerge_linked_source_records.tsv"
CHAIN = EMERGE / "provenance/liftover/hg19ToHg38.over.chain.gz"
FASTA = ROOT / "datasets/feature_sources/insilico/raw/vep/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

OUT_DIR = EMERGE / "interim"
RESCUE_OUT = OUT_DIR / "emerge_coordinate_rescue_liftover.tsv"
LINKED_OUT = OUT_DIR / "emerge_linked_source_records_grch38_rescued.tsv"
SUMMARY_OUT = OUT_DIR / "emerge_coordinate_rescue_liftover.summary.json"

MISSING = {"", ".", "NA", "N/A", "nan", "None", "Missing", "missing"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def clean(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text in MISSING:
        return ""
    if text.endswith(".0"):
        try:
            return str(int(float(text)))
        except ValueError:
            pass
    return text


def norm_chrom(value: object) -> str:
    text = clean(value)
    if text.lower().startswith("chr"):
        text = text[3:]
    return text


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1].upper()


def variant_id(chrom: str, pos: str | int, ref: str, alt: str) -> str:
    return f"{norm_chrom(chrom)}-{clean(pos)}-{clean(ref).upper()}-{clean(alt).upper()}"


def fetch_ref(fasta: pysam.FastaFile, chrom: str, pos_1based: int, ref_len: int) -> str:
    names = set(fasta.references)
    candidates = [chrom, f"chr{chrom}"]
    contig = next((c for c in candidates if c in names), "")
    if not contig or pos_1based < 1 or ref_len < 1:
        return ""
    try:
        return fasta.fetch(contig, pos_1based - 1, pos_1based - 1 + ref_len).upper()
    except Exception:
        return ""


def classify_consequence_agreement(source_consequence: str, vep_consequence: str) -> str:
    src = clean(source_consequence).lower()
    vep = clean(vep_consequence).lower()
    if not src or not vep:
        return "not_evaluated"
    if src == "missense":
        return "agree" if "missense_variant" in vep else "disagree"
    if src == "frameshift":
        return "agree" if "frameshift_variant" in vep else "disagree"
    if src == "nonsense":
        return "agree" if "stop_gained" in vep else "disagree"
    if src == "splice_region":
        return "agree" if "splice" in vep else "disagree"
    if src in {"synonymous", "silent"}:
        return "agree" if "synonymous_variant" in vep else "disagree"
    return "not_evaluated"


def liftover_one(lo: LiftOver, fasta: pysam.FastaFile, row: pd.Series) -> dict[str, str]:
    chrom = norm_chrom(row.get("chrom"))
    pos = clean(row.get("pos"))
    ref = clean(row.get("ref")).upper()
    alt = clean(row.get("alt")).upper()
    source_vid = variant_id(chrom, pos, ref, alt) if chrom and pos and ref and alt else ""

    out = {
        "source_hg19_chrom": chrom,
        "source_hg19_pos": pos,
        "source_hg19_ref": ref,
        "source_hg19_alt": alt,
        "source_hg19_variant_id": source_vid,
        "grch38_chrom": "",
        "grch38_pos": "",
        "grch38_ref": "",
        "grch38_alt": "",
        "grch38_variant_id": "",
        "liftover_status": "",
        "liftover_chain_strand": "",
        "grch38_ref_at_pos": "",
        "allele_ref_validation": "",
        "coordinate_qc_flag": "",
    }

    if not (chrom and pos and ref and alt):
        out.update(
            {
                "liftover_status": "missing_source_coordinate",
                "allele_ref_validation": "not_evaluated",
                "coordinate_qc_flag": "fail_missing_source_coordinate",
            }
        )
        return out

    try:
        pos_int = int(float(pos))
    except ValueError:
        out.update(
            {
                "liftover_status": "invalid_source_position",
                "allele_ref_validation": "not_evaluated",
                "coordinate_qc_flag": "fail_invalid_source_position",
            }
        )
        return out

    hits = lo.convert_coordinate(f"chr{chrom}", pos_int - 1)
    if not hits:
        hits = lo.convert_coordinate(chrom, pos_int - 1)
    if not hits:
        out.update(
            {
                "liftover_status": "no_mapping",
                "allele_ref_validation": "not_evaluated",
                "coordinate_qc_flag": "fail_no_mapping",
            }
        )
        return out

    # Prefer canonical autosome/sex-chrom mappings. If more than one remains,
    # keep the first but flag it.
    primary = [
        h for h in hits if norm_chrom(h[0]) in {str(i) for i in range(1, 23)} | {"X", "Y", "M", "MT"}
    ]
    chosen = primary[0] if primary else hits[0]
    new_chrom = norm_chrom(chosen[0])
    new_pos = int(chosen[1]) + 1
    strand = chosen[2]
    new_ref = ref
    new_alt = alt
    if strand == "-":
        new_ref = rc(ref)
        new_alt = rc(alt)

    ref_at_pos = fetch_ref(fasta, new_chrom, new_pos, len(new_ref))
    if ref_at_pos == new_ref:
        validation = "match"
        qc = "pass"
    elif ref_at_pos:
        validation = "mismatch"
        qc = "fail_ref_mismatch"
    else:
        validation = "not_evaluated"
        qc = "fail_reference_fetch"

    status = "lifted_single" if len(hits) == 1 else "lifted_multiple"
    out.update(
        {
            "grch38_chrom": new_chrom,
            "grch38_pos": str(new_pos),
            "grch38_ref": new_ref,
            "grch38_alt": new_alt,
            "grch38_variant_id": variant_id(new_chrom, new_pos, new_ref, new_alt),
            "liftover_status": status,
            "liftover_chain_strand": strand,
            "grch38_ref_at_pos": ref_at_pos,
            "allele_ref_validation": validation,
            "coordinate_qc_flag": qc,
        }
    )
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not CHAIN.exists():
        raise SystemExit(f"Missing liftover chain: {CHAIN}")
    if not FASTA.exists():
        raise SystemExit(f"Missing GRCh38 FASTA: {FASTA}")

    lo = LiftOver(str(CHAIN))
    fasta = pysam.FastaFile(str(FASTA))
    source = pd.read_csv(SOURCE, dtype=str, keep_default_na=False)

    records = []
    for _, row in source.iterrows():
        rescue = liftover_one(lo, fasta, row)
        source_vid = rescue["source_hg19_variant_id"]
        record = {
            "variant_uid": clean(row.get("variant_uid")),
            "variant_key": clean(row.get("variant_key")),
            "study": clean(row.get("study")),
            "source_row": clean(row.get("source_row")),
            "gene": clean(row.get("gene")),
            "target_3class": clean(row.get("target_3class")),
            "target_5class": clean(row.get("target_5class")),
            "target_classification_raw": clean(row.get("target_classification_raw")),
            "hgvs_p": clean(row.get("hgvs_p")),
            "hgvs_p_one_letter": clean(row.get("hgvs_p_one_letter")),
            "source_inferred_consequence": clean(row.get("inferred_consequence")),
            "source_cadd_phred": clean(row.get("cadd_phred")),
            "source_gnomad_max_af": clean(row.get("gnomad_max_af")),
            **rescue,
        }
        # The VEP fields are filled by a later merge step. The blank columns make
        # the intended audit schema explicit.
        record["vep_consequence_after_rescue"] = ""
        record["vep_impact_after_rescue"] = ""
        record["source_vs_vep_consequence_agreement"] = "not_evaluated"
        records.append(record)

    rescue_df = pd.DataFrame(records)
    rescue_df.to_csv(RESCUE_OUT, sep="\t", index=False)

    linked = pd.read_csv(LINKED, sep="\t", dtype=str, keep_default_na=False)
    linked = linked.merge(
        rescue_df[
            [
                "variant_uid",
                "source_hg19_variant_id",
                "grch38_variant_id",
                "grch38_chrom",
                "grch38_pos",
                "grch38_ref",
                "grch38_alt",
                "liftover_status",
                "liftover_chain_strand",
                "grch38_ref_at_pos",
                "allele_ref_validation",
                "coordinate_qc_flag",
            ]
        ],
        on="variant_uid",
        how="left",
    )
    linked["pre_rescue_resolved_variant_id"] = linked["resolved_variant_id"]
    linked["resolved_variant_id"] = linked["grch38_variant_id"].where(
        linked["coordinate_qc_flag"].eq("pass"), linked["resolved_variant_id"]
    )
    linked["source_registry_link_status"] = linked["coordinate_qc_flag"].map(
        lambda x: "linked_to_registry_grch38_liftover" if x == "pass" else "liftover_qc_failed_keep_original"
    )
    for src, dest in [
        ("grch38_chrom", "chrom"),
        ("grch38_pos", "pos"),
        ("grch38_ref", "ref"),
        ("grch38_alt", "alt"),
    ]:
        linked[dest] = linked[src].where(linked["coordinate_qc_flag"].eq("pass"), linked[dest])
    linked.to_csv(LINKED_OUT, sep="\t", index=False)

    counts = Counter()
    for col in ["target_3class", "liftover_status", "allele_ref_validation", "coordinate_qc_flag"]:
        counts.update({f"{col}:{k}": int(v) for k, v in rescue_df[col].value_counts(dropna=False).items()})
    summary = {
        "input_source": str(SOURCE.relative_to(ROOT)),
        "input_linked_source_records": str(LINKED.relative_to(ROOT)),
        "chain": str(CHAIN.relative_to(ROOT)),
        "fasta": str(FASTA.relative_to(ROOT)),
        "rescue_table": str(RESCUE_OUT.relative_to(ROOT)),
        "rescued_linked_source_records": str(LINKED_OUT.relative_to(ROOT)),
        "rows": int(len(rescue_df)),
        "unique_source_hg19_variant_ids": int(rescue_df["source_hg19_variant_id"].nunique()),
        "unique_grch38_variant_ids_pass": int(
            rescue_df.loc[rescue_df["coordinate_qc_flag"].eq("pass"), "grch38_variant_id"].nunique()
        ),
        "counts": dict(sorted(counts.items())),
        "notes": [
            "Coordinates are lifted from hg19 source coordinates to GRCh38 with UCSC hg19ToHg38 chain.",
            "coordinate_qc_flag == pass requires the lifted GRCh38 reference allele to match the local GRCh38 FASTA.",
            "Indels are retained only if their lifted POS+REF validates against GRCh38; otherwise they remain flagged for manual/remapped handling.",
            "Existing eMERGE files are not overwritten.",
        ],
    }
    SUMMARY_OUT.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
