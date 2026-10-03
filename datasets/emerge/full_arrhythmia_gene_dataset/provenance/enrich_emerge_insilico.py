#!/usr/bin/env python3
"""Retrieve MyVariant/dbNSFP in silico annotations for eMERGE validation rows.

Duplicate policy:
- Prefer eMERGE-provided values when the same field exists in the source table
  (currently CADD and gnomAD frequency).
- Fill missing fields and extra predictors from MyVariant/dbNSFP.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from urllib.parse import urlencode

from cardio_pipeline.annotation import ApiClient, MYVARIANT_API, dbnsfp_scores, extract_dbnsfp_from_myvariant


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "processed"
RESULTS = ROOT / "results"


INSILICO_FIELDS = [
    "revel",
    "cadd_phred",
    "sift",
    "polyphen2_hdiv",
    "metalr",
    "fathmm_xf",
    "alphamissense",
]


def clean(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "na", "<na>", "*"}:
        return ""
    return text


def num_or_blank(value: Any) -> str:
    text = clean(value)
    if not text:
        return ""
    try:
        return f"{float(text):.12g}"
    except ValueError:
        return ""


def first_nonempty(value: Any) -> str:
    if isinstance(value, list):
        for item in value:
            text = first_nonempty(item)
            if text:
                return text
        return ""
    if isinstance(value, dict):
        for item in value.values():
            text = first_nonempty(item)
            if text:
                return text
        return ""
    return clean(value)


def strip_hgvs_c(value: Any) -> str:
    text = first_nonempty(value)
    if not text:
        return ""
    if ":" in text:
        parts = [part for part in text.split(":") if part.startswith("c.")]
        if parts:
            return parts[-1]
        text = text.rsplit(":", 1)[-1]
    return text if text.startswith("c.") else ""


def compact_hgvs_p(value: Any) -> str:
    text = first_nonempty(value)
    if not text:
        return ""
    aa1 = {
        "Ala": "A",
        "Arg": "R",
        "Asn": "N",
        "Asp": "D",
        "Cys": "C",
        "Gln": "Q",
        "Glu": "E",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Leu": "L",
        "Lys": "K",
        "Met": "M",
        "Phe": "F",
        "Pro": "P",
        "Ser": "S",
        "Thr": "T",
        "Trp": "W",
        "Tyr": "Y",
        "Val": "V",
        "Ter": "X",
    }
    import re

    match = re.match(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|\*)$", text)
    if match:
        ref, pos, alt = match.groups()
        return f"p.{aa1.get(ref, ref)}{pos}{'X' if alt == '*' else aa1.get(alt, alt)}"
    return text


def coords_from_row(row: pd.Series) -> dict[str, Any]:
    chrom = clean(row.get("CHROM"))
    pos = clean(row.get("POS"))
    ref = clean(row.get("REF"))
    alt = clean(row.get("ALT"))
    if not (chrom and pos and ref and alt):
        return {}
    try:
        pos_int = int(float(pos))
    except ValueError:
        return {}
    return {"chrom": chrom.replace("chr", ""), "pos": pos_int, "ref": ref, "alt": alt, "coord_source": "emerge_s2"}


def hgvs_from_row(row: pd.Series) -> dict[str, str]:
    hgvs_p = clean(row.get("hgvs_p_inferred"))
    one = clean(row.get("mutName3"))
    if one and not one.startswith("p."):
        one = f"p.{one}"
    return {"hgvs_c": "", "hgvs_p": hgvs_p, "hgvs_p_one_letter": one}


def query_myvariant_hg19(
    client: ApiClient,
    coords: dict[str, Any],
    gene: str,
    hgvs: dict[str, str],
) -> dict[str, Any]:
    queries: list[tuple[str, dict[str, Any]]] = []
    if coords:
        hgvs_id = f"chr{coords['chrom']}:g.{coords['pos']}{coords['ref']}>{coords['alt']}"
        queries.append((f"{MYVARIANT_API}/variant/{hgvs_id}", {"fields": "dbnsfp", "assembly": "hg19"}))
    if gene and hgvs.get("hgvs_p"):
        queries.append((f"{MYVARIANT_API}/query", {"q": f"{gene} {hgvs['hgvs_p']}", "fields": "dbnsfp", "size": 5}))

    last_error = ""
    for url, params in queries:
        data, error = client.request_json("GET", url, params=params)
        if error:
            last_error = error
            continue
        dbnsfp = extract_dbnsfp_from_myvariant(data or {}, [gene])
        if not dbnsfp:
            continue
        scores = dbnsfp_scores(dbnsfp)
        hg38 = dbnsfp.get("hg38", {})
        if isinstance(hg38, list):
            hg38 = hg38[0] if hg38 else {}
        chrom = dbnsfp.get("chrom")
        ref = dbnsfp.get("ref")
        alt = dbnsfp.get("alt")
        if isinstance(chrom, list):
            chrom = chrom[0] if chrom else ""
        if isinstance(ref, list):
            ref = ref[0] if ref else ""
        if isinstance(alt, list):
            alt = alt[0] if alt else ""
        hgvsc = first_nonempty(dbnsfp.get("hgvsc"))
        hgvsp = first_nonempty(dbnsfp.get("hgvsp"))
        return {
            "myvariant_status": "matched" if any(value is not None for value in scores.values()) else "no_scores",
            "myvariant_query": f"{url}?{urlencode(params)}",
            "myvariant_error": "",
            "hg38_chrom": clean(chrom).replace("chr", ""),
            "hg38_pos": clean(hg38.get("start") if isinstance(hg38, dict) else ""),
            "hg38_ref": clean(ref),
            "hg38_alt": clean(alt),
            "hgvsc": hgvsc,
            "hgvsp": hgvsp,
            **scores,
        }

    return {
        "myvariant_status": "error" if last_error else "no_match",
        "myvariant_query": "; ".join(f"{url}?{urlencode(params)}" for url, params in queries),
        "myvariant_error": last_error,
    }


def update_feature_file(feature_path: Path, enriched_source: pd.DataFrame, output_path: Path) -> None:
    features = pd.read_csv(feature_path, dtype=str, keep_default_na=False)
    source = enriched_source.reset_index(drop=True)
    if len(features) != len(source):
        raise ValueError(f"Row count mismatch: {feature_path} has {len(features)} rows, source has {len(source)}")

    for field in INSILICO_FIELDS:
        if field in features.columns and f"final_{field}" in source.columns:
            features[field] = source[f"final_{field}"].map(clean)
    for field in ["gnomad_max_af", "gnomad_status", "gnomad_absent", "gnomad_max_af_source"]:
        if field in features.columns and f"final_{field}" in source.columns:
            features[field] = source[f"final_{field}"].map(clean)
    for field in ["chrom", "pos", "ref", "alt"]:
        if field in features.columns and f"final_{field}" in source.columns:
            features[field] = source[f"final_{field}"].map(clean).replace({"": "Missing"} if field != "pos" else {"": ""})
    for field in ["gene_refgene", "hgvs_c", "hgvs_p", "hgvs_p_one_letter"]:
        if field in features.columns and f"final_{field}" in source.columns:
            features[field] = source[f"final_{field}"].map(clean).replace({"": "Missing"})

    output_path.parent.mkdir(parents=True, exist_ok=True)
    features.to_csv(output_path, index=False)


def choose_value(dataset_value: Any, retrieved_value: Any) -> tuple[str, str]:
    dataset_clean = num_or_blank(dataset_value)
    retrieved_clean = num_or_blank(retrieved_value)
    if dataset_clean:
        return dataset_clean, "dataset"
    if retrieved_clean:
        return retrieved_clean, "retrieved_myvariant_dbnsfp"
    return "", "missing"


def enrich(
    source_csv: Path,
    baseline_features: Path,
    enriched_features: Path,
    output_stem: str,
    *,
    limit: int | None,
    progress_every: int,
) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    PROCESSED.mkdir(parents=True, exist_ok=True)

    source = pd.read_csv(source_csv, dtype=str, keep_default_na=False)
    if limit:
        source = source.head(limit).copy()

    client = ApiClient()
    records = []
    for idx, row in source.iterrows():
        coords = coords_from_row(row)
        hgvs = hgvs_from_row(row)
        genes = [clean(row.get("gene") or row.get("GL2"))]
        clinvar_id = clean(row.get("CV_ALLELEID"))

        myvariant = query_myvariant_hg19(client, coords, genes[0], hgvs)

        record = row.to_dict()
        record.update({f"retrieved_{k}": clean(v) for k, v in myvariant.items()})
        record["retrieved_gnomad_status"] = "not_queried_external_emerge_values_preferred"
        record["source_hg19_chrom"] = clean(row.get("CHROM"))
        record["source_hg19_pos"] = clean(row.get("POS"))
        record["source_hg19_ref"] = clean(row.get("REF"))
        record["source_hg19_alt"] = clean(row.get("ALT"))

        # Duplicate policy: if both dataset and retrieved values exist, dataset wins.
        cadd_value, cadd_source = choose_value(row.get("CADD"), myvariant.get("cadd_phred"))
        gnomad_dataset_value = num_or_blank(row.get("gnomadE_AF_max_faf99")) or num_or_blank(row.get("gnomadEG_AF"))
        gnomad_value, gnomad_source = choose_value(gnomad_dataset_value, "")

        record["final_cadd_phred"] = cadd_value
        record["final_cadd_phred_source"] = cadd_source
        record["final_gnomad_max_af"] = gnomad_value
        record["final_gnomad_max_af_source_policy"] = gnomad_source
        record["final_gnomad_status"] = "matched" if record["final_gnomad_max_af"] else "absent"
        record["final_gnomad_absent"] = "false" if record["final_gnomad_max_af"] else "true"
        record["final_gnomad_max_af_source"] = (
            "emerge_gnomadE_AF_max_faf99"
            if num_or_blank(row.get("gnomadE_AF_max_faf99"))
            else "emerge_gnomadEG_AF"
            if num_or_blank(row.get("gnomadEG_AF"))
            else "absent"
        )
        record["final_chrom"] = clean(myvariant.get("hg38_chrom")) or record["source_hg19_chrom"]
        record["final_pos"] = num_or_blank(myvariant.get("hg38_pos")) or num_or_blank(row.get("POS"))
        record["final_ref"] = clean(myvariant.get("hg38_ref")) or record["source_hg19_ref"]
        record["final_alt"] = clean(myvariant.get("hg38_alt")) or record["source_hg19_alt"]
        record["final_coordinate_source"] = "retrieved_myvariant_dbnsfp_hg38" if clean(myvariant.get("hg38_pos")) else "dataset_hg19_fallback"
        record["final_hgvs_c"] = strip_hgvs_c(myvariant.get("hgvsc"))
        record["final_gene_refgene"] = (
            f"{genes[0]}:{clean(myvariant.get('hgvsc'))}"
            if clean(myvariant.get("hgvsc"))
            else clean(row.get("GM2"))
        )
        record["final_hgvs_p"] = clean(row.get("hgvs_p_inferred")) or compact_hgvs_p(myvariant.get("hgvsp"))
        record["final_hgvs_p_one_letter"] = clean(row.get("hgvs_p_one_letter_model")) or compact_hgvs_p(myvariant.get("hgvsp"))
        record["final_hgvs_c_source"] = "retrieved_myvariant_dbnsfp" if record["final_hgvs_c"] else "missing"
        record["final_gene_refgene_source"] = "retrieved_myvariant_dbnsfp" if clean(myvariant.get("hgvsc")) else "dataset"
        record["final_hgvs_p_source"] = "dataset" if clean(row.get("hgvs_p_inferred")) else ("retrieved_myvariant_dbnsfp" if record["final_hgvs_p"] else "missing")
        for field in INSILICO_FIELDS:
            if field in {"cadd_phred"}:
                continue
            # eMERGE S2 does not ship these predictors, but keep this policy generic.
            value, source_name = choose_value(row.get(field), myvariant.get(field))
            record[f"final_{field}"] = value
            record[f"final_{field}_source"] = source_name

        records.append(record)
        if len(records) == 1 or len(records) % progress_every == 0 or len(records) == len(source):
            print(
                f"[insilico] {len(records)}/{len(source)} {genes[0]} {clean(row.get('mutName3'))} "
                f"myvariant={myvariant.get('myvariant_status')} gnomad=external_preferred",
                flush=True,
            )

    enriched_source = pd.DataFrame(records)
    enriched_source_path = PROCESSED / f"{output_stem}_insilico_enriched_source.csv"
    enriched_source.to_csv(enriched_source_path, index=False)

    audit_cols = [
        "gene",
        "mutName3",
        "hgvs_p_inferred",
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "CV_ALLELEID",
        "CV_CLNSIG",
        "CV_CLNSIG4",
        "CV_CLNREVSTAT",
        "retrieved_myvariant_status",
        "retrieved_myvariant_query",
        "retrieved_myvariant_error",
        "source_hg19_chrom",
        "source_hg19_pos",
        "source_hg19_ref",
        "source_hg19_alt",
        "retrieved_hg38_chrom",
        "retrieved_hg38_pos",
        "retrieved_hg38_ref",
        "retrieved_hg38_alt",
        "retrieved_hgvsc",
        "retrieved_hgvsp",
        "final_chrom",
        "final_pos",
        "final_ref",
        "final_alt",
        "final_coordinate_source",
        "final_gene_refgene",
        "final_gene_refgene_source",
        "final_hgvs_c",
        "final_hgvs_c_source",
        "final_hgvs_p",
        "final_hgvs_p_source",
        "CADD",
        "final_cadd_phred",
        "final_cadd_phred_source",
        "gnomadE_AF_max_faf99",
        "gnomadEG_AF",
        "final_gnomad_max_af",
        "final_gnomad_max_af_source_policy",
    ]
    for field in ["revel", "sift", "polyphen2_hdiv", "metalr", "fathmm_xf", "alphamissense"]:
        audit_cols.extend([f"retrieved_{field}", f"final_{field}", f"final_{field}_source"])
    enriched_source[[col for col in audit_cols if col in enriched_source.columns]].to_csv(
        PROCESSED / f"{output_stem}_insilico_source_priority_audit.csv",
        index=False,
    )

    update_feature_file(
        baseline_features,
        enriched_source,
        PROCESSED / f"{output_stem}_ml_baseline_features_insilico_enriched.csv",
    )
    update_feature_file(
        enriched_features,
        enriched_source,
        PROCESSED / f"{output_stem}_ml_enriched_features_insilico_enriched.csv",
    )

    summary = {
        "source_csv": str(source_csv),
        "rows": int(len(enriched_source)),
        "output_stem": output_stem,
        "myvariant_status_counts": enriched_source["retrieved_myvariant_status"].value_counts(dropna=False).to_dict(),
        "gnomad_status_counts": enriched_source["retrieved_gnomad_status"].value_counts(dropna=False).to_dict(),
        "final_non_missing": {
            field: int(enriched_source[f"final_{field}"].map(clean).astype(bool).sum())
            for field in INSILICO_FIELDS + ["gnomad_max_af"]
            if f"final_{field}" in enriched_source.columns
        },
        "final_hgvs_non_missing": {
            field: int(enriched_source[f"final_{field}"].map(clean).astype(bool).sum())
            for field in ["gene_refgene", "hgvs_c", "hgvs_p", "hgvs_p_one_letter"]
            if f"final_{field}" in enriched_source.columns
        },
        "duplicate_policy": "eMERGE-provided CADD and gnomAD values are retained when present; MyVariant/dbNSFP fills the remaining in silico predictors.",
        "clinvar_policy": "eMERGE CV_* fields are retained as shipped; CV_ALLELEID is not treated as an NCBI ClinVar Variation ID for ESummary lookup.",
    }
    (RESULTS / f"{output_stem}_insilico_enrichment_summary.json").write_text(json.dumps(summary, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scope",
        choices=[
            "candidate50",
            "candidate100",
            "candidate100_matched",
            "candidate300_matched",
            "trainratio_strict",
            "trainratio_300_relaxed",
            "full",
        ],
        default="candidate50",
    )
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--progress-every", type=int, default=50)
    args = parser.parse_args()

    if args.scope == "candidate50":
        enrich(
            PROCESSED / "emerge_candidate_50.csv",
            PROCESSED / "emerge_candidate_50_ml_baseline_features.csv",
            PROCESSED / "emerge_candidate_50_ml_enriched_features.csv",
            "emerge_candidate_50",
            limit=args.limit,
            progress_every=args.progress_every,
        )
    elif args.scope == "candidate100":
        enrich(
            PROCESSED / "emerge_candidate_100.csv",
            PROCESSED / "emerge_candidate_100_ml_baseline_features.csv",
            PROCESSED / "emerge_candidate_100_ml_enriched_features.csv",
            "emerge_candidate_100",
            limit=args.limit,
            progress_every=args.progress_every,
        )
    elif args.scope == "candidate100_matched":
        enrich(
            PROCESSED / "emerge_candidate_100_training_matched.csv",
            PROCESSED / "emerge_candidate_100_training_matched_ml_baseline_features.csv",
            PROCESSED / "emerge_candidate_100_training_matched_ml_enriched_features.csv",
            "emerge_candidate_100_training_matched",
            limit=args.limit,
            progress_every=args.progress_every,
        )
    elif args.scope == "candidate300_matched":
        enrich(
            PROCESSED / "emerge_candidate_300_training_matched.csv",
            PROCESSED / "emerge_candidate_300_training_matched_ml_baseline_features.csv",
            PROCESSED / "emerge_candidate_300_training_matched_ml_enriched_features.csv",
            "emerge_candidate_300_training_matched",
            limit=args.limit,
            progress_every=args.progress_every,
        )
    elif args.scope == "trainratio_strict":
        enrich(
            PROCESSED / "emerge_candidate_trainratio_strict.csv",
            PROCESSED / "emerge_candidate_trainratio_strict_ml_baseline_features.csv",
            PROCESSED / "emerge_candidate_trainratio_strict_ml_enriched_features.csv",
            "emerge_candidate_trainratio_strict",
            limit=args.limit,
            progress_every=args.progress_every,
        )
    elif args.scope == "trainratio_300_relaxed":
        enrich(
            PROCESSED / "emerge_candidate_trainratio_300_relaxed.csv",
            PROCESSED / "emerge_candidate_trainratio_300_relaxed_ml_baseline_features.csv",
            PROCESSED / "emerge_candidate_trainratio_300_relaxed_ml_enriched_features.csv",
            "emerge_candidate_trainratio_300_relaxed",
            limit=args.limit,
            progress_every=args.progress_every,
        )
    else:
        enrich(
            PROCESSED / "emerge_arrhythmia_variant_level_full.csv",
            PROCESSED / "emerge_arrhythmia_ml_baseline_features.csv",
            PROCESSED / "emerge_arrhythmia_ml_enriched_features.csv",
            "emerge_arrhythmia_full",
            limit=args.limit,
            progress_every=args.progress_every,
        )


if __name__ == "__main__":
    main()
