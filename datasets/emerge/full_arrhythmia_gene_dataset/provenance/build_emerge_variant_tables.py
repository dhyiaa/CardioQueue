#!/usr/bin/env python3
"""Build eMERGE external-validation seed tables from public supplements."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw"
PROCESSED = ROOT / "processed"
RESULTS = ROOT / "results"

FILE_S2 = RAW / "NIHMS1778916-supplement-File_S2_Variant_Summary.xlsx"
FILE_S3 = RAW / "NIHMS1778916-supplement-File_S3_Variants_with_disagreements_between_sequencing_centers.xlsx"
FILE_S4 = RAW / "NIHMS1778916-supplement-File_S4_In_vitro_EP_dataset.xlsx"

ARRHYTHMIA_GENES = {
    "ANK2",
    "CACNA1C",
    "KCNE1",
    "KCNE2",
    "KCNH2",
    "KCNJ2",
    "KCNQ1",
    "LMNA",
    "RYR2",
    "SCN5A",
}

CLASS_TO_NUMERIC = {"Benign": 0, "VUS": 1, "Pathogenic": 2}

IDENTIFIER_COLUMNS = [
    "variant_uid",
    "variant_key",
    "study",
    "source_row",
    "participant_id",
    "patient_match",
]

TARGET_COLUMNS = [
    "target_3class",
    "target_numeric",
    "target_5class",
    "target_classification_raw",
]

BASE_NUMERIC_FEATURES = [
    "pos",
    "gnomad_max_af",
    "revel",
    "cadd_phred",
    "sift",
    "polyphen2_hdiv",
    "metalr",
    "fathmm_xf",
    "alphamissense",
    "patient_age_at_diagnosis",
    "patient_followup_presyncope_yes_count",
    "patient_followup_presyncope_uncertain_count",
    "patient_followup_presyncope_known_count",
    "patient_followup_syncope_yes_count",
    "patient_followup_syncope_uncertain_count",
    "patient_followup_syncope_known_count",
    "patient_followup_palpitations_yes_count",
    "patient_followup_palpitations_uncertain_count",
    "patient_followup_palpitations_known_count",
    "patient_followup_chest_pain_yes_count",
    "patient_followup_chest_pain_uncertain_count",
    "patient_followup_chest_pain_known_count",
    "patient_followup_cardiac_arrest_yes_count",
    "patient_followup_cardiac_arrest_uncertain_count",
    "patient_followup_cardiac_arrest_known_count",
    "patient_followup_death_yes_count",
    "patient_followup_death_uncertain_count",
    "patient_followup_death_known_count",
    "patient_followup_other_symptom_yes_count",
    "patient_followup_other_symptom_uncertain_count",
    "patient_followup_other_symptom_known_count",
    "patient_followup_circumstances_known_count",
    "patient_qt_first",
    "patient_qt_latest",
    "patient_qt_min",
    "patient_qt_max",
    "patient_qt_count",
    "patient_lvef_first",
    "patient_lvef_latest",
    "patient_lvef_min",
    "patient_lvef_max",
    "patient_lvef_count",
]

BASE_CATEGORICAL_FEATURES = [
    "gene",
    "gene_refgene",
    "hgvs_c",
    "hgvs_p",
    "hgvs_p_one_letter",
    "inferred_consequence",
    "chrom",
    "ref",
    "alt",
    "gnomad_status",
    "gnomad_absent",
    "gnomad_max_af_source",
    "patient_sex",
    "patient_sym_presyncope",
    "patient_sym_syncope",
    "patient_sym_palpitations",
    "patient_sym_chest_pain",
    "patient_sym_cardiac_arrest",
    "patient_sym_death",
    "patient_sym_other",
    "patient_sym_circumstances",
    "patient_family_history_flags",
    "patient_followup_presyncope_any",
    "patient_followup_presyncope_latest",
    "patient_followup_syncope_any",
    "patient_followup_syncope_latest",
    "patient_followup_palpitations_any",
    "patient_followup_palpitations_latest",
    "patient_followup_chest_pain_any",
    "patient_followup_chest_pain_latest",
    "patient_followup_cardiac_arrest_any",
    "patient_followup_cardiac_arrest_latest",
    "patient_followup_death_any",
    "patient_followup_death_latest",
    "patient_followup_other_symptom_any",
    "patient_followup_other_symptom_latest",
    "patient_followup_circumstances_latest",
]

ENRICHED_NUMERIC_FEATURES = ["clinvar_stars", "clinvar_submitters"]
ENRICHED_CATEGORICAL_FEATURES = [
    "clinvar_status",
    "clinvar_classification",
    "clinvar_review_status",
    "clinvar_conflicting",
]


def clean_text(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"nan", "none", "na", "<na>"}:
        return ""
    return text


def normalize_label(value: Any) -> str:
    text = clean_text(value)
    if text in {"Pathogenic/Likely_pathogenic", "Pathogenic", "Likely_pathogenic"}:
        return "Pathogenic"
    if text in {"Benign/Likely_benign", "Benign", "Likely_benign"}:
        return "Benign"
    if text in {"VUS_conflicting", "Uncertain_significance", "Conflicting_interpretations_of_pathogenicity"}:
        return "VUS"
    return ""


def variant_key(row: pd.Series) -> str:
    return "|".join(
        clean_text(row.get(col))
        for col in ["GL2", "CHROM", "POS", "REF", "ALT"]
    )


def protein_key(row: pd.Series) -> str:
    return f"{clean_text(row.get('GL2'))}|{clean_text(row.get('mutName3'))}"


def map_consequence(value: Any) -> str:
    text = clean_text(value).lower()
    if not text:
        return "Missing"
    if "missense" in text:
        return "missense"
    if "frameshift" in text:
        return "frameshift"
    if "splice-acceptor" in text or "splice-donor" in text or "near-splice" in text:
        return "splice_region"
    if "stop-gained" in text or "stop-lost" in text:
        return "nonsense"
    if any(token in text for token in ["coding", "intron", "utr", "synonymous"]):
        return "Missing"
    return text


def infer_hgvs_p(gene: str, token: str) -> str:
    token = clean_text(token)
    if not token:
        return ""
    if token.lower().startswith("c."):
        return ""
    if token.startswith("p."):
        return token
    if re.match(r"^[A-Z][0-9]+[A-Z*]$", token):
        aa3 = {
            "A": "Ala",
            "R": "Arg",
            "N": "Asn",
            "D": "Asp",
            "C": "Cys",
            "Q": "Gln",
            "E": "Glu",
            "G": "Gly",
            "H": "His",
            "I": "Ile",
            "L": "Leu",
            "K": "Lys",
            "M": "Met",
            "F": "Phe",
            "P": "Pro",
            "S": "Ser",
            "T": "Thr",
            "W": "Trp",
            "Y": "Tyr",
            "V": "Val",
            "*": "Ter",
        }
        return f"p.{aa3.get(token[0], token[0])}{token[1:-1]}{aa3.get(token[-1], token[-1])}"
    return token


def hgvs_p_to_compact(token: str) -> str:
    token = clean_text(token)
    if not token.startswith("p."):
        return token
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
    match = re.match(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter)$", token)
    if not match:
        return token
    ref, pos, alt = match.groups()
    return f"{aa1.get(ref, ref)}{pos}{aa1.get(alt, alt)}"


def numeric_or_blank(value: Any) -> str:
    text = clean_text(value)
    if not text or text == "*":
        return ""
    try:
        return f"{float(text):.12g}"
    except ValueError:
        return ""


def categorical_or_missing(value: Any) -> str:
    text = clean_text(value)
    return text if text else "Missing"


def build_model_ready(df: pd.DataFrame, *, enriched: bool) -> pd.DataFrame:
    rows = []
    for i, row in df.reset_index(drop=True).iterrows():
        label = clean_text(row.get("reference_label_3class"))
        af = numeric_or_blank(row.get("gnomadE_AF_max_faf99")) or numeric_or_blank(row.get("gnomadEG_AF"))
        gnomad_status = "matched" if af else "absent"
        cv = clean_text(row.get("CV_CLNSIG4"))
        clinvar_status = "matched" if clean_text(row.get("CV_ALLELEID")) else "no_match"
        clinvar_conflicting = "true" if "conflicting" in cv.lower() else "false"
        out = {
            "variant_uid": f"eMERGE-III|{clean_text(row.get('GL2'))}|{clean_text(row.get('CHROM'))}:{clean_text(row.get('POS'))}:{clean_text(row.get('REF'))}>{clean_text(row.get('ALT'))}|{i + 1}",
            "variant_key": clean_text(row.get("variant_key_external")),
            "study": "eMERGE-III_Glazer",
            "source_row": str(i + 1),
            "participant_id": f"EMERGE_VARIANT_{i + 1:04d}",
            "patient_match": "False",
            "target_3class": label,
            "target_numeric": CLASS_TO_NUMERIC.get(label, ""),
            "target_5class": clean_text(row.get("finalAnnotation_postInVitro")),
            "target_classification_raw": clean_text(row.get("finalAnnotation2_postInvitro")),
            "pos": numeric_or_blank(row.get("POS")),
            "gnomad_max_af": af,
            "revel": "",
            "cadd_phred": numeric_or_blank(row.get("CADD")),
            "sift": "",
            "polyphen2_hdiv": "",
            "metalr": "",
            "fathmm_xf": "",
            "alphamissense": "",
            "gene": categorical_or_missing(row.get("GL2")),
            "gene_refgene": categorical_or_missing(row.get("GM2")),
            "hgvs_c": "Missing",
            "hgvs_p": categorical_or_missing(row.get("hgvs_p_inferred")),
            "hgvs_p_one_letter": categorical_or_missing(row.get("hgvs_p_one_letter_model")),
            "inferred_consequence": categorical_or_missing(row.get("inferred_consequence_model")),
            "chrom": categorical_or_missing(row.get("CHROM")),
            "ref": categorical_or_missing(row.get("REF")),
            "alt": categorical_or_missing(row.get("ALT")),
            "gnomad_status": gnomad_status,
            "gnomad_absent": "false" if af else "true",
            "gnomad_max_af_source": "gnomadE_AF_max_faf99" if numeric_or_blank(row.get("gnomadE_AF_max_faf99")) else ("gnomadEG_AF" if af else "absent"),
            "clinvar_stars": "",
            "clinvar_submitters": "",
            "clinvar_status": clinvar_status,
            "clinvar_classification": categorical_or_missing(row.get("CV_CLNSIG4")),
            "clinvar_review_status": categorical_or_missing(row.get("CV_CLNREVSTAT")),
            "clinvar_conflicting": clinvar_conflicting,
        }
        for col in BASE_NUMERIC_FEATURES:
            out.setdefault(col, "")
        for col in BASE_CATEGORICAL_FEATURES:
            out.setdefault(col, "Missing")
        columns = IDENTIFIER_COLUMNS + TARGET_COLUMNS + BASE_NUMERIC_FEATURES + BASE_CATEGORICAL_FEATURES
        if enriched:
            columns = columns + ENRICHED_NUMERIC_FEATURES + ENRICHED_CATEGORICAL_FEATURES
        rows.append({col: out.get(col, "") for col in columns})
    return pd.DataFrame(rows)


def collect_functional_variant_keys() -> set[str]:
    keys: set[str] = set()
    if not FILE_S4.exists():
        return keys
    xl = pd.ExcelFile(FILE_S4)
    for sheet in xl.sheet_names:
        df = pd.read_excel(FILE_S4, sheet_name=sheet, header=None, dtype=str)
        gene = ""
        if "KCNQ1" in sheet:
            gene = "KCNQ1"
        elif "SCN5A" in sheet:
            gene = "SCN5A"
        elif "KCNH2" in sheet:
            gene = "KCNH2"
        if not gene:
            continue
        first_col = df.iloc[:, 0].map(clean_text)
        for token in first_col:
            if not token or token.lower() in {"variant", "mutation", "wildtype", "wt"}:
                continue
            if any(part in token.lower() for part in ["mean", "nan"]):
                continue
            keys.add(f"{gene}|{token}")
            keys.add(f"{gene}|{hgvs_p_to_compact(token)}")
    return keys


def build_tables() -> None:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)

    s2 = pd.read_excel(FILE_S2, dtype=str)
    s2["source_study"] = "eMERGE-III Glazer Circulation"
    s2["variant_key_external"] = s2.apply(variant_key, axis=1)
    s2["protein_key_external"] = s2.apply(protein_key, axis=1)
    s2["gene"] = s2["GL2"].map(clean_text)
    s2["hgvs_p_inferred"] = [
        infer_hgvs_p(gene, token) for gene, token in zip(s2["gene"], s2["mutName3"])
    ]
    s2["hgvs_p_one_letter_model"] = s2["mutName3"].map(
        lambda x: (
            f"p.{clean_text(x)}"
            if clean_text(x) and not clean_text(x).startswith("p.") and not clean_text(x).lower().startswith("c.")
            else clean_text(x)
            if not clean_text(x).lower().startswith("c.")
            else ""
        )
    )
    s2["inferred_consequence_model"] = s2["FG2"].map(map_consequence)
    s2["reference_label_3class"] = s2["finalAnnotation2_postInvitro"].map(normalize_label)
    s2["reference_label_3class_pre_invitro"] = s2["finalAnnotation2_preInvitro"].map(normalize_label)
    s2["is_arrhythmia_gene"] = s2["gene"].isin(ARRHYTHMIA_GENES)
    s2["is_canonical"] = s2["Canonical"].map(clean_text).eq("True")
    s2["emerge_public_phenotype_row_level_available"] = False

    s3 = pd.read_excel(FILE_S3, dtype=str)
    disagreement_keys = set(s3.apply(variant_key, axis=1))
    s2["in_file_s3_sequencing_center_disagreement"] = s2["variant_key_external"].isin(disagreement_keys)

    functional_keys = collect_functional_variant_keys()
    s2["in_file_s4_functional_ep_dataset"] = s2["protein_key_external"].isin(functional_keys)

    arr = s2[
        s2["is_arrhythmia_gene"]
        & s2["reference_label_3class"].isin(["Benign", "VUS", "Pathogenic"])
    ].copy()
    arr["primary_validation_eligible"] = (
        arr["is_canonical"]
        & ~arr["in_file_s3_sequencing_center_disagreement"]
    )

    # Manuscript-friendly candidate sets. These are intentionally balanced because
    # the full public supplement is overwhelmingly VUS.
    # Prioritize functional EP variants, higher carrier count, canonical transcripts, and no disagreement.
    arr["_n_numeric"] = pd.to_numeric(arr["n"], errors="coerce").fillna(0)
    arr["_cadd_numeric"] = pd.to_numeric(arr["CADD"], errors="coerce").fillna(-1)
    arr["_n_capped"] = arr["_n_numeric"].clip(upper=10)
    arr["_priority"] = (
        arr["in_file_s4_functional_ep_dataset"].astype(int) * 1000
        + arr["primary_validation_eligible"].astype(int) * 100
        + arr["_cadd_numeric"]
        + arr["_n_capped"]
    )
    def select_candidate(quotas: dict[str, int]) -> pd.DataFrame:
        candidate_parts = []
        for label, quota in quotas.items():
            part = arr[arr["reference_label_3class"].eq(label)].copy()
            if label == "VUS":
                part = part[~part["mutName3"].isin(["D85N", "T8A"])]
            part = part.sort_values(["_priority", "_cadd_numeric", "_n_numeric"], ascending=False).head(quota)
            candidate_parts.append(part)
        return pd.concat(candidate_parts, ignore_index=True).drop(
            columns=["_n_numeric", "_n_capped", "_cadd_numeric", "_priority"]
        )

    candidate = select_candidate({"Pathogenic": 20, "VUS": 20, "Benign": 10})
    candidate100 = select_candidate({"Pathogenic": 40, "VUS": 40, "Benign": 20})

    training_like_gene_targets = {
        "Pathogenic": {"RYR2": 13, "SCN5A": 10, "CACNA1C": 7, "KCNH2": 4, "KCNQ1": 3, "LMNA": 3},
        "VUS": {"RYR2": 10, "ANK2": 8, "SCN5A": 8, "CACNA1C": 5, "KCNH2": 4, "KCNE1": 3, "KCNQ1": 1, "KCNJ2": 1},
        "Benign": {"SCN5A": 7, "KCNE1": 6, "KCNH2": 4, "KCNQ1": 1, "RYR2": 1, "ANK2": 1},
    }

    def select_training_like_100() -> pd.DataFrame:
        picked = []
        used = set()
        pool = arr[arr["primary_validation_eligible"]].copy()
        pool["_is_model_vocab_consequence"] = pool["inferred_consequence_model"].isin(
            ["missense", "frameshift", "nonsense", "splice_region", "Missing"]
        )
        pool["_training_like_priority"] = (
            pool["_is_model_vocab_consequence"].astype(int) * 10000
            + pool["in_file_s4_functional_ep_dataset"].astype(int) * 1000
            + pool["inferred_consequence_model"].eq("missense").astype(int) * 250
            + pool["_cadd_numeric"]
            + pool["_n_capped"]
        )
        for label, gene_targets in training_like_gene_targets.items():
            label_pool = pool[pool["reference_label_3class"].eq(label)].copy()
            for gene, quota in gene_targets.items():
                part = (
                    label_pool[label_pool["gene"].eq(gene)]
                    .sort_values(["_training_like_priority", "_cadd_numeric", "_n_numeric"], ascending=False)
                )
                for _, row in part.iterrows():
                    key = row["variant_key_external"]
                    if key in used:
                        continue
                    picked.append(row)
                    used.add(key)
                    if sum(1 for item in picked if item["reference_label_3class"] == label and item["gene"] == gene) >= quota:
                        break
        desired = {"Pathogenic": 40, "VUS": 40, "Benign": 20}
        for label, quota in desired.items():
            current = sum(1 for item in picked if item["reference_label_3class"] == label)
            if current >= quota:
                continue
            label_pool = pool[pool["reference_label_3class"].eq(label)].sort_values(
                ["_training_like_priority", "_cadd_numeric", "_n_numeric"], ascending=False
            )
            for _, row in label_pool.iterrows():
                key = row["variant_key_external"]
                if key in used:
                    continue
                picked.append(row)
                used.add(key)
                current += 1
                if current >= quota:
                    break
        out = pd.DataFrame(picked)
        order = {"Pathogenic": 0, "VUS": 1, "Benign": 2}
        out["_label_order"] = out["reference_label_3class"].map(order)
        out = out.sort_values(["_label_order", "gene", "_training_like_priority"], ascending=[True, True, False])
        drop_cols = [c for c in out.columns if c.startswith("_")]
        return out.drop(columns=drop_cols).reset_index(drop=True)

    candidate100_matched = select_training_like_100()

    training_gene_weights = {
        "SCN5A": 49,
        "KCNE1": 36,
        "KCNH2": 26,
        "RYR2": 16,
        "ANK2": 16,
        "CACNA1C": 14,
        "KCNQ1": 8,
        "LMNA": 3,
        "KCNJ2": 2,
        "KCNE2": 1,
    }

    def select_training_like_balanced(quotas: dict[str, int]) -> pd.DataFrame:
        picked = []
        used = set()
        pool = arr[arr["primary_validation_eligible"]].copy()
        pool["_is_model_vocab_consequence"] = pool["inferred_consequence_model"].isin(
            ["missense", "frameshift", "nonsense", "splice_region", "Missing"]
        )
        pool["_training_gene_weight"] = pool["gene"].map(training_gene_weights).fillna(0)
        pool["_training_like_priority"] = (
            pool["_is_model_vocab_consequence"].astype(int) * 10000
            + pool["inferred_consequence_model"].eq("missense").astype(int) * 2000
            + pool["in_file_s4_functional_ep_dataset"].astype(int) * 1000
            + pool["_training_gene_weight"] * 10
            + pool["_cadd_numeric"]
            + pool["_n_capped"]
        )
        for label, quota in quotas.items():
            label_pool = pool[pool["reference_label_3class"].eq(label)].sort_values(
                ["_training_like_priority", "_cadd_numeric", "_n_numeric"],
                ascending=False,
            )
            current = 0
            for _, row in label_pool.iterrows():
                key = row["variant_key_external"]
                if key in used:
                    continue
                picked.append(row)
                used.add(key)
                current += 1
                if current >= quota:
                    break
            if current < quota:
                raise ValueError(f"Could only select {current}/{quota} rows for {label}")
        out = pd.DataFrame(picked)
        order = {"Pathogenic": 0, "VUS": 1, "Benign": 2}
        out["_label_order"] = out["reference_label_3class"].map(order)
        out = out.sort_values(["_label_order", "_training_like_priority"], ascending=[True, False])
        drop_cols = [c for c in out.columns if c.startswith("_")]
        return out.drop(columns=drop_cols).reset_index(drop=True)

    candidate300_matched = select_training_like_balanced({"Pathogenic": 100, "VUS": 100, "Benign": 100})

    training_consequence_weights = {
        "Benign": {"missense": 165, "Missing": 34, "indel": 13, "splice_region": 8},
        "VUS": {"missense": 189, "frameshift": 5, "indel": 4, "nonsense": 4, "splice_region": 4, "Missing": 1},
        "Pathogenic": {"missense": 27, "frameshift": 11, "nonsense": 12, "splice_region": 4, "indel": 1},
    }

    def allocate_integer_quotas(total: int, weights: dict[str, int], available: dict[str, int]) -> dict[str, int]:
        active = {key: value for key, value in weights.items() if available.get(key, 0) > 0 and value > 0}
        if not active:
            return {}
        raw = {key: total * value / sum(active.values()) for key, value in active.items()}
        quotas = {key: min(int(raw[key]), available.get(key, 0)) for key in raw}
        while sum(quotas.values()) < total:
            candidates = [
                key for key in raw
                if quotas.get(key, 0) < available.get(key, 0)
            ]
            if not candidates:
                break
            key = max(candidates, key=lambda item: raw[item] - quotas.get(item, 0))
            quotas[key] = quotas.get(key, 0) + 1
        return quotas

    def select_training_ratio(quotas: dict[str, int], *, strict_primary_eligible: bool) -> pd.DataFrame:
        picked = []
        used = set()
        pool = arr[arr["primary_validation_eligible"]].copy() if strict_primary_eligible else arr.copy()
        pool["_is_model_vocab_consequence"] = pool["inferred_consequence_model"].isin(
            ["missense", "frameshift", "nonsense", "splice_region", "Missing", "indel"]
        )
        pool["_training_gene_weight"] = pool["gene"].map(training_gene_weights).fillna(0)
        pool["_has_hgvs_p"] = pool["hgvs_p_inferred"].map(clean_text).astype(bool)
        pool["_has_gnomad"] = (
            pool["gnomadE_AF_max_faf99"].map(numeric_or_blank).astype(bool)
            | pool["gnomadEG_AF"].map(numeric_or_blank).astype(bool)
        )
        pool["_training_ratio_priority"] = (
            pool["_is_model_vocab_consequence"].astype(int) * 10000
            + pool["_has_hgvs_p"].astype(int) * 1000
            + pool["_has_gnomad"].astype(int) * 500
            + pool["in_file_s4_functional_ep_dataset"].astype(int) * 100
            + pool["_training_gene_weight"] * 10
            + pool["_cadd_numeric"].clip(lower=0, upper=50)
            + pool["_n_capped"]
        )
        for label, quota in quotas.items():
            label_pool = pool[pool["reference_label_3class"].eq(label)].copy()
            available = label_pool["inferred_consequence_model"].value_counts().to_dict()
            consequence_quotas = allocate_integer_quotas(
                quota,
                training_consequence_weights.get(label, {"missense": 1}),
                available,
            )
            current = 0
            for consequence, consequence_quota in consequence_quotas.items():
                part = label_pool[label_pool["inferred_consequence_model"].eq(consequence)].sort_values(
                    ["_training_ratio_priority", "_cadd_numeric", "_n_numeric"],
                    ascending=False,
                )
                added = 0
                for _, row in part.iterrows():
                    key = row["variant_key_external"]
                    if key in used:
                        continue
                    picked.append(row)
                    used.add(key)
                    added += 1
                    current += 1
                    if added >= consequence_quota:
                        break
            if current < quota:
                remainder = label_pool.sort_values(
                    ["_training_ratio_priority", "_cadd_numeric", "_n_numeric"],
                    ascending=False,
                )
                for _, row in remainder.iterrows():
                    key = row["variant_key_external"]
                    if key in used:
                        continue
                    picked.append(row)
                    used.add(key)
                    current += 1
                    if current >= quota:
                        break
            if current < quota:
                raise ValueError(f"Could only select {current}/{quota} rows for {label}")
        out = pd.DataFrame(picked)
        order = {"Pathogenic": 0, "VUS": 1, "Benign": 2}
        out["_label_order"] = out["reference_label_3class"].map(order)
        out = out.sort_values(
            ["_label_order", "inferred_consequence_model", "_training_ratio_priority"],
            ascending=[True, True, False],
        )
        drop_cols = [c for c in out.columns if c.startswith("_")]
        return out.drop(columns=drop_cols).reset_index(drop=True)

    # Internal May 23 label ratio: 220 Benign, 207 VUS, 55 Pathogenic.
    # Benign availability limits the strict primary-eligible ratio-matched set.
    candidate_trainratio_strict = select_training_ratio(
        {"Benign": 119, "VUS": 112, "Pathogenic": 30},
        strict_primary_eligible=True,
    )
    candidate_trainratio_300_relaxed = select_training_ratio(
        {"Benign": 137, "VUS": 129, "Pathogenic": 34},
        strict_primary_eligible=False,
    )

    arr_out = arr.drop(columns=["_n_numeric", "_n_capped", "_cadd_numeric", "_priority"])
    s2.to_csv(PROCESSED / "emerge_s2_all_variants_extracted.csv", index=False)
    arr_out.to_csv(PROCESSED / "emerge_arrhythmia_variant_level_full.csv", index=False)
    candidate.to_csv(PROCESSED / "emerge_candidate_50.csv", index=False)
    candidate100.to_csv(PROCESSED / "emerge_candidate_100.csv", index=False)
    candidate100_matched.to_csv(PROCESSED / "emerge_candidate_100_training_matched.csv", index=False)
    candidate300_matched.to_csv(PROCESSED / "emerge_candidate_300_training_matched.csv", index=False)
    candidate_trainratio_strict.to_csv(PROCESSED / "emerge_candidate_trainratio_strict.csv", index=False)
    candidate_trainratio_300_relaxed.to_csv(PROCESSED / "emerge_candidate_trainratio_300_relaxed.csv", index=False)
    build_model_ready(arr_out, enriched=False).to_csv(
        PROCESSED / "emerge_arrhythmia_ml_baseline_features.csv", index=False
    )
    build_model_ready(arr_out, enriched=True).to_csv(
        PROCESSED / "emerge_arrhythmia_ml_enriched_features.csv", index=False
    )
    build_model_ready(candidate, enriched=False).to_csv(
        PROCESSED / "emerge_candidate_50_ml_baseline_features.csv", index=False
    )
    build_model_ready(candidate, enriched=True).to_csv(
        PROCESSED / "emerge_candidate_50_ml_enriched_features.csv", index=False
    )
    build_model_ready(candidate100, enriched=False).to_csv(
        PROCESSED / "emerge_candidate_100_ml_baseline_features.csv", index=False
    )
    build_model_ready(candidate100, enriched=True).to_csv(
        PROCESSED / "emerge_candidate_100_ml_enriched_features.csv", index=False
    )
    build_model_ready(candidate100_matched, enriched=False).to_csv(
        PROCESSED / "emerge_candidate_100_training_matched_ml_baseline_features.csv", index=False
    )
    build_model_ready(candidate100_matched, enriched=True).to_csv(
        PROCESSED / "emerge_candidate_100_training_matched_ml_enriched_features.csv", index=False
    )
    build_model_ready(candidate300_matched, enriched=False).to_csv(
        PROCESSED / "emerge_candidate_300_training_matched_ml_baseline_features.csv", index=False
    )
    build_model_ready(candidate300_matched, enriched=True).to_csv(
        PROCESSED / "emerge_candidate_300_training_matched_ml_enriched_features.csv", index=False
    )
    build_model_ready(candidate_trainratio_strict, enriched=False).to_csv(
        PROCESSED / "emerge_candidate_trainratio_strict_ml_baseline_features.csv", index=False
    )
    build_model_ready(candidate_trainratio_strict, enriched=True).to_csv(
        PROCESSED / "emerge_candidate_trainratio_strict_ml_enriched_features.csv", index=False
    )
    build_model_ready(candidate_trainratio_300_relaxed, enriched=False).to_csv(
        PROCESSED / "emerge_candidate_trainratio_300_relaxed_ml_baseline_features.csv", index=False
    )
    build_model_ready(candidate_trainratio_300_relaxed, enriched=True).to_csv(
        PROCESSED / "emerge_candidate_trainratio_300_relaxed_ml_enriched_features.csv", index=False
    )

    summary = {
        "s2_rows": int(len(s2)),
        "arrhythmia_rows": int(len(arr_out)),
        "candidate_rows": int(len(candidate)),
        "candidate100_rows": int(len(candidate100)),
        "candidate100_training_matched_rows": int(len(candidate100_matched)),
        "candidate300_training_matched_rows": int(len(candidate300_matched)),
        "candidate_trainratio_strict_rows": int(len(candidate_trainratio_strict)),
        "candidate_trainratio_300_relaxed_rows": int(len(candidate_trainratio_300_relaxed)),
        "all_gene_counts": s2["gene"].value_counts().to_dict(),
        "arrhythmia_label_counts": arr_out["reference_label_3class"].value_counts().to_dict(),
        "arrhythmia_eligible_label_counts": arr_out[arr_out["primary_validation_eligible"]][
            "reference_label_3class"
        ].value_counts().to_dict(),
        "candidate_label_counts": candidate["reference_label_3class"].value_counts().to_dict(),
        "candidate100_label_counts": candidate100["reference_label_3class"].value_counts().to_dict(),
        "candidate100_training_matched_label_counts": candidate100_matched["reference_label_3class"].value_counts().to_dict(),
        "candidate100_training_matched_gene_counts": candidate100_matched["gene"].value_counts().to_dict(),
        "candidate100_training_matched_consequence_counts": candidate100_matched["inferred_consequence_model"].value_counts().to_dict(),
        "candidate300_training_matched_label_counts": candidate300_matched["reference_label_3class"].value_counts().to_dict(),
        "candidate300_training_matched_gene_counts": candidate300_matched["gene"].value_counts().to_dict(),
        "candidate300_training_matched_consequence_counts": candidate300_matched["inferred_consequence_model"].value_counts().to_dict(),
        "candidate_trainratio_strict_label_counts": candidate_trainratio_strict["reference_label_3class"].value_counts().to_dict(),
        "candidate_trainratio_strict_gene_counts": candidate_trainratio_strict["gene"].value_counts().to_dict(),
        "candidate_trainratio_strict_consequence_counts": candidate_trainratio_strict["inferred_consequence_model"].value_counts().to_dict(),
        "candidate_trainratio_300_relaxed_label_counts": candidate_trainratio_300_relaxed["reference_label_3class"].value_counts().to_dict(),
        "candidate_trainratio_300_relaxed_gene_counts": candidate_trainratio_300_relaxed["gene"].value_counts().to_dict(),
        "candidate_trainratio_300_relaxed_consequence_counts": candidate_trainratio_300_relaxed["inferred_consequence_model"].value_counts().to_dict(),
        "file_s3_disagreement_rows_in_arrhythmia": int(arr_out["in_file_s3_sequencing_center_disagreement"].sum()),
        "file_s4_functional_rows_in_arrhythmia": int(arr_out["in_file_s4_functional_ep_dataset"].sum()),
        "notes": {
            "primary_label": "finalAnnotation2_postInvitro collapsed to Benign/VUS/Pathogenic",
            "phenotype": "No row-level phenotype values were found in File S2; File S1 is an EHR code dictionary.",
            "candidate_set": "Initial balanced 50-row set for review/scoring; not yet final for manuscript.",
        },
    }
    (RESULTS / "emerge_extraction_summary.json").write_text(json.dumps(summary, indent=2))

    rows = [
        {"file": "File S2", "column": col, "description": ""}
        for col in s2.columns
    ]
    pd.DataFrame(rows).to_csv(PROCESSED / "emerge_data_dictionary_stub.csv", index=False)


if __name__ == "__main__":
    build_tables()
