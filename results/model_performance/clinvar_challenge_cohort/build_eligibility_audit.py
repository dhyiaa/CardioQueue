#!/usr/bin/env python3
"""Audit candidates for a separate high-confidence ClinVar challenge cohort.

This script deliberately does not alter the frozen source-held-out split. ClinVar
rows remain same-ecosystem observations even when they are removed before a new
fit, so the output is an eligibility manifest rather than an external cohort.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
MATRIX = ROOT / "datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv"
CLINVAR = ROOT / "datasets/clinvar/interim/clinvar_cardiogenetics_modeling_slice.tsv"
CLINGEN_EREPO = ROOT / (
    "datasets/feature_sources/clingen/interim/variant_pathogenicity.master_registry_panel.tsv"
)
OUT_DIR = ROOT / "results/model_performance/clinvar_challenge_cohort"

CARDIAC_VCEPS = {"Cardiomyopathy VCEP", "Potassium Channel Arrhythmia VCEP"}

CARDIAC_RE = re.compile(
    r"cardio|cardiac|arrhythm|cardiomy|long qt|brugada|catecholaminergic|"
    r"aort|marfan|loeys|dilated|hypertrophic|ventricular|atrial|heart|"
    r"sudden death|channelopathy|short qt",
    re.IGNORECASE,
)
RASOPATHY_RE = re.compile(
    r"noonan|rasopath|leopard|costello|cardiofaciocutaneous", re.IGNORECASE
)


def stable_hash(value: str) -> str:
    return hashlib.sha256(f"cardioqueue-clinvar-challenge-v1|{value}".encode()).hexdigest()


def counts(frame: pd.DataFrame, columns: list[str]) -> list[dict[str, object]]:
    return frame.groupby(columns, dropna=False).size().rename("n").reset_index().to_dict("records")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    matrix = pd.read_csv(MATRIX, sep="\t", low_memory=False)
    clinvar = pd.read_csv(CLINVAR, sep="\t", low_memory=False)
    erepo = pd.read_csv(CLINGEN_EREPO, sep="\t", low_memory=False)

    matrix_fields = [
        "variant_id",
        "primary_gene",
        "model_label_3class",
        "clinvar_review_stars",
        "split_source_heldout",
        "primary_model_inclusion",
        "label_conflict_type",
        "vep_worst_consequence",
        "in_hiro",
        "in_emerge",
        "in_cardioboost",
    ]
    clinvar_fields = [
        "variant_id",
        "phenotype_list",
        "phenotype_ids",
        "review_status",
        "number_submitters",
        "last_evaluated",
        "variation_id",
        "rcv_accession",
    ]
    frame = matrix[matrix_fields].merge(
        clinvar[clinvar_fields], on="variant_id", how="left", validate="one_to_one"
    )

    erepo = erepo[erepo["expert_panel"].isin(CARDIAC_VCEPS)].copy()
    erepo["cardiac_vcep_label"] = erepo["classification"].map(
        {
            "Benign": "Benign",
            "Likely Benign": "Benign",
            "Pathogenic": "Pathogenic",
            "Likely Pathogenic": "Pathogenic",
        }
    )
    erepo = erepo[erepo["cardiac_vcep_label"].notna()].copy()
    erepo["variation_id"] = erepo["clinvar_variation_id"].astype("string")
    erepo_by_variant = (
        erepo.groupby("variation_id", dropna=False)
        .agg(
            cardiac_vcep_panel=("expert_panel", lambda s: "|".join(sorted(set(s.astype(str))))),
            cardiac_vcep_condition=("condition", lambda s: "|".join(sorted(set(s.dropna().astype(str))))),
            cardiac_vcep_label=("cardiac_vcep_label", lambda s: "|".join(sorted(set(s.astype(str))))),
            cardiac_vcep_approved_date=("approved_date", "max"),
        )
        .reset_index()
    )
    frame["variation_id"] = frame["variation_id"].astype("string")
    frame = frame.merge(erepo_by_variant, on="variation_id", how="left", validate="many_to_one")
    stars = pd.to_numeric(frame["clinvar_review_stars"], errors="coerce")
    phenotype = frame["phenotype_list"].fillna("")
    source_overlap = (
        frame[["in_hiro", "in_emerge", "in_cardioboost"]]
        .astype("string")
        .apply(lambda col: col.str.lower().isin(["true", "1", "yes", "y", "t"]))
        .any(axis=1)
    )

    frame["cardiac_phenotype_keyword"] = phenotype.str.contains(CARDIAC_RE)
    frame["rasopathy_keyword"] = phenotype.str.contains(RASOPATHY_RE)
    frame["expert_panel_status"] = stars.ge(3)
    frame["cardiac_vcep_verified"] = (
        frame["cardiac_vcep_panel"].notna()
        & frame["cardiac_vcep_label"].eq(frame["model_label_3class"])
    )
    frame["same_ecosystem_clinvar"] = True
    frame["stable_sampling_hash"] = frame["variant_id"].map(stable_hash)

    eligible = (
        frame["split_source_heldout"].isin(["train", "validation"])
        & frame["model_label_3class"].isin(["Benign", "Pathogenic"])
        & frame["primary_model_inclusion"].eq("include")
        & frame["label_conflict_type"].fillna("none").eq("none")
        & stars.ge(2)
        & (
            frame["cardiac_vcep_verified"]
            | (frame["cardiac_phenotype_keyword"] & ~frame["rasopathy_keyword"])
        )
        & ~source_overlap
    )
    candidates = frame.loc[eligible].copy()
    candidates["clinvar_review_stars"] = pd.to_numeric(
        candidates["clinvar_review_stars"], errors="coerce"
    ).astype("Int64")
    candidates["proposed_role"] = "condition-filtered challenge candidate"
    candidates.loc[candidates["expert_panel_status"], "proposed_role"] = (
        "unverified-specialty expert-panel candidate"
    )
    candidates.loc[candidates["cardiac_vcep_verified"], "proposed_role"] = (
        "cardiac-VCEP-verified seed"
    )
    candidates = candidates.sort_values(
        ["clinvar_review_stars", "model_label_3class", "primary_gene", "stable_sampling_hash"],
        ascending=[False, True, True, True],
    )

    manifest = OUT_DIR / "clinvar_cardiac_challenge_eligible.tsv"
    candidates.to_csv(manifest, sep="\t", index=False)

    expert = candidates[candidates["expert_panel_status"]]
    cardiac_vcep = candidates[candidates["cardiac_vcep_verified"]]
    two_star = candidates[~candidates["expert_panel_status"]]
    existing_external = matrix[
        matrix["split_source_heldout"].astype("string").str.startswith("external_")
        & matrix["model_label_3class"].isin(["Benign", "Pathogenic"])
    ]
    existing_labels = existing_external["model_label_3class"].value_counts().to_dict()

    # A 570-row addition with 245 P/LP and 325 B/LB would make the combined
    # descriptive benchmark exactly balanced (500/500) given the current 430.
    target = {"Benign": 325, "Pathogenic": 245}
    expert_labels = expert["model_label_3class"].value_counts().to_dict()
    cardiac_vcep_labels = cardiac_vcep["model_label_3class"].value_counts().to_dict()
    proposed_remainder = {
        label: target[label] - cardiac_vcep_labels.get(label, 0) for label in target
    }

    summary = {
        "purpose": "Eligibility audit for a separate same-ecosystem ClinVar challenge cohort",
        "frozen_primary_external_rows": int(len(existing_external)),
        "frozen_primary_external_label_counts": existing_labels,
        "eligible_rows": int(len(candidates)),
        "eligible_expert_panel_rows": int(len(expert)),
        "eligible_cardiac_vcep_verified_rows": int(len(cardiac_vcep)),
        "eligible_two_star_rows": int(len(two_star)),
        "eligible_by_star_and_label": counts(
            candidates, ["clinvar_review_stars", "model_label_3class"]
        ),
        "expert_panel_by_gene": counts(expert, ["primary_gene"]),
        "cardiac_vcep_verified_by_panel_and_label": counts(
            cardiac_vcep, ["cardiac_vcep_panel", "model_label_3class"]
        ),
        "cardiac_vcep_verified_by_gene": counts(cardiac_vcep, ["primary_gene"]),
        "provisional_570_target": target,
        "provisional_condition_filtered_remainder_after_cardiac_vcep_seed": proposed_remainder,
        "combined_1000_label_counts_if_target_met": {
            label: int(existing_labels.get(label, 0) + target[label]) for label in target
        },
        "claim_boundaries": [
            "The 430 HiRO/eMERGE/CardioBoost rows remain the primary source-held-out evaluation.",
            "ClinVar candidates are a separate same-ecosystem challenge set, not external validation.",
            "Cardiac-VCEP verification comes from the local ClinGen Evidence Repository panel identity and a concordant current binary label.",
            "Other three-star records have expert-panel review status but not locally verified cardiogenetics-panel identity; two-star records are not expert-panel records.",
            "Cardiac phenotype keyword matching is an audit screen and does not identify the reviewing expert panel affiliation.",
            "No candidate should be frozen until condition and expert-panel affiliation are checked at the assertion level.",
        ],
    }
    with (OUT_DIR / "clinvar_cardiac_challenge_eligibility.summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)

    report = f"""# ClinVar Cardiac Challenge Eligibility Audit

This audit found **{len(candidates):,}** binary, at-least-two-star ClinVar rows with either direct cardiac-VCEP verification or a cardiac phenotype keyword without an obvious RASopathy keyword, no HiRO/eMERGE/CardioBoost overlap, and an existing development assignment.

- Three-star expert-panel rows: **{len(expert):,}** ({expert_labels.get('Benign', 0)} B/LB; {expert_labels.get('Pathogenic', 0)} P/LP)
- Cardiac-VCEP-verified rows with concordant current labels: **{len(cardiac_vcep):,}** ({cardiac_vcep_labels.get('Benign', 0)} B/LB; {cardiac_vcep_labels.get('Pathogenic', 0)} P/LP)
- Two-star multiple-submitter/no-conflict rows: **{len(two_star):,}**
- Existing source-held-out cohort: **{len(existing_external):,}** ({existing_labels.get('Benign', 0)} B/LB; {existing_labels.get('Pathogenic', 0)} P/LP)

## Defensible design

Keep the 430-row HiRO/eMERGE/CardioBoost result as the primary source-held-out analysis. Build a separate 570-row ClinVar challenge cohort only after assertion-level disease-context checks. Use the locally verified cardiac-VCEP rows as a seed, then draw the remainder from condition-filtered high-confidence records using a prespecified gene-by-consequence design. A provisional target of 325 B/LB and 245 P/LP would produce a descriptive pooled benchmark of exactly 500 B/LB and 500 P/LP, but the two cohorts must also be reported separately.

## Important boundary

ClinVar review stars describe review status, not clinical specialty. The phenotype keyword filter improves relevance but cannot prove that the reviewing panel was a cardiogenetics panel. The manifest therefore remains an eligibility audit and does not change model splits.
"""
    (OUT_DIR / "CLINVAR_CARDIAC_CHALLENGE_AUDIT.md").write_text(report)


if __name__ == "__main__":
    main()
