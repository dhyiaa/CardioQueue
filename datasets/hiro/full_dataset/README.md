# HiRO / CASPER WES / VERDICT Full Dataset Handoff

This folder packages the April 2026 HiRO-derived Cardiogenetics dataset for a new project/Codex agent. It includes the raw source exports, normalized cohort files, annotation outputs, model-ready feature tables, exclusion records, and provenance scripts.

Important privacy note: the raw patient files may contain clinical export fields and dates. Keep this folder in approved local/private analysis environments only. The public GitHub snapshot intentionally excluded these files.

## What This Dataset Is

This is the internal April 2026 CASPER WES / VERDICT Cardiogenetics cohort.

Raw source scope:

- 512 classified variant rows
- 234 patient rows
- 2 studies:
  - CASPER WES
  - VERDICT

Primary model-ready scope:

- 482 variant-to-phenotype rows
- Three-class labels: `Benign`, `VUS`, `Pathogenic`
- Patient phenotype rows are matched and joined
- GC/research-team ACMG criteria are retained for reference/audit but excluded from primary ML features

## Folder Structure

### `raw_hiro_april2026/`

Original source exports:

- `ACMG AI Project_Classified Variants_07April2026(CASPER WES).csv`
  - 172 CASPER WES classified variant rows.

- `ACMG AI Project_Classified Variants_07April2026(VERDICT).csv`
  - 340 VERDICT classified variant rows.

- `ACMG_CASPER_Wes_Patient_Data_Apr_2026(in).csv`
  - 114 CASPER WES patient rows.

- `ACMG_Verdict_Patient_Data_Apr_2026(in).csv`
  - 120 VERDICT patient rows.

- `Hiro_Full_Variable_List.csv`
  - HiRO variable dictionary used to decode/interpret coded patient fields.

### `processed_normalized/`

Pipeline-normalized data:

- `normalized_variants.csv`
  - 512 rows x 20 columns.
  - One row per classified variant from CASPER WES / VERDICT.

- `normalized_patients.csv`
  - 234 rows x 67 columns.
  - One row per patient with derived phenotype summaries.
  - DOB and raw clinical dates were not exported into this normalized file.

- `analysis_cohort.csv`
  - 512 rows x 86 columns.
  - Variant rows joined to patient features where available.
  - Includes `patient_match` flag.

- `variable_dictionary_used.csv`
  - Subset of HiRO dictionary rows used by the parser.

- `preprocessing_report.json`
  - Counts, class balance, unmatched variants, duplicate rows, and criteria summaries.

### `annotations/`

External annotation outputs:

- `variant_annotations.csv`
  - 510 rows x 46 columns.
  - ClinVar, gnomAD, MyVariant/dbNSFP annotation output for preferred labeled variant rows.

- `annotated_cohort.csv`
  - 510 rows x 127 columns.
  - Preferred labeled cohort joined to external annotations.

- `annotation_cache.jsonl`
  - Cached annotation API results.

- `annotation_report.json`
  - Annotation status counts.

### `model_inputs/`

Model-ready and reference files:

- `variant_only_cohort.csv`
  - 510 rows x 127 columns.
  - Preferred labeled rows regardless of patient phenotype match.
  - Use this if the new study wants all labeled variants, including rows without matched phenotype.

- `modeling_cohort.csv`
  - 482 rows x 127 columns.
  - Preferred labeled rows with matched phenotype data.
  - This was the main cohort used for primary ML/LLM modeling.

- `ml_baseline_features.csv`
  - 482 rows x 89 columns.
  - Primary ML feature table without explicit ClinVar consensus/classification fields.

- `ml_enriched_features.csv`
  - 482 rows x 95 columns.
  - ML feature table with enriched ClinVar metadata/classification fields.

- `criterion_reference.csv`
  - 510 rows x 100 columns.
  - GC-applied ACMG criteria expanded into reference columns.
  - These are not primary ML inputs; they are reference/evaluation fields.

- `excluded_rows.csv`
  - 28 rows x 8 columns.
  - Preferred labeled rows excluded from primary modeling because no matching phenotype row was available.

- `model_input_manifest.json`
  - Main manifest for row counts, class balance, missingness, annotation coverage, feature lists, and excluded fields.

### `provenance/`

Scripts and docs copied for reproducibility:

- `provenance/scripts/build_new_cohort.py`
- `provenance/scripts/validate_new_cohort.py`
- `provenance/scripts/annotate_variants.py`
- `provenance/scripts/validate_annotations.py`
- `provenance/scripts/build_model_inputs.py`
- `provenance/scripts/validate_model_inputs.py`
- `provenance/docs/new_pipeline.md`
- `provenance/docs/project_context.md`

## Row Flow

| Step | Rows | Notes |
| --- | ---: | --- |
| Raw classified variants | 512 | 172 CASPER WES + 340 VERDICT |
| Normalized variants | 512 | Includes all source variant rows |
| Analysis cohort | 512 | Joined to patient features when available |
| Preferred labeled variant-only cohort | 510 | Drops one missing-label row and resolves exact duplicate handling |
| Primary modeling cohort | 482 | Keeps preferred labeled rows with matched phenotype data |
| Excluded from primary modeling | 28 | Labeled variants without matching phenotype row |

## Label Balance

### Full normalized variants, 512 rows

| Class | Count |
| --- | ---: |
| VUS | 227 |
| Benign | 225 |
| Pathogenic | 59 |
| Missing | 1 |

### Preferred labeled variant-only cohort, 510 rows

| Class | Count |
| --- | ---: |
| VUS | 226 |
| Benign | 225 |
| Pathogenic | 59 |

### Primary modeling cohort, 482 rows

| Class | Count |
| --- | ---: |
| Benign | 220 |
| VUS | 207 |
| Pathogenic | 55 |

## Study Balance

Primary modeling cohort:

| Study | Count |
| --- | ---: |
| VERDICT | 318 |
| CASPER WES | 164 |

Variant-only cohort:

| Study | Count |
| --- | ---: |
| VERDICT | 340 |
| CASPER WES | 170 |

## Exclusions And Issues

### Missing classification

One CASPER WES row had no final classification:

- `CASPER_WES|5395|PKP2|C.G473A:P.R158K|172`

It appears in `normalized_variants.csv` as `classification_3class = Missing`, but it is not included in the preferred labeled model-input cohort.

### Duplicate variant row

Two CASPER WES rows had the same participant/gene/variant key:

- `2085|RYR2|C.G9352A:P.G3118R`

The duplicate rows had discordant labels in the source-normalized table:

- Pathogenic
- VUS

Duplicate handling is already reflected in the preferred model-input layer.

### No matching phenotype row

There were 29 unmatched variant rows in the full 512-row analysis cohort:

- CASPER WES: 7
- VERDICT: 22

After removing the missing-label/non-preferred row, 28 preferred labeled rows remained unmatched and were excluded from the primary 482-row modeling cohort. These rows are saved in:

- `model_inputs/excluded_rows.csv`

Excluded labeled rows by class:

| Class | Count |
| --- | ---: |
| VUS | 19 |
| Benign | 5 |
| Pathogenic | 4 |

All 28 exclusions have:

- `exclusion_reasons = no_matching_phenotype_row`

## Missingness Notes

The main missingness summary is in:

- `model_inputs/model_input_manifest.json`

Important primary-modeling missingness:

| Field | Missing / 482 | Missing fraction |
| --- | ---: | ---: |
| `patient_family_history_flags` | 482 | 1.000 |
| `patient_followup_circumstances_latest` | 482 | 1.000 |
| `patient_age_at_diagnosis` | 75 | 0.156 |
| `patient_qt_latest` | 16 | 0.033 |
| `patient_qt_max` | 16 | 0.033 |
| `patient_lvef_latest` | 63 | 0.131 |
| `patient_lvef_min` | 63 | 0.131 |
| `patient_sex` | 0 | 0.000 |
| `patient_sym_syncope` | 0 | 0.000 |
| `patient_sym_cardiac_arrest` | 0 | 0.000 |

In-silico score coverage in the 482-row primary modeling cohort:

| Score | Non-missing | Missing | Missing fraction |
| --- | ---: | ---: | ---: |
| REVEL | 316 | 166 | 0.344 |
| CADD Phred | 136 | 346 | 0.718 |
| SIFT | 313 | 169 | 0.351 |
| PolyPhen-2 HDIV | 298 | 184 | 0.382 |
| MetaLR | 317 | 165 | 0.342 |
| FATHMM-XF | 323 | 159 | 0.330 |
| AlphaMissense | 317 | 165 | 0.342 |

Annotation status in the 510-row preferred labeled layer:

| Status | Count |
| --- | ---: |
| partial | 202 |
| complete | 177 |
| no_external_match | 111 |
| partial_with_error | 20 |

## Features Excluded From Primary ML Inputs

The following fields were intentionally excluded from primary ML features because they are leakage-prone, criteria/reference fields, or sensitivity-only clinical interpretation fields:

- `criteria_raw`
- `criteria_codes`
- `criteria_uncertain`
- `criteria_strengths`
- `criteria_count`
- `patient_working_diagnosis_raw`
- `patient_working_diagnosis`
- `patient_diagnosis_strength_raw`
- `patient_diagnosis_strength`
- `patient_diagnosis_feature_role`

GC-applied ACMG criteria are retained in `criterion_reference.csv` for validation and analysis, but should not be used as primary model features unless the new study explicitly defines a criteria-aware sensitivity model.

## Which File Should The Next Codex Agent Use?

Use this for the primary phenotype-aware three-class model:

- `model_inputs/ml_baseline_features.csv`

Use this if explicit ClinVar consensus/classification metadata is allowed:

- `model_inputs/ml_enriched_features.csv`

Use this if the study needs the richer cohort table before feature narrowing:

- `model_inputs/modeling_cohort.csv`

Use this if the study wants all preferred labeled variants, including rows with no phenotype match:

- `model_inputs/variant_only_cohort.csv`

Use this if the study wants the complete source-normalized variant set, including missing-label and unmatched rows:

- `processed_normalized/normalized_variants.csv`

Use raw source files only when re-running or auditing preprocessing:

- `raw_hiro_april2026/`

## Integrity Notes

Selected SHA-256 checksums at handoff:

- `model_inputs/ml_baseline_features.csv`: `5d32a044abed79d697f3f2ddee97963507e65bfaa4d8f52ffe6f094a2610c558`
- `model_inputs/ml_enriched_features.csv`: `7894586b1fc9b342d1d6c09e38b8a3ab176b14331241c9bd59344a8a11bcedf3`
- `model_inputs/modeling_cohort.csv`: `91033246189542d730d671c97ede0e028355de1e81ef89a4de7dd3dd23a8d950`
- `model_inputs/variant_only_cohort.csv`: `972fd37483ad3cd9efd884acaac6c65d666e6b008fc4a407332f22fed021525a`
- `processed_normalized/normalized_variants.csv`: `7bc2779b0544db095f27b1353d4c5db7825991f91e50eee5c5f13e4d6ab43107`

