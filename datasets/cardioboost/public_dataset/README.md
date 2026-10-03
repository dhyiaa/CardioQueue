# CardioBoost Public Dataset Handoff for Cardiogenetics

This folder packages the public CardioBoost paper datasets and the local Cardiogenetics-ready processed versions. It is intended for a future Codex agent working on the new Cardiogenetics study proposal.

Important: CardioBoost is mostly a binary variant-level dataset, not a three-class Cardiogenetics cohort. It has `Benign` and `Pathogenic` labels, but no VUS class in the processed Cardiogenetics-ready holdout set. It also does not contain the same patient-level phenotype fields as the internal Cardiogenetics/HiRO-CASPER-VERDICT dataset.

## Source

Original public CardioBoost repository:

Zhang X, Walsh R, Whiffin N, et al. Disease-specific variant pathogenicity prediction significantly improves variant interpretation in inherited cardiac conditions. Genetics in Medicine. 2021;23(1):69-79. doi:10.1038/s41436-020-00972-3

The original public repo README is copied here:

- `raw_public_cardioBoost_repo/CardioBoost_original_README.md`

## Folder Structure

### `raw_public_cardioBoost_repo/`

This is the raw public CardioBoost data folder copied from the CardioBoost manuscript repository.

It contains `.RData` files, gene lists, and preprocessing artifacts. Old model
training and prediction artifact folders were removed because this project uses
CardioBoost as a variant data source, not as an imported model.

Main subfolders:

- `raw_public_cardioBoost_repo/data/arrhythmia/`
  - Original CardioBoost arrhythmia dataset files.

- `raw_public_cardioBoost_repo/data/cardiomyopathy/`
  - Original CardioBoost cardiomyopathy dataset files.

- `raw_public_cardioBoost_repo/data/cardiomyopathy_alternative/`
  - Alternative cardiomyopathy disease/gene-specific datasets from the CardioBoost repo.

### `processed_for_cardiogenetics/`

This contains local processed versions made for the Cardiogenetics project. These are easier to use from Python/CatBoost-style workflows.

Subfolders:

- `model_inputs/`
  - CSV files that can be inspected or used as model inputs.

- `provenance/`
  - Build summaries, duplicate audits, schema comparison, and manifests.

- `results/`
  - Removed from this handoff because those files were old weighted-CatBoost
    scoring outputs, not source variant data.

## Raw Public CardioBoost Dataset Contents

The raw `.RData` objects inspected in this project include:

### Arrhythmia

| File | Object | Rows | Columns | Notes |
| --- | --- | ---: | ---: | --- |
| `arrhythmia/arm_raw_training.RData` | `train` | 532 | 85 | Raw arrhythmia training table |
| `arrhythmia/arm_train.RData` | `train` | 308 | 85 | CardioBoost arrhythmia training table |
| `arrhythmia/arm_holdout_test.RData` | `test` | 154 | 85 | CardioBoost arrhythmia holdout test table |
| `arrhythmia/arm_additional_benign_test.RData` | `arm_gnomad_rare_noexac` | 1,237 | 17 | Additional rare benign-like variants |
| `arrhythmia/arm_additional_patho_test.RData` | `arm_patho_test` | 215 | 17 | Additional pathogenic variants |

Main arrhythmia genes in the holdout file include `SCN5A`, `KCNH2`, `KCNQ1`, `CACNA1C`, `CALM2`, and `CALM1`.

### Cardiomyopathy

| File | Object | Rows | Columns | Notes |
| --- | --- | ---: | ---: | --- |
| `cardiomyopathy/cm_raw_training.RData` | `train` | 2,578 | 85 | Raw cardiomyopathy training table |
| `cardiomyopathy/cm_train.RData` | `train` | 440 | 85 | CardioBoost cardiomyopathy training table |
| `cardiomyopathy/cm_holdout_test.RData` | `test` | 218 | 85 | CardioBoost cardiomyopathy holdout test table |
| `cardiomyopathy/cm_additional_test_benign.RData` | `cm_gnomad_rare_noexac` | 2,003 | 16 | Additional rare benign-like variants |
| `cardiomyopathy/cm_additional_test_pathogenic.RData` | `cm_patho_test` | 289 | 11 | Additional pathogenic variants |

Main cardiomyopathy genes in the holdout file include `MYH7`, `MYBPC3`, `SCN5A`, `LMNA`, `TNNI3`, `TNNT2`, `PRKAG2`, and `TPM1`.

## Processed Cardiogenetics-Ready Files

The most useful files for the new Cardiogenetics project are in:

- `processed_for_cardiogenetics/model_inputs/`

### Recommended starting files

Use these first if the goal is to train/evaluate a Cardiogenetics-style model:

- `reannotated_ml_baseline_features.csv`
  - 355 rows x 89 columns.
  - Deduplicated CardioBoost holdout-overlap data reannotated into the Cardiogenetics baseline feature schema.
  - Binary labels only: 199 Pathogenic, 156 Benign.

- `reannotated_ml_enriched_features.csv`
  - 355 rows x 95 columns.
  - Same rows as the baseline file, with enriched features retained.
  - Binary labels only: 199 Pathogenic, 156 Benign.

- `reannotated_modeling_cohort.csv`
  - 355 rows x 150 columns.
  - Richer intermediate cohort table before narrowing to model feature columns.

### Other processed files

- `cardioboost_overlap_raw_rows.csv`
  - 361 row-level CardioBoost holdout rows after filtering to genes overlapping the local Cardiogenetics training gene set.

- `cardioboost_holdout_overlap_rowlevel_ml_baseline_features.csv`
  - 361 row-level feature rows.
  - Keeps exact duplicate variants that appeared in more than one CardioBoost panel.

- `cardioboost_holdout_overlap_deduplicated_ml_baseline_features.csv`
  - 355 deduplicated feature rows.
  - Older processed baseline export with CardioBoost-specific helper columns.

## Processed Dataset Counts

The clean reannotated processed set has:

| Field | Count |
| --- | ---: |
| Total deduplicated variants | 355 |
| Cardiomyopathy variants | 202 |
| Arrhythmia variants | 153 |
| Pathogenic | 199 |
| Benign | 156 |
| VUS | 0 |

Gene counts in the 355-row processed set:

| Gene | Count |
| --- | ---: |
| MYH7 | 77 |
| SCN5A | 68 |
| KCNH2 | 41 |
| MYBPC3 | 41 |
| KCNQ1 | 37 |
| CACNA1C | 22 |
| LMNA | 12 |
| TNNI3 | 12 |
| TNNT2 | 10 |
| PRKAG2 | 9 |
| TPM1 | 9 |
| CALM2 | 6 |
| GLA | 6 |
| PTPN11 | 3 |
| ACTC1 | 1 |
| PLN | 1 |

## Why This Dataset Is Binary

The processed CardioBoost holdout-overlap set uses labels compatible with the original CardioBoost pathogenicity task:

- `Benign`
- `Pathogenic`

There are no VUS rows in the processed 355-row Cardiogenetics-ready set. This differs from the internal Cardiogenetics cohort and the eMERGE full dataset, which both include VUS/uncertain labels.

For binary modeling, the processed feature files use:

- `target_3class = Benign` or `Pathogenic`
- `target_numeric = 0` for Benign
- `target_numeric = 2` for Pathogenic

The numeric value `1`, used elsewhere for VUS, is absent.

## How To Use This In A New Model

Good uses:

- External validation of a Cardiogenetics model.
- Public training signal for a binary pathogenic-vs-benign classifier.
- Public pretraining before fine-tuning/evaluating on internal Cardiogenetics data.
- A sensitivity analysis asking whether public variant-only data improve pathogenic detection.
- A separate benchmark next to eMERGE.

Risky uses:

- Do not merge directly with a three-class Benign/VUS/Pathogenic training cohort unless the modeling plan explicitly handles missing VUS labels.
- Do not treat the missing phenotype fields as true negatives or normal phenotype values.
- Do not compare performance against eMERGE without noting that CardioBoost is binary and eMERGE is heavily VUS-skewed.
- Do not use CardioBoost rows as both training and validation in the same experiment.

## Phenotype Caveat

The CardioBoost public data are variant-level pathogenicity data. They do not have the same patient-level phenotype structure used in the internal Cardiogenetics cohort.

In the Cardiogenetics-ready model files, phenotype-derived fields are therefore missing or set to `Missing`. This is intentional.

## Provenance

Useful provenance files:

- `processed_for_cardiogenetics/provenance/build_summary.json`
  - Documents the 372 combined CardioBoost holdout rows, filtering to 361 rows after local training-gene overlap, and deduplication to 355 variants.

- `processed_for_cardiogenetics/provenance/reannotated_model_input_manifest.json`
  - Documents annotation coverage, label balance, missing phenotype fields, and feature schema.

- `processed_for_cardiogenetics/provenance/reannotated_run_manifest.json`
  - Documents the reannotation run and output locations.

- `processed_for_cardiogenetics/provenance/cardioboost_exact_duplicate_audit.csv`
  - Documents variants duplicated across panels before deduplication.

## Old Results Removed

The folder originally included prior weighted CatBoost scoring outputs under:

- `processed_for_cardiogenetics/results/`

Those files were deleted for this project because they were previous model
predictions and metric summaries, not source variant data. The processed model
input CSVs and provenance files are retained.

## Suggested Rule For The Next Codex Agent

If the new study needs a public external dataset, start with:

`processed_for_cardiogenetics/model_inputs/reannotated_ml_baseline_features.csv`

If the new study needs explicit enriched annotation features, use:

`processed_for_cardiogenetics/model_inputs/reannotated_ml_enriched_features.csv`

If the new study needs original CardioBoost reproduction, inspect the raw `.RData` files under:

`raw_public_cardioBoost_repo/data/`

Keep CardioBoost separate from eMERGE in analysis tables because CardioBoost is binary and comparatively balanced, while eMERGE is three-class and VUS-heavy.
