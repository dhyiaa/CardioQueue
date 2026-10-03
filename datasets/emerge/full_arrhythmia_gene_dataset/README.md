# eMERGE Full Arrhythmia-Gene Dataset for Cardiogenetics

This folder is a handoff package for using the full public eMERGE/Glazer arrhythmia-gene variant dataset as real variant-level data in the Cardiogenetics project.

The previous external validation study used a balanced 300-variant subset. This folder keeps the full arrhythmia-gene dataset instead: 2,754 variants from the public eMERGE supplement after excluding BRCA2 and retaining arrhythmia/cardiogenetic genes.

## What Is Included

### Main full dataset

- `data/emerge_arrhythmia_variant_level_full.csv`
  - Full variant-level eMERGE arrhythmia-gene table.
  - 2,754 rows and 47 columns.
  - This is the main file to use if the new study wants the complete external eMERGE dataset.

### Model-ready feature tables

- `data/emerge_arrhythmia_ml_baseline_features.csv`
  - Full 2,754-row feature table harmonized to the May 23 Cardiogenetics/CatBoost baseline schema.
  - Patient-level phenotype fields are intentionally missing because the public eMERGE supplement does not provide row-level clinical phenotype data.

- `data/emerge_arrhythmia_ml_enriched_features.csv`
  - Full 2,754-row feature table with the enriched feature columns retained.
  - Some in-silico/HGVS fields were filled from MyVariant/dbNSFP where available.

### Provenance files

- `provenance/build_emerge_variant_tables.py`
  - Script used to build the eMERGE variant-level and model-ready tables from the public supplements.

- `provenance/enrich_emerge_insilico.py`
  - Script used to enrich missing in-silico annotation fields from MyVariant/dbNSFP.

- `provenance/emerge_extraction_summary.json`
  - Row counts, gene counts, class counts, and candidate-set summaries from the extraction run.

- The previous weighted-CatBoost scoring summary was removed because it was an
  old model-result artifact, not source variant data.

- `provenance/emerge_data_dictionary_stub.csv`
  - Source-column stub created during extraction. It is useful as a starting point, but not a polished data dictionary.

## How This Dataset Was Made

The starting point was the public eMERGE/Glazer arrhythmia-gene supplementary data:

- `NIHMS1778916-supplement-File_S2_Variant_Summary.xlsx`
  - Primary variant table.

- `NIHMS1778916-supplement-File_S3_Variants_with_disagreements_between_sequencing_centers.xlsx`
  - Used to flag variants with sequencing-center disagreements.

- `NIHMS1778916-supplement-File_S4_In_vitro_EP_dataset.xlsx`
  - Used to flag variants present in the functional electrophysiology dataset.

The extraction script loaded all 3,550 variants from File S2, then retained the arrhythmia/cardiogenetic genes overlapping the model scope:

`ANK2`, `CACNA1C`, `KCNE1`, `KCNE2`, `KCNH2`, `KCNJ2`, `KCNQ1`, `LMNA`, `RYR2`, and `SCN5A`.

`BRCA2` was excluded because it is outside the arrhythmia/cardiogenetics model scope.

The final full arrhythmia-gene table contains 2,754 variants.

## Labels

The main three-class label is:

- `reference_label_3class`

It was created by collapsing `finalAnnotation2_postInvitro` into:

- `Benign`
- `VUS`
- `Pathogenic`

The pre-in-vitro label is also retained:

- `reference_label_3class_pre_invitro`

This matters because the 300-row validation manuscript primarily discussed pre-in-vitro labels for ACMG-aligned model comparison, while post-in-vitro labels were treated as a functional-evidence stress test.

Full post-in-vitro label counts in this file:

| Label | Count |
| --- | ---: |
| VUS | 2,490 |
| Benign | 137 |
| Pathogenic | 127 |

Pre-in-vitro label counts:

| Label | Count |
| --- | ---: |
| VUS | 2,501 |
| Benign | 134 |
| Pathogenic | 119 |

## Eligibility Flags

The file includes flags that help reproduce the earlier validation logic:

- `is_arrhythmia_gene`
  - True for all rows in this full arrhythmia-gene export.

- `is_canonical`
  - True if the source supplement marked the transcript as canonical.

- `in_file_s3_sequencing_center_disagreement`
  - True if the variant appeared in the sequencing-center disagreement supplement.

- `in_file_s4_functional_ep_dataset`
  - True if the protein-level variant matched the in-vitro electrophysiology supplement.

- `primary_validation_eligible`
  - True when the row is canonical and not flagged for sequencing-center disagreement.

Primary-eligible rows:

| Group | Count |
| --- | ---: |
| Primary eligible | 2,650 |
| Not primary eligible | 104 |

Primary-eligible post-in-vitro label counts:

| Label | Count |
| --- | ---: |
| VUS | 2,409 |
| Pathogenic | 122 |
| Benign | 119 |

## Gene Counts

| Gene | Count |
| --- | ---: |
| RYR2 | 685 |
| ANK2 | 661 |
| SCN5A | 436 |
| CACNA1C | 302 |
| KCNH2 | 274 |
| KCNQ1 | 162 |
| LMNA | 130 |
| KCNE1 | 48 |
| KCNJ2 | 39 |
| KCNE2 | 17 |

## Consequence Counts

| Consequence | Count |
| --- | ---: |
| missense | 2,497 |
| Missing | 145 |
| frameshift | 58 |
| nonsense | 27 |
| splice_region | 27 |

## How This Differs From the Earlier 300-Variant Validation Set

The earlier study did not use the full 2,754-row dataset for the main external validation table. It selected a balanced 300-variant subset:

- 100 Pathogenic
- 100 VUS
- 100 Benign

That 300-row set was chosen before model scoring and was designed to give interpretable sensitivity/specificity estimates despite the public eMERGE supplement being overwhelmingly VUS. The full dataset here should be used when the new Cardiogenetics study needs the entire external eMERGE arrhythmia-gene variant pool.

## Important Caveat

The public eMERGE supplement is genotype-focused. It does not contain row-level patient phenotype values comparable to the internal HiRO/CASPER/VERDICT dataset. Therefore, phenotype features in the model-ready files were left missing rather than imputed.

## Integrity Check

SHA-256 for `data/emerge_arrhythmia_variant_level_full.csv` at handoff:

`d8c12c9f8e99b5f66d5b8f5376ffe4b9ce5956c206d61579567ec4773895e2f9`
