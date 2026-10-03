# Primary Binary CatBoost v0 Report

Generated: 2026-07-02

## Input Data

- Source matrix: `datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx_vep_hgvs_spliceai_hiro_agg.tsv`
- Split/weight table: `datasets/modeling/ready/modeling_table_with_splits_weights.tsv`
- Primary binary table: `datasets/modeling/ready/primary_binary_catboost_table.tsv`
- Total matrix rows: 86,889
- Primary binary rows: 43,201
- Binary labels: 33,978 benign, 9,223 pathogenic
- VUS rows reserved for later scoring/analysis: 42,827

## Split Strategies Built

| Split strategy | Purpose | Train | Validation | Test / held-out |
|---|---:|---:|---:|---:|
| Internal grouped split | Development, debugging, hyperparameter search | 30,240 | 6,480 | 6,481 internal test |
| Source-held-out split | External/source generalization | 36,297 | 6,406 | HiRO 70, eMERGE 203, CardioBoost 225 |
| Gene-stress split | Sparse-gene and gene-transfer stress test | 26,114 | 5,596 | sparse genes 5,894; well-represented genes 5,597 |

## Training Setup

- Model: CatBoost binary classifier
- Target: Pathogenic/likely pathogenic vs benign/likely benign
- VUS excluded from supervised binary training
- Label conflicts and invalid rows excluded from supervised binary training
- Sample weighting: `source_confidence_weight * class_weight`
- Pathogenic rows are upweighted to address class imbalance and clinical cost of false negatives.
- Leakage controls: feature selection excludes labels, review stars, source flags, source-derived label counts, patient/source IDs, variant IDs, raw source labels, and source-record fields.

## Model Artifacts

| Run | Output directory | Features | Best iteration |
|---|---|---:|---:|
| Strict internal grouped | `results/models/primary_binary_catboost_v0_strict_internal/` | 387 | 1,072 |
| Strict source-held-out | `results/models/primary_binary_catboost_v0_strict_source_heldout/` | 387 | 115 |
| Strict gene-stress | `results/models/primary_binary_catboost_v0_strict_gene_stress/` | 387 | 420 |

## Internal Grouped Results

| Set | Rows | AUROC | AUPRC | Sensitivity at 0.5 | Specificity at 0.5 | PPV at 0.5 |
|---|---:|---:|---:|---:|---:|---:|
| Train | 30,240 | 0.9998 | 0.9994 | 0.9938 | 0.9977 | 0.9917 |
| Validation | 6,480 | 0.9995 | 0.9984 | 0.9870 | 0.9965 | 0.9870 |
| Test | 6,481 | 0.9990 | 0.9972 | 0.9877 | 0.9959 | 0.9849 |

Interpretation: useful as a development benchmark, but likely optimistic because ClinVar-like variants and genes are present across train/validation/test.

## Source-Held-Out Results

| Set | Rows | AUROC | AUPRC | Sensitivity at 0.5 | Specificity at 0.5 | PPV at 0.5 |
|---|---:|---:|---:|---:|---:|---:|
| Train | 36,297 | 0.9997 | 0.9990 | 0.9918 | 0.9973 | 0.9899 |
| Validation | 6,406 | 0.9993 | 0.9972 | 0.9903 | 0.9959 | 0.9844 |
| External CardioBoost | 225 | 0.9482 | 0.9747 | 0.9441 | 0.7195 | 0.8544 |
| External eMERGE | 203 | 0.5802 | 0.6429 | 0.0088 | 0.9667 | 0.2500 |
| External HiRO | 70 | 0.9931 | 0.9904 | 0.9630 | 0.8837 | 0.8387 |

Interpretation: CardioBoost and HiRO transfer look promising. eMERGE is a major warning signal and needs source-specific investigation before paper-level claims.

## Gene-Stress Results

| Set | Rows | AUROC | AUPRC | Sensitivity at 0.5 | Specificity at 0.5 | PPV at 0.5 |
|---|---:|---:|---:|---:|---:|---:|
| Train | 26,114 | 0.9997 | 0.9990 | 0.9923 | 0.9970 | 0.9905 |
| Validation | 5,596 | 0.9996 | 0.9988 | 0.9926 | 0.9960 | 0.9874 |
| Sparse-gene stress test | 5,894 | 0.9981 | 0.9545 | 0.9850 | 0.9938 | 0.8822 |
| Well-represented test | 5,597 | 0.9989 | 0.9974 | 0.9873 | 0.9944 | 0.9822 |

Interpretation: the sparse-gene stress result is strong, but it should be validated after auditing gene membership, source composition, and feature leakage risk.

## Top Source-Held-Out Feature Families

Top features include VEP consequence/impact, gnomAD AF/popmax fields, SpliceAI scores, dbNSFP-derived SIFT/CADD/REVEL/FATHMM features, AlphaMissense/ESM scores, and protein/gene context.

## Immediate Follow-Up

1. Audit eMERGE external failures: label definitions, source IDs, coordinate mapping, gene distribution, feature missingness, and whether eMERGE rows were enriched/injected with ClinVar-derived annotations.
2. Generate prediction review tables for false negatives and false positives, especially pathogenic variants missed in eMERGE/CardioBoost/HiRO.
3. Calibrate thresholds on validation data and report clinically meaningful operating points, especially high-specificity and high-sensitivity cutoffs.
4. Train the high-confidence sensitivity model using ClinVar stars >=2 plus trusted non-ClinVar rows.
5. Score VUS rows with the primary binary model for prioritization analysis.
