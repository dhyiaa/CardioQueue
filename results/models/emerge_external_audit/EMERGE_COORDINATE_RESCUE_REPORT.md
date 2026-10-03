# eMERGE Coordinate Rescue Report

Date: 2026-07-02

## Why This Was Needed

The first source-held-out CatBoost run showed poor eMERGE external performance
(`AUROC 0.5802`, `AUPRC 0.6429`). A follow-up audit showed this should not be
interpreted as model failure: eMERGE rows were being treated as GRCh38-like
coordinates even though the source coordinates were hg19/GRCh37. That caused
dbNSFP misses, gnomAD misses, and VEP consequences that looked intronic or
intergenic for variants that should have been coding.

## Rescue Performed

Input eMERGE source records were lifted from hg19/GRCh37 to GRCh38 using the
UCSC `hg19ToHg38.over.chain.gz` chain and validated against the local GRCh38
primary assembly FASTA used by VEP.

Rescue outputs:

- `datasets/emerge/full_arrhythmia_gene_dataset/interim/emerge_coordinate_rescue_liftover.tsv`
- `datasets/emerge/full_arrhythmia_gene_dataset/interim/emerge_coordinate_rescue_liftover_with_vep_qc.tsv`
- `datasets/variant_registry/interim/clean_combined_variant_registry_with_hiro_emerge_rescues.tsv`
- `datasets/modeling/interim/final_modeling_table_local_features_hiro_emerge_rescued_gnomad_foldx_vep_hgvs_spliceai_hiro_agg.tsv`
- `datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv`
- `results/models/primary_binary_catboost_v0_strict_source_heldout_hiro_emerge_rescued/`

## Coordinate QC

| QC item | Result |
|---|---:|
| eMERGE source rows processed | 2,754 |
| Lifted single-mapping rows | 2,754 |
| GRCh38 REF allele matches | 2,754 |
| Coordinate QC pass | 2,754 |
| VEP returned after rescue | 2,754 |

## Source vs VEP Consequence QC

| Agreement class | Rows |
|---|---:|
| Broad class match | 2,578 |
| Missing source or VEP consequence | 145 |
| Partial LoF-related match | 23 |
| Broad class mismatch | 8 |

Interpretation: the coordinate rescue is strongly supported. Only 8 rows show a
broad source-vs-VEP consequence mismatch after liftOver.

## Corrected Modeling Matrix

| Item | Count |
|---|---:|
| Corrected registry/modeling rows | 85,677 |
| VEP ok | 85,660 |
| VEP missing | 17 |
| SpliceAI ok | 77,141 |
| SpliceAI missing | 8,536 |
| Primary binary training/evaluation rows | 42,990 |
| VUS scoring rows | 41,731 |

## Corrected Source-Held-Out Split

| Split | Benign | Pathogenic | Rows |
|---|---:|---:|---:|
| Train | 28,631 | 7,545 | 36,176 |
| Validation | 5,053 | 1,331 | 6,384 |
| External HiRO | 42 | 27 | 69 |
| External eMERGE | 80 | 96 | 176 |
| External CardioBoost | 53 | 132 | 185 |

Source-held-out rows are excluded from training. HiRO has priority over eMERGE,
and eMERGE has priority over CardioBoost for rows with overlapping source flags.

## Corrected Source-Held-Out CatBoost Result

| Held-out set | Rows | AUROC | AUPRC | Sensitivity at 0.5 | Specificity at 0.5 |
|---|---:|---:|---:|---:|---:|
| eMERGE | 176 | 0.9910 | 0.9914 | 0.9792 | 0.9375 |
| HiRO | 69 | 0.9877 | 0.9830 | 0.9630 | 0.8810 |
| CardioBoost | 185 | 0.9520 | 0.9812 | 0.9545 | 0.8113 |
| Validation | 6,384 | 0.9999 | 0.9996 | 0.9947 | 0.9970 |

## Interpretation

The prior eMERGE external result should be considered invalid because it was
driven by coordinate/build mismatch. After hg19-to-GRCh38 liftOver, REF allele
validation, dbNSFP/VEP/SpliceAI feature rebuilding, and source-held-out
retraining, eMERGE external performance is high and consistent with the other
trusted held-out sources.

This does not prove the final model is ready for publication. The corrected
result still needs the normal checks: leakage review, calibration, gene-level
stress testing, feature ablations, and careful reporting of source overlap.
