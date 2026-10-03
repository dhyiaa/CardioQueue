# Binary Model Results

Date: 2026-07-02

## Summary

The primary model branch is a binary CatBoost classifier:

```text
Pathogenic / Likely Pathogenic vs Benign / Likely Benign
```

This is the model branch most directly comparable to CardioBoost, because
CardioBoost was also a binary pathogenicity model with a continuous probability
score and clinical-style high-confidence thresholds:

```text
Pr >= 0.9  -> disease-causing / P-LP-like
Pr <= 0.1  -> benign / B-LB-like
0.1-0.9    -> indeterminate / VUS-like deferral
```

The binary model should remain the primary paper model. The VUS-facing analysis
should be framed as binary pathogenicity scoring plus a VUS/indeterminate
deferral zone, not as a claim that VUS is a stable biological class.

## Corrected Data Matrix

The results below use the corrected eMERGE-rescued GRCh38 feature matrix:

- `datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv`

Key counts:

| Item | Count |
|---|---:|
| Total corrected variant-level rows | 85,677 |
| Primary binary rows | 42,990 |
| Benign binary rows | 33,859 |
| Pathogenic binary rows | 9,131 |
| VUS scoring rows | 41,731 |
| VEP annotated rows | 85,660 |
| SpliceAI available rows | 77,141 |

## eMERGE Coordinate Rescue

The earlier poor eMERGE result was traced to coordinate/build mismatch. eMERGE
source coordinates were lifted from hg19/GRCh37 to GRCh38 and validated against
the same GRCh38 FASTA used by VEP.

| QC item | Result |
|---|---:|
| eMERGE source rows processed | 2,754 |
| Lifted single-mapping rows | 2,754 |
| GRCh38 REF allele matches | 2,754 |
| Coordinate QC pass | 2,754 |
| VEP returned after rescue | 2,754 |

Source-vs-VEP consequence agreement after rescue:

| Agreement class | Rows |
|---|---:|
| Broad class match | 2,578 |
| Missing source or VEP consequence | 145 |
| Partial LoF-related match | 23 |
| Broad mismatch | 8 |

Conclusion: the original eMERGE external failure should be considered invalid.
After rescue, eMERGE behaves like a strong external validation set.

## Data Splits

Three binary split strategies were used.

| Split strategy | Purpose |
|---|---|
| Internal grouped split | Baseline train/validation/test performance |
| Source-held-out split | Main external validation on HiRO, eMERGE, and CardioBoost |
| Gene-stress split | Sparse-gene vs well-represented gene stress test |

Source-held-out binary split:

| Split | Benign | Pathogenic | Rows |
|---|---:|---:|---:|
| Train | 28,631 | 7,545 | 36,176 |
| Validation | 5,053 | 1,331 | 6,384 |
| External HiRO | 42 | 27 | 69 |
| External eMERGE | 80 | 96 | 176 |
| External CardioBoost | 53 | 132 | 185 |

The HiRO external binary set is small because it is variant-level and excludes
VUS/unresolved rows. The full HiRO source-record dataset remains preserved and
is analyzed separately.

## Binary Model Performance

Strict threshold `0.5`:

| Model | Split | Rows | AUROC | AUPRC | Sensitivity | Specificity | PPV |
|---|---|---:|---:|---:|---:|---:|---:|
| Internal grouped | Test | 6,449 | 0.9990 | 0.9971 | 0.9883 | 0.9968 | 0.9883 |
| Source-held-out | External HiRO | 69 | 0.9877 | 0.9830 | 0.9630 | 0.8810 | 0.8387 |
| Source-held-out | External eMERGE | 176 | 0.9910 | 0.9914 | 0.9792 | 0.9375 | 0.9495 |
| Source-held-out | External CardioBoost | 185 | 0.9520 | 0.9812 | 0.9545 | 0.8113 | 0.9265 |
| Gene-stress | Sparse-gene stress test | 5,890 | 0.9979 | 0.9574 | 0.9848 | 0.9947 | 0.8966 |
| Gene-stress | Well-represented test | 5,565 | 0.9999 | 0.9996 | 0.9872 | 0.9976 | 0.9924 |

## CardioBoost-Comparable Deferral Performance

The closest comparison to CardioBoost is the binary model with `0.1/0.9`
clinical thresholds.

| Set | Rows | Sensitivity | Specificity | PPV | NPV | MCC | Deferral rate |
|---|---:|---:|---:|---:|---:|---:|---:|
| External combined | 430 | 0.937 | 0.943 | 0.960 | 0.977 | 0.876 | 0.119 |
| HiRO | 69 | 0.926 | 0.929 | 0.893 | 1.000 | 0.849 | 0.087 |
| eMERGE | 176 | 0.948 | 0.988 | 0.989 | 1.000 | 0.933 | 0.114 |
| CardioBoost | 185 | 0.932 | 0.887 | 0.953 | 0.903 | 0.806 | 0.135 |

External combined confusion matrix with `0.1/0.9` deferral:

| True \ Pred | B/LB | VUS/Deferred | P/LP |
|---|---:|---:|---:|
| B/LB | 127 | 38 | 10 |
| P/LP | 3 | 13 | 239 |

Interpretation: the binary model with deferral has strong sensitivity and
specificity, high PPV/NPV, and a modest deferral rate. This is the cleanest
CardioBoost-style result.

## HiRO Source-Record Binary + Deferral Analysis

The full HiRO table has 482 source records:

| HiRO source-record status | Count |
|---|---:|
| Total source records | 482 |
| With prediction | 408 |
| No prediction | 74 |

Binary model with `0.1/0.9` deferral, predicted rows only:

| True \ Pred | Benign | VUS/Deferred | Pathogenic |
|---|---:|---:|---:|
| Benign | 149 | 23 | 18 |
| VUS | 60 | 74 | 38 |
| Pathogenic | 0 | 5 | 41 |

Key interpretation:

| Metric | Value |
|---|---:|
| Exact 3-call accuracy on predicted source records | 0.647 |
| P/LP sensitivity | 41 / 46 = 0.891 |
| P/LP specificity | 0.887 |
| P/LP PPV | 0.719 |
| P/LP NPV | 1.000 |
| Pathogenic -> Benign hard errors | 0 |
| Pathogenic -> VUS deferrals | 5 |

This is clinically useful because it detects most pathogenic HiRO source
records while avoiding pathogenic-to-benign hard errors.

## ROC Curves

Generated ROC figures:

- `results/model_performance/source_heldout_hiro_emerge_rescued/roc_curves/roc_internal_grouped.png`
- `results/model_performance/source_heldout_hiro_emerge_rescued/roc_curves/roc_source_heldout.png`
- `results/model_performance/source_heldout_hiro_emerge_rescued/roc_curves/roc_external_validation_sets.png`
- `results/model_performance/source_heldout_hiro_emerge_rescued/roc_curves/roc_gene_stress.png`

External ROC AUC:

| Set | AUROC |
|---|---:|
| External HiRO | 0.9877 |
| External eMERGE | 0.9910 |
| External CardioBoost | 0.9520 |
| External combined | 0.9768 |

## Interpretation

The binary model is currently the strongest and most defensible primary model.
It is also the closest model to CardioBoost. The best paper framing is:

> We trained a cardiogenetics-specific binary pathogenic-vs-benign model and
> applied CardioBoost-style high-confidence thresholds to defer uncertain cases.

This framing avoids treating VUS as biological ground truth while still giving a
clinically useful VUS triage output.

## Key Artifacts

- `results/model_performance/source_heldout_hiro_emerge_rescued/PERFORMANCE_REVIEW.md`
- `results/model_performance/source_heldout_hiro_emerge_rescued/performance_metrics_by_subset_policy.tsv`
- `results/model_performance/source_heldout_hiro_emerge_rescued/model_level_performance_summary.tsv`
- `results/model_performance/source_heldout_hiro_emerge_rescued/roc_curves/roc_auc_summary.tsv`
- `results/model_performance/hiro_source_record_3class/hiro_source_record_predictions.tsv`
- `results/models/primary_binary_catboost_v0_strict_source_heldout_hiro_emerge_rescued/metrics.json`
