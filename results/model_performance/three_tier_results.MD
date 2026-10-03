# Three-Tier / Three-Class Model Results

Date: 2026-07-02

## Summary

This branch trains a true 3-class CatBoost model:

```text
Benign / VUS / Pathogenic
```

This is **not** the CardioBoost-equivalent model. CardioBoost was a binary
pathogenicity model with an indeterminate probability zone. This 3-class branch
is a secondary/exploratory VUS modeling arm.

The model is useful because it directly learns the VUS label and performs better
than the binary-deferral approximation for exact 3-class matching on HiRO.
However, VUS is an evidence-status label, not a stable biological class, so the
3-class branch should be interpreted cautiously.

## Training Data

Input:

- `datasets/modeling/ready/modeling_table_three_class_splits_weights_hiro_emerge_rescued.tsv`

Trainable 3-class rows:

| Class | Rows |
|---|---:|
| Benign | 33,859 |
| VUS | 41,731 |
| Pathogenic | 9,131 |
| Total | 84,721 |

Weights:

| Class | Weight |
|---|---:|
| Benign | 1.0 |
| VUS | 0.9 |
| Pathogenic | 2.5 |

Pathogenic was upweighted because it is rarer and more clinically costly to
miss. VUS was mildly downweighted because it is label-noisy and represents
classification uncertainty.

## Three-Class Splits

Source-held-out split:

| Split | Benign | VUS | Pathogenic | Rows |
|---|---:|---:|---:|---:|
| Train | 28,631 | 33,869 | 7,545 | 70,045 |
| Validation | 5,053 | 5,977 | 1,331 | 12,361 |
| External HiRO | 42 | 110 | 27 | 179 |
| External eMERGE | 80 | 1,775 | 96 | 1,951 |
| External CardioBoost | 53 | 0 | 132 | 185 |

Important: CardioBoost has no VUS rows, so it is not a valid true 3-class
external validation set. It can only be interpreted as binary-with-deferral
behavior when tested with a 3-class model.

## Internal 3-Class Performance

| Split | Rows | Accuracy | Macro-F1 | Weighted kappa | MCC | OVR macro AUROC |
|---|---:|---:|---:|---:|---:|---:|
| Validation | 12,708 | 0.9280 | 0.9015 | 0.9122 | 0.8784 | 0.9827 |
| Test | 12,709 | 0.9294 | 0.9044 | 0.9152 | 0.8808 | 0.9840 |

Internal test per-class F1:

| Class | F1 | Recall | Precision |
|---|---:|---:|---:|
| Benign | 0.9603 | 0.9683 | 0.9524 |
| VUS | 0.9280 | 0.9116 | 0.9451 |
| Pathogenic | 0.8249 | 0.8664 | 0.7872 |

## Source-Held-Out 3-Class Performance

| External set | Rows | Accuracy | Macro-F1 | Weighted kappa | MCC | OVR macro AUROC |
|---|---:|---:|---:|---:|---:|---:|
| HiRO | 179 | 0.8492 | 0.8188 | 0.7879 | 0.7160 | 0.9426 |
| eMERGE | 1,951 | 0.8893 | 0.6507 | 0.4977 | 0.4713 | 0.9242 |
| CardioBoost | 185 | 0.5946 | 0.4693 | 0.6351 | 0.4339 | NA |

HiRO is the cleanest true 3-class external validation set because it contains
Benign, VUS, and Pathogenic labels.

## HiRO Variant-Level 3-Class Result

External HiRO confusion matrix:

| True \ Pred | Benign | VUS | Pathogenic |
|---|---:|---:|---:|
| Benign | 34 | 8 | 0 |
| VUS | 8 | 100 | 2 |
| Pathogenic | 0 | 9 | 18 |

HiRO per-class performance:

| Class | F1 | Recall | Precision |
|---|---:|---:|---:|
| Benign | 0.8095 | 0.8095 | 0.8095 |
| VUS | 0.8811 | 0.9091 | 0.8547 |
| Pathogenic | 0.7660 | 0.6667 | 0.9000 |

Interpretation: the true 3-class model recognizes VUS well and makes highly
specific P/LP calls, but it sends some true pathogenic rows to VUS.

## Full HiRO Source-Record Three-Class Result

Mapped back to all 482 HiRO source records:

| HiRO source-record status | Count |
|---|---:|
| Total HiRO source records | 482 |
| With prediction | 408 |
| No prediction | 74 |

All rows:

| True \ Pred | Benign | VUS | Pathogenic | No prediction |
|---|---:|---:|---:|---:|
| Benign | 113 | 77 | 0 | 30 |
| VUS | 16 | 153 | 3 | 35 |
| Pathogenic | 0 | 16 | 30 | 9 |

Predicted rows only:

| True \ Pred | Benign | VUS | Pathogenic |
|---|---:|---:|---:|
| Benign | 113 | 77 | 0 |
| VUS | 16 | 153 | 3 |
| Pathogenic | 0 | 16 | 30 |

Key metrics on predicted HiRO source records:

| Metric | Value |
|---|---:|
| Exact 3-class accuracy | 0.7255 |
| Benign recall | 0.5947 |
| VUS recall | 0.8895 |
| Pathogenic recall | 0.6522 |
| Pathogenic -> Benign hard errors | 0 |
| Benign -> Pathogenic escalations | 0 |
| Pathogenic -> VUS errors | 16 |
| VUS -> Pathogenic escalations | 3 |

Compared with the binary-deferral model, the true 3-class model has higher exact
3-class accuracy and better VUS recognition, but lower P/LP sensitivity.

## eMERGE Three-Class Result

External eMERGE is heavily VUS-skewed:

| Class | Rows |
|---|---:|
| Benign | 80 |
| VUS | 1,775 |
| Pathogenic | 96 |

Confusion matrix:

| True \ Pred | Benign | VUS | Pathogenic |
|---|---:|---:|---:|
| Benign | 54 | 26 | 0 |
| VUS | 105 | 1,626 | 44 |
| Pathogenic | 1 | 40 | 55 |

eMERGE interpretation:

- VUS F1 is strong at 0.938.
- Benign and Pathogenic F1 are modest because the set is dominated by VUS.
- This is a valid 3-class external test, but it is highly imbalanced.

## CardioBoost Three-Class Caveat

CardioBoost has no VUS labels in the local evaluation set:

| Class | Rows |
|---|---:|
| Benign | 53 |
| VUS | 0 |
| Pathogenic | 132 |

Therefore CardioBoost is not a true 3-class external validation set. Its
confusion matrix is still useful as binary-with-deferral behavior:

| True \ Pred | Benign | VUS | Pathogenic |
|---|---:|---:|---:|
| Benign | 27 | 24 | 2 |
| Pathogenic | 4 | 45 | 83 |

The `VUS` predictions here should be interpreted as model deferrals, not as
correct/incorrect VUS classification.

## Gene-Stress 3-Class Performance

| Split | Rows | Accuracy | Macro-F1 | Weighted kappa | MCC | OVR macro AUROC |
|---|---:|---:|---:|---:|---:|---:|
| Well-represented test | 10,864 | 0.9277 | 0.9059 | 0.9163 | 0.8803 | 0.9834 |
| Sparse-gene stress test | 12,295 | 0.8737 | 0.6902 | 0.8199 | 0.7827 | 0.9574 |

Sparse-gene per-class F1:

| Class | F1 |
|---|---:|
| Benign | 0.9583 |
| VUS | 0.8696 |
| Pathogenic | 0.2427 |

Sparse pathogenic prediction remains a weak spot in the 3-class branch.

## Binary + Deferral vs True 3-Class On HiRO

Same 408 predicted HiRO source-record rows:

| Model | Exact 3-class accuracy | P/LP sensitivity | P/LP specificity | P/LP PPV | P/LP NPV |
|---|---:|---:|---:|---:|---:|
| Binary + 0.1/0.9 deferral | 0.647 | 0.891 | 0.887 | 0.719 | 1.000 |
| True 3-class model | 0.725 | 0.652 | 0.992 | 0.909 | 0.957 |

Interpretation:

- Binary + deferral is better for catching P/LP variants.
- True 3-class is better for exact Benign/VUS/Pathogenic matching and specificity.
- The binary + deferral model is more CardioBoost-comparable and should remain primary.
- The true 3-class model should be presented as secondary/exploratory VUS triage.

## Key Artifacts

- `results/model_performance/three_class_catboost_v0_review/THREE_CLASS_MODEL_REVIEW.md`
- `results/model_performance/three_class_catboost_v0_review/three_class_model_level_performance.tsv`
- `results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued/metrics.json`
- `results/models/three_class_catboost_v0_internal_hiro_emerge_rescued/metrics.json`
- `results/models/three_class_catboost_v0_gene_stress_hiro_emerge_rescued/metrics.json`
- `results/model_performance/hiro_source_record_true_three_class/hiro_source_record_true_three_class_predictions.tsv`
- `results/model_performance/hiro_source_record_true_three_class/hiro_source_record_true_three_class_confusion_all.tsv`
- `results/model_performance/hiro_source_record_true_three_class/hiro_source_record_true_three_class_confusion_predicted_only.tsv`
