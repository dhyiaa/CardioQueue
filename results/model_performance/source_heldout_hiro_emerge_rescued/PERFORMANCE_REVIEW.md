# Corrected Model Performance Review

Date: 2026-07-02

This report uses the eMERGE-rescued GRCh38 feature matrix and corrected source-held-out CatBoost models. The primary model is binary P/LP vs B/LB. The 3-call view is post-hoc: probability `<0.10` = B/LB, `>=0.90` = P/LP, otherwise VUS/defer.

## Model-Level Binary Performance at Threshold 0.5

| model | split | rows | pathogenic | benign | AUROC | AUPRC | Sensitivity@0.5 | Specificity@0.5 | PPV@0.5 | TP@0.5 | FP@0.5 | TN@0.5 | FN@0.5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Internal grouped corrected | train | 30092 | 6391 | 23701 | 0.9998 | 0.9993 | 0.9928 | 0.9973 | 0.9902 | 6345 | 63 | 23638 | 46 |
| Internal grouped corrected | validation | 6449 | 1370 | 5079 | 0.9993 | 0.9977 | 0.9854 | 0.9955 | 0.9832 | 1350 | 23 | 5056 | 20 |
| Internal grouped corrected | test | 6449 | 1370 | 5079 | 0.9990 | 0.9971 | 0.9883 | 0.9968 | 0.9883 | 1354 | 16 | 5063 | 16 |
| Source-held-out corrected | train | 36176 | 7545 | 28631 | 0.9998 | 0.9992 | 0.9932 | 0.9978 | 0.9918 | 7494 | 62 | 28569 | 51 |
| Source-held-out corrected | validation | 6384 | 1331 | 5053 | 0.9999 | 0.9996 | 0.9947 | 0.9970 | 0.9888 | 1324 | 15 | 5038 | 7 |
| Source-held-out corrected | external_cardioboost | 185 | 132 | 53 | 0.9520 | 0.9812 | 0.9545 | 0.8113 | 0.9265 | 126 | 10 | 43 | 6 |
| Source-held-out corrected | external_emerge | 176 | 96 | 80 | 0.9910 | 0.9914 | 0.9792 | 0.9375 | 0.9495 | 94 | 5 | 75 | 2 |
| Source-held-out corrected | external_hiro | 69 | 27 | 42 | 0.9877 | 0.9830 | 0.9630 | 0.8810 | 0.8387 | 26 | 5 | 37 | 1 |
| Gene-stress corrected | train | 25970 | 6208 | 19762 | 0.9997 | 0.9990 | 0.9929 | 0.9972 | 0.9910 | 6164 | 56 | 19706 | 44 |
| Gene-stress corrected | validation | 5565 | 1328 | 4237 | 0.9988 | 0.9975 | 0.9895 | 0.9950 | 0.9843 | 1314 | 21 | 4216 | 14 |
| Gene-stress corrected | sparse_gene_stress_test | 5890 | 264 | 5626 | 0.9979 | 0.9574 | 0.9848 | 0.9947 | 0.8966 | 260 | 30 | 5596 | 4 |
| Gene-stress corrected | well_represented_test | 5565 | 1331 | 4234 | 0.9999 | 0.9996 | 0.9872 | 0.9976 | 0.9924 | 1314 | 10 | 4224 | 17 |

## Requested Metrics: Source-Held-Out Model

| subset | policy | rows | exact_3class_accuracy | weighted_kappa | macro_f1 | p_lp_sensitivity | specificity_no_p_lp_escalation | ppv_p_lp_calls | npv_b_lb_calls | binary_mcc_p_call_vs_not_p | vus_deferral_rate | p_lp_positive_calls | true_positives | false_positives | p_lp_to_vus_deferrals | p_lp_to_b_lb_errors | b_lb_to_p_lp_escalations |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| validation | binary_threshold_0_5 | 6384 | 0.9966 | 0.9896 | 0.6632 | 0.9947 | 0.9970 | 0.9888 | 0.9986 | 0.9896 | 0.0000 | 1339 | 1324 | 15 | 0 | 7 | 15 |
| validation | defer_0.1_0.9 | 6384 | 0.9900 | 0.9895 | 0.6622 | 0.9872 | 0.9988 | 0.9955 | 0.9996 | 0.9891 | 0.0088 | 1320 | 1314 | 6 | 15 | 2 | 6 |
| external_all | binary_threshold_0_5 | 430 | 0.9326 | 0.8589 | 0.6196 | 0.9647 | 0.8857 | 0.9248 | 0.9451 | 0.8601 | 0.0000 | 266 | 246 | 20 | 0 | 9 | 20 |
| external_all | defer_0.1_0.9 | 430 | 0.8512 | 0.8653 | 0.5937 | 0.9373 | 0.9429 | 0.9598 | 0.9769 | 0.8758 | 0.1186 | 249 | 239 | 10 | 13 | 3 | 10 |
| external_hiro | binary_threshold_0_5 | 69 | 0.9130 | 0.8222 | 0.6072 | 0.9630 | 0.8810 | 0.8387 | 0.9737 | 0.8280 | 0.0000 | 31 | 26 | 5 | 0 | 1 | 5 |
| external_hiro | defer_0.1_0.9 | 69 | 0.8696 | 0.8604 | 0.6061 | 0.9259 | 0.9286 | 0.8929 | 1.0000 | 0.8493 | 0.0870 | 28 | 25 | 3 | 2 | 0 | 3 |
| external_emerge | binary_threshold_0_5 | 176 | 0.9602 | 0.9195 | 0.6398 | 0.9792 | 0.9375 | 0.9495 | 0.9740 | 0.9201 | 0.0000 | 99 | 94 | 5 | 0 | 2 | 5 |
| external_emerge | defer_0.1_0.9 | 176 | 0.8807 | 0.9266 | 0.6190 | 0.9479 | 0.9875 | 0.9891 | 1.0000 | 0.9325 | 0.1136 | 92 | 91 | 1 | 5 | 0 | 1 |
| external_cardioboost | binary_threshold_0_5 | 185 | 0.9135 | 0.7836 | 0.5945 | 0.9545 | 0.8113 | 0.9265 | 0.8776 | 0.7847 | 0.0000 | 136 | 126 | 10 | 0 | 6 | 10 |
| external_cardioboost | defer_0.1_0.9 | 185 | 0.8162 | 0.7666 | 0.5364 | 0.9318 | 0.8868 | 0.9535 | 0.9032 | 0.8056 | 0.1351 | 129 | 123 | 6 | 6 | 3 | 6 |

## Deferral-Policy Confusion Matrices


### external_all

| true_label | pred_Benign | pred_VUS | pred_Pathogenic |
| --- | --- | --- | --- |
| true_Benign | 127 | 38 | 10 |
| true_VUS | 0 | 0 | 0 |
| true_Pathogenic | 3 | 13 | 239 |

### external_hiro

| true_label | pred_Benign | pred_VUS | pred_Pathogenic |
| --- | --- | --- | --- |
| true_Benign | 35 | 4 | 3 |
| true_VUS | 0 | 0 | 0 |
| true_Pathogenic | 0 | 2 | 25 |

### external_emerge

| true_label | pred_Benign | pred_VUS | pred_Pathogenic |
| --- | --- | --- | --- |
| true_Benign | 64 | 15 | 1 |
| true_VUS | 0 | 0 | 0 |
| true_Pathogenic | 0 | 5 | 91 |

### external_cardioboost

| true_label | pred_Benign | pred_VUS | pred_Pathogenic |
| --- | --- | --- | --- |
| true_Benign | 28 | 19 | 6 |
| true_VUS | 0 | 0 | 0 |
| true_Pathogenic | 3 | 6 | 123 |

### validation

| true_label | pred_Benign | pred_VUS | pred_Pathogenic |
| --- | --- | --- | --- |
| true_Benign | 5006 | 41 | 6 |
| true_VUS | 0 | 0 | 0 |
| true_Pathogenic | 2 | 15 | 1314 |

## Files

- `results/model_performance/source_heldout_hiro_emerge_rescued/roc_curves/roc_internal_grouped.png`
- `results/model_performance/source_heldout_hiro_emerge_rescued/roc_curves/roc_source_heldout.png`
- `results/model_performance/source_heldout_hiro_emerge_rescued/roc_curves/roc_gene_stress.png`
- `results/model_performance/source_heldout_hiro_emerge_rescued/roc_curves/roc_external_validation_sets.png`
- `results/model_performance/source_heldout_hiro_emerge_rescued/roc_curves/roc_auc_summary.tsv`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_external_all_binary_threshold_0_5.png`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_external_all_binary_threshold_0_5.tsv`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_external_all_defer_0.1_0.9.png`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_external_all_defer_0.1_0.9.tsv`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_external_cardioboost_binary_threshold_0_5.png`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_external_cardioboost_binary_threshold_0_5.tsv`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_external_cardioboost_defer_0.1_0.9.png`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_external_cardioboost_defer_0.1_0.9.tsv`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_external_emerge_binary_threshold_0_5.png`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_external_emerge_binary_threshold_0_5.tsv`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_external_emerge_defer_0.1_0.9.png`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_external_emerge_defer_0.1_0.9.tsv`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_external_hiro_binary_threshold_0_5.png`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_external_hiro_binary_threshold_0_5.tsv`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_external_hiro_defer_0.1_0.9.png`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_external_hiro_defer_0.1_0.9.tsv`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_validation_binary_threshold_0_5.png`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_validation_binary_threshold_0_5.tsv`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_validation_defer_0.1_0.9.png`
- `results/model_performance/source_heldout_hiro_emerge_rescued/confusion_matrix_validation_defer_0.1_0.9.tsv`
- `results/model_performance/source_heldout_hiro_emerge_rescued/evaluation_manifest.json`
- `results/model_performance/source_heldout_hiro_emerge_rescued/model_level_performance_summary.tsv`
- `results/model_performance/source_heldout_hiro_emerge_rescued/performance_metrics_by_subset_policy.tsv`

## ROC AUC Curves

| Model | Split | Rows | Pathogenic | Benign | AUROC |
| --- | --- | ---: | ---: | ---: | ---: |
| Internal grouped | Validation | 6,449 | 1,370 | 5,079 | 0.9993 |
| Internal grouped | Test | 6,449 | 1,370 | 5,079 | 0.9990 |
| Source-held-out | Validation | 6,384 | 1,331 | 5,053 | 0.9999 |
| Source-held-out | External HiRO | 69 | 27 | 42 | 0.9877 |
| Source-held-out | External eMERGE | 176 | 96 | 80 | 0.9910 |
| Source-held-out | External CardioBoost | 185 | 132 | 53 | 0.9520 |
| Source-held-out | External combined | 430 | 255 | 175 | 0.9768 |
| Gene-stress | Validation | 5,565 | 1,328 | 4,237 | 0.9988 |
| Gene-stress | Well-represented test | 5,565 | 1,331 | 4,234 | 0.9999 |
| Gene-stress | Sparse-gene stress test | 5,890 | 264 | 5,626 | 0.9979 |

## CardioBoost Design Lessons

CardioBoost was designed as a disease-specific classifier for rare missense
variants in inherited cardiomyopathies and inherited arrhythmia syndromes. It
used curated pathogenic variants, presumed benign variants from 2,090 healthy
volunteers, 76 functional annotations, nested cross-validation over multiple
algorithm families, and selected AdaBoost as the final learner. It reported both
continuous ROC/PR performance and clinically interpretable high-confidence
thresholds: `Pr >= 0.9` for disease-causing, `Pr <= 0.1` for B/LB, and
intermediate scores as indeterminate/VUS-like.

Direct methodological lessons for this project:

| Lesson | Action for this project |
| --- | --- |
| Disease specificity matters | Keep cardiac-gene and disease-group framing; do not present this as a genome-wide pathogenicity model. |
| Gene/disease grouping beats tiny per-gene models | Continue reporting disease-group and sparse-gene stress performance rather than overfitting per-gene classifiers. |
| High-confidence thresholds are clinically important | Continue reporting `>=0.9`, `<0.1`, and deferral rates, not only AUROC. |
| External validation matters | Emphasize HiRO, eMERGE, and CardioBoost source-held-out testing; internal split alone is not enough. |
| Outcome/phenotype validation is powerful | Use HiRO phenotype/source-record data for a patient-event-level sensitivity analysis. |
| Missense-only scope is cleaner | Consider a missense-only primary model and a broader all-variant secondary model. |
| Feature leakage must be audited | Benchmark with and without VEP consequence, gene identity, ClinVar-derived fields, and source-derived aggregate fields. |
| PR-AUC is essential with imbalance | Continue reporting AUPRC alongside AUROC. |
