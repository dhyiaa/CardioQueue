# Three-Class Model Tuning Review

This review addresses the important methodological concern: model tuning must not be chosen because it improves HiRO, eMERGE, or CardioBoost after looking at their results. External datasets must stay external.

## Scientific Rule

Tuning decisions should be made using only the internal validation split. HiRO, eMERGE, CardioBoost, and ClinVar source-specific summaries are then used only as held-out evaluation strata.

Therefore:

- Do not choose thresholds or weights because they make HiRO look better.
- Use validation metrics to choose candidate operating points.
- Report external performance after the rule is frozen.
- Keep the original true 3-class argmax model as the cleanest baseline.
- Treat any P/LP-sensitivity-tuned model as a secondary operating point unless it improves validation and external performance consistently.

## Baseline True 3-Class Model

The original true 3-class model directly predicts `Benign`, `VUS`, or `Pathogenic` by argmax probability.

| Dataset | Rows | Accuracy | Macro-F1 | P/LP Sensitivity | P/LP Specificity | P/LP PPV |
|---|---:|---:|---:|---:|---:|---:|
| Validation | 12,361 | 0.9338 | 0.9101 | 0.8843 | 0.9713 | 0.7905 |
| HiRO source records | 325 | 0.7662 | 0.7770 | 0.7143 | 0.9931 | 0.9259 |
| eMERGE source records | 1,974 | 0.8896 | 0.6565 | 0.5816 | 0.9760 | 0.5588 |
| CardioBoost source records | 195 | 0.5949 | 0.4735 | 0.6232 | 0.9649 | 0.9773 |

Interpretation: the baseline is specific and clean, but P/LP sensitivity is modest on external source-record sets.

## Validation-Only Threshold Selection

I tested threshold rules on the existing 3-class probabilities:

```text
Call Pathogenic if prob_Pathogenic >= threshold.
Otherwise choose Benign vs VUS by whichever probability is larger.
```

Thresholds were selected using validation only. The selected thresholds did not use HiRO/eMERGE/CardioBoost.

| Validation Selection Rule | Selected Threshold | Validation Accuracy | Validation Macro-F1 | Validation P/LP Sensitivity | Validation P/LP Specificity | Validation PPV |
|---|---:|---:|---:|---:|---:|---:|
| Max validation macro-F1 | 0.700 | 0.9401 | 0.9181 | 0.8325 | 0.9850 | 0.8704 |
| Max validation weighted kappa | 0.675 | 0.9398 | 0.9181 | 0.8415 | 0.9837 | 0.8615 |
| Max validation balanced accuracy | 0.400 | 0.9278 | 0.9019 | 0.9083 | 0.9620 | 0.7426 |
| Constrained sensitivity/specificity/PPV rule | 0.425 | 0.9291 | 0.9035 | 0.9016 | 0.9643 | 0.7528 |

External source-record performance after freezing those validation-selected thresholds:

| Rule | External Set | Rows | Accuracy | Macro-F1 | P/LP Sensitivity | P/LP Specificity | P/LP PPV |
|---|---|---:|---:|---:|---:|---:|---:|
| Max validation macro-F1 | HiRO | 325 | 0.7538 | 0.7463 | 0.6000 | 0.9931 | 0.9130 |
| Max validation macro-F1 | eMERGE | 1,974 | 0.8946 | 0.6520 | 0.4796 | 0.9867 | 0.6528 |
| Max validation macro-F1 | CardioBoost | 195 | 0.4821 | 0.4300 | 0.4638 | 0.9825 | 0.9846 |
| Max validation balanced accuracy | HiRO | 325 | 0.7692 | 0.7842 | 0.7429 | 0.9931 | 0.9286 |
| Max validation balanced accuracy | eMERGE | 1,974 | 0.8820 | 0.6473 | 0.6224 | 0.9659 | 0.4880 |
| Max validation balanced accuracy | CardioBoost | 195 | 0.6308 | 0.4847 | 0.6739 | 0.9474 | 0.9688 |
| Constrained rule | HiRO | 325 | 0.7692 | 0.7842 | 0.7429 | 0.9931 | 0.9286 |
| Constrained rule | eMERGE | 1,974 | 0.8870 | 0.6556 | 0.6122 | 0.9717 | 0.5310 |
| Constrained rule | CardioBoost | 195 | 0.6154 | 0.4795 | 0.6522 | 0.9474 | 0.9677 |

Conclusion: validation-selected thresholds are scientifically valid, but they do not clearly beat the baseline in a general way. The sensitivity-oriented validation thresholds recover some P/LP sensitivity, but they reduce validation macro-F1 and do not consistently improve all external sources.

## Sequential Class-Weight Retraining

I retried weight tuning one candidate at a time, not in parallel. Each model used the same source-held-out split and was evaluated on validation first, then external sets.

Completed candidates:

| Candidate | Class Weights |
|---|---|
| `path_3_5` | Benign 1.0, VUS 0.9, Pathogenic 3.5 |
| `path_4_vus_0_8` | Benign 1.0, VUS 0.8, Pathogenic 4.0 |
| `path_5_vus_0_7` | Benign 1.0, VUS 0.7, Pathogenic 5.0 |

Argmax results:

| Candidate | Dataset | Rows | Accuracy | Macro-F1 | P/LP Sensitivity | P/LP Specificity | P/LP PPV |
|---|---|---:|---:|---:|---:|---:|---:|
| `path_3_5` | Validation | 12,361 | 0.9282 | 0.9024 | 0.9083 | 0.9625 | 0.7449 |
| `path_3_5` | HiRO | 325 | 0.7754 | 0.7887 | 0.7429 | 0.9931 | 0.9286 |
| `path_3_5` | eMERGE | 1,974 | 0.8789 | 0.6459 | 0.6224 | 0.9659 | 0.4880 |
| `path_3_5` | CardioBoost | 195 | 0.6103 | 0.4667 | 0.6667 | 0.9298 | 0.9583 |
| `path_4_vus_0_8` | Validation | 12,361 | 0.9213 | 0.8933 | 0.9226 | 0.9538 | 0.7066 |
| `path_4_vus_0_8` | HiRO | 325 | 0.7723 | 0.7864 | 0.7429 | 0.9931 | 0.9286 |
| `path_4_vus_0_8` | eMERGE | 1,974 | 0.8652 | 0.6228 | 0.6224 | 0.9611 | 0.4552 |
| `path_4_vus_0_8` | CardioBoost | 195 | 0.6205 | 0.4642 | 0.6884 | 0.8947 | 0.9406 |
| `path_5_vus_0_7` | Validation | 12,361 | 0.9072 | 0.8749 | 0.9406 | 0.9368 | 0.6424 |
| `path_5_vus_0_7` | HiRO | 325 | 0.7754 | 0.7887 | 0.7429 | 0.9931 | 0.9286 |
| `path_5_vus_0_7` | eMERGE | 1,974 | 0.8582 | 0.6238 | 0.6837 | 0.9435 | 0.3873 |
| `path_5_vus_0_7` | CardioBoost | 195 | 0.6769 | 0.4893 | 0.7536 | 0.8947 | 0.9455 |

## Interpretation

Weight tuning can improve P/LP sensitivity, but it is not a free improvement. As pathogenic weight increases:

- validation P/LP sensitivity improves;
- validation macro-F1 and accuracy decline;
- external P/LP sensitivity often improves;
- external PPV, especially in eMERGE, can decline;
- stronger weights create more pathogenic calls from VUS rows.

The most defensible sensitivity-tuned candidate so far is `path_3_5` because it improves validation P/LP sensitivity from 0.8843 to 0.9083 while keeping validation macro-F1 reasonably high at 0.9024. Stronger candidates push sensitivity higher but degrade validation macro-F1 and precision more.

## Recommendation

For the paper:

| Role | Model |
|---|---|
| Primary 3-class model | Original true 3-class argmax |
| Sensitivity-tuned 3-class model | `path_3_5` argmax |
| Exploratory high-sensitivity model | `path_5_vus_0_7` argmax, only if explicitly framed as a triage-sensitive operating point |
| Do not use as primary | Threshold chosen because it helps HiRO specifically |

The safest wording is:

> We evaluated prespecified sensitivity-oriented 3-class weighting schemes selected on internal validation. A modest pathogenic upweighting improved P/LP sensitivity across external source-record sets but traded off some validation macro-F1 and precision. Therefore, the original 3-class model remains the balanced multiclass classifier, while the modestly upweighted model is reported as a sensitivity-oriented secondary operating point.

## Output Files

- Validation-only threshold grid: `results/model_performance/three_class_scientific_tuning/validation_threshold_grid.tsv`
- Validation-selected threshold report: `results/model_performance/three_class_scientific_tuning/VALIDATION_SELECTED_THRESHOLD_REPORT.md`
- Sequential weight-tuning results: `results/model_performance/three_class_weight_tuning/three_class_weight_tuning_results.completed_candidates.tsv`
- Candidate models:
  - `results/model_performance/three_class_weight_tuning/path_3_5.cbm`
  - `results/model_performance/three_class_weight_tuning/path_4_vus_0_8.cbm`
  - `results/model_performance/three_class_weight_tuning/path_5_vus_0_7.cbm`
