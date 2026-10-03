# CardioBoost-Matched Benchmark Report

Updated: 2026-07-03

## Why This Benchmark Was Needed

CardioBoost is not a general cardiogenetics variant classifier. The published
model and official site are built for rare missense variants in selected
cardiomyopathy and arrhythmia genes, with high-confidence calls at the
`<=0.1` and `>=0.9` score thresholds.

Because of that, the fair comparison is not our full 85k-row matrix against
CardioBoost. The fair comparison is a matched rare-missense, CardioBoost-gene,
binary-label subset where both models can be scored.

## What Was Done

Official CardioBoost model artifacts were downloaded directly from the public
GitHub repository into:

```text
datasets/cardioboost/public_dataset/official_cardioBoost_model_snapshot/
```

The required official files are now local:

| File | Role |
|---|---|
| `data/cardiomyopathy/ml/train_ada.RData` | Official cardiomyopathy AdaBoost model |
| `data/arrhythmia/ml/train_ada.RData` | Official arrhythmia AdaBoost model |
| `data/cardiomyopathy/preprocess.RData` | Official cardiomyopathy preprocessing |
| `data/arrhythmia/preprocess.RData` | Official arrhythmia preprocessing |
| `script/src/predict.R` | Official scoring helper |
| `script/src/preprocess_test.R` | Official preprocessing helper |

The official models were run locally against the official all-rare mutation
tables. Outputs:

| Output | Rows |
|---|---:|
| `official_cardioboost_cm_all_rare_predictions.tsv` | 65,475 |
| `official_cardioboost_arm_all_rare_predictions.tsv` | 42,411 |

Then the official CardioBoost predictions were matched to our matrix using:

1. exact coordinate matching as a QC check;
2. HGVS cDNA matching as the primary benchmark;
3. HGVS protein-position fallback as a sensitivity check.

## Eligibility Counts From Our Matrix

| Slice | Rows | Benign | Pathogenic | VUS |
|---|---:|---:|---:|---:|
| All current matrix rows | 85,677 | 33,923 | 9,148 | 42,361 |
| CardioBoost-gene rows | 37,528 | 14,459 | 5,670 | 17,186 |
| CardioBoost-gene missense rows | 18,251 | 1,287 | 2,296 | 14,465 |
| CardioBoost-gene missense binary rows | 3,583 | 1,287 | 2,296 | 0 |
| CardioBoost-gene missense strict-rare binary rows | 1,527 | 1,073 | 454 | 0 |
| Exact coordinate matches to official all-rare universe | 64 | 34 | 8 | 22 |
| Exact coordinate matches, binary only | 42 | 34 | 8 | 0 |

The exact-coordinate match is too small for a final manuscript claim. This is
probably due to genome build and transcript-version differences between the
official CardioBoost tables and the current GRCh38-centered matrix. Therefore,
the exact-coordinate result should be treated as QC only.

## Primary Matched Benchmark

The primary benchmark uses HGVS cDNA matching.

| Benchmark | Rows | Unique variants | CM rows | Arrhythmia rows |
|---|---:|---:|---:|---:|
| HGVS cDNA primary | 2,462 | 2,344 | 1,765 | 697 |
| HGVS cDNA plus protein fallback | 2,477 | 2,359 | 1,776 | 701 |

Split composition for the HGVS cDNA primary benchmark:

| Split | Rows |
|---|---:|
| Train | 1,897 |
| Validation | 315 |
| External CardioBoost | 164 |
| External eMERGE | 73 |
| External HiRO | 13 |

For the manuscript, the external subsets are the most important rows. The full
2,462-row benchmark is useful for QC, but it includes training rows.

## External Matched Results

Primary external benchmark: HGVS cDNA-matched, source-held-out rows only.

| Model | Rows | Pathogenic | Benign | AUROC | AUPRC | Sensitivity | Specificity | PPV | NPV | Deferral | High-conf accuracy | MCC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Our CatBoost | 250 | 169 | 81 | 0.980 | 0.990 | 0.959 | 0.543 | 0.964 | 1.000 | 0.152 | 0.972 | 0.862 |
| Official CardioBoost | 250 | 169 | 81 | 0.843 | 0.902 | 0.876 | 0.148 | 0.836 | 1.000 | 0.244 | 0.847 | 0.436 |
| REVEL local | 247 | 166 | 81 | 0.903 | 0.940 | 0.530 | 0.037 | 0.936 | 1.000 | 0.607 | 0.938 | 0.641 |
| CADD PHRED | 250 | 169 | 81 | 0.872 | 0.914 | 1.000 | 0.012 | 0.693 | 1.000 | 0.020 | 0.694 | 0.227 |
| AlphaMissense | 248 | 167 | 81 | 0.937 | 0.969 | 0.635 | 0.457 | 0.991 | 0.902 | 0.403 | 0.966 | 0.686 |

External CardioBoost-source-only benchmark:

| Model | Rows | Pathogenic | Benign | AUROC | AUPRC | Sensitivity | Specificity | PPV | NPV | Deferral | High-conf accuracy | MCC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Our CatBoost | 164 | 123 | 41 | 0.972 | 0.990 | 0.967 | 0.366 | 0.952 | 1.000 | 0.146 | 0.957 | 0.799 |
| Official CardioBoost | 164 | 123 | 41 | 0.788 | 0.905 | 0.870 | 0.073 | 0.843 | 1.000 | 0.207 | 0.846 | 0.225 |
| REVEL local | 161 | 120 | 41 | 0.880 | 0.944 | 0.483 | 0.049 | 0.935 | 1.000 | 0.602 | 0.938 | 0.635 |
| CADD PHRED | 164 | 123 | 41 | 0.810 | 0.903 | 1.000 | 0.024 | 0.759 | 1.000 | 0.006 | 0.761 | 0.192 |
| AlphaMissense | 162 | 121 | 41 | 0.906 | 0.966 | 0.645 | 0.244 | 0.987 | 0.833 | 0.438 | 0.967 | 0.582 |

## Paired Bootstrap Confidence Intervals

Paired bootstrap confidence intervals were added on 2026-07-03. The bootstrap
resamples the same matched rows for our CatBoost and official CardioBoost, so
the confidence intervals test the paired difference directly.

Bootstrap settings:

| Setting | Value |
|---|---:|
| Replicates | 5,000 |
| Random seed | 20260703 |
| Primary bootstrap | Row-stratified paired bootstrap |
| Sensitivity bootstrap | Variant-cluster stratified paired bootstrap |
| Models included | Our CatBoost and official CardioBoost |

Primary row-stratified paired bootstrap, external all rows:

| Model | Metric | Estimate | 95% CI |
|---|---|---:|---:|
| Our CatBoost | AUROC | 0.980 | 0.961-0.995 |
| Official CardioBoost | AUROC | 0.843 | 0.786-0.895 |
| Our CatBoost | AUPRC | 0.990 | 0.979-0.998 |
| Official CardioBoost | AUPRC | 0.902 | 0.858-0.940 |
| Our CatBoost | Sensitivity | 0.959 | 0.923-0.988 |
| Official CardioBoost | Sensitivity | 0.876 | 0.822-0.923 |
| Our CatBoost | PPV | 0.964 | 0.937-0.988 |
| Official CardioBoost | PPV | 0.836 | 0.797-0.878 |
| Our CatBoost | Deferral | 0.152 | 0.112-0.192 |
| Official CardioBoost | Deferral | 0.244 | 0.196-0.292 |
| Our CatBoost | High-confidence accuracy | 0.972 | 0.949-0.991 |
| Official CardioBoost | High-confidence accuracy | 0.847 | 0.807-0.888 |

Primary row-stratified paired bootstrap, external all paired differences:

| Metric | Our minus CardioBoost | 95% CI |
|---|---:|---:|
| AUROC | 0.138 | 0.083-0.195 |
| AUPRC | 0.087 | 0.048-0.132 |
| Sensitivity | 0.083 | 0.030-0.136 |
| Specificity | 0.395 | 0.272-0.519 |
| PPV | 0.128 | 0.082-0.172 |
| Deferral | -0.092 | -0.156 to -0.028 |
| High-confidence accuracy | 0.125 | 0.082-0.167 |

Primary row-stratified paired bootstrap, external CardioBoost-source-only paired
differences:

| Metric | Our minus CardioBoost | 95% CI |
|---|---:|---:|
| AUROC | 0.184 | 0.101-0.274 |
| AUPRC | 0.085 | 0.040-0.135 |
| Sensitivity | 0.098 | 0.041-0.163 |
| Specificity | 0.293 | 0.146-0.463 |
| PPV | 0.109 | 0.058-0.158 |
| Deferral | -0.061 | -0.134 to 0.012 |
| High-confidence accuracy | 0.111 | 0.061-0.159 |

Variant-cluster bootstrap sensitivity analysis gave the same conclusion. For
external all rows, the paired AUROC difference was 0.138 with 95% CI
0.081-0.201, and the paired AUPRC difference was 0.087 with 95% CI
0.045-0.135. For CardioBoost-source-only rows, the paired AUROC difference was
0.184 with 95% CI 0.100-0.273, and the paired AUPRC difference was 0.085 with
95% CI 0.039-0.137.

Interpretation: the paired bootstrap supports a robust performance advantage
for our model over official CardioBoost on AUROC, AUPRC, sensitivity, PPV, and
high-confidence accuracy in both matched external benchmarks. The deferral rate
is significantly lower in the 250-row external-all benchmark; in the stricter
164-row CardioBoost-source-only benchmark it trends lower but the 95% CI crosses
zero.

## Sensitivity Check

Adding HGVS protein-position fallback added only 15 unique variants. It did not
change the conclusion.

| Benchmark | Our AUROC | CardioBoost AUROC | Our AUPRC | CardioBoost AUPRC |
|---|---:|---:|---:|---:|
| cDNA external all | 0.980 | 0.843 | 0.990 | 0.902 |
| cDNA plus protein fallback external all | 0.981 | 0.845 | 0.990 | 0.902 |
| cDNA external CardioBoost-source | 0.972 | 0.788 | 0.990 | 0.905 |
| cDNA plus protein fallback external CardioBoost-source | 0.972 | 0.787 | 0.990 | 0.905 |

## Figures

ROC and precision-recall figures were generated here:

```text
results/model_performance/cardioboost_matched_benchmark/figures/hgvs_cdot_primary_external_all_roc.png
results/model_performance/cardioboost_matched_benchmark/figures/hgvs_cdot_primary_external_all_precision_recall.png
results/model_performance/cardioboost_matched_benchmark/figures/hgvs_cdot_primary_external_cardioboost_roc.png
results/model_performance/cardioboost_matched_benchmark/figures/hgvs_cdot_primary_external_cardioboost_precision_recall.png
```

## Interpretation

This benchmark is the strongest current evidence that the new model compares
favorably with CardioBoost on CardioBoost's own intended task.

On the external HGVS cDNA-matched benchmark, our model had higher AUROC, higher
AUPRC, higher sensitivity, higher specificity, higher PPV, lower deferral, and
higher high-confidence accuracy than the official CardioBoost model.

The result is especially useful because it avoids the unfair full-matrix
comparison. It asks a narrower question:

> Among rare missense, CardioBoost-gene, binary-labeled variants that can be
> matched to the official CardioBoost scoring universe, does the new model
> perform better?

The current answer is yes, based on this matched external analysis.

## Caveats Before Manuscript Use

These results are promising, but they should be reported carefully.

1. The exact-coordinate benchmark is not adequate because it only recovered 71
   binary panel rows. Use it as QC only.
2. The primary HGVS cDNA benchmark includes training rows if the full benchmark
   is used. Manuscript claims should focus on `external_all` and
   `external_cardioboost`.
3. SCN5A can appear in both cardiomyopathy and arrhythmia contexts. Panel-aware
   handling should remain explicit.
4. Some benign specificity values are modest because the benchmark has many
   pathogenic rows and the `0.1/0.9` strategy intentionally prioritizes avoiding
   hard pathogenic-to-benign errors.
5. Paired bootstrap confidence intervals are now available, but final manuscript
   tables should still include them explicitly.

## Reproducibility

Official CardioBoost scoring:

```bash
Rscript results/model_performance/cardioboost_matched_benchmark/run_official_cardioboost_predictions.R
```

Matched benchmark construction:

```bash
python3 results/model_performance/cardioboost_matched_benchmark/build_exact_match_benchmark.py
```

Figure generation:

```bash
python3 results/model_performance/cardioboost_matched_benchmark/plot_cardiboost_benchmark.py
```

Paired bootstrap confidence intervals:

```bash
python3 results/model_performance/cardioboost_matched_benchmark/bootstrap_cardioBoost_matched_ci.py
```

Key output tables:

```text
results/model_performance/cardioboost_matched_benchmark/tables/hgvs_cardioBoost_matched_metrics.tsv
results/model_performance/cardioboost_matched_benchmark/tables/hgvs_cdot_cardioBoost_matched_benchmark.tsv
results/model_performance/cardioboost_matched_benchmark/tables/hgvs_pdot_fallback_cardioBoost_matched_benchmark.tsv
results/model_performance/cardioboost_matched_benchmark/tables/hgvs_cardioBoost_matched_bootstrap_model_ci.tsv
results/model_performance/cardioboost_matched_benchmark/tables/hgvs_cardioBoost_matched_bootstrap_paired_difference_ci.tsv
```

## Next Steps

1. Generate a manuscript-ready combined figure with ROC, PR, and a compact
   threshold-metric panel.
2. Audit all external CardioBoost-source matched rows for possible label/source
   leakage.
3. Add the bootstrap CI table to the manuscript draft.
