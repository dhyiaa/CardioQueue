# CardioBoost-Matched Comparison Plan

## Completion Status

Completed on 2026-07-02. The official CardioBoost models were run locally, a
matched HGVS cDNA benchmark was built, ROC/PR figures were generated, and the
main result report is here:

```text
results/model_performance/cardioboost_matched_benchmark/CARDIOBOOST_MATCHED_BENCHMARK_REPORT.md
```

The remaining manuscript-readiness tasks are paired bootstrap confidence
intervals, leakage/source-overlap audit, and one combined publication figure.

## Current Finding

The comparison should be improved by separating two claims:

1. **Head-to-head performance where CardioBoost is eligible.**
2. **Callability and scope advantage across the full cardiogenetics matrix.**

CardioBoost is not designed for all variants. The official manuscript and web
site restrict it to rare missense variants in selected cardiomyopathy and
arrhythmia genes, using canonical transcripts and `gnomAD_AF <= 0.1%`.

Therefore, the clean one-vs-one benchmark should not use all `86,889` rows. It
should use the subset of our data that CardioBoost can reasonably score.

## Official CardioBoost Artifacts Now Local

The official GitHub repository is reachable at commit:

```text
6d5783bf975e20e39f3f6833f16b2edb86cb3525
```

The full `git clone` disconnected, so the required official files were
downloaded directly into:

```text
datasets/cardioboost/public_dataset/official_cardioBoost_model_snapshot/
```

Current local official-model files:

| File | Role |
|---|---|
| `data/cardiomyopathy/ml/train_ada.RData` | Official cardiomyopathy AdaBoost model |
| `data/arrhythmia/ml/train_ada.RData` | Official arrhythmia AdaBoost model |
| `data/cardiomyopathy/preprocess.RData` | Official cardiomyopathy preprocessing parameters |
| `data/arrhythmia/preprocess.RData` | Official arrhythmia preprocessing parameters |
| `script/src/predict.R` | Official prediction helper |
| `script/src/preprocess_test.R` | Official preprocessing helper |
| `script/cardiomyopathy/ml/cm_prediction_all_possible_variants.R` | Official CM all-rare prediction script |
| `script/arrhythmia/ml/arm_prediction_all_possible_variants.R` | Official ARM all-rare prediction script |
| `script/requiredRpkgs.txt` | Original R package list |

The official model objects load, but prediction needs old R packages that are
not currently installed locally:

| R package | Installed locally |
|---|---|
| `mlr` | no |
| `ada` | no |
| `data.table` | no |
| `pROC` | no |
| `PRROC` | no |
| `precrec` | no |

## Current Eligibility Counts From Our Matrix

Input matrix:

```text
datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv
```

Initial benchmark-count output:

```text
results/model_performance/cardioboost_matched_benchmark/tables/cardioboost_benchmark_eligibility_counts.tsv
```

| Slice | Rows | Benign | Pathogenic | VUS | Interpretation |
|---|---:|---:|---:|---:|---|
| All corrected matrix rows | 85,677 | 33,923 | 9,148 | 42,361 | Full current model scope |
| CardioBoost gene rows | 37,528 | 14,459 | 5,670 | 17,186 | Rows in the 22 unique CardioBoost genes |
| CardioBoost gene missense rows | 18,251 | 1,287 | 2,296 | 14,465 | Broad missense slice |
| CardioBoost gene missense binary rows | 3,583 | 1,287 | 2,296 | 0 | Potential binary head-to-head pool before strict rarity/matching |
| CardioBoost gene missense strict-rare binary rows | 1,527 | 1,073 | 454 | 0 | Better matched to CardioBoost's rare-missense scope |
| Exact coordinate matches to official all-rare universe | 64 | 34 | 8 | 22 | Too low for final use; likely build/transcript mismatch |
| Exact coordinate matches, binary only | 42 | 34 | 8 | 0 | Not enough; needs liftover/HGVS matching |

The exact coordinate match to the official all-rare CardioBoost universe is too
small to use as the final benchmark. CardioBoost's original files are likely in
the original publication coordinate/transcript convention, while the current
matrix is GRCh38-centered. The next step is not to accept the 64-row match. The
next step is to build a proper CardioBoost matching layer.

## Best Benchmark Design

### Benchmark 1: Official CardioBoost-Eligible Head-To-Head

Use variants that satisfy:

- gene is in the official CardioBoost cardiomyopathy or arrhythmia list;
- VEP consequence is missense;
- binary label is Benign or Pathogenic;
- gnomAD AF/popmax is <=0.001 or confirmed absent;
- variant can be mapped to CardioBoost's canonical transcript/HGVSc/HGVSp or
  official all-rare universe.

For every matched variant, compute:

| Score | Source |
|---|---|
| Our CatBoost P(pathogenic) | Current primary binary model |
| CardioBoost pathogenicity | Official `train_ada.RData` model or released `pathogenicity` column |
| REVEL | Local dbNSFP-selected feature |
| CADD PHRED | Local dbNSFP-selected feature |
| AlphaMissense | Direct join preferred, dbNSFP fallback |

Report:

- AUROC;
- AUPRC;
- sensitivity at `0.1/0.9`;
- specificity at `0.1/0.9`;
- PPV;
- NPV;
- high-confidence classified rate;
- high-confidence accuracy;
- deferral rate;
- paired bootstrap confidence intervals.

### Benchmark 2: Source-Held-Out CardioBoost Public Dataset

This is already partly done. The current source-held-out model was tested on
CardioBoost rows that were not used for training.

Current binary + `0.1/0.9` CardioBoost-source result:

| Rows | Sensitivity | Specificity | PPV | NPV | Deferral |
|---:|---:|---:|---:|---:|---:|
| 185 | 0.932 | 0.887 | 0.953 | 0.903 | 0.135 |

This is valid external validation, but it is not the strictest published
CardioBoost head-to-head because it is not yet filtered to the exact
rare-missense benchmark and does not yet compare against CardioBoost's own score
on the same rows.

### Benchmark 3: Callability / Scope Advantage

This is separate from performance.

CardioBoost can only score rare missense variants in its supported genes and
canonical transcripts. Our model can score a broader cardiogenetics matrix,
including non-CardioBoost genes, splice variants, truncating variants, indels,
and VUS rows for prioritization.

Report:

| Dataset | Rows | Our model score available | CardioBoost eligible |
|---|---:|---:|---:|
| Full matrix | 85,677 | nearly all modelable rows | limited to rare missense in supported genes |
| Binary supervised rows | 42,990 | yes | subset only |
| VUS scoring rows | 41,731 | yes, as triage/defer score | only if rare missense and supported |
| CardioBoost-gene strict-rare binary missense | 1,527 | yes | likely yes after transcript/build matching |

This should be framed as broader applicability, not better accuracy.

## Official Web Site Role

The official CardioBoost web site is useful for spot checks and manual
validation. It is not the best primary benchmark engine because:

- it is a Shiny app;
- it is less reproducible than local model files;
- large batch automation may be fragile;
- the GitHub model/data are the publication's reproducible source.

Use the web site only to verify a handful of variants after the local official
model pipeline is working.

## Immediate Next Steps

1. Install a small R compatibility environment for `mlr` and `ada`.
2. Run official CardioBoost prediction scripts on `cm_all_rare_mutation.RData`
   and `arm_all_rare_mutation.RData`.
3. Create a transcript-aware match from our variants to CardioBoost predictions:
   coordinate liftover plus HGVSc/HGVSp fallback.
4. Build the final matched rare-missense binary benchmark table.
5. Score the same rows with our primary binary CatBoost model.
6. Add REVEL, CADD, and AlphaMissense standalone comparator columns.
7. Generate one manuscript-ready table and one ROC/PR figure.
