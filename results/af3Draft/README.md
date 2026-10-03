# AF3Draft Results

This folder stores prototype AlphaFold 3-derived outputs for the cardiogenetics project.

The `af3Draft` name is intentional. These files are not the final AF3 feature family. They are used to test whether the AF3 pathway can be integrated into the existing CatBoost model without changing the established data/modeling pipeline.

Current smoke-test report:

```text
results/model_performance/af3Draft_smoke_test/AF3DRAFT_SMOKE_TEST_REPORT.md
```

Current aggregated AF3Draft feature table:

```text
datasets/feature_sources/af3Draft/processed/af3Draft_variant_features_aggregated.tsv
```

Current augmented modeling table:

```text
datasets/modeling/interim/modeling_table_with_splits_weights_plus_af3Draft.tsv
```

Use the strict structural-only AF3Draft model as the cleaner smoke-test result. The full draft model includes provenance-like fields and should be treated as exploratory only.
