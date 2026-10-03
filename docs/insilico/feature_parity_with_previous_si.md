# In Silico Feature Parity With Previous SI

The previous SI used the following in silico predictors in the frozen evidence
layer and CatBoost feature table:

- `REVEL`
- `CADD PHRED`
- `SIFT`
- `PolyPhen-2 HDIV`
- `MetaLR`
- `FATHMM-XF`
- `AlphaMissense`

These are the minimum parity features for the current project. The current
project should include all of them, then add the newer/proposed sources:

- dbNSFP aggregate predictors: MutationTaster, MutationAssessor, FATHMM,
  PROVEAN, VEST, M-CAP, MetaSVM, Eigen, GERP++, PhyloP, PhastCons, SiPhy, and
  any relevant newer dbNSFP v5.3.1 predictors.
- Ensembl VEP consequence annotation and SpliceAI.
- CADD.
- AlphaMissense.
- ESM1b/ESM-variant.

## Current Status

Implemented and validated:

- `AlphaMissense`: local hg38 precomputed file downloaded and joined.
- `CADD PHRED`: targeted CADD GRCh38-v1.7 API workflow added and validated on
  the example ClinVar variants.

Pending for old-SI parity:

- `REVEL`
- `SIFT`
- `PolyPhen-2 HDIV`
- `MetaLR`
- `FATHMM-XF`

These should be supplied by dbNSFP. Current dbNSFP v5.3.1 is the preferred
source, but it requires institutional-email registration and an access code.
Do not use guessed S3 paths or third-party mirrors.

Pending new-proposal additions:

- dbNSFP aggregate predictor panel.
- VEP/SpliceAI.
- ESM1b/ESM-variant.

## Automated Check

Use this check on any merged in silico feature table:

```bash
python3 datasets/feature_sources/insilico/scripts/validate_insilico_requirements.py \
  --input datasets/feature_sources/insilico/interim/insilico_features_clinvar_example.tsv
```

Current expected result before dbNSFP is added:

- Present: `AlphaMissense`, `CADD_PHRED`
- Missing: `REVEL`, `SIFT`, `PolyPhen2_HDIV`, `MetaLR`, `FATHMM_XF`

After dbNSFP is added, rerun with `--strict`. The strict check should pass
before the feature table is used for modeling.

## Non-Negotiable QC

- Missing values must be real missing values, not silent API or syntax failures.
- Each source-specific join must include status and missing-reason columns.
- Coordinate-based joins must use GRCh38 unless a liftover step is documented.
- Protein-level joins must include transcript/protein mapping validation.
