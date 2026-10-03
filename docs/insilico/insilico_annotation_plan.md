# In Silico Annotation Plan

In silico tools are feature sources only. They are not independent clinical
labels and should not be used to define P/LP, VUS, or B/LB classes.

For a compact source-by-source status table, see
`docs/insilico/point_c_status_summary.md`.

## Priority Order

1. AlphaMissense hg38 precomputed scores: local file downloaded and validated.
2. dbNSFP: high-value aggregate source, but current v5.3.1 academic downloads
   require institutional-email registration and an access code.
3. CADD: precomputed GRCh38 files are large; use targeted indexed lookup or
   download only the required subset.
4. Ensembl VEP: use local cache only when storage allows; current human GRCh38
   cache is large.
5. ESM1b/ESM-variant: add only after the exact precomputed human proteome source
   and coordinate/protein mapping strategy are fixed.

## Required Feature Parity

The previous SI in silico panel requires:

- `REVEL`
- `CADD PHRED`
- `SIFT`
- `PolyPhen-2 HDIV`
- `MetaLR`
- `FATHMM-XF`
- `AlphaMissense`

Current status:

- Implemented: `AlphaMissense`, targeted `CADD PHRED`.
- Pending through dbNSFP or another verified source: `REVEL`, `SIFT`,
  `PolyPhen-2 HDIV`, `MetaLR`, `FATHMM-XF`.

The current project must include these old-SI predictors plus the new proposal
additions documented in `metadata/insilico_feature_requirements.json`.

## Local AlphaMissense Workflow

```bash
python3 datasets/feature_sources/insilico/scripts/join_alphamissense.py join \
  --input datasets/feature_sources/gnomad_v4/external/clinvar_grch38_variants.example.tsv \
  --output datasets/feature_sources/insilico/interim/alphamissense_clinvar_example.tsv

python3 datasets/feature_sources/insilico/scripts/join_alphamissense.py validate \
  --input datasets/feature_sources/insilico/interim/alphamissense_clinvar_example.tsv
```

Null handling:

- `alphamissense_status=ok`: exact hg38 `chrom-pos-ref-alt` match found.
- `alphamissense_status=not_applicable_or_not_found`: no exact AlphaMissense
  record. This is expected for non-missense variants, non-SNVs, and variants not
  represented in AlphaMissense.
- Rows marked `ok` must have `am_pathogenicity`; validation fails otherwise.

## Targeted CADD Workflow

The full CADD v1.7 GRCh38 SNV file is too large for this stage. Use the CADD API
only for small QC/test sets, because CADD explicitly says the API is experimental
and not intended for thousands or millions of variants.

```bash
python3 datasets/feature_sources/insilico/scripts/cadd_api.py doctor

python3 datasets/feature_sources/insilico/scripts/cadd_api.py annotate \
  --input datasets/feature_sources/gnomad_v4/external/clinvar_grch38_variants.example.tsv \
  --output datasets/feature_sources/insilico/interim/cadd_clinvar_example.tsv

python3 datasets/feature_sources/insilico/scripts/cadd_api.py validate \
  --input datasets/feature_sources/insilico/interim/cadd_clinvar_example.tsv
```

Null handling:

- `cadd_status=ok`: score returned with `cadd_rawscore` and `cadd_phred`.
- `cadd_status=not_found`: API returned an empty list.
- `cadd_status=not_applicable`: non-SNV for this targeted API endpoint.
- Unexpected API/JSON/network failures fail the run.

## Merge Current In Silico Features

```bash
python3 datasets/feature_sources/insilico/scripts/merge_insilico_features.py merge \
  --alphamissense datasets/feature_sources/insilico/interim/alphamissense_clinvar_example.tsv \
  --cadd datasets/feature_sources/insilico/interim/cadd_clinvar_example.tsv \
  --output datasets/feature_sources/insilico/interim/insilico_features_clinvar_example.tsv

python3 datasets/feature_sources/insilico/scripts/merge_insilico_features.py validate \
  --input datasets/feature_sources/insilico/interim/insilico_features_clinvar_example.tsv

python3 datasets/feature_sources/insilico/scripts/validate_insilico_requirements.py \
  --input datasets/feature_sources/insilico/interim/insilico_features_clinvar_example.tsv
```

## Current Size Decisions

- AlphaMissense hg38: 642,961,469 bytes compressed. Downloaded.
- CADD v1.7 GRCh38 all possible SNVs: 81G compressed. Do not download blindly.
- CADD v1.7 GRCh38 gnomAD r4.0 InDels: 1,257,151,321 bytes compressed. Feasible
  but defer until indel feature strategy is fixed.
- Ensembl VEP 116 human GRCh38 cache: 27,644,657,162 bytes compressed. Defer
  unless local VEP is chosen over REST/API or smaller targeted annotation.
- dbNSFP v5.3.1 academic branch: about 50 GB and requires institutional-email
  registration plus an access code before download links are issued.

## Required QC Before Modeling

- Every output table must include a status column and missing-reason column.
- Failed API/download/schema problems must fail the run, not create blank cells.
- Joins must use GRCh38 coordinates unless explicitly documented otherwise.
- Protein-level tools such as ESM1b require transcript/protein mapping QC before
  features are trusted.
