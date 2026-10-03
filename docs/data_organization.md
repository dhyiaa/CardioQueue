# Data Organization

This project keeps dataset sources discoverable under `datasets/`. Generated
analysis-ready tables can still go under `data/processed/` later, but the source
catalog lives in one place.

## Layout

- `datasets/`: all local datasets, API-accessible source folders, raw files,
  source-specific scripts, provenance, and source READMEs.
- `datasets/feature_sources/`: feature-injection sources such as gnomAD,
  in silico predictors, ClinGen, AlphaFold, and UniProt.
- `data/interim/`: scratch space only; should stay small.
- `data/processed/`: analysis-ready feature/label tables.
- `metadata/`: dataset manifests, provenance, source URLs, checksums, and usage
  notes.
- `scripts/`: reserved for cross-source project scripts. Source-specific scripts
  should live next to their dataset/source under `datasets/`.

## Current Sources

### Dataset 1: ClinVar

Current local file:

- `datasets/clinvar/variant_summary.txt.gz`

Validation performed on 2026-07-01:

- Source URL: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz`
- Source MD5: `f03eea5e87f0ef5f696bbc958359fa78`
- Local MD5: `f03eea5e87f0ef5f696bbc958359fa78`
- Rows: `8,991,566`
- Columns: `43`

### Dataset 2: gnomAD v4

Primary role:

- Population-frequency features for variants from any source.
- Gene-level constraint features.
- Presumed-benign proxy variants within a reviewed disease/gene panel.

Recommended access pattern:

- Use the gnomAD GraphQL API for targeted variant/gene-panel queries.
- Use bulk VCF/Hail tables only when the final gene panel and extraction logic
  justify the extra storage and compute.

Current API constants:

- API endpoint: `https://gnomad.broadinstitute.org/api`
- Dataset enum: `gnomad_r4`
- Reference genome: `GRCh38`
