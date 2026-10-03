# Point D: ClinGen Curated Evidence Plan

This document starts Point D as ClinGen curated evidence, inferred from the
current checked-in proposal. The local proposal names ClinGen Evidence Repository
as the better source for expert-curated variant classifications and evidence,
but it does not include the full original A/B/C/D source list. If Point D was
intended to be a different source, this file should be replaced before pipeline
implementation.

## Proposed Role

ClinGen should not replace ClinVar as the broad primary label source. Instead,
it should provide higher-confidence curated evidence and expert-context features
where available.

| ClinGen source area | Proposed role | Modeling use |
|---|---|---|
| ClinGen Evidence Repository / Variant Pathogenicity curations | Expert-curated variant-level classifications and evidence summaries | High-confidence metadata, label-confidence flags, possible external validation subset |
| ClinGen gene-disease validity | Gene-disease relationship strength and disease context | Gene-level features and disease-context filters |
| ClinGen dosage sensitivity | Haploinsufficiency/triplosensitivity context where relevant | Gene mechanism/context features, especially for loss-of-function interpretation |
| ClinGen VCEP specifications | Gene/disease-specific ACMG/AMP rule modifications | WE-LLM evidence rules and deterministic criterion-checking support |

## Current Access Notes

Official public ClinGen access points checked:

| Resource | URL | Observed access |
|---|---|---|
| Evidence Repository browse | `https://erepo.clinicalgenome.org/evrepo/` | Public web interface responded successfully |
| ClinGen downloads/API page | `https://search.clinicalgenome.org/kb/downloads` | Public downloads/API page responded successfully |
| Gene-disease validity browser/API-backed page | `https://search.clinicalgenome.org/kb/gene-validity` | Public browser page responded successfully and references `/api/validity` |
| Variant pathogenicity browser | `https://search.clinicalgenome.org/kb/variant-pathogenicity/all` | Linked from ClinGen downloads/navigation; needs targeted API/download check next |

Targeted implementation added:

| Script | Scope | Status |
|---|---|---|
| `datasets/feature_sources/clingen/scripts/clingen_api.py doctor` | Checks public ClinGen gene-disease validity JSON shape and gene lookup response | Added and validated |
| `datasets/feature_sources/clingen/scripts/clingen_api.py gene-validity` | Annotates a gene panel with ClinGen gene-disease validity rows | Added and validated on the 19-gene cardiogenetics seed panel |
| `datasets/feature_sources/clingen/scripts/clingen_api.py validate` | Validates required status/missing-reason fields | Added and passed |

The current script intentionally covers gene-disease validity only. Variant-level
ClinGen Evidence Repository support remains pending until an exact official
download/API endpoint and schema are confirmed.

## Candidate Features

| Feature group | Candidate fields |
|---|---|
| Variant-level ClinGen evidence | ClinGen assertion/classification, VCEP/expert panel, disease/condition, assertion date, review status, evidence summary availability |
| Gene-disease validity | gene symbol, HGNC ID, disease label, MONDO ID, inheritance, classification, expert panel, report date |
| Dosage sensitivity | haploinsufficiency score, triplosensitivity score, disease context, curation date |
| Rule guidance | applicable VCEP specification, gene/disease-specific ACMG rules, PVS1/gene mechanism notes |

## Matching Strategy

Preferred matching order:

1. ClinGen Allele Registry CA ID, if available.
2. Exact GRCh38 genomic coordinate plus reference and alternate allele.
3. Exact HGVS cDNA/protein with transcript reconciliation.
4. rsID only as a secondary cross-check, not as the primary exact-variant key.

For gene-disease validity and dosage sensitivity, match by HGNC ID when possible
and preserve disease identifiers separately. Do not treat same-gene evidence as
exact-variant evidence.

## Null and Failure Handling

| Situation | Required status |
|---|---|
| Exact variant curation found | `clingen_variant_status=ok` |
| Gene-level curation found but no exact variant curation | `clingen_variant_status=no_exact_variant_match`; preserve gene-level fields separately |
| No relevant ClinGen record found | `clingen_variant_status=not_found` |
| Coordinates/HGVS insufficient for lookup | `clingen_variant_status=not_queryable` |
| API/download/schema failure | fail the run or mark explicit error status; do not create silent nulls |

## Implementation Plan

1. Confirm exact Point D source wording from the proposal.
2. Identify official downloadable ClinGen files and API endpoints for variant
   pathogenicity and dosage sensitivity.
3. Extend `datasets/feature_sources/clingen/scripts/clingen_api.py` only after endpoints and schemas are
   verified.
4. Build a small variant-panel test query after variant-level access is fixed.
5. Save raw JSON/TSV responses under `datasets/feature_sources/clingen/raw/` or
   `datasets/feature_sources/clingen/interim/` with retrieval dates.
6. Validate that no failed API call is represented as a simple missing value.
7. Join ClinGen fields into the frozen evidence layer as metadata/features, not
   as direct replacement labels unless a separate high-confidence validation
   analysis is explicitly defined.

## Current Status

Point D is documented and partially scaffolded. Gene-disease validity annotation
is implemented as a targeted API workflow. The example cardiogenetics seed panel
returned 65 ClinGen gene-disease validity rows across 19 input genes, all with
`clingen_gene_validity_status=ok`. Variant-level ClinGen evidence is not
implemented yet.

Example workflow:

```bash
python3 datasets/feature_sources/clingen/scripts/clingen_api.py doctor

python3 datasets/feature_sources/clingen/scripts/clingen_api.py gene-validity \
  --gene-panel datasets/feature_sources/gnomad_v4/external/gene_panel.example.txt \
  --output datasets/feature_sources/clingen/interim/gene_validity.example.tsv \
  --raw-json datasets/feature_sources/clingen/raw/gene_validity.full_response.example.json

python3 datasets/feature_sources/clingen/scripts/clingen_api.py validate \
  --input datasets/feature_sources/clingen/interim/gene_validity.example.tsv
```

Example outputs:

| Output | Path | Status |
|---|---|---|
| Raw ClinGen validity API response | `datasets/feature_sources/clingen/raw/gene_validity.full_response.example.json` | Saved for reproducibility; 2.2 MB |
| Gene-panel ClinGen validity table | `datasets/feature_sources/clingen/interim/gene_validity.example.tsv` | Validated; 65 rows |
