# gnomAD v4 Dataset Plan

## Purpose

gnomAD v4 will support two parts of the study:

1. Population-frequency features for variants from ClinVar, ClinGen, or curated
   sources.
2. A presumed-benign proxy class from variants observed in the general population
   at frequencies inconsistent with severe Mendelian cardiogenetic disease.

## Current Decision

Do not download the full raw gnomAD v4 release yet. The current project need is
a clean cardiogenetics registry, not the entire population database. The
near-term strategy is:

1. Use local dbNSFP-derived gnomAD fields immediately.
2. Use direct gnomAD API only for small pilots, gene constraints, or spot checks.
3. Escalate to browser Hail Table, VCF shards, cloud query, or remote tabix only
   if we need fields that dbNSFP does not provide.

Current local gnomAD proxy table:

- `datasets/feature_sources/gnomad_v4/interim/gnomad_from_dbnsfp.clean_combined.tsv`
- 86,886 clean-registry rows.
- 42,338 rows with local gnomAD population-frequency proxy fields.
- 44,548 rows without local dbNSFP gnomAD fields.

Current direct gnomAD API outputs:

- `datasets/feature_sources/gnomad_v4/interim/gnomad_v4_gene_constraints.clean_primary.tsv`
  has 53/53 primary-panel genes found.
- `datasets/feature_sources/gnomad_v4/interim/gnomad_v4_variant_annotations.clean_registry_pilot200.tsv`
  has 144/200 variants found.

The 200-variant API pilot validates the API workflow but is too slow to use as
a one-by-one full-registry annotation strategy.

## Browser Hail Table Escalation

The direct gnomAD browser Hail Table route is now prepared for variants missing
local dbNSFP-derived gnomAD fields.

Official browser Hail Table:

```text
gs://gcp-public-data--gnomad/release/4.1/ht/browser/gnomad.browser.v4.1.sites.ht
```

Prepared local inputs:

- `datasets/feature_sources/gnomad_v4/external/gnomad_missing_from_dbnsfp.clean_combined.tsv`
- `datasets/feature_sources/gnomad_v4/external/gnomad_missing_from_dbnsfp.pilot100.tsv`

Prepared extraction script:

- `datasets/feature_sources/gnomad_v4/scripts/extract_gnomad_browser_hail_missing.py`

Current missing set:

| Group | Count |
|---|---:|
| Variants missing local dbNSFP-derived gnomAD fields | 44,548 |
| ClinVar-only | 41,771 |
| eMERGE-only | 2,740 |
| CardioBoost-only | 18 |
| HiRO-only | 10 |
| ClinVar+HiRO | 6 |
| CardioBoost+ClinVar | 3 |

Installed local tools:

- Java: Temurin 11.0.31 installed at `~/opt/jdk-11.0.31+11/Contents/Home`.
- Python/Hail: conda environment `hail-gnomad` with Python 3.11 and Hail
  0.2.138.
- Google Cloud CLI: `gcloud` 575.0.0 and `gsutil` 5.37.
- Google Cloud Storage connector:
  `datasets/feature_sources/gnomad_v4/raw/browser_hail/jars/gcs-connector-hadoop3-2.2.31-shaded.jar`.
- Local service account key for Spark/Hadoop GCS access:
  `~/.config/gcloud/gnomad-hail-local-key.json` (do not commit).

Authentication status:

- `gcloud auth login` and application-default credentials are configured for
  project `cardiogenetics-gnomad`.
- Service account `gnomad-hail-local@cardiogenetics-gnomad.iam.gserviceaccount.com`
  was created for local Spark access.

100-variant pilot result:

| Output | Rows | `ok` | `not_found` |
|---|---:|---:|---:|
| `datasets/feature_sources/gnomad_v4/interim/gnomad_browser_v4_1_missing_pilot100.tsv.bgz` | 100 | 81 | 19 |

Run command:

```bash
activate-hail-gnomad

python datasets/feature_sources/gnomad_v4/scripts/extract_gnomad_browser_hail_missing.py \
  --input datasets/feature_sources/gnomad_v4/external/gnomad_missing_from_dbnsfp.pilot100.tsv \
  --output datasets/feature_sources/gnomad_v4/interim/gnomad_browser_v4_1_missing_pilot100.tsv.bgz \
  --gcp-project YOUR_GCP_BILLING_PROJECT
```

## Recommended API Workflow

Use the API before bulk downloads:

```bash
python3 datasets/feature_sources/gnomad_v4/scripts/gnomad_api.py gene-constraints \
  --input datasets/feature_sources/gnomad_v4/external/gene_panel.example.txt \
  --output datasets/feature_sources/gnomad_v4/interim/gnomad_v4_gene_constraints.example.tsv

python3 datasets/feature_sources/gnomad_v4/scripts/gnomad_api.py validate-table \
  --input datasets/feature_sources/gnomad_v4/interim/gnomad_v4_gene_constraints.example.tsv \
  --kind genes
```

For variant annotation, provide a TSV/CSV with either `variant_id` or
`chrom,pos,ref,alt` columns:

```bash
python3 datasets/feature_sources/gnomad_v4/scripts/extract_clinvar_variants_for_gnomad.py \
  --clinvar datasets/clinvar/variant_summary.txt.gz \
  --gene-panel datasets/feature_sources/gnomad_v4/external/gene_panel.example.txt \
  --exclude-no-assertion \
  --limit 100 \
  --output datasets/feature_sources/gnomad_v4/external/clinvar_grch38_variants.example.tsv

python3 datasets/feature_sources/gnomad_v4/scripts/gnomad_api.py annotate-variants \
  --input datasets/feature_sources/gnomad_v4/external/clinvar_grch38_variants.example.tsv \
  --output datasets/feature_sources/gnomad_v4/interim/gnomad_v4_variant_annotations.tsv \
  --raw-json-dir datasets/feature_sources/gnomad_v4/interim/raw_json

python3 datasets/feature_sources/gnomad_v4/scripts/gnomad_api.py validate-table \
  --input datasets/feature_sources/gnomad_v4/interim/gnomad_v4_variant_annotations.tsv \
  --kind variants \
  --max-not-found-fraction 0.50
```

The API can rate-limit bursts of requests. Use `--sleep` and small batches while
testing, then cache outputs before scaling up.

Before a run, check that the API query syntax still matches the live gnomAD
schema:

```bash
python3 datasets/feature_sources/gnomad_v4/scripts/gnomad_api.py doctor
```

## Null and Missing-Value Rules

Null values are allowed only when their reason is explicit:

- `query_status=ok`: the API returned a variant/gene record. Missing feature
  fields should be interpreted by field context, for example absent filtering
  allele-frequency values or no genome block.
- `query_status=not_found`: the API explicitly returned a not-found response.
  This is retained as a row with `found=False`.
- Any schema error, syntax error, HTTP failure after retries, or unexpected API
  error fails the run instead of writing silent nulls.

Variant outputs include:

- `has_exome`
- `has_genome`
- `has_joint`
- `has_transcript_consequence`
- `missing_reason`

Use these columns during QC before treating blank fields as biological absence.

## Features

Variant-level features:

- Exome/genome allele count, allele number, and allele frequency.
- Joint allele count and allele number where available.
- Homozygote and hemizygote counts.
- Population maximum AF calculated from ancestry populations.
- gnomAD filtering allele frequency fields where available.
- gnomAD variant filters and flags.
- VEP consequence fields from the preferred transcript.

Gene-level features:

- pLI.
- Missense Z.
- LoF Z.
- LOEUF proxy: `oe_lof_upper`.
- Observed/expected missense and LoF metrics.

## Benign Proxy QC Draft

The presumed-benign class should not be generated until the final gene panel and
label ratios are chosen. A conservative first pass should:

- Use only reviewed disease-panel genes.
- Use GRCh38 coordinates.
- Exclude variants already labeled P/LP in ClinVar or expert-curated sources.
- Prefer PASS variants or variants without serious gnomAD filters.
- Focus on missense variants if matching the CardioBoost-style benign proxy.
- Require an AF/popmax threshold justified by disease prevalence, penetrance,
  inheritance pattern, and gene/disease mechanism.
- Track homozygote count separately rather than using it as a blind benign rule.
- Prevent leakage by separating proxy-benign variants from ClinVar B/LB labels.

## Bulk Download Decision

Do not download full gnomAD v4 VCF/Hail data until one of these is true:

- The local dbNSFP-derived gnomAD fields are missing a required feature.
- We need exhaustive benign sampling across large genomic intervals.
- We need fields not exposed by the API.

If escalation is needed, prefer targeted extraction over full local mirroring:

| Option | Use when |
|---|---|
| Browser variant Hail Table | Need browser-like exome/genome/joint frequencies and filters at scale |
| Per-chromosome VCF shards with tabix | Need direct coordinate lookup for defined intervals/variants |
| Cloud query/Hail job | Need large targeted extraction without downloading everything locally |
| Full local raw download | Last resort only, after storage and compute requirements are documented |
