# Feature Fill Status

This note tracks how feature gaps should be filled after the clean combined
registry is built.

## Current Clean Registry

The clean combined registry is:

- `datasets/variant_registry/interim/clean_combined_variant_registry.tsv`
- `datasets/variant_registry/interim/clean_combined_variant_registry_sources.tsv`
- `datasets/variant_registry/interim/clean_combined_variant_registry.summary.json`

This registry uses:

1. `datasets/clinvar/interim/clinvar_primary_after_gene_qc.tsv`
2. `datasets/hiro/full_dataset/model_inputs/ml_baseline_features.csv`
3. `datasets/emerge/full_arrhythmia_gene_dataset/data/emerge_arrhythmia_ml_baseline_features.csv`
4. `datasets/cardioboost/public_dataset/processed_for_cardiogenetics/model_inputs/reannotated_ml_baseline_features.csv`

It intentionally does not use raw, broad ClinVar directly. HiRO/eMERGE/
CardioBoost rows are preserved when coordinate-complete, and source rows without
coordinates remain counted for later HGVS/protein fallback mapping.

Current registry counts:

| Item | Count |
|---|---:|
| Deduplicated coordinate rows | 86,886 |
| Source-membership rows | 87,385 |
| ClinVar primary rows | 83,915 |
| HiRO rows with coordinates | 361 |
| HiRO rows without coordinates | 121 |
| eMERGE rows with coordinates | 2,754 |
| CardioBoost rows with coordinates | 355 |

## Local-Only Feature Registry

The first local feature pass is:

- `datasets/variant_registry/interim/clean_combined_local_feature_registry.tsv`
- `datasets/variant_registry/interim/clean_combined_local_feature_registry.summary.json`

This table joins local features only:

1. dbNSFP v5.3.1 features from
   `datasets/feature_sources/insilico/interim/dbnsfp_features_master_variant_registry.tsv`
2. Protein/UniProt/AlphaFold features from
   `datasets/feature_sources/protein_structure/interim/protein_features.tsv`

Current local feature coverage:

| Feature block | Status | Count |
|---|---|---:|
| dbNSFP exact match | `ok` | 42,338 |
| dbNSFP no exact match | `not_found` | 44,548 |
| Protein feature row ok | `ok` | 3,091 |
| Protein not mapped | `not_mapped` | 236 |
| Protein no UniProt mapping | `no_uniprot_mapping` | 2 |
| No local protein row | `not_available` | 83,557 |
| AlphaFold residue pLDDT ok | `ok` | 2,179 |
| AlphaFold fragment unresolved | `fragment_mapping_unresolved` | 912 |
| AlphaFold not queryable | `not_queryable` | 238 |

## Local Sources We Can Use Now

These can be run without waiting for external API throughput.

| Source/tool | Feature types | Current action |
|---|---|---|
| dbNSFP | REVEL, SIFT, PolyPhen, MetaLR, FATHMM-XF, CADD-like scores, AlphaMissense, ESM1b, conservation, gnomAD 4.1 fields where present | Already joined into the local feature registry; gnomAD proxy table extracted for immediate use |
| AlphaMissense local hg38 file | AlphaMissense score/class for missense variants | Local file exists; can run a full direct join if we want AlphaMissense independent of dbNSFP |
| UniProt features | Domains, regions, motifs, active/binding sites, normalized protein position | Already in `protein_features.tsv`; next step is joining and QC in final modeling table |
| AlphaFold bulk tar | Residue pLDDT and model-fragment metadata | Partially joined; next step is resolving fragmented proteins such as `TTN`, `DSP`, `RYR2` |
| DSSP/mkdssp | Secondary structure and accessibility from structures | Installed and smoke-tested; batch extraction still needed |
| FreeSASA | Solvent accessibility and buried/exposed flags | Installed and smoke-tested; batch extraction still needed |
| FoldXPro | Delta-delta-G for eligible missense variants | Installed and smoke-tested; pilot batch still needed |
| Existing source-model features | Current baseline/enriched columns in HiRO/eMERGE/CardioBoost files | Available, but should be harmonized and not confused with newly recomputed features |

## API Or External-Throughput Sources

These can run in the background, but should be staged and cached.

| Source/tool | Feature types | Current action |
|---|---|---|
| gnomAD v4 GraphQL API | AF, popmax AF, homozygote/hemizygote count, filters, consequences, gene constraints | Doctor passed; clean 53-gene constraint run and 200-variant pilot both validated; full one-by-one registry API run deferred |
| gnomAD v4 browser Hail Table | Browser-like variant frequency, filter, and annotation records for variants missing local dbNSFP gnomAD fields | 100-variant pilot completed successfully against `gs://gcp-public-data--gnomad/release/4.1/ht/browser/gnomad.browser.v4.1.sites.ht`: 81 ok, 19 not found |
| CADD API | CADD raw/PHRED for targeted variants | API script exists; use for small missing sets or decide a remote-tabix/bulk strategy |
| ClinGen APIs/downloads | Gene-disease validity, dosage, variant evidence | Gene/dosage workflows exist; variant-level full matching still pending |
| VEP REST or local VEP | Consequence/transcript/splice harmonization | Deferred unless dbNSFP/gnomAD consequences are insufficient |

## Deferred Or Decision-Dependent Sources

| Source/tool | Why deferred |
|---|---|
| Full gnomAD bulk VCF/Hail downloads | Large storage/compute; defer unless dbNSFP-derived gnomAD fields are insufficient |
| Local VEP cache and plugins | Large setup; only needed if consequence/splice harmonization requires it |
| SpliceAI plugin/precomputed source | Needs VEP/plugin decision or verified precomputed source |
| ESM1b / ESM-variant external source | dbNSFP already has ESM1b fields for many rows; independent source needs mapping QC |
| Rosetta `cartesian_ddg` | Heavier compute; best as sensitivity subset after FoldXPro pilot |
| UCSC Multiz orthologous-alignment features | Harder to reproduce; stretch goal after first model |

## Immediate Order Of Work

1. Use `clean_combined_local_feature_registry.tsv` as the first local feature
   spine.
2. Use `gnomad_from_dbnsfp.clean_combined.tsv` as the immediate local gnomAD
   population-frequency layer.
3. Merge validated gnomAD gene constraints into the local feature registry.
4. Run the prepared gnomAD browser Hail Table pilot once Google Cloud
   authentication and requester-pays billing access are configured, then scale
   to the full 44,548-row dbNSFP-missing set if the pilot validates cleanly.
5. Add DSSP/FreeSASA batch columns to the protein feature table.
6. Run a small FoldXPro missense pilot before any large delta-delta-G batch.

## gnomAD API Pilot Results

Generated inputs:

- `datasets/feature_sources/gnomad_v4/external/clean_primary_gene_panel.txt`
- `datasets/feature_sources/gnomad_v4/external/clean_combined_registry_variants.tsv`
- `datasets/feature_sources/gnomad_v4/external/clean_combined_registry_variants.pilot200.tsv`

Generated outputs:

- `datasets/feature_sources/gnomad_v4/interim/gnomad_v4_gene_constraints.clean_primary.tsv`
- `datasets/feature_sources/gnomad_v4/interim/gnomad_v4_variant_annotations.clean_registry_pilot200.tsv`

Validation results:

| gnomAD run | Rows | `ok` | `not_found` | Notes |
|---|---:|---:|---:|---|
| Clean primary gene constraints | 53 | 53 | 0 | All queried primary-panel genes returned constraint records. |
| Clean registry variant pilot | 200 | 144 | 56 | Pilot validates, but request speed is slow enough that full 86,886-row API annotation should be chunked or replaced with bulk/remote extraction. |

For the 200-variant pilot, 144 rows had joint gnomAD records, 115 had exome
records, and 94 had genome records.

## gnomAD From dbNSFP Local Layer

Immediate population-frequency features have been extracted from local dbNSFP
fields:

- `datasets/feature_sources/gnomad_v4/interim/gnomad_from_dbnsfp.clean_combined.tsv`
- `datasets/feature_sources/gnomad_v4/interim/gnomad_from_dbnsfp.clean_combined.summary.json`

| Metric | Count |
|---|---:|
| Clean-registry rows | 86,886 |
| Rows with local gnomAD proxy fields | 42,338 |
| Rows without local dbNSFP gnomAD fields | 44,548 |
| Nonmissing joint AF | 42,338 |
| Nonmissing joint popmax AF | 42,338 |
| Nonmissing joint homozygote alternate count | 42,338 |

This is now the preferred near-term gnomAD feature source. Direct gnomAD
browser/Hail/VCF extraction is reserved for missing fields, detailed filters,
transcript consequences, or exhaustive benign-proxy sampling.

## gnomAD Browser Hail Table Escalation

The direct browser Hail Table route has been prepared for the remaining variants
that were not found through local dbNSFP-derived gnomAD fields.

Prepared inputs:

- `datasets/feature_sources/gnomad_v4/external/gnomad_missing_from_dbnsfp.clean_combined.tsv`
- `datasets/feature_sources/gnomad_v4/external/gnomad_missing_from_dbnsfp.pilot100.tsv`

Prepared script:

- `datasets/feature_sources/gnomad_v4/scripts/extract_gnomad_browser_hail_missing.py`

Prepared documentation:

- `datasets/feature_sources/gnomad_v4/raw/browser_hail/README.md`
- `docs/gnomad_v4_plan.md`

Current missing-set summary:

| Group | Count |
|---|---:|
| Variants missing local dbNSFP-derived gnomAD fields | 44,548 |
| ClinVar-only | 41,771 |
| eMERGE-only | 2,740 |
| CardioBoost-only | 18 |
| HiRO-only | 10 |
| ClinVar+HiRO | 6 |
| CardioBoost+ClinVar | 3 |

The local code and core tools are ready:

- Java: Temurin 11.0.31 installed at `~/opt/jdk-11.0.31+11/Contents/Home`.
- Python/Hail: `hail-gnomad` conda environment with Python 3.11 and Hail
  0.2.138.
- Google Cloud CLI: `gcloud` 575.0.0 and `gsutil` 5.37.
- Google Cloud Storage connector:
  `datasets/feature_sources/gnomad_v4/raw/browser_hail/jars/gcs-connector-hadoop3-2.2.31-shaded.jar`.
- Local Spark service-account key:
  `~/.config/gcloud/gnomad-hail-local-key.json` (do not commit).

Pilot result:

| Output | Rows | `ok` | `not_found` |
|---|---:|---:|---:|
| `datasets/feature_sources/gnomad_v4/interim/gnomad_browser_v4_1_missing_pilot100.tsv.bgz` | 100 | 81 | 19 |

The next step is to parse the JSON payload into compact feature columns before
scaling from the 100-row pilot to the full 44,548-row missing set.
