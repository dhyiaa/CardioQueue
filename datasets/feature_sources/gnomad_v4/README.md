# gnomAD v4

Role: feature injection and benign-proxy source.

Current status:

- GraphQL API tested with dataset enum `gnomad_r4`.
- Example variant annotations and gene constraints generated.
- Clean primary gene constraints generated and validated.
- A 200-variant clean-registry API pilot generated and validated, but the API is
  too slow for a monolithic full-registry run.
- Immediate local population-frequency features are now extracted from dbNSFP
  gnomAD fields in `interim/gnomad_from_dbnsfp.clean_combined.tsv`.
- Browser Hail Table escalation is prepared for the 44,548 clean-registry
  variants missing local dbNSFP-derived gnomAD fields, and the 100-variant pilot
  completed successfully.
- Full raw VCF/Hail mirroring remains deferred; the prepared Hail route reads
  the official browser table remotely after Google authentication and
  requester-pays access are configured.

Folders:

- `raw/`: source notes.
- `external/`: gene panels and variant query inputs.
- `interim/`: example API outputs and raw JSON.
- `scripts/`: gnomAD query and ClinVar-to-gnomAD extraction scripts.

Primary near-term strategy:

1. Use `interim/gnomad_from_dbnsfp.clean_combined.tsv` now for local AF,
   popmax-AF, and homozygote-count-like features.
2. Use `interim/gnomad_v4_gene_constraints.clean_primary.tsv` for gene-level
   constraint metrics.
3. Escalate to direct gnomAD Hail/VCF extraction only for fields missing from
   dbNSFP, such as detailed filters, transcript consequences, or exhaustive
   benign-proxy sampling.

Prepared Hail escalation:

- Missing variant input:
  `external/gnomad_missing_from_dbnsfp.clean_combined.tsv`
- Pilot input:
  `external/gnomad_missing_from_dbnsfp.pilot100.tsv`
- Script:
  `scripts/extract_gnomad_browser_hail_missing.py`
- Notes:
  `raw/browser_hail/README.md`

Installed local tools:

- Java: Temurin 11.0.31 at `~/opt/jdk-11.0.31+11/Contents/Home`.
- Hail: conda environment `hail-gnomad`, Python 3.11, Hail 0.2.138.
- Google Cloud CLI: `gcloud` 575.0.0 and `gsutil` 5.37.
- GCS connector: `raw/browser_hail/jars/gcs-connector-hadoop3-2.2.31-shaded.jar`.
- Local Spark service-account key:
  `~/.config/gcloud/gnomad-hail-local-key.json` (do not commit).

Pilot output:

- `interim/gnomad_browser_v4_1_missing_pilot100.tsv.bgz`
- 100 rows queried.
- 81 rows found in the browser Hail Table.
- 19 rows not found by exact GRCh38 locus+alleles.
