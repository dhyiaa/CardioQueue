# ClinGen

Role: curated-evidence feature source.

Current status:

- Gene-disease validity API workflow validated and run for the master registry gene panel.
- Variant-level Evidence Repository gene-panel workflow validated and run for the master registry gene panel.
- Dosage sensitivity downloads validated and parsed into the final ClinGen feature table.
- Final joinable ClinGen modeling feature table generated for the current modeling skeleton.

Folders:

- `raw/`: raw API response examples.
- `external/`: source notes.
- `interim/`: generated example tables.
- `scripts/`: ClinGen API helper.

Current example outputs:

- `interim/gene_validity.example.tsv`
- `interim/gene_validity.master_registry_panel.tsv`
- `interim/variant_pathogenicity.example.tsv`
- `interim/variant_pathogenicity.master_registry_panel.tsv`
- `interim/clingen_modeling_features.tsv`
- `interim/clingen_modeling_features.summary.json`
- `raw/dosage_sensitivity/gene_dosage.csv`
- `raw/dosage_sensitivity/gene_dosage_GRCh38.tsv`

Current modeling-feature coverage:

- `interim/clingen_modeling_features.tsv`: 86,886 rows, one per modeling skeleton variant.
- Gene validity flags: 84,983 rows with a ClinGen gene-validity match; 1,903 rows not found.
- Dosage sensitivity flags: 73,792 rows with dosage data; 13,094 rows not found.
- Variant-level VCEP/ERepo evidence: 319 exact GRCh38 genomic-HGVS variant matches.

Variant-level ClinGen matching is intentionally conservative: ERepo HGVS strings
are converted to GRCh38 `chrom-pos-ref-alt` keys when possible, and only exact
matches are marked `ok`. Gene-only matches are not promoted to variant evidence.
