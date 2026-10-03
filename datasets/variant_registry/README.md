# Master Variant Registry

This folder contains the project-level variant registry.

The registry is not a modeling table by itself. It is the deduplicated variant
spine that later feature tables join onto.

Current outputs:

- `interim/master_variant_registry.tsv`: one row per GRCh38
  `chrom-pos-ref-alt` key where coordinates are available.
- `interim/master_variant_registry_sources.tsv`: one row per source membership,
  preserving source-specific row IDs, labels, and provenance.
- `interim/master_variant_registry_source_inventory.tsv`: one row per input
  source, preserving source-row counts separately from deduplicated coordinate
  counts.
- `interim/master_variant_registry.summary.json`: build summary and source
  counts.

Important conventions:

- Registry key: `chrom-pos-ref-alt`
- Preferred assembly: GRCh38
- Source flags are additive; a variant can appear in more than one source.
- HiRO rows carry `phenotype_available=true` because HiRO is the only current
  source with linked patient-level phenotype fields.
- HiRO source counts must not be read from the deduplicated coordinate registry.
  The current HiRO table has 482 private/patient-linked modeling rows. Of those,
  361 rows currently have complete coordinate keys and collapse to 224 unique
  `chrom-pos-ref-alt` variants after deduplication. The remaining 121 HiRO rows
  need HGVS/protein/transcript fallback mapping before they can enter the
  coordinate registry.
- If a HiRO variant overlaps ClinVar by coordinate, provenance still remains
  source-specific. The HiRO row keeps its private label/patient context and
  should not be described as a ClinVar observation.
- Protein-only CardioBoost V2 / DYNA rows are not included in this coordinate
  registry because they do not provide `chrom-pos-ref-alt` join keys.
