# Modeling Tables

This folder contains model-facing tables built from the cleaned coordinate
registry and validated feature layers.

The first skeleton table is intentionally variant-level: one row per
`chrom-pos-ref-alt` variant key. It preserves source-specific label provenance
so later train/test split logic can decide how to use ClinVar, HiRO, eMERGE,
and CardioBoost without losing where each label came from.

Current outputs:

- `interim/final_modeling_table_skeleton.tsv`
- `interim/final_modeling_table_skeleton.summary.json`
- `interim/dbnsfp_selected_features.tsv`
- `interim/alphamissense_direct_features.tsv`
- `interim/clingen_selected_features.tsv`
- `interim/final_modeling_table_with_clingen.tsv`
- `interim/clingen_selected_features.summary.json`
- `interim/protein_structure_selected_features.tsv`
- `interim/final_modeling_table_with_clingen_protein.tsv`
- `interim/protein_structure_selected_features.summary.json`

Important conventions:

- `genome_build` is `GRCh38`.
- `model_label_3class` is filled only when the deduplicated variant has exactly
  one label among `Benign`, `VUS`, and `Pathogenic`.
- Mixed-label variants are kept, but `clean_supervised_label` is `false` and
  `label_conflict_type` describes the conflict.
- Source-specific label columns such as `clinvar_labels`, `hiro_labels`,
  `emerge_labels`, and `cardioboost_labels` are provenance fields, not feature
  columns.
- Placeholder status columns are included for feature blocks that are not yet
  merged, such as direct gnomAD browser annotations, VEP/SpliceAI,
  DSSP/FreeSASA, and FoldX DDG.
- ClinGen is now merged in `final_modeling_table_with_clingen.tsv`. This is a
  non-destructive table; the original skeleton is unchanged.
- `clingen_selected_features.tsv` keeps compact model-ready ClinGen columns,
  binary helper flags, and trace columns for QC.
- Protein, UniProt, AlphaFold, DSSP, and FreeSASA features are now merged in
  `final_modeling_table_with_clingen_protein.tsv`. This is also
  non-destructive; `final_modeling_table_with_clingen.tsv` is unchanged.
- `protein_structure_selected_features.tsv` keeps one selected protein-context
  row per coordinate variant and records source/duplicate QC fields.
