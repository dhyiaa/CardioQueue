# AlphaFold 3 Feature Source

Role: AF3 interaction-aware complex mapping for a new CatBoost model version.

This folder is separate from the existing AlphaFold/FoldX protein-structure
feature branch. The current primary CatBoost model is preserved. AF3 outputs
will be used to build a new feature table and a new AF3-augmented model version.

Planned folders:

- `raw/`: downloaded AlphaFold Server job outputs, one folder per job.
- `interim/`: gene, sequence, partner, and completed-job manifests.
- `jobs/`: AF3 job manifests and manual submission packets.
- `processed/`: parsed AF3 residue-level and variant-level feature tables.
- `scripts/`: reproducible build scripts for audits and AF3 feature parsing.

First generated audit:

```text
results/af3/gene_panel_audit.tsv
results/af3/gene_panel_audit.summary.json
```

Current generated manifests:

```text
datasets/feature_sources/alphafold3/interim/uniprot_gene_to_sequence.tsv
datasets/feature_sources/alphafold3/interim/af3_all_entities.fasta
datasets/feature_sources/alphafold3/interim/af3_manual_cardiac_complex_whitelist.tsv
results/af3/af3_full_panel_partner_manifest.tsv
results/af3/af3_sequence_partner_manifest.summary.json
datasets/feature_sources/alphafold3/interim/evidence/uniprot_evidence_by_gene.tsv
results/af3/af3_full_panel_partner_manifest.uniprot_evidence.tsv
results/af3/af3_complex_evidence_readiness.tsv
results/af3/af3_partner_evidence_review_queue.tsv
datasets/feature_sources/alphafold3/interim/evidence/external_evidence_by_gene.tsv
results/af3/af3_full_panel_partner_manifest.external_evidence.tsv
results/af3/af3_complex_combined_evidence_readiness.tsv
results/af3/af3_partner_combined_evidence_review_queue.tsv
results/af3/af3_none_found_candidate_partners.tsv
results/af3/af3_complex_design_review.tsv
results/af3/af3_complex_groups.tsv
datasets/feature_sources/alphafold3/jobs/af3_job_manifest.tsv
datasets/feature_sources/alphafold3/jobs/manual_submission_packets/
```

Important status:

- These are review manifests, not final AlphaFold 3 server jobs.
- The sequence source is the local reviewed human UniProt FASTA, not AlphaFold 2
  structures.
- Partner rows remain `review_pending` until evidence columns are filled from
  UniProt, Reactome, IntAct, RCSB/PDB, and literature review.
- The current script does not modify existing modeling matrices or previous
  protein-structure feature tables.
- UniProt evidence has been cached locally under `raw/uniprot_evidence/`.
  This is an evidence review layer, not final AF3 submission approval.
- Reactome, RCSB/PDB, and IntAct evidence has been cached locally under
  `raw/reactome_evidence/`, `raw/rcsb_evidence/`, and `raw/intact_evidence/`.
- The combined-readiness table is the current safest file for deciding which
  complexes can move toward AF3 design review.
- Four draft manual submission packets have been created, but they remain
  `draft_needs_user_review` and `not_submitted`.

Pilot-v1 pathway:

```text
results/af3/pilot_v1/af3_pilot_v1_gene_set.tsv
results/af3/pilot_v1/af3_pilot_v1_complex_plan.tsv
results/af3/pilot_v1/af3_pilot_v1_pathway.summary.json
datasets/feature_sources/alphafold3/jobs/pilot_v1/af3_pilot_v1_submission_manifest.tsv
datasets/feature_sources/alphafold3/jobs/pilot_v1/af3_pilot_v1_completed_job_manifest_template.tsv
datasets/feature_sources/alphafold3/jobs/pilot_v1/manual_submission_packets/
datasets/feature_sources/alphafold3/raw/pilot_v1_server_outputs/
datasets/feature_sources/alphafold3/processed/pilot_v1/af3_pilot_v1_feature_schema.tsv
datasets/feature_sources/alphafold3/processed/pilot_v1/af3_pilot_v1_variant_feature_scaffold.tsv
datasets/feature_sources/alphafold3/interim/evidence/hpa_evidence_by_gene.tsv
results/af3/af3_full_panel_partner_manifest.final_evidence_stack.tsv
results/af3/pilot_v1/af3_pilot_v1_testable_genes_and_interactions.tsv
datasets/feature_sources/alphafold3/jobs/pilot_v1/server_json/
```

Pilot-v1 status:

- 23 target genes are in the first AF3 pilot path.
- 14 complex/job rows are planned: 4 existing first-wave draft packets plus
  10 second-wave review packets.
- 22,119 of 85,677 current modeling rows are in pilot-v1 genes; 22,116 already
  have protein positions available for residue mapping once AF3 outputs exist.
- The second-wave packets are review packets. They are not approved submission
  packets until the evidence checklists are manually accepted.
- Returned AlphaFold Server results should be placed under
  `raw/pilot_v1_server_outputs/` and recorded in
  `jobs/pilot_v1/af3_pilot_v1_completed_job_manifest_template.tsv`.
- HPA localization/expression columns have been added as context evidence.
  HPA is not counted as direct interaction evidence.
- AlphaFold Server JSON draft inputs have been generated in
  `jobs/pilot_v1/server_json/`. They use the `alphafoldserver` dialect and
  should be reviewed before upload.
