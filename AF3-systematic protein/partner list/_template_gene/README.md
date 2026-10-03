# Gene Template For AF3 Systematic Review

Copy this folder structure for each new gene:

```text
gene_###_GENE/
  GENE_systematic_review.md
  01_gene_context.tsv
  02_sequence_audit.tsv
  03_candidate_partner_evidence.tsv
  04_af3_design_decisions.tsv
  05_variant_mapping_audit.tsv
  09_specific_citation_sources.tsv
  af3_job_json/
  af3_job_results/
```

Required outputs:

1. Systematic evidence/decision tables.
2. Specific citation/source table with Ctrl-F terms.
3. AF3 job JSONs or a clear not-generated reason.
4. Empty `af3_job_results/` folder for Dheyaa's server outputs.
