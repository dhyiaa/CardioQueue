# Gene Panel Audit

This note audits the current source-union gene panel before modeling. It does
not freeze the final model panel yet.

## Local Source Facts

- HiRO has 57 unique gene entries in the current model input table.
- eMERGE has 10 genes.
- CardioBoost has 16 genes.
- The current curated ClinVar slice has 58 genes and is ClinVar-only. Earlier
  first-pass counts included `AKAP9`, `DMD`, and `FPGT`; those are now excluded
  upstream by the current keep panel.
- CardioBoost local gene files do not include `TTN`, `DMD`, `AKAP9`, `ANK2`,
  `RYR2`, `DSP`, `DSG2`, `PKP2`, or many other expanded-panel genes.

Generated local audit tables:

- `datasets/gene_panels/interim/source_gene_counts.tsv`
- `datasets/gene_panels/interim/gene_panel_qc_flags.tsv`
- `datasets/gene_panels/interim/gene_panel_decisions.tsv`
- `datasets/gene_panels/interim/gene_panel_decisions.summary.json`

Generated panel files:

- `datasets/gene_panels/cardiogenetics_classifier_genes.keep.txt`
- `datasets/gene_panels/cardiogenetics_classifier_genes.remove.txt`
- `datasets/gene_panels/cardiogenetics_classifier_genes.sensitivity_or_separate.txt`

## Claims Checked

| Claim | Assessment | Evidence/Action |
|---|---|---|
| The current panel is partially contaminated | Correct | Current source-union panel includes readthrough/withdrawn/typo symbols and genes whose ClinVar rows may be non-cardiac-context labels. |
| ClinVar curated counts are only ClinVar | Correct | Built only from `datasets/clinvar/variant_summary.txt.gz`; no HiRO/eMERGE/CardioBoost/gnomAD/dbNSFP rows are included. |
| `FPGT` is likely an artifact | Correct | `FPGT` is a real protein-coding gene but not a cardiogenetics gene in this context; it enters through the nearby/readthrough `FPGT-TNNI3K`/`TNNI3K` locus. Remove or quarantine from primary modeling. |
| `FPGT-TNNI3K` should simply be removed | Too blunt | It is an HGNC readthrough locus. Several local HiRO records annotated as `FPGT-TNNI3K` match `TNNI3K` ClinVar/transcript context. Map transcript-aware to `TNNI3K` when supported; otherwise quarantine. |
| `KNCH2` is suspicious | Correct | It is not a valid current panel symbol and the HiRO HGVS `c.2863C>G p.Leu955Val` matches `KCNH2` in ClinVar. Correct to `KCNH2` only when variant/transcript evidence supports it. |
| `KCNE1B` is suspicious | Correct | HGNC reports `KCNE1B` as withdrawn because it is in a false-duplication region. Collapse evidence-supported rows to `KCNE1`; otherwise quarantine. |
| `AKAP9` should be excluded from the primary model | Supported | ClinGen classifies `AKAP9`-LQTS as disputed; curated ClinVar has only 3 pathogenic rows for 3,341 total rows. |
| `ANK2` should be treated cautiously | Supported | ClinGen classifies ANK2 relationships with LQTS, Brugada, and CPVT as disputed; large VUS burden makes it risky for primary training. |
| `TTN` should be handled separately | Supported | ClinGen supports TTN for DCM, especially truncating variants, but HCM evidence is limited and TTN missense is high-noise. CardioBoost local gene files exclude `TTN`. |
| `DMD` should not be in the primary cardiac-specific model | Supported | ClinGen curation in local output is definitive for progressive muscular dystrophy, not a primary cardiogenetics label context; cardiac disease is phenotype-specific/secondary. |
| `TRPM4` is a straightforward keeper | Not safe | ClinGen classifies TRPM4-Brugada as refuted. Keep only if the target explicitly includes conduction disease with a documented evidence rule, otherwise sensitivity-only. |
| `TNNI3K` is real but sparse | Correct | ClinGen local output shows moderate evidence for DCM. Use in expanded/DCM-pooled or sensitivity analysis rather than a per-gene model. |

## Recommended Panel Policy

Use at least three panels rather than one:

1. `primary_cardioBoost_aligned_panel`: CardioBoost genes plus any explicitly
   justified core HiRO genes after ClinGen/context review.
2. `expanded_cardiogenetics_panel`: broader HCM/DCM/ARVC/channelopathy,
   phenocopy, and RASopathy genes, with source flags and gene-group pooling.
3. `sensitivity_or_quarantine_panel`: `TTN`, `DMD`, `AKAP9`, `ANK2`, `TRPM4`,
   `FPGT`, `FPGT-TNNI3K`, `KCNE1B`, `KNCH2`, and any unresolved aliases.

Primary modeling should not use the raw source-union panel until these symbols
are normalized and the final gene panel is frozen.

## Current Applied Decisions

These decisions are applied by
`datasets/gene_panels/scripts/apply_gene_panel_decisions.py`. The script does
not delete any raw rows; it creates explicit panel files and ClinVar-derived
modeling slices.

| Decision bucket | Genes | Primary model action |
|---|---|---|
| Remove contamination/quarantine | `AKAP9`, `DMD`, `FPGT`, `FPGT-TNNI3K`, `KCNE1B`, `KNCH2` | Exclude from primary panel. Current ClinVar curated slice already has zero rows for these genes; source tables may still contain alias/private rows needing normalization. |
| Separate or sensitivity-only | `ANK2`, `SOS1`, `TNNI3K`, `TRPM4`, `TTN` | Exclude from primary ClinVar training slice; keep in documented holdout/sensitivity files. |
| Conditional phenocopy/expanded scope | `EMD`, `GLA`, `HRAS`, `KRAS`, `LAMP2`, `MAP2K2`, `NRAS`, `PTPN11`, `RAF1`, `RIT1`, `SLC22A5` | Retain for now but flagged as conditional phenocopy/syndromic-overlap genes. |
| Primary keep with monitoring | `CACNB2`, `FHOD3`, `MIB1`, `SCN1B`, `SNTA1` | Retain but monitor evidence strength and feature contribution. |

Current ClinVar post-QC slices:

| Slice | Rows | Benign | Pathogenic | VUS | Notes |
|---|---:|---:|---:|---:|---|
| Primary after gene QC | 83,915 | 33,835 | 9,052 | 41,028 | Excludes TTN and Tier B sensitivity genes. |
| Holdout/sensitivity | 44,052 | 20,935 | 5,934 | 17,183 | Contains `ANK2`, `SOS1`, `TNNI3K`, `TRPM4`, and `TTN`. |
| TTN separate review | 34,804 | 17,203 | 5,765 | 11,836 | Includes a coarse `ttn_likely_truncating` flag; 2,262 rows currently flag as likely truncating/splice-like. |
| Removed contamination | 0 | 0 | 0 | 0 | Zero because the current curated ClinVar slice already excludes `AKAP9`, `DMD`, and `FPGT`; the remove file still protects future rebuilds. |
