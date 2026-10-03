# ClinVar

Primary raw source:

- `variant_summary.txt.gz`: ClinVar bulk tab-delimited variant summary.

Curated modeling slice:

- Script: `scripts/build_curated_clinvar_modeling_slice.py`
- Combined output: `interim/clinvar_cardiogenetics_modeling_slice.tsv`
- High-confidence output: `interim/clinvar_cardiogenetics_high_confidence.tsv`
- Lower-confidence output: `interim/clinvar_cardiogenetics_lower_confidence_star1.tsv`
- Summary: `interim/clinvar_cardiogenetics_modeling_slice.summary.json`

Current curation rules:

- GRCh38 only.
- Gene symbol intersects the audited cardiogenetics classifier keep panel:
  `datasets/gene_panels/cardiogenetics_classifier_genes.keep.txt`.
- Tier C likely-contamination genes are excluded from the source-union panel:
  `DMD`, `AKAP9`, `FPGT`, and `KNCH2`.
- Complete `chrom-pos-ref-alt` key from ClinVar VCF columns.
- `OriginSimple == germline`.
- Exclude rows with somatic/oncogenic ClinVar fields.
- Exclude cancer-context phenotype/name keywords.
- Exclude `Conflict_or_other` and `Other` labels from clean supervised labels.
- Keep only `Benign`, `Pathogenic`, and `VUS`.
- Split confidence tiers:
  - `high_confidence_stars_ge2`: ClinVar review stars 2 or 3.
  - `lower_confidence_star1`: ClinVar review star 1.

Current output counts:

| Slice | Rows | Benign | Pathogenic | VUS |
|---|---:|---:|---:|---:|
| Combined curated slice | 127,967 | 54,770 | 14,986 | 58,211 |
| High-confidence stars >=2 | 41,177 | 17,633 | 4,544 | 19,000 |
| Lower-confidence star 1 | 86,790 | 37,137 | 10,442 | 39,211 |

Notes:

- This is a curated ClinVar label pool, not a final train/test split.
- VUS rows are retained as a clinically meaningful third class but should be
  treated as noisier than B/P labels.
- The current modeling panel is the audited 59-gene keep list. Tier B genes
  remain conditional/review genes and can be excluded in a primary-only
  sensitivity analysis.
