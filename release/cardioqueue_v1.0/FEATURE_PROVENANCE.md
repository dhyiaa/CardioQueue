# CardioQueue v1.0 feature provenance

CardioQueue scores an annotated variant row. Its 300 predictors are retrieved from databases or calculated by reproducible code. The scorer does not ask a user to assign ACMG/AMP criteria, pathogenicity labels, segregation evidence, phenotype matches, or functional-study interpretations.

## Automatically derived groups

| Group | Source or calculation | Typical requirement |
|---|---|---|
| Variant consequence and transcript fields | Ensembl VEP | VEP cache and GRCh38 allele |
| Splicing | SpliceAI fields emitted through VEP | SpliceAI VEP plugin and files |
| Population frequency | gnomAD 4.1 browser exports and dbNSFP | Local source files or API-derived cache |
| Effect scores and conservation | dbNSFP 5.3.1a, including REVEL, CADD, SIFT, PolyPhen-2, MetaLR, ESM1b, GERP++, phyloP, and phastCons | Licensed or permitted local dbNSFP installation |
| Missense effect | AlphaMissense direct GRCh38 table and dbNSFP fallback | Local AlphaMissense table |
| Gene evidence | ClinGen-derived gene validity fields | Versioned ClinGen export |
| Protein context | Reviewed UniProt sequence and feature records | Versioned UniProt files |
| Residue geometry | AlphaFold Database v6, DSSP, and FreeSASA | Local structures and executables |
| Stability change | FoldX 5.1 | Separate FoldX license and executable |
| Missingness and summary fields | Deterministic transformations of the groups above | Repository scripts |

AlphaFold 3 is not part of CardioQueue v1.0.

## Current interface boundary

`score_variants.py` validates and scores a table that already contains every field in `model/feature_columns.json`. This release does not bundle third-party databases or the licensed FoldX executable. The manuscript repository contains the source-specific annotation and join scripts used to build the study matrix. A public workflow should pin all resource versions, normalize GRCh38 alleles, verify reference residues, and emit the frozen schema before scoring.

Missing values are accepted, but missingness patterns can differ across laboratories and can alter ranking. A deployment site should run a local retrospective validation after recreating the annotation stack.
