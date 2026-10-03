# ClinVar review-quality and structure-feature audit

## ClinVar review quality

The primary training set contained 36,176 ClinVar rows: 23,541 with one review star, 12,436 with two stars, and 199 with three stars. One-star rows were retained at weight 0.70; two-star and three-star rows received weights 1.00 and 1.20. Review status was never a predictor.

A separate model trained only on ClinVar records with at least two stars used 12,635 training rows and 2,308 validation rows. On the identical 430 source-held-out variants, its AUROC was 0.984, AUPRC 0.989, and Brier score 0.055. The primary model values were 0.978, 0.985, and 0.050. The paired AUROC difference was 0.006 (95% CI -0.001 to 0.015).

## Protein structure

The upstream source pipeline mapped 3,284 of 3,591 rows to reviewed UniProt entries. AlphaFold residue pLDDT was recovered for 2,336 rows. DSSP and FreeSASA were calculated upstream for those rows. FoldX processed 1,230 high-confidence residue mappings, returned 1,169 values, and flagged 61 reference mismatches that were excluded without imputation.

After source rows were collapsed and joined to the final 85,677-row matrix, 524 rows had AlphaFold residue pLDDT and 313 had numeric FoldX values. Neither field was present in training or validation. DSSP and FreeSASA numeric outputs were absent from the final ready matrix. Every retained structure-related feature had zero frozen-model importance. These annotations did not provide learned evidence for the frozen 338-feature baseline. The separate registry-wide enhancement is evaluated in `results/model_performance/structure_enhancement/`.
