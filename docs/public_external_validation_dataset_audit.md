# Public External Validation Dataset Audit

Date: 2026-09-25

## Decision

A stronger public external evaluation is feasible without using 570 sampled
ClinVar rows. The best design is to add two separately reported clinical
curation sources, not to combine every public variant into an unlabeled pool.

1. **SHaRe HCM 2026Q1** is the highest-value source. It contains 822 current
   binary classifications (561 P/LP and 261 B/LB) with GRCh38 coordinate IDs.
2. **Fernandez-Falgueras et al. PLOS ONE** contributes 197 binary records in
   the current CardioQueue panel (103 P/LP and 94 B/LB) after transcript-aware
   mapping.
3. The existing 430 HiRO/eMERGE/CardioBoost rows remain a separate primary
   cohort. Cohort-specific results must precede any pooled result.

After preserving precedence for 47 SHaRe variants already represented in the
existing external strata, SHaRe could add approximately 775 variants. The PLOS
source shares 43 gene/cDNA variants with SHaRe; all 43 have concordant binary
labels. Before coordinate normalization and complete cross-source
deduplication, these sources imply more than 1,300 total held-out observations
when combined with the existing 430. The final count will be lower.

## Tier 1: Suitable after required gates

### SHaRe Genomic Data Browser 2026Q1

- Public repository: https://github.com/ImperialCardioGenetics/share_browser
- Methods/preprint: https://doi.org/10.64898/2026.08.05.26359735
- Data: 2,637 GRCh38 variant rows; 822 final B/LB or P/LP rows; 26 genes.
- Binary composition: 261 B/LB and 561 P/LP.
- Curation: baseline clinical-laboratory classifications plus SHaRe
  adjudication; 899 rows in the full file are marked manually curated. The
  methods describe seven cardiovascular genetic counselors/variant curators
  and use of gene-specific ACMG/AMP specifications where available.
- Local overlap: 505/822 coordinates are already represented in the modeling
  matrix; 419 are currently train/validation rows and must be quarantined
  before a complete retraining. Forty-seven overlap existing held-out source
  strata and should retain existing source precedence.
- Strength: genuine multicenter HCM registry classifications, both labels,
  exact GRCh38 alleles, carrier/site counts, and current curation.
- Limit: predominantly HCM and sarcomeric genes; CardioBoost also used older
  SHaRe data, so SHaRe and CardioBoost are not independent evidence ecosystems.
- **Gate:** the repository states "All rights reserved" and contains no reuse
  license. Obtain written permission for publication analysis and release of a
  derived manifest before using this as a headline cohort.

### PLOS inherited-cardiovascular reinterpretation cohort

- Article: https://doi.org/10.1371/journal.pone.0297914
- Supplement: CC BY spreadsheet preserved under
  `datasets/external_validation_candidates/plos_cardiovascular_reinterpretation/`.
- Published scope: 1,425 observations from 618 index cases; the article reports
  1,131 unique variants across cardiomyopathies, channelopathies, and sudden
  cardiac death.
- Supplement: 1,261 listed rows with original class, updated class, and ACMG
  criteria. Updated labels include 176 B, 128 LB, 72 LP, 58 P, and 826 VUS.
- Current-panel binary subset: 197 unique records across 38 genes, comprising
  94 B/LB and 103 P/LP.
- Strength: independent clinical reinterpretation, broad cardiogenetic scope,
  both labels, explicit ACMG criteria, and clear CC BY reuse terms.
- Limit: gene/cDNA and protein HGVS are supplied without transcript accessions
  or genomic coordinates. Every allele requires transcript-aware normalization,
  liftover where necessary, reference validation, and manual ambiguity review.

## Tier 2: Valuable, but not independent external validation

### ClinGen cardiac VCEPs

The live ClinGen Evidence Repository contains 149 binary cardiac-VCEP
assertions: 129 from the Cardiomyopathy VCEP and 20 from the Potassium Channel
Arrhythmia VCEP. RYR2, sodium/calcium-channel, desmosomal, and congenital-heart
panels have not yet released variant assertions. These labels are authoritative
and condition-specific, but ClinGen curations are submitted to ClinVar and use
the same public evidence ecosystem. Use them as a separately named
**expert-curated challenge cohort**, not external validation.

### CardioClassifier curated variants

The public repository contains 84 variant-evidence rows and the publication
benchmarked 57 expertly curated variants. This is useful for ACMG-rule
concordance, but it is small, partly linked to ClinVar/ACGV evidence, lacks a
clean downloadable binary truth table, and overlaps the Imperial/SHaRe
ecosystem. It should not inflate the primary held-out count.

### BBI Clinical Variant Database

The BBI-CVD contains 4,770 unique clinical variants across 1,244 genes, but the
full catalog is not publicly downloadable and is not cardiology-specific.
Deidentified access can be requested from `BBICVD@UW.edu`. This is a promising
future independent source if the investigators can provide cardiac-gene rows,
five-tier laboratory classifications, coordinates, dates, and permission for
publication.

## Tier 3: One-class sensitivity cohorts only

Several public studies provide credible P/LP variants but no B/LB variants from
the same ascertainment process. They may test P/LP sensitivity but must not be
combined with benign variants from another source to estimate AUROC,
specificity, PPV, or calibration.

| Source | Public usable signal | Appropriate role |
|---|---:|---|
| Pediatric Brugada registry | 28 SCN5A P/LP variants | SCN5A sensitivity stress test |
| Early-onset atrial fibrillation cohort | 131 P/LP-positive participants; variant supplements | P/LP recovery only |
| DCM Precision Medicine Study | ACMG-classified variants from 97 test probands | DCM P/LP/VUS descriptive stress test |
| Arrhythmogenic cardiomyopathy cohort | 15 P/LP plus 66 VUS in core ACM genes | ACM sensitivity/VUS ranking |
| Pediatric cardiomyopathy cohorts | Variant supplements, mostly P/LP/VUS | Pediatric sensitivity only |

## Rejected as primary truth sets

- **ClinVar, ClinVar Miner, and CardioVar:** same label ecosystem as development.
- **Cardiac VariantFX:** aggregate burden resource, not a five-tier truth set;
  terms ask users not to publish global analyses of unpublished cohorts without
  coordination.
- **Atlas of Cardiac Genetic Variation:** useful case/control and laboratory
  evidence, but no clean licensed bulk five-tier benchmark; substantial overlap
  with LMM/OMGL and CardioBoost-era resources.
- **LOVD and locus-specific databases:** heterogeneous submitter assertions and
  inconsistent bulk access; not one prespecified ACMG-adjudicated cohort.
- **MaveDB/functional maps:** independent functional evidence, not clinical
  ACMG classifications. Suitable for mechanism-specific validation only.
- **Population cohorts such as UK Biobank/All of Us:** public papers usually
  release selected P/LP findings or aggregate counts, not complete same-source
  B/LB and P/LP truth tables.

## Required evaluation protocol

1. Freeze source releases, files, checksums, and a harmonization protocol before
   inspecting CardioQueue predictions.
2. Normalize every allele independently to GRCh38 and require reference-allele,
   transcript, cDNA, and protein concordance where available.
3. Deduplicate in this precedence order: HiRO, eMERGE, CardioBoost, SHaRe, PLOS.
   Preserve all source memberships for sensitivity analyses.
4. Remove every evaluation allele and normalized equivalent from training,
   validation, feature selection, calibration, and threshold selection.
5. Retrain the frozen model recipe once. Do not tune to any new source.
6. Report each source separately, then disease group, gene, consequence, and
   ClinVar-presence strata. A pooled estimate is secondary.
7. Report variant-level source-stratified bootstrap and gene-clustered
   bootstrap intervals. Add leave-one-source-out pooled estimates.
8. Evaluate the no-precomputed-effect-predictor model on the same rows. Public
   clinical curations often use CADD, REVEL, population frequency, ClinVar, and
   similar evidence also available to CardioQueue.
9. For SHaRe, report a ClinVar-absent subset. The current file contains 115
   binary rows with no mapped ClinVar class (95 P/LP and 20 B/LB), providing a
   harder test of dependence on ClinVar-visible variants.
10. Never call the pooled set "independent validation." Preferred wording is
    "retrospective multi-source, exact-variant-held-out evaluation."

## Reviewer-resistant claims

The enlarged evaluation can support transfer across institutions, disease
groups, genes, and variant consequences. It cannot establish prospective
clinical utility, independence from public evidence, VUS reclassification,
patient-level diagnosis, or generalization to all cardiogenetic variants.
