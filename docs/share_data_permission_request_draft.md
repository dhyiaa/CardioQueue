# Draft: SHaRe Data Permission Request

**Subject:** Request to use SHaRe 2026Q1 variant classifications for external evaluation of CardioQueue

Dear SHaRe Genomic Data Browser team,

We are developing CardioQueue, a retrospective machine-learning system for
prioritizing cardiogenetic variants for expert evidence review. The model is
trained on ClinVar-derived labels and evaluated on separately held-out clinical
and public sources. We would like to use the public
`HCM_Gen_Variant_Annotation_2026Q1.csv` release from the SHaRe Genomic Data
Browser as a separately reported HCM evaluation cohort.

Our proposed analysis would:

- use only deidentified variant-level fields already present in the public file;
- restrict evaluation to final B/LB and P/LP classifications;
- remove every SHaRe allele and normalized equivalent from model training,
  validation, feature selection, calibration, and threshold selection before
  retraining;
- deduplicate against HiRO, eMERGE, CardioBoost, and other evaluation sources;
- report SHaRe-specific performance separately from other cohorts;
- report gene-, consequence-, and ClinVar-presence sensitivity analyses;
- acknowledge SHaRe, the contributing investigators and centers, and the
  Genomic Data Browser as requested; and
- release analysis code and, if permitted, a derived variant-level evaluation
  manifest containing coordinates, gene, classification, and source provenance.

We noted that the repository states "All rights reserved" and does not include
an explicit reuse license. Could you please confirm whether the proposed
publication analysis is permitted and whether we may redistribute the derived
evaluation manifest? We would also welcome guidance on the preferred citation,
release identifier, contributor acknowledgement, and whether collaboration or
additional review by the SHaRe team would be appropriate.

No patient-level phenotype, carrier identity, or individual-level data would be
used or redistributed. We will not describe SHaRe and CardioBoost as independent
evidence ecosystems because CardioBoost used earlier SHaRe data.

Thank you for maintaining this valuable resource.

Sincerely,

[Author name and affiliation]
