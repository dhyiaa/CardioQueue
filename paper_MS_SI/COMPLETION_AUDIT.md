# BMC Medical Genomics completion audit

## Deliverables

| Item | File | Status |
|---|---|---|
| Main manuscript | `MS_BMC_MEDICAL_GENOMICS_FINAL.md` | Complete pending author fields |
| Supplementary information | `SI_BMC_MEDICAL_GENOMICS_FINAL.md` | Complete |
| Ten-paper comparator review | `literature_evidence/BMC_MEDICAL_GENOMICS_COMPARATOR_REVIEW.md` | Complete |
| Novelty positioning memo | `literature_evidence/NOVELTY_POSITIONING.md` | Complete |
| Machine-readable tables | `tables/` | Tables S1-S23 and lettered companion tables present |
| Main and supplementary figures | `figures/` | Figures 1-5 and S1-S5 present |
| Analysis environment | `environment-manuscript.txt` | Present |
| Artifact checksums | `ARTIFACT_CHECKSUMS.sha256` | 219 core, sensitivity, structure, FoldX, temporal-VUS, calibration, review-queue, and final multi-source files verified |
| Rollback checkpoint | `../backups/pre_publishability_code_audit_20260917/` | 76 archived artifacts with snapshot checksums |
| Reframing checkpoint | `../backups/pre_review_queue_reframing_20260917/` | Pre-queue manuscript package with snapshot checksums |

## Mechanical checks

| Check | Result |
|---|---:|
| Structured abstract words | 333 |
| Main manuscript words | 11,110 |
| Supplement words | 15,642 |
| Main tables | 5 |
| Main figures | 5 |
| Supplementary table files | 64 |
| Figure files | 14 |
| Markdown table-width errors | 0 |
| Missing cited artifact paths | 0 |
| Em dash or en dash | 0 |
| Prohibited stock phrases | 0 |

## Scientific checks

- Counts and metrics match the frozen rescue-aware matrix and model predictions.
- The final primary benchmark contains 775 exact-allele-deduplicated variants from five source strata, with zero exact overlap in fitting or early stopping; all 843 retained source assertions remain traceable.
- The post hoc 64-allele cross-source-concordant sensitivity reached AUROC 0.9931 (gene-clustered 95% CI, 0.9531-1.000); its limited size, 17-gene scope, and shared evidence ecosystems are disclosed.
- One discordant SHaRe/CardioBoost MYBPC3 allele is preserved in the conflict table and excluded from pooled accuracy without majority voting.
- The 5,000-sample external bootstrap reran deterministically with seed 20260916.
- The ClinVar >=2-star sensitivity model used 12,635 training and 2,308 validation rows; its paired external bootstrap used the same seed and 5,000 samples.
- CardioBoost headline intervals use variant-cluster resampling for 219 unique variants.
- Consensus filtering, source overlap, post hoc eMERGE rescue, calibration error, prevalence dependence, and possible ClinVar ecosystem circularity are disclosed.
- The frozen schema's 38 HiRO source-summary fields all had zero importance; a 300-feature ablation preserved discrimination. The public inference schema is directed to omit them.
- The frozen model's source-only structure fields had no development coverage. A post hoc registry-wide repair produced reference-validated AlphaFold, DSSP, and FreeSASA features for 2,889 training, 493 validation, and 271 held-out rows.
- The 351-feature sensitivity model improved source-held-out AUROC from 0.977 to 0.979; its gene-clustered difference was 0.0025 (95% CI, -0.00003 to 0.00502). Source heterogeneity and post hoc status are disclosed.
- In 5,890 variants from 19 genes absent from training, structure enhancement changed AUROC by 0.00064 (gene-cluster 95% CI, -0.00003 to 0.00246) and AUPRC by -0.00142 (-0.00963 to 0.00542). The absence of a clear gene-transfer gain is reported.
- A calibrator fitted only on internal validation rows reduced external Brier score by 0.00433, but its gene-cluster interval (-0.00919 to 0.00113) included no improvement. The sensitivity-specificity trade-off is disclosed.
- The held-out review queue recovered 211 of 255 P/LP variants in the first 215 reviews. Reaching 80% recovery required 207 score-ranked reviews versus a median 344 under 10,000 random orderings.
- The temporal fit used 23,851 scope-eligible January 2024 binary training rows and 4,210 validation rows; no January 2024 VUS entered fitting. All 25,737 historical VUS received one score from each of four feature sets. The primary 25,481-record queue retained 146 observed later-P/LP transitions, and still-VUS records remained outcome-unobserved rather than being labelled benign.
- In the complete primary VUS queue, the full model recovered 16/146 observed later-P/LP transitions in the first 500 reviews and 93/146 in the first 5,000. The corresponding gene-clustered 95% intervals were 5.6%-23.7% and 45.9%-90.5%. Resolved-subset and full-queue estimands are reported separately.
- Random ordering is identified as an arrival-order proxy, not observed first-come performance. No analyst-time, emergency-triage, or patient-outcome claim is made.
- The release-ready worklist contains 41,622 public-source VUS across 53 genes and excludes every HiRO-linked row. Its scores are labeled as research priorities, not classifications.
- The completed outcome-blind FoldX analysis produced 4,137/4,137 unique numeric calculations across 87 repaired fragments. The pLDDT >=70 DDG model did not improve over its eligibility-only control; the small >=50 gain over baseline did not persist against that control. DDG is not credited with the primary or geometry-enhancement result.
- The 0.1/0.9 zones are labeled exploratory and are not described as calibrated clinical probabilities.

## Author completion required

- Authors, affiliations, and corresponding-author details.
- Ethics committee, protocol number, approval date, and consent or waiver basis.
- Consent-for-publication statement.
- Repository URL or DOI, commit, licenses, and controlled-access procedure.
- Competing interests, funding, CRediT-style contributions, and acknowledgements.
- Final author verification of private-cohort descriptions and all citations.

Author-year citations were retained at the user's request. The submission system may request conversion to numbered Vancouver style.
