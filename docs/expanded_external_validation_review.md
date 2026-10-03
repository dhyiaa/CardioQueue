# Expanded External Evaluation Review

## Bottom line

The final primary benchmark contains 775 unique, reference-valid GRCh38 alleles from five source strata: 436 P/LP and 339 B/LB variants across 40 genes. It retains HiRO, eMERGE, and CardioBoost because they were source-held-out from model development, then adds only PLOS and SHaRe alleles absent from the original training and validation partitions. The resulting 300-feature CatBoost model reached AUROC 0.9841 (gene-clustered 95% CI, 0.9772-0.9890), AUPRC 0.9881 (0.9769-0.9924), sensitivity 0.9794, specificity 0.8407, and Brier score 0.0629. A post hoc sensitivity among 64 cross-source-concordant alleles reached AUROC 0.9931 (0.9531-1.000), sensitivity 1.000, and specificity 0.8929.

This is best described as a **retrospective multi-source model-held-out evaluation with exact-allele quarantine**. HiRO is model-external but institutionally internal; eMERGE, CardioBoost, PLOS, and SHaRe are institutionally external. The evaluation is materially stronger than selecting another ClinVar sample, but it is not prospective clinical validation and is not fully independent of the public ClinVar/literature evidence ecosystem.

### Were the test variants in the ClinVar training data?

**Not in the data used to fit the accepted expanded model.** Before retraining, all 1,233 candidate alleles were assigned to the external partition. The accepted model's 35,796 training rows and 6,317 early-stopping validation rows contain zero exact candidate or primary-benchmark alleles.

This does **not** mean that every benchmark allele was absent from ClinVar or had never appeared in an earlier development table. Of the 1,233 final alleles, 942 occurred in the original 85,677-row matrix and 291 were newly appended. Before the redesigned quarantine, 381 of the final benchmark alleles had been in the original training split and 77 in the original validation split. These 458 alleles were removed, the development split was rebuilt, and the model was trained again from scratch. The remaining original-matrix benchmark alleles were already source-held-out or not development-eligible.

At source level, 419 of 756 final SHaRe observations and 66 of 154 final PLOS observations had occupied the earlier train/validation splits; these source-level counts overlap and should not be summed as unique alleles. The legacy HiRO, eMERGE, and CardioBoost evaluation alleles were already in source-held-out strata. Therefore, the defensible statement is **“no exact evaluation allele was used for fitting or early stopping”**, not **“the evaluation variants were never present in ClinVar.”** Public predictor training sets and source classifications may also contain the same variants or evidence, which is why evidence-ecosystem independence is not claimed.

## What was external before expansion

The prior 430-row benchmark comprised 69 HiRO, 176 eMERGE, and 185 CardioBoost variants after source precedence. Its sources were not equivalent:

| Source | Prior rows | External status | Interpretation |
|---|---:|---|---|
| HiRO/CASPER WES/VERDICT | 69 | Model-external; institutionally internal | Unpublished local patient-linked research cohort; no HiRO allele or label informed model development. |
| eMERGE-III arrhythmia | 176 | Yes, as a public data source | External published cohort with post-functional-study labels; retrospective and connected to public evidence. |
| CardioBoost | 185 | Yes, as a public data source | External published model dataset; older SHaRe evidence and ClinVar/literature overlap prevent evidence-ecosystem independence. |

Strict GRCh38 FASTA checking retained all 69 HiRO and 176 eMERGE variants but only 174 of 185 legacy CardioBoost variants. Eleven CardioBoost coordinates failed reference-allele validation and were excluded from the new strict benchmark. One additional MYBPC3 allele (`11-47349804-C-G`) was excluded from the combined benchmark because SHaRe called it B/LB while legacy CardioBoost called it P/LP.

The historical 430-row analyses remain valid descriptions of the earlier frozen workflow. They should not be relabeled as a fully external 430-variant clinical validation set.

## New external sources

### SHaRe HCM 2026Q1

- Origin: public SHaRe HCM variant browser repository snapshot, file `HCM_Gen_Variant_Annotation_2026Q1.csv`.
- Input checksum: `f119e47ea6205d874c0dfe3b3a350617509d10a906a91fe1f998abc98eb53a45`.
- Starting binary records: 822.
- Retained after panel, coordinate, reference, and class QC: 757 (535 P/LP, 222 B/LB).
- Source-specific broad evaluation after excluding the single cross-source label conflict: 756 (535 P/LP, 221 B/LB); final primary evaluation: 337 (198 P/LP, 139 B/LB) absent from the original train/validation partitions.
- External status: institutionally external HCM registry curation, but not independent of CardioBoost or the ClinVar/literature evidence ecosystem.
- Publication constraint: the repository has no permissive reuse license and states all rights reserved. Written permission is required before publishing or redistributing a derived SHaRe manifest.

Forty-four binary rows were outside the current gene panel and 21 failed coordinate/reference requirements. No uncertain class was converted to a binary outcome. SHaRe stores its final classification separately from its ClinVar annotation and flags disagreements. Among the 337 eligible observations, 301 were marked manually curated or manually curated/check; 77 exact alleles were absent from the current local ClinVar release, 210 had nonbinary or discordant ClinVar aggregates, and 50 had concordant binary aggregates. SHaRe is therefore not a ClinVar-label export, although ClinVar and shared literature can contribute evidence.

### Fernandez-Falgueras et al. PLOS ONE 2024

- Origin: inherited cardiovascular disease reinterpretation supplement, DOI `10.1371/journal.pone.0297914`.
- License: CC BY.
- Input checksum: `79eae13f642c14ddadbef5913728745a32315f0164d1f3b6be8ded8e3dd5fbce`.
- Starting panel-matched binary records: 197.
- Exact GRCh38 alleles retained broadly: 154 (74 P/LP, 80 B/LB); final primary evaluation: 88 (24 P/LP, 64 B/LB) absent from the original train/validation partitions.
- External status: external, manually reinterpreted clinical variant series with reported ACMG criteria.

The supplement supplied gene and HGVS descriptions but not analysis-ready genomic alleles. ClinVar was used only as a coordinate resolver: a record was retained only when gene plus cDNA HGVS produced one unambiguous GRCh38 allele and the reference base matched the GRCh38 FASTA. ClinVar did not supply the outcome label. The PLOS investigators performed a 2022 ACMG/AMP reinterpretation using updated specifications and population, computational, case, segregation, functional, literature, and institutional evidence. Among the 88 eligible observations, only 10 had a concordant binary current ClinVar aggregate; 78 were conflicting, uncertain, or otherwise nonbinary. Forty-three records without a unique mapping were left out rather than inferred.

## Final 775-variant primary benchmark

| Source-specific stratum | Rows | P/LP | B/LB | Role |
|---|---:|---:|---:|---|
| HiRO | 69 | 27 | 42 | Model-external, institutionally internal |
| eMERGE | 176 | 96 | 80 | Prior source-held-out, institutionally external |
| CardioBoost, strict subset | 173 | 129 | 44 | Prior source-held-out, institutionally external |
| SHaRe eligible subset | 337 | 198 | 139 | Institutionally external; absent from original train/validation |
| PLOS eligible subset | 88 | 24 | 64 | Institutionally external; absent from original train/validation |
| Combined exact-allele-deduplicated | 775 | 436 | 339 | Primary model-held-out evaluation |

Source rows overlap and therefore exceed the 775 unique alleles. The primary manifest is `datasets/external_validation_candidates/final_model_heldout_2026/external_unique_variant_manifest.tsv`; its eligibility audit records every included and excluded source observation. The combined benchmark reached AUROC 0.9841, AUPRC 0.9881, Brier score 0.0629, sensitivity 0.9794, and specificity 0.8407. All 775 alleles were absent from the original and accepted development partitions.

### How repeated alleles are represented

The data use a two-level design. `external_source_observations.tsv` preserves every source assertion, allowing the same exact allele to have separate rows and provenance in multiple cohorts. `external_unique_variant_manifest.tsv` contains one row per exact allele for the pooled test. Concordant repeated assertions therefore count once in combined performance but remain available for source-specific results. Discordant assertions remain in the conflict audit and are excluded from the pooled primary metric rather than resolved by majority vote. Repeating the same allele, feature vector, and model score in a pooled metric would overweight well-reported variants without adding an independent test case.

The 843 retained source assertions represented 775 unique alleles. Sixty-four alleles occurred in multiple retained sources, accounting for 68 duplicate observations. Source-exclusive counts were 62 for HiRO, 172 for eMERGE, 132 for CardioBoost, 66 for PLOS, and 279 for SHaRe. The final retained assertions were class-concordant; the one detected discordant SHaRe/CardioBoost MYBPC3 allele was excluded and documented separately.

### Six-part result summary

| Analysis layer | Result | Interpretation |
|---|---|---|
| Source assertions | 843 rows: 474 P/LP and 369 B/LB | Preserves every eligible source classification and provenance |
| Unique-allele table | 775 alleles: 436 P/LP and 339 B/LB across 40 genes | Primary unit for pooled performance; no exact allele in fitting or early stopping |
| Conflict table | One MYBPC3 allele: P/LP in CardioBoost and B/LB in SHaRe | Preserved for reconciliation and excluded from pooled accuracy |
| Source-specific tests | AUROC 0.9947 HiRO, 0.9854 eMERGE, 0.9632 CardioBoost, 0.9772 PLOS, 0.9863 SHaRe | Shows performance was not confined to one contributing source |
| Combined primary test | AUROC 0.9841, AUPRC 0.9881, sensitivity 0.9794, specificity 0.8407, Brier 0.0629 | Main retrospective model-held-out estimate |
| Multi-source consensus sensitivity | 64 alleles across 17 genes; AUROC 0.9931, AUPRC 0.9949, sensitivity 1.000, specificity 0.8929, Brier 0.0465 | Post hoc support for ranking consistency when at least two retained sources agreed |

The consensus confusion matrix was TP 36, FP 3, TN 25, and FN 0. Gene-clustered 95% intervals were 0.9531-1.000 for AUROC, 0.9579-1.000 for AUPRC, 1.000-1.000 for sensitivity, and 0.6471-1.000 for specificity. The subset is small, post hoc, and not evidence-ecosystem independent, so it supports robustness rather than constituting a second external validation.

### Clinical meaning

For clinicians, the useful signal is workflow confidence rather than a new diagnostic threshold. A genetics service can use the score to order literature refresh, segregation requests, transcript review, functional-study selection, or multidisciplinary discussion. Cross-source agreement can identify records with more stable classification context, while the conflict table prevents disagreement from being hidden inside an average performance estimate. CardioQueue still does not reclassify a VUS, establish that a variant explains the patient's phenotype, or determine surveillance, treatment, or family testing.

## Broader 1,233-variant sensitivity benchmark

The source files yielded 1,330 reference-valid source observations before conflict handling. The combined manifest deduplicates exact chromosome-position-reference-alternate alleles and excludes the one discordant cross-source allele, leaving 1,233 unique variants. Source-specific analyses retain one row per source assertion, so their row counts overlap and must not be summed as unique variants.

| Source-specific stratum | Rows | P/LP | B/LB |
|---|---:|---:|---:|
| SHaRe HCM 2026Q1 | 756 | 535 | 221 |
| Fernandez-Falgueras PLOS 2024 | 154 | 74 | 80 |
| Legacy CardioBoost, strict reference-valid subset | 173 | 129 | 44 |
| Legacy eMERGE | 176 | 96 | 80 |
| Legacy HiRO | 69 | 27 | 42 |
| Combined unique benchmark | 1,233 | 801 | 432 |

Of the 1,233 variants, 942 were already represented in the prior modeling matrix and 291 were newly normalized and annotated. The complete variant list is `datasets/external_validation_candidates/expanded_external_2026/external_unique_variant_manifest.tsv`. Source-level labels, reported HGVS, mapping method, reference check, and panel status are in `external_source_observations.tsv`; excluded cross-source disagreements are in `external_label_conflicts.tsv`.

## Feature completion

The 291 new alleles were annotated with the same feature families used by the public model. Ensembl REST VEP was used because both existing local VEP installations crashed in their mixed-architecture Perl environments. All final columns were coerced to the original 300-feature schema: 203 categorical and 97 numeric predictors.

| Feature family | Covered benchmark variants | Coverage |
|---|---:|---:|
| VEP | 1,233 | 100.0% |
| SpliceAI | 1,034 | 83.9% |
| dbNSFP | 908 | 73.6% |
| gnomAD observed allele | 892 | 72.3% |
| AlphaMissense | 687 | 55.7% |
| Reviewed protein mapping | 418 | 33.9% |
| AlphaFold pLDDT | 408 | 33.1% |
| Existing FoldX result | 167 | 13.5% |

All 291 new alleles received VEP annotation; 224 received SpliceAI, 195 dbNSFP, 170 AlphaMissense and protein mapping, and 149 an observed gnomAD value. AlphaFold/DSSP/FreeSASA processing completed for all 35 required structure fragments and all 170 mapped new substitutions. FoldX was not computed for the new alleles, and its values remain explicitly unavailable rather than imputed. The expanded model therefore does not establish a FoldX DDG benefit.

## Broader sensitivity performance

The exact-allele-quarantined matrix contained 35,796 training rows, 6,317 early-stopping validation rows, and 1,233 external benchmark rows. No exact benchmark allele appeared in development. The unchanged public 300-feature recipe stopped after 314 trees.

| Evaluation | Rows | AUROC | AUPRC | Brier | Sensitivity | Specificity |
|---|---:|---:|---:|---:|---:|---:|
| Combined unique | 1,233 | 0.9841 | 0.9917 | 0.0513 | 0.9688 | 0.8657 |
| PLOS 2024 | 154 | 0.9813 | 0.9827 | 0.0540 | 0.8784 | 0.9500 |
| SHaRe 2026Q1 | 756 | 0.9858 | 0.9945 | 0.0466 | 0.9682 | 0.8733 |
| CardioBoost strict subset | 173 | 0.9632 | 0.9861 | 0.0707 | 0.9845 | 0.6818 |
| eMERGE | 176 | 0.9854 | 0.9872 | 0.0524 | 0.9792 | 0.8750 |
| HiRO | 69 | 0.9947 | 0.9926 | 0.0518 | 1.0000 | 0.8810 |
| Newly annotated variants | 291 | 0.9858 | 0.9900 | 0.0565 | 0.9740 | 0.8832 |

The combined gene-clustered 95% CIs were 0.9780-0.9882 for AUROC, 0.9790-0.9947 for AUPRC, 0.9575-0.9850 for sensitivity, and 0.8119-0.9088 for specificity. Clustering by gene avoids treating many variants in the same gene as fully independent.

## Leakage and robustness checks

Exact-allele quarantine does not remove every molecular relationship. A stricter sensitivity analysis also removed 437 development rows sharing either the same genomic site or the same VEP protein substitution with a benchmark variant. It retained 35,428 training and 6,248 validation rows and reached AUROC 0.9825, AUPRC 0.9911, Brier 0.0524, sensitivity 0.9613, and specificity 0.8657. The small decline supports robustness while showing that exact-allele separation should not be described as complete biological independence.

A second model removed 58 precomputed effect-predictor fields, including SIFT, PolyPhen, REVEL, MetaLR, FATHMM, CADD, AlphaMissense, and ESM1b families. Its 242-feature model reached AUROC 0.9873, AUPRC 0.9933, and Brier 0.0400. Relative to the primary model, the gene-clustered paired differences were +0.00317 for AUROC (95% CI, +0.00039 to +0.00601), +0.00165 for AUPRC (-0.00014 to +0.00526), and -0.01132 for Brier (-0.02070 to -0.00663). The result argues that performance is not merely inherited from these established predictor families; it does not eliminate overlap through remaining public databases or curation evidence.

## Reviewer-safe interpretation

Supported claims:

- CardioQueue generalizes across 42 evaluated cardiogenetic genes and multiple consequence classes in a 1,233-allele retrospective benchmark.
- Exact held-out alleles did not enter training or early stopping, and a stricter same-site/protein-substitution quarantine preserved high discrimination.
- The model combines a broader set of genomic, protein, and structure annotations than missense-only CardioBoost and can rank records outside rare missense variation.
- Performance persisted, and improved modestly, after removing major precomputed pathogenicity predictors.
- A post hoc 64-allele cross-source-concordant sensitivity retained high discrimination, with wide uncertainty appropriate to its size.
- The model can support evidence-review ordering, reinterpretation triage, discordance audits, and selection of variants for segregation or functional follow-up.

Claims to avoid:

- “Independent external clinical validation,” because HiRO is internal and public sources share evidence ecosystems.
- “Independent replication in 64 variants,” because the consensus subset was post hoc and contributing sources may share evidence.
- “General to all cardiogenetic mutations,” because the primary benchmark spans 40 genes rather than every gene, mechanism, ancestry, or structural-variant class.
- “Reclassifies VUS” or “assigns ACMG pathogenicity,” because the score supplies no case-level phenotype, segregation, or complete ACMG evidence.
- “FoldX improves performance,” because new FoldX values were not calculated and development coverage remains insufficient.
- “Better than all other models,” because only matched CardioBoost and routine-score comparisons support direct superiority claims.

## Reproducibility map

| Artifact | Purpose |
|---|---|
| `datasets/external_validation_candidates/build_strict_external_nondevelopment_benchmark.py` | Apply final source eligibility, exact-allele deduplication, performance estimation, and gene-clustered bootstrap. |
| `datasets/external_validation_candidates/audit_final_external_label_provenance.py` | Audit current ClinVar presence and SHaRe/PLOS label provenance without replacing source outcomes. |
| `datasets/external_validation_candidates/final_model_heldout_2026/external_unique_variant_manifest.tsv` | One row per allele in the final 775-variant primary benchmark. |
| `datasets/external_validation_candidates/final_model_heldout_2026/external_source_observations.tsv` | One row per retained source assertion. |
| `datasets/external_validation_candidates/final_model_heldout_2026/external_label_conflicts.tsv` | Cross-source binary disagreements excluded from pooled performance. |
| `datasets/external_validation_candidates/final_model_heldout_2026/eligibility_audit.tsv` | Inclusion and exclusion decision for each candidate source observation. |
| `datasets/external_validation_candidates/final_model_heldout_2026/source_accounting.tsv` | Reproducible source flow, overlap, and source-exclusive counts. |
| `datasets/external_validation_candidates/final_model_heldout_2026/source_overlap_combinations.tsv` | Exact combinations represented among overlapping alleles. |
| `results/model_performance/final_model_heldout_2026/` | Primary, source-specific, consensus, and gene-clustered results. |
| `datasets/external_validation_candidates/build_expanded_external_benchmark.py` | Normalize sources, validate GRCh38 alleles, deduplicate, quarantine, and build the expanded matrix. |
| `datasets/external_validation_candidates/expanded_external_2026/external_unique_variant_manifest.tsv` | One row per final benchmark allele. |
| `datasets/external_validation_candidates/expanded_external_2026/external_source_observations.tsv` | Every retained source assertion and mapping decision. |
| `datasets/external_validation_candidates/expanded_external_2026/new_variants_requiring_annotation.tsv` | The 291 newly added alleles. |
| `datasets/external_validation_candidates/expanded_external_2026/modeling_table_expanded_external_annotated.tsv` | Final 85,968-row annotated matrix. |
| `results/models/cardioqueue_v1_expanded_external_2026_public300/` | Accepted 300-feature model. |
| `results/model_performance/expanded_external_2026/` | Broader 1,233-variant sensitivity predictions, performance, feature coverage, and gene-clustered intervals. |
| `results/model_performance/expanded_external_2026_molecular_quarantine/` | Same-site/protein-substitution quarantine sensitivity. |
| `results/model_performance/expanded_external_2026_no_precomputed_effect_predictors/` | 242-feature predictor-family ablation. |
