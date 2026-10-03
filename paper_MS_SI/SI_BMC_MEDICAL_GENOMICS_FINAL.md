# Additional file 1: Supplementary methods and figures

# CardioQueue: development and retrospective multi-source evaluation of machine-learning ranking for cardiogenetic variant review

## Supplementary contents

1. Cohorts, outcomes, and benchmark construction
2. Predictors and model development
3. Primary multi-source evaluation
4. CardioBoost comparison and review queues
5. Robustness, calibration, structure, and FoldX
6. Temporal ClinVar VUS analysis
7. Clinical boundaries and prospective evaluation
8. Supplementary figure legends

Complete numeric tables are provided in Additional file 2. Prediction-level public outputs, schemas, scripts, environments, and checksums will be released in the study repository. HiRO participant-linked records will not be released.

## S1. Cohorts, outcomes, and benchmark construction

### S1.1 Registry and outcome definition

The registry combined normalized records relevant to inherited cardiomyopathy and arrhythmia. The original variant-level table contained 85,677 rows: 42,361 variants of uncertain significance (VUS), 33,923 benign/likely benign (B/LB), 9,148 pathogenic/likely pathogenic (P/LP), and 245 rows without an eligible supervised label. VUS were excluded from primary supervised fitting.

The binary outcome was the classification supplied by the relevant source after broad-label harmonization. Pathogenic and likely pathogenic were combined as P/LP; benign and likely benign were combined as B/LB. Conflicting or nonbinary labels were not forced into a class. Source assertions were preserved separately for provenance. Pooled performance used one row per exact allele only when retained source labels agreed.

All alleles were represented on GRCh38 and checked against the reference genome. Equivalent representations were normalized before overlap testing. Exact chromosome, position, reference, and alternate allele defined identity for development quarantine and pooled evaluation.

### S1.2 Development data

ClinVar supplied development labels. Review status, submitter identifiers, source labels, exact variant identifiers, coordinates, and aggregate classification fields were excluded from predictors. One-star records were retained with lower observation weight; a sensitivity analysis excluded them. The primary fit used 35,796 training rows and 6,317 early-stopping rows after candidate evaluation alleles were quarantined.

The gene-panel audit reviewed inherited-cardiac relevance, aliases, disease mechanism, and obvious non-cardiac label contamination. The final benchmark represented 40 genes. This does not imply that all cardiogenetic genes, structural variants, repeat expansions, deep intronic mechanisms, or complex alleles were evaluated.

### S1.3 Evaluation sources

Five sources contributed benchmark labels.

- **HiRO:** 69 exact alleles, including 27 P/LP and 42 B/LB. These alleles were absent from model development. HiRO is model-external but institutionally internal.
- **eMERGE:** 176 alleles, including 96 P/LP and 80 B/LB. eMERGE is institutionally external.
- **CardioBoost:** 173 strict reference-valid alleles, including 129 P/LP and 44 B/LB. CardioBoost is institutionally external.
- **SHaRE:** 337 of 756 candidate assertions were retained after excluding 419 alleles present in the original training or validation partitions. The retained set included 198 P/LP and 139 B/LB.
- **Fernandez-Falgueras/PLOS:** 88 of 154 candidate assertions were retained after excluding 66 alleles present in the original training or validation partitions. The retained set included 24 P/LP and 64 B/LB.

The retained sources contributed 843 assertions. After exact-allele deduplication, the primary benchmark contained 775 unique concordantly labelled alleles: 436 P/LP and 339 B/LB. Sixty-four alleles occurred in at least two sources; the resulting duplicate assertions remained available for source-specific analyses but counted once in pooled performance.

One MYBPC3 allele, GRCh38 11-47349804-C-G, was P/LP in CardioBoost and B/LB in SHaRE. It was retained in the conflict table and excluded from pooled estimates without majority voting. Source-level performance preserved each source's assertion.

### S1.4 Relationship to ClinVar

Model-held-out status means that an exact allele was absent from the model's training and early-stopping partitions. It does not prove independence from ClinVar, shared laboratories, common literature, or predictor training corpora.

Among the 337 retained SHaRE assertions, 301 were marked manually curated or manually curated/check in the source browser. Seventy-seven exact alleles were absent from the current local ClinVar release, 210 had a nonbinary or discordant ClinVar aggregate, and 50 had a concordant binary aggregate. All 88 retained PLOS alleles had a current ClinVar coordinate record because ClinVar aided coordinate resolution; only 10 had a concordant binary aggregate. Across the complete benchmark, 661 of 775 alleles had a current ClinVar coordinate record and 114 did not. These checks support distinct source adjudications but not independent evidence ecosystems.

## S2. Predictors and model development

### S2.1 Predictor groups

The CardioQueue primary model used 300 public predictors from the following groups:

- sequence consequence and transcript context;
- population frequency and observed absence;
- conservation;
- splice prediction;
- established variant-effect scores;
- gene and disease-mechanism context;
- protein sequence and amino-acid properties;
- AlphaFold residue confidence and structure availability.

Principal resources included VEP, SpliceAI, gnomAD, dbNSFP, ClinGen, UniProt, AlphaMissense, and AlphaFold. Missing annotation was retained as missing information. No absent score was converted into benign evidence.

Outcome labels, source indicators, review status, private HiRO fields, exact coordinates, variant identifiers, and source-level label summaries were excluded. The public model accepts the same schema for all sources. Exact predictor names, transformations, versions, and provenance will be distributed in the repository.

### S2.2 Model fitting and naming

The **CardioQueue primary model** is the final binary CatBoost fit evaluated on the 775-allele benchmark. It used observation weights and validation-based early stopping. Exact hyperparameters, random seeds, package versions, and the serialized model belong to the reproducibility repository because they are needed for execution but not for clinical interpretation.

The **legacy source-held-out model** is an earlier binary fit evaluated on 430 HiRO, eMERGE, and CardioBoost variants. It is retained only where queue and comparator analyses were already locked to that shared cohort.

The **temporal CardioQueue model** is a separate fit using only January 2024 ClinVar B/LB and P/LP labels. No January 2024 VUS entered that model's fitting or validation data.

An exploratory three-class model was evaluated separately. It did not replace the binary model because the binary score with an explicit deferral zone preserved greater P/LP sensitivity and had a clearer review-prioritization interpretation.

### S2.3 Evaluation and uncertainty

The primary threshold was 0.5. Exploratory review zones used scores <=0.1 for lower-priority review, >0.1 to <0.9 for deferred review, and >=0.9 for accelerated review. These thresholds were borrowed from CardioBoost for descriptive comparison and were not optimized or clinically calibrated for CardioQueue.

The principal discrimination measures were area under the receiver operating characteristic curve (AUROC) and area under the precision-recall curve (AUPRC). Threshold results included sensitivity and specificity. Calibration was assessed separately because the score is intended primarily for ranking. Primary confidence intervals used 2,000 bootstrap samples clustered by gene. Paired comparator analyses clustered by exact variant. Random queue references used 10,000 permutations.

## S3. Primary multi-source evaluation

### S3.1 Combined and source-specific results

In 775 unique alleles, the CardioQueue primary model achieved AUROC 0.9841 (gene-clustered 95% confidence interval [CI], 0.9772 to 0.9890) and AUPRC 0.9881 (0.9769 to 0.9924). At 0.5, sensitivity was 0.9794 (0.9672 to 0.9883) and specificity was 0.8407 (0.7704 to 0.8878). There were 427 true positives, 54 false positives, 285 true negatives, and 9 false negatives.

Source-specific AUROC was 0.9947 in HiRO, 0.9854 in eMERGE, 0.9632 in CardioBoost, 0.9863 in SHaRE, and 0.9772 in PLOS. Specificity ranged from 0.6818 in CardioBoost to 0.9375 in PLOS. This heterogeneity is important: the pooled result summarizes ranking across a constructed source mixture and does not guarantee identical operating characteristics in a new laboratory.

Specificity was lower than in the earlier 430-variant cohort because the expanded benchmark added new B/LB distributions, especially SHaRE and PLOS assertions. The primary model produced 24 false positives in SHaRE and 14 in CardioBoost. Stable AUROC with lower specificity indicates threshold transport and calibration shift. A laboratory should estimate local operating characteristics before using a fixed score boundary.

### S3.2 Multi-source consensus and conflicts

The post hoc consensus subset contained 64 alleles classified concordantly by at least two sources: 36 P/LP and 28 B/LB. AUROC was 0.9931 (95% CI, 0.9531 to 1.000), AUPRC was 0.9949 (0.9579 to 1.000), sensitivity was 1.000, and specificity was 0.8929. The confusion matrix contained 36 true positives, 3 false positives, 25 true negatives, and no false negatives.

This subset tests consistency where more than one source agrees; it is not a separate prospective cohort. Its interval is wide, and contributing sources may share evidence. The single source-discordant MYBPC3 allele was excluded from performance and remains an explicit reconciliation case.

### S3.3 Exploratory three-zone policy

Among 436 P/LP alleles, 405 entered accelerated review, 29 were deferred, and 2 entered lower-priority review. Among 339 B/LB alleles, 19 entered accelerated review, 115 were deferred, and 205 entered lower-priority review. The nondeferred subset covered 81.4% of benchmark alleles and had 96.7% accuracy. The accelerated zone had 95.5% positive predictive value in this source mixture, and the lower-priority zone had 99.0% negative predictive value.

These values are not expected to transport to a low-prevalence VUS backlog. The two P/LP alleles in the lower-priority zone also show why the lower zone cannot be treated as reassuring or permanently excluded from review.

The final primary model assigned 17,380 of 42,361 registry VUS to accelerated review, 15,107 to deferred review, and 9,874 to lower-priority review. The public worklist excluded all HiRO-linked records and retained 41,622 VUS: 17,258 accelerated, 14,727 deferred, and 9,637 lower-priority. VUS outcomes are unknown. These labels describe queue position only.

## S4. CardioBoost comparison and review queues

### S4.1 Executable CardioBoost comparison

The released CardioBoost artifacts were first run on their own reported holdouts. AUROC was 0.907 in 218 cardiomyopathy rows and 0.948 in 154 arrhythmia rows, supporting faithful local execution.

The strict shared comparison required matching gene and cDNA notation plus concordant protein substitution. It contained 235 observations, including 160 P/LP and 75 B/LB, representing 211 unique variants. CardioQueue AUROC was 0.984 versus 0.864 for CardioBoost. The paired difference was 0.120 (95% CI, 0.068 to 0.176). AUPRC was 0.991 versus 0.918, with a difference of 0.074 (0.037 to 0.115).

CardioQueue versus CardioBoost AUROC was 0.986 versus 0.852 among 129 cardiomyopathy observations and 0.982 versus 0.934 among 106 arrhythmia observations. A broader cDNA-only sensitivity contained 250 observations. Neither comparison reproduces CardioBoost development or its SHaRe patient-outcome analysis.

### S4.2 Fixed review budgets

The legacy 430-variant cohort contained 255 P/LP and 175 B/LB variants. The first 86 score-ranked reviews contained 86 P/LP variants. Reviewing 215 variants recovered 211 P/LP variants (82.7%) with 98.1% precision.

Recovering 80% of P/LP variants required 207 CardioQueue-ranked reviews, compared with a random-order median of 344 and an oracle minimum of 204. At 95% recovery, CardioQueue required 255 reviews, compared with a random median of 410 and an oracle minimum of 243.

At a 20% review budget, CardioQueue recall modestly exceeded VEP consequence, CADD, REVEL, direct AlphaMissense, and gnomAD rarity on each comparator's available subset. Coverage differed substantially across comparators. The small advantage over direct AlphaMissense makes broader consequence coverage a more credible operational advantage than universal superiority on shared missense variants.

Queue results depend on prevalence and case spectrum. The legacy cohort was 59.3% P/LP, far above a routine VUS backlog. Random permutation is an unprioritized reference, not observed first-come laboratory practice. No reviewer time, evidence acquisition cost, or classification-change endpoint was measured.

## S5. Robustness, calibration, structure, and FoldX

### S5.1 Label quality and leakage checks

A sensitivity fit excluded all one-star ClinVar development records, leaving 12,635 training and 2,308 validation rows. In the fixed 430-variant cohort, AUROC was 0.984 versus 0.978 for the corresponding primary-schema fit. This supports persistence of ranking under a higher-review-status label restriction; it does not remove ClinVar ecosystem overlap.

Matched ablations removed gene identity, source-record summaries, consequence/HGVS fields, and combinations of possible provenance features. Discrimination remained high, although specificity fell in the strictest ablation. Removing 58 established effect-predictor fields changed legacy-cohort AUROC by -0.00114 (95% CI, -0.00583 to 0.00399). In the broader 1,233-allele benchmark, the same ablation changed AUROC by +0.00317 (0.00039 to 0.00601). The direction varied by cohort, arguing against a single established score family explaining performance.

Earlier development artifacts included 38 HiRO-prefixed fields, all with zero feature importance. They were removed from the public schema before the primary model was finalized. No private HiRO field enters primary inference.

### S5.2 Calibration and subgroup findings

In the legacy external mixture, the calibration intercept was -0.642 and slope 0.694. The score was upward-shifted and too extreme for that source distribution. Removing class upweighting improved point calibration with a small loss of discrimination, but no weighting scheme was selected in an independent calibration cohort. CardioQueue scores should be treated as relative priorities.

The missense subset contained 321 rows and had AUROC 0.986. A small subset of 32 gnomAD-confirmed absent variants had AUROC 0.804 with a wide 95% CI of 0.609 to 0.968. This identifies a difficult sparse stratum and does not support a stable point estimate.

### S5.3 Structure geometry

Reference-checked mapping joined reviewed UniProt accessions and substitutions to AlphaFold structures. A mapping was accepted only when the reference amino acid agreed across variant annotation, UniProt, and the selected AlphaFold residue. DSSP and FreeSASA supplied secondary-structure and solvent-accessibility measures.

Registry-wide mapping produced complete AlphaFold pLDDT, DSSP, and FreeSASA values for 38,021 variants, including 2,889 legacy training, 493 validation, and 271 held-out rows. The structure-enhanced model reached AUROC 0.979 compared with 0.977 for its matched baseline in the legacy cohort. The gene-clustered difference was 0.0025 (95% CI, -0.00003 to 0.00502). In 5,890 variants from 19 genes absent from training, the difference was 0.00064 (-0.00003 to 0.00246). Structure geometry therefore showed a small source-held-out signal but no clear unseen-gene transfer gain.

### S5.4 FoldX delta-delta-G

FoldX 5.1 repaired 87 required AlphaFold fragments and returned numeric signed and absolute delta-delta-G for all 4,137 unique structure-mutation calculations. Coverage included 3,031 legacy training, 509 validation, and 717 expanded-external variant rows. In the primary 775-allele benchmark, 382 variants were eligible at pLDDT >=70, 441 at >=50, and 529 at any pLDDT.

The matched experiment separated three signals: baseline predictors, eligibility for reliable structural mapping, and eligibility plus numeric delta-delta-G. At pLDDT >=70, AUROC was 0.98393 for baseline, 0.98548 for eligibility, and 0.98298 after adding delta-delta-G. Relative to eligibility, delta-delta-G changed AUROC by -0.00233 (95% CI, -0.00390 to -0.00109). At pLDDT >=50, delta-delta-G modestly improved on baseline but did not clearly outperform eligibility. Results using all mapped residues were also not superior to eligibility.

These analyses maximize available FoldX coverage while avoiding the claim that computability itself represents stability information. They support release of delta-delta-G as a mechanistic annotation for eligible missense variants, but not its inclusion as a demonstrated source of primary-model improvement.

## S6. Temporal ClinVar VUS analysis

### S6.1 Design

Registry alleles were matched to the January 2024 ClinVar archive by exact GRCh38 allele. Broad significance was retained only when mapped terms agreed. The temporal CardioQueue model used 23,851 training and 4,210 validation variants already classified B/LB or P/LP at baseline. Every January 2024 VUS was excluded from fitting.

The model then scored all 25,481 scope-eligible January 2024 VUS once. At the July 2026 comparison snapshot, 146 were P/LP, 309 were B/LB, 24,879 remained VUS, and 147 had no unambiguous binary or VUS class. Still-VUS records were kept as outcome-unobserved and were not counted as benign.

Three sensitivity fits progressively removed established effect predictors, dynamic annotations, and other fields vulnerable to temporal circularity. The most restricted model used static molecular, gene, conservation, and consequence features. Even that model used annotations materialized after baseline and was not fully time-locked.

### S6.2 Resolved-subset and complete-queue estimands

Among the 455 VUS with an observed later binary class, the full temporal model reached AUROC 0.997 and AUPRC 0.996. This resolved subset is highly selected: variants that receive more scrutiny or new evidence are more likely to acquire a binary outcome.

The primary temporal estimand used the complete 25,481-record queue. The first 500 ranked reviews contained 16 of 146 observed later-P/LP variants; the first 5,000 contained 93. Recovering 80% required 6,644 reviews, compared with a random-order median of 20,313. Sensitivity fits retained substantial enrichment.

These results answer whether later-observed P/LP classifications were enriched toward the front of a historical VUS queue. They do not estimate the probability that any VUS will be reclassified, because most outcomes remain unobserved and reclassification is not random. A prospective claim requires archived baseline annotations, a locked model, fixed review rules, and outcomes collected after deployment.

## S7. Clinical boundaries and prospective evaluation

### S7.1 Supported use

CardioQueue may support scheduled variant-level evidence review in inherited-cardiac genetics. Potential tasks include:

- selecting aging VUS for literature and database refresh;
- identifying records for transcript or coordinate reconciliation;
- prioritizing segregation or functional follow-up when capacity is limited;
- preparing variants for multidisciplinary discussion;
- auditing high-scoring B/LB or low-scoring P/LP records for discordant evidence.

The model should be combined with a clinician-defined urgent pathway. Patient phenotype, arrhythmia burden, ventricular function, syncope, family history, and time-sensitive management questions can override molecular queue order.

### S7.2 Unsupported use

The score must not independently:

- assign an ACMG/AMP classification;
- diagnose or exclude inherited cardiac disease;
- rank patients by urgency;
- trigger treatment, surveillance, cascade testing, or reproductive decisions;
- permanently defer a low-scoring variant;
- replace phenotype, inheritance, segregation, functional, or expert evidence.

A maximum review interval and periodic unranked safety sweep are required until prospective evidence demonstrates that ranked review does not delay clinically important variants.

### S7.3 Proposed prospective study

A feasible evaluation would enroll VUS eligible for routine scheduled reinterpretation and compare CardioQueue-ranked review with random ordering or usual practice. Two cardiogenetic reviewers could independently review prespecified batches, blinded where feasible to allocation and prior model score during final ACMG/AMP adjudication.

Outcomes should include completed reviews per unit time, new evidence found, classification changes, inter-reviewer agreement, time to amended report, downstream segregation or functional studies, clinician overrides, and important variants delayed beyond a safety interval. The model, annotations, thresholds, analysis plan, and review budget should be locked before outcomes accrue. This design would test clinical workflow value; the present retrospective study does not.

## S8. Supplementary figure legends

**Figure S1. Source records and variant-level modeling rows.** Source-record counts are compared with normalized variant rows. Repeated HiRO genotypes can represent distinct participant records, so source records and exact alleles are not interchangeable units.

**Figure S2. Predictor-family coverage.** Coverage is shown for VEP, SpliceAI, gnomAD, dbNSFP, AlphaMissense, and original structure joins. Missingness is expected because several resources apply only to particular consequence classes. The completed registry-wide structure and FoldX analyses are reported in Section S5.

**Figure S3. Calibration reliability.** Predicted CardioQueue score is compared with observed P/LP fraction in validation and legacy external strata. The external intercept of -0.642 and slope of 0.694 indicate source-dependent miscalibration. The score is used for ranking, not as an absolute probability.

**Figure S4. Stratified legacy-cohort performance.** AUROC, AUPRC, sensitivity, and specificity are shown across selected gene, consequence, annotation-coverage, and rarity strata. Small strata have wide uncertainty and are hypothesis-generating.

**Figure S5. Historical-VUS cohort and bidirectional queues.** Panel A accounts for all 25,481 eligible January 2024 VUS. Panels B and C show recovery of observed later-P/LP from the high-score end and observed later-B/LB from the low-score end. Still-VUS records are outcome-unobserved. Predictor resources were not uniformly frozen to January 2024.
