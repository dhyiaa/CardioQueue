# CardioQueue: development and retrospective multi-source evaluation of machine-learning ranking for cardiogenetic variant review

## Background

Inherited cardiomyopathies and arrhythmia syndromes are genetically heterogeneous. Pathogenic variation may alter sarcomeric, desmosomal, ion-channel, calcium-handling, metabolic, or nuclear-envelope proteins, and the clinical effect depends on gene, disease mechanism, inheritance, penetrance, and phenotype. Genetic results can guide diagnosis, family screening, and prevention, but only when molecular evidence is interpreted with the patient's clinical and family context. The ACMG/AMP framework therefore combines population, computational, functional, segregation, and clinical evidence rather than treating any prediction as a classification ([Richards et al., 2015](https://doi.org/10.1038/gim.2015.30)).

Cardiogenetic testing also creates a persistent reinterpretation workload. Variants of uncertain significance (VUS) cannot direct care by themselves, yet their evidence changes as population databases, functional studies, segregation data, transcript annotations, and disease-specific knowledge mature. Laboratories and inherited-cardiac programs must decide which older records warrant renewed review when expert capacity is limited. In our registry, that problem involved 42,361 VUS.

Most computational tools address a different question. CADD, REVEL, SpliceAI, EVE, and AlphaMissense estimate general or consequence-specific molecular effects. CardioClassifier and CardioVAI assemble parts of ACMG/AMP interpretation, while GENESIS and CVD-PP focus on selected cardiac settings. CardioBoost is the closest disease-specific comparator: it predicts rare missense variation in cardiomyopathy and arrhythmia genes and has patient-outcome support from SHaRe ([Zhang et al., 2021](https://doi.org/10.1038/s41436-020-00972-3)). These tools are valuable, but discrimination alone does not tell a service how many records must be reviewed to recover a target proportion of findings.

Protein structure may add mechanistic context, particularly for missense variation. Prior cardiac studies evaluated structural and stability features in SCN5A, KCNQ1, and MYBPC3, and broader methods such as AlphScore and SIGMA use AlphaFold-derived geometry. FoldX delta-delta-G is therefore not novel by itself. Its value in a multi-gene cardiogenetic workflow must be tested against a control for whether a variant could be mapped to a reliable structure at all.

We developed CardioQueue as a review-prioritization system across missense, truncating, splice-relevant, synonymous, untranslated, and other cardiogenetic variants. It produces a continuous score and a ranked worklist; it does not assign ACMG/AMP classifications or rank patients. We evaluated the primary model in 775 model-held-out alleles from five sources, compared it directly with executable CardioBoost models, measured recovery under fixed review budgets, ranked a complete historical VUS backlog, and tested whether structure and FoldX added independent signal. The intended clinical pathway is score-informed evidence review followed by full expert adjudication (**Figure 1**).

## Methods

### Study design and intended use

This retrospective variant-level study followed TRIPOD+AI reporting principles ([Collins et al., 2024](https://www.tripod-statement.org/resources/)). The supervised outcome was a B/LB versus P/LP source classification; VUS did not enter the primary fit. CardioQueue was evaluated for molecular discrimination, transfer across source boundaries, recovery under fixed review capacity, comparison with CardioBoost, and enrichment of historical VUS that later received a binary ClinVar classification.

The intended output is a worklist for clinical laboratory scientists and genetics researchers. It may guide literature refresh, transcript review, segregation requests, functional-study selection, or multidisciplinary discussion. CardioQueue does not rank patients, establish causality, assign ACMG/AMP classifications, or support management without case-level evidence.

### Data sources

ClinVar supplied the development labels. Review status, submitted classifications, source identifiers, and exact variant identifiers were excluded from predictors. One-star records were retained with lower weight, and a sensitivity analysis restricted development to records with at least two review stars.

Five additional sources contributed evaluation labels: the institutionally internal HiRO/CASPER WES/VERDICT cohort; the institutionally external eMERGE-III arrhythmia dataset; CardioBoost's released cardiomyopathy and arrhythmia records; SHaRe HCM classifications; and the Fernandez-Falgueras inherited-cardiovascular reinterpretation cohort. All alleles were normalized to GRCh38 and reference-validated. SHaRe and PLOS outcomes came from their source adjudications, not from the model's aggregate ClinVar label. These sources can nevertheless share literature, ClinVar evidence, or prior laboratory interpretations, so model-held-out does not mean evidence-independent.

HiRO participant-linked data are not public. SHaRe-derived rows will be redistributed only with permission. Source definitions, attrition, and classification provenance are provided in Additional file 1 and the public repository.

### Gene panel and registry

Genes were reviewed for inherited cardiac relevance, mechanism, aliases, and non-cardiac label contamination. The audit and sensitivity set are reported in Additional file 1. Source records were normalized to exact GRCh38 alleles. The original registry contained 85,677 unique rows; 291 newly normalized benchmark alleles were then appended. One MYBPC3 allele with discordant CardioBoost and SHaRe labels was retained in the conflict audit and excluded from pooled performance.

### Annotation and predictors

The primary model used 300 public predictors spanning variant consequence, population frequency, conservation, established effect scores, gene context, protein sequence, and structure availability. VEP, SpliceAI, gnomAD, dbNSFP, ClinGen, UniProt, AlphaMissense, and AlphaFold supplied the principal annotation families. Outcome fields, review status, patient and source identifiers, exact variant identifiers, coordinates, source flags, private HiRO fields, and label aggregates were excluded. CatBoost handled incomplete annotation natively; no missing score was interpreted as benign evidence.

Protein substitutions were accepted only when the accession, position, and reference residue agreed across the variant annotation, reviewed UniProt sequence, and AlphaFold model. DSSP and FreeSASA described local geometry and solvent accessibility. FoldX 5.1 estimated signed and absolute delta-delta-G after structure repair. The prespecified high-confidence analysis required pLDDT >=70, with lower-confidence tiers as sensitivities. Matched comparisons tested baseline predictors, structure eligibility alone, and eligibility plus numeric delta-delta-G. The licensed FoldX executable will not be redistributed.

### Weights and split

Development records were weighted by source confidence and outcome class. The exact union of 1,233 reference-valid candidate benchmark alleles was quarantined before fitting. CardioQueue was trained on 35,796 rows, with 6,317 separate rows used for early stopping. Neither partition contained a quarantined candidate allele.

The primary benchmark retained all eligible HiRO, eMERGE, and CardioBoost alleles and only PLOS and SHaRE alleles absent from the earlier development partitions. After exact-allele deduplication, it contained 775 variants (436 P/LP and 339 B/LB) across 40 genes. HiRO was model-held-out but institutionally internal; the other four sources were institutionally external. The broader 1,233-allele quarantine was a sensitivity analysis.

Source assertions remained separate for provenance and source-specific estimates. Pooled performance counted each concordantly labelled exact allele once. Discordant assertions were reported separately and not resolved by majority vote.

### CatBoost training

We fitted a binary CatBoost classifier with observation weights and validation-based early stopping ([Prokhorenkova et al., 2018](https://doi.org/10.48550/arXiv.1706.09516)). This final fit is termed the **CardioQueue primary model**. Exact hyperparameters, seeds, schema, environment, and the trained artifact will be released with the analysis code. An earlier fit evaluated on 430 held-out variants is termed the **legacy source-held-out model** and is used only for analyses that were completed on that fixed cohort.

### Decision rules and statistics

The primary binary threshold was 0.5. An exploratory policy borrowed CardioBoost's published thresholds: scores <=0.1 indicated lower-priority review, scores >0.1 to <0.9 were deferred, and scores >=0.9 indicated accelerated review. These are workload zones, not calibrated clinical cutoffs; a score of 0.9 is not a 90% probability of pathogenicity.

We report AUROC, AUPRC, sensitivity, specificity, predictive values, Brier score, and deferral where appropriate. Primary 95% confidence intervals used 2,000 gene-clustered bootstrap samples. Paired comparator and sensitivity analyses used clustered resampling matched to their unit of analysis. Calibration was assessed with Brier score, reliability plots, and calibration intercept and slope. Full statistical definitions are in Additional file 1.

### Review-queue simulation

The fixed legacy cohort was ordered by decreasing score. We measured P/LP recovery and precision at prespecified review budgets and the number of reviews required for target recovery. Ten thousand random orderings and the theoretical oracle minimum defined unprioritized and best-possible references. VEP consequence, CADD, gnomAD rarity, REVEL, and AlphaMissense were compared only where each annotation was available.

### CardioBoost benchmark

We first reproduced CardioBoost performance on its released holdouts. The head-to-head cohort required matched gene and cDNA notation plus concordant protein substitution, yielding 235 observations from 211 unique variants. Paired differences used bootstrap resampling clustered by unique variant. A broader cDNA-only match was retained as a sensitivity analysis.

### Secondary analyses

Sensitivity analyses restricted ClinVar label quality, removed class upweighting, removed established effect-predictor families, tested potential leakage features, and evaluated protein structure and FoldX. An exploratory three-class model was also fitted but did not replace the binary ranking model.

For the temporal analysis, a separate **temporal CardioQueue model** was trained only on alleles already classified B/LB or P/LP in the January 2024 ClinVar archive. Every January 2024 VUS was excluded from fitting. The model then scored all 25,481 in-scope historical VUS. Current ClinVar status was used only as a retrospectively observed outcome; still-VUS and nonbinary records remained outcome-unobserved. Because several predictor resources postdated January 2024, this was a label-temporal analysis, not a fully time-locked prospective evaluation.

### Software and generative AI assistance

The public repository will provide the trained model, schema, software environment, analysis code, prediction-level public outputs, and checksums. OpenAI Codex assisted with drafting, language editing, literature organization, and consistency checks. The authors verified all analyses and remain responsible for the manuscript; model-generated content was not used as a study observation or predictor.

## Results

### Registry and cohorts

The registry contained 85,677 unique variant rows, including 42,361 VUS, 33,923 B/LB variants, and 9,148 P/LP variants. After appending newly normalized benchmark alleles and quarantining every candidate evaluation allele before fitting, 35,796 rows remained for training and 6,317 for early stopping (**Figure 1**).

The primary benchmark contained 775 unique, reference-valid GRCh38 alleles across 40 genes: 436 P/LP and 339 B/LB. It comprised 69 HiRO, 176 eMERGE, 173 CardioBoost, 337 SHaRE, and 88 Fernandez-Falgueras/PLOS source assertions. The SHaRE and PLOS subsets included only alleles absent from the original training and validation partitions. HiRO was model-held-out but institutionally internal; the other four sources were institutionally external.

Sixty-four alleles were concordantly classified by at least two sources. Their duplicate assertions were retained in the provenance table but each allele was counted once in pooled performance. One additional MYBPC3 allele had discordant CardioBoost and SHaRE labels and was excluded from pooled estimates without majority voting. Exact overlap between the primary benchmark and model development was zero.

### Primary model-held-out performance

The CardioQueue primary model reached AUROC 0.9841 (gene-clustered 95% CI, 0.9772 to 0.9890) and AUPRC 0.9881 (0.9769 to 0.9924). At the prespecified 0.5 threshold, sensitivity was 0.9794 (0.9672 to 0.9883) and specificity was 0.8407 (0.7704 to 0.8878). Source-specific AUROC ranged from 0.9632 for CardioBoost records to 0.9947 for HiRO (**Table 1**).

Specificity was lower than in the earlier 430-variant analysis because the 775-allele benchmark added externally adjudicated SHaRE and PLOS variants and changed the benign case mix. The model produced 54 false-positive scores among 339 B/LB alleles, including 24 in SHaRE and 14 in CardioBoost. AUROC remained high, indicating that this was chiefly a fixed-threshold transport and calibration issue rather than loss of overall ranking. The score should therefore be locally calibrated if a laboratory intends to use a fixed operational threshold.

In the post hoc 64-allele multi-source consensus subset, AUROC was 0.9931 (95% CI, 0.9531 to 1.000), sensitivity was 1.000, and specificity was 0.8929. This result supports ranking consistency when source classifications agree, but its size and shared evidence ecosystem preclude interpreting it as independent clinical validation. A broader 1,233-allele quarantine analysis, same-site/protein-substitution quarantine, removal of precomputed effect scores, and restriction to higher-review-status ClinVar development labels produced similar discrimination (Additional file 1).

### Review zones and ranked worklists

Using the exploratory 0.1 and 0.9 boundaries, 405 of 436 P/LP alleles entered accelerated review, 29 were deferred, and 2 entered lower-priority review. Among 339 B/LB alleles, 19 entered accelerated review, 115 were deferred, and 205 entered lower-priority review. The nondeferred subset covered 81.4% of the benchmark and had 96.7% accuracy; these estimates describe reference-labelled variants and do not establish clinical safety (**Figure 3**).

The primary model assigned 17,380 of 42,361 registry VUS to accelerated review, 15,107 to deferred review, and 9,874 to lower-priority review. The public research worklist excluded all HiRO-linked records and contained 41,622 VUS: 17,258 accelerated, 14,727 deferred, and 9,637 lower-priority. These zones indicate review order only. They do not revise a submitted classification, diagnose a patient, or justify management.

The fixed legacy 430-variant cohort was retained for review-budget analyses because its predictions and comparator rankings were frozen together. The first 86 reviews contained 86 P/LP variants. Recovering 80% of the 255 P/LP variants required 207 reviews, compared with a random-order median of 344 and an oracle minimum of 204. Reviewing half the cohort recovered 211 P/LP variants (82.7%) with 98.1% precision (**Figure 4A**). These values quantify enrichment within a constructed retrospective cohort, not time saved in practice.

### CardioBoost comparison

The released CardioBoost models reproduced AUROC 0.907 in their cardiomyopathy holdout and 0.948 in their arrhythmia holdout before comparison. The strict shared benchmark then retained 235 cDNA- and protein-concordant observations representing 211 unique variants. CardioQueue reached AUROC 0.984 versus 0.864 for CardioBoost, a paired difference of 0.120 (95% CI, 0.068 to 0.176). AUPRC was 0.991 versus 0.918, a difference of 0.074 (0.037 to 0.115) (**Figure 2**). This comparison is limited to harmonizable public observations and does not reproduce CardioBoost development or its patient-outcome analysis.

### Temporal VUS analysis

The temporal CardioQueue model was fitted using 23,851 training and 4,210 validation variants already classified B/LB or P/LP in the January 2024 ClinVar archive. No January 2024 VUS entered fitting. It then scored 25,481 historical VUS. By the July 2026 comparison snapshot, 146 were classified P/LP, 309 B/LB, 24,879 remained VUS, and 147 lacked an unambiguous current class.

Among all 25,481 historical VUS, the first 500 ranked reviews contained 16 of the 146 observed later-P/LP variants and the first 5,000 contained 93. Recovering 80% required 6,644 reviews, compared with a random-order median of 20,313 (**Figure 4B**). The 455 variants with an observed later binary class yielded AUROC 0.997, but this subset is selectively observed and should not be treated as a representative prospective cohort. Predictor resources were also not uniformly frozen to January 2024. The analysis therefore shows retrospective enrichment of later classifications, not prospective prediction of reclassification.

### Structure and FoldX analyses

Registry-wide, reference-checked mapping produced AlphaFold, DSSP, and FreeSASA features for 38,021 variants. A structure-enhanced sensitivity model changed AUROC from 0.977 to 0.979 in the legacy held-out cohort; the gene-clustered difference was 0.0025 (95% CI, -0.00003 to 0.00502). It did not clearly improve transfer to 19 genes absent from training.

FoldX returned numeric delta-delta-G estimates for all 4,137 eligible structure-mutation calculations. The primary 775-allele benchmark included 382 variants with FoldX results at pLDDT >=70. At this confidence threshold, adding signed and absolute delta-delta-G to structure eligibility changed AUROC by -0.00233 (95% CI, -0.00390 to -0.00109). A small improvement over the baseline at pLDDT >=50 did not clearly exceed the eligibility-only control. Thus, the study demonstrates scalable multi-gene cardiogenetic stability annotation, but FoldX did not add independent predictive value to the primary model.

## Discussion

CardioQueue separated P/LP from B/LB variants across 775 model-held-out alleles from five sources and concentrated known P/LP variants near the front of finite review queues. Its most relevant contribution is operational: it translates heterogeneous genomic, sequence, and protein annotations into a reproducible order for expert evidence review. It is not an ACMG/AMP classifier, and discrimination alone does not establish clinical utility.

The model is broader than a missense-only predictor. It covers multiple consequence classes across 40 evaluated cardiogenetic genes and integrates population, conservation, splice, established effect, protein-sequence, and structure-derived information. The direct CardioBoost comparison suggests improved ranking on the harmonizable shared subset, while CardioBoost retains patient-outcome evidence that CardioQueue does not have ([Zhang et al., 2021](https://doi.org/10.1038/s41436-020-00972-3)). CardioClassifier and CardioVAI support ACMG/AMP evidence assembly, GENESIS focuses on selected channelopathy genes, and CVD-PP predicts cardiovascular VUS pathogenicity. CardioQueue addresses a different gap: deciding which records should receive scarce expert review first.

The five-source benchmark strengthens that claim but requires careful wording. HiRO is external to model development yet internal to the contributing institution. eMERGE, CardioBoost, SHaRE, and PLOS are institutionally external. Exact held-out alleles did not enter fitting or early stopping, and previously used SHaRE/PLOS alleles were excluded. However, classifications and predictor resources can share ClinVar submissions, literature, laboratories, genes, and domains. The evaluation is therefore model-held-out and multi-source, not fully independent of the public evidence ecosystem.

For a clinical genetics service, the immediate use is scheduled backlog management. CardioQueue could rank aging VUS for literature refresh, transcript reconciliation, segregation follow-up, functional-assay selection, or multidisciplinary review. High-scoring B/LB and low-scoring P/LP records may also identify cases needing annotation or evidence reconciliation. The practical change for clinicians would be a shorter, versioned list of molecular records that deserve deeper review, while phenotype, family history, inheritance, disease mechanism, and management urgency continue to determine patient care.

Several safeguards are necessary. The score must rank variants, not patients. Urgent cases should bypass the queue. Lower-priority variants require a maximum review interval and periodic unranked safety review. No classification or management should change until a qualified laboratory completes full reinterpretation and issues an amended report. Local calibration is also important: specificity fell when the benchmark expanded from 430 to 775 variants, even though AUROC remained stable, showing that a fixed threshold does not transport automatically across source mixtures.

The historical VUS analysis is encouraging because the model was fitted without any baseline VUS and later-P/LP variants were enriched toward the front of the complete 25,481-record queue. It still does not prove that CardioQueue prospectively predicts reclassification. Only a minority of baseline VUS acquired a binary outcome, reclassification is driven by nonrandom scrutiny and evidence generation, and several predictor versions postdated baseline. A prospective study should lock the model and annotations, randomize or use a stepped-wedge design for the order of otherwise eligible VUS reviews, and compare time to completed review, evidence yield, classification change, reviewer time, and safety overrides.

Protein structure is a secondary mechanistic contribution. AlphaFold geometry could be computed at scale, and FoldX produced delta-delta-G for 4,137 cardiogenetic substitutions. Few cardiogenetic studies have tested stability estimates across this many genes with reference-residue checks, confidence tiers, and an eligibility-only control. The negative result matters: numeric delta-delta-G did not improve the high-confidence model beyond knowing that a variant mapped to a reliable structure. FoldX should therefore remain an interpretable annotation for selected missense review, not a claimed source of primary model performance.

Limitations include retrospective source selection, a development set derived from ClinVar, incomplete independence among evaluation sources, uncertain ancestry representation, sparse strata, and exclusion of structural variants and many noncoding mechanisms. The primary cohort's P/LP prevalence differs from a clinical VUS backlog, so precision and workload estimates will not transport directly. Calibration was imperfect, and no prospective laboratory workflow, reviewer-time endpoint, classification-change endpoint, or patient outcome was measured. The 40-gene benchmark supports broad cardiogenetic applicability within its observed scope, not universality across all inherited cardiac disease genes or variant types.

The next decisive test is prospective clinical workflow evaluation. Two cardiogenetic reviewers could compare a CardioQueue-ranked arm with random or usual-practice ordering, using blinded adjudication and prespecified review budgets. Such a study would establish whether prioritization increases actionable evidence yield without delaying important variants. It would also support local threshold selection and reveal whether the model complements clinician judgment or restates information already used in review.

## Conclusions

CardioQueue showed high discrimination across 775 model-held-out cardiogenetic variants and enriched known or later-observed P/LP variants near the front of retrospective review queues. Its clinical role is to prioritize records for expert reinterpretation, not to classify variants or direct patient care. Prospective evaluation must now test whether a locked, locally calibrated queue improves review yield, timeliness, and consistency without delaying clinically important variants.

## List of abbreviations

ACMG: American College of Medical Genetics and Genomics  
AMP: Association for Molecular Pathology  
AUPRC: area under the precision-recall curve  
AUROC: area under the receiver operating characteristic curve  
B/LB: benign/likely benign  
CI: confidence interval  
GRCh: Genome Reference Consortium human assembly  
HGVS: Human Genome Variation Society  
NPV: negative predictive value  
P/LP: pathogenic/likely pathogenic  
PPV: positive predictive value  
VUS: variant of uncertain significance  
WES: whole-exome sequencing

## Declarations

### Ethics approval and consent to participate

[REQUIRED BEFORE SUBMISSION: Insert the research ethics committee name, protocol/reference number, approval date, and consent or waiver basis for HiRO/CASPER WES/VERDICT records. State the data-governance approval and confirm compliance with the Declaration of Helsinki.]

### Consent for publication

[REQUIRED BEFORE SUBMISSION: State whether identifiable individual data are present. If none, use: "Not applicable. The manuscript contains no identifiable individual-level information."]

### Availability of data and materials

Public resources are available from their cited sources under their access terms. At publication, a public GitHub repository and an archived release will contain the trained model, inference code, frozen feature schema, a research-only public-source VUS worklist, cleaned non-identifying tables permitted for redistribution, analysis scripts, manifests, environment, and checksums. The worklist excludes all HiRO-linked rows and does not provide clinical classifications. SHaRe-derived row-level data will be distributed only if written permission permits; otherwise, the repository will provide acquisition instructions, checksums, and scripts that authorized users can run locally. Insert the final repository URL, archive DOI, commit, code license, and model license before submission. HiRO/CASPER WES/VERDICT records contain participant-linked clinical information and are not publicly available. Access conditions are [INSERT CONTROLLED-ACCESS PROCESS, GOVERNANCE CONTACT, AND ELIGIBILITY CONDITIONS].

### Competing interests

[REQUIRED: Declare all interests using author initials, or state: "The authors declare that they have no competing interests."]

### Funding

[REQUIRED: List every funding source and grant number, and describe the funders' roles. If there was no external funding, state that explicitly.]

### Authors' contributions

[REQUIRED: Provide author initials and contributions for conceptualization, data curation, methodology, software, analysis, investigation, visualization, supervision, funding, writing, and final approval.]

### Acknowledgements

[REQUIRED: Name contributors who do not meet authorship criteria and confirm permission, or state "Not applicable."]

### Authors' information

Not applicable.

## Additional files

**Additional file 1 (.pdf): Supplementary methods and figures.** Cohort construction, predictor groups, sensitivity analyses, calibration, structure and FoldX analyses, review-queue methods, temporal VUS analysis, and supplementary figure legends.

**Additional file 2 (.xlsx): Supplementary tables.** Source accounting, primary and source-specific performance, CardioBoost comparison, robustness analyses, structure and FoldX results, review-queue outcomes, temporal VUS results, and gene-panel quality control.

## Table

**Table 1. Performance of the CardioQueue primary model in model-held-out alleles.** Pooled performance counts each exact concordantly labelled allele once. Source rows retain one assertion per source and therefore overlap. HiRO is institutionally internal; all other named sources are institutionally external.

| Cohort | n (P/LP; B/LB) | AUROC | AUPRC | Sensitivity | Specificity |
|---|---:|---:|---:|---:|---:|
| Combined unique alleles | 775 (436; 339) | 0.984 | 0.988 | 0.979 | 0.841 |
| HiRO | 69 (27; 42) | 0.995 | 0.993 | 1.000 | 0.881 |
| eMERGE | 176 (96; 80) | 0.985 | 0.987 | 0.979 | 0.875 |
| CardioBoost | 173 (129; 44) | 0.963 | 0.986 | 0.984 | 0.682 |
| SHaRE | 337 (198; 139) | 0.986 | 0.992 | 0.985 | 0.827 |
| Fernandez-Falgueras/PLOS | 88 (24; 64) | 0.977 | 0.964 | 0.917 | 0.938 |

## Figure legends

**Figure 1. CardioQueue study workflow and interpretation boundary.** The numbered path separates registry assembly, exact-allele quarantine, model fitting, evaluation in 775 held-out alleles, robustness analyses, and retrospective translation to review order. The lower band shows the required sequence from model score to expert ACMG/AMP adjudication and, only when independently supported, amended classification or management. Urgent or management-sensitive cases must bypass score order.

**Figure 2. Matched comparison with public CardioBoost.** Receiver operating characteristic and precision-recall curves compare 235 cDNA- and protein-concordant observations, representing 211 unique variants, scored by both models. This subset does not recreate the CardioBoost development cohort or its patient-outcome analysis.

**Figure 3. Three-zone review analysis.** Panel A shows primary-model score distributions for 339 B/LB and 436 P/LP model-held-out alleles. Panel B gives zone counts at the exploratory boundaries of <=0.1, >0.1 to <0.9, and >=0.9. Panel C applies the same zones to all 42,361 registry VUS and overlapping source subsets, which must not be summed. VUS outcomes are unknown. The zones set review order only and cannot independently change classification or care.

**Figure 4. Cumulative recovery under score-ranked review.** Panel A shows recovery of 255 P/LP variants in the fixed 430-variant legacy queue; the band is a 95% bootstrap interval and the gray reference summarizes 10,000 random orderings. Panel B ranks all 25,481 eligible January 2024 ClinVar VUS and shows recovery of the 146 variants observed as P/LP by July 2026. Still-VUS records are outcome-unobserved. Predictor resources were not uniformly frozen to baseline, so Panel B is a retrospective label-temporal analysis, not prospective validation.

## References

Adam F, Fluri M, Scherz A, Rabaglio M. Occurrence of variants of unknown clinical significance in genetic testing for hereditary breast and ovarian cancer syndrome and Lynch syndrome: a literature review and analytical observational retrospective cohort study. *BMC Med Genomics*. 2023;16:7. doi:10.1186/s12920-023-01437-7.

Al-Mahrami N, Albalushi A, Al Hattali F, et al. CardioVar: a machine learning framework for pathogenicity prediction of cardiomyopathy genetic variants. *Bioinformatics Advances*. 2026;6:vbag135. doi:10.1093/bioadv/vbag135.

Alirezaie N, Kernohan KD, Hartley T, Majewski J, Hocking TD. ClinPred: prediction tool to identify disease-relevant nonsynonymous single-nucleotide variants. *Am J Hum Genet*. 2018;103:474-483. doi:10.1016/j.ajhg.2018.08.005.

Anderson D, Lassmann T. An expanded phenotype centric benchmark of variant prioritisation tools. *Hum Mutat*. 2022;43:539-546. doi:10.1002/humu.24362.

Arbustini E, Behr ER, Carrier L, et al. Interpretation and actionability of genetic variants in cardiomyopathies: a position statement from the European Society of Cardiology Council on cardiovascular genomics. *Eur Heart J*. 2022;43:1901-1916. doi:10.1093/eurheartj/ehab895.

Cannon S, Williams M, Gunning AC, Wright CF. Evaluation of in silico pathogenicity prediction tools for the classification of small in-frame indels. *BMC Med Genomics*. 2023;16:36. doi:10.1186/s12920-023-01454-6.

Chao K, Wilson M, Goodrich J, gnomAD Production Team. gnomAD v4.1. gnomAD browser. 2024. Available from: https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1. Accessed 16 September 2026.

Cheng J, Novati G, Pan J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science*. 2023;381:eadg7492. doi:10.1126/science.adg7492.

Collins GS, Moons KGM, Dhiman P, et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. *BMJ*. 2024;385:e078378. doi:10.1136/bmj-2023-078378.

Dahary D, Golan Y, Mazor Y, et al. Genome analysis and knowledge-driven variant interpretation with TGex. *BMC Med Genomics*. 2019;12:200. doi:10.1186/s12920-019-0647-8.

Draelos RL, et al. GENESIS: gene-specific machine learning models for variants of uncertain significance found in catecholaminergic polymorphic ventricular tachycardia and long QT syndrome-associated genes. *Circ Arrhythm Electrophysiol*. 2022;15:e010326. doi:10.1161/CIRCEP.121.010326.

Delgado J, Radusky LG, Cianferoni D, Serrano L. FoldX 5.0: working with RNA, small molecules and a new graphical interface. *Bioinformatics*. 2019;35:4168-4169. doi:10.1093/bioinformatics/btz184.

Frazer J, Notin P, Dias M, et al. Disease variant prediction with deep generative models of evolutionary data. *Nature*. 2021;599:91-95. doi:10.1038/s41586-021-04043-8.

dbNSFP. Releases: dbNSFP v5.3.1. 2026. Available from: https://www.dbnsfp.org/releases/. Accessed 16 September 2026.

EMBL-EBI. AlphaFold Database version 6 release notes. 2025. Available from: https://www.ebi.ac.uk/pdbe/news/alphafold-database-release-notes. Accessed 16 September 2026.

Fernandez-Falgueras A, et al. The importance of variant reinterpretation in inherited cardiovascular diseases: establishing the optimal timeframe. *PLoS One*. 2024;19:e0297914. doi:10.1371/journal.pone.0297914.

Glazer AM, Davogustto G, Shaffer CM, et al. Arrhythmia variant associations and reclassifications in the eMERGE-III sequencing study. *Circulation*. 2022;145:877-891. doi:10.1161/CIRCULATIONAHA.121.055562.

Ioannidis NM, Rothstein JH, Pejaver V, et al. REVEL: an ensemble method for predicting the pathogenicity of rare missense variants. *Am J Hum Genet*. 2016;99:877-885. doi:10.1016/j.ajhg.2016.08.016.

Jagadeesh KA, Wenger AM, Berger MJ, et al. M-CAP eliminates a majority of variants of uncertain significance in clinical exomes at high sensitivity. *Nat Genet*. 2016;48:1581-1586. doi:10.1038/ng.3703.

Jaganathan K, Kyriazopoulou Panagiotopoulou S, McRae JF, et al. Predicting splicing from primary sequence with deep learning. *Cell*. 2019;176:535-548.e24. doi:10.1016/j.cell.2018.12.015.

Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction with AlphaFold. *Nature*. 2021;596:583-589. doi:10.1038/s41586-021-03819-2.

Li Q, Wang K. InterVar: clinical interpretation of genetic variants by the 2015 ACMG-AMP guidelines. *Am J Hum Genet*. 2017;100:267-280. doi:10.1016/j.ajhg.2017.01.004.

Kabsch W, Sander C. Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. *Biopolymers*. 1983;22:2577-2637. doi:10.1002/bip.360221211.

Karczewski KJ, Francioli LC, Tiao G, et al. The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature*. 2020;581:434-443. doi:10.1038/s41586-020-2308-7.

Kircher M, Witten DM, Jain P, O'Roak BJ, Cooper GM, Shendure J. A general framework for estimating the relative pathogenicity of human genetic variants. *Nat Genet*. 2014;46:310-315. doi:10.1038/ng.2892.

Kroncke BM, Mendenhall J, Smith DK, et al. Protein structure aids predicting functional perturbation of missense variants in SCN5A and KCNQ1. *Comput Struct Biotechnol J*. 2019;17:206-214. doi:10.1016/j.csbj.2019.01.008.

Kwong A, Ho CYS, Shin VY, Au CH, Chan TL, Ma ESK. How does re-classification of variants of unknown significance (VUS) impact the management of patients at risk for hereditary breast cancer? *BMC Med Genomics*. 2022;15:122. doi:10.1186/s12920-022-01270-4.

Landrum MJ, Lee JM, Benson M, et al. ClinVar: improving access to variant interpretations and supporting evidence. *Nucleic Acids Res*. 2018;46:D1062-D1067. doi:10.1093/nar/gkx1153.

Landstrom AP, Chahal AA, Ackerman MJ, et al. Interpreting incidentally identified variants in genes associated with heritable cardiovascular disease: a scientific statement from the American Heart Association. *Circ Genom Precis Med*. 2023;16:e000092. doi:10.1161/HCG.0000000000000092.

Liu X, Li C, Mou C, Dong Y, Tu Y. dbNSFP v4: a comprehensive database of transcript-specific functional predictions and annotations for human nonsynonymous and splice-site SNVs. *Genome Med*. 2020;12:103. doi:10.1186/s13073-020-00803-9.

Manshaei R, DeLong S, Andric V, et al. GeneTerpret: a customizable multilayer approach to genomic variant prioritization and interpretation. *BMC Med Genomics*. 2022;15:31. doi:10.1186/s12920-022-01166-3.

McLaren W, Gil L, Hunt SE, et al. The Ensembl Variant Effect Predictor. *Genome Biol*. 2016;17:122. doi:10.1186/s13059-016-0974-4.

Mitternacht S. FreeSASA: an open source C library for solvent accessible surface area calculations. *F1000Res*. 2016;5:189. doi:10.12688/f1000research.7931.1.

Nicora G, Limongelli I, Gambelli P, et al. CardioVAI: an automatic implementation of ACMG-AMP variant interpretation guidelines in the diagnosis of cardiovascular diseases. *Hum Mutat*. 2018;39:1835-1846. doi:10.1002/humu.23665.

Prokhorenkova L, Gusev G, Vorobev A, Dorogush AV, Gulin A. CatBoost: unbiased boosting with categorical features. In: *Advances in Neural Information Processing Systems 31*. 2018:6638-6648.

Qi H, Zhang H, Zhao Y, et al. MVP predicts the pathogenicity of missense variants by deep learning. *Nat Commun*. 2021;12:374. doi:10.1038/s41467-020-20847-0.

Ramaker ME, Abdulrahim JW, Corey KM, et al. Cardiovascular Disease Pathogenicity Predictor (CVD-PP): a tissue-specific in silico tool for discriminating pathogenicity of variants of unknown significance in cardiovascular disease genes. *Circ Genom Precis Med*. 2024;17:e004464. doi:10.1161/CIRCGEN.123.004464.

Rehm HL, Berg JS, Brooks LD, et al. ClinGen, the Clinical Genome Resource. *N Engl J Med*. 2015;372:2235-2242. doi:10.1056/NEJMsr1406261.

Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. *Genet Med*. 2015;17:405-424. doi:10.1038/gim.2015.30.

Schmidt A, Roner S, Mai K, Klinkhammer H, Kircher M, Ludwig KU. Predicting the pathogenicity of missense variants using features derived from AlphaFold2. *Bioinformatics*. 2023;39:btad280. doi:10.1093/bioinformatics/btad280.

Schymkowitz J, Borg J, Stricher F, et al. The FoldX web server: an online force field. *Nucleic Acids Res*. 2005;33:W382-W388. doi:10.1093/nar/gki387.

Tian Y, Pesaran T, Chamberlin A, et al. REVEL and BayesDel outperform other in silico meta-predictors for clinical variant classification. *Sci Rep*. 2019;9:12752. doi:10.1038/s41598-019-49224-8.

Suay-Corredera C, Pricolo MR, Herrero-Galan E, et al. Protein haploinsufficiency drivers identify MYBPC3 variants that cause hypertrophic cardiomyopathy. *J Biol Chem*. 2021;297:100854. doi:10.1016/j.jbc.2021.100854.

UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids Res*. 2025;53:D609-D617. doi:10.1093/nar/gkae1010.

UniProt Consortium. UniProt release 2026_02. 2026. Available from: https://www.uniprot.org/help/2026-06-10-release. Accessed 16 September 2026.

Van Calster B, McLernon DJ, van Smeden M, Wynants L, Steyerberg EW. Calibration: the Achilles heel of predictive analytics. *BMC Med*. 2019;17:230. doi:10.1186/s12916-019-1466-7.

Varadi M, Anyango S, Deshpande M, et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. *Nucleic Acids Res*. 2022;50:D439-D444. doi:10.1093/nar/gkab1061.

Whiffin N, Walsh R, Govind R, et al. CardioClassifier: disease- and gene-specific computational decision support for clinical genome interpretation. *Genet Med*. 2018;20:1246-1254. doi:10.1038/gim.2017.258.

Wu Y, Li R, Sun S, Weile J, Roth FP. Improved pathogenicity prediction for rare human missense variants. *Am J Hum Genet*. 2021;108:1891-1906. doi:10.1016/j.ajhg.2021.08.012.

Zhao H, et al. SIGMA leverages protein structural information to predict the pathogenicity of missense variants. *Cell Rep Methods*. 2024;4:100687. doi:10.1016/j.crmeth.2023.100687.

Zhang X, Walsh R, Whiffin N, et al. Disease-specific variant pathogenicity prediction significantly improves variant interpretation in inherited cardiac conditions. *Genet Med*. 2021;23:69-79. doi:10.1038/s41436-020-00972-3.
