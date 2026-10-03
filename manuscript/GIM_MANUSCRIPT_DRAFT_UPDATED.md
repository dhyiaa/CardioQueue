# A Source-Held-Out CatBoost Model for Cardiogenetic Variant Prioritization With Explicit VUS Deferral

## Article Type

Original research article for *Genetics in Medicine*.

## Abstract

### Purpose

Cardiogenetic testing often identifies rare variants that remain difficult to classify. We built a cardiogenetics-specific machine-learning model for pathogenic versus benign variant prioritization with explicit deferral of uncertain calls.

### Methods

We integrated ClinVar, HiRO/CASPER WES/VERDICT, eMERGE-III, and public CardioBoost variants into a GRCh38 registry. The analysis-ready matrix contained 85,677 rows: 33,923 benign/likely benign, 9,148 pathogenic/likely pathogenic, and 42,361 VUS. The primary model was a binary CatBoost classifier trained on 42,990 benign or pathogenic rows. VUS rows were scored after training. We used source-held-out external validation, matched benchmarking against official CardioBoost, leakage ablation, calibration analysis, and stratified performance review.

### Results

On 430 source-held-out external binary rows, the model achieved AUROC 0.977, AUPRC 0.985, sensitivity 0.965, specificity 0.886, PPV 0.925, and Brier score 0.053. With 0.1/0.9 high-confidence thresholds, sensitivity was 0.937, specificity 0.943, PPV 0.960, NPV 0.977, and deferral 0.119. On a 250-row matched CardioBoost benchmark, our model outperformed official CardioBoost: AUROC 0.980 versus 0.843 and high-confidence accuracy 0.972 versus 0.847.

### Conclusion

A feature-rich source-held-out cardiogenetics CatBoost model improved high-confidence prioritization compared with public CardioBoost on matched rows while preserving a clinically explicit deferral zone.

## Introduction

Inherited cardiomyopathy and arrhythmia testing has made rare-variant discovery routine, but interpretation remains uneven. ACMG/AMP criteria provide the clinical framework, yet many variants lack segregation, functional, population, or disease-specific evidence [1]. ClinVar is valuable because it aggregates submitted classifications, review status, and disease assertions, but it also contains heterogeneous evidence levels, conflicting assertions, and disease-context ambiguity [2].

This problem is especially visible for missense variation. General deleteriousness tools such as CADD, REVEL, AlphaMissense, and ESM-based scores capture useful molecular signals, but they are not designed only for inherited cardiac disease [3-6]. CardioBoost addressed that gap by training disease-specific binary classifiers for rare missense variants in inherited cardiomyopathy and arrhythmia genes, then applying high-confidence thresholds: scores at least 0.9 were considered disease-causing, scores at most 0.1 were considered benign or likely benign, and intermediate scores were treated as indeterminate [7].

That design choice is important. CardioBoost did not train a true Benign/VUS/Pathogenic model. It trained a binary pathogenicity model and used a probability interval as a deferral zone. We followed the same clinical logic because VUS is not a stable biological class. A VUS means that current evidence is insufficient. Some VUS later become benign, some become pathogenic, and many remain unresolved.

Here we report a larger cardiogenetics-specific CatBoost model built from public and local variant sources. The model uses modern annotation resources, including gnomAD, dbNSFP, VEP, SpliceAI, AlphaMissense, ClinGen, UniProt, AlphaFold-derived protein features, and FoldX [3-6,8-15]. We tested the model under source-held-out validation, then directly compared it with the official public CardioBoost models on the same external matched rows.

The central question was simple: can a larger, feature-rich cardiogenetics model improve high-confidence variant prioritization while still deferring uncertain cases?

## Materials and Methods

### Study Design

We designed a retrospective variant-level modeling study with external source-held-out validation. The primary endpoint was binary classification of pathogenic/likely pathogenic versus benign/likely benign variants. VUS rows were excluded from supervised binary training and used for post-training triage.

Three analyses were prespecified for interpretation:

1. A primary binary CatBoost model with standard 0.5 threshold and CardioBoost-style 0.1/0.9 deferral thresholds.
2. A matched public CardioBoost benchmark using rows that both models could score.
3. Secondary analyses, including true three-class modeling, all-label VUS triage, leakage ablation, calibration, and stratified performance.

### Data Sources

ClinVar was the main public label source. We used the local ClinVar bulk file after germline filtering, conflict handling, and cardiogenetics gene-panel QC. Clinical significance was mapped into benign/likely benign, VUS, and pathogenic/likely pathogenic classes. Review stars and conflict status were retained for weighting and sensitivity analysis, but direct ClinVar confidence metadata was excluded from the primary predictor set.

HiRO/CASPER WES/VERDICT contributed 482 private patient-linked source records with genotype, phenotype, and ACMG-style adjudication fields. These records were not collapsed as meaningless duplicates. Repeated variants may represent distinct patient events, phenotypes, or source classifications. In the current modelable data, 408 HiRO source records linked to model predictions; 74 remain unresolved after coordinate and HGVS rescue.

eMERGE contributed 2,754 arrhythmia-gene source records from the eMERGE-III arrhythmia study [16]. An early poor eMERGE result was traced to coordinate build mismatch. We lifted the eMERGE coordinates from GRCh37 to GRCh38, checked reference alleles against the same FASTA used for VEP, and reran annotation. All 2,754 eMERGE rows passed single-mapping and GRCh38 reference checks.

CardioBoost contributed 355 public source records and official public model artifacts from the CardioBoost GitHub repository [7]. The local CardioBoost source is binary in the current processed evaluation set and has no VUS labels. It was used for external testing and public-model benchmarking.

### Gene Panel Review

We did not use a raw union of all genes found in the input files. Genes were reviewed for cardiac relevance and contamination risk. AKAP9, DMD, FPGT, FPGT-TNNI3K, KCNE1B, and KNCH2-like artifacts were excluded from the primary model. TTN was handled separately because truncating variants are important in dilated cardiomyopathy, while TTN missense variation is large, VUS-heavy, and mechanistically different from the rare missense scope of CardioBoost. Broader syndrome genes, sparse genes, and mixed-scope genes were retained only with explicit primary, expanded, monitor, sensitivity, or special-handling flags.

### Variant Registry and Label Handling

The analysis-ready table contained 85,677 rows and 420 columns. Rows were keyed by GRCh38 chromosome, position, reference allele, and alternate allele. This table was derived from a larger 86,889-row annotated intermediate matrix; all manuscript-facing results use the 85,677-row analysis-ready matrix.

The matrix contained 33,923 benign/likely benign rows, 42,361 VUS rows, 9,148 pathogenic/likely pathogenic rows, and 245 missing, conflicting, or non-supervised labels. The primary binary slice contained 42,990 rows: 33,859 benign and 9,131 pathogenic. The exploratory three-class slice contained 84,721 rows: 33,859 benign, 41,731 VUS, and 9,131 pathogenic. Full source and label counts are provided in Supplementary Sections S1 and S3.

### Feature Engineering

The feature matrix integrated variant identity, molecular consequence, splicing, population frequency, in silico pathogenicity, conservation, gene-disease validity, protein context, and source-record evidence features. Coverage was high for VEP consequence and HGVS fields (85,660 of 85,677 rows, 99.98%), SpliceAI (77,141 rows, 90.04%), and gnomAD final allele frequency (62,937 rows, 73.46%). dbNSFP covered 43,774 rows (51.09%), AlphaMissense direct joins covered 38,023 rows (44.38%), and FoldX DDG covered 313 rows (0.37%). Detailed feature coverage is provided in Supplementary Section S4.

dbNSFP contributed missense-focused in silico and conservation features, including SIFT, PolyPhen-2, REVEL, MetaLR, FATHMM-XF, CADD, GERP, PhyloP, and PhastCons [4]. VEP supplied consequence and HGVS fields [8]. SpliceAI supplied splice donor and acceptor gain/loss scores [9]. gnomAD-derived features included allele frequency, population maximum allele frequency, observed status, and homozygote count [10]. AlphaMissense was joined directly when possible and through dbNSFP when available [5]. ClinGen supplied gene-disease validity and limited variant-level evidence flags [11].

Protein features included UniProt domains and sites, AlphaFold residue confidence, DSSP secondary structure, FreeSASA solvent accessibility, and FoldX stability estimates for eligible variants [12-15]. These features were intentionally sparse because they require reliable protein-position mapping and, for FoldX, compatible structure and residue context.

### Model Training

The primary learner was CatBoost, chosen because it handles mixed feature types, nonlinear effects, missingness, and tabular biomedical data well [17]. The primary model used 338 selected non-leakage features after excluding labels, raw source classifications, patient identifiers, direct variant identifiers, and review metadata as predictors.

Sample weights combined label confidence and class imbalance. Pathogenic rows were upweighted because they were less common and because missing a pathogenic variant is clinically more serious. ClinVar one-star rows were not discarded from the primary model; they were retained with lower confidence contribution. Label conflicts and VUS rows were excluded from supervised binary training.

### Split Strategy

We used three split strategies: an internal grouped split for development, a source-held-out split for the main external validation, and a gene-stress split for sparse-gene evaluation. The paper-facing source-held-out split held out HiRO, eMERGE, and CardioBoost source strata for external validation. Rows were grouped by variant identity so the same GRCh38 chrom-pos-ref-alt variant could not appear in both training and test under different source names. The source-held-out binary split contained 36,176 training rows, 6,384 validation rows, and 430 external rows. Full split composition is provided in Supplementary Section S5.

### Matched CardioBoost Benchmark

The official public CardioBoost cardiomyopathy and arrhythmia AdaBoost models were run locally using the official preprocessing and prediction scripts [7]. Official predictions were matched to our model predictions by HGVS cDNA for the primary benchmark. Exact coordinate matching was used as QC but recovered too few rows for the main comparison, likely due to genome-build and transcript-version differences.

The matched benchmark was restricted to source-held-out rows with binary labels that both models could score. This created a 250-row external benchmark and a stricter 164-row CardioBoost-source-only benchmark. We used paired bootstrap resampling with 5,000 replicates to estimate confidence intervals for model differences.

### Leakage and Ablation Audit

To address possible shortcut learning, we retrained the primary binary model after removing: direct gene identity, VEP consequence and HGVS features, ClinVar metadata, source-record aggregates, HGVS/transcript identifiers, and all of these together in a strict low-leakage stress test. Each ablation used the same split and CatBoost recipe as the primary model.

## Results

### Source-Held-Out Binary Performance

At the standard 0.5 threshold, the primary model performed strongly across all source-held-out external sets.

**Table 1. Source-held-out external validation of the primary binary model**

| External set | Rows | AUROC | AUPRC | Sensitivity at 0.5 | Specificity at 0.5 | PPV at 0.5 | Sensitivity with deferral | Specificity with deferral | PPV with deferral | Deferral |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Combined external | 430 | 0.977 | 0.985 | 0.965 | 0.886 | 0.925 | 0.937 | 0.943 | 0.960 | 0.119 |
| HiRO | 69 | 0.988 | 0.983 | 0.963 | 0.881 | 0.839 | 0.926 | 0.929 | 0.893 | 0.087 |
| eMERGE | 176 | 0.991 | 0.991 | 0.979 | 0.938 | 0.949 | 0.948 | 0.988 | 0.989 | 0.114 |
| CardioBoost | 185 | 0.952 | 0.981 | 0.955 | 0.811 | 0.926 | 0.932 | 0.887 | 0.953 | 0.135 |

The external combined deferral confusion matrix and source-specific confusion matrices are provided in Supplementary Section S6.

### Same-Row Benchmark Against Official CardioBoost

The matched public-model benchmark directly compared our CatBoost predictions with official CardioBoost predictions on the same source-held-out rows. On 250 external matched rows, our model had higher AUROC, AUPRC, sensitivity, specificity, PPV, and high-confidence accuracy, with a lower deferral rate.

**Table 2. Same-row external benchmark against official CardioBoost**

| Metric | Our CatBoost | Official CardioBoost | Difference | 95% CI |
|---|---:|---:|---:|---:|
| AUROC | 0.980 | 0.843 | 0.138 | 0.083 to 0.195 |
| AUPRC | 0.990 | 0.902 | 0.087 | 0.048 to 0.132 |
| Sensitivity | 0.959 | 0.876 | 0.083 | 0.030 to 0.136 |
| Specificity | 0.543 | 0.148 | 0.395 | 0.272 to 0.519 |
| PPV | 0.964 | 0.836 | 0.128 | 0.082 to 0.172 |
| Deferral | 0.152 | 0.244 | -0.092 | -0.156 to -0.028 |
| High-confidence accuracy | 0.972 | 0.847 | 0.125 | 0.082 to 0.167 |

The stricter 164-row CardioBoost-source-only benchmark gave the same broad conclusion: AUROC 0.972 for our model versus 0.788 for official CardioBoost, AUPRC 0.990 versus 0.905, sensitivity 0.967 versus 0.870, and high-confidence accuracy 0.957 versus 0.846.

### Calibration and Stratified Performance

The validation Brier score was 0.003. On external data, Brier scores were 0.053 combined, 0.068 for HiRO, 0.030 for eMERGE, and 0.069 for CardioBoost. Mean predicted pathogenic probability was mildly higher than observed pathogenic fraction on external data, suggesting mild overconfidence. The model is therefore best framed as a high-confidence prioritization tool, not as a fully calibrated absolute risk model.

Stratified analyses showed strong performance in arrhythmia genes, cardiomyopathy genes, missense variants, dbNSFP-matched variants, and gnomAD-observed variants. The gnomAD-confirmed-absent stratum was small and weaker, so it should be interpreted as a cautionary subgroup. Detailed calibration and stratified performance tables are provided in Supplementary Section S11.

### All-Label VUS Triage

We applied the binary model to benign, VUS, and pathogenic rows using the same three zones: benign-like, deferred, and pathogenic-like. This analysis is not a true three-class classifier. It asks which VUS resemble benign or pathogenic training examples.

Across 85,677 model rows, exact three-zone agreement was 0.670, P/LP sensitivity was 0.981, and VUS deferral was 0.178. In HiRO source records with predictions, no pathogenic source record was called benign-like: 41 of 46 pathogenic rows were called pathogenic-like and 5 were deferred.

The high number of VUS-to-pathogenic-like calls in ClinVar and eMERGE should not be treated as errors by itself. It identifies variants whose current uncertain label conflicts with molecular evidence learned from binary P/LP versus B/LB examples.
Full all-label triage tables, source-record counts, and hard-error counts are provided in Supplementary Section S7.

### True Three-Class Model

The true Benign/VUS/Pathogenic CatBoost model was trained as a secondary experiment. It improved exact VUS label matching in HiRO, but had lower pathogenic sensitivity than the binary deferral model.

On HiRO source records with predictions, the true three-class model achieved exact accuracy 0.7255, pathogenic recall 0.6522, and VUS recall 0.8895. It made no pathogenic-to-benign hard errors and no benign-to-pathogenic escalations. However, 16 pathogenic HiRO source records were sent to VUS. This supports keeping the binary deferral model as the primary clinical triage model and treating the true three-class model as an exploratory evidence-status model.

### Leakage and Ablation Audit

The model remained strong after removal of reviewer-sensitive feature families. External-all AUROC and AUPRC were 0.977 and 0.984 for the baseline. Removing gene identity gave AUROC 0.977 and AUPRC 0.984. Removing source-record aggregates gave AUROC 0.977 and AUPRC 0.983. Removing VEP consequence and HGVS-style annotations gave AUROC 0.973 and AUPRC 0.980. A strict low-leakage model that removed gene identity, consequence/HGVS/transcript identifiers, ClinVar metadata, and source-record aggregates still achieved AUROC 0.974 and AUPRC 0.981, although specificity and PPV decreased.

These results do not prove absence of all leakage, but they argue against one obvious shortcut driving the external results.

## Discussion

This study supports a focused claim. A cardiogenetics-specific CatBoost model trained on a larger integrated dataset, with modern annotation features and source-held-out validation, improved high-confidence variant prioritization compared with official public CardioBoost on matched external rows.

The strongest evidence is not the internal split. Internal performance was near perfect, which is expected in a large feature-rich variant table and is not persuasive by itself. The meaningful evidence is the source-held-out validation and the matched public-model benchmark. On the 430-row external binary set, AUROC was 0.977 and AUPRC was 0.985. On the 250-row matched CardioBoost benchmark, our model improved AUROC by 0.138 and high-confidence accuracy by 0.125, with paired bootstrap confidence intervals that did not cross zero.

The model also preserved a clinically useful deferral zone. This matters because forcing every variant into pathogenic or benign categories is not how clinical interpretation works. The 0.1/0.9 thresholds produced a combined external deferral rate of 0.119, PPV 0.960, and NPV 0.977. In HiRO source records, there were no pathogenic-to-benign hard errors under the all-label deferral analysis.

The comparison to CardioBoost should still be written with care. We did not reproduce the private CardioBoost training cohort or original development environment. We ran the public CardioBoost model artifacts locally and compared both models on the same matched external rows. That is a fair public-model comparison, not a full replication of the original CardioBoost study.

The three-class model tells a different story. It better matches VUS labels but misses more pathogenic cases by sending them to VUS. Because VUS is an evidence state, we consider the binary model with deferral the main clinical model. The three-class model remains useful for secondary analysis and for studying uncertainty patterns.

Several limitations remain. First, ClinVar labels are submitted clinical assertions and can encode ascertainment, submitter behavior, and evidence heterogeneity [2]. Second, HiRO is patient-linked and valuable but small after binary variant-level filtering. Third, CardioBoost matching by HGVS cDNA was necessary because exact coordinate matching recovered too few rows. Fourth, calibration was mildly overconfident on external data, so raw probabilities should not be presented as absolute risk.

In summary, this model is mature enough for a manuscript centered on binary cardiogenetic pathogenicity prioritization with explicit deferral, external source-held-out validation, and direct public-model benchmarking.

## Figure Placeholders

The main manuscript uses five display items, consistent with the GIM original research limit: three figures and two tables. Detailed tables are moved to the Supplementary Information.

**Figure 1. Study design and registry construction.** Flow from ClinVar, HiRO, eMERGE, and CardioBoost to the GRCh38 registry, analysis-ready matrix, binary training slice, VUS scoring pool, feature families, and source-held-out validation.

**Figure 2. Primary model performance and calibration.** ROC and precision-recall curves for validation, HiRO, eMERGE, CardioBoost, and combined external validation, with a calibration reliability panel.

**Figure 3. Matched CardioBoost benchmark and clinical triage.** ROC and precision-recall curves comparing our model with official CardioBoost on the same external matched rows, paired bootstrap CI summary, and three-zone VUS triage distribution.

## Data and Code Availability

ClinVar, gnomAD, dbNSFP, VEP, SpliceAI, AlphaMissense, UniProt, ClinGen, AlphaFold, and CardioBoost public resources are available from their respective repositories subject to their licenses. HiRO/CASPER WES/VERDICT source records contain local patient-linked information and are not public in raw form. The analysis code, derived non-identifying tables, model reports, and figure-generation scripts are stored in the project repository. Final public release should exclude private HiRO identifiers and any restricted-license raw files.

## Ethics Statement

HiRO/CASPER WES/VERDICT records should be reported under the applicable local ethics and data-governance approvals. Public resources were used according to their stated licenses and access terms.

## References

1. Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants. *Genet Med*. 2015;17:405-424.
2. Landrum MJ, Lee JM, Benson M, et al. ClinVar: improving access to variant interpretations and supporting evidence. *Nucleic Acids Res*. 2018;46:D1062-D1067.
3. Ioannidis NM, Rothstein JH, Pejaver V, et al. REVEL: an ensemble method for predicting the pathogenicity of rare missense variants. *Am J Hum Genet*. 2016;99:877-885.
4. Liu X, Li C, Mou C, Dong Y, Tu Y. dbNSFP v4: a comprehensive database of transcript-specific functional predictions and annotations. *Genome Med*. 2020;12:103.
5. Cheng J, Novati G, Pan J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science*. 2023;381:eadg7492.
6. Meier J, Rao R, Verkuil R, et al. Language models enable zero-shot prediction of the effects of mutations on protein function. *NeurIPS*. 2021.
7. Zhang S, Walsh R, Whiffin N, et al. Disease-specific variant pathogenicity prediction significantly improves variant interpretation in inherited cardiac conditions. *Genet Med*. 2021.
8. McLaren W, Gil L, Hunt SE, et al. The Ensembl Variant Effect Predictor. *Genome Biol*. 2016;17:122.
9. Jaganathan K, Kyriazopoulou Panagiotopoulou S, McRae JF, et al. Predicting splicing from primary sequence with deep learning. *Cell*. 2019;176:535-548.
10. Karczewski KJ, Francioli LC, Tiao G, et al. The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature*. 2020;581:434-443.
11. Rehm HL, Berg JS, Brooks LD, et al. ClinGen: the Clinical Genome Resource. *N Engl J Med*. 2015;372:2235-2242.
12. UniProt Consortium. UniProt: the Universal Protein Knowledgebase. *Nucleic Acids Res*. 2025.
13. Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction with AlphaFold. *Nature*. 2021;596:583-589.
14. Varadi M, Anyango S, Deshpande M, et al. AlphaFold Protein Structure Database. *Nucleic Acids Res*. 2022;50:D439-D444.
15. Schymkowitz J, Borg J, Stricher F, et al. The FoldX web server: an online force field. *Nucleic Acids Res*. 2005;33:W382-W388.
16. Glazer AM, Wada Y, Li B, et al. Arrhythmia variant associations and reclassifications in the eMERGE-III sequencing study. *Circulation*. 2022.
17. Prokhorenkova L, Gusev G, Vorobev A, Dorogush AV, Gulin A. CatBoost: unbiased boosting with categorical features. *NeurIPS*. 2018.
18. Kircher M, Witten DM, Jain P, et al. A general framework for estimating the relative pathogenicity of human genetic variants. *Nat Genet*. 2014;46:310-315.
19. Freeman PJ, Hart RK, Gretton LJ, et al. VariantValidator: accurate validation, mapping, and formatting of sequence variation descriptions. *Hum Mutat*. 2018;39:61-68.
