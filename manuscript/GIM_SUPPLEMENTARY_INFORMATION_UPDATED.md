# Supplementary Information

# A Source-Held-Out CatBoost Model for Cardiogenetic Variant Prioritization With Explicit VUS Deferral

## Supplementary Contents

S1. Data source integration and registry construction  
S2. Gene-panel quality control  
S3. HiRO, eMERGE, and CardioBoost source-record handling  
S4. Feature engineering and coverage  
S5. Model training, weighting, and split assignment  
S6. Primary binary model performance  
S7. All-label VUS triage  
S8. True three-class model  
S9. Matched public CardioBoost benchmark  
S10. Leakage and ablation audit  
S11. Calibration and stratified performance  
S12. Reproducibility files and scripts  

## How This Supplement Links to the Main Manuscript

The main manuscript is kept within the GIM original research display-item limit by using three main figures and two main tables. Detailed source counts, feature coverage, split composition, secondary models, calibration, stratified performance, and benchmarking details are placed here.

| Main manuscript item | What it shows | Supporting SI sections |
|---|---|---|
| Figure 1 | Registry construction, source flow, feature families, and source-held-out design | S1-S5 |
| Figure 2 | Primary model ROC/PR performance and calibration | S6 and S11 |
| Figure 3 | Matched CardioBoost benchmark and clinical triage | S7 and S9 |
| Table 1 | Source-held-out external validation of the primary binary model | S5-S6 |
| Table 2 | Same-row benchmark against official CardioBoost with paired bootstrap CIs | S9 |

Suggested supplementary figures:

| Supplementary figure | Content | Source artifact |
|---|---|---|
| Figure S1 | Source records versus matrix rows | `results/data_characteristics/figures/source_records_vs_matrix_rows.png` |
| Figure S2 | Feature-group coverage | `results/manuscript_figures_tables/figures/feature_group_coverage_manuscript.png` |
| Figure S3 | Calibration reliability plot | `results/manuscript_figures_tables/figures/primary_binary_calibration_reliability.png` |
| Figure S4 | Stratified external performance heatmap | `results/manuscript_figures_tables/figures/external_all_stratified_performance_heatmap.png` |
| Figure S5 | VUS three-zone triage distribution | `results/manuscript_figures_tables/figures/vus_three_zone_triage_distribution.png` |

Suggested supplementary tables:

| Supplementary table | Content | SI section |
|---|---|---|
| Table S1 | Matrix row-count reconciliation and label distribution | S1 |
| Table S2 | Gene-panel QC and excluded/sensitivity genes | S2 |
| Table S3 | HiRO, eMERGE, and CardioBoost source-record handling | S3 |
| Table S4 | Feature family coverage | S4 |
| Table S5 | Split composition and weighting rules | S5 |
| Table S6 | Full primary binary model performance and confusion matrices | S6 |
| Table S7 | All-label VUS triage and hard-error counts | S7 |
| Table S8 | True three-class model results | S8 |
| Table S9 | CardioBoost matched benchmark details and comparators | S9 |
| Table S10 | Leakage and ablation audit | S10 |
| Table S11 | Calibration and stratified performance | S11 |

## S1. Data Source Integration and Registry Construction

The project used a GRCh38-centered registry keyed by chromosome, position, reference allele, and alternate allele. The registry was built to keep one row per resolved variant while preserving source-record tables for datasets where repeated variants may carry different source evidence.

Two row counts appear in project outputs:

| Matrix | Rows | Columns | Use |
|---|---:|---:|---|
| Annotated intermediate feature matrix | 86,889 | 463 | GRCh38 registry after feature annotation, before final ready-table cleanup |
| Analysis-ready modeling matrix | 85,677 | 420 | Final matrix used for splits, weights, model training, validation, calibration, and all-label triage |

All manuscript-facing model results use the 85,677-row analysis-ready matrix.

The current analysis-ready label distribution is:

| Label | Rows |
|---|---:|
| Benign / likely benign | 33,923 |
| VUS | 42,361 |
| Pathogenic / likely pathogenic | 9,148 |
| Missing, conflict, or non-supervised label | 245 |

The primary binary supervised slice contains:

| Class | Rows |
|---|---:|
| Benign / likely benign | 33,859 |
| Pathogenic / likely pathogenic | 9,131 |
| Total | 42,990 |

The exploratory three-class slice contains:

| Class | Rows |
|---|---:|
| Benign / likely benign | 33,859 |
| VUS | 41,731 |
| Pathogenic / likely pathogenic | 9,131 |
| Total | 84,721 |

## S2. Gene-Panel Quality Control

The gene panel was not accepted as a raw union of all source genes. We audited genes for cardiac relevance, disease mechanism, label quality, sparse pathogenic evidence, and possible source artifacts.

Primary concerns included:

| Gene or gene group | Issue | Handling |
|---|---|---|
| AKAP9 | Weak or disputed long-QT evidence in current ClinGen-style interpretation; almost no P/LP signal in the local ClinVar slice | Excluded from primary model |
| DMD | Primary neuromuscular disease gene; cardiac disease is often secondary to Duchenne or Becker muscular dystrophy | Excluded from primary cardiogenetics model |
| FPGT and FPGT-TNNI3K | Read-through or locus artifact risk | Excluded or quarantined |
| KNCH2-like artifacts | Typographic or alias issue for KCNH2 | Excluded after audit |
| TTN | Important for truncating DCM variants, but missense variation is massive and VUS-heavy | Separate or special handling |
| ANK2, SOS1, TNNI3K, TRPM4 | Real or possible cardiac relevance but sparse or mixed evidence | Sensitivity, monitor, or subgroup handling |

This audit was needed because ClinVar gene slices can pull in variants submitted for non-cardiac disease contexts. The final primary model therefore used explicit inclusion flags and source-specific QC, not simple gene-name matching.

## S3. HiRO, eMERGE, and CardioBoost Source-Record Handling

### S3.1 HiRO

HiRO/CASPER WES/VERDICT contributed 482 private source records. These records include genotype, phenotype, and ACMG-style adjudication fields. We preserved source records because repeated genotypes can correspond to different patients, phenotypes, or classification events.

Current HiRO linkage:

| HiRO status | Count |
|---|---:|
| Total source records | 482 |
| Linked source records with model prediction | 408 |
| Unresolved source records | 74 |
| Resolved variant-level rows | 240 |

The 240 variant rows are not the full value of HiRO. The source-record analysis remains important because the same resolved variant can appear in multiple patient-level records.

HiRO source-record deferral confusion for the primary binary model:

| True label | Benign-like | Deferred | Pathogenic-like |
|---|---:|---:|---:|
| Benign | 149 | 23 | 18 |
| VUS | 60 | 74 | 38 |
| Pathogenic | 0 | 5 | 41 |

Key point: among predicted HiRO source records, no pathogenic source record was called benign-like.

### S3.2 eMERGE

eMERGE contributed 2,754 arrhythmia-gene rows. An initial poor external result was traced to coordinate mismatch. We then lifted eMERGE coordinates from GRCh37 to GRCh38, checked the reference allele against the GRCh38 FASTA, and reran VEP annotation.

eMERGE coordinate rescue:

| QC item | Result |
|---|---:|
| Source rows processed | 2,754 |
| Single-mapping lifted rows | 2,754 |
| GRCh38 REF allele matches | 2,754 |
| Coordinate QC pass | 2,754 |
| VEP returned after rescue | 2,754 |

Source-versus-VEP consequence agreement:

| Agreement class | Rows |
|---|---:|
| Broad class match | 2,578 |
| Missing source or VEP consequence | 145 |
| Partial LoF-related match | 23 |
| Broad mismatch | 8 |

The earlier eMERGE failure should therefore be treated as invalid. The corrected eMERGE set is a valid external source stratum.

### S3.3 CardioBoost

CardioBoost contributed 355 public source records and 353 linked variant rows. Many CardioBoost rows overlap ClinVar coordinates, so source-held-out splitting and matched benchmarking were required. CardioBoost has no VUS labels in the processed local evaluation table, so it is not a true three-class external validation set.

## S4. Feature Engineering and Coverage

The final ready matrix used feature families that capture different evidence types.

| Feature family | Main contribution | Coverage in 85,677-row matrix |
|---|---|---:|
| VEP consequence and HGVS | Transcript consequence, protein change, molecular impact | 99.98% |
| SpliceAI | Splice donor and acceptor gain/loss scores | 90.04% |
| gnomAD | Allele frequency, popmax AF, observed status, homozygote count | 73.46% final AF |
| dbNSFP | Missense predictors and conservation scores | 51.09% |
| AlphaMissense | Protein-level missense effect score | 44.38% direct |
| ClinGen | Gene-disease validity and limited variant-level evidence | High for gene-level, sparse for variant-level |
| Protein structure | UniProt domains, AlphaFold pLDDT, DSSP, FreeSASA | Sparse, position-dependent |
| FoldX | Stability change for compatible missense variants | 0.37% in ready matrix |

Feature coverage:

| Feature family | Covered rows | Coverage |
|---|---:|---:|
| VEP | 85,660 / 85,677 | 99.98% |
| SpliceAI | 77,141 / 85,677 | 90.04% |
| gnomAD final AF | 62,937 / 85,677 | 73.46% |
| gnomAD observed | 55,387 / 85,677 | 64.65% |
| dbNSFP | 43,774 / 85,677 | 51.09% |
| AlphaMissense direct | 38,023 / 85,677 | 44.38% |
| FoldX DDG | 313 / 85,677 | 0.37% |

dbNSFP does not cover every row because many variants are non-missense, indels, splice-region variants, or otherwise not represented in the dbNSFP missense-oriented records. This is expected and was handled with missingness-aware CatBoost features.

## S5. Model Training, Weighting, and Split Assignment

The primary model was binary: P/LP versus B/LB. This was selected as the main model because VUS is an evidence-status label. VUS rows were scored after training but not used as positive or negative ground truth in the primary supervised task.

Weights combined confidence and class imbalance:

| Row type | Weighting logic |
|---|---|
| Pathogenic / likely pathogenic | Upweighted because the class is smaller and false negatives are clinically costly |
| ClinVar 2-star or 3-star | Higher confidence than 1-star |
| ClinVar 1-star | Retained, lower confidence contribution |
| HiRO patient-linked records | High-trust source, reported separately at source-record level |
| eMERGE and CardioBoost | Source-specific handling, leakage-aware splitting |
| Conflicts | Excluded from supervised binary training |
| VUS | Excluded from supervised binary training, scored afterward |

Three split strategies were created:

| Split | Purpose |
|---|---|
| Internal grouped split | Development, debugging, and hyperparameter search |
| Source-held-out split | Main manuscript validation |
| Gene-stress split | Sparse-gene and gene-transfer stress testing |

Source-held-out binary split:

| Split | Benign | Pathogenic | Rows |
|---|---:|---:|---:|
| Train | 28,631 | 7,545 | 36,176 |
| Validation | 5,053 | 1,331 | 6,384 |
| External HiRO | 42 | 27 | 69 |
| External eMERGE | 80 | 96 | 176 |
| External CardioBoost | 53 | 132 | 185 |

## S6. Primary Binary Model Performance

Standard threshold 0.5:

| Model / split | Rows | AUROC | AUPRC | Sensitivity | Specificity | PPV |
|---|---:|---:|---:|---:|---:|---:|
| Internal grouped test | 6,449 | 0.9990 | 0.9971 | 0.9883 | 0.9968 | 0.9883 |
| Source-held-out HiRO | 69 | 0.9877 | 0.9830 | 0.9630 | 0.8810 | 0.8387 |
| Source-held-out eMERGE | 176 | 0.9910 | 0.9914 | 0.9792 | 0.9375 | 0.9495 |
| Source-held-out CardioBoost | 185 | 0.9520 | 0.9812 | 0.9545 | 0.8113 | 0.9265 |
| Gene-stress sparse-gene test | 5,890 | 0.9979 | 0.9574 | 0.9848 | 0.9947 | 0.8966 |
| Gene-stress well-represented test | 5,565 | 0.9999 | 0.9996 | 0.9872 | 0.9976 | 0.9924 |

CardioBoost-style 0.1/0.9 deferral thresholds:

| Set | Rows | Sensitivity | Specificity | PPV | NPV | MCC | Deferral |
|---|---:|---:|---:|---:|---:|---:|---:|
| External combined | 430 | 0.937 | 0.943 | 0.960 | 0.977 | 0.876 | 0.119 |
| HiRO | 69 | 0.926 | 0.929 | 0.893 | 1.000 | 0.849 | 0.087 |
| eMERGE | 176 | 0.948 | 0.988 | 0.989 | 1.000 | 0.933 | 0.114 |
| CardioBoost | 185 | 0.932 | 0.887 | 0.953 | 0.903 | 0.806 | 0.135 |

External combined confusion matrix under deferral:

| True label | Benign-like | Deferred | Pathogenic-like |
|---|---:|---:|---:|
| Benign | 127 | 38 | 10 |
| Pathogenic | 3 | 13 | 239 |

## S7. All-Label VUS Triage

The binary model was applied to benign, VUS, and pathogenic rows using the same probability zones:

| Probability | Output |
|---:|---|
| <= 0.1 | Benign-like |
| 0.1 to 0.9 | Deferred / uncertain |
| >= 0.9 | Pathogenic-like |

This is a clinical triage analysis, not a true three-class classifier.

| Dataset | Rows with prediction | True benign | True VUS | True pathogenic | Exact three-zone agreement | P/LP sensitivity | VUS deferral |
|---|---:|---:|---:|---:|---:|---:|---:|
| All variant rows | 85,677 | 33,923 | 42,361 | 9,148 | 0.670 | 0.981 | 0.178 |
| ClinVar variant rows | 83,915 | 33,796 | 40,872 | 9,038 | 0.674 | 0.982 | 0.174 |
| HiRO variant rows | 240 | 43 | 139 | 31 | 0.568 | 0.903 | 0.325 |
| eMERGE variant rows | 2,754 | 130 | 2,394 | 111 | 0.498 | 0.946 | 0.427 |
| CardioBoost variant rows | 353 | 57 | 0 | 137 | 0.820 | 0.934 | 0.261 |
| HiRO source records | 408 | 190 | 172 | 46 | 0.647 | 0.891 | 0.250 |
| eMERGE source records | 2,754 | 137 | 2,490 | 127 | 0.497 | 0.945 | 0.427 |
| CardioBoost source records | 355 | 156 | 0 | 199 | 0.642 | 0.869 | 0.259 |

Important hard-error counts:

| Dataset | P/LP to benign-like | Benign to pathogenic-like | VUS to pathogenic-like |
|---|---:|---:|---:|
| All variant rows | 22 | 43 | 18,020 |
| ClinVar variant rows | 19 | 36 | 17,666 |
| HiRO source records | 0 | 18 | 38 |
| eMERGE source records | 0 | 1 | 738 |
| CardioBoost source records | 5 | 30 | 0 |

VUS-to-pathogenic-like calls should be interpreted as prioritization signals, not ground-truth errors.

## S8. True Three-Class Model

A secondary CatBoost model was trained to directly predict Benign, VUS, and Pathogenic.

Training data:

| Class | Rows |
|---|---:|
| Benign | 33,859 |
| VUS | 41,731 |
| Pathogenic | 9,131 |
| Total | 84,721 |

Class weights:

| Class | Weight |
|---|---:|
| Benign | 1.0 |
| VUS | 0.9 |
| Pathogenic | 2.5 |

Source-held-out three-class split:

| Split | Benign | VUS | Pathogenic | Rows |
|---|---:|---:|---:|---:|
| Train | 28,631 | 33,869 | 7,545 | 70,045 |
| Validation | 5,053 | 5,977 | 1,331 | 12,361 |
| External HiRO | 42 | 110 | 27 | 179 |
| External eMERGE | 80 | 1,775 | 96 | 1,951 |
| External CardioBoost | 53 | 0 | 132 | 185 |

External results:

| External set | Rows | Accuracy | Macro-F1 | Weighted kappa | MCC | OVR macro AUROC |
|---|---:|---:|---:|---:|---:|---:|
| HiRO | 179 | 0.8492 | 0.8188 | 0.7879 | 0.7160 | 0.9426 |
| eMERGE | 1,951 | 0.8893 | 0.6507 | 0.4977 | 0.4713 | 0.9242 |
| CardioBoost | 185 | 0.5946 | 0.4693 | 0.6351 | 0.4339 | NA |

CardioBoost has no VUS rows, so it is not a valid true three-class validation set.

HiRO source-record three-class confusion, predicted rows only:

| True label | Benign | VUS | Pathogenic |
|---|---:|---:|---:|
| Benign | 113 | 77 | 0 |
| VUS | 16 | 153 | 3 |
| Pathogenic | 0 | 16 | 30 |

The true three-class model improves exact VUS label matching but has lower P/LP sensitivity than the binary deferral model. Therefore, it is reported as secondary.

## S9. Matched Public CardioBoost Benchmark

The official public CardioBoost model artifacts were downloaded from the public GitHub repository and run locally with the official preprocessing and prediction scripts.

Local official files included:

| File | Role |
|---|---|
| `data/cardiomyopathy/ml/train_ada.RData` | Official cardiomyopathy AdaBoost model |
| `data/arrhythmia/ml/train_ada.RData` | Official arrhythmia AdaBoost model |
| `data/cardiomyopathy/preprocess.RData` | Official cardiomyopathy preprocessing |
| `data/arrhythmia/preprocess.RData` | Official arrhythmia preprocessing |
| `script/src/predict.R` | Official scoring helper |
| `script/src/preprocess_test.R` | Official preprocessing helper |

Official CardioBoost prediction outputs:

| Output | Rows |
|---|---:|
| Official cardiomyopathy all-rare predictions | 65,475 |
| Official arrhythmia all-rare predictions | 42,411 |

Eligibility counts from our matrix:

| Slice | Rows | Benign | Pathogenic | VUS |
|---|---:|---:|---:|---:|
| All current matrix rows | 85,677 | 33,923 | 9,148 | 42,361 |
| CardioBoost-gene rows | 37,528 | 14,459 | 5,670 | 17,186 |
| CardioBoost-gene missense rows | 18,251 | 1,287 | 2,296 | 14,465 |
| CardioBoost-gene missense binary rows | 3,583 | 1,287 | 2,296 | 0 |
| CardioBoost-gene missense strict-rare binary rows | 1,527 | 1,073 | 454 | 0 |
| Exact coordinate matches to official all-rare universe | 64 | 34 | 8 | 22 |
| Exact coordinate matches, binary only | 42 | 34 | 8 | 0 |

Exact coordinate matching was too small for the main benchmark. The primary benchmark used HGVS cDNA matching, with coordinate matching retained as QC.

Primary external HGVS cDNA-matched benchmark:

| Model | Rows | Pathogenic | Benign | AUROC | AUPRC | Sensitivity | Specificity | PPV | NPV | Deferral | High-confidence accuracy |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Our CatBoost | 250 | 169 | 81 | 0.980 | 0.990 | 0.959 | 0.543 | 0.964 | 1.000 | 0.152 | 0.972 |
| Official CardioBoost | 250 | 169 | 81 | 0.843 | 0.902 | 0.876 | 0.148 | 0.836 | 1.000 | 0.244 | 0.847 |
| REVEL local | 247 | 166 | 81 | 0.903 | 0.940 | 0.530 | 0.037 | 0.936 | 1.000 | 0.607 | 0.938 |
| CADD PHRED | 250 | 169 | 81 | 0.872 | 0.914 | 1.000 | 0.012 | 0.693 | 1.000 | 0.020 | 0.694 |
| AlphaMissense | 248 | 167 | 81 | 0.937 | 0.969 | 0.635 | 0.457 | 0.991 | 0.902 | 0.403 | 0.966 |

Paired bootstrap differences, our model minus official CardioBoost:

| Metric | Difference | 95% CI |
|---|---:|---:|
| AUROC | 0.138 | 0.083 to 0.195 |
| AUPRC | 0.087 | 0.048 to 0.132 |
| Sensitivity | 0.083 | 0.030 to 0.136 |
| Specificity | 0.395 | 0.272 to 0.519 |
| PPV | 0.128 | 0.082 to 0.172 |
| Deferral | -0.092 | -0.156 to -0.028 |
| High-confidence accuracy | 0.125 | 0.082 to 0.167 |

This benchmark is a same-row public-model comparison. It is not a reconstruction of the private CardioBoost development cohort.

## S10. Leakage and Ablation Audit

The leakage audit retrained the primary binary model after removing sensitive feature families.

| Ablation | Features | Dropped | Description |
|---|---:|---:|---|
| Baseline | 338 | 0 | Current primary feature-selection rules |
| No gene identity | 328 | 10 | Removes direct gene and gene-level identity features |
| No VEP consequence/HGVS | 302 | 36 | Removes VEP consequence, impact, transcript, and HGVS-like annotations |
| No ClinVar metadata | 338 | 0 | Confirms residual ClinVar metadata was already excluded |
| No source-record aggregates | 300 | 38 | Removes HiRO/source aggregate features |
| No HGVS/transcript identifiers | 316 | 22 | Removes HGVS, transcript, protein accession, rsID/CAID, and identifier-like fields |
| Strict low-leakage | 246 | 92 | Removes all high-risk families together |

External-all performance:

| Ablation | AUROC | AUPRC | Sensitivity | Specificity | PPV | Deferral | High-confidence accuracy |
|---|---:|---:|---:|---:|---:|---:|---:|
| Baseline | 0.977 | 0.984 | 0.973 | 0.897 | 0.932 | 0.098 | 0.966 |
| No gene identity | 0.977 | 0.984 | 0.973 | 0.886 | 0.925 | 0.100 | 0.966 |
| No VEP consequence/HGVS | 0.973 | 0.980 | 0.957 | 0.886 | 0.924 | 0.160 | 0.967 |
| No ClinVar metadata | 0.977 | 0.984 | 0.973 | 0.897 | 0.932 | 0.098 | 0.966 |
| No source-record aggregates | 0.977 | 0.983 | 0.973 | 0.891 | 0.929 | 0.102 | 0.966 |
| No HGVS/transcript identifiers | 0.978 | 0.985 | 0.973 | 0.914 | 0.943 | 0.093 | 0.964 |
| Strict low-leakage | 0.974 | 0.981 | 0.973 | 0.811 | 0.883 | 0.212 | 0.973 |

The audit supports the conclusion that the primary results are not explained by one obvious leakage-prone feature family.

## S11. Calibration and Stratified Performance

Calibration summary:

| Subset | Rows | Brier | Mean predicted P/LP probability | Observed P/LP fraction | Mean error |
|---|---:|---:|---:|---:|---:|
| Validation | 6,384 | 0.003 | 0.211 | 0.208 | 0.003 |
| External all | 430 | 0.053 | 0.627 | 0.593 | 0.034 |
| HiRO | 69 | 0.068 | 0.444 | 0.391 | 0.052 |
| eMERGE | 176 | 0.030 | 0.575 | 0.545 | 0.029 |
| CardioBoost | 185 | 0.069 | 0.745 | 0.714 | 0.031 |

The external Brier scores and mean errors suggest mild overconfidence. The model should be presented as a prioritization model unless validation-only calibration is added and externally tested.

Stratified external-all results:

| Stratum | Rows | AUROC | AUPRC | Sensitivity | Specificity | PPV |
|---|---:|---:|---:|---:|---:|---:|
| Missense | 315 | 0.987 | 0.989 | 0.977 | 0.887 | 0.914 |
| Stop gained | 32 | 0.867 | 0.990 | 1.000 | 0.500 | 0.968 |
| dbNSFP matched | 348 | 0.985 | 0.989 | 0.981 | 0.896 | 0.938 |
| dbNSFP missing | 82 | 0.949 | 0.962 | 0.878 | 0.854 | 0.857 |
| Arrhythmia genes | 265 | 0.985 | 0.988 | 0.975 | 0.907 | 0.939 |
| Cardiomyopathy genes | 151 | 0.965 | 0.980 | 0.946 | 0.847 | 0.906 |
| gnomAD observed | 326 | 0.984 | 0.984 | 0.977 | 0.901 | 0.919 |
| gnomAD not joined | 72 | 0.997 | 1.000 | 0.985 | 1.000 | 1.000 |
| gnomAD confirmed absent | 32 | 0.787 | 0.839 | 0.692 | 0.737 | 0.643 |
| Well-represented genes | 406 | 0.976 | 0.985 | 0.963 | 0.882 | 0.925 |
| Sparse P/LP genes | 24 | 0.986 | 0.981 | 1.000 | 0.929 | 0.909 |

The gnomAD-confirmed-absent stratum is small. It is useful as a cautionary result but should not be overinterpreted.

## S12. Reproducibility Files and Scripts

Main local artifacts:

| Artifact | Path |
|---|---|
| Current manuscript progress | `PROJECT_PROPOSAL-PROGRESS.md` |
| Primary manuscript draft | `MANUSCRIPT_DRAFT_CURRENT.md` |
| Updated GIM draft | `manuscript/GIM_MANUSCRIPT_DRAFT_UPDATED.md` |
| Updated SI draft | `manuscript/GIM_SUPPLEMENTARY_INFORMATION_UPDATED.md` |
| Primary binary results | `results/model_performance/binary_Results.MD` |
| Three-class results | `results/model_performance/three_tier_results.MD` |
| Source-held-out review | `results/model_performance/source_heldout_hiro_emerge_rescued/PERFORMANCE_REVIEW.md` |
| CardioBoost benchmark | `results/model_performance/cardioboost_matched_benchmark/CARDIOBOOST_MATCHED_BENCHMARK_REPORT.md` |
| Leakage audit | `results/model_performance/leakage_ablation_audit/LEAKAGE_ABLATION_AUDIT_REPORT.md` |
| Data characteristics | `results/data_characteristics/DATA_CHARACTERISTICS_REPORT.md` |
| Manuscript figures and tables | `results/manuscript_figures_tables/` |

Important generated figure tables:

| Output | Path |
|---|---|
| External validation table | `results/manuscript_figures_tables/tables/external_validation_performance.tsv` |
| Calibration table | `results/manuscript_figures_tables/tables/calibration_brier_summary.tsv` |
| Stratified performance | `results/manuscript_figures_tables/tables/stratified_performance.tsv` |
| All-label triage summary | `results/manuscript_figures_tables/tables/primary_binary_three_zone_all_labels_summary.tsv` |
| Feature coverage | `results/manuscript_figures_tables/tables/feature_group_coverage_manuscript.tsv` |

CardioBoost benchmark commands:

```bash
Rscript results/model_performance/cardioboost_matched_benchmark/run_official_cardioboost_predictions.R
python3 results/model_performance/cardioboost_matched_benchmark/build_exact_match_benchmark.py
python3 results/model_performance/cardioboost_matched_benchmark/plot_cardiboost_benchmark.py
```
