# A Source-Held-Out CatBoost Model for Cardiogenetic Variant Prioritization

## Current Readiness Statement

The project is now strong enough for a manuscript built around one main claim:
a large, feature-rich, cardiogenetics-specific binary CatBoost model with an
explicit VUS deferral zone generalizes across source-held-out HiRO, eMERGE, and
CardioBoost strata, and outperforms the official public CardioBoost model on a
same-row, CardioBoost-eligible matched benchmark.

This claim is intentionally narrow. It does not say we reproduced the full
private CardioBoost development study. It says we ran the public CardioBoost
models locally, matched their scores to the same external rows scored by our
model, and compared both models under the same high-confidence threshold
framework.

Calibration analysis, Brier scores, reliability plots, stratified performance
tables, leakage/ablation audits, all-label clinical triage tables, and paired
bootstrap confidence intervals for the CardioBoost comparison have now been
generated. The remaining work is manuscript polish, optional calibration
testing, final multipanel figure layout, and deeper rescue of the unresolved
HiRO source records.

## Count Reconciliation

Two row counts appear in the project because they describe different analysis
stages.

| Matrix | Rows | Columns | Meaning |
|---|---:|---:|---|
| Annotated intermediate feature matrix | 86,889 | 463 | Large GRCh38 registry after feature annotation, before final ready-table cleanup |
| Analysis-ready modeling matrix | 85,677 | 420 | Current matrix used for final splits, weights, external validation, calibration, stratified analyses, and all-label triage |

All current manuscript-facing results use the 85,677-row analysis-ready matrix.
The 86,889-row matrix should be described only as an intermediate annotated
registry.

## Abstract

Variant interpretation remains a major bottleneck in inherited cardiomyopathy
and arrhythmia genetics. Existing disease-specific tools, including
CardioBoost, focus mainly on rare missense variants and use a binary pathogenic
versus benign framework with an indeterminate score zone for uncertain cases
(Zhang et al., 2021). However, current public and local resources now support a
larger cardiogenetics-specific model with broader feature coverage and multiple
external validation strata.

We built a GRCh38-centered cardiogenetics registry integrating ClinVar,
HiRO/CASPER WES/VERDICT, eMERGE-III arrhythmia variants, and the public
CardioBoost dataset. The current analysis-ready matrix contains 85,677 unique
variant rows. It includes 33,923 benign/likely benign rows, 9,148
pathogenic/likely pathogenic rows, 42,361 VUS rows, and 245 missing, conflicting,
or non-supervised rows. The primary binary training slice contains 42,990 rows:
33,859 benign and 9,131 pathogenic. VUS rows were excluded from supervised
binary training and scored afterward.

The primary model is a binary CatBoost classifier trained to separate
pathogenic/likely pathogenic from benign/likely benign variants. A
CardioBoost-style three-zone rule was then applied to the predicted pathogenic
probability: `<=0.1` benign-like, `0.1-0.9` deferred, and `>=0.9`
pathogenic-like. On source-held-out binary external validation, the model
achieved AUROC 0.977 and AUPRC 0.985 across 430 external rows. At the standard
0.5 threshold, sensitivity was 0.965, specificity 0.886, and PPV 0.925. Under
the 0.1/0.9 high-confidence rule, the external deferral rate was 0.119 and
high-confidence accuracy was 0.966.

On the same-row CardioBoost-matched benchmark, our model outperformed the
official public CardioBoost model on 250 CardioBoost-eligible external rows:
AUROC 0.980 versus 0.843, AUPRC 0.990 versus 0.902, sensitivity 0.959 versus
0.876, PPV 0.964 versus 0.836, deferral 0.152 versus 0.244, and high-confidence
accuracy 0.972 versus 0.847. Paired bootstrap confidence intervals supported
these differences. The paired AUROC difference was 0.138 with 95% CI
0.083-0.195, and the paired high-confidence accuracy difference was 0.125 with
95% CI 0.082-0.167.

These results support a manuscript centered on binary cardiogenetic
pathogenicity prediction with explicit VUS deferral, source-held-out validation,
and direct public-model benchmarking against CardioBoost.

## Introduction

Clinical genetic testing for inherited cardiomyopathies and arrhythmia
syndromes frequently returns variants whose interpretation is difficult. This is
especially true for rare missense variants and variants of uncertain
significance. ACMG/AMP guidelines provide a structured framework, but the
evidence required for confident classification is incomplete for many variants
(Richards et al., 2015). ClinVar provides broad public access to submitted
variant classifications, but review status, submitter agreement, and disease
context vary across records (Landrum et al., 2018).

CardioBoost is the most direct prior comparator for this project. It is a
disease-specific classifier for rare missense variants in inherited
cardiomyopathy and inherited arrhythmia genes. It estimates binary pathogenicity
probability, then uses two clinical thresholds: `>=0.9` for disease-causing,
`<=0.1` for benign or likely benign, and the middle range as indeterminate
(Zhang et al., 2021). This matters because CardioBoost is not a true
Benign/VUS/Pathogenic classifier. Its VUS-like output is a deferral zone.

The current project extends that framework in three practical ways. First, it
uses a larger integrated registry. Second, it adds modern annotations, including
gnomAD v4.1, dbNSFP, AlphaMissense, VEP/HGVS, SpliceAI, ClinGen context, and
selected protein-structure features from AlphaFold, UniProt, DSSP, FreeSASA,
and FoldX. Third, it tests generalization on source-held-out HiRO, eMERGE, and
CardioBoost strata instead of relying only on an internal random split.

The main question is concrete: can a larger cardiogenetics-specific model
improve high-confidence variant prioritization while preserving a clinically
safe uncertain zone?

## Methods

### Data Sources

ClinVar was used as the main public label source after germline filtering,
conflict handling, and gene-panel QC. The current analysis includes 83,915
ClinVar variant rows with model predictions. ClinVar labels were mapped to
benign/likely benign, VUS, and pathogenic/likely pathogenic. Review stars were
retained for confidence weighting and sensitivity analysis, but not used as
biological model features.

HiRO/CASPER WES/VERDICT contributed 482 private source records with genotype,
phenotype fields, and ACMG-style adjudication. These records are not treated as
simple duplicate variants. Repeated genotypes can represent different patient
events, phenotypes, or classifications. After coordinate rescue, 408 source
records link to the model registry and 74 remain unresolved. These linked
records correspond to 240 variant-level matrix rows.

eMERGE contributed 2,754 arrhythmia-gene source records. An initial weak
external result was traced to coordinate build mismatch. After GRCh37 to GRCh38
rescue, all 2,754 rows lifted with single mappings, passed GRCh38 REF allele
checks, and received VEP annotations. eMERGE is therefore treated as a real
external variant source, not a supplemental annotation set.

CardioBoost contributed 355 public source records and 353 linked variant rows.
This local CardioBoost table is binary in the current processed cohort, with no
VUS rows. It is useful for external testing and for public-model benchmarking,
but leakage control is essential because many CardioBoost variants overlap with
ClinVar.

### Gene Panel QC

The model does not use a raw source-union gene panel. Genes with likely
contamination, weak cardiac evidence, or unclear scope were removed or
quarantined. AKAP9, DMD, FPGT, FPGT-TNNI3K, KCNE1B, and KNCH2 were excluded
from the primary panel. ANK2, SOS1, TNNI3K, TRPM4, and TTN were assigned to
sensitivity or special handling. TTN is handled carefully because truncating
variants are important in dilated cardiomyopathy, while TTN missense variation
is large, VUS-heavy, and biologically different from the rare missense task that
CardioBoost addressed.

### Feature Engineering

The analysis-ready matrix contains 85,677 unique GRCh38 variant rows and 420
columns. The primary CatBoost model used 338 selected non-leakage feature
columns after excluding labels, raw source labels, review metadata as
predictors, source IDs, patient identifiers, and direct variant identifiers.

Feature coverage across the ready matrix was:

| Feature family | Covered rows | Coverage |
|---|---:|---:|
| VEP consequence/HGVS | 85,660 / 85,677 | 99.98% |
| SpliceAI | 77,141 / 85,677 | 90.04% |
| gnomAD final allele frequency | 62,937 / 85,677 | 73.46% |
| gnomAD observed status | 55,387 / 85,677 | 64.65% |
| dbNSFP | 43,774 / 85,677 | 51.09% |
| AlphaMissense direct | 38,023 / 85,677 | 44.38% |
| FoldX DDG | 313 / 85,677 | 0.37% |

This coverage pattern is important. The model is powered mainly by VEP,
SpliceAI, gnomAD, dbNSFP, AlphaMissense, conservation, and ClinGen-related
features. Protein-structure and FoldX features are valuable for mapped missense
variants, but they are sparse in the full matrix and should be framed as
high-value subfeatures, not as universal coverage.

### Model Formulations

The primary model is binary: pathogenic/likely pathogenic versus benign/likely
benign. VUS rows were excluded from supervised binary training because VUS is an
evidence status, not a stable biological class. A VUS can later become benign or
pathogenic when new evidence appears. Training the primary model to learn VUS
directly risks teaching it to predict evidence gaps.

Sample weights combined source confidence and class imbalance. Pathogenic rows
were upweighted because false negatives are clinically more costly and because
pathogenic rows are less common than benign rows in the binary training slice.
The final primary binary training slice contained 42,990 rows.

A secondary true three-class CatBoost model was also trained on Benign, VUS, and
Pathogenic labels. It is useful for VUS behavior analysis and direct label
matching, but it is not the correct direct comparator to CardioBoost.

### Splits and Leakage Control

Three split strategies were implemented.

| Split | Purpose |
|---|---|
| Internal grouped split | Development, debugging, and hyperparameter search |
| Source-held-out split | Main external validation on HiRO, eMERGE, and CardioBoost |
| Gene-stress split | Sparse-gene and gene-transfer stress testing |

The source-held-out split is the paper-facing split. It trains on non-held-out
rows, validates internally, and tests externally on HiRO, eMERGE, and
CardioBoost. Splits are grouped by variant identity so the same
chrom-pos-ref-alt variant cannot appear on both sides under different source
names.

A leakage-oriented ablation audit retrained the model after removing sensitive
feature families: gene identity, ClinVar metadata, source-record aggregates,
VEP/consequence/HGVS fields, HGVS/transcript identifiers, and all of these
together in a strict low-leakage stress test.

## Results

### Dataset Size and Source Composition

The final ready matrix contains 85,677 unique variant rows. The supervised
labels are:

| Label | Rows |
|---|---:|
| Benign / likely benign | 33,923 |
| VUS | 42,361 |
| Pathogenic / likely pathogenic | 9,148 |
| Missing, conflict, or other non-supervised label | 245 |

The primary binary training pool contains 42,990 rows: 33,859 benign and 9,131
pathogenic. The exploratory three-class pool contains 84,721 rows: 33,859
benign, 41,731 VUS, and 9,131 pathogenic.

| Source | Source records | Linked/scored source records | Variant rows with prediction |
|---|---:|---:|---:|
| ClinVar | 83,915 | 83,915 | 83,915 |
| HiRO | 482 | 408 | 240 |
| eMERGE | 2,754 | 2,754 | 2,754 |
| CardioBoost | 355 | 355 | 353 |

HiRO is reported at both source-record and variant levels because repeated
genotypes may represent distinct patient-level evidence. This is why the
variant-level HiRO count is smaller than the source-record count.

### Primary Binary Model

At a standard 0.5 threshold, the source-held-out binary model performed well on
all external sets:

| External set | Rows | AUROC | AUPRC | Sensitivity | Specificity | PPV | Brier |
|---|---:|---:|---:|---:|---:|---:|---:|
| Combined external | 430 | 0.977 | 0.985 | 0.965 | 0.886 | 0.925 | 0.053 |
| HiRO | 69 | 0.988 | 0.983 | 0.963 | 0.881 | 0.839 | 0.068 |
| eMERGE | 176 | 0.991 | 0.991 | 0.979 | 0.938 | 0.949 | 0.030 |
| CardioBoost | 185 | 0.952 | 0.981 | 0.955 | 0.811 | 0.926 | 0.069 |

Using the CardioBoost-style `0.1/0.9` high-confidence thresholds, external
performance remained strong:

| External set | Rows | Sensitivity | Specificity | PPV | NPV | Deferral rate |
|---|---:|---:|---:|---:|---:|---:|
| Combined external | 430 | 0.937 | 0.943 | 0.960 | 0.977 | 0.119 |
| HiRO | 69 | 0.926 | 0.929 | 0.893 | 1.000 | 0.087 |
| eMERGE | 176 | 0.948 | 0.988 | 0.989 | 1.000 | 0.114 |
| CardioBoost | 185 | 0.932 | 0.887 | 0.953 | 0.903 | 0.135 |

These results are the strongest evidence that the model generalizes. They also
make the model directly comparable in concept to CardioBoost because both use a
binary pathogenicity score and an indeterminate zone.

### Calibration

Calibration was evaluated with Brier score and reliability curves.

| Subset | Rows | Brier | Mean predicted P/LP probability | Observed P/LP fraction | Mean error |
|---|---:|---:|---:|---:|---:|
| Validation | 6,384 | 0.003 | 0.211 | 0.208 | 0.003 |
| External all | 430 | 0.053 | 0.627 | 0.593 | 0.034 |
| HiRO | 69 | 0.068 | 0.444 | 0.391 | 0.052 |
| eMERGE | 176 | 0.030 | 0.575 | 0.545 | 0.029 |
| CardioBoost | 185 | 0.069 | 0.745 | 0.714 | 0.031 |

The model is very well calibrated on internal validation. It is mildly
overconfident on external data, especially in the small HiRO set. This supports
using the current probabilities as ranking and triage scores. If the final
paper makes strong claims about absolute probability calibration, validation-set
Platt or isotonic calibration should be tested before final submission.

### Stratified Performance

Stratified external-all performance was generated by gene group, sparse-gene
status, consequence group, dbNSFP coverage, and gnomAD status.

| Stratum | Rows | AUROC | AUPRC | Sensitivity | Specificity | PPV |
|---|---:|---:|---:|---:|---:|---:|
| Arrhythmia genes | 265 | 0.985 | 0.988 | 0.975 | 0.907 | 0.939 |
| Cardiomyopathy genes | 151 | 0.965 | 0.980 | 0.946 | 0.847 | 0.906 |
| Missense variants | 315 | 0.987 | 0.989 | 0.977 | 0.887 | 0.914 |
| dbNSFP matched | 348 | 0.985 | 0.989 | 0.981 | 0.896 | 0.938 |
| dbNSFP missing | 82 | 0.949 | 0.962 | 0.878 | 0.854 | 0.857 |
| gnomAD observed | 326 | 0.984 | 0.984 | 0.977 | 0.901 | 0.919 |
| gnomAD confirmed absent | 32 | 0.787 | 0.839 | 0.692 | 0.737 | 0.643 |
| Sparse P/LP genes | 24 | 0.986 | 0.981 | 1.000 | 0.929 | 0.909 |

Performance is strongest in missense, dbNSFP-matched, and gnomAD-observed
variants. dbNSFP-missing rows still perform reasonably but are weaker. The
gnomAD-confirmed-absent stratum is small and performs worse, so it should be
reported as a high-caution descriptive subgroup.

### All-Label Clinical Triage

The binary model was also applied to all Benign, VUS, and Pathogenic rows using
the same three zones. This is not a true three-class classifier. It is a triage
analysis that asks which VUS look benign-like, which remain uncertain, and
which look pathogenic-like.

| Dataset | Rows | True B | True VUS | True P | Exact 3-zone agreement | P/LP sensitivity | Specificity vs non-P/LP | P/LP PPV | VUS deferral |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| All variant rows | 85,677 | 33,923 | 42,361 | 9,148 | 0.670 | 0.981 | 0.763 | 0.332 | 0.178 |
| ClinVar variant rows | 83,915 | 33,796 | 40,872 | 9,038 | 0.674 | 0.982 | 0.763 | 0.334 | 0.174 |
| HiRO variant rows | 240 | 43 | 139 | 31 | 0.568 | 0.903 | 0.813 | 0.452 | 0.325 |
| eMERGE variant rows | 2,754 | 130 | 2,394 | 111 | 0.498 | 0.946 | 0.715 | 0.127 | 0.427 |
| CardioBoost variant rows | 353 | 57 | 0 | 137 | 0.820 | 0.934 | 0.895 | 0.955 | 0.261 |
| HiRO source records | 408 / 482 | 190 | 172 | 46 | 0.647 | 0.891 | 0.845 | 0.423 | 0.250 |

The all-label table is intentionally harsher than the binary endpoint because
it treats VUS deferral as the expected output. Many VUS rows are pushed into
benign-like or pathogenic-like zones. That should be framed as prioritization,
not as ordinary classification error.

The most clinically important hard-error pattern is favorable. Across all
variant rows, 22 of 9,148 P/LP rows were called benign-like. In HiRO source
records, 0 of 46 pathogenic records were called benign-like.

### Secondary True Three-Class Model

The true three-class model was trained on 84,721 rows: 33,859 benign, 41,731
VUS, and 9,131 pathogenic. It directly predicts Benign, VUS, or Pathogenic.

The three-class model had stronger exact label matching for VUS but lower
pathogenic sensitivity than the binary deferral model. On HiRO source records,
the binary deferral model had exact three-zone agreement 0.647 and P/LP
sensitivity 0.891. The true three-class model had exact three-class accuracy
0.766 and P/LP sensitivity 0.714. Therefore, the three-class model is useful as
a secondary VUS behavior analysis, but the binary deferral model remains the
primary clinical triage model.

### Leakage and Ablation Audit

Reviewer-sensitive feature families were removed one at a time and together.
External-all AUROC/AUPRC remained stable:

| Ablation | External rows | AUROC | AUPRC | Sensitivity | Specificity | PPV | Deferral |
|---|---:|---:|---:|---:|---:|---:|---:|
| Baseline current | 430 | 0.977 | 0.984 | 0.973 | 0.897 | 0.932 | 0.098 |
| No gene identity | 430 | 0.977 | 0.984 | 0.973 | 0.886 | 0.925 | 0.100 |
| No VEP/consequence/HGVS | 430 | 0.973 | 0.980 | 0.957 | 0.886 | 0.924 | 0.160 |
| No ClinVar metadata | 430 | 0.977 | 0.984 | 0.973 | 0.897 | 0.932 | 0.098 |
| No source-record aggregates | 430 | 0.977 | 0.983 | 0.973 | 0.891 | 0.929 | 0.102 |
| Strict low-leakage stress test | 430 | 0.974 | 0.981 | 0.973 | 0.811 | 0.883 | 0.212 |

The audit does not show evidence that the main external result is driven by one
obvious leakage-prone feature family. The strict low-leakage model should be
reported as a robustness stress test, not as the deployable model, because it
raises the deferral rate and reduces specificity and PPV.

## Comparison With CardioBoost

CardioBoost reported PR-AUC 0.91 for cardiomyopathy and 0.96 for inherited
arrhythmia. With its 0.1/0.9 thresholds, it reported sensitivity 69.5% and
specificity 56.0% for cardiomyopathy, and sensitivity 83.3% and specificity
78.6% for arrhythmia. Its indeterminate rates were 29.8% and 11.7% (Zhang et
al., 2021).

Our source-held-out binary model with the same threshold concept achieved
combined external sensitivity 93.7%, specificity 94.3%, PPV 96.0%, NPV 97.7%,
and deferral rate 11.9%. On the held-out CardioBoost source stratum, it
achieved sensitivity 93.2%, specificity 88.7%, PPV 95.3%, NPV 90.3%, and
deferral rate 13.5%.

The cleaner direct comparison is the same-row public-model benchmark. We ran
the official public CardioBoost cardiomyopathy and arrhythmia AdaBoost model
objects locally on the official public all-rare mutation tables. The resulting
CardioBoost scores were matched to our CatBoost predictions by gene and HGVS
cDNA. Exact coordinate matching recovered too few binary rows for the primary
comparison, likely due to build and transcript-version differences, so exact
matching was used as QC only.

The primary matched benchmark contained 250 source-held-out rows: 164
CardioBoost-source rows, 73 eMERGE-source rows, and 13 HiRO-source rows. On
these same rows, our CatBoost model achieved AUROC 0.980, AUPRC 0.990,
sensitivity 0.959, specificity 0.543, PPV 0.964, NPV 1.000, deferral 0.152,
and high-confidence accuracy 0.972. Official CardioBoost achieved AUROC 0.843,
AUPRC 0.902, sensitivity 0.876, specificity 0.148, PPV 0.836, NPV 1.000,
deferral 0.244, and high-confidence accuracy 0.847.

Paired bootstrap confidence intervals supported the difference. In the 250-row
benchmark, our model minus official CardioBoost was:

| Metric | Difference | 95% CI |
|---|---:|---:|
| AUROC | 0.138 | 0.083-0.195 |
| AUPRC | 0.087 | 0.048-0.132 |
| Sensitivity | 0.083 | 0.030-0.136 |
| PPV | 0.128 | 0.082-0.172 |
| Deferral | -0.092 | -0.156 to -0.028 |
| High-confidence accuracy | 0.125 | 0.082-0.167 |

In the stricter 164-row CardioBoost-source-only subset, our model achieved
AUROC 0.972 and AUPRC 0.990, compared with AUROC 0.788 and AUPRC 0.905 for
official CardioBoost. The paired AUROC difference was 0.184 with 95% CI
0.101-0.274.

This benchmark supports a public-model same-row comparison. It should not be
described as a reconstruction of CardioBoost's original private development or
clinical outcome cohorts.

## Discussion

The strongest result is not the internal validation score. Internal validation
is useful, but it can be optimistic when ClinVar-like variants and similar genes
appear across train and test. The stronger evidence is source-held-out
generalization to HiRO, eMERGE, and CardioBoost after careful coordinate rescue.

The model also addresses a key clinical problem: how to handle VUS without
pretending VUS is a stable biological class. The primary binary model with a
deferral zone is closer to clinical triage and closer to CardioBoost's design.
The all-label analysis then becomes useful because it shows which VUS are
benign-like, which remain uncertain, and which are pathogenic-like.

HiRO is small but important. Its value comes from patient-linked source records
and expert adjudication, not from large N. The source-record analysis shows that
the binary deferral model did not hard-call any predicted HiRO pathogenic record
as benign-like. That is a clinically meaningful safety signal, although 74 HiRO
records still need coordinate rescue before the full 482-record table can be
scored.

The CardioBoost comparison is now much stronger than a broad narrative
comparison. We tested both models on the same CardioBoost-eligible rows and used
paired bootstrap confidence intervals. The result favors our model across AUROC,
AUPRC, sensitivity, PPV, deferral, and high-confidence accuracy.

The main caution is probability calibration. Ranking and high-confidence triage
are strong, but external Brier scores show mild overconfidence. This is fixable
if absolute risk calibration becomes a manuscript claim. For now, the safest
language is that the model provides a pathogenicity prioritization score with
explicit high-confidence and deferred zones.

## Limitations

First, ClinVar remains the dominant label source. It is essential but
heterogeneous. Review status was used for weighting and sensitivity analysis,
not as a direct biological feature.

Second, dbNSFP and AlphaMissense coverage is incomplete across the full matrix.
This is expected because many variants are not missense-like SNVs. CatBoost can
handle missingness, but performance must be stratified by coverage.

Third, protein-structure and FoldX features are sparse. They should be framed as
mapped-missense subfeatures until protein-position mapping improves.

Fourth, the CardioBoost benchmark is a public-model same-row comparison, not a
full reproduction of CardioBoost's private original study.

Fifth, the true three-class model may learn evidence status as much as biology.
It remains secondary.

## Work Remaining Before Submission

1. Add bootstrap confidence intervals for source-held-out external validation
   and key stratified results, especially small strata such as HiRO and gnomAD
   confirmed-absent variants.
2. Test optional validation-only calibration with Platt scaling or isotonic
   regression if final claims depend on absolute probabilities.
3. Assemble final multipanel manuscript figures from the generated ROC/PR,
   calibration, source-flow, feature-coverage, all-label triage, and stratified
   performance outputs.
4. Continue HiRO unresolved-record rescue. The current audit identified 74
   unlinked source records and two high-confidence manual candidates, but the
   remaining rows need HGVS normalization or manual review.
5. Finalize reviewer-safe wording for the CardioBoost comparison.

## Conclusion

The project has reached the point where the primary model is scientifically
credible. The best current manuscript is not a generic three-class classifier.
It is a cardiogenetics-specific binary CatBoost model with CardioBoost-style
deferral, large-scale public and private source integration, source-held-out
validation, calibration and stratified reporting, leakage stress tests, and a
same-row public-model benchmark against official CardioBoost.

The clean claim is:

> A larger, feature-richer, externally validated cardiogenetics model improved
> high-confidence pathogenicity prioritization on CardioBoost-eligible matched
> external rows while preserving an explicit VUS deferral zone.

## References

Cheng J, Novati G, Pan J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. Science. 2023;381:eadg7492.

Glazer AM, Davogustto G, Shaffer CM, Vanoye CG, et al. Arrhythmia variant associations and reclassifications in the eMERGE-III sequencing study. Circulation. 2022;145:877-891.

Ioannidis NM, Rothstein JH, Pejaver V, et al. REVEL: an ensemble method for predicting the pathogenicity of rare missense variants. American Journal of Human Genetics. 2016;99:877-885.

Jaganathan K, Kyriazopoulou Panagiotopoulou S, McRae JF, et al. Predicting splicing from primary sequence with deep learning. Cell. 2019;176:535-548.e24.

Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction with AlphaFold. Nature. 2021;596:583-589.

Karczewski KJ, Francioli LC, Tiao G, et al. The mutational constraint spectrum quantified from variation in 141,456 humans. Nature. 2020;581:434-443.

Kircher M, Witten DM, Jain P, O'Roak BJ, Cooper GM, Shendure J. A general framework for estimating the relative pathogenicity of human genetic variants. Nature Genetics. 2014;46:310-315.

Landrum MJ, Lee JM, Benson M, et al. ClinVar: improving access to variant interpretations and supporting evidence. Nucleic Acids Research. 2018;46:D1062-D1067.

Liu X, Li C, Mou C, Dong Y, Tu Y. dbNSFP v4: a comprehensive database of transcript-specific functional predictions and annotations for human nonsynonymous and splice-site SNVs. Genome Medicine. 2020;12:103.

McLaren W, Gil L, Hunt SE, et al. The Ensembl Variant Effect Predictor. Genome Biology. 2016;17:122.

Prokhorenkova L, Gusev G, Vorobev A, Dorogush AV, Gulin A. CatBoost: unbiased boosting with categorical features. Advances in Neural Information Processing Systems. 2018.

Rehm HL, Berg JS, Brooks LD, et al. ClinGen, the Clinical Genome Resource. New England Journal of Medicine. 2015;372:2235-2242.

Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of ACMG and AMP. Genetics in Medicine. 2015;17:405-424.

Schymkowitz J, Borg J, Stricher F, Nys R, Rousseau F, Serrano L. The FoldX web server: an online force field. Nucleic Acids Research. 2005;33:W382-W388.

The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. Nucleic Acids Research. 2025;53:D609-D617.

Varadi M, Anyango S, Deshpande M, et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. Nucleic Acids Research. 2022;50:D439-D444.

Zhang X, Walsh R, Whiffin N, et al. Disease-specific variant pathogenicity prediction significantly improves variant interpretation in inherited cardiac conditions. Genetics in Medicine. 2021;23:69-79. doi:10.1038/s41436-020-00972-3.
