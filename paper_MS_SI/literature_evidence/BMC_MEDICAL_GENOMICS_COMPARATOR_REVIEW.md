# BMC Medical Genomics comparator review

## Purpose

Ten full-text articles from *BMC Medical Genomics* were reviewed to guide the manuscript's scope, structure, reporting, and claim boundaries. The set covers variant pathogenicity prediction, genomic prioritization, machine learning, workflow design, cardiogenetics, and VUS reclassification.

## Articles reviewed

| Article | Relevant design lesson |
|---|---|
| Cannon et al. (2023) | Separate performance on familiar variants from variants absent from common development resources. |
| Manshaei et al. (2022) | Distinguish gene validity, variant evidence, and case-level causality. |
| Dahary et al. (2019) | Show the complete annotation-to-review workflow and preserve traceable evidence. |
| Shi et al. (2019) | Define features and comparator methods precisely for a machine-learning predictor. |
| Rao et al. (2018) | Frame ranking around a concrete diagnostic bottleneck and evaluate named alternatives. |
| Anzar et al. (2019) | Compare an ensemble with individual tools on identical observations. |
| Zheng et al. (2020) | Report sensitivity, specificity, and discrimination across explicit validation conditions. |
| Schon et al. (2021) | Connect computational filtering to a clinical genomics workflow without overstating diagnosis. |
| Kim et al. (2024) | Anchor cardiogenetic interpretation in cohort counts, genes, domains, and classification categories. |
| Kwong et al. (2022) | Treat VUS as an evidence state and distinguish reclassification from model prioritization. |

## Decisions applied to the manuscript

The abstract reports cohort counts, discrimination, error-sensitive metrics, deferral, and uncertainty intervals. Figure 1 presents source integration, filtering, features, and held-out evaluation as one workflow. The Results separate frozen source-held-out performance from the HGVS-matchable CardioBoost comparison. VUS scores are called review-priority signals, never reclassifications. Methods name data releases, genome build, split precedence, weights, feature counts, hyperparameters, thresholds, software, and bootstrap settings. The Discussion states that source separation cannot remove shared ClinVar evidence or overlap in predictor training corpora.

## Full references

Cannon S, Williams M, Gunning AC, Wright CF. Evaluation of in silico pathogenicity prediction tools for the classification of small in-frame indels. *BMC Med Genomics*. 2023;16:36. doi:10.1186/s12920-023-01454-6.

Manshaei R, DeLong S, Andric V, et al. GeneTerpret: a customizable multilayer approach to genomic variant prioritization and interpretation. *BMC Med Genomics*. 2022;15:31. doi:10.1186/s12920-022-01166-3.

Dahary D, Golan Y, Mazor Y, et al. Genome analysis and knowledge-driven variant interpretation with TGex. *BMC Med Genomics*. 2019;12:200. doi:10.1186/s12920-019-0647-8.

Shi F, Yao Y, Bin Y, Zheng CH, Xia J. Computational identification of deleterious synonymous variants in human genomes using a feature-based approach. *BMC Med Genomics*. 2019;12(Suppl 1):12. doi:10.1186/s12920-018-0455-6.

Rao A, Vg S, Joseph T, Kotte S, Sivadasan N, Srinivasan R. Phenotype-driven gene prioritization for rare diseases using graph convolution on heterogeneous networks. *BMC Med Genomics*. 2018;11:57. doi:10.1186/s12920-018-0372-8.

Anzar I, Sverchkova A, Stratford R, Clancy T. NeoMutate: an ensemble machine learning framework for the prediction of somatic mutations in cancer. *BMC Med Genomics*. 2019;12:63. doi:10.1186/s12920-019-0508-5.

Zheng T, Zhu X, Zhang X, et al. A machine learning framework for genotyping the structural variations with copy number variant. *BMC Med Genomics*. 2020;13(Suppl 6):79. doi:10.1186/s12920-020-00733-w.

Schon U, Holzer A, Laner A, et al. HPO-driven virtual gene panel: a new efficient approach in molecular autopsy of sudden unexplained death. *BMC Med Genomics*. 2021;14:94. doi:10.1186/s12920-021-00946-7.

Kim OH, Kim J, Kim Y, et al. Exploring novel MYH7 gene variants using in silico analyses in Korean patients with cardiomyopathy. *BMC Med Genomics*. 2024;17:225. doi:10.1186/s12920-024-02000-8.

Kwong A, Ho CYS, Shin VY, Au CH, Chan TL, Ma ESK. How does re-classification of variants of unknown significance impact the management of patients at risk for hereditary breast cancer? *BMC Med Genomics*. 2022;15:122. doi:10.1186/s12920-022-01270-4.
