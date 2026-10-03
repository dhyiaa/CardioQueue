# Novelty positioning for the cardiogenetic CatBoost study

## Claims supported by the analysis

1. The study assembles 338 predictors across an 85,677-row cardiogenetic registry that includes multiple consequence classes, not missense variants alone.
2. Exact variants linked to HiRO, eMERGE, or CardioBoost are excluded before fitting and reported in separate source-held-out strata.
3. Public CardioBoost artifacts are executed on the same 250 matchable observations, with paired uncertainty clustered across 219 unique variants.
4. The model combines forced-binary evaluation with an explicit 0.1 to 0.9 deferral interval, post-training VUS prioritization, and a separate three-class analysis.
5. A ClinVar review-quality sensitivity fit removes one-star development records and tests the identical 430 external variants.
6. A registry-wide extension resolves overlapping AlphaFold v6 fragments, verifies reference residues against UniProt and PDB sequence, and adds DSSP and FreeSASA values for 3,382 development and 271 held-out variants.
7. The structure-enhanced sensitivity model shows a small ranking gain, with source heterogeneity and gene-clustered uncertainty reported. A separate outcome-blind FoldX analysis completed 4,137 unique calculations and used an eligibility-only control; generic DDG did not add independent value.
8. A held-out worklist simulation translates discrimination into finite review capacity. Recovering 80% of 255 P/LP variants required 208 score-ranked reviews versus a median 344 under 10,000 random orderings.
9. Worklist concentration is checked within each held-out source and against official CardioBoost scores on the same 250 matchable observations.
10. The planned release will join the trained models, frozen schemas, inference and analysis code, redistributable cleaned tables, environment, manifests, and checksums.

## Claims not supported

- Do not call this the first cardiogenetic machine-learning pathogenicity model. CardioBoost predates it.
- Do not call it the largest cardiogenetic model. CardioVar analyzed 117,582 ClinVar variants from 250 genes.
- Do not call it the first structure-aware pathogenicity model. AlphScore and SIGMA already used AlphaFold2-derived features.
- Do not state that FoldX improved the classifier. The pLDDT >=70 DDG model was worse than eligibility alone, the >=50 gain over baseline did not clearly exceed eligibility alone, and all-pLDDT inclusion added noise.
- Do not describe the small AlphaFold, DSSP, and FreeSASA sensitivity result as a universal structure benefit. Overall AUROC uncertainty crossed zero under gene-cluster resampling, coverage was missense-selected, and effects differed by source.
- Do not state that the primary model trained only on high-review ClinVar records. It retained 23,541 one-star training rows at weight 0.70.
- Do not describe the held-out sets as prospective or clinically independent validation.
- Do not describe random ordering as observed first-come performance. No arrival timestamps were available.
- Do not claim that 39.5% fewer variant reviews means 39.5% less analyst time. Review complexity and duration were not measured.
- Do not describe the score as emergency, patient-risk, or treatment triage. The study has variant labels, not patient outcomes or time-to-event endpoints.

## The central answer to "so what?"

The study tests whether a broad cardiogenetic score transfers across provenance boundaries and whether that ranking can organize a finite evidence-review queue. In the held-out simulation, reviewing 48.4% of records recovered 80% of known P/LP variants, while an unprioritized queue required a median 80.0% of records. Source-specific error, calibration drift, label quality, comparator eligibility, structural-feature attrition, and deferral remain visible. The next study can freeze this worklist, recalibrate it locally, and directly measure analyst time, evidence generation, classification change, and delayed actionable findings.

## Literature anchors

- Zhang et al. CardioBoost, *Genetics in Medicine* (2021): https://doi.org/10.1038/s41436-020-00972-3
- Al-Mahrami et al. CardioVar, *Bioinformatics Advances* (2026): https://doi.org/10.1093/bioadv/vbag135
- Schmidt et al. AlphScore, *Bioinformatics* (2023): https://doi.org/10.1093/bioinformatics/btad280
- Zhao et al. SIGMA, *Cell Genomics* (2024): https://pmc.ncbi.nlm.nih.gov/articles/PMC10831939/
