# Improving BMC Paper Progress

## Objective

Strengthen the BMC Medical Genomics submission around its clinically bounded contribution: a reproducible, multi-source cardiogenetic *review-ranking* model, not an autonomous pathogenicity classifier or patient-risk model. This log is the restart point for any subsequent Codex session.

## Status

Started: 2026-09-24

Primary manuscript: `paper_MS_SI/MS_BMC_MEDICAL_GENOMICS_FINAL.md`

Supplement: `paper_MS_SI/SI_BMC_MEDICAL_GENOMICS_FINAL.md`

Release package: `release/cardioqueue_v1.0/`

## Ten Review Workstreams

| Workstream | Role and reason | Status | Evidence / decision |
|---|---|---|---|
| 1. Manuscript architecture | Read the full BMC manuscript and identify gaps in clinical narrative, novelty, and claim boundaries. | Complete | The manuscript is 431 lines and already appropriately frames CardioQueue as review prioritization. Strengthen clinician-facing translation, not clinical-performance claims. |
| 2. Supplement audit | Read the full SI for concordance, methods, limitations, and reproducibility. | Complete | The SI is 1,116 lines. It records key limits: frozen baseline lacked numeric structure data in development; FoldX DDG did not enter the learned baseline or structure-enhanced model. |
| 3. Active modeling-code audit | Inspect active training, queue, structure-enhancement, and release inference code. | Complete | Primary pipeline separates external sources, uses a 300-feature public schema, and ships strict inference checks. Queue metrics are correctly score-ranking metrics, not analyst-time outcomes. |
| 4. Data/provenance audit | Review source, feature, split, and coverage documentation. | Complete | The analysis-ready registry contains 85,677 rows and 420 columns; the frozen model uses 300 public database-derived features. Exact held-out source variants were excluded before fitting. |
| 5. Structure-method audit | Evaluate AlphaFold/DSSP/FreeSASA/FoldX claims against actual coverage and sensitivity analyses. | Complete | Registry-wide AlphaFold v6 geometry, DSSP, and FreeSASA were added in a 351-feature sensitivity fit, with small source-held-out ranking gains. FoldX exists for 313 final rows and none in development, so it must remain feasibility/provenance evidence and future work, not a learned-performance claim. |
| 6. Clinical-workflow review | Answer what a clinician should take from this and how practice might change. | Complete | The appropriate clinical role is a laboratory/genetics-service worklist for evidence review, segregation/functional follow-up, and multidisciplinary preparation. It must not determine patient urgency, diagnosis, treatment, or cascade testing without expert interpretation. |
| 7. Cardiac-specific comparator search | Literature workstream covering CardioBoost, CardioClassifier, CardioVAI, GENESIS, CVD-PP, CardioVar, and gene-specific cardiac models. | Complete | The comparison must distinguish ACMG decision support, gene-specific missense predictors, disease-focused broad predictors, and CardioQueue's queue endpoint. |
| 8. General predictor comparator search | Literature workstream covering CADD, M-CAP, REVEL, ClinPred, MVP, VARITY, EVE, AlphaMissense, and related benchmark literature. | Complete | General predictors mainly predict variant effect/pathogenicity and usually focus on missense changes; they are not directly evaluated as evidence-review queue systems across multiple cardiogenetic consequence classes. |
| 9. Structure comparator search | Literature workstream covering AlphScore, SIGMA, FoldX/AlphaFold clinical use, and structure-informed evidence studies. | Complete | Structure-informed prediction is established outside cardiogenetics and in individual cardiac genes. Novelty must be claimed as integration and audit in a multi-gene cardiogenetic review workflow, not first use of structure or DDG. |
| 10. Clinical-standards and final consistency review | Compare framing with ACMG/AMP and cardiovascular guidance, then reconcile manuscript/SI language. | Complete | AHA and ESC guidance confirms expert multidisciplinary interpretation, phenotype and family evidence remain essential. Corrected the completion audit's obsolete 208 to the reproduced value of 207 reviews at 80% recovery. |

## Literature Workstreams (five search passes)

1. Cardiac-specific variant models and protein structure literature.
2. PubMed search for disease-specific models, clinical prioritization, and AlphaFold/FoldX use.
3. General and structure-informed models: VARITY, SIGMA, AlphScore, MVP, EVE, CardioClassifier, CardioVar.
4. Clinical guideline and broad predictor search: CADD, ClinPred, ACMG/AMP, AHA, ESC, InterVar.
5. Benchmark and predictor-threshold search: M-CAP, BayesDel, REVEL, structure-informed PP3/BP4, DS-MVP.

The comparator matrix was added as SI Section S21 and covers 22 primary/review papers. Sources should be rechecked before submission for journal formatting and newly released versions.

## Facts That Must Remain Accurate

- Registry: 85,677 analysis-ready GRCh38 rows; 33,923 B/LB, 42,361 VUS, 9,148 P/LP, plus 245 non-supervised/conflicting/missing-label rows.
- Frozen release: CatBoost with 300 public database-derived features; held-out sources were excluded before fitting.
- Primary external evaluation: 430 rows, 255 P/LP and 175 B/LB; AUROC 0.978 and AUPRC 0.985 in the manuscript's frozen v1.0 result.
- Queue result: 80% recovery requires 207 reviews, versus a random-order median 344; this measures retrospective score concentration only, not analyst time saved.
- CardioBoost comparison: 235 strict cDNA/protein-concordant observations, 211 unique variants; official artifacts passed positive controls before comparison.
- Structure enhancement: reference-validated AlphaFold v6 geometry, DSSP, and FreeSASA were available for 2,889 training, 493 validation, and 271 held-out rows. The enhanced model gained 0.0025 external AUROC, with gene-clustered 95% CI -0.00003 to 0.00502.
- FoldX: 313 final matrix DDG values, none in training/validation; do not say FoldX drove headline or enhanced-model performance.
- Scores are not calibrated pathogenicity probabilities in a new laboratory and do not replace ACMG/AMP, phenotype, family segregation, functional evidence, or clinical judgment.

## Planned Edits

- [x] Add clear clinician-facing workflow, boundaries, and prospective end points to the main manuscript.
- [x] Expand Introduction and Discussion with a calibrated comparator taxonomy and citations.
- [x] Add an SI comparator evidence matrix covering at least 20 papers/models.
- [x] Add an SI clinical-translation section distinguishing current use from future opportunities.
- [x] Harmonize the obsolete 208 review count in the completion audit to 207.
- [x] Re-run static consistency checks after edits.

## Changes Made This Session

- Updated `paper_MS_SI/MS_BMC_MEDICAL_GENOMICS_FINAL.md`: reframed the Introduction and Discussion around a clinician-facing evidence-review workflow; added cardiac, general-predictor, structure-model, AHA, and ESC context; stated the multi-consequence scope and its limits; clarified the FoldX boundary; and expanded references.
- Updated `paper_MS_SI/SI_BMC_MEDICAL_GENOMICS_FINAL.md`: added S20 Clinical Translation and Future Evaluation plus S21 Comparator Evidence Matrix (22 models/papers).
- Updated `release/cardioqueue_v1.0/README.md`: added intended workflow and explicit non-use boundaries.
- Updated `paper_MS_SI/COMPLETION_AUDIT.md`: corrected the 80% queue result from the obsolete 208 to 207 and updated table inventory wording.

## Verification Completed

- Markdown table-width scan of the main manuscript and SI: no mismatches.
- Unicode dash scan of the main manuscript and SI: none found.
- Reference-presence scan: every newly cited lead author has a corresponding reference entry in the main manuscript.
- Queue consistency scan: manuscript, SI, generated review-queue tables, and audit now use 207 score-ranked reviews for 80% P/LP recovery. Numeric occurrences of 208 that remain are unrelated values such as calibration slopes or row indices.
- Current word counts: main manuscript 8,573; SI 11,130; progress log 1,033.
- No model retraining was run: all performance claims retain their existing frozen artifacts and generated result tables. The active training, review-queue, structure-enhancement, and release-inference scripts were inspected for consistency with the written claims.

## Restart Instructions

1. Read this file, then `paper_MS_SI/MS_BMC_MEDICAL_GENOMICS_FINAL.md`, `paper_MS_SI/SI_BMC_MEDICAL_GENOMICS_FINAL.md`, and `paper_MS_SI/COMPLETION_AUDIT.md`.
2. Preserve the scientific boundaries above, especially FoldX coverage and retrospective-only claims.
3. Finish the planned edits, check citations against primary sources, then update this log with files changed and validation results.

## Second-Pass Clinical Utility and FoldX Audit (2026-09-24)

### Why CardioQueue is more than a reorderer

CardioQueue is a cardiogenetic variant-effect model whose most defensible present deployment is an evidence-review queue. It performs known B/LB versus P/LP discrimination, transfers to held-out sources, supports an explicit deferral zone, and retrospectively enriches later-resolved VUS near the top of a reanalysis worklist. The queue is the clinical action layer, not the whole scientific contribution. It must not be described as independently reclassifying a VUS because ACMG/AMP classification still requires phenotype, segregation, allelic, functional, and other evidence that the current model does not receive.

### FoldX denominator and corrected experiment

- The 85,677-row registry is not an appropriate FoldX denominator. FoldX stability delta-delta-G is defined for amino-acid substitutions and is not biologically applicable to synonymous, splice-only/intronic, truncating, UTR, and other non-substitution records.
- Registry-wide structure mapping identified 38,023 candidate amino-acid-substitution rows and mapped 38,021 across 57 UniProt accessions and 105 AlphaFold v6 fragments. This broad registry set includes VUS and rows outside the supervised primary analysis.
- In the frozen primary binary cohort, the outcome-blind corrected preflight selected 2,951 rows representing 2,860 unique substitutions across 82 fragments and 117 mutation batches. Split row counts are 2,365 train, 391 validation, 17 HiRO, 55 eMERGE, and 123 CardioBoost. Eligibility requires successful canonical mapping, reference-residue agreement, and local AlphaFold pLDDT >=70.
- The original matrix contained 313 numeric FoldX values, all outside development, so FoldX could not contribute to the learned baseline.
- A legacy FoldXPro attempt covered only 140/2,860 unique eligible calculations (123 training, 20 validation, and zero usable external values; duplicate substitutions can cross row splits). Tool/process incompatibilities made that run unsuitable for performance inference.
- The corrected FoldX 5.1 pipeline repairs each fragment once, batches BuildModel mutations, and passed an 11/11 smoke test. Script: `datasets/feature_sources/protein_structure/scripts/run_registry_foldx51_repaired_ddg.py`.
- A full corrected run was resumed in exec session `94404`. Sixty-seven repaired PDB files were cached. Eight large MYH7/RYR2 fragments remained under active repair for the full 20-minute compute window, after which the run was interrupted rather than silently excluding them. No final DDG table or performance result was produced. Rerunning the command will reuse completed repairs but retry interrupted fragments. Expected output: `datasets/feature_sources/protein_structure/interim/registry_foldx51_repaired_ddg.tsv`.

After completion, run in order:

1. `python3 datasets/feature_sources/protein_structure/scripts/build_structure_foldx_augmentation.py`
2. Train a DDG-only source-heldout model using `registry_foldx_only_features.tsv` and output `results/models/primary_binary_catboost_foldx_only_source_heldout`.
3. Train a structure-plus-DDG source-heldout model using `registry_structure_foldx_features.tsv` and output `results/models/primary_binary_catboost_structure_foldx_source_heldout`.
4. `python3 results/model_performance/structure_foldx_enhancement/evaluate_structure_foldx_enhancement.py`
5. Add SI Section S17 only after reviewing coverage, paired bootstrap intervals, and feature importance. Do not imply benefit if the result is null or adverse.

### Clinically useful present roles

- Scheduled VUS reanalysis watchlist: decide which old VUS records receive evidence review first.
- Evidence-acquisition triage: prioritize segregation studies, phenotype reconciliation, literature review, or functional testing, without treating the score as evidence itself.
- Inherited-cardiac multidisciplinary meeting preparation: surface technically plausible variants for expert discussion.
- Laboratory/database quality control: flag high-scoring B/LB or low-scoring P/LP records for transcript, allele, provenance, and annotation checks.
- Capacity planning: set transparent review budgets and a deferral band instead of silently forcing every record into a binary action.
- Mechanism-aware research queue: combine sequence, regional constraint, structural context, and eligible FoldX stability effects to nominate variants for experimental follow-up.

Potential future roles requiring new data and prospective validation include patient-specific phenotype-gated ranking, automated report-update triggers, gene-panel curation surveillance, active learning from adjudicated reviews, and trial/registry recruitment support. None is a current validated use.

### High-value, relatively easy additions

1. Finish the corrected all-eligible FoldX comparison. This directly resolves a visible incompleteness and provides a clean DDG-only versus geometry-only versus combined ablation.
2. Produce an external error atlas for the 9 false negatives and 16 false positives at the frozen 0.5 threshold, stratified by source, gene, consequence, and feature coverage. This is more clinically informative than another global AUROC.
3. Run a focused no-precomputed-meta-predictor ablation excluding CADD, REVEL, AlphaMissense, SIFT, PolyPhen, MetaLR, ESM1b, and similar scores. This tests whether the system adds value beyond repackaging existing predictors.
4. Test rank stability over multiple seeds or bootstrap refits, reporting top-k overlap and rank correlation. A clinical queue must be reproducible, not merely accurate on average.
5. Generate concise SHAP/evidence cards for representative high-priority VUS and external errors, with a warning that explanations describe model behavior rather than ACMG evidence.

Implementation started in this pass:

- Added and ran `results/model_performance/external_error_atlas/analyze_external_errors.py` against the promoted 300-feature `cardioqueue_v1_source_heldout` artifact. At the frozen 0.5 threshold, the 430 held-out rows contain 9 false negatives and 16 false positives: CardioBoost 6/8, eMERGE 2/4, and HiRO 1/4, respectively. Eleven false positives are VEP missense variants. Three false negatives are annotated as intergenic and one as intronic, making transcript/annotation reconciliation a specific audit target. Seventeen of 25 discordances have structure features; 12 had legacy numeric FoldX values. Treat these as model/label discordances, not proven source-label errors.
- Added a focused `no_precomputed_effect_predictors` ablation and an `--only` option to `results/model_performance/leakage_ablation_audit/run_primary_binary_leakage_ablation.py`. The ablation removes SIFT, PolyPhen-2, REVEL, MetaLR, FATHMM-XF, CADD, AlphaMissense, and ESM-1b fields.
- Extended `datasets/modeling/scripts/train_primary_binary_catboost.py` with repeatable `--exclude-regex` and `--thread-count` options so the same ablation can be run directly against the promoted 300-feature public model recipe. A two-thread test was interrupted to leave CPU capacity for FoldX; the subsequent full-machine run completed in `results/models/cardioqueue_v1_no_precomputed_effect_predictors_source_heldout`.
- Completed the full public-schema no-precomputed-effect-predictor fit after FoldX was stopped. It retained 242 predictors. On 430 held-out rows, AUROC was 0.976863 versus 0.978006 for CardioQueue v1.0; AUPRC was 0.984762 versus 0.984757; Brier score was 0.053963 versus 0.049953; and the 0.5 confusion matrix was 246 TP, 18 FP, 157 TN, and 9 FN. Paired 5,000-bootstrap ablation-minus-baseline differences were -0.001143 AUROC (95% CI -0.005827 to 0.003989), 0.000004 AUPRC (-0.003613 to 0.004321), and 0.004010 Brier (-0.004957 to 0.013436). This argues against simple repackaging of those eight predictor families while not eliminating broader public-evidence circularity.
- Added `results/model_performance/no_precomputed_effect_predictor_ablation/evaluate_ablation.py`, machine-readable result tables, and SI Tables S10b/S10c. Added an explicit S17 documenting the FoldX preflight and unfinished compute status instead of implying an uncompleted sensitivity result.

### Publishability judgment

The work is plausibly publishable in *BMC Medical Genomics* as a retrospective methods/resource paper if it is framed as multi-source, multi-consequence cardiogenetic evidence-review prioritization with reproducible artifacts. Its strongest differentiators are the 85,677-row registry, explicit source-heldout evaluation, direct CardioBoost artifact comparison, temporal VUS analysis, abstention/worklist endpoint, deployable release checks, and outcome-blind structure/FoldX audit. It is not yet a clinical-utility or autonomous-classification study. Principal reviewer risks remain label/predictor circularity, consensus selection of easier cases, limited external sample size, retrospective temporal leakage, weak calibration, and absence of measured analyst time, prospective yield, or patient outcomes.

### Second-pass files changed and verification

- Revised `paper_MS_SI/MS_BMC_MEDICAL_GENOMICS_FINAL.md` with evidence-acquisition triage, VUS surveillance, quality-control and future-use language; the 242-feature predictor ablation and paired intervals; and explicit wording that corrected DDG integration remains incomplete.
- Revised `paper_MS_SI/SI_BMC_MEDICAL_GENOMICS_FINAL.md` with S10 public-schema ablation results and a new S17 corrected FoldX preflight/compute-status section.
- Added `paper_MS_SI/tables/Table_S10b_no_effect_predictor_ablation.tsv`, `Table_S10c_no_effect_predictor_paired_bootstrap.tsv`, and `Table_S17_foldx_preflight_eligibility.tsv`.
- Added the external error atlas and paired effect-predictor ablation evaluation scripts and outputs under `results/model_performance/`.
- Trained `results/models/cardioqueue_v1_no_precomputed_effect_predictors_source_heldout`; this supersedes the earlier session note that no retraining had been run.
- `py_compile` passed for all scripts changed in this pass. Markdown table-width checks found no mismatches. Current word counts are 8,880 for the main manuscript and 11,603 for the SI. No FoldX or model-training process remains running.

## GeneTerpret Full-Paper Transfer Review (2026-09-24)

Source read in full: Manshaei et al., *GeneTerpret: a customizable multilayer approach to genomic variant prioritization and interpretation*, BMC Medical Genomics 2022;15:31. The supplied text included the complete Background, Implementation, Results, Discussion, Conclusions, references, and supplementary-file descriptions.

### Transferable strengths adopted

- GeneTerpret explicitly separates gene-disease validity, variant pathogenicity, and patient-level causality. CardioQueue now uses that same conceptual boundary in the Introduction, Methods, Discussion, and a redesigned Figure 1: CardioQueue models molecular variant effect/review priority, curated evidence constrains gene validity, and patient phenotype/family evidence remains outside the model.
- GeneTerpret validates modules against different reference sources instead of using one metric for the entire platform. CardioQueue's Methods now describe a five-question validation ladder: discrimination, source transfer, fixed-capacity recovery, executable comparator concordance, and historical VUS enrichment, with protein/leakage stress tests kept separate.
- GeneTerpret reports whether expert-selected variants appear in the top 10 or 50. CardioQueue already has a stronger fixed-budget framework; Figure 5 is now a two-panel top-k recovery figure showing both the 430-row source-heldout queue and the 457-row temporal VUS queue with clinically readable annotations.
- GeneTerpret discusses discordant variants individually. CardioQueue's 25 external model/label discordances are now summarized in the main Results and detailed in SI S6, with variant-level Table S6b and source-level Table S6c.
- GeneTerpret attempts to measure analyst time directly. SI S20 now specifies a stronger randomized, counterbalanced, multi-reviewer crossover study with analyst minutes as the primary endpoint, evidence/reclassification yield as effectiveness outcomes, delayed actionable findings as safety outcomes, mixed-effects analysis, frozen model/protocol, and blinded adjudication.

### Weaknesses not copied

- GeneTerpret's 97.2% gene-validity result is essentially recovery of 695/715 positive ClinGen associations and does not establish specificity against unsupported genes.
- Its 83.5% DECIPHER statement does not provide the full discrimination, calibration, or uncertainty reporting expected for CardioQueue.
- Patient-level validation is small: 10 families and 20 unrelated probands, with only five reported variants of interest in the latter cohort.
- The analyst-time assessment was internal and apparently non-randomized, without a crossover design, blinded outcome adjudication, uncertainty intervals, or safety endpoint.
- User-adjustable weights and thresholds improve flexibility but complicate reproducibility. CardioQueue retains a frozen schema, score rule, and provenance manifests for evaluation.

### Files changed in this transfer pass

- Added `results/manuscript_figures_tables/generate_cardioqueue_workflow.py` and generated `paper_MS_SI/figures/Figure_1_cardioqueue_workflow.png` while retaining the former filename as a compatibility alias.
- Extended `results/model_performance/review_queue_simulation/evaluate_review_queue.py` so Figure 5 presents source-heldout and temporal-VUS recovery side by side.
- Extended `results/model_performance/external_error_atlas/analyze_external_errors.py` to emit Tables S6b and S6c.
- Revised the main manuscript's conceptual framing, validation rationale, discordance Results, comparator argument, and Figure 1/Figure 5 legends.
- Revised SI S6, S20, the comparator matrix, and the figure inventory.

Transfer-pass verification: all three changed Python scripts passed `py_compile`; both manuscripts have zero Markdown table-width mismatches and zero replacement characters. Current word counts are 9,108 for the main manuscript and 11,933 for the SI. Figure 1 rendered at 3,930 x 2,129 pixels and Figure 5 at 3,329 x 1,409 pixels; both were visually inspected for readable text, complete curves, and non-overlapping annotations.

## ClinVar Challenge-Cohort Design Audit (2026-09-25)

### Decision

The proposed additional 570 ClinVar variants can strengthen robustness testing, but they must not be described as an expansion of the independent/source-held-out cohort. The defensible hierarchy is:

1. Primary source-held-out evaluation: the existing 430 HiRO, eMERGE, and CardioBoost rows.
2. Secondary same-ecosystem ClinVar challenge: 570 prespecified variants removed before a complete retraining.
3. Exploratory mixed-provenance benchmark: 1,000 rows, always accompanied by the two separate cohort results.

A larger round number does not compensate for shared ClinVar evidence. The added cohort is valuable for label quality, gene/consequence breadth, calibration stress, and queue robustness, not for a claim of 1,000 independent external variants.

### Local eligibility result

Added and ran `results/model_performance/clinvar_challenge_cohort/build_eligibility_audit.py`. It joins the rescue-aware modeling table to the curated ClinVar phenotype fields and requires an existing train/validation assignment, binary consensus label, primary-model inclusion, at least two stars, a cardiac phenotype keyword, no obvious RASopathy keyword, and no HiRO/eMERGE/CardioBoost overlap.

- Eligible candidates: 12,341.
- Two-star multiple-submitter/no-conflict candidates: 12,228 (9,985 B/LB; 2,243 P/LP).
- Three-star expert-panel candidates after cardiac-VCEP or phenotype screening: 113 (53 B/LB; 60 P/LP).
- Local ClinGen Evidence Repository matching directly verifies 92 eligible, current-label-concordant cardiac-VCEP rows: 83 Cardiomyopathy VCEP and 9 Potassium Channel Arrhythmia VCEP. The broader local ERepo contains 118 binary cardiac-VCEP assertions, but external overlap, current eligibility, phenotype handling, and label concordance reduce the usable challenge seed.
- A provisional addition of 325 B/LB and 245 P/LP would make the descriptive mixed 1,000-row benchmark exactly balanced at 500/500. After including the 92 verified cardiac-VCEP rows, this would require 273 B/LB and 205 P/LP condition-filtered rows. This is a planning target, not a frozen sample.

Generated artifacts:

- `results/model_performance/clinvar_challenge_cohort/clinvar_cardiac_challenge_eligible.tsv`
- `results/model_performance/clinvar_challenge_cohort/clinvar_cardiac_challenge_eligibility.summary.json`
- `results/model_performance/clinvar_challenge_cohort/CLINVAR_CARDIAC_CHALLENGE_AUDIT.md`

### Required before selection and retraining

- Verify disease context and expert-panel affiliation at the ClinVar assertion/ClinGen Evidence Repository level. A cardiac-panel gene and a three-star aggregate record do not by themselves prove cardiogenetics-specific adjudication.
- Freeze a dated 570-row manifest before inspecting model predictions. Record variant ID, label, condition, review status, panel/submitter provenance, gene, consequence, selection stratum, inclusion probability, and random seed.
- Include all eligible condition-verified cardiac expert-panel variants where possible, then draw the two-star remainder with prespecified label, gene, and consequence stratification. Cap dominance by any one gene and avoid a benign-synonymous-heavy challenge set.
- Add a distinct `challenge_clinvar` split value. Never name it `external_clinvar`, because multiple scripts interpret `external_*` as independent source holding.
- Remove all selected alleles from training, validation, tuning, calibration, and feature selection before retraining. Directly scoring a model that already trained on these rows is invalid.
- Report unweighted constructed-cohort metrics and design-weighted prevalence-sensitive metrics. PPV, NPV, accuracy, and AUPRC from a balanced sample do not estimate routine clinical prevalence.
- Report variant-level stratified and gene-clustered bootstrap intervals; retain source-stratified intervals for the original 430.

### Rerun dependency chain after the manifest is frozen

1. Rebuild the rescue-aware source-heldout split with the challenge quarantine.
2. Retrain the primary, high-star, source-weight-only, no-precomputed-predictor, structure, and any completed DDG sensitivity models with unchanged recipes.
3. Regenerate discrimination, calibration, deferral, review-queue, subgroup, error-atlas, and paired-ablation results separately for the 430 and 570 cohorts.
4. Regenerate the CardioBoost matched comparison because primary-model predictions will change after retraining.
5. Regenerate affected figures/tables and then update manuscript/SI numbers. Do not pool first and report subcohorts second.

### Current status

The eligibility audit passed `py_compile` and reproduced the counts above. No split was changed, no 570-row sample was frozen, no model was retrained, and no manuscript performance number was changed in this pass. This pause is intentional: assertion-level provenance must be verified before selection to avoid converting a methodological improvement into a mislabeled same-source test.

Pipeline preparation was completed without changing the frozen outputs. `datasets/modeling/scripts/build_training_slices_and_splits.py` now accepts an optional `--clinvar-challenge-manifest`, validates unique and eligible IDs, rejects overlap with HiRO/eMERGE/CardioBoost, and creates `split_source_heldout_clinvar_challenge` with the distinct value `challenge_clinvar`. A default-mode regression test reproduced all 85,677 existing split assignments exactly. A 10-row smoke manifest produced 10 challenge rows, retained all 430 external rows, and rebuilt train/validation as 36,167/6,383. The production option must not be used until the final 570-row manifest is frozen.

## Public External-Dataset Search (2026-09-25)

### Search organization

Ten independent search assignments were launched across cardiomyopathy registries, inherited arrhythmia cohorts, ClinGen VCEPs, locus-specific databases, clinical decision-support datasets, hospital cohorts, population-screening studies, sudden-death/pediatric cohorts, public data repositories, and hostile-review methodology. The subagent service reached its usage limit during execution: one ClinGen workstream completed, while the remaining searches were completed directly using primary articles, supplements, repositories, and local file inspection. Do not imply that ten agent reports completed successfully.

### Main result

The strongest route is not a 570-row ClinVar sample. Two public clinical-curation sources can produce a much larger and more defensible exact-variant-held-out evaluation:

1. SHaRe HCM 2026Q1: 822 binary variants, including 561 P/LP and 261 B/LB across 26 genes, with GRCh38 coordinate IDs and contemporary HCM-focused curation.
2. Fernandez-Falgueras et al. PLOS ONE inherited-cardiovascular reinterpretation cohort: 197 binary records in the current CardioQueue panel, including 103 P/LP and 94 B/LB across 38 genes, with applied ACMG criteria.

The two sources share 43 gene/cDNA variants, and all 43 have concordant binary labels (27 P/LP and 16 B/LB). The existing 430 should remain separately reported. After preserving existing-source precedence, SHaRe alone could add approximately 775 rows, taking the exact-variant-held-out evaluation above 1,200 before adding nonduplicate PLOS records. Final counts require coordinate normalization and full deduplication.

### SHaRe audit

- Public repository: `https://github.com/ImperialCardioGenetics/share_browser`.
- Source file inspected: `HCM_Gen_Variant_Annotation_2026Q1.csv`, 2,637 rows and 36 fields.
- Final classifications: 561 P/LP, 261 B/LB, and 837 VUS; other rows are non-HCM VUS or missing.
- Exact local-coordinate overlap among binary rows: 505/822. Of these, 351 are currently training, 68 validation, 42 external CardioBoost, 5 external HiRO, and 39 not eligible. Therefore, 419 current development variants must be quarantined before a complete retraining.
- Forty-seven binary SHaRe rows overlap an existing external stratum. Preserve current source precedence rather than counting them twice.
- The full binary SHaRe set contains 115 variants without a mapped ClinVar class (95 P/LP, 20 B/LB). Report this as a difficult ClinVar-absent sensitivity subset.
- Major limitation: the repository states "All rights reserved" and has no reuse license. Written permission from the SHaRe investigators is required before publication analysis or redistribution of a derived manifest. Contacts in the repository include Sophie Hespe for genetic data and Pantazis Theotokis for the browser.
- CardioBoost used older SHaRe cardiomyopathy data, so SHaRe and CardioBoost are not independent evidence ecosystems even after exact deduplication.

### PLOS reinterpretation-cohort audit

- Article DOI: `10.1371/journal.pone.0297914`; article license CC BY.
- Official supplement downloaded to `datasets/external_validation_candidates/plos_cardiovascular_reinterpretation/pone.0297914.s001.xlsx`.
- SHA-256: `79eae13f642c14ddadbef5913728745a32315f0164d1f3b6be8ded8e3dd5fbce`.
- Spreadsheet: 1,261 listed rows with original class, updated class, and ACMG criteria. Updated classes comprise 176 B, 128 LB, 72 LP, 58 P, 826 VUS, and one blank.
- Current CardioQueue panel: 197 binary records, 94 B/LB and 103 P/LP across 38 genes.
- Limitation: only gene/cDNA and protein HGVS are supplied. Transcript-aware conversion to GRCh38, reference validation, and manual ambiguity review are mandatory before freezing a manifest.

### Other source decisions

- Live ClinGen ERepo cardiac VCEPs currently provide 149 binary assertions: 129 cardiomyopathy and 20 potassium-channel arrhythmia. Use as an expert-curated challenge, not independent external validation, because ClinGen curations are submitted to ClinVar.
- CardioClassifier provides a small expert-curated/rule-concordance resource but overlaps ClinVar, ACGV, and the Imperial/SHaRe ecosystem.
- BBI-CVD has 4,770 unique clinical variants but requires investigator access and is not cardiology-specific. Requesting cardiac-gene rows is a strong future option.
- Brugada, early-atrial-fibrillation, DCM, arrhythmogenic-cardiomyopathy, and pediatric cohorts mostly publish P/LP and VUS only. Use as P/LP-sensitivity or VUS-ranking stress tests, never as binary AUROC cohorts paired with benign variants from another source.
- Cardiac VariantFX, ACGV, LOVD, CardioVar, population cohorts, and functional assay repositories fail at least one primary requirement: five-tier same-source truth, independence, bulk access, reuse permission, or balanced labels.

### New documentation

- Added `docs/public_external_validation_dataset_audit.md` with the candidate hierarchy, rejection reasons, links, reviewer risks, and required protocol.
- Added `datasets/external_validation_candidates/plos_cardiovascular_reinterpretation/README.md` with source, license, checksum, counts, and mapping limitations.
- Added `datasets/external_validation_candidates/audit_public_external_candidates.py`. With the public SHaRe 2026Q1 CSV it reproduced 822 SHaRe binary rows, 197 PLOS panel-binary rows, 43 preliminary gene/cDNA overlaps, zero binary-label discordances, and a preliminary 976-row nonduplicate union. Output: `results/model_performance/public_external_dataset_audit.summary.json`.
- Added `docs/share_data_permission_request_draft.md` so the rights question is resolved before publication analysis or redistribution.

### Next actions before modeling

1. Obtain written SHaRe permission.
2. Build a transcript-aware PLOS coordinate resolver and produce an ambiguity report; do not infer transcripts silently.
3. Freeze source versions and a precedence-aware union manifest.
4. Remove all normalized equivalent alleles from development and retrain the unchanged recipe once.
5. Report current 430, SHaRe, PLOS, ClinGen challenge, and pooled results separately. Use source-stratified and gene-clustered bootstrap intervals plus a no-precomputed-score ablation.

No new model was trained and no manuscript performance claim was changed during this search pass.

## Live ClinGen Cardiac VCEP Search (2026-09-25)

The live ClinGen Evidence Repository summary API was queried by VCEP affiliation ID. Counts below are current published variant-condition assertions returned by the API; B/LB and P/LP are the binary classes relevant to CardioQueue.

| VCEP | Affiliation | Status | Published assertions | Binary B/LB or P/LP | Published genes |
|---|---:|---|---:|---:|---|
| Cardiomyopathy | 50002 | Step 4 approved Apr. 2017 | 216 | 129 (50 B, 9 LB, 35 LP, 35 P) | MYH7, MYBPC3, ACTC1, MYL2, MYL3, TNNI3, TNNT2, TPM1 |
| Potassium Channel Arrhythmia | 50108 | Step 4 approved May 2025 | 37 | 20 (4 B, 4 LB, 8 LP, 4 P) | KCNQ1 |
| RYR2 | 50169 | Step 1 completed Aug. 2025 | 0 | 0 | None yet |
| Sodium and Calcium Channel Arrhythmia | 50160 | Step 1 completed Mar. 2025 | 0 | 0 | None yet; planned initial focus is SCN5A |
| Desmosomal Cardiomyopathy | 50148 | Step 1 completed May 2023 | 0 | 0 | None yet; planned scope is PKP2 and DSP |
| Congenital Heart Disease | 50130 | Step 2 completed Oct. 2025; pilot in progress | 0 | 0 | None yet; planned scope is GATA4, MYH6, NKX2-5, NR2F2 |

The two publishing panels therefore provide only 149 currently downloadable binary assertions across these six requested VCEPs. API results had one unique CAID per binary assertion. The ERepo homepage reported a temporary VCI-to-ERepo update interruption from 2026-09-23 through 2026-09-30, so the counts are the latest publicly obtainable snapshot but may omit newly approved records awaiting synchronization.

Primary access points:

- Bulk CSV/TSV: `https://search.clinicalgenome.org/kb/downloads`
- ERepo API documentation: `https://genboree.org/gitlab/clingen/erepo/erepo-records-table/-/wikis/Evidence-Repository-API`
- Query template: `https://erepo.clinicalgenome.org/evrepo/api/summary/classifications?columns=affiliationId&values=50002&matchTypes=exact&pgSize=500&pg=1&matchMode=and`

External-validation conclusion: these are high-quality, condition-specific expert assertions, but they are not an independent evidence ecosystem. ClinGen explicitly reviews ClinVar data and submits approved VCEP classifications back to ClinVar. CardioQueue development labels came from ClinVar, so exact overlap and earlier ClinVar appearances must be audited. Existing or retrospectively withheld VCEP variants should be described as a ClinGen expert-curated/same-ecosystem challenge, not external validation. A defensible temporal evaluation would require a previously locked model and variants entirely absent from the development ClinVar snapshot, with no post-cutoff ClinVar-derived feature updates. The current six VCEPs cannot supply approximately 500 external variants.

## Expanded 1,233-Variant Evaluation Completed (2026-09-25)

### Objective and agent continuity

This pass executed the previously recommended design rather than sampling another 500 ClinVar variants. Earlier agent/search workstreams identified and audited SHaRe, PLOS, ClinGen, and alternative cardiac datasets. Implementation, annotation, retraining, evaluation, and manuscript synchronization were completed in the main workspace so every decision traces to one frozen set of files. Do not imply ten completed agent reports: the subagent service previously reached its usage limit, as documented above.

### External-source review

The earlier 430-row benchmark had three different externality levels: HiRO (69) is an internal unpublished source-held-out dataset; eMERGE (176) is a public institutionally external source; and CardioBoost (185) is public and external but connected to older SHaRe and the broader ClinVar/literature ecosystem. Strict GRCh38 reference validation retained HiRO and eMERGE but only 174 CardioBoost alleles; 11 CardioBoost coordinates failed reference-base validation. The old queue and comparator results remain historical frozen analyses.

New sources:

- SHaRe HCM 2026Q1: 822 starting binary rows; 757 reference-valid in-panel variants (535 P/LP, 222 B/LB). One cross-source disagreement was excluded from source evaluation, leaving 756 observations. Institutionally external, but not evidence-ecosystem independent. Written reuse permission remains required.
- Fernandez-Falgueras et al. PLOS ONE 2024: 197 panel-binary records; 154 mapped uniquely and reference-valid (74 P/LP, 80 B/LB). Forty-three ambiguous records were not inferred. PLOS supplied labels; ClinVar was used only to resolve coordinates. CC BY.

One conflict was excluded: `11-47349804-C-G` in `MYBPC3`, called B/LB by SHaRe and P/LP by legacy CardioBoost.

Final benchmark: 1,233 unique GRCh38 alleles, 801 P/LP and 432 B/LB across 42 genes. Source-specific analyses contain 1,328 observations because alleles can occur in multiple sources. Of the unique alleles, 942 already existed in the matrix and 291 required annotation.

Primary manifests:

- `datasets/external_validation_candidates/expanded_external_2026/external_unique_variant_manifest.tsv`
- `datasets/external_validation_candidates/expanded_external_2026/external_source_observations.tsv`
- `datasets/external_validation_candidates/expanded_external_2026/external_label_conflicts.tsv`
- `datasets/external_validation_candidates/expanded_external_2026/new_variants_requiring_annotation.tsv`
- `datasets/external_validation_candidates/expanded_external_2026/build_summary.json`

### Cohort builder and annotation

Added `datasets/external_validation_candidates/build_expanded_external_benchmark.py`. It normalizes source rows, validates exact alleles against GRCh38, applies the panel, preserves source assertions, excludes disagreements, deduplicates, removes exact benchmark alleles from development, and regenerates a seeded 85:15 split. The 85,968-row matrix contained 35,796 training, 6,317 validation, and 1,233 `external_expanded` rows, with zero exact benchmark alleles in development. Seed: 20260925.

Added `annotate_vep_rest.py` and `merge_expanded_external_features.py`. Local VEP 116 crashed because the existing environment mixed ARM and x86 components. A clean ARM64 environment at `~/opt/miniforge3/envs/vep-116-arm64` also segfaulted during Perl compilation, so annotation moved transparently to the official Ensembl REST endpoint. All 291 new alleles received VEP output.

New-allele coverage: dbNSFP 195/291, AlphaMissense 170/291, SpliceAI 224/291, reviewed protein mapping 170/291, and observed gnomAD 149/291. Structure processing validated all 170 mapped substitutions and completed DSSP/FreeSASA for all 35 required AlphaFold fragments. FoldX was not run for the 291 new variants; values remain explicitly unavailable.

Final matrix: `datasets/external_validation_candidates/expanded_external_2026/modeling_table_expanded_external_annotated.tsv`, 85,968 rows by 420 columns. All 300 model fields retained original typing (203 categorical, 97 numeric). Benchmark coverage: VEP 100.0%, SpliceAI 83.9%, dbNSFP 73.6%, gnomAD 72.3%, AlphaMissense 55.7%, protein mapping 33.9%, AlphaFold pLDDT 33.1%, and existing FoldX 13.5%.

### Models and results

The first diagnostic fit at `results/models/cardioqueue_v1_expanded_external_2026` had 338 features because HiRO aggregate fields were admitted. It is discarded and must not be reported.

Accepted model: `results/models/cardioqueue_v1_expanded_external_2026_public300`. It stopped after 314 trees. Across 1,233 alleles: AUROC 0.984094 (gene-clustered 95% CI, 0.977960-0.988231), AUPRC 0.991677 (0.978968-0.994651), Brier 0.051293 (0.041224-0.071000), sensitivity 0.968789 (0.957539-0.985041), and specificity 0.865741 (0.811897-0.908821). Confusion matrix: TP 776, FP 58, TN 374, FN 25. Calibration intercept -0.7298 and slope 0.8515.

Source AUROCs: PLOS 0.9813, SHaRe 0.9858, strict CardioBoost 0.9632, eMERGE 0.9854, and HiRO 0.9947. The 291 newly annotated alleles reached AUROC 0.9858 and AUPRC 0.9900.

Evaluation files are under `results/model_performance/expanded_external_2026/`. The evaluator is `datasets/external_validation_candidates/evaluate_expanded_external_benchmark.py`.

### Sensitivity analyses

Added `build_molecular_quarantine_sensitivity.py`. It removed 437 additional development rows: 383 shared a genomic site, 87 shared an exact VEP protein substitution, and 33 met both criteria. The stricter matrix contained 35,428 training and 6,248 validation rows. Model `cardioqueue_v1_expanded_external_2026_molecular_quarantine_public300` reached AUROC 0.982487, AUPRC 0.991082, Brier 0.052362, sensitivity 0.961298, and specificity 0.865741.

The 242-feature no-precomputed-predictor model excluded 58 SIFT, PolyPhen, REVEL, MetaLR, FATHMM, CADD, AlphaMissense, ESM1b, and HiRO-family fields. Model `cardioqueue_v1_expanded_external_2026_no_precomputed_effect_predictors` reached AUROC 0.987261, AUPRC 0.993327, Brier 0.039975, sensitivity 0.967541, and specificity 0.898148.

Added `compare_expanded_models.py`. Five thousand paired gene-clustered samples gave AUROC difference +0.003167 (95% CI, +0.000393 to +0.006011), AUPRC difference +0.001650 (-0.000142 to +0.005264), and Brier difference -0.011318 (-0.020697 to -0.006631). Established scores did not drive the expanded result, although public-evidence circularity remains possible.

### Manuscript and handoff

Added `docs/expanded_external_validation_review.md`. Updated the main manuscript abstract, methods, primary Results, leakage analyses, Discussion, limitations, conclusion, and availability. Added SI section S22 and four machine-readable S22 tables for performance, clustered intervals, feature coverage, and paired ablation results.

Before submission:

1. Obtain written SHaRe permission; do not redistribute row-level SHaRe data before permission.
2. Use “retrospective multi-source source-held-out evaluation,” not “independent external clinical validation.”
3. Do not claim generalization to all cardiogenetic mutations; the benchmark covers 42 genes and not every mechanism, ancestry, CNV, or structural variant.
4. Do not claim FoldX improves performance; new-variant FoldX was not run and development DDG coverage remains insufficient.
5. Keep the earlier 430-row queue analysis as a legacy operational experiment; queue budgets were not rerun on the new 1,233-row cohort in this pass.

### Clarification: ClinVar overlap versus model holdout

The accepted expanded model has **zero exact benchmark alleles** in its training or early-stopping validation partitions. This required moving 381 final benchmark alleles out of the earlier ClinVar-derived training split and 77 out of the earlier validation split, then retraining from scratch. Of the 1,233 benchmark alleles, 942 were present somewhere in the original matrix and 291 were newly appended. Therefore use the claim “no exact evaluation allele informed fitting or early stopping.” Do not claim that all evaluation variants were absent from ClinVar or had never appeared in a development table. At source level, 419/756 SHaRe and 66/154 PLOS observations had occupied the prior train/validation partitions; legacy HiRO, eMERGE, and CardioBoost variants were already source-held-out. Public predictor training corpora and classification evidence can still overlap conceptually or directly.

## Final Five-Source Model-Held-Out Design (2026-09-25)

### Agent and direct-audit organization

The user requested ten agents. The subagent service had a hard thread limit from prior work: six focused audits completed, covering SHaRe provenance, PLOS reinterpretation, CardioBoost/SHaRe relationships, eMERGE provenance, HiRO terminology, and exact local cohort accounting. Four remaining workstreams were completed directly using the primary papers, repository source, local ClinVar snapshot, manifests, and manuscript audit. Do not claim ten successfully completed agents.

### Final eligibility rule

The primary model-held-out benchmark now contains:

- all 69 strict HiRO observations because they were source-held-out from model development;
- all 176 strict eMERGE observations;
- all 173 strict CardioBoost observations;
- 337/756 SHaRe observations whose exact alleles were absent from the original training and validation partitions;
- 88/154 PLOS observations whose exact alleles were absent from the original training and validation partitions.

HiRO is model-external but institutionally internal. The other four sources are institutionally external. Exact cross-source alleles were counted once and the discordant `MYBPC3` allele `11-47349804-C-G` remained excluded.

The five sources contributed 843 observations. Sixty-four variants occurred in more than one retained source combination; after exact-allele deduplication, the final benchmark contains 775 variants, 436 P/LP and 339 B/LB, across 40 genes. Exact overlap with both original and accepted model training/validation partitions is zero.

### Classification-provenance audit

PLOS is an external institutional reinterpretation, not a ClinVar-label export. The Fernandez-Falgueras team reinterpreted variants in 2022 using ACMG/AMP criteria, updated criterion specifications, gnomAD, computational evidence, literature, functional and segregation evidence, affected cases, and institutional patient data. The paper does not name ClinVar as the source of updated labels. ClinVar was used in this project only to resolve PLOS gene+cDNA descriptions to GRCh38 coordinates. All 88 retained PLOS alleles consequently have a current ClinVar coordinate record, but only 10 have a concordant binary aggregate; 78 are conflicting, uncertain, or otherwise nonbinary.

SHaRe is a mixed-provenance, SHaRe-led expert reclassification source rather than a ClinVar copy. Baseline calls came from participating clinical laboratories, followed by systematic SHaRe re-evaluation. The 2026Q1 release stores SHaRe and ClinVar classifications in separate fields and flags disagreement. Among 337 retained observations, 299 are marked manually curated and two manually curated/check. Seventy-seven are absent from the current local ClinVar release; 210 have a discordant or nonbinary aggregate; only 50 have a concordant binary aggregate. ClinVar was one possible evidence source alongside gnomAD, literature, functional evidence, segregation, and internal SHaRe cases, so complete evidence independence is not claimed.

CardioBoost's historical 129-variant SHaRe replication cohort was not its training or primary holdout. The current 2026Q1 SHaRe release nevertheless overlaps the local strict CardioBoost subset at 42 exact alleles; 41 labels agree and one discordance was excluded. Source-specific results remain separate and the combined benchmark is deduplicated.

eMERGE classifications were produced through sequencing-center ACMG/AMP interpretation, ClinVar-informed manual review for part of the gene set, and electrophysiology for a small subset. It is external to model development but not fully independent of public evidence.

Across the complete 775-variant benchmark, 661 exact alleles occur in the current local ClinVar snapshot and 114 do not. ClinVar presence was not the exclusion rule because a current ClinVar record does not show that its label entered model fitting. The primary safeguard is absence of the exact allele from original and accepted development partitions, combined with use of the source-specific outcome label.

### Final performance

The accepted 300-feature model reached:

- AUROC 0.984107; gene-clustered 95% CI 0.977170-0.989027.
- AUPRC 0.988132; 95% CI 0.976910-0.992414.
- Brier score 0.062939; 95% CI 0.048967-0.082540.
- Sensitivity 0.979358; 95% CI 0.967164-0.988272.
- Specificity 0.840708; 95% CI 0.770403-0.887760.
- Confusion matrix: TP 427, FP 54, TN 285, FN 9.

Source AUROCs: HiRO 0.9947, eMERGE 0.9854, CardioBoost 0.9632, PLOS 0.9772, and SHaRe 0.9863. The broader 1,233-variant exact-quarantine analysis remains a sensitivity, not the primary benchmark.

### New artifacts

- `datasets/external_validation_candidates/build_strict_external_nondevelopment_benchmark.py`
- `datasets/external_validation_candidates/audit_final_external_label_provenance.py`
- `datasets/external_validation_candidates/final_model_heldout_2026/external_unique_variant_manifest.tsv`
- `datasets/external_validation_candidates/final_model_heldout_2026/external_source_observations.tsv`
- `datasets/external_validation_candidates/final_model_heldout_2026/eligibility_audit.tsv`
- `datasets/external_validation_candidates/final_model_heldout_2026/label_provenance_audit.json`
- `results/model_performance/final_model_heldout_2026/`
- `paper_MS_SI/tables/Table_S23a_final_model_heldout_performance.tsv`
- `paper_MS_SI/tables/Table_S23b_final_model_heldout_gene_bootstrap.tsv`

### Concise manuscript wording

> We evaluated the locked model in 775 unique, label-concordant GRCh38 variants from five source-held-out cohorts. HiRO was institutionally internal but entirely external to model development; eMERGE, CardioBoost, PLOS, and SHaRe were institutionally external. To prevent exact-allele leakage, we excluded 66 PLOS and 419 SHaRe observations whose alleles had occupied the original training or early-stopping validation partitions, then deduplicated overlaps across sources. PLOS outcomes came from a published ACMG/AMP reinterpretation and SHaRe outcomes from registry-specific re-evaluation rather than from the model's ClinVar development label. No exact evaluation allele informed fitting or model selection, although ClinVar, literature, and other evidence resources can overlap across classification systems. We therefore describe this as a retrospective multi-source model-held-out evaluation, not fully evidence-independent or prospective clinical validation.

### Final verification

The final scripts compiled and the benchmark was rebuilt with seed 20260925 and 2,000 gene-clustered bootstrap samples. An independent manifest check confirmed 775 unique variants (436 P/LP, 339 B/LB) across 40 genes, 843 retained source observations, 64 alleles represented in multiple sources, and 68 duplicate source observations removed by exact-allele deduplication. None of the 775 variants occupied an original training or validation partition, and all 775 were assigned to `external_expanded` rather than training or validation in the accepted retrained matrix. The confusion matrix was reproduced as TP 427, FP 54, TN 285, and FN 9. The PLOS/SHaRe ClinVar-provenance audit was rerun against the current local snapshot and reproduced 661 benchmark alleles present and 114 absent.

### Repeated-allele analysis rule

The external data now have an explicit two-level interpretation. `external_source_observations.tsv` keeps one row for each source assertion, including repeated exact alleles, so source provenance and agreement remain auditable. `external_unique_variant_manifest.tsv` keeps one row per exact GRCh38 allele for combined performance. Concordant repeated classifications contribute once to pooled metrics; otherwise the identical molecular record and score would be given extra statistical weight. Discordant source labels remain in the conflict audit and are excluded from the primary pooled endpoint rather than resolved post hoc. Source-specific analyses may each include the allele once, but those estimates are not pooled as independent observations. The summary accounting is in `paper_MS_SI/tables/Table_S23c_final_model_heldout_source_accounting.tsv`.

## Multi-Source Consensus Integration (2026-09-26)

The final builder now emits all six requested analysis layers: source assertions, unique alleles, conflicts, source-specific estimates, the combined primary result, and a post hoc multi-source consensus sensitivity. New generated artifacts are `external_label_conflicts.tsv`, `source_accounting.tsv`, `source_overlap_combinations.tsv`, and `multi_source_consensus_predictions.tsv` under the final data/results directories. Tables S23a-S23c were regenerated directly from these outputs.

The one conflict is `MYBPC3 11-47349804-C-G`, classified P/LP by CardioBoost and B/LB by SHaRe. Both assertions remain auditable, but the allele is excluded from pooled accuracy without majority voting. The 64 consensus alleles have concordant labels in at least two retained sources, include 36 P/LP and 28 B/LB variants across 17 genes, and contribute one prediction each.

Consensus performance: AUROC 0.993056 (gene-clustered 95% CI, 0.953106-1.000000), AUPRC 0.994877 (0.957889-1.000000), Brier score 0.046536 (0.012819-0.087072), sensitivity 1.000000 (1.000000-1.000000), specificity 0.892857 (0.647059-1.000000), PPV 0.923077, and NPV 1.000000. Confusion matrix: TP 36, FP 3, TN 25, FN 0. This is explicitly reported as a post hoc supportive sensitivity, not independent replication.

Clinical interpretation added throughout the manuscript: the model can help a genetics service order evidence review, segregation follow-up, transcript reconciliation, functional-study consideration, and multidisciplinary discussion. Agreement across sources increases confidence that ranking is not tied to one classification source, while the conflict table identifies variants requiring reconciliation. Neither pooled nor consensus performance changes diagnosis or management without patient phenotype, inheritance, segregation, mechanism, and complete ACMG/AMP review.

The manuscript package audit and README now include the final S23 outputs. `ARTIFACT_CHECKSUMS.sha256` was expanded and refreshed against the current 80-file package after verification detected 25 stale hashes from earlier regenerated artifacts; all 80 entries now pass checksum validation. The structured abstract was tightened to 334 words while retaining the five-source primary result, post hoc consensus result, and clinical interpretation.

## Maximized FoldX / Delta-Delta-G Development Analysis (started 2026-09-26; completed 2026-09-27)

### Why FoldX was absent from development rows

The earlier FoldX source workflow began from 3,591 HiRO/eMERGE/CardioBoost records rather than the accepted training/validation matrix. Its 313 numeric values were therefore source-cohort values (19 HiRO, 4 eMERGE, 144 CardioBoost, and 146 ineligible/nonprimary rows), with zero development rows. A later registry attempt used the old frozen split and an unreliable FoldXPro workflow; it produced only 123 numeric training and 20 validation results. This was a workflow-scope/completion problem, not evidence that FoldX was impossible for development variants.

### Outcome-blind maximal mapping audit

Added `datasets/feature_sources/protein_structure/scripts/build_maximized_foldx_mapping.py`. It maps the accepted 85,968-row annotated matrix in this precedence order: direct AlphaMissense protein mapping; uniquely aligned dbNSFP UniProt/position fields; then exact VEP amino-acid substitution on the gene's dominant reviewed canonical UniProt accession. Every accepted mapping must be a single amino-acid substitution, agree with the reviewed UniProt reference sequence, fall within an AlphaFold v6 fragment, and agree with that fragment's reference residue. Labels are carried only for reporting and are never used in mapping or selection.

Primary-binary coverage:

| Split | Cohort rows | Any mapped pLDDT | pLDDT >=50 | pLDDT >=70 | pLDDT >=90 |
|---|---:|---:|---:|---:|---:|
| Train | 35,796 | 3,031 | 2,639 | 2,404 | 1,574 |
| Validation | 6,317 | 509 | 450 | 416 | 288 |
| Expanded external | 1,233 | 717 | 620 | 526 | 233 |

Mapping methods among training rows: 2,734 direct AlphaMissense, 141 aligned dbNSFP rescue, 16 exact flattened VEP/canonical rescue, and 140 transcript-level VEP rescues with one unique amino-acid consequence matching the reviewed canonical sequence. Of 35,796 training rows, 31,180 are not exact missense substitutions and therefore do not have a biologically defined single-substitution FoldX DDG; 1,585 are VEP-missense rows without one unambiguous reviewed-protein mapping after transcript-level rescue; 3,031 pass mapping. Among mapped rows, 392 have pLDDT below 50 and 235 have pLDDT 50-69. The prespecified analysis hierarchy is pLDDT >=70 primary, pLDDT >=50 sensitivity, and all mapped residues exploratory.

The final 775-allele model-held-out benchmark contains 529 computable mappings at any pLDDT, 441 at pLDDT >=50, and 382 at pLDDT >=70.

Generated mapping artifacts:

- `datasets/feature_sources/protein_structure/interim/maximized_foldx_mapping_audit.tsv`
- `datasets/feature_sources/protein_structure/interim/maximized_foldx_structure_features.tsv`
- `datasets/feature_sources/protein_structure/interim/maximized_foldx_mapping.summary.json`

### Full FoldX 5.1 calculation and model-comparison design

The corrected method runs FoldX 5.1 `RepairPDB` once per AlphaFold v6 fragment and batched `BuildModel` calculations thereafter. A resumed pLDDT >=70 repair run is active; at this checkpoint, 80/82 old-scope fragments were repaired and the MYH7 and DSP large-fragment repairs were still computing. The expanded all-mapped preflight contains 4,257 variant rows, 4,137 unique structure/mutation calculations, 87 fragments, and 139 batches. Do not report DDG performance until that run and the expanded run complete.

The FoldX runner cache was corrected to verify exact calculation keys and mutations before accepting a cached batch; equal row count alone is no longer sufficient. Its selection statement now records the actual pLDDT threshold rather than hard-coding 70.

Added `build_maximized_foldx_augmentations.py` to generate matched eligibility-only and signed/absolute DDG augmentation tables at pLDDT >=70, >=50, and any mapped pLDDT after calculations finish. Added `evaluate_maximized_foldx_models.py` for whole-validation, covered-validation, all-775-held-out, and covered-held-out metrics with paired gene-cluster bootstrap differences.

The required comparison is: (1) baseline retrained after excluding legacy `foldx_ddg_*` fields; (2) FoldX-eligibility-only control; and (3) eligibility plus numeric signed/absolute DDG. Because high-confidence FoldX coverage is enriched for pathogenic missense variants (training pLDDT >=70: 1,987 P/LP and 322 B/LB), improvement over the baseline alone is insufficient. A defensible incremental-DDG claim requires the numeric DDG model to improve over the eligibility-only control, especially among the 380 high-confidence covered held-out alleles. A null result must be reported as such; inclusion alone is not novelty.

### Remaining steps

1. Finish the active corrected FoldX repair/build run; inspect every failed/timeout/parse status.
2. Run the expanded all-mapped calculation using the maximized audit and reuse only mutation-verified caches.
3. Build pLDDT-tier augmentation files.
4. Retrain matched baseline, eligibility-only, and DDG models with identical splits, weights, seed, and CatBoost settings.
5. Run whole-cohort and covered-cohort paired evaluation, then decide whether the manuscript supports incremental value, feasibility-only reporting, or a null ablation result.
6. Only after final results, update the MS/SI, tables, progress record, and package checksums.

### Terminal handoff (2026-09-27)

The corrected old-scope pLDDT >=70 run completed successfully: 2,951 variant rows, 2,860 unique calculations, 82/82 repaired fragments, and 2,860/2,860 numeric DDG results with no failed, timeout, or parse statuses. These results and 117 completed batch caches are retained under `interim/registry_foldx51_repaired/`.

The expanded driver was started and then intentionally interrupted at the user's request to conserve Codex credits. No FoldX process remains active. The five expanded-only fragments (`P0DP23`, `P0DP24`, `Q13936`, and RYR2 fragments F4/F15) remain pending. Resume with `datasets/feature_sources/protein_structure/scripts/run_maximized_foldx_analysis.sh`; it will reuse successful calculations by exact structure/mutation key and then run augmentations, seven matched models, and all three tier evaluations.

### Completed maximized FoldX results (2026-09-27)

The user resumed the terminal workflow. All 87 required AlphaFold fragments repaired successfully. FoldX returned numeric values for all 4,137 unique structure-mutation calculations and 4,257 primary-binary variant rows: 3,031 training, 509 validation, and 717 expanded external. There were no repair, build, timeout, or parse failures. The final 775-allele benchmark contains 382 covered alleles at pLDDT >=70, 441 at >=50, and 529 at any pLDDT.

Seven matched models completed: one clean baseline excluding legacy `foldx_ddg_*` fields, plus eligibility-only and numeric-DDG models at pLDDT >=70, >=50, and any mapped pLDDT. Corrected tier-specific evaluations used 2,000 paired gene-cluster bootstrap samples.

Primary pLDDT >=70 result across all 775 held-out alleles:

- Baseline: AUROC 0.98393, AUPRC 0.98787, Brier 0.06229.
- Eligibility only: AUROC 0.98548, AUPRC 0.98871, Brier 0.05743.
- Eligibility + numeric DDG: AUROC 0.98298, AUPRC 0.98666, Brier 0.06427.
- DDG minus eligibility: AUROC -0.00233 (95% CI -0.00390 to -0.00109), AUPRC -0.00203 (-0.00491 to -0.00046), Brier +0.00704 (+0.00460 to +0.01041). Numeric DDG was worse than eligibility alone.

pLDDT >=50 result across all 775:

- Eligibility + DDG improved modestly over baseline: AUROC +0.00089 (+0.00019 to +0.00171), AUPRC +0.00079 (+0.00014 to +0.00201), Brier -0.00315 (-0.00503 to -0.00101).
- It did not clearly improve over eligibility alone: all corresponding intervals included zero.
- Among 441 covered held-out alleles, DDG versus baseline changed AUROC by +0.00235 (+0.00044 to +0.00513) and Brier by -0.00552 (-0.00933 to -0.00190), but again did not beat eligibility alone.

Using every mapped residue did not improve over baseline and worsened Brier versus eligibility alone by +0.00249 (+0.00045 to +0.00417). The scientific conclusion is therefore a completed negative incremental-value ablation: registry-scale, reference-validated FoldX integration is feasible, but generic DDG does not add independent predictive value beyond mapping eligibility and the existing CardioQueue predictors. Do not claim FoldX improved CardioQueue or make a first-use claim. The defensible novelty is multi-gene cardiogenetic scale, complete reference checking, long-protein fragment handling, confidence tiers, exact external quarantine, and an eligibility-only control that prevents missingness/coverage from masquerading as stability signal.

New/updated artifacts:

- `datasets/feature_sources/protein_structure/interim/maximized_foldx51_ddg.tsv`
- `datasets/feature_sources/protein_structure/interim/maximized_foldx51_ddg.summary.json`
- `results/models/maximized_foldx_2026/`
- `results/model_performance/maximized_foldx_enhancement/{ge70,ge50,any}/`
- `paper_MS_SI/tables/Table_S17_foldx_coverage_eligibility.tsv`
- `paper_MS_SI/tables/Table_S17b_foldx_tier_metrics.tsv`
- `paper_MS_SI/tables/Table_S17c_foldx_paired_bootstrap.tsv`
- `paper_MS_SI/tables/Table_S17d_foldx_feature_importance.tsv`

The main manuscript, SI S17, README, completion audit, and novelty memo were updated to replace the unfinished-run language with the completed results and clinical interpretation.

### Final FoldX package verification (2026-09-27)

Removed the obsolete high-confidence preflight table after replacing it with the complete confidence-tier tables. Corrected the FoldX runner's cached-progress display so cached rows cannot make the submitted-job counter appear off by one; this display-only correction does not alter any calculation or result.

The reproducibility manifest now covers 168 artifacts, including maximized mapping and calculation outputs, augmentation tables, all seven matched models, all three tier evaluations, scripts, and Tables S17-S17d. All 168 checksums verified. Python compilation, shell syntax, manuscript Markdown table widths, and direct assertions against the FoldX summary, coverage table, and held-out metrics passed. The structured abstract is 315 words; the final manuscript and SI contain 10,838 and 15,154 whitespace-delimited words, respectively. No FoldX, model-training, or evaluation process remains active.

Final interpretation: FoldX was successfully scaled to 4,137 unique reference-validated calculations, including 3,031 training variants. It did not show independent predictive benefit over matched mapping eligibility and therefore was not added to the production CardioQueue model. This completed negative ablation is reported transparently and supports a narrower novelty claim about auditable, multi-gene cardiogenetic structure integration rather than improved classification from generic delta-delta-G.

## Complete Historical-VUS Temporal Queue Upgrade (2026-09-27)

The prior S19 analysis conditioned the queue on 457 January 2024 VUS already known to have a later binary classification. Although its separate temporal model excluded those VUS from fitting, the 117-review workload and 100% precision statements were not representative of a real backlog because 25,280 still-VUS or currently nonbinary records were absent from the denominator. The ordinary latest-label CardioQueue model was confirmed unsuitable for this test: 448/457 resolved VUS occupied its ordinary training or validation partitions.

Five existing agent threads were assigned independent cohort-integrity, leakage, statistical, clinical, and figure/manuscript audits. Four completed and were integrated. The fifth cohort-integrity thread hit its usage limit; its assertions were independently covered by the statistical/figure agents and direct local tests.

The temporal builder now provides an in-scope split. Only January 2024 B/LB or P/LP variants within the prespecified model scope enter fitting: 23,851 training and 4,210 validation rows. Every January 2024 VUS remains excluded. A separate inference pass scores all 25,737 historical VUS, including the primary 25,481 scope-eligible records and 256 outside-scope sensitivity records. The primary fixed-horizon accounting is 146 observed later P/LP, 309 later B/LB, 24,879 still VUS, and 147 without an unambiguous current class. Still-VUS records are never coded as benign.

Four temporal models were trained with identical in-scope splits: full public 300 features; 242 features after removing eight common effect-predictor families; 191 low-circularity features after additionally removing SpliceAI, availability indicators, and current scope/ClinGen fields; and 19 static molecular features. Resolved-subset AUROCs were 0.9970, 0.9963, 0.9936, and 0.9892, respectively. These are conditional metrics, not full-queue accuracy or prospective validation.

In the complete 25,481-record primary queue, the full model recovered 16/146 observed later-P/LP transitions in the first 500 reviews (11.0%; gene-clustered 95% CI 5.6%-23.7%) and 93/146 in the first 5,000 (63.7%; 45.9%-90.5%). Eighty-percent recovery required 6,644 reviews versus a random-order median of 20,313 (95% simulation interval 18,550-21,831). The low-circularity model recovered 24 and 94 at the same budgets. Yield is explicitly labelled a lower bound because unresolved records have unknown outcomes.

Added a two-ended evidence-refresh analysis: high scores nominate possible P/LP evidence review and low scores nominate possible B/LB evidence review. Added 5,000 gene-cluster bootstrap replicates, 10,000 random orderings, absolute review budgets, score-decile accounting, the complete prediction file, a redesigned Figure 5 temporal panel, and new Figure S6.

Primary new artifacts:

- `results/model_performance/temporal_clinvar_2024/evaluate_all_historical_vus.py`
- `results/model_performance/temporal_clinvar_2024/all_historical_vus_predictions.tsv`
- `results/models/cardioqueue_temporal_clinvar_2024_in_scope_{public300,no_effect_predictors,low_circularity,static19}/`
- `paper_MS_SI/tables/Table_S19_all_vus_transition_counts.tsv`
- `paper_MS_SI/tables/Table_S19b_all_vus_queue_budgets.tsv`
- `paper_MS_SI/tables/Table_S19c_all_vus_queue_targets.tsv`
- `paper_MS_SI/tables/Table_S19d_resolved_vus_model_metrics.tsv`
- `paper_MS_SI/tables/Table_S19e_all_vus_gene_bootstrap.tsv`
- `paper_MS_SI/tables/Table_S19f_all_vus_score_deciles.tsv`
- `paper_MS_SI/tables/Table_S19g_all_vus_two_ended_queue.tsv`
- `paper_MS_SI/figures/Figure_5_review_queue_cumulative_gain.png`
- `paper_MS_SI/figures/Figure_S6_all_historical_vus_bidirectional_gain.png`

The main manuscript, SI S19, contents, legends, README, and completion audit were updated. The old resolved-only S19 queue tables were removed. The current comparison file is identified by its 1 July 2026 UTC timestamp and SHA-256. The reproducibility manifest now contains 219 artifacts.
