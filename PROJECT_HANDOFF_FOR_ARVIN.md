# CardioQueue project handoff for Arvin and Codex

## Start here

This repository contains the code, documentation, frozen public release, manuscript package, aggregate results, and provenance needed to finish the CardioQueue paper. The full working directory on the originating computer was approximately 164 GB. Large raw/intermediate datasets, local environments, licensed software, AF3 bulk outputs, redundant trained models, and participant-linked HiRO data are not in Git.

The current manuscript is a strong working draft with verified analyses, but it is **not yet ready for submission**. Arvin's latest editorial comments in `paper_MS_SI/Arvin Feedback.txt` have not all been applied, and author-controlled declarations remain unresolved.

Recommended reading order for a fresh Codex session:

1. This file.
2. `paper_MS_SI/Arvin Feedback.txt`.
3. `paper_MS_SI/MS_BMC_MEDICAL_GENOMICS_FINAL.md`.
4. `paper_MS_SI/SI_BMC_MEDICAL_GENOMICS_FINAL.md`.
5. `paper_MS_SI/COMPLETION_AUDIT.md`.
6. `paper_MS_SI/README.md` and `paper_MS_SI/Manuscript Writing skills.txt`.
7. Only then open the specific tables, figures, reports, manifests, or scripts needed for the task.

## Project in one paragraph

CardioQueue is a CatBoost-based research model that ranks inherited-cardiac-disease variants for expert evidence review. It combines consequence, population-frequency, conservation, established effect-prediction, gene, protein-sequence, and structure-related annotations. The primary supervised outcome is benign/likely benign (B/LB) versus pathogenic/likely pathogenic (P/LP); VUS were excluded from primary fitting and later ranked as a review backlog. The output is a continuous review-priority score. It is not an ACMG/AMP classification, a patient-priority score, a diagnosis, or a treatment recommendation.

## What was done

### 1. Sources and registry

- ClinVar supplied the main development labels and the historical VUS analysis.
- HiRO/CASPER WES/VERDICT supplied an institutionally internal, model-held-out cohort. These participant-linked records are private.
- eMERGE-III, released CardioBoost records, SHaRe HCM classifications, and the Fernandez-Falgueras/PLOS inherited-cardiovascular reinterpretation cohort supplied additional evaluation assertions.
- Exact variants were normalized to GRCh38 and reference-checked. Cross-source duplicates and discordant labels were tracked rather than silently majority-voted.
- The original registry contained 85,677 unique rows: 42,361 VUS, 33,923 B/LB, 9,148 P/LP, and 245 other/missing labels.

### 2. Annotation and feature engineering

The pipeline brought together VEP, SpliceAI, gnomAD, dbNSFP, ClinGen, UniProt, AlphaMissense, AlphaFold, DSSP, FreeSASA, and selected protein-context features. Outcome fields, exact identifiers, source flags, participant/source-record identifiers, review metadata, private HiRO variables, and label aggregates were excluded from the public 300-feature inference schema.

Protein substitutions were mapped only when accession, position, and reference residue were consistent across annotation, UniProt, and AlphaFold. FoldX 5.1 delta-delta-G calculations were analyzed separately because FoldX is licensed and because structure availability can itself carry signal.

### 3. Model development

- The **CardioQueue primary model** is the current 300-public-feature binary CatBoost model.
- A union of 1,233 candidate evaluation alleles was quarantined before fitting.
- The final primary fit used 35,796 training rows and 6,317 early-stopping/validation rows.
- The main model-held-out benchmark contained 775 exact-allele-deduplicated variants (436 P/LP and 339 B/LB) across 40 genes and five source strata.
- The **legacy source-held-out model** used the earlier fixed 430-variant benchmark and remains only for analyses whose predictions and comparator rankings were frozen together.
- The **temporal CardioQueue model** was trained on January 2024 binary ClinVar labels and scored January 2024 VUS. It is a separate label-temporal analysis.
- A true three-class model exists only as an exploratory secondary analysis. It is not the main model.

### 4. Evaluation and sensitivity analyses

The project evaluated pooled and source-specific AUROC/AUPRC, fixed-threshold performance, exploratory 0.1/0.9 review zones, gene-clustered bootstrap intervals, calibration, leakage ablations, higher-review-status labels, removal of precomputed effect predictors, gene-transfer stress tests, structure geometry, FoldX, a direct executable CardioBoost comparison, review-budget simulations, and retrospective temporal enrichment of later ClinVar classifications.

The direct CardioBoost comparison uses harmonizable shared public observations. It does not reproduce CardioBoost's private development or patient-outcome analyses.

## Current headline results

Use the current main manuscript and `paper_MS_SI/tables/` as the authoritative paper-facing sources. Do not substitute older values from `PROJECT_PROPOSAL-PROGRESS.md` or older result bundles when they differ.

- Primary 775-allele benchmark: AUROC 0.9841 (gene-clustered 95% CI 0.9772-0.9890) and AUPRC 0.9881 (0.9769-0.9924).
- At threshold 0.5: sensitivity 0.9794 (0.9672-0.9883) and specificity 0.8407 (0.7704-0.8878).
- Exploratory 0.1/0.9 zones covered 81.4% of the benchmark without deferral and had 96.7% accuracy in that nondeferred subset. These are workload zones, not clinical probability cutoffs.
- Strict shared CardioBoost comparison: 235 observations representing 211 unique variants; AUROC 0.984 for CardioQueue versus 0.864 for CardioBoost, paired difference 0.120 (95% CI 0.068-0.176).
- Public research worklist: 41,622 VUS after excluding all HiRO-linked rows.
- Temporal analysis: 25,481 January 2024 VUS scored; 146 had a later observed P/LP classification by the comparison snapshot. This is retrospective enrichment, not prospective prediction.
- Structure geometry produced only a small sensitivity-model change. FoldX did not add independent value at the prespecified pLDDT >=70 tier and should not be credited with the primary performance.

## Scientific framing that must be preserved

- CardioQueue prioritizes **variants for evidence review**, not patients.
- It does not assign ACMG/AMP classifications or replace expert adjudication.
- VUS is an evidence category, not a stable biological class. The primary model is binary; applying review zones to VUS is prioritization, not reclassification.
- Scores and the 0.1/0.9 zones are not calibrated clinical probabilities.
- Model-held-out does not guarantee evidence independence because sources may share ClinVar, literature, laboratory, or consortium evidence.
- HiRO is held out from model development but is institutionally internal. eMERGE, CardioBoost, SHaRe, and PLOS are institutionally external.
- The SHaRe/PLOS and 64-allele consensus analyses have important reuse, overlap, size, and shared-evidence limitations.
- The temporal analysis is label-temporal, not fully time-locked, because some predictor resources postdate January 2024.
- The review-queue results describe enrichment in constructed retrospective cohorts, not proven time savings, clinical safety, or patient benefit.

## Manuscript status and Arvin's requested direction

The files named `FINAL` are the latest assembled versions, not author-approved submission files. Arvin requested a less technical and more clinically readable paper. The next editing pass should:

- remove code-like workflow language from the main text and retain only scientifically necessary methods;
- strengthen the clinical cardiomyopathy/genetics background;
- use one consistent name: **CardioQueue primary model**, with the legacy and temporal models named only when genuinely needed;
- remove unnecessary implementation details, class-weight numbers, seeds, tree counts, and computer-language explanations from the prose while keeping reproducibility in code and repository records;
- simplify tables, removing metrics such as Brier score from clinical-facing tables unless essential to the argument;
- reduce the supplementary burden substantially; 23 supplementary tables is not an acceptable final presentation;
- consolidate supplementary tables into a single XLSX and figures into a single captioned PDF if that remains the author decision;
- place tables and figures at the end of the manuscript and give each figure one clear results-section story;
- replace crowded or placeholder-like flow charts and harmonize plot styles/legends;
- keep references as verified inline links if Zotero-compatible linked text is desired;
- do not add an abstract, keywords, author block, or other front matter until Arvin asks for them;
- use direct, adult, non-vague scientific language.

Before submission, authors must provide or approve authors/affiliations, corresponding-author details, ethics and consent/waiver language, funding, competing interests, contributions, acknowledgements, controlled-access procedure, licenses, final repository URL/commit, and archive DOI. Never invent these.

## Repository map

| Path | Purpose | Share status |
|---|---|---|
| `paper_MS_SI/` | Current manuscript, SI, feedback, figures, tables, evidence notes, checksums | Included; primary writing workspace |
| `release/cardioqueue_v1.0/` | Compact public model, 300-feature schema, inference script, provenance, and public VUS worklist | Included; intended public release |
| `datasets/**/scripts/` and provenance scripts | Data acquisition, normalization, feature engineering, and registry code | Included where not data-bearing |
| `results/` | Analysis scripts, aggregate reports, figures, metrics, and many local row-level/model products | Code and aggregate artifacts included; heavy/row-level products excluded |
| `docs/` | Design decisions, data organization, audits, and implementation plans | Included |
| `metadata/` | Source manifests and feature requirements | Included |
| `AF3-systematic protein/` | Systematic AF3 design reviews and job definitions | Reviews/definitions included; server outputs excluded |
| `PROJECT_PROPOSAL-PROGRESS.md` | Long chronological development log | Included but not authoritative over newer manuscript artifacts |
| `improving BMC paper progress.MD` | Later paper-improvement and expanded-evaluation log | Included as historical decision context |

## What is deliberately absent from Git

- the roughly 150 GB local `datasets/` payload: downloaded public sources, raw files, intermediates, ready matrices, archives, VEP/dbNSFP/SpliceAI/AlphaFold assets, and most row-level tables;
- participant-linked or row-level HiRO/CASPER WES/VERDICT data;
- SHaRe-derived row-level data unless written redistribution permission is obtained;
- the local `.tools/` software environments and package caches;
- licensed FoldX executables and local FoldX working files;
- AF3/AlphaFold server outputs, MSAs, bulky structures, and full-data JSON files;
- redundant experimental CatBoost model binaries and row-level prediction tables under `results/`;
- local backups, temporary files, logs, and editor metadata.

The clone is sufficient to edit/finalize the paper, inspect code and aggregate evidence, run the compact public inference package, and reproduce lightweight figure/table generation when its inputs are included. It is **not** sufficient for a clean-room rebuild of every analysis.

## Environment and running the public release

The analysis environment snapshot is `paper_MS_SI/environment-manuscript.txt`. The core recorded versions include Python-compatible CatBoost 1.2.10, pandas 2.3.3, NumPy 2.0.2, scikit-learn 1.6.1, statsmodels 0.14.6, SciPy 1.13.1, and Matplotlib 3.9.4. Structure analyses additionally used mkdssp 4.6.1 and FreeSASA 2.2.1.

For the compact public scorer:

```bash
cd release/cardioqueue_v1.0
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python score_variants.py annotated_variants.tsv cardioqueue_scores.tsv
```

The input must contain `variant_id` plus every ordered field in `model/feature_columns.json`. See `release/cardioqueue_v1.0/README.md` and `FEATURE_PROVENANCE.md`. Do not interpret the output as a clinical classification or calibrated pathogenicity probability.

## Reacquiring data for a full rebuild

Public inputs must be downloaded again from their original providers under their current terms. Start with `metadata/datasets.json`, `datasets/README.md`, `docs/data_organization.md`, and the source-specific scripts/READMEs. Verify checksums and versions; several resources are version-sensitive or restrict redistribution. In particular, dbNSFP and AlphaMissense have their own licenses, FoldX requires a separate license, and VEP/ClinVar/gnomAD snapshots can change.

Private HiRO data require the institution's approved access route and must stay outside Git. A full external benchmark rebuild may also require written SHaRe permission or local execution against an authorized copy. If data cannot be obtained, do not claim that every analysis was independently reproduced; use the frozen aggregate tables, manifests, checksums, and audit trail for manuscript verification.

## Recommended next Codex task

Ask Codex to revise the manuscript from `paper_MS_SI/Arvin Feedback.txt` while preserving every supported result and limitation. A useful prompt is:

> Read `AGENTS.md`, `PROJECT_HANDOFF_FOR_ARVIN.md`, `paper_MS_SI/Arvin Feedback.txt`, the current MS/SI, and the completion audit. Propose and then implement a clinically readable manuscript revision that applies Arvin's feedback. Verify all retained numbers against the current paper tables/manifests, keep author-controlled fields explicit, do not invent citations or analyses, and do not expose private HiRO or unlicensed data.

After editorial agreement, generate the final journal-format artifacts, render and inspect every figure/table/document, run consistency checks, update repository/commit/license fields, and create an immutable archived release/DOI.

