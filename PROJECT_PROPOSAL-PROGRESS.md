# Project Proposal

## Paper-Facing Current Status, 2026-07-03

This section supersedes older exploratory notes below. Earlier parts of this
progress document are intentionally preserved because they show how decisions
were made, but the paper-facing counts and claims should use the values in this
section unless a later dated section explicitly replaces them.

### Main Claim

The strongest current manuscript claim is:

> A larger, feature-richer, externally validated cardiogenetics CatBoost model
> improves high-confidence pathogenicity prioritization on CardioBoost-eligible
> matched external rows while preserving an explicit VUS deferral zone.

The primary model is not a true three-class model. It is a binary
Pathogenic/Likely Pathogenic versus Benign/Likely Benign model with a
CardioBoost-style deferral zone:

| Probability | Output |
|---:|---|
| `<=0.1` | Benign-like |
| `0.1-0.9` | Deferred / uncertain |
| `>=0.9` | Pathogenic-like |

This is the reviewer-safe framing because VUS is an evidence status, not a
stable biological class. The true three-class model remains a secondary
analysis for VUS behavior and label-matching.

### Count Reconciliation

Two row counts appear in the project because they refer to different stages.

| Matrix | Rows | Columns | Use |
|---|---:|---:|---|
| Annotated intermediate feature matrix | `86,889` | `463` | Large GRCh38 registry after feature annotation, before final ready-table cleanup |
| Analysis-ready modeling matrix | `85,677` | `420` | Current matrix used for final splits, weights, external validation, calibration, stratified analyses, and all-label triage |

Use `85,677` for manuscript-facing current results. Use `86,889` only when
describing the older annotated intermediate matrix.

Current analysis-ready label counts:

| Label | Rows |
|---|---:|
| Benign / likely benign | `33,923` |
| VUS | `42,361` |
| Pathogenic / likely pathogenic | `9,148` |
| Missing, conflict, or other non-supervised label | `245` |

Training slice counts:

| Slice | Rows | Composition |
|---|---:|---|
| Primary binary training slice | `42,990` | `33,859` benign, `9,131` pathogenic |
| Exploratory three-class slice | `84,721` | `33,859` benign, `41,731` VUS, `9,131` pathogenic |
| VUS scoring pool | `42,361` | VUS rows scored after binary training |

### Completed Analyses

| Analysis | Status | Output |
|---|---|---|
| Primary binary CatBoost source-held-out model | Done | `results/models/primary_binary_catboost_v0_strict_source_heldout_hiro_emerge_rescued/` |
| All-label binary three-zone triage | Done | `results/model_performance/primary_binary_three_zone_all_sources/` |
| True three-class CatBoost model | Done, secondary | `results/model_performance/true_three_class_all_sources/` |
| Three-class weight/threshold tuning | Done, no external-test tuning used | `results/model_performance/three_class_threshold_tuning/` |
| CardioBoost same-row public-model benchmark | Done | `results/model_performance/cardioboost_matched_benchmark/` |
| CardioBoost paired bootstrap confidence intervals | Done | `results/model_performance/cardioboost_matched_benchmark/` |
| Leakage and ablation audit | Done | `results/model_performance/leakage_ablation_audit/` |
| Calibration analysis and reliability plot | Done | `results/manuscript_figures_tables/` |
| Stratified performance tables | Done | `results/manuscript_figures_tables/tables/stratified_performance.tsv` |
| Feature coverage and source-flow figures | Done | `results/manuscript_figures_tables/figures/` |
| HiRO unresolved-record audit | Done | `results/data_characteristics/hiro_unresolved_record_audit/` |

### Primary External Validation Result

Source-held-out binary validation at the standard `0.5` threshold:

| External set | Rows | AUROC | AUPRC | Sensitivity | Specificity | PPV | Brier |
|---|---:|---:|---:|---:|---:|---:|---:|
| Combined external | `430` | `0.977` | `0.985` | `0.965` | `0.886` | `0.925` | `0.053` |
| HiRO | `69` | `0.988` | `0.983` | `0.963` | `0.881` | `0.839` | `0.068` |
| eMERGE | `176` | `0.991` | `0.991` | `0.979` | `0.938` | `0.949` | `0.030` |
| CardioBoost | `185` | `0.952` | `0.981` | `0.955` | `0.811` | `0.926` | `0.069` |

CardioBoost-style `0.1/0.9` high-confidence threshold result:

| External set | Rows | Sensitivity | Specificity | PPV | NPV | Deferral |
|---|---:|---:|---:|---:|---:|---:|
| Combined external | `430` | `0.937` | `0.943` | `0.960` | `0.977` | `0.119` |
| HiRO | `69` | `0.926` | `0.929` | `0.893` | `1.000` | `0.087` |
| eMERGE | `176` | `0.948` | `0.988` | `0.989` | `1.000` | `0.114` |
| CardioBoost | `185` | `0.932` | `0.887` | `0.953` | `0.903` | `0.135` |

### CardioBoost Matched Benchmark

The official public CardioBoost cardiomyopathy and arrhythmia AdaBoost model
objects were run locally on the official public all-rare mutation tables.
Official CardioBoost scores were matched to our source-held-out CatBoost
predictions by gene and HGVS cDNA. Exact coordinate matching was used as QC but
not as the primary benchmark because it recovered too few binary rows, likely
due to build and transcript-version differences.

This is a fair public-model same-row benchmark. It is not a reconstruction of
CardioBoost's original private training and clinical outcome cohorts.

External all, `250` matched rows:

| Model | AUROC | AUPRC | Sensitivity | Specificity | PPV | NPV | Deferral | High-conf accuracy |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Our CatBoost | `0.980` | `0.990` | `0.959` | `0.543` | `0.964` | `1.000` | `0.152` | `0.972` |
| Official CardioBoost | `0.843` | `0.902` | `0.876` | `0.148` | `0.836` | `1.000` | `0.244` | `0.847` |

Paired bootstrap differences, our model minus official CardioBoost:

| Metric | Difference | 95% CI |
|---|---:|---:|
| AUROC | `0.138` | `0.083-0.195` |
| AUPRC | `0.087` | `0.048-0.132` |
| Sensitivity | `0.083` | `0.030-0.136` |
| PPV | `0.128` | `0.082-0.172` |
| Deferral | `-0.092` | `-0.156 to -0.028` |
| High-conf accuracy | `0.125` | `0.082-0.167` |

External CardioBoost-source-only, `164` rows:

| Model | AUROC | AUPRC | Sensitivity | Specificity | PPV | NPV | Deferral | High-conf accuracy |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Our CatBoost | `0.972` | `0.990` | `0.967` | `0.366` | `0.952` | `1.000` | `0.146` | `0.957` |
| Official CardioBoost | `0.788` | `0.905` | `0.870` | `0.073` | `0.843` | `1.000` | `0.207` | `0.846` |

### Calibration and Stratified Performance

Calibration analysis, Brier scores, and reliability plots are complete.

| Subset | Rows | Brier | Mean predicted P/LP probability | Observed P/LP fraction | Mean error |
|---|---:|---:|---:|---:|---:|
| Validation | `6,384` | `0.003` | `0.211` | `0.208` | `0.003` |
| External all | `430` | `0.053` | `0.627` | `0.593` | `0.034` |
| HiRO | `69` | `0.068` | `0.444` | `0.391` | `0.052` |
| eMERGE | `176` | `0.030` | `0.575` | `0.545` | `0.029` |
| CardioBoost | `185` | `0.069` | `0.745` | `0.714` | `0.031` |

Interpretation: ranking and triage performance are strong, but the model is
mildly overconfident on external data. If the final manuscript claims calibrated
absolute probabilities, run validation-only Platt or isotonic calibration and
evaluate it externally.

Current feature coverage:

| Feature group | Covered rows | Coverage |
|---|---:|---:|
| VEP | `85,660 / 85,677` | `99.98%` |
| SpliceAI | `77,141 / 85,677` | `90.04%` |
| gnomAD final AF | `62,937 / 85,677` | `73.46%` |
| gnomAD observed | `55,387 / 85,677` | `64.65%` |
| dbNSFP | `43,774 / 85,677` | `51.09%` |
| AlphaMissense direct | `38,023 / 85,677` | `44.38%` |
| FoldX DDG | `313 / 85,677` | `0.37%` |

Key stratified external-all results:

| Stratum | Rows | AUROC | AUPRC | Sensitivity | Specificity | PPV |
|---|---:|---:|---:|---:|---:|---:|
| Arrhythmia genes | `265` | `0.985` | `0.988` | `0.975` | `0.907` | `0.939` |
| Cardiomyopathy genes | `151` | `0.965` | `0.980` | `0.946` | `0.847` | `0.906` |
| Missense consequence | `315` | `0.987` | `0.989` | `0.977` | `0.887` | `0.914` |
| dbNSFP matched | `348` | `0.985` | `0.989` | `0.981` | `0.896` | `0.938` |
| dbNSFP missing | `82` | `0.949` | `0.962` | `0.878` | `0.854` | `0.857` |
| gnomAD observed | `326` | `0.984` | `0.984` | `0.977` | `0.901` | `0.919` |
| gnomAD confirmed absent | `32` | `0.787` | `0.839` | `0.692` | `0.737` | `0.643` |
| Sparse P/LP genes | `24` | `0.986` | `0.981` | `1.000` | `0.929` | `0.909` |

The gnomAD-confirmed-absent stratum is small and weaker. Report it as a
high-caution descriptive subgroup, not as a failed main endpoint.

### All-Label Clinical Triage

The binary model was applied to all Benign, VUS, and Pathogenic rows using the
same three probability zones. This is not a true three-class classifier.

| Dataset | Rows with prediction | Benign | VUS | Pathogenic | Exact 3-zone agreement | P/LP sensitivity | Specificity vs non-P/LP | P/LP PPV | VUS deferral |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| All variant rows | `85,677` | `33,923` | `42,361` | `9,148` | `0.670` | `0.981` | `0.763` | `0.332` | `0.178` |
| ClinVar variant rows | `83,915` | `33,796` | `40,872` | `9,038` | `0.674` | `0.982` | `0.763` | `0.334` | `0.174` |
| HiRO variant rows | `240` | `43` | `139` | `31` | `0.568` | `0.903` | `0.813` | `0.452` | `0.325` |
| eMERGE variant rows | `2,754` | `130` | `2,394` | `111` | `0.498` | `0.946` | `0.715` | `0.127` | `0.427` |
| CardioBoost variant rows | `353` | `57` | `0` | `137` | `0.820` | `0.934` | `0.895` | `0.955` | `0.261` |
| HiRO source records | `408 / 482` | `190` | `172` | `46` | `0.647` | `0.891` | `0.845` | `0.423` | `0.250` |
| eMERGE source records | `2,754` | `137` | `2,490` | `127` | `0.497` | `0.945` | `0.719` | `0.140` | `0.427` |
| CardioBoost source records | `355` | `156` | `0` | `199` | `0.642` | `0.869` | `0.808` | `0.852` | `0.259` |

Important hard-error counts:

| Dataset | P/LP-to-benign errors | Benign-to-P/LP escalations | VUS-to-P/LP prioritizations |
|---|---:|---:|---:|
| All variant rows | `22` | `43` | `18,020` |
| ClinVar variant rows | `19` | `36` | `17,666` |
| HiRO source records | `0` | `18` | `38` |
| eMERGE source records | `0` | `1` | `738` |
| CardioBoost source records | `5` | `30` | `0` |

The VUS-to-P/LP number is large because the model was trained to identify
pathogenic-like evidence, not to learn VUS as a biological class. This is useful
for prioritization, but it requires future reclassification, expert review, or
case-level validation before clinical claims.

### Remaining Work Before Submission

| Priority | Task | Why |
|---:|---|---|
| 1 | Add bootstrap confidence intervals for source-held-out external validation and key stratified results | Needed for final tables, especially small strata |
| 2 | Test optional validation-only probability calibration | Useful only if final claims depend on calibrated probabilities |
| 3 | Assemble final multipanel manuscript figures | Combine ROC/PR, calibration, source flow, feature coverage, triage, and stratified results |
| 4 | Continue HiRO unresolved-record rescue | Current model scores `408 / 482`; remaining `74` need HGVS normalization or manual review |
| 5 | Finalize strict wording for CardioBoost comparison | Prevents overclaiming versus the private CardioBoost study |

## ClinVar as the Primary Label Source for a CatBoost Variant Classification Project

Yes — but not in the clean way you probably want.

ClinVar usually provides the final clinical classification for a variant, but it does not always provide a clean, structured ACMG table for every entry.

## What ClinVar reliably provides

ClinVar provides the following information that is useful for the model:

- Clinical significance: such as Pathogenic / Likely Pathogenic, VUS, or Benign / Likely Benign
- Review status: including whether criteria were provided, whether multiple submitters were involved, or whether an expert panel reviewed the classification
- Condition or phenotype information, although this often needs cleaning and normalization
- Submitter-level assertions
- Conflict status
- Indirect evidence of whether assertion criteria were provided

ClinVar explicitly states that it does not compute conclusions itself; instead, it reports classifications submitted by external submitters, while conflicts and uncertainty are shown clearly.

## What “criteria provided” means

When ClinVar says “criteria provided,” it means that the submitter supplied assertion criteria or a method for how the variant was classified. For one-star review status, ClinVar defines this as a single submitter with criteria provided.

However, this does not always mean that the variant-level ACMG codes are available. For example, a submission may only say something like “classified according to ACMG/AMP guidelines” without providing a structured list of specific ACMG rules such as PM2 + PP3 + PS4.

## Where ACMG codes may appear

In some cases, submitters include ACMG codes in the comment or evidence text. Some ClinVar pages explicitly state that the following ACMG criteria were applied and list codes such as PM1, PM2, PM5, PP3, and PP5.

But this is not standardized enough to rely on as a complete structured ACMG dataset.

## Best approach for this study

### Use ClinVar for labels

ClinVar should be used as a source of real-world submitted clinical-significance labels:

- Pathogenic / Likely Pathogenic → P/LP
- Uncertain significance → VUS
- Benign / Likely Benign → B/LB

### Use review status as confidence

Review status should be used as a proxy for confidence and reliability:

- No assertion criteria provided → avoid or use only for sensitivity analysis
- Criteria provided, single submitter → usable, but lower confidence
- Criteria provided, multiple submitters, no conflicts → stronger
- Reviewed by expert panel → high-confidence / gold-standard

### Do not assume ClinVar contains complete ACMG evidence tables

For ACMG criterion-level data, better sources include:

1. ClinGen Evidence Repository, which is better suited for expert-curated variant classifications and evidence.
2. ClinVar submitter comments and XML, which can be mined but are incomplete and noisy.
3. The WE-LLM pipeline, which can extract or propose ACMG criteria from ClinVar, ClinGen, literature, gnomAD, and functional evidence.
4. Manual adjudication by a cardiogeneticist or expert reviewer for selected variants, which would serve as a gold-standard subset for difficult cases.

## Bottom line

ClinVar contains ACMG-based classifications and sometimes ACMG criteria, but it is not a complete structured ACMG-code database.

For this project, the intended framing is:

> ClinVar will be used as the primary public source of submitted clinical-significance labels and review-status metadata. Variant-level ACMG/AMP evidence codes will not be assumed to be uniformly available from ClinVar; instead, criterion-level evidence will be extracted where available from submitter comments and XML and supplemented using ClinGen Evidence Repository records, source-grounded WE-LLM evidence curation, and expert review for selected variants.

## Project intent

The goal of this project is to build a CatBoost-based model that learns from clinically meaningful variant features and labels derived from publicly available classification sources. The model will be trained to predict clinical significance categories while accounting for the fact that labels are noisy, partially observed, and sometimes conflicting.

The central design choice is to use ClinVar as the primary label source, while treating ACMG evidence as incomplete and potentially supplementary rather than uniformly available. This approach is realistic for a public-data project and reflects the actual structure of available evidence in the field.

## Protein-structure and stability feature plan

Protein-structure features will be used as model features, not as labels. The
first-pass structure layer will use reviewed UniProt mapping, AlphaFold residue
confidence, UniProt domain/site annotations, DSSP secondary structure, and
FreeSASA solvent accessibility. These features are useful because many
cardiogenetic missense variants act through disruption of conserved domains,
buried structural cores, channel regions, sarcomeric motifs, or other protein
contexts that are not fully captured by sequence-only predictors.

Delta-delta-G features will be added only for missense variants after the
variant-to-protein mapping is stable. The recommended study design is:

1. Use FoldX as the primary first-pass delta-delta-G engine if academic
   licensing is approved. It is fast, local, PDB-based, and practical for
   thousands of variants across AlphaFold structures.
2. Use Rosetta `cartesian_ddg` as a smaller sensitivity/validation analysis,
   not as the default whole-cohort engine. It is more computationally expensive
   but useful for checking whether FoldX-derived stability effects are robust
   in high-priority genes or variants.
3. Treat DynaMut2, PremPS, DDGun, and similar web/server tools as secondary
   comparison or spot-check resources unless their bulk/API terms and
   reproducibility are explicitly confirmed.

All stability outputs must include `ddg_method`, `ddg_status`, and
`ddg_missing_reason`. Delta-delta-G values from different engines must not be
merged as a single numeric feature without preserving the method identity,
because tools differ in calibration, input assumptions, and sign conventions.

## Dataset integration, leakage control, and sparse-gene strategy

The study will use a unified variant-level feature table built around a fixed
genome assembly and explicit source flags. GRCh38 is the preferred coordinate
system because ClinVar, gnomAD v4, current dbNSFP, AlphaMissense, and the local
protein-structure feature workflow are all easiest to maintain on GRCh38. Any
legacy GRCh37-derived source, including CardioBoost-derived coordinates, must be
lifted over or otherwise mapped with a documented status field before it is
joined to the main table.

The first implementation milestone is a master variant registry: one row per
`chrom-pos-ref-alt` coordinate key, with source flags and a companion
source-membership table. This registry is the spine for all later feature joins
and train/test split control. It is intentionally not the final modeling table;
it preserves provenance first, then downstream scripts can apply stricter label,
gene-panel, review-status, and class-balance rules.

### Dataset roles

| Dataset/source | Join key | What it adds | Primary role |
|---|---|---|---|
| ClinVar bulk | `chrom`, `pos`, `ref`, `alt`, gene | Public P/LP, VUS, B/LB labels; review status; submitter/conflict metadata | Main public label source after confidence filtering |
| gnomAD v4 | `chrom`, `pos`, `ref`, `alt` and gene | Allele frequency, popmax AF, homozygote/hemizygote counts, constraint metrics, benign-proxy candidates | Feature source and presumed-benign source |
| dbNSFP | `chrom`, `pos`, `ref`, `alt` | REVEL, SIFT, PolyPhen, MetaLR, FATHMM-XF, conservation, constraint, and related in silico predictors | Main in silico feature join |
| AlphaMissense | gene/protein position/alternate amino acid and genomic key when available | Missense pathogenicity score | Feature source |
| AlphaFold/UniProt/DSSP/FreeSASA/FoldXPro | UniProt accession and residue position | pLDDT, protein position, domains, motifs, topology, secondary structure, solvent accessibility, and delta-delta-G | Protein-structure feature family |
| ClinGen | gene, disease, HGVS/variant identifiers where available | Expert curation flags, gene-disease validity, dosage, VCEP classification context | High-confidence evidence feature and sensitivity stratum |
| HiRO / CASPER WES / VERDICT | Existing cleaned local table plus genomic/protein identifiers | Expert labels and linked patient-level phenotype fields | Highest-trust local stratum; phenotype-aware analysis only |
| eMERGE / Glazer arrhythmia | Existing cleaned local table plus variant identifiers | Real public variant-level labels, including many VUS | Core public variant source or held-out validation stratum depending on split design |
| CardioBoost coordinate-rich public handoff | HGVS/genomic identity | Public B/P labels and precomputed/reannotated feature fields | Binary B/P training source or held-out CardioBoost-comparison stratum |
| CardioBoost V2 / DYNA Zenodo | Protein-sequence fields | Protein-sequence B/P/VUS examples, not coordinate-rich | Auxiliary sequence-level analysis; not a replacement for coordinate joins |

HiRO is not valuable because it adds a large number of variants. It is valuable
because it is the only current source with linked patient-level phenotype fields
and blinded expert re-adjudication. Those variants should be treated as the
highest-trust stratum and analyzed with a `phenotype_available` flag. Phenotype
features should not be mixed into the main public-only model unless there is a
separate phenotype-aware experiment, because ClinVar, CardioBoost, eMERGE, and
gnomAD-derived rows do not carry comparable patient-level phenotype fields.
For bookkeeping, HiRO should be counted as private/patient-linked source rows
first, then mapped into the deduplicated coordinate registry second; repeated
variants and rows lacking complete `chrom-pos-ref-alt` keys must not be mistaken
for missing HiRO observations. If a HiRO row overlaps ClinVar by genomic
coordinate, that overlap should be treated as a shared variant identity only;
the HiRO label, patient context, and adjudication remain private HiRO
provenance.

### Feature annotation pipeline for new variants

Every ClinVar or gnomAD-derived variant that is not already represented in a
processed source table must be passed through a reproducible annotation
pipeline. The minimum feature pipeline is:

1. Normalize coordinates and alleles on GRCh38.
2. Deduplicate by `chrom-pos-ref-alt` with source-specific label provenance.
3. Join gnomAD v4 frequency and constraint features.
4. Join dbNSFP in silico predictors after the local file passes checksum QC.
5. Join AlphaMissense for missense variants.
6. Add VEP or equivalent consequence/transcript annotations where needed,
   especially for splice and consequence harmonization.
7. Add protein features from UniProt and AlphaFold for variants with reliable
   protein HGVS/residue mapping.
8. Add DSSP/FreeSASA structural context and FoldXPro delta-delta-G only for
   eligible missense variants with trustworthy protein mapping.

The old CardioBoost feature space should be used as a reproducibility guide,
not as a hard requirement. The most important near-term feature set is a
competitive public-data feature table combining dbNSFP, gnomAD, AlphaMissense,
ClinGen, and protein-structure features. UCSC Multiz/orthologous-alignment
features are a stretch goal because they are harder to reproduce and should not
block the main model.

### Leakage-aware analysis design

CardioBoost-derived variants require explicit leakage control. If the model is
trained on CardioBoost variants, it cannot be presented as an independent test
against CardioBoost's published performance. Therefore the study should define
two analyses:

1. Primary comparison analysis: hold out CardioBoost-derived coordinate-rich
   variants as an external binary B/P test stratum, train on non-CardioBoost
   public/internal sources, and evaluate generalization.
2. Expanded training analysis: include CardioBoost variants in training and
   evaluate on separate held-out ClinVar/eMERGE/HiRO-derived splits with
   source-stratified reporting.

All train/test splits must be deduplicated by genomic variant identity and, for
protein-only auxiliary datasets, by gene/protein change where genomic identity
is unavailable. Rows derived from the same underlying variant must not appear on
both sides of a split under different source names.

### Class balance and VUS sampling

**Updated July 2 training-slice decision:** the first serious model should be
trained as a binary pathogenic-vs-benign classifier, not as a primary 3-class
Benign/Pathogenic/VUS model. The 3-class model should still be tried, but as an
exploratory/sensitivity analysis rather than the main result.

The key reason is conceptual: **VUS is not a true biological class**. A VUS
often means "not enough evidence yet," not "intermediate pathogenicity." Some
VUS are future benign variants, some are future pathogenic variants, and many
are simply under-evidenced. If the model is trained to predict VUS as a third
class, it may learn human uncertainty, submitter behavior, and evidence gaps
rather than variant effect. The stronger paper story is to train a calibrated
binary P/LP-vs-B/LB model and then apply it to VUS rows to prioritize variants
that look more pathogenic-like or benign-like.

After conservative basic filters:

- include primary rows only;
- remove label conflicts;
- remove `ref == alt` no-op rows;
- keep clean labels only;

the current HiRO-aggregation matrix gives the following sample-size tradeoff:

| Training slice | Rows | Benign | Pathogenic | VUS |
|---|---:|---:|---:|---:|
| Clean 3-class available | `86,028` | `33,978` | `9,223` | `42,827` |
| Clean binary only | `43,201` | `33,978` | `9,223` | `0` |
| Clean VUS only | `42,827` | `0` | `0` | `42,827` |
| ClinVar clean binary | `42,862` | `33,811` | `9,051` | `0` |
| ClinVar stars >=2 binary | `15,190` | `11,965` | `3,225` | `0` |
| ClinVar star 1 binary | `27,672` | `21,846` | `5,826` | `0` |
| Non-ClinVar source binary | `339` | `167` | `172` | `0` |
| HiRO binary | `70` | `43` | `27` | `0` |
| eMERGE binary | `203` | `90` | `113` | `0` |
| CardioBoost binary | `225` | `82` | `143` | `0` |

Thus, training only on high-confidence ClinVar stars `>=2` would reduce the
binary training pool from about `43.2k` rows to about `15.2k` rows. That is
still enough for a sensitivity model, but it is much smaller and less broad.
The primary model should therefore include all clean binary rows while using
sample weights to reflect source/confidence quality.

Recommended model formulations:

| Model | Training data | Role |
|---|---|---|
| Primary model | Binary P/LP vs B/LB, all clean rows, weighted by confidence | Main classifier |
| High-confidence sensitivity model | Binary P/LP vs B/LB, ClinVar stars >=2 plus trusted sources | Robustness check |
| Exploratory 3-class model | Benign / Pathogenic / VUS | Secondary experiment only |
| VUS scoring analysis | Apply binary model to VUS rows | Main VUS-facing clinical analysis |

For ClinVar review stars, do not simply discard all 1-star rows in the primary
model. Use weights instead, and exclude review stars from the biological feature
set so they do not leak label confidence into prediction.

Example weighting concept:

| Row type | Suggested weight concept |
|---|---:|
| HiRO expert/private patient-linked rows | High |
| ClinVar 3-star rows | High |
| ClinVar 2-star rows | Medium-high |
| ClinVar 1-star rows | Lower |
| CardioBoost/eMERGE rows | Source-specific and split-controlled |
| Label conflicts | Exclude from supervised training |
| VUS rows | Exclude from binary training; score later |

The first training slice should therefore be:

```text
binary labels only: Benign vs Pathogenic
primary_model_inclusion == include
label_conflict_type == none
ref != alt
exclude VUS from supervised binary training
include ClinVar star 1 rows, but with lower sample weight
```

The high-confidence sensitivity model should be:

```text
same binary task
ClinVar stars >=2 only, plus trusted non-ClinVar rows
```

Then compare the all-clean weighted model against the high-confidence-only
model. If performance and calibration are similar, the all-clean weighted model
is preferred because it uses broader evidence and more genes. If the
high-confidence-only model is more stable, it can become the conservative
primary result.

VUS rows should be handled as follows:

| VUS use | Recommendation |
|---|---|
| Training target | Not primary |
| Model scoring | Yes; very important |
| Calibration/evidence prioritization | Yes |
| Case studies | Yes |
| Future reclassification analysis | Yes |
| 3-class experiment | Exploratory only |

The preferred manuscript framing is:

> We trained a binary pathogenic-vs-benign classifier using high-confidence and
> weighted lower-confidence labels, then applied it to VUS to prioritize
> variants likely to behave like pathogenic or benign variants.

That is more defensible than claiming the primary model "predicts VUS," because
VUS is an evidence status, not a stable biological endpoint.

The next modeling artifact should be a training-slice table with columns such
as:

```text
training_slice_primary_binary
training_slice_high_confidence_binary
training_slice_three_class_exploratory
sample_weight
split_group
exclude_reason
```

The older note below remains useful for a possible exploratory 3-class/VUS
analysis, but it is no longer the recommended primary model formulation.

The expected bottleneck is not raw VUS availability. ClinVar and the DYNA
protein-level files provide many more VUS examples than can be used naively.
The modeling dataset should therefore use stratified VUS sampling rather than
letting the VUS class dominate training.

Recommended sampling rule:

- Build P/LP, B/LB, and VUS pools per gene and consequence class.
- Exclude conflicting ClinVar labels from the primary training set.
- Keep high-review-status P/LP and B/LB variants preferentially.
- Sample VUS at a fixed ratio, such as at most 2x to 3x the number of P/LP
  variants per gene, while preserving gene and consequence diversity.
- Keep an unsampled VUS pool for calibration, sensitivity analysis, and
  out-of-distribution evaluation.

This makes the synthetic-data arm scientifically sharper: synthetic
augmentation should not be framed as merely fixing a small total-N problem. It
should be tested as a way to improve performance in sparse genes and sparse
mechanistic subgroups where real P/LP examples remain limited.

### Sparse-gene strategy

For genes with fewer than about 30 real P/LP variants after all real sources are
combined, the analysis will use three complementary strategies:

| Strategy | Use | Rationale |
|---|---|---|
| Disease/mechanism pooling | Pool related sparse genes, such as sarcomere HCM genes, when biology and VCEP logic support it | Lets the model learn shared mechanisms without pretending each sparse gene has enough independent signal |
| Sparse-gene reporting flag | Report metrics separately for sparse genes vs. well-represented genes | Directly tests whether model performance degrades in low-N genes |
| Targeted synthetic augmentation | Add synthetic/evidence-derived examples specifically for sparse genes or sparse mechanisms | Tests whether augmentation closes the sparse-gene performance gap rather than only improving aggregate metrics |

Sparse-gene performance should be a prespecified secondary analysis. If
augmentation helps primarily in sparse genes, that is a stronger scientific
finding than a small aggregate improvement in the full dataset.

## Current project status, implementation log, and what is left to do

This section is a durable handoff for future Codex sessions and collaborators.
It records the current workspace state after the dataset organization,
feature-engineering, source-record preservation, and background annotation work.
The important design point is that this project now keeps two related but
different data layers:

1. A variant-level modeling layer: one row per `chrom-pos-ref-alt`, used for
   CatBoost features, leakage control, and train/test splitting.
2. Source-record preservation layers: one row per original source record for
   HiRO, eMERGE, and CardioBoost, used to preserve patient/source-level
   classifications, phenotype fields, and provenance.

The variant-level modeling matrix must not be interpreted as the raw dataset
size for sources such as HiRO. HiRO in particular contains repeated genotypes
across different patient/source records, and those records can carry different
phenotype signatures and sometimes different classifications.

### Snapshot date and current best model-facing files

Current local status snapshot: 2026-07-02 in the local workspace.

Important note for future readers: the July 2 updates below supersede older
July 1 counts later in this section where they conflict. The older details are
kept because they document the path taken, but the current model-facing matrix
is now the gnomAD/FoldX/VEP/HGVS/SpliceAI matrix with HiRO source-record
aggregation features.

| Purpose | Current file |
|---|---|
| Clean variant registry after HiRO coordinate rescue promotion | `datasets/variant_registry/interim/clean_combined_variant_registry_with_hiro_rescues.tsv` |
| Source-membership table after HiRO rescue promotion | `datasets/variant_registry/interim/clean_combined_variant_registry_sources_with_hiro_rescues.tsv` |
| Frozen local-feature matrix before HiRO rescue application | `datasets/modeling/interim/final_modeling_table_local_features_frozen.tsv` |
| Best current local-feature matrix, rescue-aware | `datasets/modeling/interim/final_modeling_table_local_features_frozen_with_hiro_rescues.tsv` |
| Current best population/stability/consequence matrix before HGVS/SpliceAI | `datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx_vep.tsv` |
| VEP + HGVS + SpliceAI matrix before HiRO aggregation | `datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx_vep_hgvs_spliceai.tsv` |
| Current best model-facing matrix with HiRO source-record aggregation | `datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx_vep_hgvs_spliceai_hiro_agg.tsv` |
| HiRO source-record aggregation feature table | `datasets/modeling/interim/hiro_source_record_aggregation_features.tsv` |
| HiRO source-record-level modeling/evaluation table | `datasets/modeling/interim/hiro_source_record_level_table.tsv` |
| HiRO full source-record preservation table | `datasets/hiro/full_dataset/interim/hiro_linked_source_records.tsv` |
| eMERGE full source-record preservation table | `datasets/emerge/full_arrhythmia_gene_dataset/data/emerge_linked_source_records.tsv` |
| CardioBoost full source-record preservation table | `datasets/cardioboost/public_dataset/processed_for_cardiogenetics/model_inputs/cardioboost_linked_source_records.tsv` |

### July 2 implementation update: gnomAD, FoldX, VEP, HGVS, and SpliceAI

This update records the current state after the gnomAD browser extraction,
FoldX full batch, VEP installation/debugging, VEP consequence annotation, and
the start of the HGVS/SpliceAI VEP rerun.

#### Current best model-facing matrix

The current best completed matrix is:

`datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx_vep.tsv`

| Item | Current value |
|---|---:|
| Rows | `86,889` |
| Columns | `403` |
| Duplicate `variant_id`s | `0` |
| gnomAD final status: observed | `54,173` |
| gnomAD final status: confirmed absent | `17,275` |
| gnomAD final status: not joined | `15,441` |
| gnomAD AF nonmissing | `64,472` |
| FoldX selected variants joined | `1,151` |
| FoldX selected `ok` variants | `1,099` in the merged matrix |
| FoldX selected `ref_mismatch` variants | `52` |
| VEP annotated variants | `86,872` |
| VEP missing variants | `17` |

The immediate predecessor matrix
`datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx.tsv`
merged gnomAD browser/dbNSFP population features and FoldX DDG features. The
current VEP matrix adds consequence, transcript, protein-position, impact,
canonical/MANE, and variant-class features.

#### gnomAD population features

Direct gnomAD browser/Hail extraction is now complete and parsed. It was used
only for variants missing useful dbNSFP-derived gnomAD coverage, rather than
downloading all raw gnomAD.

| gnomAD item | Count / status |
|---|---:|
| Browser/Hail variants extracted | `44,548 / 44,548` |
| Browser status `ok` | `27,273` |
| Browser status `not_found` | `17,275` |
| Parsed browser feature columns | `103` |
| Matrix rows with final gnomAD status `observed` | `54,173` |
| Matrix rows with final gnomAD status `confirmed_absent` | `17,275` |
| Matrix rows with final gnomAD status `not_joined` | `15,441` |

Safe interpretation rule: `confirmed_absent` means the exact variant was queried
and not found in the gnomAD browser output; these rows can receive AF=0 only
with `gnomad_confirmed_absent_flag=true`. `not_joined` is not equivalent to
AF=0 and must remain missing/unknown.

#### FoldX DDG features

The FoldX primary batch completed and was merged into the modeling matrix with
explicit status/missingness flags.

| FoldX item | Count |
|---|---:|
| Raw batch rows | `1,230` |
| Selected variant-level FoldX rows | `1,152` |
| Selected `ok` variant-level rows | `1,100` |
| Selected `ref_mismatch` variant-level rows | `52` |
| `ok` DDG values in current merged matrix | `1,099` |
| FoldX duplicate variant groups | `17` |
| FoldX discordant duplicate groups | `2` |

FoldX mismatch QC is complete:

| Mismatch QC item | Count / interpretation |
|---|---|
| Raw `ref_mismatch` rows | `61` |
| Selected `ref_mismatch` variants | `52` |
| Largest mismatch cluster | `CACNA1C` / UniProt `Q13936`, `37` raw rows and `36` selected variants |
| Source distribution | eMERGE `32`, HiRO `14`, CardioBoost `6` selected mismatch variants |
| Primary rule | Treat all `ref_mismatch` DDG values as missing; keep `foldx_ddg_ref_mismatch_flag=true`; do not impute DDG |

The likely cause of most mismatches is transcript/isoform or protein-position
incompatibility, especially in CACNA1C. These should be remapped with a
canonical transcript/protein workflow before any rerun, not silently rescued.

#### VEP consequence annotation

VEP setup required several engineering fixes:

1. The downloaded Ensembl VEP cache was successfully unpacked locally.
2. The native macOS arm64 conda VEP build repeatedly segfaulted before
   annotation.
3. Rosetta 2 was installed, and a separate `vep-116-x86` conda environment was
   created using `osx-64` packages.
4. Direct `arch -x86_64` VEP execution worked under Rosetta.
5. A consequence-first VEP run was completed without HGVS or SpliceAI.

Completed VEP consequence output:

| VEP item | Count / file |
|---|---|
| Raw VEP transcript/consequence rows | `2,605,524` |
| Unique variants returned by VEP | `86,872` |
| Feature rows after collapse | `86,889` |
| VEP feature columns merged | `27` |
| Completed VEP matrix | `datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx_vep.tsv` |

Top current VEP consequence counts in the completed matrix:

| Consequence | Variants |
|---|---:|
| `missense_variant` | `40,280` |
| `synonymous_variant` | `17,080` |
| `intron_variant` | `11,088` |
| `frameshift_variant` | `4,276` |
| `splice_region_variant` | `3,192` |
| `stop_gained` | `2,345` |

The `17` VEP-missing rows have been QC'd. All are ClinVar rows where
`ref == alt`, such as `A>A`, `G>G`, or `AGG>AGG`. These are no-op/non-variant
records, not true VEP failures. Audit file:

`datasets/feature_sources/insilico/interim/vep/vep_missing_variant_qc.tsv`

Primary handling rule: keep source provenance, but flag/exclude these no-op
alleles from variant-effect annotation and primary modeling unless there is a
specific reason to keep them in a source-record-only audit.

#### HGVS and SpliceAI rerun in progress

A richer VEP run is currently in progress in detached screen session
`vep_hgvs_spliceai_x86`.

Purpose:

1. Add `HGVSc` and `HGVSp` from VEP.
2. Add SpliceAI SNV scores from the local Ensembl MANE GRCh38 SNV file.
3. Preserve the previous completed consequence-only matrix by writing to a new
   output suffix.

Current in-progress target:

`datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx_vep_hgvs_spliceai.tsv`

Current state at the time of this proposal update:

| Item | Status |
|---|---|
| GRCh38 primary assembly FASTA | Downloading in background; expected compressed size `881,964,081` bytes |
| FASTA download behavior | Resumable via `curl --continue-at -` |
| FASTA destination | `datasets/feature_sources/insilico/raw/vep/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz` |
| SpliceAI SNV file | Local, `27G`, with `.tbi` index |
| SpliceAI plugin | Locally patched to allow SNV-only use because the public Ensembl GRCh38 MANE file is SNV-only |
| Monitor command | `tail -f datasets/feature_sources/insilico/raw/vep/logs/vep_cache_install_annotate.log` |

If this screen stops before producing the final HGVS/SpliceAI matrix, rerun:

```bash
screen -dmS vep_hgvs_spliceai_x86 bash -lc 'VEP_ENV_NAME=vep-116-x86 VEP_CONDA_SUBDIR=osx-64 VEP_ENABLE_HGVS=1 VEP_ENABLE_SPLICEAI=1 VEP_RUN_SUFFIX=vep_hgvs_spliceai ./datasets/feature_sources/insilico/scripts/run_vep_cache_install_annotate.sh'
```

The rerun should resume the FASTA download if it was interrupted, then
decompress and index the FASTA, smoke-test VEP, run full VEP, parse features,
and merge the HGVS/SpliceAI matrix.

### Major work completed so far

1. Organized the project around a `datasets/` directory with separate areas for
   raw, interim, processed, feature-source, and modeling outputs.
2. Read and interpreted the study proposal, then aligned the dataset plan around
   ClinVar, gnomAD, dbNSFP, AlphaMissense, ClinGen, HiRO, eMERGE, CardioBoost,
   and protein-structure features.
3. Downloaded and validated ClinVar bulk `variant_summary.txt.gz`; curated it
   into a germline cardiogenetics slice and then a primary gene-QC slice.
4. Audited the cardiogenetics gene panel and separated primary genes,
   sensitivity genes, TTN-separate handling, and likely contamination.
5. Built a clean combined variant registry from ClinVar, HiRO, eMERGE, and
   CardioBoost using `chrom-pos-ref-alt` as the coordinate identity.
6. Built source-membership tables so each variant-level row retains source
   flags and source-label provenance.
7. Downloaded dbNSFP v5.3.1a GRCh38, validated the clean download, installed
   indexed-query support through the local Python/Hail environment, and joined
   selected dbNSFP features.
8. Joined local dbNSFP-derived gnomAD fields and created a plan to use direct
   gnomAD browser/Hail only for variants missing useful dbNSFP gnomAD coverage.
9. Downloaded/validated AlphaMissense and built a direct AlphaMissense join
   independent of dbNSFP.
10. Queried/processed ClinGen gene validity, dosage, and variant-level evidence
    into selected model-facing ClinGen features.
11. Built UniProt/AlphaFold/protein-position parsing features for local sources.
12. Installed and ran DSSP/mkdssp and FreeSASA for eligible AlphaFold residues.
13. Installed/smoke-tested FoldXPro, ran a pilot, then started the full primary
    FoldX batch for eligible high-confidence missense variants.
14. Started VEP/SpliceAI cache download in a resumable background screen.
15. Started gnomAD browser/Hail chunk extraction in a resumable background
    screen after installing Java, Hail, gcloud, and Google authentication.
16. Promoted high-confidence HiRO coordinate rescues without overwriting the
    original registry.
17. Created linked source-record tables for HiRO, eMERGE, and CardioBoost to
    avoid losing patient/source-level information when modeling rows are
    deduplicated by genotype.

### Dataset status

| Source | Current local status | Current count | Primary role | Remaining work |
|---|---|---:|---|---|
| ClinVar bulk | Local, validated, curated, gene-QC filtered | `83,915` primary rows after gene-panel QC | Main public B/LB, VUS, P/LP label source | Preserve confidence fields; keep holdout/sensitivity/TTN strata separate |
| HiRO / CASPER WES / VERDICT | Local, source-record table preserved, coordinate rescues promoted | `482` raw HiRO source records; `408` linked to promoted registry; `73` unresolved; `1` rescued but not promoted | Highest-trust private/patient-linked phenotype and adjudication stratum | Resolve the remaining unresolved rows; build safe aggregate phenotype/source features |
| eMERGE / Glazer arrhythmia | Local, profiled, linked source-record table built | `2,754` source records; `2,754` unique linked variants | Public arrhythmia-gene source, VUS-heavy validation/training stratum | Decide training vs held-out/source-stratified evaluation role |
| CardioBoost coordinate-rich public handoff | Local, processed, linked source-record table built | `355` source records; `353` unique variants | Public binary B/P cardiogenetics source | Use with leakage-aware splits; do not use as independent benchmark if trained on it |
| CardioBoost V2 / DYNA Zenodo | Local and processed | `21,472` protein-sequence rows | Auxiliary sequence-level analysis | Keep as extra only; it lacks coordinate completeness and does not replace original CardioBoost |
| gnomAD v4 | Local dbNSFP-derived gnomAD features done; direct browser/Hail extraction running | `12 / 45` Hail chunks complete in current run | Population frequency, popmax AF, homozygote count, benign proxy support | Finish/repair chunk extraction; merge browser fields for dbNSFP-missing variants |
| dbNSFP v5.3.1a | Local, validated, selected features built and merged into frozen local matrix | `42,338` dbNSFP `ok` rows in current rescue-aware matrix | Main in silico feature block | Rerun feature joins for the 3 newly added rescue-aware variants later |
| AlphaMissense | Local direct join done | `36,773` direct exact matches | Independent missense pathogenicity feature | Rerun for newly appended rescue-aware rows later |
| ClinGen | Gene validity, dosage, and variant evidence selected/merged | `84,983` gene-validity `ok` rows; `319` exact variant evidence rows in earlier selected table | Expert evidence and context flags | Use as high-confidence evidence feature, not a broad replacement label source |
| Protein structure | UniProt, AlphaFold, DSSP, FreeSASA selected/merged | `3,089` protein-feature `ok` rows in current rescue-aware matrix | Protein context, domains/sites, pLDDT, secondary structure, solvent accessibility | Extend mapping only after primary matrix QC; merge FoldX when finished |

### Registry and modeling matrix status

The registry is the provenance spine. The modeling matrix is the feature table.
They should remain related but not confused.

| Artifact | Rows | Notes |
|---|---:|---|
| `clean_combined_variant_registry.tsv` | `86,886` | Pre-HiRO-rescue clean coordinate registry |
| `clean_combined_variant_registry_with_hiro_rescues.tsv` | `86,889` | Rescue-aware registry; 3 new unique variant rows added |
| `clean_combined_variant_registry_sources_with_hiro_rescues.tsv` | `87,434` | Source-membership table; preserves source labels/provenance |
| `final_modeling_table_local_features_frozen.tsv` | `86,886` | ClinGen + protein/structure + dbNSFP selected + direct AlphaMissense |
| `final_modeling_table_local_features_frozen_with_hiro_rescues.tsv` | `86,889` | Current best local-feature matrix; duplicate variant IDs checked: `0` |

Current rescue-aware matrix feature coverage:

| Feature/status field | Count |
|---|---:|
| Total rows | `86,889` |
| Duplicate `variant_id`s | `0` |
| Rows with `in_hiro=true` / private-linked | `240` variant-level rows |
| dbNSFP `ok` | `42,338` |
| dbNSFP `not_found` | `44,548` |
| dbNSFP `not_joined` | `3` newly appended HiRO rescue-aware rows |
| Direct AlphaMissense `ok` | `36,773` |
| Direct AlphaMissense `not_applicable_or_not_found` | `50,113` |
| Direct AlphaMissense `not_joined` | `3` newly appended HiRO rescue-aware rows |
| ClinGen gene validity `ok` | `84,983` |
| ClinGen gene validity `not_found` | `1,903` |
| Protein feature `ok` | `3,089` |
| Protein feature `not_available` | `83,562` |

### Source-record preservation and duplicate-genotype audit

This audit was added after identifying an important issue: HiRO has genotype,
ACMG/adjudication, and phenotype information. Repeated genotypes should be
collapsed for variant-level feature learning and leakage control, but the
underlying patient/source records must be preserved.

| Dataset | Raw/source records | Unique resolved variants | Duplicate genotype groups | Discordant labels among duplicates | Preservation table |
|---|---:|---:|---:|---:|---|
| HiRO | `482` | `239` resolved variant IDs in linked table | `46` | `4` | `datasets/hiro/full_dataset/interim/hiro_linked_source_records.tsv` |
| eMERGE | `2,754` | `2,754` | `0` | `0` | `datasets/emerge/full_arrhythmia_gene_dataset/data/emerge_linked_source_records.tsv` |
| CardioBoost | `355` | `353` | `2` | `0` | `datasets/cardioboost/public_dataset/processed_for_cardiogenetics/model_inputs/cardioboost_linked_source_records.tsv` |
| ClinVar primary slice | `83,915` | `83,915` | `0` | `0` | ClinVar curated table is already variant-level for the current slice |

HiRO source-record detail:

| HiRO status | Count |
|---|---:|
| Raw HiRO records | `482` |
| Coordinate-complete before rescue | `359` |
| Coordinate-incomplete before rescue | `123` |
| High-confidence coordinate rescues promoted | `49` |
| HiRO records linked to promoted registry in linked table | `408` |
| Unresolved no coordinates | `73` |
| Rescued but not promoted | `1` |
| Duplicate genotype groups with variable phenotype signature | `28` |
| Duplicate genotype groups with discordant 3-class labels | `4` |

Interpretation: `in_hiro=true` on 240 variant-level rows does not mean only 240
HiRO observations exist. It means 240 unique coordinate variants in the current
matrix have at least one HiRO source record attached. The full HiRO evidence
remains in `hiro_linked_source_records.tsv`.

### HiRO coordinate rescue status

The first HiRO rescue pass used HGVS-like cDNA/protein terms against ClinVar
variant summary. Only high-confidence matches were promoted.

| Rescue category | Count |
|---|---:|
| HiRO rows reviewed for rescue | `123` |
| High-confidence rescues | `50` candidates |
| High-confidence rescues promoted | `49` |
| Existing registry variants enriched | `46` |
| New unique variant rows added | `3` |
| Review candidates not promoted | `3` |
| Not rescued | `70` |

The promoted rows are auditable in:
`datasets/hiro/full_dataset/interim/hiro_coordinate_rescues_promoted.tsv`.

### Gene-panel QC status

The current gene-panel QC layer separates primary modeling genes from
contamination, alias/quarantine genes, TTN-separate handling, and
sensitivity-only genes.

| Bucket | Genes | Current action |
|---|---|---|
| Remove or quarantine | `AKAP9`, `DMD`, `FPGT`, `FPGT-TNNI3K`, `KCNE1B`, `KNCH2` | Exclude from primary panel |
| Sensitivity or separate analysis | `ANK2`, `SOS1`, `TNNI3K`, `TRPM4`, `TTN` | Exclude from the primary ClinVar training slice or handle separately |
| Conditional phenocopy or expanded-scope genes | `EMD`, `GLA`, `HRAS`, `KRAS`, `LAMP2`, `MAP2K2`, `NRAS`, `PTPN11`, `RAF1`, `RIT1`, `SLC22A5` | Retain with conditional/syndromic-overlap flags |
| Keep but monitor | `CACNB2`, `FHOD3`, `MIB1`, `SCN1B`, `SNTA1` | Retain but monitor evidence strength and feature contribution |

Current post-QC ClinVar slices:

| Slice | Rows | Benign | Pathogenic | VUS |
|---|---:|---:|---:|---:|
| Primary after gene QC | `83,915` | `33,835` | `9,052` | `41,028` |
| Holdout/sensitivity | `44,052` | `20,935` | `5,934` | `17,183` |
| TTN separate review | `34,804` | `17,203` | `5,765` | `11,836` |

### In silico feature status

| Feature source | Current status | What is done | What remains |
|---|---|---|---|
| dbNSFP v5.3.1a | Downloaded, MD5-validated, selected-feature table built, merged into local frozen matrix | Selected columns include REVEL, SIFT, PolyPhen2 HDIV, MetaLR, FATHMM-XF, CADD, conservation, ESM1b, AlphaMissense, and gnomAD/dbNSFP popmax fields | Rerun selected-feature build after rescue-aware registry becomes the default source |
| Direct AlphaMissense | Local direct join done and merged | `36,773` exact matches in current matrix | Rerun for newly appended rescued variants and later final registry |
| CADD | Available through dbNSFP for covered variants; targeted API path documented | CADD-like fields are present from dbNSFP for covered variants | Decide whether separate CADD API/bulk adds enough beyond dbNSFP |
| VEP | Cache download running in background | Cache file reached about `18.30 / 27.64 GB` (`66.2%`) at latest check | Finish download, install/test VEP, annotate consequence/transcript/splice features |
| SpliceAI | Pending through VEP/plugin or verified precomputed source | Requirement documented | Merge after VEP/cache setup |
| ESM1b / ESM-variant | Present through dbNSFP for covered variants | dbNSFP ESM1b fields selected | Confirm whether a separate direct ESM source is worth adding later |

dbNSFP coverage in the current matrix is not expected to be 100%. The `ok`
matches are mostly SNV/protein-coding/scored variants. Missing rows include
indels, MNVs, noncoding variants, synonymous/UTR-like rows, and variants outside
dbNSFP's scored scope. This is not currently treated as a failed join.

### gnomAD status

The current strategy is not to download full raw gnomAD unless necessary. The
project uses dbNSFP-derived gnomAD fields immediately and direct gnomAD
browser/Hail extraction only for variants missing dbNSFP coverage.

Current gnomAD direct extraction status:

| Item | Status |
|---|---|
| Local dbNSFP-derived gnomAD fields | Done for dbNSFP-covered variants |
| gnomAD gene constraint table | Available locally for the 53-gene panel |
| Browser/Hail chunk input | `45` chunks total |
| Completed browser/Hail chunks | `12 / 45` |
| Current running chunk at last check | `chunk_0012` |
| Current concern | Hail process appears long-running/stalled on `chunk_0012`; if it does not complete, switch to smaller chunks or one-chunk defensive extraction |

Remaining gnomAD tasks:

1. Finish or repair direct gnomAD browser/Hail extraction.
2. Combine chunk outputs into one gnomAD browser feature table.
3. Parse useful browser fields, including AF, popmax AF, homozygote/hemizygote
   counts, filters, and transcript/consequence fields if available.
4. Merge browser-derived fields with dbNSFP-derived gnomAD fields.
5. Define presumed-benign proxy rules separately from ClinVar B/LB labels.

### ClinGen status

ClinGen is now a mostly completed local feature block for gene-level context,
with limited but valuable variant-level evidence.

| ClinGen feature block | Current status |
|---|---|
| Gene-disease validity | Selected features built and merged |
| Dosage sensitivity | Selected features built and merged |
| Variant-level / VCEP evidence | Exact matching completed where available; limited coverage |
| Selected feature table | `datasets/modeling/interim/clingen_selected_features.tsv` |
| Merged matrix with ClinGen | Already included in local frozen matrix |

Interpretation: ClinGen should be used as expert evidence/context and
sensitivity-stratum metadata, not as a broad replacement source of labels.
Variant-level exact evidence is high-confidence but sparse.

### Protein, AlphaFold, DSSP, FreeSASA, and FoldX status

| Component | Current status | Remaining work |
|---|---|---|
| UniProt reviewed human FASTA/features | Local | Keep version pinned and mapping decisions documented |
| AlphaFold human structures | Local | Improve fragmented-protein mapping for large proteins later |
| Protein-position parser | Built for HiRO/eMERGE/CardioBoost HGVS-like fields | Extend only after primary matrix QC |
| UniProt domains/sites/motifs | Built and merged | Keep functional flags and domain/site names |
| AlphaFold pLDDT | Built and merged where mappable | Improve unresolved fragment mapping later |
| DSSP/mkdssp | Installed, batch features built, merged | None urgent |
| FreeSASA | Installed, batch features built, merged | None urgent |
| FoldXPro smoke test | Done | None |
| FoldXPro pilot | Done: 9 variants, 7 ok, 2 ref/isoform mismatches | Use QC rules |
| FoldX full primary batch | Running | Merge only after completion/QC |

Latest FoldX batch status at the time of this update:

| FoldX item | Count |
|---|---:|
| Primary eligible variants selected | `1,230` |
| Batch rows attempted so far | `286` |
| `ok` DDG results | `244` |
| `ref_mismatch` | `42` |

The mismatch count is important. FoldX should be merged with explicit
`ddg_status`, `ddg_missing_reason`, and eligibility flags. Do not silently
impute failed/mismatched DDG values.

### Background jobs currently running

As of the July 2 update, the earlier FoldX, gnomAD, and VEP-cache jobs have
completed or been superseded by merged outputs. The only current live screen job
expected at the time of this update is the richer VEP HGVS/SpliceAI rerun.

| Screen/job | Current status | Monitor |
|---|---|---|
| `vep_hgvs_spliceai_x86` | Running; downloading/indexing GRCh38 FASTA, then rerunning VEP with HGVS + SpliceAI SNV | `datasets/feature_sources/insilico/raw/vep/logs/vep_cache_install_annotate.log` |

Use `screen -ls` to check live sessions. If the HGVS/SpliceAI VEP rerun stops,
rerun the documented `screen -dmS vep_hgvs_spliceai_x86 ...` command in the
July 2 VEP section above.

### Current priority next steps

| Priority | Task | Status | Why |
|---:|---|---|---|
| 1 | Let HGVS/SpliceAI VEP rerun finish | Running in `vep_hgvs_spliceai_x86` | Adds `HGVSc`, `HGVSp`, and SpliceAI SNV features on top of completed VEP consequence fields |
| 2 | Validate and merge the HGVS/SpliceAI VEP rerun | Waiting on current screen job | Confirm row counts, SpliceAI coverage, HGVS coverage, and no accidental overwrite of the completed VEP matrix |
| 3 | Generate a full QC report for the current best completed matrix | Ready now using `final_modeling_table_local_features_gnomad_foldx_vep.tsv` | Shows model-ready rows, label balance, feature missingness, source-record preservation, no-op allele flags, leakage/duplicate checks, and gene-stratified counts |
| 4 | Decide whether the first baseline model uses consequence-only VEP or waits for HGVS/SpliceAI | Pending current rerun | Consequence-only VEP is already usable; HGVS/SpliceAI may improve splice/HGVS harmonization but should not block all QC |
| 5 | Build safe source-record aggregation features | Needed before modeling | Add variant-level summaries such as source record count, source label discordance, HiRO phenotype evidence availability, and private/source-linked flags |
| 6 | Define final modeling slices and split keys | Needed before modeling | Freeze primary/sensitivity/TTN/no-op allele exclusions, source-stratified splits, and leakage groups |
| 7 | Build first baseline CatBoost-ready table | Next modeling milestone | Convert the current feature matrix into a clean train/validation/test table with selected features, labels, missingness indicators, and split assignments |
| 8 | Optional: improve protein/HGVS mapping and rerun FoldX mismatches | Deferred | Only needed if structural features are to be a primary claim; current mismatch flags are safe for baseline modeling |
| 9 | Optional: separate direct ESM/other in silico sources beyond dbNSFP | Deferred | dbNSFP already provides ESM1b and many in silico predictors for covered rows |

### Work that can be deferred until after the first model

These are valuable but should not block a first baseline model:

1. Separate direct ESM1b/ESM-variant source beyond dbNSFP.
2. Rosetta `cartesian_ddg` sensitivity runs.
3. UCSC Multiz or custom orthologous-alignment features.
4. More aggressive gnomAD benign-proxy sampling.
5. Phenotype-aware HiRO-only modeling experiments.
6. Full-scale protein mapping expansion to all ClinVar/gnomAD rows beyond the
   currently parsed local source records.

### Main risks to control

| Risk | Control |
|---|---|
| Source leakage across train/test splits | Split by variant identity and source membership; do not allow the same coordinate variant on both sides |
| Source-record collapse hiding patient evidence | Preserve linked source-record tables; aggregate explicitly only when intended |
| ClinVar label noise | Use review stars, conflict filters, source labels, and sensitivity analyses |
| VUS class domination | Use stratified VUS sampling rather than all available VUS rows |
| TTN overwhelming the dataset | Keep TTN separate or consequence-restricted |
| Non-cardiac/syndromic contamination | Preserve gene-panel QC buckets and conditional phenocopy flags |
| Silent feature missingness | Require explicit status and missing-reason columns for every feature source |
| Mixing phenotype and non-phenotype rows | Keep phenotype-aware HiRO analyses separate from public variant-only models |
| Protein mapping errors | Require UniProt, residue, wild-type amino-acid, AlphaFold-fragment, and FoldX reference-AA QC before structure/DDG features |
| Treating gnomAD benign proxies as labels | Keep gnomAD proxy benign candidates separate from ClinVar/eMERGE/CardioBoost/HiRO labels |


thinking about the study and ideas to consider and reivew of the current status: 

These are exactly the right questions to ask before you commit to the full synthetic arm. Let me go through each one carefully.

---

## 1. With 83,915 ClinVar rows (post-QC), do you even need synthetic data?

The honest answer is: **probably not for overall performance, but potentially yes for specific subproblems.**

Here is why the answer is nuanced rather than a simple no.

**What 83,915 rows actually means for training.** After your VUS subsample (2–3x P/LP per gene), your training set will be roughly:

- P/LP: ~9,000 (your post-QC count)
- B/LB: ~33,000 (you have more than enough, will also subsample to balance)
- VUS: ~18,000–27,000 (2–3x P/LP)
- **Practical training N: ~40,000–50,000 variants**

That is 35–44x CardioBoost's training set. At that scale, the marginal benefit of adding synthetic data to the overall model is likely small — you are well past the point where raw N is your binding constraint for the aggregate metrics (macro-F1, AUROC, AUPRC). Adding synthetic data at this scale risks being a null result for the primary endpoint, which is fine scientifically but weakens the paper's main claim.

**Where synthetic data still has a defensible, testable role** is exactly what your proposal already identifies: the sparse-gene problem. Look at your per-gene P/LP counts. After removing TTN, DMD, AKAP9, and the quarantine genes, you will have genes like:

- PLN: maybe 15–25 P/LP
- ACTC1: maybe 20–35 P/LP
- MYL2: maybe 25–40 P/LP
- MYL3: maybe 15–30 P/LP
- CALM1/CALM2/CALM3: each fewer than 30 P/LP
- TRDN: maybe 20–40 P/LP
- SMAD3/TGFBR1/TGFBR2: each fewer than 50 P/LP

For these genes, a model trained only on real data will have wide confidence intervals and unreliable per-gene performance. Synthetic augmentation targeted specifically at sparse genes — not at the whole dataset — is a much sharper, more interesting scientific claim: "we show that LLM-generated synthetic evidence augmentation selectively improves classification in data-sparse genes without degrading performance in data-rich genes." That is publishable, testable, and not something any prior study has done.

**Revised recommendation for the synthetic arm:** rather than augmenting the full 40,000–50,000-variant training set, define the synthetic augmentation arm as targeted enrichment for genes with fewer than 50 real P/LP variants. Generate synthetic evidence-layer cases for those specific genes, test whether per-gene metrics improve in the sparse group specifically, and report the full-dataset metrics as unchanged (or slightly improved) as a safety check showing the augmentation doesn't hurt the well-represented genes. This is a more defensible, more interesting claim than "augmentation improves overall performance" on a dataset where you already have 35x more data than the prior benchmark.

---

## 2. Would this model be way better than CardioBoost?

Almost certainly yes on the metrics that matter, for reasons that are now documented in your own data. But let me be precise about what "better" means and where the gains will come from.

**Where you will clearly outperform CardioBoost:**

P/LP sensitivity should be higher because the primary model will be a binary
P/LP-vs-benign classifier trained on a much larger, better annotated dataset
than CardioBoost. Their reported TPR of 69.5% (CM) and 83.3% (IAS) at the 90%
confidence threshold was achieved with only 238 and 168 pathogenic training
variants respectively. With more than 9,000 P/LP rows before final split
decisions, plus modern annotation and source weighting, the model has
dramatically more signal to learn from. The VUS-facing improvement should come
from applying this calibrated binary model to VUS rows, not from treating VUS as
a primary biological class.

Gene coverage is broader — you cover ARVC, CPVT, and aortopathy genes that CardioBoost explicitly excluded. On those gene families, CardioBoost has no published performance at all, so you are not "beating" them — you are doing something they never attempted.

Feature richness is substantially greater. AlphaMissense alone, which postdates CardioBoost, is consistently shown to outperform REVEL and M-CAP (CardioBoost's best-performing in silico features) on most published benchmarks. Adding AlphaMissense as a feature should give measurable PR-AUC improvement independent of sample size.

**Where the comparison is more nuanced:**

CardioBoost's pathogenic-CM class came from 9,007 expert-diagnosed clinical patients. Those variants are high-trust, clinically validated pathogenic calls. Your ClinVar P/LP class at ≥1-star has inherent label noise — some single-submitter assertions will be wrong. Your ≥2-star subset (5,179 P/LP in the high-confidence slice) is more comparable to CardioBoost's quality, but that is only about 3x their CM training P/LP count, not 35x. This means your headline performance gains will be real but the label-quality caveat must be stated clearly.

**The honest performance expectation:** on your own held-out test set, you will almost certainly show PR-AUC ≥ 0.92–0.95 across disease groups, better calibration (lower Brier score), and dramatically lower indeterminate-classification rates than CardioBoost. On CardioBoost's own held-out variants (if used as an independent test stratum), you should match or exceed their numbers on arrhythmia genes and likely match on cardiomyopathy genes.

---

## 3. Is CardioBoost the best/most recent model? Full literature review of cardiogenetic variant classifiers

No — CardioBoost is not the most recent, but it is the most widely cited disease-specific cardiac classifier and the one most directly comparable to your design. Here is the complete landscape:

### Published disease-specific cardiogenetic variant classifiers

**APPRAISE / Ruklisa et al. (Genome Medicine, 2015)**
Bayesian syndrome- and gene-specific probability model. MYH7-focused initially, then extended. N: small, gene-specific. No VUS class. Predates gnomAD. Primarily a statistical model rather than ML. Still cited as the conceptual precursor to CardioBoost.

**PolyPhen-HCM (Jordan et al., AJHG, 2011)**
HCM-specific adaptation of PolyPhen, MYH7 and MYBPC3 only. Tiny training set. No VUS class. Historically important but superseded.

**PathoPredictor (Evans et al., Genome Research, 2019)**
Disease-specific ML classifier for cardiac genes. Random forest. ~600 training variants. Cardiomyopathy genes only. No published VUS class performance. Less widely cited than CardioBoost, partially because the training data and code are not as cleanly public.

**CardioBoost (Zhang et al., Genetics in Medicine, 2021)**
Current field benchmark. AdaBoost. 1,147 variants total. 22 genes. Binary B/P only (no VUS class in training). PR-AUC 0.91 (CM) and 0.96 (IAS). Fully public GitHub + Zenodo. This is the model you are directly benchmarking against.

**CardioClassifier (Tayoun et al., Genetics in Medicine, 2019)**
Not a pure ML classifier — it is a rule-based ACMG/AMP implementation for cardiac genes, not a trained model. It automates criterion application rather than learning from training data. Not directly comparable to your model but worth mentioning as an alternative approach in your Introduction.

**DYNA (Cedars-Sinai / Zenodo 2024)**
Used CardioBoost's public dataset plus ClinVar-extracted VUS panels, trained a protein-language-model-based classifier. Published as a preprint/conference paper. Uses the Zenodo record you already have. Binary B/P on the CardioBoost gene panel, with VUS scoring as a secondary output. This is the most recent direct competitor and the one you must benchmark against alongside CardioBoost.

**ClinPred, BayesDel, VEST4 (genome-wide, not cardiac-specific)**
These are genome-wide tools that happen to be applied to cardiac variants in some studies. Not disease-specific classifiers. CardioBoost already outperforms them; you should replicate that comparison.

**REVEL, M-CAP, PrimateAI (genome-wide in silico tools)**
These are the genome-wide benchmarks CardioBoost used as comparators. You should include them as comparators too, since any reviewer will ask "how do you compare to the best genome-wide tools."

**AlphaMissense as a classifier (Cheng et al., Science, 2023)**
AlphaMissense is not a disease-specific classifier — it's a genome-wide pathogenicity predictor. But it outperforms REVEL and M-CAP on most benchmarks and is now used as a comparator in any variant pathogenicity paper. You should include it as a standalone comparator (not just as a feature) to show your trained model adds value on top of AlphaMissense alone.

**The gap your study fills:** no published study has combined (1) a large
weighted cardiogenetics-specific binary pathogenic-vs-benign classifier with
explicit VUS scoring/prioritization, (2) the full spectrum of cardiac disease
groups including ARVC and aortopathy, (3) a training set exceeding 40,000 clean
binary cardiogenetic variants before final split decisions, (4)
AlphaFold/AlphaMissense/ESM1b/FoldX-aware features, and (5) external validation
on multiple independent cohorts. DYNA comes closest but is protein-language-model
based (not tabular CatBoost), uses CardioBoost-derived data, and does not build
from the same large-scale ClinVar/eMERGE/HiRO/CardioBoost integrated feature
matrix.

---

## 4. Is this paper publishable, and where?

**Yes, clearly publishable.** The real-data-scale-up arm alone — without any synthetic augmentation — is a complete, novel, methodologically strong paper. The reasons:

The sample size advance is the single largest in the field's history. Going from 1,147 (CardioBoost) to 40,000–50,000 training variants with full public reproducibility is a major methodological contribution that any reviewer in this space will recognize immediately.

The feature set advance is real and citable. AlphaFold/AlphaMissense/ESM1b/FoldX ΔΔG are all post-2021 tools that no prior cardiogenetic classifier has used. The ΔΔG structural stability angle in particular is novel for this task.

The gene panel advance covers disease groups no prior classifier has handled in a systematic, multi-gene, ML framework (ARVC, CPVT, aortopathy).

The VUS-facing analysis is a direct clinical improvement over CardioBoost's
binary-only output. The primary model should remain binary P/LP vs. B/LB, but
it should be calibrated and applied to VUS rows to estimate whether each VUS
looks more pathogenic-like or benign-like. This is closer to the clinical
question than treating VUS as a stable third biological class.

Full public reproducibility is increasingly required by top journals and is a genuine differentiator from all prior work in this space.

**Target journals, ranked:**

1. **Genetics in Medicine** — natural home given CardioBoost published there, direct lineage, same audience. Impact factor ~11. This should be your primary target.
2. **npj Genomic Medicine** — broad genomics audience, open access, good for methods-heavy papers with a clinical application angle.
3. **Circulation: Genomic and Precision Medicine** — American Heart Association journal, cardiac-audience reach, appropriate for a study with direct clinical cardiogenetics application.
4. **JAMIA** — if framing leans more toward the informatics/ML-methods side.
5. **European Heart Journal — Digital Health** — good fit if you want a European cardiology audience.

---

## 5. Things to pay close attention to

**The VUS class is your biggest methodological challenge.** ClinVar VUS labels are extremely heterogeneous — some are VUS because of genuinely uncertain evidence, some are VUS because they were classified years ago with incomplete data, some are VUS in one submitter and P/LP in another. Your model will learn "VUS" as a catch-all label for uncertainty, which is not the same as learning the biological properties that make a variant uncertain. You need to be very explicit in your Methods about what VUS means in your training set, and your Discussion should address whether the model is learning "variant features associated with classificatory uncertainty" rather than "variant features associated with biological pathogenic mechanism." These are related but not identical.

**Your dbNSFP coverage gap is a real issue.** You have 42,338 dbNSFP matches out of 86,889 total registry rows — about 49% coverage. That means roughly half your variants have missing values for REVEL, SIFT, PolyPhen, and other core in silico features. CatBoost handles missing values natively, which is one reason you chose it, but you need to report feature missingness explicitly and check whether model performance differs between variants with full dbNSFP coverage versus those without. If the model performs well only on the 49% with full in silico annotation, that is a significant limitation.

**The protein structure coverage is currently very low.** Your current matrix shows only 3,089 protein feature "ok" rows out of 86,889 — about 3.6% coverage. This is expected given that protein-position parsing requires complete HGVS-to-protein-residue mapping, which is hard to do reliably at scale from ClinVar's heterogeneous notation. But it means AlphaFold/DSSP/FreeSASA/FoldX features will be missing for ~96% of your training variants. These features will likely have little influence on the model's learned weights because of the sparse coverage, even though they are methodologically interesting. You should either invest in improving protein-position mapping coverage substantially before the first model, or frame the structural features as a "proof of concept for future integration" rather than claiming them as a primary contribution of this model. FoldX on 244 variants (your current batch) is far too sparse to affect a model trained on 40,000+ variants.

**Label leakage between ClinVar features and ClinVar labels.** Since ClinVar is both your label source and the source of some features (review stars, submitter count, conflict metadata), you must be careful. Review stars should not be a training feature for a model where the label came from ClinVar — a 2-star variant has more submitters agreeing on its label, so review stars would partially encode label reliability rather than variant biology. Exclude review stars from the feature set or only use them as a stratum variable for sensitivity analysis.

**The FoldX mismatch rate of 42/286 (14.7%) is worth investigating.** Reference amino acid mismatches between the AlphaFold structure and the ClinVar variant's reference residue usually indicate isoform differences or transcript annotation mismatches. These are not random — they cluster in specific genes or transcript annotation situations. Before merging FoldX ΔΔG values, you need to understand whether the 14.7% mismatch variants are a biased subset (e.g., concentrated in specific genes or consequence types) that would introduce systematic error if silently excluded.

**The CardioBoost benchmark comparison requires careful framing.** If you train on any CardioBoost variants (even from the Zenodo/DYNA repackaged files), you cannot claim your held-out test performance on CardioBoost's original holdout set is an independent benchmark. Your current proposal correctly handles this by defining two separate analyses, but make sure the primary published comparison uses the clean split: trained on non-CardioBoost data, tested on CardioBoost variants as an independent external validation stratum.

**TTN is sitting in a separate holdout stratum with 34,804 rows and 5,765 P/LP.** Do not leave this decision deferred too long. If you eventually include TTN truncating variants (frameshift/nonsense/splice), those 5,765 P/LP rows are a major addition — but the model then needs to learn TTN-specific mechanism (haploinsufficiency via truncation, not missense disruption), which is biologically quite different from sarcomere missense pathogenicity. A clean approach: train the primary model without TTN entirely, then train a sensitivity model that includes TTN truncating-only variants, and report whether TTN inclusion improves or degrades performance on the held-out non-TTN test set. This directly answers whether TTN belongs in a unified cardiac variant classifier or needs its own model.

---

## 2026-07-02 Project Status Update: Data Engineering and Model Results

This section records the current state of the project after dataset integration,
feature annotation, eMERGE coordinate rescue, binary modeling, and exploratory
3-class modeling.

### Current Corrected Modeling Matrix

The current primary matrix is:

```text
datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv
```

Key counts:

| Item | Count |
|---|---:|
| Total corrected variant-level rows | 85,677 |
| Primary binary rows | 42,990 |
| Benign binary rows | 33,859 |
| Pathogenic binary rows | 9,131 |
| VUS scoring rows | 41,731 |
| VEP annotated rows | 85,660 |
| SpliceAI available rows | 77,141 |

The true 3-class matrix is:

```text
datasets/modeling/ready/modeling_table_three_class_splits_weights_hiro_emerge_rescued.tsv
```

3-class trainable rows:

| Class | Rows |
|---|---:|
| Benign | 33,859 |
| VUS | 41,731 |
| Pathogenic | 9,131 |
| Total | 84,721 |

### Major Data Engineering Completed

| Area | Status |
|---|---|
| ClinVar | Clean cardiogenetics slice created; conflicts/quarantine genes handled |
| HiRO | 482 source records preserved; 409 resolved/linkable source records; 73 still unresolved |
| eMERGE | hg19/GRCh37 to GRCh38 coordinate rescue completed and validated |
| CardioBoost | Public dataset integrated and used as held-out external source |
| dbNSFP | Downloaded, validated, joined, selected feature table created |
| gnomAD | dbNSFP-derived + browser/Hail-derived fields merged |
| AlphaMissense | Direct join completed |
| VEP/HGVS/SpliceAI | Local VEP cache installed/tested; annotation completed |
| ClinGen | Gene/variant features merged into modeling matrix |
| Protein features | UniProt/AlphaFold/DSSP/FreeSASA features generated where mappings exist |
| FoldX | Full primary batch merged where eligible; mismatches flagged |

### eMERGE Coordinate Rescue

The initial poor eMERGE result was diagnosed as a coordinate/build issue rather
than true model failure. eMERGE source coordinates were lifted from hg19/GRCh37
to GRCh38 and validated against the local GRCh38 FASTA.

| QC item | Result |
|---|---:|
| eMERGE source rows processed | 2,754 |
| Lifted single-mapping rows | 2,754 |
| GRCh38 REF allele matches | 2,754 |
| Coordinate QC pass | 2,754 |
| VEP returned after rescue | 2,754 |

After rescue, eMERGE external performance improved from invalid/poor to strong.

### Data Splitting Strategy

Three evaluation strategies are now implemented:

| Split | Purpose |
|---|---|
| Internal grouped split | Baseline model performance and debugging |
| Source-held-out split | Main external validation on HiRO, eMERGE, and CardioBoost |
| Gene-stress split | Sparse-gene vs well-represented gene generalization |

The source-held-out split is the most important paper-facing split because it
tests independent external sources rather than random ClinVar-like mixing.

### Primary Binary Model Results

The primary model remains the binary CatBoost classifier:

```text
Pathogenic / Likely Pathogenic vs Benign / Likely Benign
```

This is the closest model to CardioBoost, because CardioBoost was also a binary
pathogenicity model with high-confidence thresholds rather than a true VUS-trained
3-class classifier.

Strict threshold `0.5`:

| Evaluation | Rows | AUROC | AUPRC | Sensitivity | Specificity | PPV |
|---|---:|---:|---:|---:|---:|---:|
| Internal grouped test | 6,449 | 0.9990 | 0.9971 | 0.9883 | 0.9968 | 0.9883 |
| External HiRO | 69 | 0.9877 | 0.9830 | 0.9630 | 0.8810 | 0.8387 |
| External eMERGE | 176 | 0.9910 | 0.9914 | 0.9792 | 0.9375 | 0.9495 |
| External CardioBoost | 185 | 0.9520 | 0.9812 | 0.9545 | 0.8113 | 0.9265 |
| Sparse-gene stress test | 5,890 | 0.9979 | 0.9574 | 0.9848 | 0.9947 | 0.8966 |

CardioBoost-comparable `0.1/0.9` deferral thresholds:

| External set | Rows | Sensitivity | Specificity | PPV | NPV | MCC | Deferral rate |
|---|---:|---:|---:|---:|---:|---:|---:|
| External combined | 430 | 0.937 | 0.943 | 0.960 | 0.977 | 0.876 | 0.119 |
| HiRO | 69 | 0.926 | 0.929 | 0.893 | 1.000 | 0.849 | 0.087 |
| eMERGE | 176 | 0.948 | 0.988 | 0.989 | 1.000 | 0.933 | 0.114 |
| CardioBoost | 185 | 0.932 | 0.887 | 0.953 | 0.903 | 0.806 | 0.135 |

Interpretation: the binary model with CardioBoost-style deferral is the strongest
and most defensible main model. It has high sensitivity, high specificity, high
PPV/NPV, and modest deferral rates.

### Full HiRO Source-Record Binary + Deferral Result

The HiRO source-record table contains all 482 records:

| HiRO source-record status | Count |
|---|---:|
| Total source records | 482 |
| With prediction | 408 |
| No prediction | 74 |

Binary model with `0.1/0.9` deferral, predicted rows only:

| True \ Pred | Benign | VUS/Deferred | Pathogenic |
|---|---:|---:|---:|
| Benign | 149 | 23 | 18 |
| VUS | 60 | 74 | 38 |
| Pathogenic | 0 | 5 | 41 |

Key signal: `41/46` predicted HiRO pathogenic source records were called
Pathogenic, `5/46` were deferred to VUS, and `0/46` were hard-called Benign.

### Exploratory True 3-Class Model Results

A separate true 3-class CatBoost branch was trained:

```text
Benign / VUS / Pathogenic
```

This model is not directly CardioBoost-comparable. It is a secondary VUS modeling
arm.

Internal 3-class test:

| Metric | Value |
|---|---:|
| Rows | 12,709 |
| Accuracy | 0.9294 |
| Macro-F1 | 0.9044 |
| Weighted kappa | 0.9152 |
| MCC | 0.8808 |
| OVR macro AUROC | 0.9840 |

Source-held-out HiRO 3-class:

| Metric | Value |
|---|---:|
| Rows | 179 |
| Accuracy | 0.8492 |
| Macro-F1 | 0.8188 |
| Weighted kappa | 0.7879 |
| MCC | 0.7160 |
| OVR macro AUROC | 0.9426 |

HiRO variant-level 3-class confusion matrix:

| True \ Pred | Benign | VUS | Pathogenic |
|---|---:|---:|---:|
| Benign | 34 | 8 | 0 |
| VUS | 8 | 100 | 2 |
| Pathogenic | 0 | 9 | 18 |

Full HiRO source-record true 3-class result, predicted rows only:

| True \ Pred | Benign | VUS | Pathogenic |
|---|---:|---:|---:|
| Benign | 113 | 77 | 0 |
| VUS | 16 | 153 | 3 |
| Pathogenic | 0 | 16 | 30 |

The true 3-class model improves exact 3-class accuracy and VUS recognition, but
it is more conservative for pathogenic calls than the binary-deferral model.

### Binary + Deferral vs True 3-Class On HiRO

Same 408 predicted HiRO source-record rows:

| Model | Exact 3-class accuracy | P/LP sensitivity | P/LP specificity | P/LP PPV | P/LP NPV |
|---|---:|---:|---:|---:|---:|
| Binary + 0.1/0.9 deferral | 0.647 | 0.891 | 0.887 | 0.719 | 1.000 |
| True 3-class model | 0.725 | 0.652 | 0.992 | 0.909 | 0.957 |

Interpretation:

- Binary + deferral is better for catching P/LP variants.
- True 3-class is better for exact Benign/VUS/Pathogenic matching and P/LP specificity.
- The binary + deferral model is more clinically balanced and more comparable to CardioBoost.
- The true 3-class model is valuable as a secondary exploratory VUS triage model.

### Current Paper Framing

The recommended paper framing is:

1. Primary model: binary cardiogenetics-specific P/LP vs B/LB CatBoost classifier.
2. Clinical interpretation: CardioBoost-style `0.1/0.9` high-confidence thresholds with VUS/indeterminate deferral.
3. External validation: HiRO, eMERGE, CardioBoost source-held-out testing.
4. Secondary analysis: true 3-class Benign/VUS/Pathogenic model for VUS triage.
5. Stress testing: sparse-gene performance and HiRO source-record/patient-event-level analysis.

### Current Result Files

Binary results:

```text
results/model_performance/binary_Results.MD
results/model_performance/source_heldout_hiro_emerge_rescued/PERFORMANCE_REVIEW.md
```

Three-class results:

```text
results/model_performance/three_tier_results.MD
results/model_performance/three_class_catboost_v0_review/THREE_CLASS_MODEL_REVIEW.md
```

Key model artifacts:

```text
results/models/primary_binary_catboost_v0_strict_source_heldout_hiro_emerge_rescued/
results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued/
results/model_performance/hiro_source_record_3class/
results/model_performance/hiro_source_record_true_three_class/
```

### Historical Remaining Work List, Now Largely Completed

| Task | Why |
|---|---|
| Feature leakage audit | Completed in `results/model_performance/leakage_ablation_audit/` |
| Calibration curves | Completed in `results/manuscript_figures_tables/` |
| Direct comparator benchmarks | Completed for official CardioBoost, REVEL, CADD, AlphaMissense where available |
| Missense-only model | Still optional; closest apples-to-apples analysis if we want a dedicated CardioBoost-scope model |
| HiRO source-record/patient-event analysis | Partly completed; `408 / 482` source records scored and `74` unresolved records audited |
| VUS reclassification or expert-review analysis | Still future work; best way to validate VUS-facing prioritization |
| Sparse-gene improvement | Still future work; current sparse-gene results are descriptive because external N is small |

## 2026-07-02 Continued Progress: CardioBoost Comparison, 3-Class Tuning, and Interpretation

This continuation records the conclusions after reviewing the CardioBoost paper
carefully and comparing it against the current local results. The goal is to
keep the project scientifically defensible and avoid overstating results.

### What CardioBoost Actually Did

CardioBoost is not a true 3-class Benign/VUS/Pathogenic classifier. It is a
binary pathogenicity probability model for rare missense variants in inherited
cardiac disease genes. Its clinical-style output uses two high-confidence
thresholds:

| Probability | CardioBoost interpretation |
|---:|---|
| `Pr >= 0.9` | Disease-causing |
| `Pr <= 0.1` | Benign / likely benign |
| `0.1 < Pr < 0.9` | Indeterminate / VUS-like |

Therefore, the correct direct comparison is not our true 3-class model. The
closest comparison is the current primary binary CatBoost model with the same
`0.1/0.9` high-confidence thresholds.

Important CardioBoost benchmark numbers from the paper:

| Metric | CardioBoost cardiomyopathy | CardioBoost arrhythmia |
|---|---:|---:|
| PR-AUC | `0.91` | `0.96` |
| Overall accuracy with `0.1/0.9` thresholds | `63.3%` | `81.2%` |
| High-confidence classified | `70.2%` | `88.3%` |
| High-confidence accuracy | `90.2%` | `91.9%` |
| Indeterminate rate | `29.8%` | `11.7%` |
| TPR / sensitivity | `69.5%` | `83.3%` |
| TNR / specificity | `56.0%` | `78.6%` |
| PPV | `86.3%` | `90.9%` |
| NPV | `96.6%` | `93.2%` |

CardioBoost also has several important design constraints:

- rare missense variants only;
- disease-specific cardiomyopathy and arrhythmia models;
- two binary classifiers rather than a true VUS model;
- training data in the hundreds to low thousands, not tens of thousands;
- high-quality disease-specific pathogenic sources and healthy-volunteer benign
  controls;
- clinical association/outcome validation in SHaRe, which we have not yet
  replicated.

### Our CardioBoost-Style Binary Result

The current primary binary CatBoost model is:

```text
P/LP vs B/LB
```

with CardioBoost-style thresholds:

```text
P(pathogenic) >= 0.9  -> P/LP-like
P(pathogenic) <= 0.1  -> B/LB-like
0.1-0.9               -> VUS / indeterminate / deferred
```

External combined binary result:

| Metric | Current model |
|---|---:|
| AUROC | `0.9768` |
| Sensitivity | `0.937` |
| Specificity | `0.943` |
| PPV | `0.960` |
| NPV | `0.977` |
| Deferral / indeterminate rate | `0.119` |
| High-confidence classified | `0.881` |

External source-specific binary + deferral results:

| External set | Rows | Sensitivity | Specificity | PPV | NPV | Deferral rate |
|---|---:|---:|---:|---:|---:|---:|
| HiRO | `69` | `0.926` | `0.929` | `0.893` | `1.000` | `0.087` |
| eMERGE | `176` | `0.948` | `0.988` | `0.989` | `1.000` | `0.114` |
| CardioBoost | `185` | `0.932` | `0.887` | `0.953` | `0.903` | `0.135` |

This was the strongest comparison before the official CardioBoost matched
benchmark was completed. It still matters because it shows the primary model's
source-held-out behavior across HiRO, eMERGE, and CardioBoost. However, the
cleaner direct comparison is now the completed same-row CardioBoost-matched
benchmark documented below.

The defensible wording is:

> Using a CardioBoost-style binary pathogenicity framework with `0.1/0.9`
> high-confidence thresholds, the larger cardiogenetics CatBoost model achieved
> high external sensitivity, specificity, PPV, NPV, and a low deferral rate,
> comparing favorably with published CardioBoost benchmark ranges.

The wording to avoid for this broad source-held-out table alone is:

> We definitively outperform CardioBoost.

That stronger claim required a CardioBoost-matched benchmark, which has now
been completed in the later section. The design requirement was:

```text
rare missense only
CardioBoost genes only
binary labels only
CardioBoost held out from training
report PR-AUC, ROC-AUC, sensitivity, specificity, PPV, NPV,
deferral rate, high-confidence accuracy, and high-confidence call rate
```

### Why Binary Is The Main Comparison And 3-Class Is Secondary

The binary model is the correct comparison to CardioBoost because CardioBoost
does not learn VUS as a class. Its VUS-like category is an indeterminate
probability zone.

The true 3-class model answers a different question:

```text
Can a model learn Benign vs VUS vs Pathogenic labels directly?
```

This is scientifically interesting but harder to defend as the primary endpoint
because VUS is an evidence state, not a stable biological state. A VUS can later
become benign or pathogenic as evidence accumulates.

Current model roles:

| Role | Model |
|---|---|
| Primary CardioBoost-comparable model | Binary CatBoost with `0.1/0.9` deferral |
| Primary true 3-class model | Original 3-class CatBoost argmax |
| Secondary sensitivity-oriented 3-class model | `path_3_5` weighted 3-class CatBoost |
| VUS-facing clinical analysis | Apply binary model to VUS rows and rank as benign-like, unresolved, or pathogenic-like |

### True 3-Class Model Status

The original true 3-class model directly predicts:

```text
Benign / VUS / Pathogenic
```

Current source-record external performance:

| External set | Rows with prediction | Accuracy | Macro-F1 | P/LP sensitivity | P/LP specificity | P/LP PPV |
|---|---:|---:|---:|---:|---:|---:|
| HiRO | `325` | `0.7662` | `0.7770` | `0.7143` | `0.9931` | `0.9259` |
| eMERGE | `1,974` | `0.8896` | `0.6565` | `0.5816` | `0.9760` | `0.5588` |
| CardioBoost | `195` | `0.5949` | `0.4735` | `0.6232` | `0.9649` | `0.9773` |

Interpretation:

- the 3-class model is more specific and cleaner for P/LP calls;
- it recognizes VUS better than binary-deferral does;
- it misses more P/LP by sending them to VUS;
- it should not replace the binary-deferral model as the primary clinical
  triage model.

### Three-Class Tuning Experiments

A methodological concern came up: it would be scientifically wrong to tune
thresholds or weights because they improve HiRO, eMERGE, or CardioBoost after
looking at those results. That would invite a reviewer to say the model was
fit to the external test set.

The corrected tuning rule is:

1. tune thresholds or weights using the internal validation split only;
2. freeze the selected rule/model;
3. evaluate HiRO, eMERGE, CardioBoost, and ClinVar only after freezing.

#### Threshold Tuning

For the true 3-class model, thresholds were selected using validation only. The
rule tested was:

```text
Call Pathogenic if prob_Pathogenic >= threshold.
Otherwise choose Benign vs VUS by whichever probability is larger.
```

Validation-selected thresholds:

| Validation selection rule | Selected threshold | Validation accuracy | Validation macro-F1 | Validation P/LP sensitivity | Validation P/LP specificity | Validation PPV |
|---|---:|---:|---:|---:|---:|---:|
| Max validation macro-F1 | `0.700` | `0.9401` | `0.9181` | `0.8325` | `0.9850` | `0.8704` |
| Max validation weighted kappa | `0.675` | `0.9398` | `0.9181` | `0.8415` | `0.9837` | `0.8615` |
| Max validation balanced accuracy | `0.400` | `0.9278` | `0.9019` | `0.9083` | `0.9620` | `0.7426` |
| Constrained sensitivity/specificity/PPV rule | `0.425` | `0.9291` | `0.9035` | `0.9016` | `0.9643` | `0.7528` |

Conclusion: validation-selected thresholding is scientifically valid, but it
does not clearly beat the original 3-class argmax model in a general way. More
sensitive thresholds recover more P/LP but reduce macro-F1 and/or PPV.

#### Weight Tuning

Class-weight tuning was retried sequentially, one model at a time, not in
parallel. The completed candidates were:

| Candidate | Class weights |
|---|---|
| `path_3_5` | Benign `1.0`, VUS `0.9`, Pathogenic `3.5` |
| `path_4_vus_0_8` | Benign `1.0`, VUS `0.8`, Pathogenic `4.0` |
| `path_5_vus_0_7` | Benign `1.0`, VUS `0.7`, Pathogenic `5.0` |

Argmax results:

| Candidate | Dataset | Accuracy | Macro-F1 | P/LP sensitivity | P/LP specificity | P/LP PPV |
|---|---|---:|---:|---:|---:|---:|
| Original 3-class | Validation | `0.9338` | `0.9101` | `0.8843` | `0.9713` | `0.7905` |
| `path_3_5` | Validation | `0.9282` | `0.9024` | `0.9083` | `0.9625` | `0.7449` |
| `path_4_vus_0_8` | Validation | `0.9213` | `0.8933` | `0.9226` | `0.9538` | `0.7066` |
| `path_5_vus_0_7` | Validation | `0.9072` | `0.8749` | `0.9406` | `0.9368` | `0.6424` |

External source-record P/LP sensitivity:

| Candidate | HiRO | eMERGE | CardioBoost |
|---|---:|---:|---:|
| Original 3-class | `0.7143` | `0.5816` | `0.6232` |
| `path_3_5` | `0.7429` | `0.6224` | `0.6667` |
| `path_4_vus_0_8` | `0.7429` | `0.6224` | `0.6884` |
| `path_5_vus_0_7` | `0.7429` | `0.6837` | `0.7536` |

Interpretation:

- weight tuning can improve P/LP sensitivity;
- it is not a free improvement;
- as pathogenic weight increases, validation macro-F1 and PPV decline;
- the original model remains the best balanced 3-class model;
- `path_3_5` is the most defensible sensitivity-oriented 3-class candidate so
  far, because it improves validation P/LP sensitivity while keeping macro-F1
  reasonably high.

The conclusion is not "tuning failed." The conclusion is:

> We did not find a tuning adjustment that clearly improves the 3-class model
> globally. We found a sensitivity-oriented alternative that catches more P/LP
> variants but trades off overall multiclass performance and precision.

### What We Have Actually Improved Over CardioBoost

The project appears to improve on CardioBoost in several ways, but each claim
needs the right denominator and comparison. This section was originally written
before the same-row CardioBoost benchmark was completed; it is now retained as
interpretive context, with the completed benchmark documented immediately below.

| Area | Current interpretation |
|---|---|
| Binary pathogenicity performance | Strong; compares favorably to CardioBoost-style high-confidence metrics |
| Sample size | Much larger than CardioBoost; tens of thousands of rows vs hundreds/low thousands |
| Feature breadth | Broader: dbNSFP, gnomAD v4/browser, AlphaMissense, VEP/SpliceAI, ClinGen, protein features, FoldX |
| Gene breadth | Broader cardiogenetics panel, including genes outside CardioBoost's published scope |
| VUS handling | More explicit: VUS scoring and true 3-class exploratory modeling |
| HiRO validation | Adds private phenotype-linked source-record testing that CardioBoost did not have |
| Clinical outcome validation | Not yet done; CardioBoost has SHaRe outcome association, so we cannot claim superiority there |

The cleanest broad source-held-out claim:

> Compared with CardioBoost's binary high-confidence framework, the current
> model achieves high external sensitivity and specificity with a low
> indeterminate rate in a much larger, feature-richer cardiogenetics dataset.

The stronger claim now supported by the completed same-row benchmark and paired
bootstrap confidence intervals:

> Our source-held-out CatBoost model outperformed the official public
> CardioBoost model on the same CardioBoost-eligible matched external rows.

The language still needs care. This is a public-model matched benchmark, not a
full reconstruction of CardioBoost's private original development dataset.

### CardioBoost-Matched Subset: Planned Design Now Completed

The benchmark design was:

| Filter / rule | Reason |
|---|---|
| Rare missense variants only | CardioBoost is rare-missense only |
| CardioBoost genes only | Avoid comparing on genes CardioBoost never attempted |
| Binary B/LB vs P/LP labels only | CardioBoost is binary |
| Hold out CardioBoost-derived variants | Avoid leakage |
| Use `0.1/0.9` thresholds | Match CardioBoost's reported clinical threshold framework |
| Report PR-AUC and ROC-AUC | Match CardioBoost global discrimination metrics |
| Report high-confidence classified rate | Match CardioBoost's call-rate metric |
| Report high-confidence accuracy | Match CardioBoost's clinical-threshold accuracy |
| Report sensitivity, specificity, PPV, NPV | Match Table 1 style |
| Compare REVEL, CADD, AlphaMissense, and CardioBoost score where available | Reviewer-expected baseline |

This analysis has now been completed and is the main "are we better than
CardioBoost on the same public-model matched rows?" answer.

### Current Best Overall Scientific Framing

The strongest framing is:

1. **Primary endpoint:** binary P/LP-vs-B/LB prediction with
   CardioBoost-style `0.1/0.9` deferral.
2. **Primary comparison:** CardioBoost-style high-confidence metrics and PR-AUC
   on a matched rare-missense subset.
3. **External validation:** HiRO, eMERGE, and CardioBoost source-held-out
   strata.
4. **VUS analysis:** apply the binary model to VUS as a prioritization/ranking
   task rather than claiming VUS is a biological class.
5. **Secondary model:** true 3-class CatBoost model for direct
   Benign/VUS/Pathogenic label prediction.
6. **Sensitivity analysis:** `path_3_5` weighted 3-class model as a
   P/LP-sensitive operating point, not the primary 3-class model.
7. **Future clinical validation:** VUS expert re-review, HiRO phenotype-linked
   analysis, and any available outcome association.

### Updated Result Files

New/updated result files from this continuation:

```text
results/model_performance/primary_binary_three_zone_all_sources/PRIMARY_BINARY_THREE_ZONE_ALL_SOURCES.md
results/model_performance/primary_binary_three_zone_binary_only/PRIMARY_BINARY_BINARY_ONLY_REPORT.md
results/model_performance/true_three_class_all_sources/TRUE_THREE_CLASS_ALL_SOURCES_REPORT.md
results/model_performance/three_class_scientific_tuning/THREE_CLASS_TUNING_REVIEW.md
results/model_performance/three_class_scientific_tuning/VALIDATION_SELECTED_THRESHOLD_REPORT.md
results/model_performance/three_class_weight_tuning/three_class_weight_tuning_results.completed_candidates.tsv
```

Candidate sensitivity-tuned 3-class models:

```text
results/model_performance/three_class_weight_tuning/path_3_5.cbm
results/model_performance/three_class_weight_tuning/path_4_vus_0_8.cbm
results/model_performance/three_class_weight_tuning/path_5_vus_0_7.cbm
```

### Updated Immediate Next Steps From This Earlier Stage

| Priority | Task | Why |
|---:|---|---|
| 1 | Build CardioBoost-matched rare-missense binary benchmark | Completed below |
| 2 | Add comparator scores: CardioBoost, REVEL, CADD, AlphaMissense where available | Completed below |
| 3 | Generate PR-AUC/ROC-AUC and high-confidence metrics on matched subset | Completed below |
| 4 | Keep original 3-class model as balanced secondary model | Best overall 3-class performance |
| 5 | Report `path_3_5` only as sensitivity-oriented secondary 3-class model | Better P/LP sensitivity but tradeoff in macro-F1/PPV |
| 6 | Add CIs and leakage audit before final manuscript claim | Prevents overclaiming |

## CardioBoost-Matched Benchmark Completed

Updated: 2026-07-02

We completed the first real head-to-head comparison against the official
CardioBoost model.

Official CardioBoost model artifacts were downloaded locally from the public
repository into:

```text
datasets/cardioboost/public_dataset/official_cardioBoost_model_snapshot/
```

The official cardiomyopathy and arrhythmia AdaBoost models were run locally on
their all-rare mutation tables. This produced:

| Official output | Rows |
|---|---:|
| CardioBoost cardiomyopathy all-rare predictions | 65,475 |
| CardioBoost arrhythmia all-rare predictions | 42,411 |

The fair comparison is not the full 85k-row matrix, because CardioBoost is only
intended for rare missense variants in its supported genes. Therefore, we built
a matched benchmark using CardioBoost-gene, missense, binary-labeled variants
and compared:

- our primary binary CatBoost model;
- official CardioBoost;
- REVEL;
- CADD PHRED;
- AlphaMissense.

The exact-coordinate match was too small to use as a final claim:

| Exact-coordinate slice | Rows | Benign | Pathogenic | VUS |
|---|---:|---:|---:|---:|
| Matches to official all-rare universe | 64 | 34 | 8 | 22 |
| Binary only | 42 | 34 | 8 | 0 |

This likely reflects genome-build and transcript-version differences, so it is
treated as QC only.

The primary matched benchmark uses HGVS cDNA matching:

| Benchmark | Rows | Unique variants | CM rows | Arrhythmia rows |
|---|---:|---:|---:|---:|
| HGVS cDNA primary | 2,462 | 2,344 | 1,765 | 697 |
| HGVS cDNA plus protein fallback | 2,477 | 2,359 | 1,776 | 701 |

### What This Benchmark Represents

This benchmark is a same-row comparison between our CatBoost model and the
official CardioBoost model. The official CardioBoost cardiomyopathy and
arrhythmia AdaBoost models were run locally using the public CardioBoost model
objects and official public all-rare mutation tables. The resulting official
CardioBoost pathogenicity scores were then matched to our own CatBoost
predictions.

Therefore, for the reported matched rows:

- our CatBoost model and official CardioBoost were scored on the same variants;
- the same binary labels were used for metric calculation;
- the same CardioBoost-style `0.1/0.9` deferral thresholds were applied;
- REVEL, CADD PHRED, AlphaMissense, and M-CAP were also included as comparator
  scores when available.

The 250-row external benchmark is not a reconstruction of the complete original
CardioBoost paper dataset. CardioBoost used additional private/curated training
data that are not fully available in the public release. That is not a problem
for this comparison, because the goal here is not to reproduce every original
CardioBoost training variant. The goal is to ask whether both models perform
differently on the same CardioBoost-eligible variants that are present in our
current integrated matrix and can be matched to the official CardioBoost
scoring universe.

The 250-row external HGVS cDNA benchmark consists of source-held-out rows:

| External source stratum | Rows |
|---|---:|
| CardioBoost-source held out | 164 |
| eMERGE-source held out | 73 |
| HiRO-source held out | 13 |
| Total | 250 |

The stricter CardioBoost-source-only benchmark uses only the 164 held-out
CardioBoost-source rows. This is the cleanest direct CardioBoost-source
comparison. The 250-row version is broader because it also includes HiRO and
eMERGE external rows that are CardioBoost-eligible and match the official
CardioBoost scoring universe.

Exact coordinate matching recovered too few binary rows to support a manuscript
claim, likely because of build and transcript-version differences between the
official CardioBoost tables and our GRCh38-centered registry. For that reason,
exact-coordinate matching is treated as QC only, and gene + HGVS cDNA matching
is the primary benchmark. HGVS protein-position matching was used only as a
small sensitivity check.

Important fairness note: our model held out CardioBoost-source rows during
source-held-out training. The official CardioBoost model may have seen some
public CardioBoost-source variants during its own development, because those
data come from the CardioBoost release. If anything, this should favor the
official CardioBoost model on those rows, not our model. This makes the
CardioBoost-source-only matched result especially useful, but still not the
same as reproducing CardioBoost's private original training and validation
cohorts.

The correct manuscript wording is:

> We compared our source-held-out CatBoost model against the official public
> CardioBoost model on the same CardioBoost-eligible variants matched by gene
> and HGVS cDNA. The benchmark included 250 external rows, including 164
> CardioBoost-source held-out rows. Both models were scored on the same variants
> using the same binary labels and the same `0.1/0.9` deferral thresholds.

The wording to avoid is:

> We reproduced the full CardioBoost paper dataset.

We did not reproduce the full private CardioBoost dataset, and we should not
claim that. The fair claim is a matched, same-row, public-model benchmark.

For manuscript claims, the source-held-out external rows are the key rows, not
the full benchmark that includes training rows.

External HGVS cDNA-matched benchmark:

| Model | Rows | AUROC | AUPRC | Sensitivity | Specificity | PPV | NPV | Deferral | High-conf accuracy |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Our CatBoost | 250 | 0.980 | 0.990 | 0.959 | 0.543 | 0.964 | 1.000 | 0.152 | 0.972 |
| Official CardioBoost | 250 | 0.843 | 0.902 | 0.876 | 0.148 | 0.836 | 1.000 | 0.244 | 0.847 |
| REVEL local | 247 | 0.903 | 0.940 | 0.530 | 0.037 | 0.936 | 1.000 | 0.607 | 0.938 |
| CADD PHRED | 250 | 0.872 | 0.914 | 1.000 | 0.012 | 0.693 | 1.000 | 0.020 | 0.694 |
| AlphaMissense | 248 | 0.937 | 0.969 | 0.635 | 0.457 | 0.991 | 0.902 | 0.403 | 0.966 |

External CardioBoost-source-only benchmark:

| Model | Rows | AUROC | AUPRC | Sensitivity | Specificity | PPV | NPV | Deferral | High-conf accuracy |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Our CatBoost | 164 | 0.972 | 0.990 | 0.967 | 0.366 | 0.952 | 1.000 | 0.146 | 0.957 |
| Official CardioBoost | 164 | 0.788 | 0.905 | 0.870 | 0.073 | 0.843 | 1.000 | 0.207 | 0.846 |
| REVEL local | 161 | 0.880 | 0.944 | 0.483 | 0.049 | 0.935 | 1.000 | 0.602 | 0.938 |
| CADD PHRED | 164 | 0.810 | 0.903 | 1.000 | 0.024 | 0.759 | 1.000 | 0.006 | 0.761 |
| AlphaMissense | 162 | 0.906 | 0.966 | 0.645 | 0.244 | 0.987 | 0.833 | 0.438 | 0.967 |

Paired bootstrap confidence intervals were added on 2026-07-03 using 5,000
paired bootstrap replicates. The primary CI uses a row-stratified paired
bootstrap; a variant-cluster stratified bootstrap was also run as a sensitivity
analysis because some variants appear in more than one CardioBoost panel
context.

External all paired differences, row-stratified bootstrap:

| Metric | Our minus CardioBoost | 95% CI |
|---|---:|---:|
| AUROC | 0.138 | 0.083-0.195 |
| AUPRC | 0.087 | 0.048-0.132 |
| Sensitivity | 0.083 | 0.030-0.136 |
| Specificity | 0.395 | 0.272-0.519 |
| PPV | 0.128 | 0.082-0.172 |
| Deferral | -0.092 | -0.156 to -0.028 |
| High-confidence accuracy | 0.125 | 0.082-0.167 |

External CardioBoost-source-only paired differences, row-stratified bootstrap:

| Metric | Our minus CardioBoost | 95% CI |
|---|---:|---:|
| AUROC | 0.184 | 0.101-0.274 |
| AUPRC | 0.085 | 0.040-0.135 |
| Sensitivity | 0.098 | 0.041-0.163 |
| Specificity | 0.293 | 0.146-0.463 |
| PPV | 0.109 | 0.058-0.158 |
| Deferral | -0.061 | -0.134 to 0.012 |
| High-confidence accuracy | 0.111 | 0.061-0.159 |

The variant-cluster bootstrap gave the same conclusion for the key ranking
metrics. External all AUROC difference was 0.138 with 95% CI 0.081-0.201 and
AUPRC difference was 0.087 with 95% CI 0.045-0.135. CardioBoost-source-only
AUROC difference was 0.184 with 95% CI 0.100-0.273 and AUPRC difference was
0.085 with 95% CI 0.039-0.137.

Interpretation: the CI analysis supports a robust advantage for our model over
official CardioBoost on AUROC, AUPRC, sensitivity, PPV, and high-confidence
accuracy in both matched external benchmarks. The deferral rate is significantly
lower in the 250-row external-all benchmark. In the stricter 164-row
CardioBoost-source-only benchmark, deferral trends lower but the CI crosses
zero.

Interpretation:

- this is the strongest current evidence that our model compares favorably with
  CardioBoost on CardioBoost's own intended rare-missense binary task;
- the result should be described as a matched external benchmark, not a
  full-matrix comparison;
- the exact-coordinate result should not be used as a main claim;
- paired bootstrap confidence intervals are now available and should be included
  in final manuscript tables.

Main report:

```text
results/model_performance/cardioboost_matched_benchmark/CARDIOBOOST_MATCHED_BENCHMARK_REPORT.md
```

Figures:

```text
results/model_performance/cardioboost_matched_benchmark/figures/hgvs_cdot_primary_external_all_roc.png
results/model_performance/cardioboost_matched_benchmark/figures/hgvs_cdot_primary_external_all_precision_recall.png
results/model_performance/cardioboost_matched_benchmark/figures/hgvs_cdot_primary_external_cardioboost_roc.png
results/model_performance/cardioboost_matched_benchmark/figures/hgvs_cdot_primary_external_cardioboost_precision_recall.png
```

## Leakage and Ablation Audit, 2026-07-03

We completed a stricter leakage/ablation audit for the primary binary CatBoost
model. The audit retrained matched models using the same source-held-out split,
training slice, sample weights, and CatBoost settings, while removing
reviewer-sensitive feature families one at a time.

This was not a post-hoc deletion of columns at prediction time. Each ablation
was retrained from scratch.

Audit location:

```text
results/model_performance/leakage_ablation_audit/
```

Main report:

```text
results/model_performance/leakage_ablation_audit/LEAKAGE_ABLATION_AUDIT_REPORT.md
```

Feature-family ablations:

| Ablation | Features Used | Dropped From Baseline | Purpose |
|---|---:|---:|---|
| baseline_current | 338 | 0 | Current primary feature-selection rules |
| no_gene_identity | 328 | 10 | Remove direct gene/gene-level identity features |
| no_vep_consequence_hgvs | 302 | 36 | Remove VEP/gnomAD consequence, impact, transcript, and HGVS-like annotations while keeping SpliceAI scores |
| no_clinvar_metadata | 338 | 0 | Check whether residual ClinVar/review/confidence fields remain after the default leakage filter |
| no_source_record_aggregates | 300 | 38 | Remove HiRO/source-record aggregate and phenotype/evidence fields |
| no_hgvs_transcript_identifiers | 316 | 22 | Remove HGVS, transcript, protein accession, rsID/CAID, and identifier-like fields |
| strict_low_leakage | 246 | 92 | Remove gene identity, consequence/HGVS/transcript identifiers, ClinVar metadata, and source-record aggregates together |

External-all ablation results:

| Ablation | Rows | AUROC | AUPRC | Sensitivity @0.5 | Specificity @0.5 | PPV @0.5 | Deferral | High-conf accuracy |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline_current | 430 | 0.977 | 0.984 | 0.973 | 0.897 | 0.932 | 0.098 | 0.966 |
| no_gene_identity | 430 | 0.977 | 0.984 | 0.973 | 0.886 | 0.925 | 0.100 | 0.966 |
| no_vep_consequence_hgvs | 430 | 0.973 | 0.980 | 0.957 | 0.886 | 0.924 | 0.160 | 0.967 |
| no_clinvar_metadata | 430 | 0.977 | 0.984 | 0.973 | 0.897 | 0.932 | 0.098 | 0.966 |
| no_source_record_aggregates | 430 | 0.977 | 0.983 | 0.973 | 0.891 | 0.929 | 0.102 | 0.966 |
| no_hgvs_transcript_identifiers | 430 | 0.978 | 0.985 | 0.973 | 0.914 | 0.943 | 0.093 | 0.964 |
| strict_low_leakage | 430 | 0.974 | 0.981 | 0.973 | 0.811 | 0.883 | 0.212 | 0.973 |

Interpretation:

- removing gene identity did not collapse external performance, which argues
  against simple gene memorization;
- removing residual ClinVar metadata changed nothing because the default
  leakage filter was already excluding those columns;
- removing source-record aggregates barely changed external-all performance,
  so the main model is not dependent on HiRO phenotype/source-record shortcuts;
- removing VEP/consequence/HGVS features caused only a small AUROC/AUPRC drop,
  although deferral increased;
- the strict low-leakage stress test remained strong for AUROC and AUPRC, but
  specificity and PPV fell and deferral increased, especially on CardioBoost.

The audit supports the claim that the primary model is not driven by one obvious
leakage-prone feature family. For the manuscript, the current primary model can
remain the feature-rich model, while `strict_low_leakage` should be reported as
a robustness stress test, not as the deployable model.

Figure:

```text
results/model_performance/leakage_ablation_audit/external_all_ablation_metrics.png
```

## Manuscript Figures, Calibration, and Stratified Performance, 2026-07-03

We generated a manuscript-facing result bundle for the primary binary CatBoost
model with `0.1/0.9` deferral. The bundle uses the source-held-out model
predictions and the rescue-aware ready matrix.

Bundle location:

```text
results/manuscript_figures_tables/
```

Main report:

```text
results/manuscript_figures_tables/MANUSCRIPT_RESULTS_BUNDLE_REPORT.md
```

Generated figures:

```text
results/manuscript_figures_tables/figures/primary_binary_roc_pr_panel.png
results/manuscript_figures_tables/figures/primary_binary_calibration_reliability.png
results/manuscript_figures_tables/figures/source_flow_modeling_rows.png
results/manuscript_figures_tables/figures/feature_group_coverage_manuscript.png
results/manuscript_figures_tables/figures/external_all_stratified_performance_heatmap.png
```

Generated tables:

```text
results/manuscript_figures_tables/tables/external_validation_performance.tsv
results/manuscript_figures_tables/tables/calibration_brier_summary.tsv
results/manuscript_figures_tables/tables/calibration_curve_points.tsv
results/manuscript_figures_tables/tables/source_flow_modeling_rows.tsv
results/manuscript_figures_tables/tables/feature_group_coverage_manuscript.tsv
results/manuscript_figures_tables/tables/stratified_performance.tsv
```

External validation summary:

| Subset | Rows | P/LP | B/LB | AUROC | AUPRC | Brier | Sensitivity | Specificity | PPV | Deferral | High-conf accuracy |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Validation | 6,384 | 1,331 | 5,053 | 1.000 | 1.000 | 0.003 | 0.995 | 0.997 | 0.989 | 0.009 | 0.999 |
| External all | 430 | 255 | 175 | 0.977 | 0.985 | 0.053 | 0.965 | 0.886 | 0.925 | 0.119 | 0.966 |
| HiRO | 69 | 27 | 42 | 0.988 | 0.983 | 0.068 | 0.963 | 0.881 | 0.839 | 0.087 | 0.952 |
| eMERGE | 176 | 96 | 80 | 0.991 | 0.991 | 0.030 | 0.979 | 0.938 | 0.949 | 0.114 | 0.994 |
| CardioBoost | 185 | 132 | 53 | 0.952 | 0.981 | 0.069 | 0.955 | 0.811 | 0.926 | 0.135 | 0.944 |

All-label clinical triage summary:

This second table applies the same binary model to all Benign, VUS, and
Pathogenic rows using three zones:

```text
<=0.1      Benign-like
0.1-0.9    Deferred / VUS-zone
>=0.9      Pathogenic-like
```

This is not a true three-class classifier. It is a clinical triage analysis.
The exact 3-zone accuracy is therefore informative but not the primary
performance metric, because VUS is an evidence category rather than a stable
biological class.

| Dataset | Rows with prediction | Benign | VUS | Pathogenic | Pred benign-like | Pred deferred | Pred pathogenic-like | Exact 3-zone accuracy | P/LP sensitivity | Specificity vs non-P/LP | P/LP PPV | VUS deferral rate |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| All variant rows | 85,677 | 33,923 | 42,361 | 9,148 | 43,292 | 15,249 | 27,136 | 0.670 | 0.981 | 0.763 | 0.332 | 0.178 |
| ClinVar variant rows | 83,915 | 33,796 | 40,872 | 9,038 | 42,668 | 14,574 | 26,673 | 0.674 | 0.982 | 0.763 | 0.334 | 0.174 |
| HiRO variant rows | 240 | 43 | 139 | 31 | 95 | 78 | 67 | 0.568 | 0.903 | 0.813 | 0.452 | 0.325 |
| eMERGE variant rows | 2,754 | 130 | 2,394 | 111 | 719 | 1,176 | 859 | 0.498 | 0.946 | 0.715 | 0.127 | 0.427 |
| CardioBoost variant rows | 353 | 57 | 0 | 137 | 60 | 92 | 201 | 0.820 | 0.934 | 0.895 | 0.955 | 0.261 |
| HiRO source records | 408 / 482 | 190 | 172 | 46 | 209 | 102 | 97 | 0.647 | 0.891 | 0.845 | 0.423 | 0.250 |
| eMERGE source records | 2,754 | 137 | 2,490 | 127 | 719 | 1,176 | 859 | 0.497 | 0.945 | 0.719 | 0.140 | 0.427 |
| CardioBoost source records | 355 | 156 | 0 | 199 | 60 | 92 | 203 | 0.642 | 0.869 | 0.808 | 0.852 | 0.259 |

Interpretation: the all-label triage analysis is much harsher than the binary
endpoint because VUS rows are treated as if the "correct" output should be
deferred. Many ClinVar and eMERGE VUS rows are pushed into benign-like or
pathogenic-like zones. This should not be framed as model error in the same way
as P/LP-to-benign or B/LB-to-pathogenic errors. It is better framed as VUS
prioritization: which VUS look benign-like, which remain uncertain, and which
look pathogenic-like.

Clinically important hard-error counts:

| Dataset | P/LP-to-benign errors | Benign-to-P/LP escalations | VUS-to-P/LP prioritizations |
|---|---:|---:|---:|
| All variant rows | 22 | 43 | 18,020 |
| ClinVar variant rows | 19 | 36 | 17,666 |
| HiRO source records | 0 | 18 | 38 |
| eMERGE source records | 0 | 1 | 738 |
| CardioBoost source records | 5 | 30 | 0 |

The VUS-to-P/LP number is large because the model was trained to separate
pathogenic-like from benign-like evidence, not to learn VUS as a biological
class. This is useful for prioritization, but it must be validated with future
reclassification, expert review, or case-level evidence before making clinical
claims.

Calibration summary:

| Subset | Rows | Brier | Mean predicted P/LP probability | Observed P/LP fraction | Mean error |
|---|---:|---:|---:|---:|---:|
| Validation | 6,384 | 0.003 | 0.211 | 0.208 | 0.003 |
| External all | 430 | 0.053 | 0.627 | 0.593 | 0.034 |
| HiRO | 69 | 0.068 | 0.444 | 0.391 | 0.052 |
| eMERGE | 176 | 0.030 | 0.575 | 0.545 | 0.029 |
| CardioBoost | 185 | 0.069 | 0.745 | 0.714 | 0.031 |

Interpretation: the model is very well calibrated on the internal validation
split. It is mildly overconfident on external data, with the largest mean error
in HiRO. This is not surprising because HiRO is small and patient-linked. The
final manuscript should include calibration plots and Brier score, and it may
be worth testing Platt/isotonic calibration later if we want calibrated
probabilities instead of ranking/triage scores.

Feature coverage across the final ready matrix:

| Feature group | Covered rows | Coverage |
|---|---:|---:|
| VEP | 85,660 / 85,677 | 99.98% |
| SpliceAI | 77,141 / 85,677 | 90.04% |
| gnomAD final AF | 62,937 / 85,677 | 73.46% |
| gnomAD observed | 55,387 / 85,677 | 64.65% |
| dbNSFP | 43,774 / 85,677 | 51.09% |
| AlphaMissense direct | 38,023 / 85,677 | 44.38% |
| FoldX DDG | 313 / 85,677 | 0.37% |

Key external-all stratified results:

| Stratum | Rows | AUROC | AUPRC | Sensitivity | Specificity | PPV |
|---|---:|---:|---:|---:|---:|---:|
| Arrhythmia genes | 265 | 0.985 | 0.988 | 0.975 | 0.907 | 0.939 |
| Cardiomyopathy genes | 151 | 0.965 | 0.980 | 0.946 | 0.847 | 0.906 |
| Missense consequence | 315 | 0.987 | 0.989 | 0.977 | 0.887 | 0.914 |
| dbNSFP matched | 348 | 0.985 | 0.989 | 0.981 | 0.896 | 0.938 |
| dbNSFP missing | 82 | 0.949 | 0.962 | 0.878 | 0.854 | 0.857 |
| gnomAD observed | 326 | 0.984 | 0.984 | 0.977 | 0.901 | 0.919 |
| gnomAD confirmed absent | 32 | 0.787 | 0.839 | 0.692 | 0.737 | 0.643 |
| Sparse P/LP genes | 24 | 0.986 | 0.981 | 1.000 | 0.929 | 0.909 |

Interpretation: performance is strongest in dbNSFP-matched, missense, and
gnomAD-observed variants. dbNSFP-missing variants still perform reasonably but
are weaker. The gnomAD-confirmed-absent stratum is small and performs worse, so
it should be treated as an uncertainty/high-caution subgroup. Sparse-gene
external rows look good, but the sample size is only 24, so this should be
reported descriptively.

Updated next steps:

| Priority | Task | Why |
|---:|---|---|
| 1 | Add bootstrap CIs for external validation and key stratified results | Needed for final paper tables, especially small external strata |
| 2 | Test optional probability calibration on validation only | Could improve external probability calibration without changing ranking |
| 3 | Generate final multi-panel manuscript figure layout | Combine ROC/PR, calibration, feature coverage, and source flow |
| 4 | Improve unresolved HiRO coordinate mapping | More HiRO source records can be scored |
| 5 | Decide final manuscript framing for strict low-leakage and stratified subgroup results | Robustness/supporting analyses, not the main deployable model |

## HiRO Unresolved-Record Audit, 2026-07-03

We completed a row-level audit of the remaining HiRO source records that do not
link to a promoted variant-registry row.

Report:

```text
results/data_characteristics/hiro_unresolved_record_audit/HIRO_UNRESOLVED_74_AUDIT_REPORT.md
```

Main tables:

```text
results/data_characteristics/hiro_unresolved_record_audit/tables/hiro_unresolved_74_record_audit.tsv
results/data_characteristics/hiro_unresolved_record_audit/tables/hiro_unresolved_review_candidates.tsv
results/data_characteristics/hiro_unresolved_record_audit/tables/blocker_by_label.tsv
results/data_characteristics/hiro_unresolved_record_audit/tables/blocker_by_consequence.tsv
results/data_characteristics/hiro_unresolved_record_audit/tables/recommended_actions.tsv
```

Current HiRO source-record status:

| Category | Source records |
|---|---:|
| Total HiRO source records | 482 |
| Linked to registry/modelable variant ID | 408 |
| Unresolved no coordinates | 73 |
| Rescued but not promoted | 1 |
| Total unresolved/not promoted | 74 |

Unresolved/not-promoted labels:

| Label | Rows |
|---|---:|
| Benign | 30 |
| Pathogenic | 9 |
| VUS | 35 |

Main blockers:

| Blocker | Rows | Interpretation |
|---|---:|---|
| Needs HGVS normalization, missense | 26 | cDNA/protein text exists, but no safe GRCh38 VCF key yet |
| Needs HGVS normalization, indel/splice/LoF | 22 | frameshift/indel/splice/nonsense HGVS needs transcript-aware normalization |
| No variant text or HGVS | 15 | current processed table does not contain enough variant detail; need original report/raw source |
| Gene symbol/readthrough cleanup | 5 | e.g. FPGT/FPGT-TNNI3K/readthrough/alias issues before mapping |
| ClinVar review candidate not promoted | 3 | candidate exists but evidence was too weak for automatic promotion |
| Partial coordinate/allele | 1 | some coordinate fields present but incomplete |
| Manual review | 1 | does not fit simple categories |
| Rescued coordinate not promoted | 1 | prior rescue candidate did not enter promoted registry |

Review candidates found by the new audit:

| HiRO record | Label | Candidate | Why not auto-promoted |
|---|---|---|---|
| FHOD3 c.646G>A / p.V216I | VUS | `18-36611991-G-A` | unique internal matrix gene+cDNA candidate, but transcript context should be checked |
| FLNC c.3301G>A / p.G1101S | VUS | `7-128845054-G-A` | unique internal matrix gene+cDNA candidate, but transcript context should be checked |
| DSC2 c.2187G>A / p.Ala733Thr | VUS | `18-31070779-C-T` | protein-level/isoform-dependent match to an existing HiRO variant |
| RYR2 c.6022+5 / likely c.6022+5G>A | VUS | ClinVar `222789` | source text is malformed/truncated, manual review needed |
| CASQ2 p.Asp383dup | VUS | ClinVar `179599` | repeat/duplication representation differs; manual review needed |
| AKAP9 p.Arg1609Gln | VUS | ClinVar `190484` | AKAP9 is already a questionable/excluded gene; manual review only |

Important interpretation:

- The unresolved HiRO rows are not accidental duplicate drops. They remain
  preserved as source records but cannot yet be joined to a safe variant-level
  feature row.
- The largest blocker is transcript-aware HGVS-to-genomic normalization.
- Offline VEP cannot parse HGVS input, so the offline VEP cache is not enough
  for this rescue step.
- The next technical rescue path is online VEP HGVS, VariantValidator,
  Mutalyzer, or a local transcript-aware HGVS normalizer.
- The two strongest immediate candidates are both VUS rows, so promoting them
  would improve HiRO VUS scoring coverage but would not change the current
  binary HiRO P/LP vs B/LB external performance.

Decision for now: do not change the registry/model matrix from this audit
alone. Use the generated review-candidate table for manual transcript review,
then promote only records with a unique and transcript-consistent GRCh38
chrom-pos-ref-alt.
