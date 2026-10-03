# New Data Pipeline

This project layer is the clean replacement for the original Colab monolith.
It starts by normalizing the April 2026 CASPER WES and VERDICT files in
`new_data/`.

For a short handoff summary of the current ML results, preferred triage model,
and LLM run order, see `docs/current_status_and_next_steps.md`.

## Inputs

- `ACMG AI Project_Classified Variants_07April2026(CASPER WES).csv`
- `ACMG AI Project_Classified Variants_07April2026(VERDICT).csv`
- `ACMG_CASPER_Wes_Patient_Data_Apr_2026(in).csv`
- `ACMG_Verdict_Patient_Data_Apr_2026(in).csv`

## Run Preprocessing

```bash
python3 scripts/build_new_cohort.py
python3 scripts/validate_new_cohort.py
```

Outputs are written to `outputs/new_data/`:

- `normalized_variants.csv`: one row per classified variant.
- `normalized_patients.csv`: one row per patient with derived phenotype features.
- `analysis_cohort.csv`: variant rows joined to patient features.
- `variable_dictionary_used.csv`: HiRO dictionary rows used by this parser.
- `preprocessing_report.json`: counts, class balance, unmatched variants, and criteria summary.

## Design Notes

- CASPER WES and VERDICT are kept as separate study labels.
- Variant rows without matching patient data are retained and flagged with
  `patient_match=false`.
- Classifications are normalized into both a five-ish class field and the main
  three-class outcome (`Pathogenic`, `VUS`, `Benign`).
- CASPER WES is preferred over VERDICT only for exact duplicate
  `participant_id + gene + variant` keys, reflecting Brianna's note that CASPER
  WES was reviewed more recently.
- DOB and raw clinical dates are not exported in the normalized files.
- Diagnosis and diagnosis strength are exported but marked as
  sensitivity-only/potential-leakage fields for downstream ML/LLM work.
- Follow-up wide columns are summarized into `any`, `latest`, `yes_count`,
  `uncertain_count`, and `known_count` features.
- Coded categories are read from `Hiro_Full_Variable_List.csv` when available.

## External Annotation

To parse variants and plan external queries without making network calls:

```bash
python3 scripts/annotate_variants.py --dry-run --limit 10
```

To run live annotation with resume/cache support:

```bash
python3 scripts/annotate_variants.py --preferred-labeled-only
python3 scripts/validate_annotations.py
```

Outputs are written to `outputs/annotations/`:

- `variant_annotations.csv`: one annotation row per processed variant.
- `annotated_cohort.csv`: normalized cohort joined to annotation fields.
- `annotation_report.json`: status counts and cache statistics.
- `annotation_cache.jsonl`: append-only cache keyed by normalized gene and
  variant text.

External annotation statuses distinguish `matched`, `no_match`, `absent`,
`skipped`, and `error`. A gnomAD `absent` status is represented with
`gnomad_max_af=0.0` and `gnomad_absent=true`; API errors leave AF blank.

## Frozen Model Inputs

After annotation completes, build the ML/LLM-ready input layer:

```bash
python3 scripts/build_model_inputs.py
python3 scripts/validate_model_inputs.py
```

Outputs are written to `outputs/model_inputs/`:

- `modeling_cohort.csv`: preferred, labeled rows with matched phenotype data.
- `variant_only_cohort.csv`: preferred, labeled rows regardless of phenotype
  match.
- `ml_baseline_features.csv`: primary ML features without ClinVar consensus,
  diagnosis fields, or GC-applied ACMG criteria.
- `ml_enriched_features.csv`: baseline ML features plus structured ClinVar
  consensus/metadata fields.
- `llm_baseline_evidence.jsonl`: structured LLM evidence blocks with safe
  ClinVar metadata but no explicit ClinVar classification cue.
- `llm_enriched_evidence.jsonl`: same evidence blocks with explicit
  `clinvar_label_hint`.
- `criterion_reference.csv`: GC-applied ACMG criteria expanded into
  criterion-level booleans and strength modifiers for reasoning validation.
- `model_input_manifest.json`: row counts, class balance, status coverage,
  missingness, feature lists, and excluded leakage/sensitivity fields.

GC-applied ACMG criteria are intentionally kept out of the primary ML feature
tables. They are used as reference outputs for criterion-level validation and
can be used later in an explicitly labelled sensitivity/ablation analysis.

## CatBoost ML

Run the manuscript-aligned CatBoost baseline and enriched models:

```bash
python3 scripts/run_ml.py
python3 scripts/validate_ml_outputs.py
```

The default configuration uses 5-fold `StratifiedGroupKFold`, grouped by
participant identifier, and fixed CatBoost hyperparameters matching the draft
methods: `iterations=500`, `learning_rate=0.03`, `depth=6`, `random_seed=42`.

Run the leakage audit after building model inputs and ML outputs:

```bash
python3 scripts/audit_model_leakage.py
```

This writes `outputs/qc/model_leakage_audit.json`. Direct use of target,
identifier, GC-criteria, or diagnosis columns and patient fold leakage are
treated as hard failures; high-cardinality HGVS fields are reported as softer
variant-identity warnings.

Explore sensitivity/specificity tradeoffs for pathogenic probability
thresholds:

```bash
python3 scripts/analyze_ml_thresholds.py
```

This writes `outputs/ml_threshold/pathogenic_threshold_sweep.csv` and
`outputs/ml_threshold/pathogenic_threshold_summary.csv`. The default CatBoost
argmax output is conservative; thresholded outputs can be used for triage
operating-point analysis without retraining.

By default, this also renders a formal thresholded triage confusion matrix at
`p_pathogenic >= 0.20`:

- `outputs/ml_threshold/catboost_threshold_0p2_triage_predictions.csv`
- `outputs/ml_threshold/catboost_threshold_0p2_triage_metrics.csv`
- `outputs/ml_threshold/catboost_baseline_threshold_0p2_confusion_matrix.csv`
- `outputs/ml_threshold/catboost_enriched_threshold_0p2_confusion_matrix.csv`
- `outputs/figures/ml/ml_threshold_0p2_triage_confusion_matrices.pdf`
- `outputs/ml_threshold/preferred_ml_triage_model.json`

The preferred ML triage operating point is the enriched CatBoost model with
`p_pathogenic >= 0.20`. The no-exact-variant-identity model below is retained
only as a robustness/sensitivity check, not as the preferred model.

Run the strict exact-variant-identity sensitivity model:

```bash
python3 scripts/run_ml_sensitivity_no_variant_identity.py
```

This removes `pos`, `chrom`, `ref`, `alt`, `gene_refgene`, `hgvs_c`, `hgvs_p`,
and `hgvs_p_one_letter`, then reruns baseline and enriched CatBoost under the
same grouped cross-validation design. Outputs are written to
`outputs/ml_sensitivity/no_exact_variant_identity/`.

Outputs are written to `outputs/ml/`:

- `catboost_baseline_oof_predictions.csv` and
  `catboost_enriched_oof_predictions.csv`: fold-held-out predictions.
- `catboost_oof_predictions.csv`: combined baseline/enriched predictions and
  pathogenic probabilities, ready for hybrid and triage analyses.
- `catboost_metrics.csv` and `catboost_summary.json`: sensitivity,
  specificity, PPV, NPV, MCC, weighted kappa, macro-F1, VUS deferral,
  over-diagnosis, AUC, bootstrap CIs, and exact McNemar comparison.
- `*_confusion_matrix.csv`, `*_fold_metrics.csv`, and
  `*_feature_importance.csv`: audit artifacts for model behavior.

Create manuscript-style ML figures:

```bash
python3 scripts/make_ml_figures.py
python3 scripts/validate_ml_figures.py
```

Figure outputs are written to `outputs/figures/ml/` as both PNG and PDF:

- `ml_confusion_matrices`: side-by-side baseline/enriched three-class
  confusion matrices.
- `ml_metrics_comparison`: sensitivity, specificity, PPV, NPV, MCC, weighted
  kappa, and macro-F1.
- `ml_roc_curves`: binary pathogenic-vs-nonpathogenic ROC curves.
- `ml_precision_recall_curves`: binary pathogenic-vs-nonpathogenic precision
  recall curves.
- `ml_feature_importance_top20`: CatBoost feature importance.
- `ml_shap_pathogenic_top20`: native CatBoost SHAP summary for the pathogenic
  class.
- `ml_shap_pathogenic_values.csv`: mean absolute SHAP values for audit.

## LLM Prompt Preparation

Prepare LLM prompt files without making API calls:

```bash
python3 scripts/prepare_llm_prompts.py
python3 scripts/validate_llm_prompts.py
```

Outputs are written to `outputs/llm/`:

- `primary_baseline_prompts.jsonl`: primary-agent prompts using baseline
  evidence and safe ClinVar metadata only.
- `primary_enriched_prompts.jsonl`: primary-agent prompts with explicit
  ClinVar consensus cue.
- `llm_reference_key.jsonl`: separate answer-key/reference labels used later
  for evaluation joins. These are intentionally not stored inside the API
  prompt-job files.
- `llm_prompt_manifest.json`: prompt counts, schemas, and an explicit
  `llm_api_calls_made=false` flag.

These scripts are intentionally non-spending. They generate prompt jobs and
schemas only; they do not call OpenAI or any other LLM API.

When credits are available, the primary-agent runner is available but dry-runs
by default:

```bash
python3 scripts/run_llm_primary.py --prompt-file outputs/llm/primary_baseline_prompts.jsonl
```

It will not call the API unless `--execute-api` is explicitly supplied.

Paid runs require `OPENAI_API_KEY` in the shell environment or in a project-root
`.env` file. The runner preflights this setting before starting API execution.

The LLM runner is resume-safe. Completed rows are stored one JSON object per
line in the output file, flushed after each variant. If a run is interrupted,
re-running the same command skips successful variants and retries failures by
default:

```bash
python3 scripts/run_llm_primary.py \
  --prompt-file outputs/llm/primary_baseline_prompts.jsonl \
  --output outputs/llm/primary_baseline_outputs.jsonl \
  --execute-api
```

Use the enriched prompt file and output path for Track B:

```bash
python3 scripts/run_llm_primary.py \
  --prompt-file outputs/llm/primary_enriched_prompts.jsonl \
  --output outputs/llm/primary_enriched_outputs.jsonl \
  --execute-api
```

Useful options:

- `--limit 3`: smoke-test only a few rows.
- `--no-retry-failures`: skip previously failed rows instead of retrying them.
- `--no-resume`: start a fresh output file.

Monitor a long run from another terminal:

```bash
python3 scripts/monitor_llm_progress.py \
  --prompt-file outputs/llm/primary_baseline_prompts.jsonl \
  --output outputs/llm/primary_baseline_outputs.jsonl \
  --watch \
  --interval 60
```

Re-run local validation on an existing output file without making API calls:

```bash
python3 scripts/revalidate_llm_outputs.py \
  --input outputs/llm/primary_baseline_outputs.jsonl
```

The monitor and revalidator count latest status by `variant_uid|track`, so
duplicate raw JSONL rows from failed attempts and successful retries are
preserved as an audit trail without being double-counted as final outputs.

Evaluate completed primary LLM outputs offline:

```bash
python3 scripts/evaluate_llm_primary.py \
  --outputs outputs/llm/primary_baseline_outputs.jsonl \
  --reference outputs/llm/llm_reference_key.jsonl \
  --track baseline \
  --output-dir outputs/llm_analysis \
  --figure-dir outputs/figures/llm
```

This writes predictions, validation, metrics, confusion matrices, criterion
summaries, error summaries, and figures without making API calls.

Use the same evaluator for enriched outputs:

```bash
python3 scripts/evaluate_llm_primary.py \
  --outputs outputs/llm/primary_enriched_outputs.jsonl \
  --reference outputs/llm/llm_reference_key.jsonl \
  --track enriched \
  --output-dir outputs/llm_analysis \
  --figure-dir outputs/figures/llm
```

## Adversarial LLM Pathway

The adversarial pathway mirrors the earlier monolithic-code design but is
resume-safe and dry-run-safe in the modular project. It contains a Critic agent
and an Arbitrator agent:

- Critic audits the primary JSON for population-frequency thresholds,
  in-silico PP3 thresholds, and final-label synthesis logic.
- Arbitrator runs only if the critic rejects/modifies the primary output or
  recommends a different label.
- Output is append-only JSONL, flushed after each variant.
- No API calls occur unless `--execute-api` is passed.

Prepare critic prompt files without spending credits:

```bash
python3 scripts/prepare_adversarial_prompts.py \
  --track baseline \
  --evidence outputs/model_inputs/llm_baseline_evidence.jsonl \
  --primary-output outputs/llm/primary_baseline_outputs.jsonl \
  --output outputs/llm/adversarial_baseline_critic_prompts.jsonl

python3 scripts/prepare_adversarial_prompts.py \
  --track enriched \
  --evidence outputs/model_inputs/llm_enriched_evidence.jsonl \
  --primary-output outputs/llm/primary_enriched_outputs.jsonl \
  --output outputs/llm/adversarial_enriched_critic_prompts.jsonl
```

Dry-run the adversarial runner before spending credits:

```bash
python3 scripts/run_llm_adversarial.py \
  --track baseline \
  --evidence outputs/model_inputs/llm_baseline_evidence.jsonl \
  --primary-output outputs/llm/primary_baseline_outputs.jsonl \
  --output outputs/llm/adversarial_baseline_outputs.jsonl
```

Paid adversarial smoke test, only when ready:

```bash
python3 scripts/run_llm_adversarial.py \
  --track baseline \
  --evidence outputs/model_inputs/llm_baseline_evidence.jsonl \
  --primary-output outputs/llm/primary_baseline_outputs.jsonl \
  --output outputs/llm/adversarial_baseline_outputs.jsonl \
  --limit 3 \
  --execute-api
```

If interrupted, rerun the same command. Successful adversarial rows will be
skipped and failures retried by default.

Monitor adversarial progress without making API calls:

```bash
python3 scripts/monitor_adversarial_progress.py \
  --output outputs/llm/adversarial_baseline_outputs.jsonl \
  --expected-jobs 482 \
  --watch \
  --interval 60
```

Evaluate completed adversarial outputs offline:

```bash
python3 scripts/evaluate_llm_adversarial.py \
  --outputs outputs/llm/adversarial_baseline_outputs.jsonl \
  --reference outputs/llm/llm_reference_key.jsonl \
  --track baseline \
  --output-dir outputs/llm_analysis \
  --figure-dir outputs/figures/llm
```

Build the current ML-vs-LLM comparison table/figures:

```bash
python3 scripts/make_model_comparison.py
```

This combines `outputs/ml/catboost_oof_predictions.csv`,
`outputs/llm_analysis/primary_baseline_predictions.csv`, optional
`outputs/llm_analysis/primary_enriched_predictions.csv` if present, optional
`outputs/llm_analysis/adversarial_baseline_predictions.csv` if present, and
optional `outputs/llm_analysis/adversarial_enriched_predictions.csv` if
present, plus the preferred ML threshold predictions at `outputs/ml_threshold/`.
It writes comparison tables to `outputs/model_comparison/` and figures to
`outputs/figures/comparison/`.

Build the hybrid diagnostic confirmation and triage workflow analyses:

```bash
python3 scripts/make_hybrid_analysis.py
```

This is an offline step and makes no API calls. It reimplements the prior
manuscript/SI hybrid logic on the current cohort:

- Diagnostic hybrid: default CatBoost diagnostic predictions are used as the
  base; ML-pathogenic calls require LLM pathogenic confirmation, otherwise the
  hybrid output becomes VUS for human tie-breaking.
- Triage tournament: ML-only comparators are included alongside ML+LLM queues.
  For ML+LLM queues, consensus pathogenic variants are ranked first, ML-only
  pathogenic discordances second, and all other variants lower, with ML
  pathogenic probability used for within-zone ordering.

Key outputs:

- `outputs/hybrid/hybrid_master_predictions.csv`
- `outputs/hybrid/diagnostic_hybrid_metrics.csv`
- `outputs/hybrid/diagnostic_hybrid_metrics.md`
- `outputs/hybrid/diagnostic_hybrid_predictions_long.csv`
- `outputs/hybrid/diagnostic_hybrid_mcnemar.csv`
- `outputs/hybrid/triage_tournament_results.csv`
- `outputs/hybrid/triage_tournament_results_manuscript.csv`
- `outputs/hybrid/triage_winning_queue.csv`
- `outputs/hybrid/triage_zone_summary.csv`
- `outputs/hybrid/hybrid_analysis_manifest.json`
- `outputs/figures/hybrid/diagnostic_hybrid_leaderboard.pdf`
- `outputs/figures/hybrid/triage_tournament_velocity.pdf`
- `outputs/figures/hybrid/triage_cumulative_yield.pdf`
