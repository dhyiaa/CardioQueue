# Project Context Notes

These notes summarize the working assumptions from the March-April 2026 email
chain and meeting discussion. They intentionally do not store CWL credentials,
Zoom passcodes, or other operational access details.

## Data Context

- The original 361-variant analysis should be framed carefully because the
  reference labels may reflect commercial/lab interpretations at time of
  reporting, including older records, rather than a uniform contemporary
  cardiogenetics re-interpretation.
- The April 2026 update contains 512 classified variants from two studies:
  CASPER WES and VERDICT.
- CASPER WES classifications should be preferred for exact duplicate variants
  between CASPER WES and VERDICT because CASPER WES was reviewed more recently.
- The new classified-variant files include GC/research-team applied ACMG
  criteria, which should replace ClinVar-derived criterion activation wherever
  possible.
- Modified ACMG criterion strengths are represented with underscore notation,
  such as `PVS1_Moderate`.
- Doug's patient files are pivoted to one row per patient. Repeated ECG, echo,
  and follow-up events appear as wide columns such as `QT_1` through `QT_25`.
- The HiRO variable dictionary file is the key for interpreting coded fields.

## Modeling Context

- Diagnosis and diagnosis strength should be retained for audit and sensitivity
  analyses but excluded from primary ML/LLM feature sets unless explicitly
  testing a sensitivity model. They may encode downstream clinical
  interpretation and can create leakage.
- Safe primary feature candidates include variant identity, GC-applied ACMG
  criteria, objective symptoms, family history flags, ECG summaries, echo
  summaries, sex, and age.
- ClinVar should not be treated as the source of truth for ACMG criteria in the
  new pipeline because Brianna flagged errors in ClinVar algorithmic criterion
  extraction.
- ClinGen gene-specific ACMG rules may matter for later prompt/model refinement,
  especially strength modifications like MYH7 PVS1 being moderate rather than
  very strong.

## Compute Context

- ARC/Sockeye or related GPU access may become useful for large-scale reruns,
  but the preprocessing layer is intentionally local and lightweight.

## Current Analysis Status

- The April 2026 ML pipeline has been rerun end to end on 482
  variant-to-phenotype rows.
- Default CatBoost argmax is retained as the diagnostic ML comparator.
- Enriched CatBoost with `p_pathogenic >= 0.20` is the preferred ML triage
  operating point. It is a review-prioritization signal, not autonomous final
  diagnosis.
- LLM prompt jobs are prepared and validated for baseline and enriched primary
  LLM runs. The answer-key reference labels are stored separately from API
  prompt-job files.
- The next major project step is to run primary LLM baseline/enriched outputs
  when credits and time are available, then evaluate LLM and hybrid ML-LLM
  workflows.

See `docs/current_status_and_next_steps.md` for the concise handoff summary and
exact run commands.
