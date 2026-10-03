# CardioQueue v1.0

CardioQueue ranks variants in inherited cardiac disease genes for expert evidence review. It does not assign ACMG/AMP classifications and must not be used as a stand-alone basis for diagnosis, treatment, family testing, or emergency triage.

## Included

- `model/cardioqueue_v1.0.cbm`: frozen CatBoost model.
- `model/feature_columns.json`: ordered 300-feature schema and categorical fields.
- `model/training_metrics.json`: development and retrospective evaluation summary.
- `model/feature_importance.tsv`: CatBoost feature importance.
- `worklists/public_source_vus_research_worklist.tsv`: ranked 41,622-row research worklist with every HiRO-linked row excluded.
- `score_variants.py`: batch inference with strict schema checks.
- `FEATURE_PROVENANCE.md`: source, calculation, and deployment requirements for each feature family.
- `requirements.txt`: minimal inference environment.

The public model excludes 38 HiRO-prefixed fields. In retrospective source-held-out evaluation, its AUROC was 0.978 when the three held-out sources were pooled. Performance, calibration, and review yield can change with laboratory ascertainment, ancestry, gene mix, annotation versions, and local prevalence.

## Intended workflow

Use the score after a laboratory's usual variant-quality control and annotation. It can order molecular records for literature refresh, transcript reconciliation, segregation or functional follow-up, database-discordance review, and multidisciplinary preparation. The intended outcome is a transparent evidence-review queue, not an ACMG/AMP classification, diagnosis, patient-urgency score, or treatment recommendation.

## Input

Provide a tab-delimited table containing `variant_id` and every field listed in `model/feature_columns.json`. Values may be missing. Categorical missing values are converted to `__MISSING__`; numeric fields are parsed as numbers.

```bash
python score_variants.py annotated_variants.tsv cardioqueue_scores.tsv
```

The output contains `variant_id`, `cardioqueue_version`, and `review_priority_score`. Higher scores place a variant earlier in the review queue. The score is not a calibrated probability of pathogenicity in a new clinical population.

## Reproducibility and data

The manuscript repository will provide analysis scripts, frozen schemas, checksums, redistributable derived tables, and environment records. ClinVar, ClinGen, gnomAD, dbNSFP, AlphaMissense, UniProt, AlphaFold, DSSP, FreeSASA, and CardioBoost remain subject to their source terms and versioning. FoldX must be obtained separately under its own license and is not redistributed.

## Citation

Citation metadata will be added after manuscript acceptance. Cite the archived release DOI and manuscript together once available.
