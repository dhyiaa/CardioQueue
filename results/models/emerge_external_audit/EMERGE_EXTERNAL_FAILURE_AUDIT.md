# eMERGE External Validation Failure Audit

Generated: 2026-07-02

## Question

The first strict source-held-out CatBoost model performed poorly on the eMERGE external split:

| External set | Rows | AUROC | AUPRC |
|---|---:|---:|---:|
| CardioBoost | 225 | 0.9482 | 0.9747 |
| HiRO | 70 | 0.9931 | 0.9904 |
| eMERGE | 203 | 0.5802 | 0.6429 |

This audit checks whether eMERGE was accidentally used in training and why performance collapsed.

## Was eMERGE Used In Training?

No for the source-held-out run.

The split builder assigns clean binary rows with `in_emerge == True` to `external_emerge`, unless the same row is also HiRO, which has higher source priority. The eMERGE external rows in this run are source-pure:

| in_clinvar | in_hiro | in_emerge | in_cardioboost | rows |
|---|---|---|---|---:|
| False | False | True | False | 203 |

So this result is not caused by eMERGE leakage into the source-held-out training set.

## Main Finding

The poor eMERGE result is most likely a coordinate/assembly/variant-representation problem, not a true model-generalization result.

Evidence:

| Feature/annotation check | eMERGE external | HiRO external | CardioBoost external |
|---|---:|---:|---:|
| dbNSFP `ok` | 0 / 203 | 47 / 70 | 207 / 225 |
| dbNSFP `not_found` | 203 / 203 | 20 / 70 | 18 / 225 |
| AlphaMissense direct `ok` | 0 / 203 | 34 / 70 | 188 / 225 |
| SpliceAI `ok` | 17 / 203 | 63 / 70 | 209 / 225 |
| gnomAD confirmed absent | 201 / 203 | 10 / 70 | 17 / 225 |
| VEP MODIFIER impact | 195 / 203 | 2 / 70 | 17 / 225 |
| VEP missense fraction | 1.5% | 61.4% | 90.2% |

The original eMERGE source table says the binary eMERGE variants are mostly coding or splice-relevant:

| eMERGE original inferred consequence | Benign | Pathogenic |
|---|---:|---:|
| missense | 117 | 44 |
| frameshift | 0 | 45 |
| nonsense | 0 | 20 |
| splice_region | 2 | 16 |
| Missing | 18 | 2 |

But after our GRCh38 feature annotation pipeline, the held-out eMERGE variants are mostly VEP intronic/intergenic/MODIFIER. That is biologically inconsistent with the source table and explains the model behavior.

## Likely Cause

The eMERGE source coordinates appear to have been handled as if they were GRCh38 in the current registry/model matrix, but the eMERGE provenance code treats the source coordinates as hg19/GRCh37 and attempts MyVariant/dbNSFP rescue to hg38.

The relevant provenance code has fields such as:

- `source_hg19_chrom`
- `source_hg19_pos`
- `source_hg19_ref`
- `source_hg19_alt`
- `retrieved_hg38_*`
- `final_coordinate_source`

However, the local final eMERGE handoff file does not preserve `final_coordinate_source`, so rows that failed MyVariant rescue may have silently fallen back to source hg19 coordinates. Those fallback coordinates then entered the global GRCh38 matrix and were annotated against GRCh38 VEP/dbNSFP/gnomAD, causing false intronic/intergenic calls and total dbNSFP misses.

## Why The Model Scored eMERGE Low

The strict model learned strong signals from VEP consequence/impact, gnomAD frequency, SpliceAI, dbNSFP, AlphaMissense/ESM, and protein features.

For eMERGE external rows, many of those signals were missing or benign-like because the coordinates likely pointed to the wrong GRCh38 loci:

- all eMERGE external rows had dbNSFP `not_found`
- all AlphaMissense direct scores were missing
- most SpliceAI values were missing
- almost every row was marked gnomAD-confirmed absent
- VEP mostly called them intronic/intergenic/MODIFIER

So the low probabilities are expected under the current, likely-wrong annotation state.

## Files Written By This Audit

- `results/models/emerge_external_audit/emerge_external_audit.summary.json`
- `results/models/emerge_external_audit/external_emerge_rows_with_predictions.tsv`
- `results/models/emerge_external_audit/external_emerge_by_gene_label_prob.tsv`
- `results/models/emerge_external_audit/external_emerge_false_negatives_threshold_0_5.tsv`
- `results/models/emerge_external_audit/external_emerge_false_positives_threshold_0_5.tsv`
- `results/models/emerge_external_audit/external_emerge_source_overlap.tsv`

## Recommended Fix Before Interpreting eMERGE

1. Quarantine current eMERGE external model result as `invalid until coordinate rescue`.
2. Rebuild eMERGE coordinates with explicit provenance:
   - keep original hg19 fields
   - derive GRCh38 via a reproducible liftOver chain or trusted MyVariant/dbNSFP mapping
   - store `emerge_coordinate_source`, `emerge_liftover_status`, and `emerge_coordinate_qc_flag`
3. Validate rescued eMERGE rows before model training:
   - VEP consequence should agree with source `inferred_consequence_model` for most coding/splice variants
   - dbNSFP should recover missense/nonsense scores where applicable
   - AlphaMissense should map for missense rows
   - gnomAD should not be `confirmed_absent` for nearly all benign rows when source gnomAD AF is nonzero
4. Rebuild registry, feature joins, splits, and source-held-out model.
5. Only then report eMERGE as an external validation set.

## Interim Interpretation

The current eMERGE poor result should not be interpreted as a true model failure. It is best interpreted as a successful QC alarm: the source-held-out split exposed that eMERGE is not currently harmonized to the same coordinate/annotation space as ClinVar, HiRO, and CardioBoost.
