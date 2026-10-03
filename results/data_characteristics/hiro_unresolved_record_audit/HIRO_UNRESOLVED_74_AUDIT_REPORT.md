# HiRO Unresolved Source-Record Audit

This audit explains the HiRO records that still do not link to a promoted variant-registry row. It preserves every HiRO source record and does not modify the registry or modeling matrix.

## Summary

- HiRO source records: 482
- Linked to registry: 408
- Unresolved or not promoted: 74
- Unique internal gene+cDNA candidates: 2
- Unique internal protein-only candidates: 1
- ClinVar review candidates from the previous rescue pass: 3

## Status By Label

| hiro_registry_link_status | Benign | Pathogenic | VUS |
| --- | --- | --- | --- |
| rescued_but_not_in_promoted_registry | 0 | 0 | 1 |
| unresolved_no_coordinates | 30 | 9 | 34 |

## Main Blockers

| blocker_category | rows |
| --- | --- |
| needs_hgvs_normalization_missense | 26 |
| needs_hgvs_normalization_indel_splice_lof | 22 |
| no_variant_text_or_hgvs | 15 |
| gene_symbol_or_readthrough_cleanup | 5 |
| clinvar_review_candidate_not_promoted | 3 |
| partial_coordinate_or_allele | 1 |
| needs_manual_review | 1 |
| rescued_coordinate_not_promoted | 1 |

## Recommended Actions

| recommended_next_action | rows |
| --- | --- |
| Run online VEP HGVS, VariantValidator, Mutalyzer, or a local transcript-aware HGVS normalizer. | 45 |
| Cannot rescue from current processed table; return to raw source/report. | 15 |
| Clean gene alias/readthrough first, then run HGVS normalization. | 5 |
| Manual ClinVar/HGVS review; do not auto-promote. | 3 |
| Review for promotion: unique internal matrix gene+cDNA candidate. | 2 |
| Manual review. | 2 |
| Recover missing allele fields from original VCF/report or normalize HGVS. | 1 |
| Manual review only: unique protein-level candidate may be isoform-dependent. | 1 |

## Review Candidates

| variant_uid | target_3class | gene | hgvs_c | hgvs_p | inferred_consequence | matrix_cdna_candidate_variant_ids | matrix_protein_candidate_variant_ids | rescue_clinvar_name | recommended_next_action |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASPER_WES|263|FHOD3|C.G646A:P.V216I|67 | VUS | FHOD3 | c.646G>A | p.V216I | missense | 18-36611991-G-A |  | nan | Review for promotion: unique internal matrix gene+cDNA candidate. |
| CASPER_WES|6087|FLNC|C.G3301A:P.G1101S|76 | VUS | FLNC | c.3301G>A | p.G1101S | missense | 7-128845054-G-A |  | nan | Review for promotion: unique internal matrix gene+cDNA candidate. |
| VERDICT|5924|RYR2|C.60225G>A|30 | VUS | RyR2 | c.6022+5 | Missing | splice_region |  |  | NM_001035.3(RYR2):c.6022+5G>A | Manual ClinVar/HGVS review; do not auto-promote. |
| VERDICT|1787|CASQ2|C.1144_1149DUPGAT>P.ASP383DUP|289 | VUS | CASQ2 | c.1144_1149dupGAT> | p.Asp383dup | indel |  |  | NM_001232.4(CASQ2):c.1132GAT[7] (p.Asp383dup) | Manual ClinVar/HGVS review; do not auto-promote. |
| VERDICT|3247|DSC2|C.2187G>AP.ALA733THR->2197|299 | VUS | DSC2 | c.2187G>A | p.Ala733Thr | missense |  | 18-31070779-C-T | nan | Manual review only: unique protein-level candidate may be isoform-dependent. |
| VERDICT|5486|AKAP9|C.4825_4826DELAGINSCA>P.ARG1609GLN|325 | VUS | AKAP9 | c.4825_4826delAGinsCA> | p.Arg1609Gln | indel |  |  | NM_005751.5(AKAP9):c.4825_4826delinsCA (p.Arg1609Gln) | Manual ClinVar/HGVS review; do not auto-promote. |

## Interpretation

Most unresolved rows have HGVS-like variant text but no safe genomic VCF key. The largest blocker is transcript-aware HGVS normalization, especially for missense, indel, frameshift, and splice records. Offline VEP cannot parse HGVS input, so the remaining rescue work needs online VEP, VariantValidator, Mutalyzer, or a local transcript-aware HGVS normalizer.

The current audit finds two unique internal gene+cDNA candidates and one protein-only internal candidate. These should be reviewed before promotion because transcript choice changes cDNA numbering for some genes. The three older ClinVar review candidates remain manual-review candidates and should not be automatically promoted.

## Files

- `tables/hiro_unresolved_74_record_audit.tsv`
- `tables/hiro_unresolved_review_candidates.tsv`
- `tables/blocker_by_label.tsv`
- `tables/blocker_by_consequence.tsv`
- `hiro_unresolved_74_audit.summary.json`
