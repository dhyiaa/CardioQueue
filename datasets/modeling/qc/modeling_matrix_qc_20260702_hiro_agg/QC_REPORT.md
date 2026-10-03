# Modeling Matrix QC Report

Generated: `2026-07-02T15:20:09`
Matrix: `datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx_vep_hgvs_spliceai_hiro_agg.tsv`

## Executive Summary

- Rows: `86,889`
- Columns: `463`
- Duplicate `variant_id` rows: `0`
- Duplicate `coordinate_key` rows: `0`
- `ref == alt` no-op rows: `17`
- Binary supervised rows, before final split decisions: `43,282`
- Primary-include binary supervised rows: `43,217`

## Label Balance

| model_label_3class | rows | percent |
| --- | --- | --- |
| VUS | 43458 | 50.016 |
| Benign | 34042 | 39.179 |
| Pathogenic | 9240 | 10.634 |
| <missing> | 149 | 0.171 |

## Source Coverage

| sources | rows | percent |
| --- | --- | --- |
| clinvar | 83546 | 96.153 |
| emerge | 2752 | 3.167 |
| cardioboost\|clinvar | 234 | 0.269 |
| clinvar\|hiro | 131 | 0.151 |
| cardioboost | 115 | 0.132 |
| hiro | 105 | 0.121 |
| cardioboost\|clinvar\|hiro | 2 | 0.002 |
| clinvar\|emerge | 2 | 0.002 |
| cardioboost\|hiro | 2 | 0.002 |

## Source Record Preservation

| source | path | status | source_rows | distinct_source_records | linked_source_rows | unresolved_source_rows | unique_resolved_variants |
| --- | --- | --- | --- | --- | --- | --- | --- |
| HiRO | datasets/hiro/full_dataset/interim/hiro_linked_source_records.tsv | ok | 482 | 482 | 408 | 73 | 239 |
| eMERGE | datasets/emerge/full_arrhythmia_gene_dataset/data/emerge_linked_source_records.tsv | ok | 2754 | 2754 | 2754 | 0 | 2754 |
| CardioBoost | datasets/cardioboost/public_dataset/processed_for_cardiogenetics/model_inputs/cardioboost_linked_source_records.tsv | ok | 355 | 355 | 355 | 0 | 353 |

## Split Roles And Scope

| default_split_role | rows | percent |
| --- | --- | --- |
| training_candidate_unassigned | 83546 | 96.153 |
| emerge_candidate_validation_or_source_stratum | 2092 | 2.408 |
| excluded_gene_scope | 695 | 0.8 |
| cardioboost_binary_source_stratum | 225 | 0.259 |
| hiro_high_trust_stratum | 182 | 0.209 |
| label_conflict_review | 149 | 0.171 |

| modeling_scope | rows | percent |
| --- | --- | --- |
| primary | 73536 | 84.632 |
| expanded_scope | 9254 | 10.65 |
| primary_monitor | 3403 | 3.916 |
| sensitivity_only | 680 | 0.783 |
| exclude_or_quarantine | 13 | 0.015 |
| ttn_separate | 3 | 0.003 |

## Label Conflicts

| label_conflict_type | rows | percent |
| --- | --- | --- |
| none | 86740 | 99.829 |
| benign_vus_conflict | 86 | 0.099 |
| pathogenic_vus_conflict | 61 | 0.07 |
| benign_pathogenic_conflict | 2 | 0.002 |

Top conflict rows are written to `label_conflict_rows.tsv`.

## Feature Group Coverage

| feature_group | columns_present | status_columns | value_columns | rows_with_positive_status | positive_status_percent | rows_with_any_value | value_coverage_percent |
| --- | --- | --- | --- | --- | --- | --- | --- |
| dbNSFP | 12 | 1 | 11 | 42338 | 48.727 | 42338 | 48.727 |
| AlphaMissense direct | 3 | 1 | 1 | 36773 | 42.322 | 36773 | 42.322 |
| gnomAD | 5 | 2 | 3 | 71448 | 82.229 | 71448 | 82.229 |
| ClinGen | 5 | 3 | 2 | 84983 | 97.806 | 84983 | 97.806 |
| Protein structure | 9 | 5 | 3 | 3089 | 3.555 | 2179 | 2.508 |
| VEP | 6 | 2 | 4 | 86872 | 99.98 | 86872 | 99.98 |
| SpliceAI | 5 | 1 | 4 | 75954 | 87.415 | 75954 | 87.415 |

## Core Feature Missingness

| column | present | nonmissing | missing | missing_percent | unique_nonmissing |
| --- | --- | --- | --- | --- | --- |
| alphafold_plddt_status | True | 86889 | 0 | 0.0 | 4 |
| alphafold_residue_plddt | True | 2179 | 84710 | 97.492 | 988 |
| alphamissense_dbnsfp_score | True | 41877 | 45012 | 51.804 | 35708 |
| alphamissense_direct_class | True | 36773 | 50116 | 57.678 | 3 |
| alphamissense_direct_score | True | 36773 | 50116 | 57.678 | 8590 |
| alphamissense_direct_status | True | 86889 | 0 | 0.0 | 3 |
| cadd_phred | True | 42338 | 44551 | 51.273 | 3024 |
| clingen_dosage_status | True | 86886 | 3 | 0.003 | 2 |
| clingen_gene_validity_max_classification | True | 84983 | 1906 | 2.194 | 5 |
| clingen_gene_validity_status | True | 86889 | 0 | 0.0 | 3 |
| clingen_variant_evidence_classification | True | 319 | 86570 | 99.633 | 5 |
| clingen_variant_evidence_status | True | 86889 | 0 | 0.0 | 3 |
| dbnsfp_status | True | 86889 | 0 | 0.0 | 3 |
| dssp_secondary_structure_class | True | 2178 | 84711 | 97.493 | 6 |
| dssp_status | True | 86889 | 0 | 0.0 | 4 |
| esm1b_score | True | 42075 | 44814 | 51.576 | 40158 |
| fathmm_xf_coding_score | True | 38636 | 48253 | 55.534 | 35922 |
| foldx_ddg_kcal_mol | True | 1099 | 85790 | 98.735 | 1064 |
| foldx_ddg_status | True | 86889 | 0 | 0.0 | 3 |
| freesasa_relative | True | 2178 | 84711 | 97.493 | 1637 |
| freesasa_status | True | 86889 | 0 | 0.0 | 4 |
| gerp_rs | True | 42334 | 44555 | 51.278 | 2252 |
| gnomad_browser_status | True | 86889 | 0 | 0.0 | 3 |
| gnomad_final_af | True | 64472 | 22417 | 25.8 | 28863 |
| gnomad_final_homozygote_count | True | 64477 | 22412 | 25.794 | 2026 |
| gnomad_final_popmax_af | True | 71443 | 15446 | 17.777 | 27565 |
| gnomad_final_status | True | 86889 | 0 | 0.0 | 3 |
| metalr_score | True | 38319 | 48570 | 55.899 | 9536 |
| phastcons100way_vertebrate | True | 42338 | 44551 | 51.273 | 962 |
| phylop100way_vertebrate | True | 42338 | 44551 | 51.273 | 9346 |
| polyphen2_hdiv_score | True | 41709 | 45180 | 51.997 | 18144 |
| protein_feature_status | True | 86889 | 0 | 0.0 | 4 |
| revel_score | True | 41719 | 45170 | 51.986 | 24403 |
| sift_score | True | 41734 | 45155 | 51.969 | 19641 |
| spliceai_status | True | 86889 | 0 | 0.0 | 2 |
| vep_SpliceAI_pred_DS_AG | True | 75954 | 10935 | 12.585 | 101 |
| vep_SpliceAI_pred_DS_AL | True | 75954 | 10935 | 12.585 | 101 |
| vep_SpliceAI_pred_DS_DG | True | 75954 | 10935 | 12.585 | 101 |
| vep_SpliceAI_pred_DS_DL | True | 75954 | 10935 | 12.585 | 101 |
| vep_annotation_status | True | 86889 | 0 | 0.0 | 2 |
| vep_hgvsc | True | 86872 | 17 | 0.02 | 85549 |
| vep_hgvsp | True | 86872 | 17 | 0.02 | 61209 |
| vep_impact | True | 86872 | 17 | 0.02 | 4 |
| vep_status | True | 86889 | 0 | 0.0 | 2 |
| vep_worst_consequence | True | 86872 | 17 | 0.02 | 23 |

## Top Genes By Label Count

| primary_gene | Benign | Pathogenic | VUS | nan | labeled_3class_rows | supervised_binary_rows | pathogenic_fraction_binary | sparse_pathogenic_lt30 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RYR2 | 4009 | 270 | 5553 | 2 | 9832 | 4279 | 0.0631 | False |
| FLNC | 2361 | 608 | 2702 | 0 | 5671 | 2969 | 0.2048 | False |
| DSP | 1849 | 817 | 2665 | 2 | 5331 | 2666 | 0.3065 | False |
| MYH7 | 1910 | 423 | 2869 | 26 | 5202 | 2333 | 0.1813 | False |
| SCN5A | 1625 | 554 | 2684 | 26 | 4863 | 2179 | 0.2542 | False |
| MYBPC3 | 1306 | 1041 | 1872 | 17 | 4219 | 2347 | 0.4435 | False |
| CACNA1C | 1920 | 100 | 1666 | 12 | 3686 | 2020 | 0.0495 | False |
| KCNH2 | 1193 | 715 | 1731 | 19 | 3639 | 1908 | 0.3747 | False |
| KCNQ1 | 1395 | 568 | 928 | 10 | 2891 | 1963 | 0.2894 | False |
| HCN4 | 910 | 20 | 1240 | 0 | 2170 | 930 | 0.0215 | True |
| LMNA | 625 | 535 | 987 | 6 | 2147 | 1160 | 0.4612 | False |
| PKP2 | 684 | 341 | 1062 | 2 | 2087 | 1025 | 0.3327 | False |
| DSG2 | 683 | 129 | 1222 | 0 | 2034 | 812 | 0.1589 | False |
| RBM20 | 863 | 103 | 1034 | 0 | 2000 | 966 | 0.1066 | False |
| DSC2 | 649 | 98 | 1104 | 2 | 1851 | 747 | 0.1312 | False |
| GLA | 301 | 1099 | 425 | 1 | 1825 | 1400 | 0.785 | False |
| TRDN | 819 | 55 | 590 | 0 | 1464 | 874 | 0.0629 | False |
| PRKAG2 | 657 | 26 | 656 | 3 | 1339 | 683 | 0.0381 | True |
| JUP | 585 | 36 | 710 | 0 | 1331 | 621 | 0.058 | False |
| BAG3 | 376 | 154 | 733 | 0 | 1263 | 530 | 0.2906 | False |
| RAF1 | 574 | 51 | 627 | 0 | 1252 | 625 | 0.0816 | False |
| DES | 409 | 122 | 673 | 0 | 1204 | 531 | 0.2298 | False |
| SLC22A5 | 517 | 246 | 419 | 1 | 1182 | 763 | 0.3224 | False |
| PTPN11 | 481 | 136 | 467 | 1 | 1084 | 617 | 0.2204 | False |
| CACNB2 | 529 | 1 | 529 | 0 | 1059 | 530 | 0.0019 | True |
| TMEM43 | 450 | 2 | 594 | 1 | 1046 | 452 | 0.0044 | True |
| TNNT2 | 449 | 68 | 498 | 2 | 1015 | 517 | 0.1315 | False |
| TPM1 | 441 | 52 | 458 | 7 | 951 | 493 | 0.1055 | False |
| MAP2K2 | 498 | 19 | 428 | 0 | 945 | 517 | 0.0368 | True |
| ACTC1 | 375 | 22 | 432 | 1 | 829 | 397 | 0.0554 | True |
| CASQ2 | 402 | 78 | 327 | 0 | 807 | 480 | 0.1625 | False |
| TNNI3 | 291 | 82 | 392 | 4 | 765 | 373 | 0.2198 | False |
| LAMP2 | 321 | 132 | 282 | 0 | 735 | 453 | 0.2914 | False |
| ANK2 | 48 | 13 | 613 | 0 | 674 | 61 | 0.2131 | True |
| MIB1 | 265 | 27 | 350 | 0 | 642 | 292 | 0.0925 | True |
| SCN1B | 234 | 37 | 365 | 0 | 636 | 271 | 0.1365 | False |
| HRAS | 331 | 29 | 274 | 0 | 634 | 360 | 0.0806 | True |
| KCNJ2 | 204 | 63 | 363 | 0 | 630 | 267 | 0.236 | False |
| SNTA1 | 278 | 0 | 336 | 0 | 614 | 278 | 0.0 | True |
| MYL2 | 252 | 15 | 321 | 0 | 588 | 267 | 0.0562 | True |

## Sparse Pathogenic Genes

| primary_gene | Benign | Pathogenic | VUS | nan | labeled_3class_rows | supervised_binary_rows | pathogenic_fraction_binary | sparse_pathogenic_lt30 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SNTA1 | 278 | 0 | 336 | 0 | 614 | 278 | 0.0 | True |
| CACNB2 | 529 | 1 | 529 | 0 | 1059 | 530 | 0.0019 | True |
| FPGT-TNNI3K | 0 | 1 | 3 | 0 | 4 | 1 | 1.0 | True |
| DMD | 0 | 1 | 2 | 0 | 3 | 1 | 1.0 | True |
| TMEM43 | 450 | 2 | 594 | 1 | 1046 | 452 | 0.0044 | True |
| TTN | 0 | 2 | 1 | 0 | 3 | 2 | 1.0 | True |
| KCNE2 | 54 | 3 | 97 | 0 | 154 | 57 | 0.0526 | True |
| MYL3 | 163 | 5 | 295 | 0 | 463 | 168 | 0.0298 | True |
| FHOD3 | 173 | 6 | 272 | 1 | 451 | 179 | 0.0335 | True |
| CALM3 | 143 | 7 | 50 | 0 | 200 | 150 | 0.0467 | True |
| PLN | 43 | 11 | 88 | 1 | 142 | 54 | 0.2037 | True |
| NRAS | 121 | 12 | 154 | 0 | 287 | 133 | 0.0902 | True |
| ANK2 | 48 | 13 | 613 | 0 | 674 | 61 | 0.2131 | True |
| MYL2 | 252 | 15 | 321 | 0 | 588 | 267 | 0.0562 | True |
| CALM1 | 123 | 18 | 38 | 0 | 179 | 141 | 0.1277 | True |
| MAP2K2 | 498 | 19 | 428 | 0 | 945 | 517 | 0.0368 | True |
| HCN4 | 910 | 20 | 1240 | 0 | 2170 | 930 | 0.0215 | True |
| CALM2 | 136 | 20 | 66 | 0 | 222 | 156 | 0.1282 | True |
| ACTC1 | 375 | 22 | 432 | 1 | 829 | 397 | 0.0554 | True |
| KCNE1 | 128 | 23 | 211 | 1 | 362 | 151 | 0.1523 | True |
| PRKAG2 | 657 | 26 | 656 | 3 | 1339 | 683 | 0.0381 | True |
| MIB1 | 265 | 27 | 350 | 0 | 642 | 292 | 0.0925 | True |
| HRAS | 331 | 29 | 274 | 0 | 634 | 360 | 0.0806 | True |

## Main QC Flags

- Rows with `label_conflict_type != none` should be excluded from clean supervised training or reviewed.
- `ref == alt` rows explain the VEP-missing rows and should not be treated as annotation failures.
- HiRO/eMERGE/CardioBoost source records are preserved separately from the variant-level matrix.
- `gnomad_final_status == confirmed_absent` is distinct from `not_joined`; only confirmed absent can support an AF=0 indicator.
- FoldX `ref_mismatch` rows should keep DDG missing and use the mismatch flag rather than imputation.

## Output Files

- `all_status_missingness.tsv`
- `clean_supervised_label_counts.tsv`
- `clingen_gene_validity_status_counts.tsv`
- `clingen_variant_evidence_status_counts.tsv`
- `core_feature_missingness.tsv`
- `default_split_role_counts.tsv`
- `dssp_status_counts.tsv`
- `duplicate_variant_id_rows.tsv`
- `feature_group_coverage.tsv`
- `foldx_status_counts.tsv`
- `freesasa_status_counts.tsv`
- `gene_label_counts.tsv`
- `gene_panel_decision_counts.tsv`
- `gnomad_final_status_counts.tsv`
- `label_conflict_counts.tsv`
- `label_conflict_rows.tsv`
- `label_counts.tsv`
- `label_counts_by_default_split_role.tsv`
- `label_counts_by_gene_panel_decisions.tsv`
- `label_counts_by_modeling_scope.tsv`
- `label_counts_by_primary_gene.tsv`
- `label_counts_by_sources.tsv`
- `linked_source_record_duplicate_variants.tsv`
- `linked_source_record_summary.tsv`
- `modeling_scope_counts.tsv`
- `multi_source_record_variant_rows.tsv`
- `protein_feature_status_counts.tsv`
- `ref_alt_equal_rows.tsv`
- `source_overlap_counts.tsv`
- `sparse_pathogenic_genes.tsv`
- `spliceai_status_counts.tsv`
- `vep_status_counts.tsv`
- `summary.json`
