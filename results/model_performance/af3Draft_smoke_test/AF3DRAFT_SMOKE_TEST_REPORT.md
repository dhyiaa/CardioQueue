# AF3Draft Smoke Test

This is a prototype analysis. The current AlphaFold 3 layer is named **AF3Draft** because the partner list and AF3 job coverage are not final.

## What Was Tested

- Base model: current primary weighted binary CatBoost model.
- Smoke-test models: same CatBoost setup and same `split_source_heldout` split, with only added `af3Draft_` features.
- AF3 job IDs and complex names were kept for QC but excluded from training features.
- A stricter AF3Draft model also excludes AF3 status/count/partner-gene provenance fields and keeps only structure-like columns.
- This test is not a final claim. It asks whether the draft AF3 pathway has enough signal to justify scaling tomorrow.

## AF3Draft Coverage

| scope | rows | genes | rows_with_any_af3Draft_output | rows_with_strong_or_moderate_af3Draft_pair | rows_with_strong_af3Draft_pair | pct_all_rows | pct_binary_rows | pct_vus_rows |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all_modeling_rows | 85677 | 61 | 32143 | 13884 | 7220 | 37.52% | 32.30% | 32.78% |
| primary_binary_rows | 42990 | 53 | 15557 | 6762 | 3479 | 36.19% | 15.73% | 15.96% |
| vus_rows | 42361 | 61 | 16461 | 7059 | 3684 | 38.86% | 16.42% | 16.66% |

More precise structural coverage:

| scope | rows | trusted AF3Draft pair context | trusted pair plus residue feature | trusted pair plus 5A interface | trusted pair plus 8A interface |
|---|---:|---:|---:|---:|---:|
| all_modeling_rows | 85,677 | 13,884 | 11,048 | 1,723 | 3,075 |
| primary_binary_rows | 42,990 | 6,762 | 4,650 | 732 | 1,276 |
| vus_rows | 42,361 | 7,059 | 6,339 | 986 | 1,792 |

The trusted-pair count means the variant's gene has a strong/moderate AF3Draft interaction context. The residue-feature and interface columns are stricter variant-level structural coverage.

## Performance Comparison

| model | subset | scope | rows | auroc | auprc | sensitivity | specificity | ppv | npv |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| baseline_primary_binary_existing | external_cardioboost | all_rows | 185 | 0.9520 | 0.9812 | 0.9545 | 0.8113 | 0.9265 | 0.8776 |
| baseline_primary_binary_existing | external_cardioboost | af3Draft_trusted_pair_rows | 35 | 0.9400 | 0.9703 | 1.0000 | 0.7000 | 0.8929 | 1.0000 |
| baseline_primary_binary_existing | external_emerge | all_rows | 176 | 0.9910 | 0.9914 | 0.9792 | 0.9375 | 0.9495 | 0.9740 |
| baseline_primary_binary_existing | external_emerge | af3Draft_trusted_pair_rows | 42 | 0.9886 | 0.9922 | 0.9565 | 1.0000 | 1.0000 | 0.9500 |
| baseline_primary_binary_existing | external_hiro | all_rows | 69 | 0.9877 | 0.9830 | 0.9630 | 0.8810 | 0.8387 | 0.9737 |
| baseline_primary_binary_existing | external_hiro | af3Draft_trusted_pair_rows | 21 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| baseline_primary_binary_existing | train | all_rows | 36176 | 0.9998 | 0.9992 | 0.9932 | 0.9978 | 0.9918 | 0.9982 |
| baseline_primary_binary_existing | train | af3Draft_trusted_pair_rows | 5665 | 0.9999 | 0.9997 | 0.9932 | 0.9991 | 0.9970 | 0.9979 |
| baseline_primary_binary_existing | validation | all_rows | 6384 | 0.9999 | 0.9996 | 0.9947 | 0.9970 | 0.9888 | 0.9986 |
| baseline_primary_binary_existing | validation | af3Draft_trusted_pair_rows | 999 | 1.0000 | 0.9999 | 0.9915 | 0.9987 | 0.9957 | 0.9974 |
| af3Draft_plus_primary_binary_fullDraft | train | all_rows | 36176 | 0.9998 | 0.9994 | 0.9938 | 0.9982 | 0.9932 | 0.9984 |
| af3Draft_plus_primary_binary_fullDraft | train | af3Draft_trusted_pair_rows | 5665 | 0.9999 | 0.9998 | 0.9932 | 0.9991 | 0.9970 | 0.9979 |
| af3Draft_plus_primary_binary_fullDraft | validation | all_rows | 6384 | 0.9999 | 0.9996 | 0.9962 | 0.9970 | 0.9888 | 0.9990 |
| af3Draft_plus_primary_binary_fullDraft | validation | af3Draft_trusted_pair_rows | 999 | 1.0000 | 0.9999 | 0.9957 | 0.9974 | 0.9915 | 0.9987 |
| af3Draft_plus_primary_binary_fullDraft | external_cardioboost | all_rows | 185 | 0.9565 | 0.9837 | 0.9545 | 0.8113 | 0.9265 | 0.8776 |
| af3Draft_plus_primary_binary_fullDraft | external_cardioboost | af3Draft_trusted_pair_rows | 35 | 0.9640 | 0.9850 | 1.0000 | 0.7000 | 0.8929 | 1.0000 |
| af3Draft_plus_primary_binary_fullDraft | external_emerge | all_rows | 176 | 0.9887 | 0.9873 | 0.9792 | 0.9500 | 0.9592 | 0.9744 |
| af3Draft_plus_primary_binary_fullDraft | external_emerge | af3Draft_trusted_pair_rows | 42 | 0.9840 | 0.9899 | 0.9565 | 0.9474 | 0.9565 | 0.9474 |
| af3Draft_plus_primary_binary_fullDraft | external_hiro | all_rows | 69 | 0.9885 | 0.9826 | 0.9630 | 0.8810 | 0.8387 | 0.9737 |
| af3Draft_plus_primary_binary_fullDraft | external_hiro | af3Draft_trusted_pair_rows | 21 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| af3Draft_plus_primary_binary_strictStructural | train | all_rows | 36176 | 0.9998 | 0.9994 | 0.9938 | 0.9979 | 0.9922 | 0.9984 |
| af3Draft_plus_primary_binary_strictStructural | train | af3Draft_trusted_pair_rows | 5665 | 0.9999 | 0.9998 | 0.9947 | 0.9991 | 0.9970 | 0.9984 |
| af3Draft_plus_primary_binary_strictStructural | validation | all_rows | 6384 | 0.9999 | 0.9996 | 0.9955 | 0.9970 | 0.9888 | 0.9988 |
| af3Draft_plus_primary_binary_strictStructural | validation | af3Draft_trusted_pair_rows | 999 | 1.0000 | 0.9999 | 0.9957 | 0.9987 | 0.9957 | 0.9987 |
| af3Draft_plus_primary_binary_strictStructural | external_cardioboost | all_rows | 185 | 0.9548 | 0.9828 | 0.9545 | 0.8113 | 0.9265 | 0.8776 |
| af3Draft_plus_primary_binary_strictStructural | external_cardioboost | af3Draft_trusted_pair_rows | 35 | 0.9560 | 0.9807 | 1.0000 | 0.7000 | 0.8929 | 1.0000 |
| af3Draft_plus_primary_binary_strictStructural | external_emerge | all_rows | 176 | 0.9885 | 0.9868 | 0.9792 | 0.9375 | 0.9495 | 0.9740 |
| af3Draft_plus_primary_binary_strictStructural | external_emerge | af3Draft_trusted_pair_rows | 42 | 0.9886 | 0.9922 | 0.9565 | 0.9474 | 0.9565 | 0.9474 |
| af3Draft_plus_primary_binary_strictStructural | external_hiro | all_rows | 69 | 0.9903 | 0.9859 | 0.9630 | 0.8810 | 0.8387 | 0.9737 |
| af3Draft_plus_primary_binary_strictStructural | external_hiro | af3Draft_trusted_pair_rows | 21 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |

## AF3Draft Features Used

### af3Draft_plus_primary_binary_fullDraft

```text
af3Draft_status
af3Draft_n_long_rows
af3Draft_n_jobs
af3Draft_n_complexes
af3Draft_n_chain_gene_contexts
af3Draft_has_residue_feature
af3Draft_best_pair_support_rank
af3Draft_best_pair_support_class
af3Draft_has_strong_or_moderate_pair
af3Draft_has_strong_pair
af3Draft_best_residue_plddt
af3Draft_mean_residue_plddt
af3Draft_min_distance_to_partner_angstrom
af3Draft_within_partner_interface_5A
af3Draft_within_partner_interface_8A
af3Draft_best_pair_iptm
af3Draft_best_pair_min_pae
af3Draft_best_pair_max_contact_probability
af3Draft_nearest_partner_genes
```
### af3Draft_plus_primary_binary_strictStructural

```text
af3Draft_has_residue_feature
af3Draft_best_pair_support_rank
af3Draft_has_strong_or_moderate_pair
af3Draft_has_strong_pair
af3Draft_best_residue_plddt
af3Draft_mean_residue_plddt
af3Draft_min_distance_to_partner_angstrom
af3Draft_within_partner_interface_5A
af3Draft_within_partner_interface_8A
af3Draft_best_pair_iptm
af3Draft_best_pair_min_pae
af3Draft_best_pair_max_contact_probability
```

## Interpretation

The most important comparison is not only global performance, because current AF3Draft coverage is small. The key smoke-test readout is whether AF3Draft-covered rows improve or stay stable when these features are added. If global metrics are unchanged, that is expected at this stage because most modelable variants have no AF3Draft output yet.
