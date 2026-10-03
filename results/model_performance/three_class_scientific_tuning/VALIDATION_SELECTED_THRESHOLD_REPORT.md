# Scientifically Controlled 3-Class Threshold Selection

Thresholds are selected using the internal validation split only. External HiRO, eMERGE, CardioBoost, and ClinVar summaries below are held-out evaluations of those validation-selected rules.

Decision rule: call `Pathogenic` if `prob_Pathogenic >= threshold`; otherwise choose `Benign` vs `VUS` by the larger probability.

## max_validation_macro_f1

| Dataset | path_threshold | rows | accuracy | macro_f1 | weighted_kappa | plp_sensitivity | plp_specificity | plp_ppv | plp_npv | binary_mcc_plp | pred_pathogenic | plp_to_vus | plp_to_benign | benign_to_plp | vus_to_plp |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| validation_selection_set | 0.7000 | 12361 | 0.9401 | 0.9181 | 0.9283 | 0.8325 | 0.9850 | 0.8704 | 0.9799 | 0.8337 | 1273 | 220 | 3 | 0 | 165 |
| external_hiro_source_records | 0.7000 | 325 | 0.7538 | 0.7463 | 0.6993 | 0.6000 | 0.9931 | 0.9130 | 0.9536 | 0.7170 | 23 | 14 | 0 | 0 | 2 |
| external_emerge_source_records | 0.7000 | 1974 | 0.8946 | 0.6520 | 0.4957 | 0.4796 | 0.9867 | 0.6528 | 0.9732 | 0.5402 | 72 | 50 | 1 | 0 | 25 |
| external_cardioboost_source_records | 0.7000 | 195 | 0.4821 | 0.4300 | 0.5676 | 0.4638 | 0.9825 | 0.9846 | 0.4308 | 0.4305 | 65 | 70 | 4 | 1 | 0 |
| external_hiro_variant_rows | 0.7000 | 179 | 0.8324 | 0.7870 | 0.7577 | 0.5556 | 0.9868 | 0.8824 | 0.9259 | 0.6621 | 17 | 12 | 0 | 0 | 2 |
| external_emerge_variant_rows | 0.7000 | 1974 | 0.8946 | 0.6520 | 0.4957 | 0.4796 | 0.9867 | 0.6528 | 0.9732 | 0.5402 | 72 | 50 | 1 | 0 | 25 |
| external_cardioboost_variant_rows | 0.7000 | 194 | 0.4794 | 0.4287 | 0.5660 | 0.4599 | 0.9825 | 0.9844 | 0.4308 | 0.4285 | 64 | 70 | 4 | 1 | 0 |
| clinvar_variant_rows | 0.7000 | 83688 | 0.9454 | 0.9295 | 0.9331 | 0.8700 | 0.9872 | 0.8916 | 0.9843 | 0.8665 | 8819 | 1129 | 46 | 11 | 945 |

## max_validation_weighted_kappa

| Dataset | path_threshold | rows | accuracy | macro_f1 | weighted_kappa | plp_sensitivity | plp_specificity | plp_ppv | plp_npv | binary_mcc_plp | pred_pathogenic | plp_to_vus | plp_to_benign | benign_to_plp | vus_to_plp |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| validation_selection_set | 0.6750 | 12361 | 0.9398 | 0.9181 | 0.9283 | 0.8415 | 0.9837 | 0.8615 | 0.9809 | 0.8338 | 1300 | 208 | 3 | 0 | 180 |
| external_hiro_source_records | 0.6750 | 325 | 0.7538 | 0.7463 | 0.6993 | 0.6000 | 0.9931 | 0.9130 | 0.9536 | 0.7170 | 23 | 14 | 0 | 0 | 2 |
| external_emerge_source_records | 0.6750 | 1974 | 0.8941 | 0.6508 | 0.4945 | 0.4796 | 0.9861 | 0.6438 | 0.9732 | 0.5361 | 73 | 50 | 1 | 0 | 26 |
| external_cardioboost_source_records | 0.6750 | 195 | 0.4974 | 0.4366 | 0.5795 | 0.4855 | 0.9825 | 0.9853 | 0.4409 | 0.4466 | 68 | 67 | 4 | 1 | 0 |
| external_hiro_variant_rows | 0.6750 | 179 | 0.8324 | 0.7870 | 0.7577 | 0.5556 | 0.9868 | 0.8824 | 0.9259 | 0.6621 | 17 | 12 | 0 | 0 | 2 |
| external_emerge_variant_rows | 0.6750 | 1974 | 0.8941 | 0.6508 | 0.4945 | 0.4796 | 0.9861 | 0.6438 | 0.9732 | 0.5361 | 73 | 50 | 1 | 0 | 26 |
| external_cardioboost_variant_rows | 0.6750 | 194 | 0.4948 | 0.4355 | 0.5780 | 0.4818 | 0.9825 | 0.9851 | 0.4409 | 0.4447 | 67 | 67 | 4 | 1 | 0 |
| clinvar_variant_rows | 0.6750 | 83688 | 0.9452 | 0.9295 | 0.9331 | 0.8777 | 0.9861 | 0.8840 | 0.9852 | 0.8665 | 8974 | 1059 | 46 | 11 | 1030 |

## max_validation_balanced_accuracy

| Dataset | path_threshold | rows | accuracy | macro_f1 | weighted_kappa | plp_sensitivity | plp_specificity | plp_ppv | plp_npv | binary_mcc_plp | pred_pathogenic | plp_to_vus | plp_to_benign | benign_to_plp | vus_to_plp |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| validation_selection_set | 0.4000 | 12361 | 0.9278 | 0.9019 | 0.9178 | 0.9083 | 0.9620 | 0.7426 | 0.9886 | 0.7978 | 1628 | 119 | 3 | 2 | 417 |
| external_hiro_source_records | 0.4000 | 325 | 0.7692 | 0.7842 | 0.7274 | 0.7429 | 0.9931 | 0.9286 | 0.9697 | 0.8131 | 28 | 9 | 0 | 0 | 2 |
| external_emerge_source_records | 0.4000 | 1974 | 0.8820 | 0.6473 | 0.4985 | 0.6224 | 0.9659 | 0.4880 | 0.9800 | 0.5247 | 125 | 36 | 1 | 0 | 64 |
| external_cardioboost_source_records | 0.4000 | 195 | 0.6308 | 0.4847 | 0.6600 | 0.6739 | 0.9474 | 0.9688 | 0.5455 | 0.5652 | 96 | 41 | 4 | 3 | 0 |
| external_hiro_variant_rows | 0.4000 | 179 | 0.8547 | 0.8287 | 0.7976 | 0.7037 | 0.9868 | 0.9048 | 0.9494 | 0.7680 | 21 | 8 | 0 | 0 | 2 |
| external_emerge_variant_rows | 0.4000 | 1974 | 0.8820 | 0.6473 | 0.4985 | 0.6224 | 0.9659 | 0.4880 | 0.9800 | 0.5247 | 125 | 36 | 1 | 0 | 64 |
| external_cardioboost_variant_rows | 0.4000 | 194 | 0.6289 | 0.4841 | 0.6590 | 0.6715 | 0.9474 | 0.9684 | 0.5455 | 0.5639 | 95 | 41 | 4 | 3 | 0 |
| clinvar_variant_rows | 0.4000 | 83688 | 0.9352 | 0.9155 | 0.9245 | 0.9449 | 0.9666 | 0.7741 | 0.9931 | 0.8363 | 11032 | 452 | 46 | 16 | 2476 |

## max_validation_plp_mcc

| Dataset | path_threshold | rows | accuracy | macro_f1 | weighted_kappa | plp_sensitivity | plp_specificity | plp_ppv | plp_npv | binary_mcc_plp | pred_pathogenic | plp_to_vus | plp_to_benign | benign_to_plp | vus_to_plp |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| validation_selection_set | 0.6750 | 12361 | 0.9398 | 0.9181 | 0.9283 | 0.8415 | 0.9837 | 0.8615 | 0.9809 | 0.8338 | 1300 | 208 | 3 | 0 | 180 |
| external_hiro_source_records | 0.6750 | 325 | 0.7538 | 0.7463 | 0.6993 | 0.6000 | 0.9931 | 0.9130 | 0.9536 | 0.7170 | 23 | 14 | 0 | 0 | 2 |
| external_emerge_source_records | 0.6750 | 1974 | 0.8941 | 0.6508 | 0.4945 | 0.4796 | 0.9861 | 0.6438 | 0.9732 | 0.5361 | 73 | 50 | 1 | 0 | 26 |
| external_cardioboost_source_records | 0.6750 | 195 | 0.4974 | 0.4366 | 0.5795 | 0.4855 | 0.9825 | 0.9853 | 0.4409 | 0.4466 | 68 | 67 | 4 | 1 | 0 |
| external_hiro_variant_rows | 0.6750 | 179 | 0.8324 | 0.7870 | 0.7577 | 0.5556 | 0.9868 | 0.8824 | 0.9259 | 0.6621 | 17 | 12 | 0 | 0 | 2 |
| external_emerge_variant_rows | 0.6750 | 1974 | 0.8941 | 0.6508 | 0.4945 | 0.4796 | 0.9861 | 0.6438 | 0.9732 | 0.5361 | 73 | 50 | 1 | 0 | 26 |
| external_cardioboost_variant_rows | 0.6750 | 194 | 0.4948 | 0.4355 | 0.5780 | 0.4818 | 0.9825 | 0.9851 | 0.4409 | 0.4447 | 67 | 67 | 4 | 1 | 0 |
| clinvar_variant_rows | 0.6750 | 83688 | 0.9452 | 0.9295 | 0.9331 | 0.8777 | 0.9861 | 0.8840 | 0.9852 | 0.8665 | 8974 | 1059 | 46 | 11 | 1030 |

## validation_constrained_sens90_spec90_ppv60_max_macro_f1

| Dataset | path_threshold | rows | accuracy | macro_f1 | weighted_kappa | plp_sensitivity | plp_specificity | plp_ppv | plp_npv | binary_mcc_plp | pred_pathogenic | plp_to_vus | plp_to_benign | benign_to_plp | vus_to_plp |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| validation_selection_set | 0.4250 | 12361 | 0.9291 | 0.9035 | 0.9188 | 0.9016 | 0.9643 | 0.7528 | 0.9878 | 0.8008 | 1594 | 128 | 3 | 2 | 392 |
| external_hiro_source_records | 0.4250 | 325 | 0.7692 | 0.7842 | 0.7274 | 0.7429 | 0.9931 | 0.9286 | 0.9697 | 0.8131 | 28 | 9 | 0 | 0 | 2 |
| external_emerge_source_records | 0.4250 | 1974 | 0.8870 | 0.6555 | 0.5074 | 0.6122 | 0.9717 | 0.5310 | 0.9796 | 0.5460 | 113 | 37 | 1 | 0 | 53 |
| external_cardioboost_source_records | 0.4250 | 195 | 0.6154 | 0.4795 | 0.6483 | 0.6522 | 0.9474 | 0.9677 | 0.5294 | 0.5460 | 93 | 44 | 4 | 3 | 0 |
| external_hiro_variant_rows | 0.4250 | 179 | 0.8547 | 0.8287 | 0.7976 | 0.7037 | 0.9868 | 0.9048 | 0.9494 | 0.7680 | 21 | 8 | 0 | 0 | 2 |
| external_emerge_variant_rows | 0.4250 | 1974 | 0.8870 | 0.6555 | 0.5074 | 0.6122 | 0.9717 | 0.5310 | 0.9796 | 0.5460 | 113 | 37 | 1 | 0 | 53 |
| external_cardioboost_variant_rows | 0.4250 | 194 | 0.6134 | 0.4789 | 0.6473 | 0.6496 | 0.9474 | 0.9674 | 0.5294 | 0.5446 | 92 | 44 | 4 | 3 | 0 |
| clinvar_variant_rows | 0.4250 | 83688 | 0.9370 | 0.9183 | 0.9263 | 0.9410 | 0.9691 | 0.7869 | 0.9927 | 0.8424 | 10808 | 487 | 46 | 16 | 2287 |
