# Systematic AF3 V2 Gene 001-015 CatBoost Test

This is an additive V2 experiment using systematic AlphaFold 3 outputs from genes 1-15 only. The established baseline tables and models were not modified.

## AF3 Interaction Threshold Policy

- Strong interface: pair ipTM >= 0.70, pair/min PAE <= 5, max contact probability >= 0.60.
- Exploratory supported interface: pair ipTM >= 0.60, pair/min PAE <= 8, max contact probability >= 0.45.
- The main model uses structural columns only. AF3 job IDs, partner names, status text, and count/provenance columns are excluded from training.
- Unsupported jobs are retained in QC files but are not treated as positive interaction evidence.

## Coverage

| scope | rows | genes | rows_with_any_systematicAF3_output | rows_with_exploratory_or_strong_pair | rows_with_strong_pair | pct_any_output | pct_exploratory_or_strong_pair | pct_strong_pair |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all_modeling_rows | 85677 | 61 | 40764 | 22899 | 14136 | 47.58% | 26.73% | 16.50% |
| primary_binary_rows | 42990 | 53 | 16595 | 9109 | 5474 | 38.60% | 21.19% | 12.73% |
| vus_rows | 42361 | 61 | 23977 | 13650 | 8574 | 56.60% | 32.22% | 20.24% |

## Performance

| model | subset | scope | rows | auroc | auprc | sensitivity | specificity | ppv | npv |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| baseline_primary_binary_existing | external_cardioboost | all_rows | 185 | 0.9520 | 0.9812 | 0.9545 | 0.8113 | 0.9265 | 0.8776 |
| baseline_primary_binary_existing | external_cardioboost | systematicAF3_any_output_rows | 127 | 0.9392 | 0.9819 | 0.9381 | 0.8333 | 0.9479 | 0.8065 |
| baseline_primary_binary_existing | external_cardioboost | systematicAF3_trusted_pair_rows | 85 | 0.9239 | 0.9826 | 0.9275 | 0.8125 | 0.9552 | 0.7222 |
| baseline_primary_binary_existing | external_cardioboost | systematicAF3_strong_pair_rows | 51 | 0.9415 | 0.9846 | 0.9756 | 0.7000 | 0.9302 | 0.8750 |
| baseline_primary_binary_existing | external_emerge | all_rows | 176 | 0.9910 | 0.9914 | 0.9792 | 0.9375 | 0.9495 | 0.9740 |
| baseline_primary_binary_existing | external_emerge | systematicAF3_any_output_rows | 154 | 0.9976 | 0.9982 | 0.9880 | 0.9577 | 0.9647 | 0.9855 |
| baseline_primary_binary_existing | external_emerge | systematicAF3_trusted_pair_rows | 115 | 0.9976 | 0.9980 | 0.9831 | 0.9643 | 0.9667 | 0.9818 |
| baseline_primary_binary_existing | external_emerge | systematicAF3_strong_pair_rows | 100 | 0.9971 | 0.9981 | 0.9828 | 0.9762 | 0.9828 | 0.9762 |
| baseline_primary_binary_existing | external_hiro | all_rows | 69 | 0.9877 | 0.9830 | 0.9630 | 0.8810 | 0.8387 | 0.9737 |
| baseline_primary_binary_existing | external_hiro | systematicAF3_any_output_rows | 49 | 0.9892 | 0.9822 | 1.0000 | 0.9032 | 0.8571 | 1.0000 |
| baseline_primary_binary_existing | external_hiro | systematicAF3_trusted_pair_rows | 25 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| baseline_primary_binary_existing | external_hiro | systematicAF3_strong_pair_rows | 12 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| baseline_primary_binary_existing | train | all_rows | 36176 | 0.9998 | 0.9992 | 0.9932 | 0.9978 | 0.9918 | 0.9982 |
| baseline_primary_binary_existing | train | systematicAF3_any_output_rows | 13899 | 0.9998 | 0.9993 | 0.9965 | 0.9979 | 0.9943 | 0.9987 |
| baseline_primary_binary_existing | train | systematicAF3_trusted_pair_rows | 7607 | 0.9999 | 0.9998 | 0.9977 | 0.9976 | 0.9920 | 0.9993 |
| baseline_primary_binary_existing | train | systematicAF3_strong_pair_rows | 4541 | 1.0000 | 0.9999 | 0.9962 | 0.9980 | 0.9933 | 0.9989 |
| baseline_primary_binary_existing | validation | all_rows | 6384 | 0.9999 | 0.9996 | 0.9947 | 0.9970 | 0.9888 | 0.9986 |
| baseline_primary_binary_existing | validation | systematicAF3_any_output_rows | 2366 | 0.9999 | 0.9998 | 0.9953 | 0.9960 | 0.9890 | 0.9983 |
| baseline_primary_binary_existing | validation | systematicAF3_trusted_pair_rows | 1277 | 0.9999 | 0.9998 | 0.9900 | 0.9959 | 0.9868 | 0.9969 |
| baseline_primary_binary_existing | validation | systematicAF3_strong_pair_rows | 770 | 0.9998 | 0.9995 | 0.9839 | 0.9932 | 0.9786 | 0.9949 |
| systematicAF3_v2_primary_binary_structuralOnly | train | all_rows | 36176 | 0.9998 | 0.9993 | 0.9931 | 0.9978 | 0.9918 | 0.9982 |
| systematicAF3_v2_primary_binary_structuralOnly | train | systematicAF3_any_output_rows | 13899 | 0.9998 | 0.9993 | 0.9967 | 0.9980 | 0.9946 | 0.9988 |
| systematicAF3_v2_primary_binary_structuralOnly | train | systematicAF3_trusted_pair_rows | 7607 | 0.9999 | 0.9998 | 0.9977 | 0.9978 | 0.9926 | 0.9993 |
| systematicAF3_v2_primary_binary_structuralOnly | train | systematicAF3_strong_pair_rows | 4541 | 1.0000 | 0.9999 | 0.9962 | 0.9983 | 0.9943 | 0.9989 |
| systematicAF3_v2_primary_binary_structuralOnly | validation | all_rows | 6384 | 0.9999 | 0.9996 | 0.9955 | 0.9968 | 0.9881 | 0.9988 |
| systematicAF3_v2_primary_binary_structuralOnly | validation | systematicAF3_any_output_rows | 2366 | 0.9999 | 0.9998 | 0.9969 | 0.9954 | 0.9875 | 0.9988 |
| systematicAF3_v2_primary_binary_structuralOnly | validation | systematicAF3_trusted_pair_rows | 1277 | 0.9999 | 0.9998 | 0.9934 | 0.9949 | 0.9836 | 0.9979 |
| systematicAF3_v2_primary_binary_structuralOnly | validation | systematicAF3_strong_pair_rows | 770 | 0.9998 | 0.9995 | 0.9892 | 0.9914 | 0.9735 | 0.9966 |
| systematicAF3_v2_primary_binary_structuralOnly | external_cardioboost | all_rows | 185 | 0.9554 | 0.9835 | 0.9621 | 0.7925 | 0.9203 | 0.8936 |
| systematicAF3_v2_primary_binary_structuralOnly | external_cardioboost | systematicAF3_any_output_rows | 127 | 0.9447 | 0.9844 | 0.9485 | 0.8333 | 0.9485 | 0.8333 |
| systematicAF3_v2_primary_binary_structuralOnly | external_cardioboost | systematicAF3_trusted_pair_rows | 85 | 0.9266 | 0.9839 | 0.9420 | 0.8125 | 0.9559 | 0.7647 |
| systematicAF3_v2_primary_binary_structuralOnly | external_cardioboost | systematicAF3_strong_pair_rows | 51 | 0.9512 | 0.9873 | 1.0000 | 0.7000 | 0.9318 | 1.0000 |
| systematicAF3_v2_primary_binary_structuralOnly | external_emerge | all_rows | 176 | 0.9900 | 0.9896 | 0.9792 | 0.9375 | 0.9495 | 0.9740 |
| systematicAF3_v2_primary_binary_structuralOnly | external_emerge | systematicAF3_any_output_rows | 154 | 0.9973 | 0.9979 | 0.9880 | 0.9577 | 0.9647 | 0.9855 |
| systematicAF3_v2_primary_binary_structuralOnly | external_emerge | systematicAF3_trusted_pair_rows | 115 | 0.9976 | 0.9980 | 0.9831 | 0.9464 | 0.9508 | 0.9815 |
| systematicAF3_v2_primary_binary_structuralOnly | external_emerge | systematicAF3_strong_pair_rows | 100 | 0.9971 | 0.9981 | 0.9828 | 0.9286 | 0.9500 | 0.9750 |
| systematicAF3_v2_primary_binary_structuralOnly | external_hiro | all_rows | 69 | 0.9868 | 0.9809 | 0.9630 | 0.8810 | 0.8387 | 0.9737 |
| systematicAF3_v2_primary_binary_structuralOnly | external_hiro | systematicAF3_any_output_rows | 49 | 0.9875 | 0.9785 | 1.0000 | 0.9032 | 0.8571 | 1.0000 |
| systematicAF3_v2_primary_binary_structuralOnly | external_hiro | systematicAF3_trusted_pair_rows | 25 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| systematicAF3_v2_primary_binary_structuralOnly | external_hiro | systematicAF3_strong_pair_rows | 12 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | train | all_rows | 36176 | 0.9998 | 0.9993 | 0.9931 | 0.9980 | 0.9925 | 0.9982 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | train | systematicAF3_any_output_rows | 13899 | 0.9999 | 0.9994 | 0.9967 | 0.9979 | 0.9943 | 0.9988 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | train | systematicAF3_trusted_pair_rows | 7607 | 1.0000 | 0.9998 | 0.9977 | 0.9978 | 0.9926 | 0.9993 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | train | systematicAF3_strong_pair_rows | 4541 | 1.0000 | 0.9999 | 0.9962 | 0.9983 | 0.9943 | 0.9989 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | validation | all_rows | 6384 | 0.9999 | 0.9996 | 0.9947 | 0.9970 | 0.9888 | 0.9986 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | validation | systematicAF3_any_output_rows | 2366 | 0.9999 | 0.9998 | 0.9953 | 0.9960 | 0.9890 | 0.9983 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | validation | systematicAF3_trusted_pair_rows | 1277 | 0.9999 | 0.9997 | 0.9900 | 0.9959 | 0.9868 | 0.9969 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | validation | systematicAF3_strong_pair_rows | 770 | 0.9998 | 0.9995 | 0.9839 | 0.9932 | 0.9786 | 0.9949 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | external_cardioboost | all_rows | 185 | 0.9504 | 0.9809 | 0.9545 | 0.8113 | 0.9265 | 0.8776 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | external_cardioboost | systematicAF3_any_output_rows | 127 | 0.9381 | 0.9819 | 0.9381 | 0.8333 | 0.9479 | 0.8065 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | external_cardioboost | systematicAF3_trusted_pair_rows | 85 | 0.9221 | 0.9827 | 0.9275 | 0.8125 | 0.9552 | 0.7222 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | external_cardioboost | systematicAF3_strong_pair_rows | 51 | 0.9463 | 0.9859 | 0.9756 | 0.7000 | 0.9302 | 0.8750 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | external_emerge | all_rows | 176 | 0.9892 | 0.9890 | 0.9688 | 0.9500 | 0.9588 | 0.9620 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | external_emerge | systematicAF3_any_output_rows | 154 | 0.9966 | 0.9975 | 0.9759 | 0.9718 | 0.9759 | 0.9718 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | external_emerge | systematicAF3_trusted_pair_rows | 115 | 0.9970 | 0.9975 | 0.9831 | 0.9643 | 0.9667 | 0.9818 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | external_emerge | systematicAF3_strong_pair_rows | 100 | 0.9963 | 0.9977 | 0.9828 | 0.9762 | 0.9828 | 0.9762 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | external_hiro | all_rows | 69 | 0.9894 | 0.9846 | 0.9630 | 0.8810 | 0.8387 | 0.9737 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | external_hiro | systematicAF3_any_output_rows | 49 | 0.9892 | 0.9819 | 1.0000 | 0.9032 | 0.8571 | 1.0000 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | external_hiro | systematicAF3_trusted_pair_rows | 25 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| systematicAF3_v2_primary_binary_structuralPlusFlags | external_hiro | systematicAF3_strong_pair_rows | 12 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |

## Pair QC Preview

| target_gene | job_id | gene_a | gene_b | pair_support_class | pair_iptm | pair_pae_min_summary | max_contact_probability |
| --- | --- | --- | --- | --- | --- | --- | --- |
| FLNC | AF3_SYS_G002_FLNC_TAIL_2593_2725_HOMODIMER_DRAFT | FLNC | FLNC | exploratory_supported_interface | 0.680 | 0.930 | 0.970 |
| DSC2 | AF3_SYS_G015_DSC2_CYTO_716_901_JUP_FULL_LENGTH_DRAFT | DSC2 | JUP | exploratory_supported_interface | 0.680 | 1.150 | 0.990 |
| MYH7 | AF3_SYS_G004_MYH7_NECK_760_840_MYL2_MYL3_DRAFT | MYL2 | MYL3 | exploratory_supported_interface | 0.660 | 2.270 | 0.460 |
| CACNA1C | AF3_SYS_G007_CACNA1C_CACNB2_CALM1_CA4_FULL_LENGTH_DRAFT | CACNA1C | CALM1 | exploratory_supported_interface | 0.630 | 6.510 | 0.700 |
| PKP2 | AF3_SYS_G011_PKP2_DSP_NTERM_1_584_DRAFT | PKP2 | DSP | exploratory_supported_interface | 0.630 | 2.230 | 0.920 |
| MYH7 | AF3_SYS_G004_MYH7_MOTOR_85_778_ACTC1_DRAFT | MYH7 | ACTC1 | exploratory_supported_interface | 0.600 | 7.110 | 0.820 |
| DSP | AF3_SYS_G003_DSP_NTERM_1_584_PKP2_JUP_DRAFT | DSP | PKP2 | moderate_possible_interface | 0.590 | 2.560 | 0.930 |
| DSP | AF3_SYS_G003_DSP_NTERM_1_584_PKP2_JUP_DRAFT | DSP | JUP | moderate_possible_interface | 0.590 | 3.330 | 0.930 |
| KCNQ1 | AF3_SYS_G009_KCNQ1_KCNE1_CALM1_CA16_TETRAMER_DRAFT | KCNQ1 | KCNQ1 | moderate_possible_interface | 0.590 | 1.440 | 1.000 |
| KCNQ1 | AF3_SYS_G009_KCNQ1_KCNE1_CALM1_CA16_TETRAMER_DRAFT | KCNQ1 | KCNQ1 | moderate_possible_interface | 0.580 | 1.480 | 1.000 |
| KCNQ1 | AF3_SYS_G009_KCNQ1_KCNE1_CALM1_CA16_TETRAMER_DRAFT | KCNQ1 | KCNQ1 | moderate_possible_interface | 0.580 | 1.430 | 1.000 |
| KCNQ1 | AF3_SYS_G009_KCNQ1_KCNE1_CALM1_CA16_TETRAMER_DRAFT | KCNQ1 | KCNQ1 | moderate_possible_interface | 0.580 | 1.450 | 1.000 |
| MYBPC3 | AF3_SYS_G006_MYBPC3_NTERM_1_452_ACTC1_DRAFT | MYBPC3 | ACTC1 | moderate_possible_interface | 0.570 | 4.120 | 0.560 |
| KCNQ1 | AF3_SYS_G009_KCNQ1_KCNE1_CALM1_CA16_TETRAMER_DRAFT | KCNQ1 | KCNQ1 | moderate_possible_interface | 0.570 | 1.870 | 0.930 |
| KCNQ1 | AF3_SYS_G009_KCNQ1_KCNE1_CALM1_CA16_TETRAMER_DRAFT | KCNQ1 | KCNQ1 | moderate_possible_interface | 0.570 | 1.890 | 0.930 |
| PKP2 | AF3_SYS_G011_PKP2_DSP_NTERM_JUP_DRAFT | PKP2 | DSP | moderate_possible_interface | 0.570 | 3.100 | 0.910 |
| DSP | AF3_SYS_G003_DSP_NTERM_1_584_PKP2_JUP_DRAFT | PKP2 | JUP | moderate_possible_interface | 0.560 | 3.720 | 0.860 |
| KCNQ1 | AF3_SYS_G009_KCNQ1_KCNE1_CALM1_CA16_TETRAMER_DRAFT | CALM1 | CALM1 | moderate_possible_interface | 0.560 | 4.710 | 0.930 |
| KCNQ1 | AF3_SYS_G009_KCNQ1_KCNE1_CALM1_CA16_TETRAMER_DRAFT | CALM1 | CALM1 | moderate_possible_interface | 0.550 | 4.770 | 0.920 |
| KCNQ1 | AF3_SYS_G009_KCNQ1_KCNE1_CALM1_CA16_TETRAMER_DRAFT | CALM1 | CALM1 | moderate_possible_interface | 0.550 | 5.120 | 0.920 |
| KCNQ1 | AF3_SYS_G009_KCNQ1_KCNE1_CALM1_CA16_TETRAMER_DRAFT | CALM1 | CALM1 | moderate_possible_interface | 0.550 | 4.870 | 0.920 |
| PKP2 | AF3_SYS_G011_PKP2_DSP_NTERM_JUP_DRAFT | DSP | JUP | moderate_possible_interface | 0.550 | 4.380 | 0.870 |
| KCNH2 | AF3_SYS_G008_KCNH2_HOMOTETRAMER_FULL_LENGTH_DRAFT | KCNH2 | KCNH2 | moderate_possible_interface | 0.540 | 1.760 | 0.970 |
| KCNH2 | AF3_SYS_G008_KCNH2_HOMOTETRAMER_FULL_LENGTH_DRAFT | KCNH2 | KCNH2 | moderate_possible_interface | 0.540 | 1.790 | 0.970 |
| KCNH2 | AF3_SYS_G008_KCNH2_HOMOTETRAMER_FULL_LENGTH_DRAFT | KCNH2 | KCNH2 | moderate_possible_interface | 0.540 | 1.800 | 0.970 |
| KCNH2 | AF3_SYS_G008_KCNH2_HOMOTETRAMER_FULL_LENGTH_DRAFT | KCNH2 | KCNH2 | moderate_possible_interface | 0.530 | 1.860 | 0.970 |
| CACNA1C | AF3_SYS_G007_CACNA1C_CACNB2_CALM1_CA4_FULL_LENGTH_DRAFT | CACNA1C | CACNB2 | moderate_possible_interface | 0.520 | 1.780 | 0.960 |
| CACNA1C | AF3_SYS_G007_CACNA1C_CACNB2_FULL_LENGTH_DRAFT | CACNA1C | CACNB2 | moderate_possible_interface | 0.520 | 1.940 | 0.960 |
| RYR2 | AF3_SYS_G001_RYR2_CaMBD1_1940_1965_CALM1_CA4_DRAFT | RYR2 | RYR2 | moderate_possible_interface | 0.500 | 5.360 | 0.780 |
| KCNH2 | AF3_SYS_G008_KCNH2_HOMOTETRAMER_FULL_LENGTH_DRAFT | KCNH2 | KCNH2 | moderate_possible_interface | 0.500 | 2.720 | 0.710 |

## AF3 Features Used

### systematicAF3_v2_primary_binary_structuralOnly

```text
systematicAF3_best_residue_plddt
systematicAF3_mean_residue_plddt
systematicAF3_min_distance_to_partner_angstrom
systematicAF3_within_partner_interface_5A
systematicAF3_within_partner_interface_8A
systematicAF3_best_pair_iptm
systematicAF3_best_pair_min_pae
systematicAF3_best_pair_max_contact_probability
```

### systematicAF3_v2_primary_binary_structuralPlusFlags

```text
systematicAF3_best_pair_support_rank
systematicAF3_has_exploratory_or_strong_pair
systematicAF3_has_strong_pair
systematicAF3_best_residue_plddt
systematicAF3_mean_residue_plddt
systematicAF3_min_distance_to_partner_angstrom
systematicAF3_within_partner_interface_5A
systematicAF3_within_partner_interface_8A
systematicAF3_best_pair_iptm
systematicAF3_best_pair_min_pae
systematicAF3_best_pair_max_contact_probability
```

## Outputs

- Aggregated AF3 variant features: `datasets/feature_sources/alphafold3/processed/systematic_v2_gene001_015/systematicAF3_v2_variant_features_aggregated.tsv`
- Long AF3 variant features: `datasets/feature_sources/alphafold3/processed/systematic_v2_gene001_015/systematicAF3_v2_variant_features_long.tsv`
- Merged modeling table: `datasets/modeling/interim/modeling_table_with_splits_weights_plus_systematicAF3_v2_gene001_015.tsv`
- Metrics: `results/model_performance/systematic_AF3_v2_gene001_015/systematic_AF3_v2_metrics.tsv`
