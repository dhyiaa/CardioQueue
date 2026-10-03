# Primary Binary CatBoost On Binary-Labeled Rows Only

Truth labels are restricted to `Benign` and `Pathogenic`; true VUS rows are excluded from this report.

Prediction rule remains:

- `P(pathogenic) <= 0.1`: Benign-like
- `0.1 < P(pathogenic) < 0.9`: VUS/defer
- `P(pathogenic) >= 0.9`: Pathogenic-like

Metrics are shown two ways: first counting deferrals as non-pathogenic for binary P/LP-vs-rest metrics, and second on only called rows after excluding deferrals.

## Summary Metrics

| Dataset | binary_rows_with_prediction | true_benign | true_pathogenic | pred_benign | pred_vus_defer | pred_pathogenic | deferral_rate | sensitivity_pathogenic_counting_defers_as_negative | specificity_counting_defers_as_non_pathogenic | ppv_pathogenic | npv_counting_defers_as_non_pathogenic | mcc_counting_defers_as_non_pathogenic | called_rows_excluding_defers | called_accuracy_excluding_defers | called_sensitivity_excluding_defers | called_specificity_excluding_defers | called_ppv_excluding_defers | called_npv_excluding_defers | called_mcc_excluding_defers | pathogenic_to_vus_deferrals | pathogenic_to_benign_errors | benign_to_pathogenic_escalations |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| All variant rows | 43071 | 33923 | 9148 | 33589 | 464 | 9018 | 0.0108 | 0.9811 | 0.9987 | 0.9952 | 0.9949 | 0.9850 | 42607 | 0.9985 | 0.9976 | 0.9987 | 0.9952 | 0.9993 | 0.9954 | 151 | 22 | 43 |
| ClinVar variant rows | 42834 | 33796 | 9038 | 33495 | 426 | 8913 | 0.0099 | 0.9822 | 0.9989 | 0.9960 | 0.9953 | 0.9862 | 42408 | 0.9987 | 0.9979 | 0.9989 | 0.9960 | 0.9994 | 0.9961 | 142 | 19 | 36 |
| HiRO variant rows | 74 | 43 | 31 | 36 | 7 | 31 | 0.0946 | 0.9032 | 0.9302 | 0.9032 | 0.9302 | 0.8335 | 67 | 0.9552 | 1.0000 | 0.9231 | 0.9032 | 1.0000 | 0.9131 | 3 | 0 | 3 |
| eMERGE variant rows | 241 | 130 | 111 | 107 | 28 | 106 | 0.1162 | 0.9459 | 0.9923 | 0.9906 | 0.9556 | 0.9422 | 213 | 0.9953 | 1.0000 | 0.9907 | 0.9906 | 1.0000 | 0.9907 | 6 | 0 | 1 |
| CardioBoost variant rows | 194 | 57 | 137 | 34 | 26 | 134 | 0.1340 | 0.9343 | 0.8947 | 0.9552 | 0.8500 | 0.8170 | 168 | 0.9464 | 0.9771 | 0.8378 | 0.9552 | 0.9118 | 0.8406 | 6 | 3 | 6 |
| HiRO source records | 236 | 190 | 46 | 149 | 28 | 59 | 0.1186 | 0.8913 | 0.9053 | 0.6949 | 0.9718 | 0.7287 | 208 | 0.9135 | 1.0000 | 0.8922 | 0.6949 | 1.0000 | 0.7874 | 5 | 0 | 18 |
| eMERGE source records | 264 | 137 | 127 | 108 | 35 | 121 | 0.1326 | 0.9449 | 0.9927 | 0.9917 | 0.9510 | 0.9402 | 229 | 0.9956 | 1.0000 | 0.9908 | 0.9917 | 1.0000 | 0.9913 | 7 | 0 | 1 |
| CardioBoost source records | 355 | 156 | 199 | 60 | 92 | 203 | 0.2592 | 0.8693 | 0.8077 | 0.8522 | 0.8289 | 0.6791 | 263 | 0.8669 | 0.9719 | 0.6471 | 0.8522 | 0.9167 | 0.6899 | 21 | 5 | 30 |

## Confusion Matrix Files

- All variant rows: `results/model_performance/primary_binary_three_zone_binary_only/all_variant_rows_binary_only_confusion.tsv`
- ClinVar variant rows: `results/model_performance/primary_binary_three_zone_binary_only/clinvar_variant_rows_binary_only_confusion.tsv`
- HiRO variant rows: `results/model_performance/primary_binary_three_zone_binary_only/hiro_variant_rows_binary_only_confusion.tsv`
- eMERGE variant rows: `results/model_performance/primary_binary_three_zone_binary_only/emerge_variant_rows_binary_only_confusion.tsv`
- CardioBoost variant rows: `results/model_performance/primary_binary_three_zone_binary_only/cardioboost_variant_rows_binary_only_confusion.tsv`
- HiRO source records: `results/model_performance/primary_binary_three_zone_binary_only/hiro_source_records_binary_only_confusion.tsv`
- eMERGE source records: `results/model_performance/primary_binary_three_zone_binary_only/emerge_source_records_binary_only_confusion.tsv`
- CardioBoost source records: `results/model_performance/primary_binary_three_zone_binary_only/cardioboost_source_records_binary_only_confusion.tsv`
