# CardioQueue v1.0 Three-Zone Report

Model: `results/models/cardioqueue_v1_source_heldout/primary_binary_catboost_model.cbm`

Three-zone rule:

- `P(pathogenic) <= 0.1`: Benign-like
- `0.1 < P(pathogenic) < 0.9`: VUS/defer
- `P(pathogenic) >= 0.9`: Pathogenic-like

Important: this is a binary model with a deferral zone. The VUS/defer output means insufficient binary confidence, not a separately learned VUS biology class.

## Summary Metrics

| Dataset | rows_total | rows_with_prediction | rows_without_prediction | true_benign | true_vus | true_pathogenic | pred_benign | pred_vus_defer | pred_pathogenic | exact_3class_accuracy | weighted_kappa | macro_f1 | plp_sensitivity | specificity_non_plp | ppv_plp | npv_non_plp | binary_mcc_plp_vs_rest | vus_deferral_rate | plp_positive_calls | true_positives | false_positives | plp_to_vus_deferrals | plp_to_benign_errors | benign_to_plp_escalations | vus_to_plp_escalations |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| All variant rows | 85677 | 85677 | 0 | 33923 | 42361 | 9148 | 44161 | 14332 | 27184 | 0.6606 | 0.7206 | 0.6155 | 0.9820 | 0.7627 | 0.3316 | 0.9972 | 0.4948 | 0.1673 | 27184 | 8983 | 18103 | 147 | 18 | 39 | 18064 |
| ClinVar variant rows | 83915 | 83915 | 0 | 33796 | 40872 | 9038 | 43479 | 13710 | 26726 | 0.6652 | 0.7264 | 0.6178 | 0.9832 | 0.7623 | 0.3337 | 0.9973 | 0.4968 | 0.1634 | 26726 | 8886 | 17746 | 138 | 14 | 35 | 17711 |
| HiRO variant rows | 240 | 240 | 0 | 43 | 139 | 31 | 103 | 72 | 65 | 0.5634 | 0.5566 | 0.5723 | 0.9032 | 0.8242 | 0.4667 | 0.9804 | 0.5703 | 0.3000 | 65 | 28 | 32 | 3 | 0 | 2 | 30 |
| eMERGE variant rows | 2754 | 2754 | 0 | 130 | 2394 | 111 | 810 | 1097 | 847 | 0.4706 | 0.2317 | 0.3545 | 0.9459 | 0.7199 | 0.1293 | 0.9967 | 0.2897 | 0.3983 | 847 | 105 | 707 | 5 | 1 | 1 | 706 |
| CardioBoost variant rows | 353 | 353 | 0 | 57 | 0 | 137 | 70 | 86 | 197 | 0.8299 | 0.8113 | 0.5602 | 0.9197 | 0.9298 | 0.9692 | 0.8281 | 0.8230 | 0.2436 | 197 | 126 | 4 | 8 | 3 | 4 | 0 |
| HiRO source records | 482 | 408 | 74 | 190 | 172 | 46 | 219 | 94 | 95 | 0.6422 | 0.5675 | 0.6145 | 0.8913 | 0.8508 | 0.4316 | 0.9840 | 0.5554 | 0.2304 | 95 | 41 | 54 | 5 | 0 | 18 | 36 |
| eMERGE source records | 2754 | 2754 | 0 | 137 | 2490 | 127 | 810 | 1097 | 847 | 0.4695 | 0.2364 | 0.3566 | 0.9449 | 0.7233 | 0.1417 | 0.9963 | 0.3037 | 0.3983 | 847 | 120 | 727 | 6 | 1 | 1 | 726 |
| CardioBoost source records | 355 | 355 | 0 | 156 | 0 | 199 | 70 | 86 | 199 | 0.6648 | 0.6322 | 0.4782 | 0.8593 | 0.8205 | 0.8593 | 0.8205 | 0.6798 | 0.2423 | 199 | 171 | 28 | 23 | 5 | 28 | 0 |

## Confusion Matrix Files

- All variant rows: `results/model_performance/primary_binary_three_zone_all_sources/all_variant_rows_confusion.tsv`
- ClinVar variant rows: `results/model_performance/primary_binary_three_zone_all_sources/clinvar_variant_rows_confusion.tsv`
- HiRO variant rows: `results/model_performance/primary_binary_three_zone_all_sources/hiro_variant_rows_confusion.tsv`
- eMERGE variant rows: `results/model_performance/primary_binary_three_zone_all_sources/emerge_variant_rows_confusion.tsv`
- CardioBoost variant rows: `results/model_performance/primary_binary_three_zone_all_sources/cardioboost_variant_rows_confusion.tsv`
- HiRO source records: `results/model_performance/primary_binary_three_zone_all_sources/hiro_source_records_confusion.tsv`
- eMERGE source records: `results/model_performance/primary_binary_three_zone_all_sources/emerge_source_records_confusion.tsv`
- CardioBoost source records: `results/model_performance/primary_binary_three_zone_all_sources/cardioboost_source_records_confusion.tsv`
