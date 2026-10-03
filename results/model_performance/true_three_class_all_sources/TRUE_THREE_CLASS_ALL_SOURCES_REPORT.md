# True 3-Class CatBoost All-Sources Report

Model: `results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued/three_class_catboost_model.cbm`

Prediction classes are learned directly by the model:

- `Benign`
- `VUS`
- `Pathogenic`

Important: this is the exploratory true 3-class model, not the primary CardioBoost-comparable binary + deferral model.

## Summary Metrics

| Dataset | rows_total | rows_with_prediction | rows_without_prediction | true_benign | true_vus | true_pathogenic | pred_benign | pred_vus | pred_pathogenic | exact_3class_accuracy | weighted_kappa | macro_f1 | plp_sensitivity | specificity_non_plp | ppv_plp | npv_non_plp | binary_mcc_plp_vs_rest | vus_prediction_rate | plp_positive_calls | true_positives | false_positives | plp_to_vus_deferrals | plp_to_benign_errors | benign_to_plp_escalations | vus_to_plp_escalations |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| All variant rows | 84721 | 84721 | 0 | 33859 | 41731 | 9131 | 34439 | 40034 | 10248 | 0.9396 | 0.9278 | 0.9229 | 0.9226 | 0.9759 | 0.8220 | 0.9905 | 0.8544 | 0.4725 | 10248 | 8424 | 1824 | 657 | 50 | 16 | 1808 |
| ClinVar variant rows | 83688 | 83688 | 0 | 33780 | 40870 | 9038 | 34288 | 39218 | 10182 | 0.9417 | 0.9307 | 0.9254 | 0.9282 | 0.9760 | 0.8239 | 0.9912 | 0.8585 | 0.4686 | 10182 | 8389 | 1793 | 603 | 46 | 14 | 1779 |
| HiRO variant rows | 179 | 179 | 0 | 42 | 110 | 27 | 42 | 117 | 20 | 0.8492 | 0.7879 | 0.8188 | 0.6667 | 0.9868 | 0.9000 | 0.9434 | 0.7424 | 0.6536 | 20 | 18 | 2 | 9 | 0 | 0 | 2 |
| eMERGE variant rows | 1974 | 1974 | 0 | 83 | 1793 | 98 | 164 | 1708 | 102 | 0.8896 | 0.5066 | 0.6565 | 0.5816 | 0.9760 | 0.5588 | 0.9781 | 0.5472 | 0.8652 | 102 | 57 | 45 | 40 | 1 | 0 | 45 |
| CardioBoost variant rows | 194 | 194 | 0 | 57 | 0 | 137 | 34 | 73 | 87 | 0.5928 | 0.6424 | 0.4728 | 0.6204 | 0.9649 | 0.9770 | 0.5140 | 0.5361 | 0.3763 | 87 | 85 | 2 | 48 | 4 | 2 | 0 |
| HiRO source records | 482 | 325 | 157 | 167 | 123 | 35 | 119 | 179 | 27 | 0.7662 | 0.7219 | 0.7770 | 0.7143 | 0.9931 | 0.9259 | 0.9664 | 0.7945 | 0.5508 | 27 | 25 | 2 | 10 | 0 | 0 | 2 |
| eMERGE source records | 2754 | 1974 | 780 | 83 | 1793 | 98 | 164 | 1708 | 102 | 0.8896 | 0.5066 | 0.6565 | 0.5816 | 0.9760 | 0.5588 | 0.9781 | 0.5472 | 0.8652 | 102 | 57 | 45 | 40 | 1 | 0 | 45 |
| CardioBoost source records | 355 | 195 | 160 | 57 | 0 | 138 | 34 | 73 | 88 | 0.5949 | 0.6436 | 0.4735 | 0.6232 | 0.9649 | 0.9773 | 0.5140 | 0.5375 | 0.3744 | 88 | 86 | 2 | 48 | 4 | 2 | 0 |

## Confusion Matrix Files

- All variant rows: `results/model_performance/true_three_class_all_sources/all_variant_rows_confusion.tsv`
- ClinVar variant rows: `results/model_performance/true_three_class_all_sources/clinvar_variant_rows_confusion.tsv`
- HiRO variant rows: `results/model_performance/true_three_class_all_sources/hiro_variant_rows_confusion.tsv`
- eMERGE variant rows: `results/model_performance/true_three_class_all_sources/emerge_variant_rows_confusion.tsv`
- CardioBoost variant rows: `results/model_performance/true_three_class_all_sources/cardioboost_variant_rows_confusion.tsv`
- HiRO source records: `results/model_performance/true_three_class_all_sources/hiro_source_records_confusion.tsv`
- eMERGE source records: `results/model_performance/true_three_class_all_sources/emerge_source_records_confusion.tsv`
- CardioBoost source records: `results/model_performance/true_three_class_all_sources/cardioboost_source_records_confusion.tsv`
