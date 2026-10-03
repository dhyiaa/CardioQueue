# 3-Class CatBoost Threshold Adjustment Report

Adjustment: call `Pathogenic` whenever `prob_Pathogenic >= 0.2`; otherwise choose `Benign` vs `VUS` by the larger of those two probabilities.

This keeps the true 3-class model but shifts the P/LP operating point toward higher sensitivity.

## Summary Metrics

| Dataset | rows | true_benign | true_vus | true_pathogenic | pred_benign | pred_vus | pred_pathogenic | accuracy | macro_f1 | weighted_kappa | plp_sensitivity | plp_specificity | plp_ppv | plp_npv | binary_mcc | plp_to_vus | plp_to_benign | benign_to_plp | vus_to_plp |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| All variant rows | 84721 | 33859 | 41731 | 9131 | 34417 | 35737 | 14567 | 0.8996 | 0.8673 | 0.8912 | 0.9724 | 0.9248 | 0.6095 | 0.9964 | 0.7373 | 210 | 42 | 40 | 5648 |
| ClinVar variant rows | 83688 | 33780 | 40870 | 9038 | 34266 | 35008 | 14414 | 0.9016 | 0.8693 | 0.8940 | 0.9757 | 0.9250 | 0.6118 | 0.9968 | 0.7404 | 182 | 38 | 34 | 5562 |
| HiRO variant rows | 179 | 42 | 110 | 27 | 42 | 109 | 28 | 0.8492 | 0.8287 | 0.8024 | 0.8148 | 0.9605 | 0.7857 | 0.9669 | 0.7639 | 5 | 0 | 0 | 6 |
| eMERGE variant rows | 1974 | 83 | 1793 | 98 | 164 | 1567 | 243 | 0.8435 | 0.6152 | 0.4581 | 0.8265 | 0.9136 | 0.3333 | 0.9902 | 0.4893 | 16 | 1 | 2 | 160 |
| CardioBoost variant rows | 194 | 57 | 0 | 137 | 32 | 44 | 118 | 0.7113 | 0.4973 | 0.6631 | 0.8029 | 0.8596 | 0.9322 | 0.6447 | 0.6183 | 23 | 4 | 8 | 0 |
| HiRO source records | 325 | 167 | 123 | 35 | 119 | 169 | 37 | 0.7723 | 0.7946 | 0.7460 | 0.8857 | 0.9793 | 0.8378 | 0.9861 | 0.8442 | 4 | 0 | 0 | 6 |
| eMERGE source records | 1974 | 83 | 1793 | 98 | 164 | 1567 | 243 | 0.8435 | 0.6152 | 0.4581 | 0.8265 | 0.9136 | 0.3333 | 0.9902 | 0.4893 | 16 | 1 | 2 | 160 |
| CardioBoost source records | 195 | 57 | 0 | 138 | 32 | 44 | 119 | 0.7128 | 0.4977 | 0.6639 | 0.8043 | 0.8596 | 0.9328 | 0.6447 | 0.6192 | 23 | 4 | 8 | 0 |

## Confusion Matrix Files

- All variant rows: `results/model_performance/three_class_threshold_tuning/all_variant_rows_threshold_0_2_confusion.tsv`
- ClinVar variant rows: `results/model_performance/three_class_threshold_tuning/clinvar_variant_rows_threshold_0_2_confusion.tsv`
- HiRO variant rows: `results/model_performance/three_class_threshold_tuning/hiro_variant_rows_threshold_0_2_confusion.tsv`
- eMERGE variant rows: `results/model_performance/three_class_threshold_tuning/emerge_variant_rows_threshold_0_2_confusion.tsv`
- CardioBoost variant rows: `results/model_performance/three_class_threshold_tuning/cardioboost_variant_rows_threshold_0_2_confusion.tsv`
- HiRO source records: `results/model_performance/three_class_threshold_tuning/hiro_source_records_threshold_0_2_confusion.tsv`
- eMERGE source records: `results/model_performance/three_class_threshold_tuning/emerge_source_records_threshold_0_2_confusion.tsv`
- CardioBoost source records: `results/model_performance/three_class_threshold_tuning/cardioboost_source_records_threshold_0_2_confusion.tsv`
