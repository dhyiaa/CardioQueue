# Manuscript Figures, Calibration, and Stratified Performance

This bundle summarizes the primary binary CatBoost P/LP vs B/LB model with 0.1/0.9 deferral. The figures and tables are generated from the source-held-out model predictions and the rescue-aware ready matrix.

## External Validation

| subset | rows | pathogenic | benign | auroc | auprc | brier | sensitivity | specificity | ppv | deferral_rate_0_1_0_9 | high_conf_accuracy_0_1_0_9 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| validation | 6384 | 1331 | 5053 | 1.000 | 1.000 | 0.003 | 0.995 | 0.997 | 0.989 | 0.008 | 0.999 |
| external_all | 430 | 255 | 175 | 0.978 | 0.985 | 0.050 | 0.965 | 0.909 | 0.939 | 0.116 | 0.971 |
| external_hiro | 69 | 27 | 42 | 0.989 | 0.983 | 0.062 | 0.963 | 0.905 | 0.867 | 0.072 | 0.969 |
| external_emerge | 176 | 96 | 80 | 0.989 | 0.989 | 0.029 | 0.979 | 0.950 | 0.959 | 0.114 | 0.987 |
| external_cardioboost | 185 | 132 | 53 | 0.956 | 0.983 | 0.066 | 0.955 | 0.849 | 0.940 | 0.135 | 0.956 |

## All-Label Clinical Triage

This table applies the same binary model to Benign, VUS, and Pathogenic rows using three probability zones: `<=0.1` benign-like, `0.1-0.9` deferred/VUS-zone, and `>=0.9` pathogenic-like. This is not a true three-class classifier. It is a clinical triage analysis.

| dataset | rows_with_prediction | true_benign | true_vus | true_pathogenic | pred_benign | pred_vus_defer | pred_pathogenic | exact_3class_accuracy | plp_sensitivity | specificity_non_plp | ppv_plp | vus_deferral_rate | plp_to_benign_errors | benign_to_plp_escalations | vus_to_plp_escalations |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| All variant rows | 85677 | 33923 | 42361 | 9148 | 44161 | 14332 | 27184 | 0.661 | 0.982 | 0.763 | 0.332 | 0.167 | 18 | 39 | 18064 |
| ClinVar variant rows | 83915 | 33796 | 40872 | 9038 | 43479 | 13710 | 26726 | 0.665 | 0.983 | 0.762 | 0.334 | 0.163 | 14 | 35 | 17711 |
| HiRO variant rows | 240 | 43 | 139 | 31 | 103 | 72 | 65 | 0.563 | 0.903 | 0.824 | 0.467 | 0.300 | 0 | 2 | 30 |
| eMERGE variant rows | 2754 | 130 | 2394 | 111 | 810 | 1097 | 847 | 0.471 | 0.946 | 0.720 | 0.129 | 0.398 | 1 | 1 | 706 |
| CardioBoost variant rows | 353 | 57 | 0 | 137 | 70 | 86 | 197 | 0.830 | 0.920 | 0.930 | 0.969 | 0.244 | 3 | 4 | 0 |
| HiRO source records | 408 | 190 | 172 | 46 | 219 | 94 | 95 | 0.642 | 0.891 | 0.851 | 0.432 | 0.230 | 0 | 18 | 36 |
| eMERGE source records | 2754 | 137 | 2490 | 127 | 810 | 1097 | 847 | 0.469 | 0.945 | 0.723 | 0.142 | 0.398 | 1 | 1 | 726 |
| CardioBoost source records | 355 | 156 | 0 | 199 | 70 | 86 | 199 | 0.665 | 0.859 | 0.821 | 0.859 | 0.242 | 5 | 28 | 0 |

## Calibration

| subset | rows | brier | mean_predicted_probability | observed_pathogenic_fraction | calibration_intercept_like_mean_error |
| --- | --- | --- | --- | --- | --- |
| validation | 6384 | 0.003 | 0.211 | 0.208 | 0.003 |
| external_all | 430 | 0.050 | 0.621 | 0.593 | 0.028 |
| external_hiro | 69 | 0.062 | 0.435 | 0.391 | 0.044 |
| external_emerge | 176 | 0.029 | 0.570 | 0.545 | 0.024 |
| external_cardioboost | 185 | 0.066 | 0.740 | 0.714 | 0.026 |

## Feature Coverage

| feature_group | covered_rows | total_rows | coverage_fraction |
| --- | --- | --- | --- |
| VEP | 85660 | 85677 | 1.000 |
| SpliceAI | 77141 | 85677 | 0.900 |
| gnomAD final AF | 62937 | 85677 | 0.735 |
| gnomAD observed | 55387 | 85677 | 0.646 |
| dbNSFP | 43774 | 85677 | 0.511 |
| AlphaMissense direct | 38023 | 85677 | 0.444 |
| FoldX DDG | 313 | 85677 | 0.004 |

## Key Stratified External-All Results

| stratum_type | stratum | rows | pathogenic | benign | auroc | auprc | sensitivity | specificity | ppv |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| consequence_group | missense | 315 | 173 | 142 | 0.988 | 0.990 | 0.977 | 0.915 | 0.934 |
| consequence_group | stop_gained | 32 | 30 | 2 | 0.883 | 0.992 | 1.000 | 0.500 | 0.968 |
| dbnsfp_coverage_group | dbNSFP_matched | 348 | 214 | 134 | 0.985 | 0.989 | 0.981 | 0.918 | 0.950 |
| dbnsfp_coverage_group | dbNSFP_missing | 82 | 41 | 41 | 0.954 | 0.965 | 0.878 | 0.878 | 0.878 |
| gene_group | Arrhythmia | 265 | 157 | 108 | 0.984 | 0.987 | 0.975 | 0.917 | 0.944 |
| gene_group | Cardiomyopathy | 151 | 92 | 59 | 0.968 | 0.982 | 0.946 | 0.898 | 0.935 |
| gnomad_status_group | gnomAD_observed | 326 | 175 | 151 | 0.984 | 0.983 | 0.977 | 0.921 | 0.934 |
| gnomad_status_group | gnomAD_not_joined | 72 | 67 | 5 | 0.997 | 1.000 | 0.985 | 1.000 | 1.000 |
| gnomad_status_group | gnomAD_confirmed_absent | 32 | 13 | 19 | 0.804 | 0.851 | 0.692 | 0.789 | 0.692 |
| sparse_gene_group | well_represented_PLP_ge30 | 406 | 245 | 161 | 0.977 | 0.985 | 0.963 | 0.901 | 0.937 |
| sparse_gene_group | sparse_PLP_lt30 | 24 | 10 | 14 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

## Generated Files

- `figures/primary_binary_roc_pr_panel.png`
- `figures/primary_binary_calibration_reliability.png`
- `figures/source_flow_modeling_rows.png`
- `figures/feature_group_coverage_manuscript.png`
- `figures/external_all_stratified_performance_heatmap.png`
- `figures/vus_three_zone_triage_distribution.png`
- `tables/external_validation_performance.tsv`
- `tables/primary_binary_three_zone_all_labels_summary.tsv`
- `tables/primary_binary_three_zone_all_labels_confusions_long.tsv`
- `tables/calibration_brier_summary.tsv`
- `tables/calibration_curve_points.tsv`
- `tables/source_flow_modeling_rows.tsv`
- `tables/feature_group_coverage_manuscript.tsv`
- `tables/stratified_performance.tsv`
