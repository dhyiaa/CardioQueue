# Systematic AF3 V2 All-Rows Three-Zone Evaluation

Binary model applied to all B/VUS/P rows using <=0.1 benign-like, 0.1-0.9 VUS/defer, >=0.9 pathogenic-like. VUS/defer is a triage zone, not a learned VUS class.

| dataset | rows | true_benign | true_vus | true_pathogenic | pred_benign_like | pred_vus_defer | pred_pathogenic_like | exact_3zone_accuracy | vus_deferral_rate | p_lp_sensitivity | specificity_vs_benign | ppv_pathogenic_like_vs_binary |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| All variant rows | 85432 | 33923 | 42361 | 9148 | 43897 | 14386 | 27149 | 0.662 | 0.330 | 0.983 | 0.999 | 0.996 |
| ClinVar variant rows | 83706 | 33796 | 40872 | 9038 | 43255 | 13765 | 26686 | 0.667 | 0.328 | 0.984 | 0.999 | 0.996 |
| HiRO variant rows | 213 | 43 | 139 | 31 | 92 | 59 | 62 | 0.559 | 0.381 | 0.935 | 0.953 | 0.935 |
| eMERGE variant rows | 2635 | 130 | 2394 | 111 | 740 | 1085 | 810 | 0.483 | 0.442 | 0.946 | 0.992 | 0.991 |
| CardioBoost variant rows | 194 | 57 | 0 | 137 | 36 | 25 | 133 | 0.830 |  | 0.934 | 0.912 | 0.962 |
| Systematic AF3 any-output rows | 248 | 57 | 72 | 119 | 73 | 45 | 130 | 0.714 | 0.375 | 0.924 | 0.912 | 0.957 |
| Systematic AF3 trusted-pair rows | 152 | 33 | 37 | 82 | 38 | 25 | 89 | 0.737 | 0.378 | 0.902 | 0.909 | 0.961 |
| Systematic AF3 strong-pair rows | 94 | 21 | 24 | 49 | 20 | 14 | 60 | 0.755 | 0.375 | 0.959 | 0.857 | 0.940 |
