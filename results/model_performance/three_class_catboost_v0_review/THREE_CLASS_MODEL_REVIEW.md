# Three-Class CatBoost V0 Review

Date: 2026-07-02

This is a separate exploratory 3-class branch: Benign / VUS / Pathogenic. It does not replace the primary binary model. VUS is an evidence-status label, so results are useful for triage but should be interpreted more cautiously than binary P/LP vs B/LB performance.

## Data Used

- Total corrected matrix rows: `85677`
- Trainable 3-class rows: `84721`
- Label counts: Benign `33859`, VUS `41731`, Pathogenic `9131`
- Weights: source confidence x class weight; class weights are Benign `1.0`, VUS `0.9`, Pathogenic `2.5`.

## Source-Held-Out Split

| split | Benign | VUS | Pathogenic | rows | eligible_rows |
| --- | --- | --- | --- | --- | --- |
| external_cardioboost | 53 | 0 | 132 | 185 | 185 |
| external_emerge | 80 | 1775 | 96 | 1951 | 1951 |
| external_hiro | 42 | 110 | 27 | 179 | 179 |
| not_eligible | 64 | 630 | 17 | 956 | 0 |
| train | 28631 | 33869 | 7545 | 70045 | 70045 |
| validation | 5053 | 5977 | 1331 | 12361 | 12361 |

## Model-Level Results

| model | split | rows | benign | vus | pathogenic | accuracy | macro_f1 | weighted_f1 | weighted_kappa | multiclass_mcc | ovr_macro_auroc | benign_f1 | vus_f1 | pathogenic_f1 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Internal grouped 3-class | train | 59304 | 23701 | 29211 | 6392 | 0.9389 | 0.9223 | 0.9394 | 0.9273 | 0.8975 | 0.9889 | 0.9608 | 0.9375 | 0.8687 |
| Internal grouped 3-class | validation | 12708 | 5077 | 6254 | 1377 | 0.9280 | 0.9015 | 0.9283 | 0.9122 | 0.8784 | 0.9827 | 0.9602 | 0.9269 | 0.8175 |
| Internal grouped 3-class | test | 12709 | 5081 | 6266 | 1362 | 0.9294 | 0.9044 | 0.9299 | 0.9152 | 0.8808 | 0.9840 | 0.9603 | 0.9280 | 0.8249 |
| Source-held-out 3-class | train | 70045 | 28631 | 33869 | 7545 | 0.9432 | 0.9283 | 0.9436 | 0.9326 | 0.9048 | 0.9899 | 0.9633 | 0.9408 | 0.8809 |
| Source-held-out 3-class | validation | 12361 | 5053 | 5977 | 1331 | 0.9338 | 0.9101 | 0.9343 | 0.9232 | 0.8888 | 0.9866 | 0.9650 | 0.9306 | 0.8348 |
| Source-held-out 3-class | external_cardioboost | 185 | 53 | 0 | 132 | 0.5946 | 0.4693 | 0.7300 | 0.6351 | 0.4339 |  | 0.6429 | 0.0000 | 0.7650 |
| Source-held-out 3-class | external_emerge | 1951 | 80 | 1775 | 96 | 0.8893 | 0.6507 | 0.8996 | 0.4977 | 0.4713 | 0.9242 | 0.4500 | 0.9380 | 0.5641 |
| Source-held-out 3-class | external_hiro | 179 | 42 | 110 | 27 | 0.8492 | 0.8188 | 0.8469 | 0.7879 | 0.7160 | 0.9426 | 0.8095 | 0.8811 | 0.7660 |
| Gene-stress 3-class | train | 50698 | 19764 | 24726 | 6208 | 0.9370 | 0.9222 | 0.9375 | 0.9275 | 0.8963 | 0.9892 | 0.9622 | 0.9348 | 0.8697 |
| Gene-stress 3-class | validation | 10864 | 4233 | 5299 | 1332 | 0.9260 | 0.9054 | 0.9264 | 0.9152 | 0.8773 | 0.9826 | 0.9595 | 0.9234 | 0.8333 |
| Gene-stress 3-class | sparse_gene_stress_test | 12295 | 5626 | 6405 | 264 | 0.8737 | 0.6902 | 0.8967 | 0.8199 | 0.7827 | 0.9574 | 0.9583 | 0.8696 | 0.2427 |
| Gene-stress 3-class | well_represented_test | 10864 | 4236 | 5301 | 1327 | 0.9277 | 0.9059 | 0.9283 | 0.9163 | 0.8803 | 0.9834 | 0.9627 | 0.9256 | 0.8296 |

## Source-Held-Out Confusion Matrices


### external_hiro

| true_label | pred_Benign | pred_VUS | pred_Pathogenic |
| --- | --- | --- | --- |
| true_Benign | 34 | 8 | 0 |
| true_VUS | 8 | 100 | 2 |
| true_Pathogenic | 0 | 9 | 18 |

### external_emerge

| true_label | pred_Benign | pred_VUS | pred_Pathogenic |
| --- | --- | --- | --- |
| true_Benign | 54 | 26 | 0 |
| true_VUS | 105 | 1626 | 44 |
| true_Pathogenic | 1 | 40 | 55 |

### external_cardioboost

| true_label | pred_Benign | pred_VUS | pred_Pathogenic |
| --- | --- | --- | --- |
| true_Benign | 27 | 24 | 2 |
| true_VUS | 0 | 0 | 0 |
| true_Pathogenic | 4 | 45 | 83 |

### validation

| true_label | pred_Benign | pred_VUS | pred_Pathogenic |
| --- | --- | --- | --- |
| true_Benign | 4908 | 144 | 1 |
| true_VUS | 208 | 5458 | 311 |
| true_Pathogenic | 3 | 151 | 1177 |

## HiRO Full Source-Record Evaluation

- Total HiRO source records: `482`
- With predictions: `408`
- No prediction: `74`
- Exact 3-class accuracy on predicted source records: `0.7255`
- Benign recall: `0.5947`
- VUS recall: `0.8895`
- Pathogenic recall: `0.6522`
- Pathogenic-to-benign errors: `0`
- Benign-to-pathogenic escalations: `0`

### HiRO source records, all rows

| true_label | pred_Benign | pred_VUS | pred_Pathogenic | pred_No_prediction |
| --- | --- | --- | --- | --- |
| true_Benign | 113 | 77 | 0 | 30 |
| true_VUS | 16 | 153 | 3 | 35 |
| true_Pathogenic | 0 | 16 | 30 | 9 |

### HiRO source records with predictions only

| true_label | pred_Benign | pred_VUS | pred_Pathogenic |
| --- | --- | --- | --- |
| true_Benign | 113 | 77 | 0 |
| true_VUS | 16 | 153 | 3 |
| true_Pathogenic | 0 | 16 | 30 |

## Interpretation

- The 3-class model learns VUS much better than the binary deferral approximation, especially on HiRO.
- HiRO external variant-level performance is promising: accuracy `0.8492`, macro-F1 `0.8188`, VUS F1 `0.8811`, Pathogenic F1 `0.7660`.
- eMERGE external is dominated by VUS rows; VUS F1 is high, but Benign and Pathogenic F1 are modest.
- CardioBoost has no VUS rows, so 3-class metrics there are not directly comparable; many CardioBoost P/LP rows are predicted VUS rather than Benign.
- Sparse-gene pathogenic performance remains the weak spot: sparse-gene Pathogenic F1 is low because pathogenic sparse-gene rows are few and model precision drops.
- The primary paper model should probably remain binary, while this 3-class model can be presented as exploratory VUS triage unless further validation improves it.

## Artifacts

- `results/model_performance/hiro_source_record_true_three_class/hiro_source_record_true_three_class_confusion_all.png`
- `results/model_performance/hiro_source_record_true_three_class/hiro_source_record_true_three_class_confusion_all.tsv`
- `results/model_performance/hiro_source_record_true_three_class/hiro_source_record_true_three_class_confusion_predicted_only.png`
- `results/model_performance/hiro_source_record_true_three_class/hiro_source_record_true_three_class_confusion_predicted_only.tsv`
- `results/model_performance/hiro_source_record_true_three_class/hiro_source_record_true_three_class_predictions.tsv`
- `results/model_performance/hiro_source_record_true_three_class/hiro_source_record_true_three_class_summary.json`
- `results/model_performance/hiro_source_record_true_three_class/hiro_variant_level_true_three_class_predictions.tsv`
- `results/model_performance/three_class_catboost_v0_review/three_class_model_level_performance.tsv`
- `results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued/confusion_matrix_split_three_class_source_heldout_external_cardioboost.tsv`
- `results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued/confusion_matrix_split_three_class_source_heldout_external_emerge.tsv`
- `results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued/confusion_matrix_split_three_class_source_heldout_external_hiro.tsv`
- `results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued/confusion_matrix_split_three_class_source_heldout_train.tsv`
- `results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued/confusion_matrix_split_three_class_source_heldout_validation.tsv`
- `results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued/feature_columns.json`
- `results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued/feature_importance.tsv`
- `results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued/metrics.json`
- `results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued/three_class_catboost_model.cbm`
- `results/models/three_class_catboost_v0_source_heldout_hiro_emerge_rescued/three_class_predictions_split_three_class_source_heldout.tsv`
