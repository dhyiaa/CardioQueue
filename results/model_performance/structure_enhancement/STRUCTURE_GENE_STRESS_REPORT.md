# Structure enhancement under gene holdout

The gene-stress test contained 5,890 variants from 19 genes absent from training. Structure fields covered 282 rows across 17 genes.

Across all held-out-gene rows, baseline AUROC was 0.99794 and structure-enhanced AUROC was 0.99858. The gene-clustered difference was 0.00064 (95% CI -0.00003 to 0.00246). AUPRC changed from 0.95743 to 0.95602.

Among structure-covered rows, AUROC changed from 0.98142 to 0.98072. The gene-clustered difference was -0.00070 (95% CI -0.00304 to 0.00332). The analysis does not establish improved discrimination for unseen genes.

Predictions can change on structure-missing rows because the enhanced model is refitted in full. Those changes cannot be attributed to direct structural measurements.
