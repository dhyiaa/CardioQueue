# Leakage and Ablation Audit

This audit retrained the primary binary CatBoost model after removing reviewer-sensitive feature families. All runs used the same source-held-out split, sample weights, CatBoost recipe, and binary P/LP vs B/LB training slice.

The purpose is to test whether performance depends on obvious shortcuts such as gene identity, ClinVar/review metadata, VEP consequence/HGVS annotations, or source-record aggregate fields.

## Ablation Definitions

| Ablation | Features | Dropped From Baseline | Description |
|---|---:|---:|---|
| baseline_current | 338 | 0 | Current primary feature-selection rules. |
| no_gene_identity | 328 | 10 | Remove direct and gene-level identity features, including gene symbols and gene validity flags. |
| no_vep_consequence_hgvs | 302 | 36 | Remove VEP/gnomAD consequence, impact, transcript, and HGVS-like annotations while keeping SpliceAI scores. |
| no_clinvar_metadata | 338 | 0 | Remove any residual ClinVar/review/confidence metadata beyond the default leakage filter. |
| no_source_record_aggregates | 300 | 38 | Remove source-record aggregates and HiRO phenotype/source evidence fields. |
| no_hgvs_transcript_identifiers | 316 | 22 | Remove HGVS, transcript, protein-accession, rsID/CAID, and direct identifier-like fields. |
| strict_low_leakage | 246 | 92 | Remove gene identity, consequence/HGVS/transcript identifiers, ClinVar metadata, and source-record aggregates together. |

## Main Performance Summary

| Ablation | Subset | Rows | AUROC | AUPRC | Sens @0.5 | Spec @0.5 | PPV @0.5 | Deferral | High-conf Accuracy |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline_current | validation | 6384 | 1.000 | 1.000 | 0.996 | 0.997 | 0.990 | 0.009 | 0.999 |
| baseline_current | external_hiro | 69 | 0.991 | 0.987 | 0.963 | 0.905 | 0.867 | 0.072 | 0.969 |
| baseline_current | external_emerge | 176 | 0.990 | 0.990 | 0.990 | 0.938 | 0.950 | 0.085 | 0.988 |
| baseline_current | external_cardioboost | 185 | 0.949 | 0.980 | 0.962 | 0.830 | 0.934 | 0.119 | 0.945 |
| baseline_current | external_all | 430 | 0.977 | 0.984 | 0.973 | 0.897 | 0.932 | 0.098 | 0.966 |
| no_gene_identity | validation | 6384 | 1.000 | 1.000 | 0.995 | 0.997 | 0.988 | 0.008 | 0.999 |
| no_gene_identity | external_hiro | 69 | 0.989 | 0.984 | 0.963 | 0.881 | 0.839 | 0.072 | 0.953 |
| no_gene_identity | external_emerge | 176 | 0.991 | 0.991 | 0.990 | 0.925 | 0.941 | 0.085 | 0.988 |
| no_gene_identity | external_cardioboost | 185 | 0.953 | 0.981 | 0.962 | 0.830 | 0.934 | 0.124 | 0.951 |
| no_gene_identity | external_all | 430 | 0.977 | 0.984 | 0.973 | 0.886 | 0.925 | 0.100 | 0.966 |
| no_vep_consequence_hgvs | validation | 6384 | 0.999 | 0.998 | 0.997 | 0.991 | 0.969 | 0.029 | 0.998 |
| no_vep_consequence_hgvs | external_hiro | 69 | 0.989 | 0.980 | 1.000 | 0.929 | 0.900 | 0.116 | 0.984 |
| no_vep_consequence_hgvs | external_emerge | 176 | 0.980 | 0.979 | 0.948 | 0.963 | 0.968 | 0.176 | 0.972 |
| no_vep_consequence_hgvs | external_cardioboost | 185 | 0.965 | 0.986 | 0.955 | 0.736 | 0.900 | 0.162 | 0.955 |
| no_vep_consequence_hgvs | external_all | 430 | 0.973 | 0.980 | 0.957 | 0.886 | 0.924 | 0.160 | 0.967 |
| no_clinvar_metadata | validation | 6384 | 1.000 | 1.000 | 0.996 | 0.997 | 0.990 | 0.009 | 0.999 |
| no_clinvar_metadata | external_hiro | 69 | 0.991 | 0.987 | 0.963 | 0.905 | 0.867 | 0.072 | 0.969 |
| no_clinvar_metadata | external_emerge | 176 | 0.990 | 0.990 | 0.990 | 0.938 | 0.950 | 0.085 | 0.988 |
| no_clinvar_metadata | external_cardioboost | 185 | 0.949 | 0.980 | 0.962 | 0.830 | 0.934 | 0.119 | 0.945 |
| no_clinvar_metadata | external_all | 430 | 0.977 | 0.984 | 0.973 | 0.897 | 0.932 | 0.098 | 0.966 |
| no_source_record_aggregates | validation | 6384 | 1.000 | 1.000 | 0.996 | 0.997 | 0.989 | 0.008 | 0.999 |
| no_source_record_aggregates | external_hiro | 69 | 0.989 | 0.984 | 0.963 | 0.905 | 0.867 | 0.087 | 0.968 |
| no_source_record_aggregates | external_emerge | 176 | 0.991 | 0.991 | 0.990 | 0.925 | 0.941 | 0.097 | 0.987 |
| no_source_record_aggregates | external_cardioboost | 185 | 0.950 | 0.979 | 0.962 | 0.830 | 0.934 | 0.114 | 0.945 |
| no_source_record_aggregates | external_all | 430 | 0.977 | 0.983 | 0.973 | 0.891 | 0.929 | 0.102 | 0.966 |
| no_hgvs_transcript_identifiers | validation | 6384 | 1.000 | 1.000 | 0.995 | 0.998 | 0.991 | 0.009 | 0.999 |
| no_hgvs_transcript_identifiers | external_hiro | 69 | 0.988 | 0.982 | 0.963 | 0.905 | 0.867 | 0.087 | 0.968 |
| no_hgvs_transcript_identifiers | external_emerge | 176 | 0.993 | 0.993 | 0.990 | 0.963 | 0.969 | 0.097 | 0.987 |
| no_hgvs_transcript_identifiers | external_cardioboost | 185 | 0.952 | 0.981 | 0.962 | 0.849 | 0.941 | 0.092 | 0.940 |
| no_hgvs_transcript_identifiers | external_all | 430 | 0.978 | 0.985 | 0.973 | 0.914 | 0.943 | 0.093 | 0.964 |
| strict_low_leakage | validation | 6384 | 0.999 | 0.993 | 0.993 | 0.984 | 0.941 | 0.040 | 0.998 |
| strict_low_leakage | external_hiro | 69 | 0.987 | 0.974 | 1.000 | 0.905 | 0.871 | 0.130 | 0.967 |
| strict_low_leakage | external_emerge | 176 | 0.970 | 0.973 | 0.948 | 0.912 | 0.929 | 0.188 | 0.979 |
| strict_low_leakage | external_cardioboost | 185 | 0.972 | 0.989 | 0.985 | 0.585 | 0.855 | 0.265 | 0.971 |
| strict_low_leakage | external_all | 430 | 0.974 | 0.981 | 0.973 | 0.811 | 0.883 | 0.212 | 0.973 |

## Interpretation

The baseline external-all AUROC was 0.977 and AUPRC was 0.984. Compare each ablation against that row. Large drops after removing one feature family would indicate possible dependence on that family.
The strict_low_leakage run is intentionally harsh. It should not be treated as the final clinical model, but it is useful as a reviewer-facing stress test.

Key findings:

- Removing direct gene identity did not reduce external-all AUROC/AUPRC. Sensitivity was unchanged, while specificity fell slightly from 0.897 to 0.886.
- Removing residual ClinVar metadata changed nothing. This confirms that the default training feature filter had already removed these columns.
- Removing source-record aggregates, including HiRO phenotype/evidence fields, left external-all AUROC/AUPRC essentially unchanged. This supports the claim that the primary variant model is not being driven by patient/source-record shortcuts.
- Removing VEP/consequence/HGVS-style annotations caused only a small external-all drop: AUROC 0.977 to 0.973 and AUPRC 0.984 to 0.980. Deferral increased from 0.098 to 0.160.
- Removing HGVS/transcript/direct identifier-like fields did not hurt the model. External-all AUROC/AUPRC were 0.978/0.985.
- The strict low-leakage model remained strong for ranking: AUROC 0.974 and AUPRC 0.981. However, specificity fell to 0.811, PPV fell to 0.883, and deferral increased to 0.212. This means the full feature-rich model is better calibrated for high-confidence calling, while the strict model is a useful stress test.

Overall, the audit does not show evidence that the primary external performance is explained by a single leakage-prone feature family. The main manuscript can use the feature-rich primary model, with this audit reported as a robustness analysis.

Saved machine-readable outputs:

- `ablation_metrics.tsv`
- `ablation_feature_set_summary.tsv`
- `external_all_ablation_deltas_vs_baseline.tsv`
- `ablation_specs.json`
- `predictions_<ablation>.tsv`
- `external_all_ablation_metrics.png`
