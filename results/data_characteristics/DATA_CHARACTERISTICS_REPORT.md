# Data Characteristics Report

This report summarizes the current final cardiogenetics modeling matrix and its source/feature coverage.

## Current Matrix

| item | path | rows | columns | size_mb | notes |
| --- | --- | --- | --- | --- | --- |
| Final modeling matrix | datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx_vep_hgvs_spliceai_hiro_agg.tsv | 86889 | 463 | 322.76 | One row per unique GRCh38 chrom-pos-ref-alt variant. |
| QC report | datasets/modeling/qc/modeling_matrix_qc_20260702_hiro_agg/QC_REPORT.md |  |  | 0.0114 | Full missingness, label, source, duplicate, and leakage QC. |

## Source Composition

| source | source_records | linked_source_records | unresolved_source_records | unique_resolved_source_variants | variant_rows_in_final_matrix | matrix_rows_also_in_clinvar | matrix_rows_not_in_clinvar | source_specific_review_stars_nonmissing |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ClinVar | 83915 | 83915 | 0 | 83915 | 83915 | 83915 | 0 | 83915 |
| HiRO | 482 | 408 | 73 | 239 | 240 | 133 | 107 | 0 |
| eMERGE | 2754 | 2754 | 0 | 2754 | 2754 | 2 | 2752 | 0 |
| CardioBoost | 355 | 355 | 0 | 353 | 353 | 236 | 117 | 0 |

## Label Balance

| model_label_3class | rows | percent |
| --- | --- | --- |
| VUS | 43458 | 50.016 |
| Benign | 34042 | 39.179 |
| Pathogenic | 9240 | 10.634 |
| <missing> | 149 | 0.171 |

## Label Balance By Source Flag

| source | variant_rows_in_matrix | benign | pathogenic | vus | missing_or_conflict | binary_supervised_rows | pathogenic_fraction_binary |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ClinVar | 83915 | 33827 | 9051 | 40892 | 145 | 42878 | 0.2111 |
| HiRO | 240 | 44 | 31 | 141 | 24 | 75 | 0.4133 |
| eMERGE | 2754 | 137 | 126 | 2490 | 1 | 263 | 0.4791 |
| CardioBoost | 353 | 82 | 143 | 0 | 128 | 225 | 0.6356 |

## Feature Family Coverage

| feature_group | features_checked | rows_with_any_value | coverage_percent |
| --- | --- | --- | --- |
| Genotype / variant identity | 6 | 86889 | 100.0 |
| VEP consequence / HGVS | 4 | 86872 | 99.98 |
| SpliceAI | 4 | 75954 | 87.415 |
| gnomAD population frequency | 3 | 71448 | 82.229 |
| dbNSFP in silico | 6 | 42338 | 48.727 |
| Conservation | 3 | 42338 | 48.727 |
| AlphaMissense / ESM | 3 | 42173 | 48.537 |
| ClinGen | 2 | 84983 | 97.806 |
| Protein domains / AlphaFold | 3 | 3093 | 3.56 |
| DSSP / FreeSASA | 2 | 2178 | 2.507 |
| FoldX DDG | 1 | 1099 | 1.265 |
| HiRO source-record evidence | 6 | 238 | 0.274 |

## Key Subscore Coverage

| feature_group | feature | nonmissing_rows | missing_rows | coverage_percent | meaning |
| --- | --- | --- | --- | --- | --- |
| Genotype / variant identity | variant_id | 86889 | 0 | 100.0 | GRCh38 chrom-pos-ref-alt identity; primary join key. |
| Genotype / variant identity | chrom | 86889 | 0 | 100.0 |  |
| Genotype / variant identity | pos | 86889 | 0 | 100.0 |  |
| Genotype / variant identity | ref | 86889 | 0 | 100.0 |  |
| Genotype / variant identity | alt | 86889 | 0 | 100.0 |  |
| Genotype / variant identity | primary_gene | 86889 | 0 | 100.0 |  |
| VEP consequence / HGVS | vep_worst_consequence | 86872 | 17 | 99.98 | Predicted molecular consequence, used for consequence-aware learning. |
| VEP consequence / HGVS | vep_impact | 86872 | 17 | 99.98 |  |
| VEP consequence / HGVS | vep_hgvsc | 86872 | 17 | 99.98 | Transcript-level HGVS from VEP. |
| VEP consequence / HGVS | vep_hgvsp | 86872 | 17 | 99.98 | Protein-level HGVS from VEP. |
| SpliceAI | vep_SpliceAI_pred_DS_AG | 75954 | 10935 | 87.415 |  |
| SpliceAI | vep_SpliceAI_pred_DS_AL | 75954 | 10935 | 87.415 |  |
| SpliceAI | vep_SpliceAI_pred_DS_DG | 75954 | 10935 | 87.415 |  |
| SpliceAI | vep_SpliceAI_pred_DS_DL | 75954 | 10935 | 87.415 |  |
| gnomAD population frequency | gnomad_final_af | 64472 | 22417 | 74.2 | Overall population allele frequency; strong benign/prior-probability feature. |
| gnomAD population frequency | gnomad_final_popmax_af | 71443 | 15446 | 82.223 | Maximum population allele frequency; useful for ancestry-aware benign evidence. |
| gnomAD population frequency | gnomad_final_homozygote_count | 64477 | 22412 | 74.206 | Homozygote observation count; helps flag variants inconsistent with severe dominant disease. |
| dbNSFP in silico | sift_score | 41734 | 45155 | 48.031 | Sequence conservation/tolerance predictor; lower values imply deleteriousness. |
| dbNSFP in silico | polyphen2_hdiv_score | 41709 | 45180 | 48.003 | Protein impact predictor trained for Mendelian disease discrimination. |
| dbNSFP in silico | revel_score | 41719 | 45170 | 48.014 | Ensemble missense pathogenicity predictor. |
| dbNSFP in silico | metalr_score | 38319 | 48570 | 44.101 | Ensemble in silico pathogenicity score. |
| dbNSFP in silico | fathmm_xf_coding_score | 38636 | 48253 | 44.466 | Coding variant pathogenicity prediction. |
| dbNSFP in silico | cadd_phred | 42338 | 44551 | 48.727 | Genome-wide deleteriousness score. |
| Conservation | gerp_rs | 42334 | 44555 | 48.722 | Evolutionary constraint/conservation. |
| Conservation | phylop100way_vertebrate | 42338 | 44551 | 48.727 | PhyloP conservation across vertebrates. |
| Conservation | phastcons100way_vertebrate | 42338 | 44551 | 48.727 | PhastCons conservation across vertebrates. |
| AlphaMissense / ESM | alphamissense_direct_score | 36773 | 50116 | 42.322 | AlphaMissense missense pathogenicity score from direct table join. |
| AlphaMissense / ESM | alphamissense_dbnsfp_score | 41877 | 45012 | 48.196 | AlphaMissense score via dbNSFP. |
| AlphaMissense / ESM | esm1b_score | 42075 | 44814 | 48.424 | Protein language model variant-effect score via dbNSFP. |
| ClinGen | clingen_gene_validity_max_classification | 84983 | 1906 | 97.806 | Strength of ClinGen gene-disease validity evidence. |
| ClinGen | clingen_variant_evidence_classification | 319 | 86570 | 0.367 | Variant-level expert/VCEP evidence when matched. |
| Protein domains / AlphaFold | uniprot_accession | 3091 | 83798 | 3.557 | Protein accession used for protein-position features. |
| Protein domains / AlphaFold | protein_position | 3093 | 83796 | 3.56 | Mapped residue position. |
| Protein domains / AlphaFold | alphafold_residue_plddt | 2179 | 84710 | 2.508 | AlphaFold per-residue confidence. |
| DSSP / FreeSASA | dssp_secondary_structure_class | 2178 | 84711 | 2.507 | Secondary structure class from AlphaFold/DSSP. |
| DSSP / FreeSASA | freesasa_relative | 2178 | 84711 | 2.507 | Relative solvent accessibility. |
| FoldX DDG | foldx_ddg_kcal_mol | 1099 | 85790 | 1.265 | Predicted stability change from FoldX. |
| HiRO source-record evidence | hiro_source_record_count | 238 | 86651 | 0.274 | Number of HiRO patient/source classification records mapped to this variant. |
| HiRO source-record evidence | hiro_unique_patient_count | 238 | 86651 | 0.274 | Number of unique HiRO patients carrying this variant. |
| HiRO source-record evidence | hiro_pathogenic_record_count | 37 | 86852 | 0.043 | Number of HiRO source records classified as pathogenic. |
| HiRO source-record evidence | hiro_benign_record_count | 51 | 86838 | 0.059 | Number of HiRO source records classified as benign. |
| HiRO source-record evidence | hiro_vus_record_count | 154 | 86735 | 0.177 | Number of HiRO source records classified as VUS. |
| HiRO source-record evidence | hiro_label_discordant_flag | 4 | 86885 | 0.005 | Whether HiRO source records mapped to the variant disagree across 3-class labels. |

## Genotype / Consequence Characteristics

### Ref/Alt Variant Type
| variant_type_by_ref_alt | rows | percent |
| --- | --- | --- |
| SNV | 78660 | 90.529 |
| Deletion | 5187 | 5.97 |
| Insertion/duplication | 2587 | 2.977 |
| MNV | 455 | 0.524 |

### Top VEP Consequences
| vep_worst_consequence | rows | percent |
| --- | --- | --- |
| missense_variant | 40280 | 46.358 |
| synonymous_variant | 17080 | 19.657 |
| intron_variant | 11088 | 12.761 |
| frameshift_variant | 4276 | 4.921 |
| splice_region_variant | 3192 | 3.674 |
| stop_gained | 2345 | 2.699 |
| non_coding_transcript_exon_variant | 2123 | 2.443 |
| 3_prime_UTR_variant | 1726 | 1.986 |
| splice_donor_variant | 912 | 1.05 |
| splice_acceptor_variant | 888 | 1.022 |
| 5_prime_UTR_variant | 841 | 0.968 |
| inframe_deletion | 721 | 0.83 |
| intergenic_variant | 552 | 0.635 |
| inframe_insertion | 305 | 0.351 |
| stop_lost | 221 | 0.254 |
| start_lost | 134 | 0.154 |
| upstream_gene_variant | 54 | 0.062 |
| downstream_gene_variant | 47 | 0.054 |
| protein_altering_variant | 42 | 0.048 |
| stop_retained_variant | 35 | 0.04 |
| <missing> | 17 | 0.02 |
| coding_sequence_variant | 6 | 0.007 |
| incomplete_terminal_codon_variant | 2 | 0.002 |
| transcript_ablation | 2 | 0.002 |

## gnomAD Population Coverage

| gnomad_final_status | rows | percent |
| --- | --- | --- |
| observed | 54173 | 62.347 |
| confirmed_absent | 17275 | 19.882 |
| not_joined | 15441 | 17.771 |

## ClinVar Confidence

| slice | clinvar_review_stars | rows | percent |
| --- | --- | --- | --- |
| all_clinvar_rows | 1.0 | 54410 | 64.839 |
| all_clinvar_rows | 2.0 | 29137 | 34.722 |
| all_clinvar_rows | 3.0 | 368 | 0.439 |
| primary_include_clinvar_rows | 1.0 | 54409 | 64.839 |
| primary_include_clinvar_rows | 2.0 | 29137 | 34.722 |
| primary_include_clinvar_rows | 3.0 | 368 | 0.439 |

## HiRO / eMERGE / CardioBoost Source Identity Audit

| dataset | local_source_identity | source_records | variant_rows_in_matrix | overlap_with_clinvar_variant_rows | source_specific_review_stars_nonmissing | interpretation |
| --- | --- | --- | --- | --- | --- | --- |
| HiRO | Internal CASPER WES / VERDICT source records with patient phenotype and ACMG-style classification fields. | 482 | 240 | 133 | 0 | Real private HiRO data. Some variants overlap public ClinVar by coordinate, but HiRO labels/phenotypes remain distinct source-record evidence. |
| eMERGE | Public eMERGE/Glazer arrhythmia-gene supplement; real cohort variant data from 21,846 participants, genotype-focused in local files. | 2754 | 2754 | 2 | 0 | Real eMERGE source variants. Labels may use ClinVar/manual ACMG evidence per manuscript, but local rows are not ClinVar-imported rows. |
| CardioBoost | Public CardioBoost coordinate-rich variant dataset reprocessed for this study. | 355 | 353 | 236 | 0 | Real public CardioBoost labels/features. Many rows overlap ClinVar coordinates and must be source/leakage controlled. |

Interpretation:

- HiRO is real internal CASPER WES / VERDICT source data with patient-linked phenotype fields and private source classifications. ClinVar overlap means some exact coordinates are also in ClinVar, not that the HiRO records are ClinVar rows.
- eMERGE is real public cohort/supplement variant data from the Glazer/eMERGE arrhythmia study. The local supplement is genotype-focused; the manuscript states variants were classified using a mix of sequencing-center ACMG/AMP review, ClinVar annotations, manual review, and in vitro evidence for selected VUS. That means its labels can use ClinVar evidence, but the rows are not simply imported ClinVar rows.
- CardioBoost is real public CardioBoost coordinate-rich variant data. It has substantial ClinVar coordinate overlap and therefore needs leakage-aware splitting.

## Figures

- `figures/feature_group_coverage.png`
- `figures/gnomad_status_counts.png`
- `figures/label_distribution_by_source.png`
- `figures/protein_structure_status.png`
- `figures/source_records_vs_matrix_rows.png`
- `figures/source_records_vs_matrix_rows_non_clinvar_zoom.png`
- `figures/top_genes_label_distribution.png`
- `figures/vep_consequence_distribution.png`

## Tables

- `tables/clinvar_confidence.tsv`
- `tables/data_inventory.tsv`
- `tables/feature_group_coverage.tsv`
- `tables/feature_subscore_coverage.tsv`
- `tables/gene_label_counts.tsv`
- `tables/gnomad_status_counts.tsv`
- `tables/hiro_482_to_240_flow.tsv`
- `tables/hiro_emerge_cardioboost_source_identity_audit.tsv`
- `tables/hiro_matrix_set_mismatch_audit.tsv`
- `tables/hiro_resolved_variant_source_record_counts.tsv`
- `tables/label_by_source.tsv`
- `tables/label_counts.tsv`
- `tables/qc_feature_group_coverage_status_based.tsv`
- `tables/source_overlap.tsv`
- `tables/source_record_summary.tsv`
- `tables/sparse_pathogenic_genes.tsv`
- `tables/variant_type_counts.tsv`
- `tables/vep_consequence_counts.tsv`

## Short Scientific Takeaway

The dataset is broad and feature-rich for a cardiogenetics variant model: 86,889 unique variant rows, 43,282 binary supervised rows, 412 columns, and strong coverage for VEP/SpliceAI/gnomAD/ClinGen. In silico missense scores cover about half the full variant matrix because many variants are not dbNSFP-scored missense-like variants. Protein-structure features are intentionally sparse and should be treated as optional high-value features for mapped missense/protein-changing variants rather than universal inputs.
