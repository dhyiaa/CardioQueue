# Dataset Source Status

This file is the project-level tracker for the proposal data sources. It records
what is already local, what is API-ready, what is pending, and what QC is needed
before a source can be used for modeling.

## Source Summary

| Proposal point | Source | Role in project | Current status | Primary local documentation |
|---|---|---|---|---|
| A | ClinVar bulk tab-delimited download | Primary public submitted clinical-significance labels and ClinVar metadata | Downloaded, validated, and curated modeling slice built | `datasets/clinvar/README.md` |
| B | gnomAD v4 | Population frequency, popmax AF, homozygote counts, constraint metrics, benign-class proxy | API-ready for targeted panels; example annotations validated | `docs/gnomad_v4_plan.md` |
| C | In silico annotation tools | Feature engineering only; not a variant source or label source | Partially implemented | `docs/insilico/point_c_status_summary.md` |
| D | ClinGen curated evidence, inferred from current proposal text | Expert-curated gene/variant evidence and higher-confidence criterion context | Gene-disease validity API workflow added and example validated; exact Point D wording still needs confirmation | `docs/clingen/point_d_clingen_plan.md` |
| E | Protein structure, AlphaFold, and stability features | Protein structural context, residue confidence, and mutation stability features | AlphaFold and UniProt local; DSSP/FreeSASA installed; FoldXPro installed and one-variant DDG smoke test passed | `docs/protein_structure/point_e_protein_structure_plan.md` |
| Core variant source | eMERGE/Glazer arrhythmia-gene variant dataset | Real variant-level training/evaluation data | Local interim handoff package added and profiled | `datasets/emerge/full_arrhythmia_gene_dataset/README.md` |
| Core variant source | CardioBoost public dataset | Public variant-level training/evaluation data, mostly binary B/P | Local processed package added and profiled; old scoring artifacts removed | `datasets/cardioboost/public_dataset/README.md` |
| Protein-sequence source | CardioBoost V2 / DYNA Zenodo dataset | CardioBoost/ClinVar protein sequence pairs for CM/ARM plus DYNA splicing files | Downloaded from Zenodo, MD5-validated, processed | `datasets/cardioboost/public_dataset/v2_dyna_zenodo_sequence_dataset/README.md` |
| Core variant source | HiRO / CASPER WES / VERDICT | Internal/high-value cardiogenetics variant data with phenotype fields | Local full dataset added and profiled | `datasets/hiro/full_dataset/README.md` |

## Core Variant Training Sources

## Master Variant Registry

| Component | Status | Notes |
|---|---|---|
| Coordinate registry | Built | `datasets/variant_registry/interim/master_variant_registry.tsv`; 161,108 GRCh38-oriented `chrom-pos-ref-alt` rows |
| Source membership table | Built | `datasets/variant_registry/interim/master_variant_registry_sources.tsv`; 161,802 source-membership rows |
| Source flags | Built | Includes `in_clinvar`, `in_hiro`, `in_emerge`, `in_cardioboost`, and `phenotype_available` |
| Source inventory | Built | `datasets/variant_registry/interim/master_variant_registry_source_inventory.tsv`; separates source-row counts from deduplicated coordinate counts |
| Gene panel | Audited keep/remove panel added | Source-union registry panel has 63 genes; modeling keep panel has 59 genes after removing Tier C contamination genes `DMD`, `AKAP9`, `FPGT`, and `KNCH2` |
| Main caveat | Open | HiRO has 482 private/patient-linked modeling rows, but only 361 currently have complete coordinate keys and these collapse to 224 unique coordinate variants; the remaining 121 rows need protein/HGVS fallback mapping |

| Source | Local status | Label structure | Current usable count | Key caveat |
|---|---|---|---:|---|
| ClinVar | Downloaded, checksum-validated, and curated | Benign, VUS, Pathogenic after excluding conflict/other labels from the clean slice | 127,967 curated variants on the audited 59-gene modeling panel; 41,177 high-confidence stars >=2 and 86,790 star-1 lower-confidence rows | Public submitted labels remain noisier than HiRO/expert adjudication; Tier B genes should be sensitivity-analysis candidates. |
| eMERGE/Glazer | Local and profiled | Benign, VUS, Pathogenic | 2,754 variants; 2,650 primary-validation eligible | Public supplement lacks row-level patient phenotype fields. |
| HiRO / CASPER WES / VERDICT | Local and profiled | Benign, VUS, Pathogenic plus phenotype fields | 482 private/patient-linked modeling rows; 361 with complete coordinate keys; 224 deduplicated coordinate variants | Internal/private; do not collapse the source-row count to the coordinate-registry count when describing HiRO. Coordinate overlap with ClinVar is not the same as ClinVar provenance. |
| CardioBoost public | Local and profiled | Benign, Pathogenic; no VUS in processed cohort | 355 deduplicated processed variants | Binary dataset; keep source flag and avoid pretending missing VUS means benign/pathogenic certainty. |
| CardioBoost V2 / DYNA Zenodo | Local and processed | Benign, Pathogenic, VUS protein sequence pairs | 21,472 protein-sequence rows | Not coordinate-rich; use for sequence modeling or as a protein-level auxiliary dataset. |

## Feature Injection Sources

| Source | Feature role | Current status |
|---|---|---|
| gnomAD v4 | AF, popmax, homozygote/hemizygote counts, constraint metrics, benign-proxy sampling | API-ready; examples validated. |
| In silico tools | CADD, REVEL, SIFT, PolyPhen, MetaLR, FATHMM-XF, AlphaMissense, etc. | Partial; AlphaMissense local, CADD API-ready, dbNSFP GRCh38 BGZF downloaded/MD5-validated and full registry dbNSFP feature table built; final merged feature table pending. |
| ClinGen | Gene/variant expert evidence, gene-disease validity, dosage sensitivity | API/local partial; gene validity, dosage, and variant-pathogenicity examples validated. |
| Protein/AlphaFold/UniProt | Residue confidence, structural context, domains, stability features | Partial; AlphaFold and UniProt local, unified protein table generated, DSSP/FreeSASA installed, FoldXPro smoke-tested; batch DDG and residue-level DSSP/FreeSASA merge pending. |

## Point C Snapshot

| In silico source | Status | Notes |
|---|---|---|
| AlphaMissense hg38 | Downloaded and validated | Official hg38 file is local; example join passed. |
| CADD v1.7 GRCh38 | API-ready for small targeted sets | Bulk all-SNV file is too large to download blindly at this stage. |
| dbNSFP v5.3.1 academic branch | Downloaded, checksum-validated, and joined | Required for old-SI parity predictors: REVEL, SIFT, PolyPhen-2 HDIV, MetaLR, FATHMM-XF. Current `dbNSFP5.3.1a_grch38.gz` MD5 matches the official checksum: `afff632f13325c687477d64ca0c63b1a`. Full registry join produced 161,108 rows with 79,630 exact dbNSFP matches. |
| Ensembl VEP 116 GRCh38 | Deferred | Local cache is large; choose local VEP only if needed. |
| SpliceAI | Pending VEP/plugin decision | Add through VEP plugin or another verified precomputed source. |
| ESM1b / ESM-variant | Pending source/mapping confirmation | Needs exact public source and transcript/protein mapping QC before use. |

Current Point C is pipeline-ready and has a validated dbNSFP registry feature
layer. It is not fully modeling-complete until selected dbNSFP columns are
merged/renamed into the final modeling feature table and strict feature
requirements pass.

## Point D Assumption

The checked-in `PROJECT_PROPOSAL.md` does not contain the full A/B/C/D list. It
does, however, identify the ClinGen Evidence Repository as the better source for
expert-curated variant classifications and evidence. Therefore Point D has been
started as an inferred ClinGen curated-evidence source.

If the intended Point D is different, replace `docs/clingen/point_d_clingen_plan.md`
with the correct source-specific plan before implementation.

Current Point D implementation status:

| ClinGen component | Status |
|---|---|
| Gene lookup API | Checked with `MYH7` and `SCN5A` |
| Gene-disease validity API | Doctor passed; endpoint returned 3,640 rows |
| Cardiogenetics seed-panel validity output | Validated; 19 input genes, 65 output rows, 0 missing genes |
| Variant-level Evidence Repository annotation | Implemented for gene panels; 19-gene seed panel produced 189 rows, 179 ok rows |
| Dosage sensitivity annotation | Downloads implemented; gene dosage CSV and GRCh38 TSV downloaded |

## Point E Snapshot

| Component | Status | Notes |
|---|---|---|
| UniProt reviewed-human mapping | Implemented for gene panels | Required before AlphaFold joins |
| AlphaFold DB gene-panel metadata | Validated | 19 seed genes tested: 16 `ok`, 3 `not_found` |
| AlphaFold human bulk tar | Downloaded and inspected | 5,177,506,304 bytes; 47,172 compressed PDB/CIF model files |
| AlphaFold bulk seed-panel coverage | Validated | 19 seed genes tested: 19 `ok`; covers `TTN`, `DSP`, and `RYR2` through segmented models |
| UniProt reviewed human FASTA | Downloaded | 20,431 reviewed human sequences |
| UniProt reviewed human feature TSV | Downloaded | 20,431 reviewed human protein records plus header |
| AlphaFold residue pLDDT | Example validated | Requires `uniprot_accession` plus `protein_position` |
| Delta-delta-G | Smoke-tested | FoldXPro Mac Silicon runs locally; ACTC1 M134T BuildModel smoke test produced `ddg_kcal_mol=-0.168657`. Rosetta cartesian_ddg remains a smaller sensitivity subset option |
| DSSP/FreeSASA structural context | Installed/smoke-tested | `mkdssp 4.6.1` and FreeSASA Python API installed in `.tools/envs/protein-structure`; MYH7 AlphaFold smoke test passed |
| Domain/functional-site features | Pending | Requires UniProt/Pfam/InterPro extraction |

## eMERGE Snapshot

| Component | Status | Notes |
|---|---|---|
| Full arrhythmia-gene variant table | Local and checksum-verified | 2,754 rows, 47 columns; SHA-256 matches the handoff README. |
| Baseline model-ready table | Local | 2,754 rows, 89 columns; harmonized to the previous Cardiogenetics/CatBoost baseline schema. |
| Enriched model-ready table | Local | 2,754 rows, 95 columns; retains additional enrichment columns. |
| Labels | Local | Post-in-vitro labels: 2,490 VUS, 137 Benign, 127 Pathogenic. |
| Primary validation eligibility | Local | 2,650 eligible rows; 104 excluded by canonical/disagreement logic. |
| In silico completeness in model-ready tables | Partial | `gnomad_max_af` is 100%; `cadd_phred` is 95.5%; REVEL/SIFT/PolyPhen/MetaLR/FATHMM-XF/AlphaMissense are currently empty in these tables. |
| Recommended use | Ready as a real variant-level source after schema review | Keep a source flag and use leakage-aware splits/deduplication across ClinVar, HiRO, eMERGE, and CardioBoost. |

## CardioBoost Public Snapshot

| Component | Status | Notes |
|---|---|---|
| Raw public CardioBoost data | Local | Arrhythmia and cardiomyopathy `.RData` variant datasets retained. |
| Old model/scoring artifacts | Removed | Processed `results/`, raw `ml/`, and raw `prediction/` folders were deleted. |
| Clean reannotated baseline table | Local | 355 rows, 89 columns; 199 Pathogenic, 156 Benign. |
| Clean reannotated enriched table | Local | 355 rows, 95 columns with in silico fields populated for most variants. |
| Label structure | Binary | No VUS rows in the processed CardioBoost cohort. |
| Feature completeness in reannotated files | Strong but not complete | `revel` 100%, `metalr` 100%, `sift` 99.7%, `polyphen2_hdiv` 98.3%, `gnomad_max_af` 95.8%, `alphamissense` 94.1%, `fathmm_xf` 93.2%, `cadd_phred` 92.7%. |
| Recommended use | Real public variant-level source | Best for binary B/P training, public pretraining, source-stratified experiments, or careful integration into a three-class model with explicit missing-VUS handling. |

## CardioBoost V2 / DYNA Zenodo Snapshot

| Component | Status | Notes |
|---|---|---|
| Zenodo record | Downloaded and MD5-validated | DOI `10.5281/zenodo.13397296`; license CC-BY-4.0. |
| Raw files | Local | 12 CSV files, including CardioBoost CM/ARM train/test, ClinVar protein B/P and VUS panels, and DYNA splicing files. |
| Processed protein table | Local | 21,472 rows, 16 columns. |
| Protein clinical-label counts | Local | 18,619 VUS, 1,702 Benign, 1,151 Pathogenic. |
| Original CardioBoost protein rows | Local | 1,120 rows: 608 Pathogenic, 512 Benign. |
| ClinVar protein rows | Local | 20,352 rows: 18,619 VUS, 1,190 Benign, 543 Pathogenic. |
| Key limitation | Important | Sequence pairs only; no gene symbol or `chrom-pos-ref-alt` join key. |
| Recommended use | Auxiliary sequence-level source | Useful for protein-sequence modeling/pretraining and VUS sequence-pair experiments, not a direct replacement for coordinate-based CardioBoost/ClinVar tables. |

## QC Standard Across Sources

Every source-specific output table must include:

- source status column;
- missing-reason column;
- source version or retrieval date;
- genome build or coordinate system when relevant;
- join key used for matching;
- validation command or checksum when available.

Silent nulls are not acceptable. A null feature must mean a documented biological
or source-coverage absence, not an API failure, syntax error, schema mismatch, or
unvalidated join.
