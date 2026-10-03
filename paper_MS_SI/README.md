# BMC Medical Genomics manuscript package

## Primary files

| File | Role |
|---|---|
| `MS_BMC_MEDICAL_GENOMICS_FINAL.md` | Research article formatted for *BMC Medical Genomics* |
| `SI_BMC_MEDICAL_GENOMICS_FINAL.md` | Additional file 1, expanded methods and results |
| `COMPLETION_AUDIT.md` | Verification record and author-only submission fields |
| `environment-manuscript.txt` | Package versions used to regenerate manuscript statistics |
| `ARTIFACT_CHECKSUMS.sha256` | Checksums for 219 frozen core, sensitivity, structure, FoldX, temporal-VUS, calibration, review-queue, and final multi-source artifacts |
| `literature_evidence/BMC_MEDICAL_GENOMICS_COMPARATOR_REVIEW.md` | Review of 10 related articles from the target journal |
| `literature_evidence/NOVELTY_POSITIONING.md` | Supported novelty claims and claims to avoid |
| `../results/model_performance/clinvar_structure_audit/CLINVAR_STRUCTURE_AUDIT.md` | ClinVar review-quality and protein-structure provenance audit |
| `../results/model_performance/structure_enhancement/STRUCTURE_ENHANCEMENT_REPORT.md` | Registry-wide AlphaFold, DSSP, and FreeSASA sensitivity analysis |
| `../results/model_performance/structure_enhancement/STRUCTURE_GENE_STRESS_REPORT.md` | Structure-enhanced evaluation on genes absent from training |
| `../results/model_performance/maximized_foldx_enhancement/` | Completed FoldX coverage, matched eligibility controls, confidence-tier metrics, and paired gene-clustered analyses |
| `../results/model_performance/calibration_sensitivity/` | Internal-validation calibration sensitivity and paired external evaluation |
| `../results/model_performance/review_queue_simulation/REVIEW_QUEUE_REPORT.md` | Held-out score-ranked worklist simulation and limitations |
| `../results/model_performance/review_queue_simulation/public_source_vus_research_worklist.tsv` | Research-only public-source VUS queue with HiRO-linked rows excluded |
| `../results/model_performance/final_model_heldout_2026/` | Final 775-allele primary, source-specific, conflict-aware, and 64-allele consensus results |
| `../datasets/external_validation_candidates/final_model_heldout_2026/` | Source assertions, unique manifest, conflicts, eligibility, provenance, and overlap accounting |
| `tables/` | Machine-readable supplementary tables |
| `figures/` | Main and supplementary figure files |

The older `MS_GIM_DRAFT.md` and `SI_GIM_DRAFT.md` files are retained as historical drafts. They are not the submission version.

The pre-enhancement BMC package, primary model, sensitivity model, and edited scripts are preserved under `backups/pre_publishability_code_audit_20260917/` with checksums for 76 archived artifacts.

The package immediately before review-queue reframing is preserved under `backups/pre_review_queue_reframing_20260917/` with a separate checksum manifest.

## Submission boundary

The scientific and editorial package is complete. Bracketed fields in the main manuscript require facts held by the authors. They must be resolved before submission; no ethics, funding, authorship, or data-access details were inferred.
