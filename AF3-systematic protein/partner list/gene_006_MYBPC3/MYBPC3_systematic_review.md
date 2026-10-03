# Gene 006: MYBPC3 AF3 Systematic Review

Date: 2026-07-06

## Current Decision

MYBPC3 is gene #6 in the systematic AF3 queue. It is a high-priority cardiomyopathy gene with 4,236 included variant rows and 1,041 pathogenic rows, but it should be modeled as domain-level sarcomere interactions, not as one giant full-length thick-filament job.

Reason:

```text
MYBPC3 canonical protein length = 1,274 aa
MYBPC3 is modular: Ig-like and fibronectin-type domains, an N-terminal regulatory region, a C-terminal anchoring region, and phosphorylation/disordered context.
```

The top-3 packet focuses on three mechanisms:

1. N-terminal MYBPC3 regulatory region with MYH7 S2.
2. N-terminal MYBPC3 regulatory region with cardiac actin as an F-actin proxy.
3. C-terminal MYBPC3 thick-filament anchoring region with MYH7 S2 as an exploratory test.

## Local Dataset Context

| Metric | Value |
|---|---:|
| Included MYBPC3 rows from priority/ready table | 4,236 |
| Binary rows | 2,347 |
| VUS rows | 1,872 |
| Pathogenic rows | 1,041 |
| Benign rows | 1,306 |
| Missing label rows in ready model table | 17 |
| Rows with VEP protein position | 4,235 |
| Rows with dbNSFP protein position | 2,024 |
| Rows with older `protein_position` | 48 |

Important mapping decision: MYBPC3 AF3 feature mapping should use VEP/UniProt Q14896 reconciliation, not the older sparse `protein_position` column.

## Sequence Audit

| Field | Value |
|---|---|
| MYBPC3 UniProt accession | Q14896 |
| Reviewed entry | yes |
| Canonical protein length | 1,274 aa |
| Baseline sequence policy | reviewed UniProt canonical MYBPC3 amino-acid sequence |
| PTM/isoform policy | no phosphorylation/PTM states in baseline; note this as a limitation because MYBPC3 regulation is phosphorylation-sensitive |
| Full-length AF3 | held for now; domain-level jobs are more interpretable |

## Candidate Partner Summary

| Partner | Priority | Why |
|---|---|---|
| MYH7 S2 | top-3 | MYBPC3 binds MHC/MYH7, and IntAct supports MYBPC3-MYH7 physical association |
| ACTC1 | top-3 | MYBPC3 binds F-actin/native thin filaments, and IntAct supports MYBPC3-ACTC1 physical association |
| MYH7 S2 with MYBPC3 C-terminus | top-3 exploratory | Tests a separate C-terminal thick-filament anchoring context |
| MYL2/MYL3 | hold | Important for MYH7, but not direct MYBPC3-first biology |
| Expanded thin filament | hold | Biologically relevant but too large/redundant for this top-3 packet |
| Noisy IntAct hits | hold | Not promoted without cardiac-specific support |

Detailed partner evidence table:

```text
03_candidate_partner_evidence.tsv
```

## Proposed AF3 Designs

| Design | Residues | Partner | Status |
|---|---:|---|---|
| Nterm-MYH7-S2 | MYBPC3 1-452 | MYH7 838-963 | top-3 draft JSON ready |
| Nterm-ACTC1 | MYBPC3 1-452 | ACTC1 full length | top-3 draft JSON ready |
| Cterm-MYH7-S2 | MYBPC3 971-1274 | MYH7 838-963 | top-3 draft JSON ready; exploratory |
| Full-MYBPC3-MYH7-S2 | MYBPC3 1-1274 | MYH7 838-963 | held, not generated |
| Expanded thin filament | MYBPC3 1-452 | ACTC1/TPM1/TNNI3/TNNT2 | held, not generated |

Detailed design table:

```text
04_af3_design_decisions.tsv
```

## Why Not Full-Length First

A full-length MYBPC3 job could be submitted later, but it would be harder to interpret. MYBPC3 has multiple Ig/Fn-like domains, a disordered region, and phosphorylation-regulated behavior. A single full-length model may create low-confidence domain packing that looks impressive but is not useful as a variant feature.

The domain-level packet is more reviewer-safe: each job has a clear biological question and a clear feature interpretation.

## Source Links For Manual Review

- UniProt MYBPC3: https://www.uniprot.org/uniprotkb/Q14896/entry
- UniProt MYH7: https://www.uniprot.org/uniprotkb/P12883/entry
- UniProt ACTC1: https://www.uniprot.org/uniprotkb/P68032/entry
- Reactome striated muscle contraction: https://reactome.org/content/detail/R-HSA-390522
- MYBPC3-ACTC1 IntAct/PubMed support: https://pubmed.ncbi.nlm.nih.gov/21569246/
- MYBPC3-MYH7/ACTC1 IntAct/PubMed support: https://pubmed.ncbi.nlm.nih.gov/23980194/
- MYH7-MYBPC3 IntAct/PubMed support: https://pubmed.ncbi.nlm.nih.gov/30623132/
- MYH7 S2 boundary support: https://pubmed.ncbi.nlm.nih.gov/17095604/

Specific Ctrl-F source table:

```text
09_specific_citation_sources.tsv
```

## AF3 Draft Job JSONs

Top-3 upload file:

```text
af3_job_json/AF3_SYS_gene_006_MYBPC3_TOP3_upload_alphafoldserver.json
```

Individual JSONs:

```text
af3_job_json/AF3_SYS_G006_MYBPC3_NTERM_1_452_MYH7_S2_838_963_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G006_MYBPC3_NTERM_1_452_ACTC1_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G006_MYBPC3_CTERM_971_1274_MYH7_S2_838_963_DRAFT.alphafoldserver.json
```

Place AlphaFold Server outputs here:

```text
af3_job_results/
```

## What To Watch When AF3 Results Come Back

For MYBPC3, do not promote a job just because it folds. Promote only if pair confidence, interface PAE, and contact evidence support a biologically plausible MYBPC3-partner interface. The ACTC1 job should be described as an actin proxy, not a complete F-actin/thin-filament model. The C-terminal job is exploratory and should be treated carefully.
