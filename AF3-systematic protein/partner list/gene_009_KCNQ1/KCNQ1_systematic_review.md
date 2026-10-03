# Gene 009: KCNQ1 AF3 Systematic Review

Date: 2026-07-06

## Current Decision

KCNQ1 is gene #9 in the systematic AF3 queue. It has 2,792 included variant rows, 1,916 binary rows, and 852 VUS rows in the current rescue-aware modeling table. It is a high-priority long-QT/cardiac repolarization gene and has one of the clearest cardiology-specific partner contexts: KCNQ1 with KCNE1 in the IKs channel.

The first-pass packet uses three jobs:

1. Tetrameric KCNQ1-KCNE1-CALM1-Ca IKs regulatory complex.
2. KCNQ1-KCNE1 full-length pair control.
3. KCNQ1-CALM1-Ca full-length pair control.

This packet intentionally keeps KCNE1 included for AF3, even though KCNE1 had mixed include/quarantine status in other modeling tables. That is because you explicitly asked to include KCNE1 in the AF3 partner workflow, and biologically it is the central IKs auxiliary subunit.

## Local Dataset Context

| Metric | Value |
|---|---:|
| Included KCNQ1 rows | 2,792 |
| Binary rows | 1,916 |
| VUS rows | 852 |
| Pathogenic rows | 534 |
| Benign rows | 1,382 |
| Missing-label rows | 24 |
| HiRO-overlap rows | 7 |

KCNQ1 has enough benign and pathogenic labels for later gene-stratified analysis. The older `protein_position` field is sparse, so final AF3 mapping should use VEP/UniProt/dbNSFP reconciliation.

## Sequence Audit

| Field | Value |
|---|---|
| KCNQ1 UniProt accession | P51787 |
| Reviewed entry | yes |
| Canonical protein length | 676 aa |
| KCNE1 UniProt accession | P15382 |
| KCNE1 protein length | 129 aa |
| CALM1 UniProt accession | P0DP23 |
| CALM1 protein length | 149 aa |
| Baseline sequence policy | reviewed UniProt canonical amino-acid sequence |
| PTM/isoform policy | no membrane, voltage-state, PIP2, phosphorylation, or isoform-2 modulation in baseline |

## Candidate Partner Summary

| Partner | Priority | Why |
|---|---|---|
| KCNE1 | top-3 lead | Core IKs beta subunit; strongest cardiac-specific KCNQ1 partner |
| CALM1 + Ca2+ | top-3 | KCNQ1 has annotated calmodulin interaction regions and CALM regulates channel activity |
| KCNQ1-KCNE1-CALM1-Ca2+ | top-3 | Combined IKs regulatory context |
| KCNE2 | hold | KCNQ1-KCNE2 exists but is not the primary IKs cardiogenetics context |
| AKAP9 | hold/excluded-gene context | KCNQ1 has AKAP9 region, but AKAP9 is not suitable for the primary cardiac classifier |
| PIP2 | hold ligand review | Important channel ligand, but needs separate ligand/membrane setup |

Detailed partner evidence table:

```text
03_candidate_partner_evidence.tsv
```

## Proposed AF3 Designs

| Design | Residues | Partner | Status |
|---|---:|---|---|
| Tetrameric IKs regulatory complex | KCNQ1 1-676 x4 | KCNE1 x4 + CALM1 x4 + 16 Ca2+ | top-3 draft JSON ready |
| KCNQ1-KCNE1 pair control | KCNQ1 1-676 | KCNE1 1-129 | top-3 draft JSON ready |
| KCNQ1-CALM1 calcium pair control | KCNQ1 1-676 | CALM1 + 4 Ca2+ | top-3 draft JSON ready |
| Compact CALM-region sensitivity | KCNQ1 360-390 and/or 505-535 | CALM1 + 4 Ca2+ | held optional |
| Compact C-terminal KCNE1 sensitivity | KCNQ1 530-580 | KCNE1 109-129 | held optional |
| KCNQ1-KCNE2 auxiliary sensitivity | KCNQ1 1-676 | KCNE2 1-123 | held |

Detailed design table:

```text
04_af3_design_decisions.tsv
```

## Why The Tetrameric Job Comes First

The older pilot-v2 jobs were useful but simplified: one KCNQ1, one KCNE1, one CALM1, and one calcium ion. The systematic packet upgrades the lead job to a tetrameric IKs draft because KCNQ1 is a channel alpha subunit and IKs is not truly a single-chain pair in biology.

However, the tetrameric job is still a draft. It lacks membrane, voltage state, PIP2, and phosphorylation context. That is why the pair-control jobs are included. If the tetramer is messy, the pair controls may be easier to interpret.

## Source Links For Manual Review

- UniProt KCNQ1: https://www.uniprot.org/uniprotkb/P51787/entry
- UniProt KCNE1: https://www.uniprot.org/uniprotkb/P15382/entry
- UniProt CALM1: https://www.uniprot.org/uniprotkb/P0DP23/entry
- KCNQ1/KCNE1 IKs support: https://pubmed.ncbi.nlm.nih.gov/20533308/ ; https://pubmed.ncbi.nlm.nih.gov/8900283/ ; https://pubmed.ncbi.nlm.nih.gov/9108097/ ; https://pubmed.ncbi.nlm.nih.gov/9312006/
- KCNQ1/CALM support: https://pubmed.ncbi.nlm.nih.gov/25441029/
- KCNQ1/KCNE1 C-terminal interaction support: https://pubmed.ncbi.nlm.nih.gov/25037568/
- KCNQ1 PIP2 structural/ligand context: https://pubmed.ncbi.nlm.nih.gov/31883792/
- Reactome: https://reactome.org/
- RCSB/PDB: https://www.rcsb.org/

Specific Ctrl-F source table:

```text
09_specific_citation_sources.tsv
```

## AF3 Draft Job JSONs

Top-3 upload file:

```text
af3_job_json/AF3_SYS_gene_009_KCNQ1_TOP3_upload_alphafoldserver.json
```

Individual JSONs:

```text
af3_job_json/AF3_SYS_G009_KCNQ1_KCNE1_CALM1_CA16_TETRAMER_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G009_KCNQ1_KCNE1_FULL_LENGTH_PAIR_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G009_KCNQ1_CALM1_CA4_FULL_LENGTH_PAIR_DRAFT.alphafoldserver.json
```

Place AlphaFold Server outputs here:

```text
af3_job_results/
```

## What To Watch When AF3 Results Come Back

For KCNQ1, do not promote the tetrameric job solely because it looks channel-like. Use pairwise ipTM, inter-chain PAE, contact confidence, and variant proximity to the KCNE1 or CALM interfaces. If the tetrameric job is weak, the simpler KCNQ1-KCNE1 and KCNQ1-CALM1 pair jobs become the safer feature sources. PIP2 and AKAP9 remain deferred.
