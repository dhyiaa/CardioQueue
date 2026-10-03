# Gene 011: PKP2 AF3 Systematic Review

Date: 2026-07-06

## Current Decision

PKP2 is gene #11 in the systematic AF3 queue. It has 2,089 included variant rows, 1,025 binary rows, and 1,062 VUS rows in the current rescue-aware modeling table. It is a core arrhythmogenic cardiomyopathy/desmosomal gene with enough pathogenic rows for meaningful subgroup review.

The first-pass packet uses three jobs:

1. PKP2 + DSP N-terminal 1-584 + JUP reduced desmosomal plaque complex.
2. PKP2 + JUP full-length pair control.
3. PKP2 + DSP N-terminal 1-584 pair control.

This packet intentionally avoids the full desmosome core and full-length DSP. The full desmosome core is too large and biologically complicated for a first-pass PKP2 packet. The DSP N-terminal 1-584 region is used because UniProt annotates it as the interaction region with JUP and PKP2.

## Local Dataset Context

| Metric | Value |
|---|---:|
| Included PKP2 rows | 2,089 |
| Binary rows | 1,025 |
| VUS rows | 1,062 |
| Pathogenic rows | 341 |
| Benign rows | 684 |

PKP2 is model-relevant and biologically important. Later AF3 features should focus on whether a variant lands near the PKP2-DSP/JUP plaque interface, not on broad epithelial or generic signaling interactions.

## Sequence Audit

| Field | Value |
|---|---|
| PKP2 UniProt accession | Q99959 |
| Reviewed entry | yes |
| Canonical protein length | 881 aa |
| DSP UniProt accession | P15924 |
| DSP region used | 1-584 |
| JUP UniProt accession | P14923 |
| JUP protein length | 745 aa |
| Baseline sequence policy | reviewed UniProt canonical amino-acid sequence |
| PTM/isoform policy | no phosphorylation, methylation, isoform, membrane, or desmosomal mechanical state in baseline |

## Candidate Partner Summary

| Partner | Priority | Why |
|---|---|---|
| DSP N-terminal 1-584 + JUP | top-3 lead | Best reduced desmosomal plaque context for PKP2 |
| JUP | top-3 | Direct PKP2-JUP plaque interaction, clean pair control |
| DSP N-terminal 1-584 | top-3 | Direct PKP2-DSP interface proxy with annotated DSP boundary |
| DSC2/DSG2 cytoplasmic tails | hold | Plausible cadherin-plaque bridge, but better handled in DSC2/DSG2 packets |
| SCN5A/ANK3/GJA1 | hold | Interesting arrhythmia/conduction context, but local evidence is by similarity |
| Generic IntAct proteomics hits | exclude | Too noisy and not cardiac-desmosome specific |

Detailed partner evidence table:

```text
03_candidate_partner_evidence.tsv
```

## Proposed AF3 Designs

| Design | Residues | Partner | Status |
|---|---:|---|---|
| Reduced desmosomal plaque | PKP2 1-881 | DSP 1-584 + JUP 1-745 | top-3 draft JSON ready |
| PKP2-JUP pair | PKP2 1-881 | JUP 1-745 | top-3 draft JSON ready |
| PKP2-DSP pair | PKP2 1-881 | DSP 1-584 | top-3 draft JSON ready |
| PKP2-DSC2 cytoplasmic tail | PKP2 1-881 | DSC2 716-901 | held |
| PKP2-DSG2 cytoplasmic tail | PKP2 1-881 | DSG2 635-1118 | held |
| PKP2-SCN5A/GJA1/ANK3 | boundary needed | conduction partners | held |

Detailed design table:

```text
04_af3_design_decisions.tsv
```

## Why The Reduced Plaque Job Comes First

PKP2 is not mainly a soluble enzyme or a simple one-partner protein. It sits in the desmosomal plaque. The most relevant first-pass biology is therefore a reduced plaque assembly with DSP and JUP. Full-length DSP is 2,871 aa and would dominate the job. Using DSP 1-584 keeps the annotated PKP2/JUP interaction region while keeping the job interpretable.

The pair-control jobs are included because they help separate whether a strong interface comes from PKP2-JUP, PKP2-DSP, or only the combined trimer context.

## Source Links For Manual Review

- UniProt PKP2: https://www.uniprot.org/uniprotkb/Q99959/entry
- UniProt DSP: https://www.uniprot.org/uniprotkb/P15924/entry
- UniProt JUP: https://www.uniprot.org/uniprotkb/P14923/entry
- HPA PKP2: https://www.proteinatlas.org/PKP2
- PKP2 binding/function paper: https://pubmed.ncbi.nlm.nih.gov/11790773/
- PKP2 desmosome localization paper: https://pubmed.ncbi.nlm.nih.gov/22781308/
- PKP2/ARVC desmosome assembly paper: https://pubmed.ncbi.nlm.nih.gov/19533476/
- PKP2-DSP IntAct/local evidence reference: https://pubmed.ncbi.nlm.nih.gov/25225338/
- DSP/JUP ARVC desmosome reference: https://pubmed.ncbi.nlm.nih.gov/16917092/

Specific Ctrl-F source table:

```text
09_specific_citation_sources.tsv
```

## AF3 Draft Job JSONs

Top-3 upload file:

```text
af3_job_json/AF3_SYS_gene_011_PKP2_TOP3_upload_alphafoldserver.json
```

Individual JSONs:

```text
af3_job_json/AF3_SYS_G011_PKP2_DSP_NTERM_JUP_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G011_PKP2_JUP_FULL_LENGTH_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G011_PKP2_DSP_NTERM_1_584_DRAFT.alphafoldserver.json
```

Place AlphaFold Server outputs here:

```text
af3_job_results/
```

## What To Watch When AF3 Results Come Back

For PKP2, the useful signal will be specific interface confidence: pairwise ipTM, inter-chain PAE, contact confidence, and whether mapped variants sit near the PKP2-DSP or PKP2-JUP interface. Do not use generic global confidence alone. If the trimer is weak but one pair-control is strong, use the pair-control result as the safer AF3 feature source.
