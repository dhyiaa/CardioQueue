# Gene 015: DSC2 AF3 Systematic Review

Date: 2026-07-06

## Current Decision

DSC2 is gene #15 in the systematic AF3 queue. It has 1,853 included variant rows, 747 binary rows, and 1,104 VUS rows in the current rescue-aware modeling table. DSC2 is a primary desmosomal cadherin gene for arrhythmogenic cardiomyopathy. The first-pass AF3 packet should therefore focus on cardiac desmosome biology: extracellular cadherin adhesion and cytoplasmic plaque attachment.

The first-pass packet uses three jobs:

1. DSC2 extracellular 53-690 plus DSG2 extracellular 50-608.
2. DSC2 cytoplasmic tail 716-901 plus PKP2 full-length.
3. DSC2 cytoplasmic tail 716-901 plus JUP full-length.

This packet does not repeat the old full-length DSG2-DSC2 job by default. That prior AF3 output was not interface-usable, so the new packet separates extracellular adhesion from plaque-tail interactions.

## Local Dataset Context

| Metric | Value |
|---|---:|
| Included DSC2 rows | 1,853 |
| Binary rows | 747 |
| VUS rows | 1,104 |
| Pathogenic rows | 98 |
| Benign rows | 649 |

DSC2 has a large VUS burden. If AF3 produces a confident interface for either the ectodomain or plaque-tail jobs, those outputs may become useful DSC2-specific VUS prioritization features.

## Sequence Audit

| Field | Value |
|---|---|
| DSC2 UniProt accession | Q02487 |
| Reviewed entry | yes |
| Canonical protein length | 901 aa |
| DSC2 extracellular proxy used | 53-690 |
| DSC2 cytoplasmic-tail boundary used | 716-901 |
| DSG2 UniProt accession | Q14126 |
| DSG2 extracellular proxy used | 50-608 |
| PKP2 UniProt accession | Q99959 |
| JUP UniProt accession | P14923 |
| Baseline sequence policy | reviewed UniProt canonical amino-acid sequence, with documented domain slices |
| PTM/proteoform policy | no calcium occupancy, glycosylation, phosphorylation, membrane, or desmosome mechanical state in baseline |

The extracellular cadherin job is a protein-only approximation. Real desmosomal cadherins are calcium-dependent, glycosylated, and membrane anchored. A weak AF3 result should be treated as no usable interface, not as proof that the interaction is biologically absent.

## Candidate Partner Summary

| Partner | Priority | Why |
|---|---|---|
| DSG2 extracellular domain | top-3 compact rescue | Core desmosomal cadherin partner; prior full-length pair was weak |
| PKP2 | top-3 | Direct DSC2 plaque partner with IntAct/UniProt support |
| JUP | top-3 | Core desmosomal plaque partner with curated DSC2 support |
| Full-length DSC2-DSG2 | hold | Already tested and not interface-usable |
| DSP | hold | Relevant but larger and covered by DSP/PKP2 packets |
| GJA1/DSG3/broad candidates | hold/exclude | Lower priority or less cardiac-specific for first-pass DSC2 AF3 |

Detailed partner evidence table:

```text
03_candidate_partner_evidence.tsv
```

## Proposed AF3 Designs

| Design | Residues | Partner | Status |
|---|---:|---|---|
| DSC2-DSG2 extracellular cadherin rescue | DSC2 53-690 | DSG2 50-608 | top-3 draft JSON ready |
| DSC2 cytoplasmic-tail/PKP2 plaque context | DSC2 716-901 | PKP2 1-881 | top-3 draft JSON ready |
| DSC2 cytoplasmic-tail/JUP plaque context | DSC2 716-901 | JUP 1-745 | top-3 draft JSON ready |
| Full-length DSC2-DSG2 repeat | DSC2 1-901 | DSG2 1-1118 | held, already weak |
| DSC2 cytoplasmic-tail/DSP context | DSC2 716-901 | DSP 1-584 | held |

Detailed design table:

```text
04_af3_design_decisions.tsv
```

## Why These Jobs Come First

DSC2 has two clinically relevant surfaces: extracellular desmosomal adhesion and intracellular plaque anchoring. The old full-length DSC2/DSG2 run mixed both surfaces with membrane and unmodeled cadherin chemistry, and it did not pass interface QC. The new jobs split the biology into smaller, interpretable tests.

The DSC2-DSG2 ectodomain job asks whether compact extracellular cadherin regions produce a usable adhesive interface. The PKP2 and JUP jobs ask whether DSC2 tail variants can be interpreted in relation to desmosomal plaque partners.

## Source Links For Manual Review

- UniProt DSC2: https://www.uniprot.org/uniprotkb/Q02487/entry
- UniProt DSG2: https://www.uniprot.org/uniprotkb/Q14126/entry
- UniProt PKP2: https://www.uniprot.org/uniprotkb/Q99959/entry
- UniProt JUP: https://www.uniprot.org/uniprotkb/P14923/entry
- UniProt DSP: https://www.uniprot.org/uniprotkb/P15924/entry
- HPA DSC2: https://www.proteinatlas.org/DSC2
- DSC2/desmosome partner reference: https://pubmed.ncbi.nlm.nih.gov/21062920/
- PKP2/desmosome localization reference: https://pubmed.ncbi.nlm.nih.gov/22781308/
- Local prior AF3 DSG2-DSC2 QC: results/af3/current_completed_af3_summary/completed_af3_pair_qc.tsv

Specific Ctrl-F source table:

```text
09_specific_citation_sources.tsv
```

## AF3 Draft Job JSONs

Top-3 upload file:

```text
af3_job_json/AF3_SYS_gene_015_DSC2_TOP3_upload_alphafoldserver.json
```

Individual JSONs:

```text
af3_job_json/AF3_SYS_G015_DSC2_ECTO_53_690_DSG2_ECTO_50_608_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G015_DSC2_CYTO_716_901_PKP2_FULL_LENGTH_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G015_DSC2_CYTO_716_901_JUP_FULL_LENGTH_DRAFT.alphafoldserver.json
```

Place AlphaFold Server outputs here:

```text
af3_job_results/
```

## What To Watch When AF3 Results Come Back

For DSC2, use pairwise interface evidence: ipTM, inter-chain PAE, contact probability, and whether mapped DSC2 variants land near a confident interface. Do not use global ranking alone. For the extracellular job, remember that missing calcium, glycosylation, and membrane orientation are real limitations.
