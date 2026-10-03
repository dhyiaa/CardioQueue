# Gene 013: DSG2 AF3 Systematic Review

Date: 2026-07-06

## Current Decision

DSG2 is gene #13 in the systematic AF3 queue. It has 2,034 included variant rows, 812 binary rows, and 1,222 VUS rows in the current rescue-aware modeling table. DSG2 is a primary desmosomal cadherin gene for arrhythmogenic cardiomyopathy, so the first-pass AF3 jobs should focus on cardiac desmosome biology, not generic epithelial signaling or viral receptor biology.

The first-pass packet uses three jobs:

1. DSG2 extracellular 50-608 plus DSC2 extracellular 53-690.
2. DSG2 extracellular 50-608 homodimer.
3. DSG2 cytoplasmic tail 635-1118 plus PKP2 full-length.

This packet does not repeat the old full-length DSG2-DSC2 job by default. That prior AF3 output was not interface-usable, so the new packet splits DSG2 into extracellular cadherin adhesion and cytoplasmic plaque-context designs.

## Local Dataset Context

| Metric | Value |
|---|---:|
| Included DSG2 rows | 2,034 |
| Binary rows | 812 |
| VUS rows | 1,222 |
| Pathogenic rows | 129 |
| Benign rows | 683 |

DSG2 has a large VUS burden. If AF3 produces a strong interface for a clean domain-level design, it may be useful for prioritizing DSG2 VUS by whether they land near a confident cadherin adhesion or plaque interface.

## Sequence Audit

| Field | Value |
|---|---|
| DSG2 UniProt accession | Q14126 |
| Reviewed entry | yes |
| Canonical protein length | 1,118 aa |
| DSG2 extracellular proxy used | 50-608 |
| DSG2 cytoplasmic-tail boundary used | 635-1118 |
| DSC2 UniProt accession | Q02487 |
| DSC2 extracellular proxy used | 53-690 |
| PKP2 UniProt accession | Q99959 |
| Baseline sequence policy | reviewed UniProt canonical amino-acid sequence, with documented domain slices |
| PTM/proteoform policy | no calcium occupancy, glycosylation, phosphorylation, membrane, or desmosome mechanical state in baseline |

The extracellular cadherin designs are protein-only approximations. Real desmosomal cadherins are calcium-dependent, glycosylated, membrane-anchored proteins. Therefore, these jobs are useful only if pairwise AF3 confidence is strong and local geometry is plausible.

## Candidate Partner Summary

| Partner | Priority | Why |
|---|---|---|
| DSC2 extracellular domain | top-3 compact rescue | Core desmosomal cadherin partner; prior full-length pair was weak |
| DSG2 extracellular homodimer | top-3 | UniProt-supported homodimerization; clean DSG2 self-adhesion comparator |
| PKP2 | top-3 | Direct desmosomal plaque partner for DSG2 cytoplasmic-tail context |
| Full-length DSG2-DSC2 | hold | Already tested and not interface-usable |
| JUP/DSP plaque context | hold | Relevant but less direct for DSG2 than PKP2 and overlaps other packets |
| CTNNB1 | hold | Broader junction/EMT biology, not top cardiogenetics desmosome context |
| Adenovirus fiber proteins | exclude | Infection biology, not inherited cardiogenetics context |

Detailed partner evidence table:

```text
03_candidate_partner_evidence.tsv
```

## Proposed AF3 Designs

| Design | Residues | Partner | Status |
|---|---:|---|---|
| DSG2-DSC2 extracellular cadherin rescue | DSG2 50-608 | DSC2 53-690 | top-3 draft JSON ready |
| DSG2 extracellular homodimer | DSG2 50-608 x2 | DSG2 self | top-3 draft JSON ready |
| DSG2 cytoplasmic-tail/PKP2 plaque context | DSG2 635-1118 | PKP2 1-881 | top-3 draft JSON ready |
| Full-length DSG2-DSC2 repeat | DSG2 1-1118 | DSC2 1-901 | held, already weak |
| DSG2-JUP/DSP plaque context | DSG2 635-1118 | JUP or DSP | held |

Detailed design table:

```text
04_af3_design_decisions.tsv
```

## Why These Jobs Come First

DSG2 has two biologically distinct surfaces that matter for cardiogenetics: extracellular desmosomal adhesion and intracellular plaque anchoring. The prior full-length DSG2-DSC2 job mixed signal peptide, membrane, extracellular, and cytoplasmic regions in one model and was weak. The new top-3 packet separates those questions.

The first job asks whether a compact DSG2/DSC2 extracellular cadherin pair gives a usable adhesion interface. The second asks whether DSG2 homodimerization is more stable than the heterotypic pair. The third tests the intracellular plaque context with PKP2, the most directly supported DSG2 plaque partner in the local evidence stack.

## Source Links For Manual Review

- UniProt DSG2: https://www.uniprot.org/uniprotkb/Q14126/entry
- UniProt DSC2: https://www.uniprot.org/uniprotkb/Q02487/entry
- UniProt PKP2: https://www.uniprot.org/uniprotkb/Q99959/entry
- UniProt JUP: https://www.uniprot.org/uniprotkb/P14923/entry
- UniProt DSP: https://www.uniprot.org/uniprotkb/P15924/entry
- HPA DSG2: https://www.proteinatlas.org/DSG2
- DSG2/PKP2 interaction reference: https://pubmed.ncbi.nlm.nih.gov/11790773/
- DSC2/desmosome partner reference: https://pubmed.ncbi.nlm.nih.gov/21062920/
- DSG2 homodimer reference: https://pubmed.ncbi.nlm.nih.gov/31845994/
- DSG2/CTNNB1 held-partner reference: https://pubmed.ncbi.nlm.nih.gov/29910125/
- Local prior AF3 DSG2-DSC2 QC: results/af3/current_completed_af3_summary/completed_af3_pair_qc.tsv

Specific Ctrl-F source table:

```text
09_specific_citation_sources.tsv
```

## AF3 Draft Job JSONs

Top-3 upload file:

```text
af3_job_json/AF3_SYS_gene_013_DSG2_TOP3_upload_alphafoldserver.json
```

Individual JSONs:

```text
af3_job_json/AF3_SYS_G013_DSG2_ECTO_50_608_DSC2_ECTO_53_690_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G013_DSG2_ECTO_50_608_HOMODIMER_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G013_DSG2_CYTO_635_1118_PKP2_FULL_LENGTH_DRAFT.alphafoldserver.json
```

Place AlphaFold Server outputs here:

```text
af3_job_results/
```

## What To Watch When AF3 Results Come Back

For DSG2, do not use global ranking alone. Check pairwise ipTM, inter-chain PAE, contact probability, and whether mapped DSG2 variants land near a confident interface. Because cadherin biology depends on calcium, glycosylation, and membrane geometry, weak AF3 outputs should be treated as no usable interface, not as negative biology.
