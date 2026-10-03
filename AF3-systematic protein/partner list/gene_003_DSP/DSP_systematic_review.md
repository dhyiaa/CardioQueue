# Gene 003: DSP AF3 Systematic Review

Date: 2026-07-06

## Current Decision

DSP is gene #3 in the systematic AF3 queue. It is a high-burden cardiogenetics gene, but full-length desmoplakin should not be the default design for every partner.

Reason:

```text
DSP canonical protein length = 2,871 aa
```

DSP is modular: an N-terminal plaque region binds JUP/PKP proteins, a long central rod is mostly structural, and a C-terminal plakin-repeat/globular region anchors intermediate filaments. The systematic AF3 design therefore uses two biology-focused contexts:

1. N-terminal desmosome plaque interactions with JUP and PKP2.
2. C-terminal intermediate-filament context with DES/desmin.

## Local Dataset Context

| Metric | Value |
|---|---:|
| Included DSP rows from priority file | 5333 |
| Current matrix DSP rows | 5333 |
| Binary rows | 2666 |
| VUS rows | 2665 |
| Pathogenic rows | 817 |
| Benign rows | 1849 |
| Current rows with `protein_position` | 17 |
| Current rows with VEP protein position-like value | 4751 |
| Current rows with dbNSFP protein position raw | 2914 |

Important caveat: VEP gives numeric-like DSP protein positions for 4751 rows, but 582 rows are non-numeric or dash values. The AF3 mapper should explicitly flag those before claiming variant-level coverage.

## Sequence Audit

| Field | Value |
|---|---|
| UniProt accession | P15924 |
| Reviewed entry | yes |
| Canonical protein length | 2,871 aa |
| Baseline sequence policy | reviewed UniProt canonical desmoplakin DPI/P15924-1 amino-acid sequence |
| PTM/isoform policy | no PTMs or alternative isoforms in baseline; separate sensitivity jobs only if justified |
| Full-length AF3 | not recommended as first systematic design |

## Candidate Partner Summary

| Partner | Priority | Why |
|---|---|---|
| JUP/plakoglobin | primary | UniProt and PubMed support DSP N-terminal interaction; core desmosome plaque biology |
| PKP2/plakophilin-2 | primary | UniProt and PubMed support DSP/PKP2 interaction; ARVC-relevant desmosome plaque partner |
| PKP2 + JUP | primary context | Prior local AF3 pilot showed trimer context supported DSP-JUP better than full-length pair-only DSP-JUP |
| DES/desmin | primary | Cardiac intermediate-filament partner proxy for DSP C-terminal IF-binding region |
| DSC2 | secondary | Desmosome context, but not promoted as direct DSP design in this gene packet |
| DSG2 | secondary | Desmosome cadherin context, better handled in DSG2/DSC2 packet |

Detailed partner evidence table:

```text
03_candidate_partner_evidence.tsv
```

## Proposed AF3 Designs

| Design | Residues | Partner | Status |
|---|---:|---|---|
| NTERM-JUP | 1-584 | JUP | draft JSON ready |
| NTERM-PKP2 | 1-584 | PKP2 | draft JSON ready |
| NTERM-PKP2-JUP | 1-584 | PKP2 + JUP | primary draft JSON ready |
| CTERM-DES | 1946-2871 | DES | primary draft JSON ready |
| DomainB-DES | 2244-2446 | DES | draft JSON ready, compact IF-region test |
| DomainC-DES | 2609-2822 | DES | sensitivity draft JSON ready |
| Full DSP-JUP | 1-2871 | JUP | not recommended; prior full-length pair weak |
| DSG2/DSC2 cadherin context | boundary needed | DSG2/DSC2 | deferred to cadherin/desmosome packet |

Detailed design table:

```text
04_af3_design_decisions.tsv
```

## Prior AF3 Pilot Lesson Incorporated

The earlier AF3Draft pilot tested full-length DSP-related desmosome contexts. The standalone full-length DSP-JUP pair was not interface-usable, while the PKP2-DSP-JUP trimer produced moderate DSP-JUP support. For this systematic packet, we avoid repeating the full-length pair-only design and instead create smaller DSP N-terminal plaque jobs.

## Source Links For Manual Review

- UniProt DSP: https://www.uniprot.org/uniprotkb/P15924/entry
- UniProt JUP: https://www.uniprot.org/uniprotkb/P14923/entry
- UniProt PKP2: https://www.uniprot.org/uniprotkb/Q99959/entry
- UniProt DES: https://www.uniprot.org/uniprotkb/P17661/entry
- DSP-JUP N-terminal paper: https://pubmed.ncbi.nlm.nih.gov/9348293/
- DSP ARVC/JUP interaction paper: https://pubmed.ncbi.nlm.nih.gov/16917092/
- DSP intermediate-filament disease paper: https://pubmed.ncbi.nlm.nih.gov/11063735/

Specific source/citation table for Ctrl-F review:

```text
09_specific_citation_sources.tsv
```

## AF3 Draft Job JSONs

Combined upload candidate:

```text
af3_job_json/AF3_SYS_gene_003_DSP_all_draft_jobs_alphafoldserver.json
```

Manifest:

```text
af3_job_json/AF3_SYS_gene_003_DSP_json_manifest.tsv
```

Jobs intentionally not generated yet:

```text
af3_job_json/AF3_SYS_gene_003_DSP_not_generated_jobs.tsv
```

| AF3 Job | Entities | Upload Status | Evidence IDs |
|---|---|---|---|
| `AF3_SYS_G003_DSP_NTERM_1_584_JUP_DRAFT` | DSP 1-584 + JUP | draft review needed | `SRC_UNIPROT_DSP_JUP_PKP2_REGION`, `SRC_DSP_JUP_PMID9348293` |
| `AF3_SYS_G003_DSP_NTERM_1_584_PKP2_DRAFT` | DSP 1-584 + PKP2 | draft review needed | `SRC_UNIPROT_DSP_JUP_PKP2_REGION`, `SRC_PKP2_DSP_PMID11790773` |
| `AF3_SYS_G003_DSP_NTERM_1_584_PKP2_JUP_DRAFT` | DSP 1-584 + PKP2 + JUP | primary draft review needed | `SRC_PRIOR_AF3_DSP_PKP2_JUP`, `SRC_UNIPROT_DSP_JUP_PKP2_REGION` |
| `AF3_SYS_G003_DSP_CTERM_1946_2871_DES_DRAFT` | DSP 1946-2871 + DES | primary draft review needed | `SRC_UNIPROT_DSP_IF`, `SRC_DSP_IF_PMID30354334` |
| `AF3_SYS_G003_DSP_DOMAIN_B_2244_2446_DES_DRAFT` | DSP 2244-2446 + DES | draft review needed | `SRC_UNIPROT_DSP_IF`, `SRC_DSP_IF_PMID30354334` |
| `AF3_SYS_G003_DSP_DOMAIN_C_2609_2822_DES_DRAFT` | DSP 2609-2822 + DES | sensitivity draft review needed | `SRC_UNIPROT_DSP_IF`, `SRC_DSP_IF_PMID30354334` |

## What To Watch When AF3 Results Come Back

For DSP plaque jobs, a useful result should show the DSP N-terminal fragment contacting JUP and/or PKP2 with acceptable pair ipTM, low interface PAE, and contacts near the expected DSP 1-584 region. For DES jobs, interpret cautiously because desmin biology is filamentous; a protein-only DES chain is a useful proxy but not a full intermediate-filament lattice.

## Next Action For DSP

1. Upload the combined JSON or selected individual JSONs after manual review.
2. Place AlphaFold Server output folders under:

```text
af3_job_results/
```

3. After output is added, run pair-level ipTM/PAE/contact QC and compare against the prior pilot-v2 full-length desmosome jobs.
4. Build the DSP AF3 feature mapper using numeric-like VEP/dbNSFP protein positions and flag dash, range, truncating, and isoform-mismatch rows.
