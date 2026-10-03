# Gene 007: CACNA1C AF3 Systematic Review

Date: 2026-07-06

## Current Decision

CACNA1C is gene #7 in the systematic AF3 queue. It has 3,557 included variant rows, 1,987 binary rows, and 1,536 VUS rows in the current rescue-aware modeling table. The AF3 packet is high priority because CACNA1C is a core L-type calcium-channel gene and because the current variant burden is large enough that even a partial AF3 feature layer could affect many rows.

The first-pass packet uses three jobs only:

1. Full-length CACNA1C with CACNB2, CALM1, and four calcium ions.
2. Full-length CACNA1C with CACNB2 only.
3. CACNA1C IQ/C-terminal calmodulin-binding region 1648-1705 with CALM1 and four calcium ions.

This keeps the submission practical while still separating the beta-subunit signal from the calmodulin/IQ signal.

## Local Dataset Context

| Metric | Value |
|---|---:|
| Included CACNA1C rows | 3,557 |
| Binary rows | 1,987 |
| VUS rows | 1,536 |
| Pathogenic rows | 97 |
| Benign rows | 1,890 |
| Missing-label rows | 34 |
| HiRO-overlap rows | 10 |

Source composition is mostly ClinVar and eMERGE, with small CardioBoost and HiRO overlap. That means CACNA1C AF3 features should be treated as a broad public-data structure layer, not as a HiRO-specific signal.

## Sequence Audit

| Field | Value |
|---|---|
| CACNA1C UniProt accession | Q13936 |
| Reviewed entry | yes |
| Canonical protein length | 2,221 aa |
| CACNB2 UniProt accession | Q08289 |
| CACNB2 protein length | 660 aa |
| CALM1 UniProt accession | P0DP23 |
| CALM1 protein length | 149 aa |
| Baseline sequence policy | reviewed UniProt canonical amino-acid sequence |
| PTM/isoform policy | no phosphorylation, membrane, lipid, voltage-state, or splice-isoform states in baseline |

Important caveat: this is not the literal fully modified channel in a cardiomyocyte membrane. It is a systematic, reproducible protein-sequence baseline. The AF3 output must be judged by interface confidence and biological plausibility before any feature is promoted.

## Candidate Partner Summary

| Partner | Priority | Why |
|---|---|---|
| CACNB2 | top-3 | Core L-type calcium-channel beta-2 subunit; direct channel-regulatory context |
| CALM1 + Ca2+ | top-3 | Calmodulin mediates calcium-dependent inactivation and CACNA1C has annotated CaM-binding regions |
| CACNB2 + CALM1 + Ca2+ | top-3 | Tests combined channel-regulatory context in one full-length draft |
| CACNA2D/gamma family | hold | Biologically relevant but exact sequence/isoform design is not ready |
| STAC2 | hold | Annotated region exists, but not first cardiac cardiogenetics mechanism here |
| AID-CACNB2 compact job | hold optional | Strong domain rationale, but generate after full-length alpha/beta result if needed |

Detailed partner evidence table:

```text
03_candidate_partner_evidence.tsv
```

## Proposed AF3 Designs

| Design | Residues | Partner | Status |
|---|---:|---|---|
| Full channel regulatory complex | CACNA1C 1-2221 | CACNB2 + CALM1 + 4 Ca2+ | top-3 draft JSON ready |
| Alpha-beta pair control | CACNA1C 1-2221 | CACNB2 | top-3 draft JSON ready |
| IQ-CALM compact job | CACNA1C 1648-1705 | CALM1 + 4 Ca2+ | top-3 draft JSON ready |
| AID-beta compact sensitivity | CACNA1C 414-458 | CACNB2 114-460 | held optional |
| Expanded auxiliary channel context | boundary/sequence needed | CACNA2D/gamma/STAC candidates | not JSON-ready yet |

Detailed design table:

```text
04_af3_design_decisions.tsv
```

## Why Full-Length And Domain-Level Both Matter

CACNA1C is different from compact sarcomere proteins. A full-length channel job can expose multiple cytosolic regions at once, but it may also produce confusing domain placement because the AlphaFold Server job does not include a membrane or voltage state. The IQ-calmodulin peptide job is much smaller and easier to interpret, but it cannot capture long-range channel context.

So the packet intentionally includes both. The full-length jobs are broad screens. The IQ job is the clean mechanistic test.

## Source Links For Manual Review

- UniProt CACNA1C: https://www.uniprot.org/uniprotkb/Q13936/entry
- UniProt CACNB2: https://www.uniprot.org/uniprotkb/Q08289/entry
- UniProt CALM1: https://www.uniprot.org/uniprotkb/P0DP23/entry
- UniProt CALM2: https://www.uniprot.org/uniprotkb/P0DP24/entry
- CACNA1C AID support: https://pubmed.ncbi.nlm.nih.gov/15141227/
- CACNA1C IQ/CALM support: https://pubmed.ncbi.nlm.nih.gov/16299511/ ; https://pubmed.ncbi.nlm.nih.gov/16338416/ ; https://pubmed.ncbi.nlm.nih.gov/19279214/ ; https://pubmed.ncbi.nlm.nih.gov/20953164/
- N-terminal CALM-binding support: https://pubmed.ncbi.nlm.nih.gov/22518098/
- Calmodulin/CACNA1C calcium-dependent inactivation support: https://pubmed.ncbi.nlm.nih.gov/26969752/
- L-type calcium-channel regulation support: https://pubmed.ncbi.nlm.nih.gov/31454269/
- CACNB2 membrane-targeting/function support: https://pubmed.ncbi.nlm.nih.gov/17525370/ ; https://pubmed.ncbi.nlm.nih.gov/36424916/
- Reactome: https://reactome.org/
- RCSB/PDB: https://www.rcsb.org/

Specific Ctrl-F source table:

```text
09_specific_citation_sources.tsv
```

## AF3 Draft Job JSONs

Top-3 upload file:

```text
af3_job_json/AF3_SYS_gene_007_CACNA1C_TOP3_upload_alphafoldserver.json
```

Individual JSONs:

```text
af3_job_json/AF3_SYS_G007_CACNA1C_CACNB2_CALM1_CA4_FULL_LENGTH_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G007_CACNA1C_CACNB2_FULL_LENGTH_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G007_CACNA1C_IQ_1648_1705_CALM1_CA4_DRAFT.alphafoldserver.json
```

Place AlphaFold Server outputs here:

```text
af3_job_results/
```

## What To Watch When AF3 Results Come Back

For CACNA1C, do not promote a full-length job just because it returns a folded channel-like object. We need interface-specific evidence: pairwise ipTM, inter-chain PAE, contact confidence, and whether disease-relevant variants map near a plausible CACNB2 or CALM1-facing surface. If full-length alpha/beta is weak or messy, the next sensitivity job should be the compact AID-CACNB2 design documented in the held table.
