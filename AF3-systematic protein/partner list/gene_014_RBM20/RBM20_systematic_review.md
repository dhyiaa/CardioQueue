# Gene 014: RBM20 AF3 Systematic Review

Date: 2026-07-06

## Current Decision

RBM20 is gene #14 in the systematic AF3 queue. It has 2,000 included variant rows, 966 binary rows, and 1,034 VUS rows in the current rescue-aware modeling table. RBM20 is a real cardiogenetics gene, especially for dilated cardiomyopathy, but it is not a normal protein-protein interaction target. Its main mechanism is RNA binding and regulation of cardiac pre-mRNA splicing.

Because of that, this packet is intentionally conservative. I did not promote noisy IntAct hits into protein partners. The first-pass JSONs are review-grade RBM20 monomer/domain jobs only:

1. RBM20 390-660, covering the U1-type zinc finger, RRM, RS region, and the clinically important RS hotspot.
2. RBM20 500-610, focused on the RRM domain.
3. RBM20 1140-1227, focused on the C-terminal matrin-type zinc finger and truncation-sensitive C-terminal region.

These jobs may help with domain-level local structure/confidence context. They should not be described as protein-protein interaction evidence.

## Local Dataset Context

| Metric | Value |
|---|---:|
| Included RBM20 rows | 2,000 |
| Binary rows | 966 |
| VUS rows | 1,034 |
| Pathogenic rows | 103 |
| Benign rows | 863 |

RBM20 has many VUS rows, but the AF3 layer should not pretend to solve RNA-splicing biology with a generic PPI model. For now, the safest feature layer is domain membership, disordered-region flags, PTM/hotspot flags, and maybe AF3 local confidence if these monomer/domain jobs produce interpretable structures.

## Sequence Audit

| Field | Value |
|---|---|
| RBM20 UniProt accession | Q5T481 |
| Reviewed entry | yes |
| Canonical protein length | 1,227 aa |
| U1-type zinc finger | 409-443 |
| RRM domain | 518-593 |
| RS region | 628-655 |
| C-terminal matrin-type zinc finger | 1161-1192 |
| Baseline sequence policy | reviewed UniProt canonical amino-acid sequence, with documented domain slices |
| PTM/proteoform policy | no phosphorylation, RNA, snRNP, granule state, or spliceosome assembly state in baseline |

RBM20 phosphorylation is not a minor detail. UniProt notes that phosphorylation regulates localization, including Ser-635 and Ser-637 in the RS region. Therefore, a static unmodified protein structure can miss the clinically relevant mechanism.

## Candidate Partner Summary

| Partner | Priority | Why |
|---|---|---|
| RBM20 390-660 monomer/domain context | top-3 review-grade | Captures U1 zinc finger, RRM, RS region, and disease hotspot |
| RBM20 RRM 500-610 | top-3 review-grade | Focused RNA-recognition motif context |
| RBM20 C-terminal zinc finger 1140-1227 | top-3 review-grade | Captures matrin-type zinc finger and C-terminal truncation-sensitive region |
| RNA 5-UCUU-3 motif | held | Mechanistically central, but RNA design/schema needs review |
| U1/U2 snRNP components | held | UniProt-supported association, but exact AF3 entity set not curated |
| NTAQ1/LNX1/EGLN3 and other IntAct hits | excluded | Noisy/high-throughput and not cardiac-splicing-specific enough |

Detailed partner evidence table:

```text
03_candidate_partner_evidence.tsv
```

## Proposed AF3 Designs

| Design | Residues | Partner | Status |
|---|---:|---|---|
| RBM20 U1 zinc finger/RRM/RS hotspot | 390-660 | none, monomer | review-grade JSON ready |
| RBM20 RRM | 500-610 | none, monomer | review-grade JSON ready |
| RBM20 C-terminal matrin-type zinc finger | 1140-1227 | none, monomer | review-grade JSON ready |
| RBM20 RRM/RS + RNA 5-UCUU-3 | boundary/RNA needed | RNA motif | held |
| RBM20 + U1/U2 snRNP components | entity set needed | snRNP proteins/RNA | held |

Detailed design table:

```text
04_af3_design_decisions.tsv
```

## Why This Is Not A Standard Partner Packet

RBM20 regulates mRNA splicing of cardiac targets such as TTN, CACNA1C, CAMK2D, and PDLIM5. UniProt specifically says RBM20 binds a 5-UCUU-3 intronic motif and associates with components of U1/U2 snRNP complexes. That biology points toward RNA and spliceosome context, not a simple stable protein-protein complex.

The local IntAct table contains candidates such as NTAQ1 and LNX1, but these are mainly high-throughput/two-hybrid or proximity hits. They are not strong enough to become publication-grade AF3 partners for a cardiogenetics classifier.

## Source Links For Manual Review

- UniProt RBM20: https://www.uniprot.org/uniprotkb/Q5T481/entry
- RBM20 familial DCM discovery: https://pubmed.ncbi.nlm.nih.gov/19712804/
- RBM20 regulates titin splicing: https://pubmed.ncbi.nlm.nih.gov/22466703/
- RBM20 cardiac pre-mRNA processing / UCUU motif: https://pubmed.ncbi.nlm.nih.gov/24960161/
- RBM20 RSRSP phosphorylation and localization: https://pubmed.ncbi.nlm.nih.gov/29895960/
- RBM20 nuclear localization mutations: https://pubmed.ncbi.nlm.nih.gov/32840935/
- RBM20 granules / cardiomyopathy: https://pubmed.ncbi.nlm.nih.gov/33188278/
- RBM20 gain-of-function splicing/granule redistribution: https://pubmed.ncbi.nlm.nih.gov/34732726/
- RBM20 phosphorylation/RNA-interactome: https://pubmed.ncbi.nlm.nih.gov/35427468/
- Local IntAct RBM20 evidence: datasets/feature_sources/alphafold3/raw/intact_evidence/Q5T481.tab25.txt

Specific Ctrl-F source table:

```text
09_specific_citation_sources.tsv
```

## AF3 Draft Job JSONs

Top-3 review-grade upload file:

```text
af3_job_json/AF3_SYS_gene_014_RBM20_TOP3_upload_alphafoldserver.json
```

Individual JSONs:

```text
af3_job_json/AF3_SYS_G014_RBM20_ZNF_RRM_RS_390_660_MONOMER_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G014_RBM20_RRM_500_610_MONOMER_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G014_RBM20_CTERM_ZNF_1140_1227_MONOMER_DRAFT.alphafoldserver.json
```

Place AlphaFold Server outputs here:

```text
af3_job_results/
```

## What To Watch When AF3 Results Come Back

For RBM20, useful outputs would be local confidence in the RRM/zinc-finger regions and whether disease-hotspot residues fall in structured or disordered context. Do not use these outputs as interface features. If we later build an RNA-aware RBM20 AF3 design, that should be a separate version with explicit RNA sequence, U1/U2 snRNP entity choices, and a written justification.
