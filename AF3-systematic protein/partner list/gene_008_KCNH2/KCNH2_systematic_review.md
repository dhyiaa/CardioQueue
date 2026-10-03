# Gene 008: KCNH2 AF3 Systematic Review

Date: 2026-07-06

## Current Decision

KCNH2 is gene #8 in the systematic AF3 queue. It has 3,478 included variant rows, 1,870 binary rows, and 1,575 VUS rows in the current rescue-aware modeling table. It is a high-priority LQTS/cardiac repolarization gene, but it needs a careful AF3 strategy because the first KCNH2-KCNE2 AF3 run was not interface-usable.

The first-pass packet uses three jobs:

1. Full-length KCNH2 homotetramer.
2. KCNH2 PAS/PAC region 1-150 with KCNH2 CNBHD region 720-870.
3. KCNH2 with KCNE2 as a biologically supported but low-priority comparator.

The lead design is the KCNH2 homotetramer, not KCNE2, because the existing local KCNE2 AF3 output had weak interface confidence.

## Local Dataset Context

| Metric | Value |
|---|---:|
| Included KCNH2 rows | 3,478 |
| Binary rows | 1,870 |
| VUS rows | 1,575 |
| Pathogenic rows | 701 |
| Benign rows | 1,169 |
| Missing-label rows | 33 |
| HiRO-overlap rows | 10 |

KCNH2 has enough pathogenic and benign labels to be useful for later gene-stratified analysis. The current older `protein_position` field is sparse for KCNH2, so AF3 mapping should use the better VEP/UniProt/dbNSFP reconciliation pathway before final feature generation.

## Prior AF3 Evidence

The previous KCNH2-KCNE2 server run is documented locally as `AF3_PILOTV1_REVIEW_KCNH2_KCNE2`. Its best model had:

| Metric | Value |
|---|---:|
| Ranking score | 0.58 |
| ipTM | 0.31 |
| pTM | 0.40 |
| Best pair PAE minimum | about 12.5 |
| Feature call | not interface-usable |

So KCNE2 remains biologically relevant, but we should not use the old KCNE2 interface as a confident AF3 feature.

## Sequence Audit

| Field | Value |
|---|---|
| KCNH2 UniProt accession | Q12809 |
| Reviewed entry | yes |
| Canonical protein length | 1,159 aa |
| KCNE2 UniProt accession | Q9Y6J6 |
| KCNE2 protein length | 123 aa |
| Baseline sequence policy | reviewed UniProt canonical amino-acid sequence |
| PTM/isoform policy | no phosphorylation, glycosylation, membrane, voltage-state, or isoform A/B-USO states in baseline |

## Candidate Partner Summary

| Partner | Priority | Why |
|---|---|---|
| KCNH2 homotetramer | top-3 lead | Direct hERG/Kv11.1 channel assembly context and replacement for weak KCNE2 output |
| KCNH2 PAS/PAC-CNBHD | top-3 | Compact regulatory-domain design using annotated KCNH2 PAS/PAC and cNMP-binding regions |
| KCNE2 | top-3 comparator | Biologically supported IKr auxiliary context, but prior AF3 interface was weak |
| KCNE1 | hold | Possible KCNH2 association, but stronger first-pass role is KCNQ1/IKs |
| RNF207/DNAJB12/DNAJB14/NEDD4L/NDFIP | hold | Trafficking/degradation biology, not clean first-pass AF3 interface targets |

Detailed partner evidence table:

```text
03_candidate_partner_evidence.tsv
```

## Proposed AF3 Designs

| Design | Residues | Partner | Status |
|---|---:|---|---|
| hERG homotetramer | KCNH2 1-1159 x4 | KCNH2 homomeric copies | top-3 draft JSON ready |
| PAS/PAC-CNBHD regulatory pair | KCNH2 1-150 | KCNH2 720-870 | top-3 draft JSON ready |
| KCNH2-KCNE2 comparator | KCNH2 1-1159 | KCNE2 1-123 | top-3 comparator JSON ready; prior weak |
| Channel-core tetramer sensitivity | KCNH2 398-870 x4 | KCNH2 homomeric copies | held optional |
| KCNH2-KCNE1 auxiliary comparator | KCNH2 1-1159 | KCNE1 1-129 | held |

Detailed design table:

```text
04_af3_design_decisions.tsv
```

## Why The Homotetramer Comes First

KCNH2 is a channel. The most defensible first-pass structural question is whether variants lie near a hERG channel assembly/interface context. The KCNE2 relationship matters biologically, but the first AF3 run did not support a usable interface. Therefore, the systematic packet keeps KCNE2 visible but does not let it dominate the feature strategy.

The PAS/PAC-CNBHD job is included because it is compact and mechanistically interpretable. It asks a different question from the full homotetramer: whether variants fall near a regulatory-domain packing surface inside KCNH2.

## Source Links For Manual Review

- UniProt KCNH2: https://www.uniprot.org/uniprotkb/Q12809/entry
- UniProt KCNE2: https://www.uniprot.org/uniprotkb/Q9Y6J6/entry
- UniProt KCNE1: https://www.uniprot.org/uniprotkb/P15382/entry
- Reactome: https://reactome.org/
- RCSB/PDB: https://www.rcsb.org/
- KCNH2/KCNE association support: https://pubmed.ncbi.nlm.nih.gov/9230439/ ; https://pubmed.ncbi.nlm.nih.gov/10219239/
- KCNE2 IKr support: https://pubmed.ncbi.nlm.nih.gov/10219239/ ; https://pubmed.ncbi.nlm.nih.gov/12185453/
- KCNH2 trafficking/regulatory examples: https://pubmed.ncbi.nlm.nih.gov/25281747/ ; https://pubmed.ncbi.nlm.nih.gov/27916661/ ; https://pubmed.ncbi.nlm.nih.gov/26363003/

Specific Ctrl-F source table:

```text
09_specific_citation_sources.tsv
```

## AF3 Draft Job JSONs

Top-3 upload file:

```text
af3_job_json/AF3_SYS_gene_008_KCNH2_TOP3_upload_alphafoldserver.json
```

Individual JSONs:

```text
af3_job_json/AF3_SYS_G008_KCNH2_HOMOTETRAMER_FULL_LENGTH_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G008_KCNH2_PAS_PAC_1_150_CNBHD_720_870_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G008_KCNH2_KCNE2_FULL_LENGTH_COMPARATOR_DRAFT.alphafoldserver.json
```

Place AlphaFold Server outputs here:

```text
af3_job_results/
```

## What To Watch When AF3 Results Come Back

For KCNH2, promote features only if the homotetramer or PAS/CNBHD job has strong or moderate interface support. If the homotetramer is dominated by disorder or unclear membrane placement, use the held channel-core tetramer sensitivity design next. For KCNE2, require a major improvement over the prior weak result before using any KCNE2 distance/contact features.
