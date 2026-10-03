# Gene 010: HCN4 AF3 Systematic Review

Date: 2026-07-06

## Current Decision

HCN4 is gene #10 in the systematic AF3 queue. It has 2,170 included variant rows, 930 binary rows, and 1,240 VUS rows in the current rescue-aware modeling table. It is a high-burden but sparse-pathogenic cardiogenetics gene: only 20 pathogenic rows are present, while most labeled rows are benign or VUS.

The first-pass packet uses three jobs:

1. Full-length HCN4 homotetramer.
2. HCN4 C-linker/CNBD tetramer with cAMP ligand draft.
3. HCN4 compact channel-core/C-linker/CNBD tetramer.

This is intentionally a channel-assembly/regulatory packet, not a broad protein-partner packet. HCN4 has weaker or less cardiac-specific interaction candidates, but the strongest publication-safe biological context is HCN4 channel tetramerization and cyclic nucleotide regulation.

## Local Dataset Context

| Metric | Value |
|---|---:|
| Included HCN4 rows | 2,170 |
| Binary rows | 930 |
| VUS rows | 1,240 |
| Pathogenic rows | 20 |
| Benign rows | 910 |
| Missing-label rows | 0 |
| HiRO-overlap rows | present in HiRO source records, mainly VUS |

HCN4 should be treated as a sparse-P/LP gene in downstream subgroup analysis. AF3 features may be most useful for VUS prioritization and mechanism review, not for gene-specific supervised learning.

## Sequence Audit

| Field | Value |
|---|---|
| HCN4 UniProt accession | Q9Y3Q4 |
| Reviewed entry | yes |
| Canonical protein length | 1,203 aa |
| Protein name | Potassium/sodium hyperpolarization-activated cyclic nucleotide-gated channel 4 |
| Baseline sequence policy | reviewed UniProt canonical amino-acid sequence |
| PTM/isoform policy | no membrane, voltage-state, phosphorylation, or S-palmitoylation in baseline |
| Ligand policy | cAMP/CNBD job is draft and must pass server/QC review |

## Candidate Partner Summary

| Partner | Priority | Why |
|---|---|---|
| HCN4 homotetramer | top-3 lead | Curated channel assembly state and strongest direct AF3 context |
| HCN4 C-linker/CNBD + cAMP | top-3 | Direct cyclic nucleotide regulation; strong IntAct/PDB evidence |
| HCN4 channel-core tetramer | top-3 | Smaller fallback if full-length tetramer is rejected or messy |
| HCN1/HCN2/HCN3 heteromers | hold | Possible but less cardiac-specific for first-pass HCN4 model |
| FYN/ABL1 | hold | Peptide-array/signaling evidence, not primary channel-complex context |
| PEX5L/IRAG1/IRAG2 | hold | UniProt by-similarity notes need targeted manual review |

Detailed partner evidence table:

```text
03_candidate_partner_evidence.tsv
```

## Proposed AF3 Designs

| Design | Residues | Partner | Status |
|---|---:|---|---|
| HCN4 full-length homotetramer | HCN4 1-1203 x4 | HCN4 x4 | top-3 draft JSON ready, near size ceiling |
| HCN4 C-linker/CNBD+cAMP tetramer | HCN4 518-724 x4 | 4 cAMP ligands as CCD_CMP | top-3 draft JSON ready, ligand QC needed |
| HCN4 channel-core/C-linker/CNBD tetramer | HCN4 209-724 x4 | HCN4 x4 | top-3 draft JSON ready |
| HCN4-HCN1/2/3 heteromer | full-length/domain TBD | HCN1/HCN2/HCN3 | held |
| HCN4 kinase/signaling association | boundary needed | FYN/ABL1 | held |

Detailed design table:

```text
04_af3_design_decisions.tsv
```

## Why The Homotetramer Comes First

HCN4 is an ion channel, and UniProt/IntAct support HCN4 homotetramerization. For a cardiogenetics classifier, this is cleaner than broad interaction mining because the homotetramer directly matches the cardiac pacemaker channel mechanism. The full-length tetramer is large, so a compact channel-core/CNBD tetramer is included as a practical fallback.

The cAMP job is biologically important because HCN4 is cyclic nucleotide gated and the CNBD is a major regulatory region. However, this job has an extra technical caveat: the server must accept `CCD_CMP` as cAMP. If not, generate a protein-only CNBD tetramer and document the ligand failure.

## Source Links For Manual Review

- UniProt HCN4: https://www.uniprot.org/uniprotkb/Q9Y3Q4/entry
- Reactome HCN channels: https://reactome.org/content/detail/R-HSA-1296061
- RCSB 4HBN: https://www.rcsb.org/structure/4HBN
- HPA HCN4: https://www.proteinatlas.org/HCN4
- HCN4 human heart channel characterization: https://pubmed.ncbi.nlm.nih.gov/10228147/
- HCN channel heteromerization: https://pubmed.ncbi.nlm.nih.gov/12928435/
- HCN4 cAMP-dependent gating: https://pubmed.ncbi.nlm.nih.gov/20829353/
- HCN4 tetramerization dynamics: https://pubmed.ncbi.nlm.nih.gov/22006928/
- HCN4 cyclic dinucleotide/cAMP pocket context: https://pubmed.ncbi.nlm.nih.gov/24776929/
- cAMP/CNBD ligand binding work: https://pubmed.ncbi.nlm.nih.gov/24605759/

Specific Ctrl-F source table:

```text
09_specific_citation_sources.tsv
```

## AF3 Draft Job JSONs

Top-3 upload file:

```text
af3_job_json/AF3_SYS_gene_010_HCN4_TOP3_upload_alphafoldserver.json
```

Individual JSONs:

```text
af3_job_json/AF3_SYS_G010_HCN4_HOMOTETRAMER_FULL_LENGTH_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G010_HCN4_CLINKER_CNBD_518_724_CAMP4_TETRAMER_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G010_HCN4_CORE_209_724_TETRAMER_DRAFT.alphafoldserver.json
```

Place AlphaFold Server outputs here:

```text
af3_job_results/
```

## What To Watch When AF3 Results Come Back

For HCN4, do not promote the full-length tetramer just because it is biologically complete. Use pairwise ipTM, inter-chain PAE, contact confidence, and variant proximity to the tetramer interface, pore/C-linker, and CNBD/cAMP region. If the full-length job is weak, the compact channel-core and CNBD jobs may be better feature sources. If the cAMP ligand is rejected or poorly placed, keep ligand features missing and use protein-only CNBD geometry.
