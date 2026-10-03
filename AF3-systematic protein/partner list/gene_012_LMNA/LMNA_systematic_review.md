# Gene 012: LMNA AF3 Systematic Review

Date: 2026-07-06

## Current Decision

LMNA is gene #12 in the systematic AF3 queue. It has 2,073 included variant rows, 1,153 binary rows, and 912 VUS rows in the current rescue-aware modeling table. It is a primary cardiogenetics gene for cardiolaminopathy, but its biology is not a simple soluble protein-partner system. Lamin A/C forms nuclear lamina assemblies, has proteoform processing, and interacts with inner nuclear membrane proteins.

The first-pass packet uses three jobs:

1. Mature Lamin-A/C 1-646 homodimer.
2. LMNA tail 384-646 plus emerin/EMD nucleoplasmic region 1-222.
3. LMNA IF rod 31-387 homodimer.

This packet does not repeat the old full-length LMNA-EMD job by default. That previous AF3 output was not interface-usable, so the new packet is a compact rescue plus lamin self-assembly packet.

## Local Dataset Context

| Metric | Value |
|---|---:|
| Included LMNA rows | 2,073 |
| Binary rows | 1,153 |
| VUS rows | 912 |
| Pathogenic rows | 529 |
| Benign rows | 624 |

LMNA is model-relevant and biologically important. Later AF3 features should focus on whether a variant lands near a confident lamin self-assembly interface or a confident LMNA-emerin compact interface. Full-length weak LMNA-EMD predictions should not be converted into positive interaction features.

## Sequence Audit

| Field | Value |
|---|---|
| LMNA UniProt accession | P02545 |
| Reviewed entry | yes |
| Canonical protein length | 664 aa |
| Mature Lamin-A/C boundary used | 1-646 |
| LMNA IF rod boundary used | 31-387 |
| LMNA tail boundary used | 384-646 |
| EMD UniProt accession | P50402 |
| EMD nucleoplasmic boundary used | 1-222 |
| Baseline sequence policy | reviewed UniProt canonical amino-acid sequence, with documented mature/domain slices |
| PTM/proteoform policy | no farnesylation, proteolytic processing chemistry, phosphorylation, methylation, membrane, or isoform-specific state in baseline |

The canonical LMNA sequence is prelamin-A/C. For the default assembly job, I used the mature Lamin-A/C 1-646 boundary to avoid modeling the prelamin-A C-terminal processing/farnesylation region as if it were a normal unmodified tail.

## Candidate Partner Summary

| Partner | Priority | Why |
|---|---|---|
| LMNA mature homodimer | top-3 | Core lamin assembly biology with strong UniProt/literature support |
| EMD nucleoplasmic 1-222 | top-3 compact rescue | Direct curated nuclear-envelope partner; prior full-length LMNA-EMD was weak |
| LMNA IF rod homodimer | top-3 | Structured rod/polymerization proxy for variants in the major IF domain |
| Full-length LMNA-EMD | hold | Already tested and not interface-usable |
| TMEM43-EMD-LMNA | hold | Already tested and not interface-usable; membrane/topology context likely missing |
| SUN1/SUN2/SYNE2/LEMD2/TMPO/BANF1/NARF/TMEM201 | hold | Real candidates, but need separate cardiac-context and boundary review |

Detailed partner evidence table:

```text
03_candidate_partner_evidence.tsv
```

## Proposed AF3 Designs

| Design | Residues | Partner | Status |
|---|---:|---|---|
| LMNA mature homodimer | LMNA 1-646 x2 | LMNA self | top-3 draft JSON ready |
| LMNA tail-EMD compact rescue | LMNA 384-646 | EMD 1-222 | top-3 draft JSON ready |
| LMNA IF rod homodimer | LMNA 31-387 x2 | LMNA self | top-3 draft JSON ready |
| Full-length LMNA-EMD repeat | LMNA 1-664 | EMD 1-254 | held, already weak |
| TMEM43-EMD-LMNA trimer | full-length chains | TMEM43 + EMD | held, already weak |

Detailed design table:

```text
04_af3_design_decisions.tsv
```

## Why These Jobs Come First

LMNA is a nuclear lamina protein. The most defensible first question is whether variants sit near lamin assembly regions that AF3 can model with useful confidence. That is why the mature LMNA homodimer and IF rod homodimer are prioritized.

The LMNA-EMD interaction is biologically important, but our previous full-length LMNA-EMD job was weak: pair ipTM 0.10, minimum inter-chain PAE 30.33, and maximum contact probability 0.05. Therefore, the new LMNA-EMD design uses only LMNA tail 384-646 and EMD nucleoplasmic 1-222. This does not prove the real cell uses that isolated pair; it is a testable compact rescue design. If it is still weak, we should say LMNA-EMD did not yield usable AF3 interface features in this setup.

## Source Links For Manual Review

- UniProt LMNA: https://www.uniprot.org/uniprotkb/P02545/entry
- UniProt EMD: https://www.uniprot.org/uniprotkb/P50402/entry
- UniProt TMEM43: https://www.uniprot.org/uniprotkb/Q9BTV4/entry
- HPA LMNA: https://www.proteinatlas.org/LMNA
- HPA EMD: https://www.proteinatlas.org/EMD
- LMNA homodimer/polymerization: https://pubmed.ncbi.nlm.nih.gov/15476822/
- LMNA lamina assembly/protofilament reference: https://pubmed.ncbi.nlm.nih.gov/31434876/
- LMNA homodimer reference: https://pubmed.ncbi.nlm.nih.gov/33706103/
- LMNA-emerin interaction reference: https://pubmed.ncbi.nlm.nih.gov/12475961/
- EMD/prelamin localization reference: https://pubmed.ncbi.nlm.nih.gov/19323649/
- TMEM43/emerin reference: https://pubmed.ncbi.nlm.nih.gov/18230648/
- Local prior AF3 LMNA-EMD QC: results/af3/current_completed_af3_summary/completed_af3_pair_qc.tsv

Specific Ctrl-F source table:

```text
09_specific_citation_sources.tsv
```

## AF3 Draft Job JSONs

Top-3 upload file:

```text
af3_job_json/AF3_SYS_gene_012_LMNA_TOP3_upload_alphafoldserver.json
```

Individual JSONs:

```text
af3_job_json/AF3_SYS_G012_LMNA_MATURE_1_646_HOMODIMER_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G012_LMNA_TAIL_384_646_EMD_NUC_1_222_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G012_LMNA_IF_ROD_31_387_HOMODIMER_DRAFT.alphafoldserver.json
```

Place AlphaFold Server outputs here:

```text
af3_job_results/
```

## What To Watch When AF3 Results Come Back

For LMNA, the useful signal is not global ranking alone. We need pairwise interface evidence: ipTM, inter-chain PAE, contact probability, and whether mapped LMNA variants sit near a confident interface. If the compact LMNA-EMD job remains weak, then LMNA AF3 features should come mainly from lamin assembly designs or be marked unavailable for EMD-interface context.
