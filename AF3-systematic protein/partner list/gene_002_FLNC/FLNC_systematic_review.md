# Gene 002: FLNC AF3 Systematic Review

Date: 2026-07-06

## Current Decision

FLNC is gene #2 in the systematic AF3 queue and one of the highest-burden genes in the current modeling matrix. It should **not** be submitted as one full-length multi-partner AF3 job.

Reason:

```text
FLNC canonical protein length = 2,725 aa
```

FLNC is a modular actin/Z-disc scaffold with an N-terminal actin-binding domain, 24 filamin repeats, hinge regions, and a C-terminal self-association/dimerization region. The reviewer-safe design is therefore domain-level AF3 jobs, not one giant complex.

## Local Dataset Context

| Metric | Value |
|---|---:|
| Included FLNC rows from priority file | 5671 |
| Current matrix FLNC rows | 5671 |
| Binary rows | 2969 |
| VUS rows | 2702 |
| Pathogenic rows | 608 |
| Benign rows | 2361 |
| Current rows with `protein_position` | 6 |
| Current rows with VEP protein position | 5670 |
| Current rows with dbNSFP protein position raw | 2811 |

Important caveat: the older `protein_position` column only covers 6 FLNC rows, but VEP protein positions are present for 5,670 rows. So the AF3 feature mapper should use and QC the VEP-derived protein positions for FLNC, then reconcile them against UniProt Q14315 before variant-level AF3 features are claimed.

## Sequence Audit

| Field | Value |
|---|---|
| UniProt accession | Q14315 |
| Reviewed entry | yes |
| Canonical protein length | 2,725 aa |
| Baseline sequence policy | reviewed UniProt canonical amino-acid sequence |
| PTM/isoform policy | no PTMs or alternative isoforms in baseline; separate sensitivity jobs only if justified |
| Full-length AF3 | not recommended as first design |

## Candidate Partner Summary

| Partner | Priority | Why |
|---|---|---|
| FLNC homodimer/tail | primary | C-terminal self-association/dimerization is compact, disease-linked, and structurally interpretable |
| ACTC1/actin | primary candidate | FLNC has an N-terminal actin-binding domain and cardiac actin context is biologically relevant |
| XIRP1 | primary candidate | UniProt maps XIRP1 interaction to FLNC repeat 20 |
| PGM5 | primary candidate | UniProt maps PGM5 interaction to FLNC repeats 18-21; PGM5/Xin/FLNC Z-disc biology is relevant |
| INPPL1/SHIP2 | secondary | Curated FLNC interaction region exists, but cardiac specificity is weaker than actin/Z-disc/dimer jobs |
| MYOT | deferred | Supported Z-disc partner, but exact AF3 fragment boundary needs review |
| ACTN2 | deferred | Relevant Z-disc context, but direct AF3 design was not promoted from general co-localization alone |

Detailed partner evidence table:

```text
03_candidate_partner_evidence.tsv
```

## Proposed AF3 Designs

| Design | Residues | Partner | Status |
|---|---:|---|---|
| ABD | 1-259 | ACTC1 | draft JSON ready, simplified actin monomer |
| Tail homodimer | 2593-2725 | FLNC copy | primary candidate draft JSON ready |
| Repeat20-XIRP1 | 2244-2306 | XIRP1 | draft JSON ready, partner boundary needs review |
| R18-21-PGM5 | 1947-2401 | PGM5 | draft JSON ready, PGM5 boundary needs review |
| R22-24-INPPL1 | 2403-2724 | INPPL1 | secondary draft JSON ready |
| R24 compact homodimer | 2630-2724 | FLNC copy | sensitivity draft JSON ready |
| MYOT Z-disc context | boundary needed | MYOT | not JSON-ready yet |
| ACTN2 Z-disc context | boundary needed | ACTN2 | not JSON-ready yet |

Detailed design table:

```text
04_af3_design_decisions.tsv
```

## Source Links For Manual Review

- UniProt FLNC: https://www.uniprot.org/uniprotkb/Q14315/entry
- Human Protein Atlas FLNC: https://www.proteinatlas.org/ENSG00000128591-FLNC
- FLNC Z-disc review: https://pmc.ncbi.nlm.nih.gov/articles/PMC7216277/
- FLNC-actin cardiac paper: https://pubmed.ncbi.nlm.nih.gov/37492967/
- PGM5/Xin/FLNC Z-disc paper: https://pubmed.ncbi.nlm.nih.gov/24963132/

Specific source/citation table for Ctrl-F review:

```text
09_specific_citation_sources.tsv
```

## AF3 Draft Job JSONs

Combined upload candidate:

```text
af3_job_json/AF3_SYS_gene_002_FLNC_all_draft_jobs_alphafoldserver.json
```

Manifest:

```text
af3_job_json/AF3_SYS_gene_002_FLNC_json_manifest.tsv
```

Jobs intentionally not generated yet:

```text
af3_job_json/AF3_SYS_gene_002_FLNC_not_generated_jobs.tsv
```

| AF3 Job | Entities | Upload Status | Evidence IDs |
|---|---|---|---|
| `AF3_SYS_G002_FLNC_ABD_1_259_ACTC1_DRAFT` | FLNC 1-259 + ACTC1 | draft review needed | `SRC_UNIPROT_FLNC_ACTIN_BINDING`, `SRC_FLNC_CARDIAC_ACTIN_PMID37492967`, `SRC_UNIPROT_ACTC1` |
| `AF3_SYS_G002_FLNC_TAIL_2593_2725_HOMODIMER_DRAFT` | FLNC 2593-2725 x2 | draft review needed | `SRC_UNIPROT_FLNC_SELF_ASSOCIATION_TAIL`, `SRC_UNIPROT_FLNC_DIMER_PMID12525170` |
| `AF3_SYS_G002_FLNC_R20_2244_2306_XIRP1_DRAFT` | FLNC 2244-2306 + XIRP1 full length | draft review needed | `SRC_UNIPROT_FLNC_XIRP1_R20`, `SRC_UNIPROT_XIRP1` |
| `AF3_SYS_G002_FLNC_R18_21_1947_2401_PGM5_DRAFT` | FLNC 1947-2401 + PGM5 full length | draft review needed | `SRC_UNIPROT_FLNC_PGM5_R18_21`, `SRC_UNIPROT_PGM5` |
| `AF3_SYS_G002_FLNC_R22_24_2403_2724_INPPL1_DRAFT` | FLNC 2403-2724 + INPPL1 full length | secondary draft review needed | `SRC_UNIPROT_FLNC_INPPL1_R22_24`, `SRC_UNIPROT_INPPL1` |
| `AF3_SYS_G002_FLNC_R24_2630_2724_HOMODIMER_COMPACT_DRAFT` | FLNC 2630-2724 x2 | sensitivity draft review needed | `SRC_UNIPROT_FLNC_SELF_ASSOCIATION_TAIL` |

## What To Watch When AF3 Results Come Back

For FLNC, a useful AF3 result should show a plausible interface in the expected FLNC fragment. Stronger candidates should have high pair ipTM, low interface PAE, and residue contacts close to the curated FLNC interaction region. Jobs using full XIRP1, PGM5, or INPPL1 should be treated as exploratory until the partner-side binding fragment is narrowed.

## Next Action For FLNC

1. Upload the combined JSON or individual JSONs to AlphaFold Server after manual review.
2. Place server output folders under:

```text
af3_job_results/
```

3. After output is added, run the same QC path used for RYR2: pair-level ipTM/PAE/contact review, chain-level pLDDT, and feature-decision table.
4. Build the FLNC AF3 feature mapper from VEP-derived protein positions first, then reconcile against UniProt Q14315 and flag any transcript/isoform mismatches.
