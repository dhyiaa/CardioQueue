# Gene 001: RYR2 AF3 Systematic Review

Date: 2026-07-04

## Current Decision

RYR2 is the highest-burden gene in the included matrix, but it should **not** be submitted as a full-length AlphaFold 3 complex first.

Reason:

```text
RYR2 canonical protein length = 4,967 aa
```

A full-length RYR2 tetramer or full RYR2 plus partner complex is too large and not the cleanest biological design for AlphaFold Server. RYR2 should be handled as domain-level or structure-informed jobs.

## Local Dataset Context

| Metric | Value |
|---|---:|
| Included RYR2 rows | 9,465 |
| Binary rows | 4,255 |
| VUS rows | 5,193 |
| Pathogenic rows | 261 |
| Benign rows | 3,994 |
| Current rows with protein position | 18 |

Important caveat: only 18 current RYR2 rows have direct `protein_position` in the ready matrix. Before AF3 feature extraction is useful for RYR2, we need a better RYR2 protein-position mapping pass from VEP/dbNSFP/HGVS.

## Sequence Audit

| Field | Value |
|---|---|
| UniProt accession | Q92736 |
| Reviewed entry | yes |
| Canonical protein length | 4,967 aa |
| Baseline full-length AF3 | not recommended |
| Recommended strategy | domain-level or experimentally structure-informed jobs |

Baseline AF3 sequence should remain the reviewed UniProt canonical RYR2 sequence unless a domain/PTM/isoform-specific sensitivity job is explicitly defined.

## Candidate Partner Summary

| Partner | Priority | Why |
|---|---|---|
| RYR2 homomer/tetramer | biological context, not first AF3 job | RYR2 functions as a tetrameric SR calcium-release channel |
| FKBP1B | high priority | Direct regulatory partner with UniProt/IntAct/Reactome/RCSB support |
| CALM1/CALM2/CALM3 | high priority | Calmodulin regulation is central and structure-supported |
| CASQ2/TRDN | secondary domain context | Luminal calcium-release complex context, but RYR2 domain boundary must be defined |
| Ca2+ | conditional | Include only when specific Ca/CaM or calcium-release design requires it |
| S100A1 | deferred | Plausible regulator, needs citation/evidence review before submission |
| PKA/AKAP6/phosphatases | deferred | PTM/regulatory-state biology, not baseline unmodified AF3 |

Detailed partner evidence table:

```text
03_candidate_partner_evidence.tsv
```

## Proposed AF3 Designs

| Design | Status | Decision |
|---|---|---|
| RYR2 domain + FKBP1B | structure contact regions identified; design choice still needed | primary |
| RYR2 CaM-binding domain + CALM1 + optional Ca2+ | domain boundary needed | primary |
| CASQ2 + TRDN +/- RYR2 luminal domain | candidate RYR2 LL1 boundary identified | secondary |
| Full RYR2 tetramer | defer | not practical for AF3 Server |

Detailed design table:

```text
04_af3_design_decisions.tsv
```

## Source Links For Manual Review

- UniProt RYR2: https://www.uniprot.org/uniprotkb/Q92736/entry
- Reactome RYR2 entity: https://reactome.org/content/detail/R-HSA-937283
- Reactome RYR calcium-release reaction: https://reactome.org/content/detail/R-HSA-2855020
- NCBI Gene RYR2: https://www.ncbi.nlm.nih.gov/gene/6262
- RCSB 7U9T: https://www.rcsb.org/structure/7U9T
- CaM-bound RyR2 CaMBD structure article: https://pmc.ncbi.nlm.nih.gov/articles/PMC8211408/

Specific source/citation table for Ctrl-F review:

```text
09_specific_citation_sources.tsv
```

Rule for this systematic AF3 branch:

```text
Every gene packet must include source IDs, specific evidence claims, URLs/PMIDs/PDB IDs,
and Ctrl-F terms so the evidence can be manually checked.
```

## Domain Boundary Review

RYR2 can now move forward for calmodulin-domain jobs, but not for full-length or FKBP1B jobs yet.

| Domain design | RYR2 residues | Partner | Status |
|---|---|---|---|
| CaMBD1 | 1940-1965 | CALM1 +/- Ca2+ | candidate |
| CaMBD2 | 3580-3611 | CALM1 +/- Ca2+ | primary candidate |
| CaMBD3 | 4246-4275 | CALM1 + Ca2+ | primary candidate, structure-supported |
| FKBP1B binding context | boundary needed | FKBP1B | not JSON-ready yet |
| Luminal CASQ2/TRDN context | boundary needed | CASQ2/TRDN | not JSON-ready yet |

Detailed boundary table:

```text
06_domain_boundary_review.tsv
```

Recommended next RYR2 action: make draft AF3 JSON for CaMBD2-CALM1-Ca2+ and CaMBD3-CALM1-Ca2+ after extracting the exact UniProt canonical fragments and adding flanking residues if needed. Keep FKBP1B and CASQ2/TRDN RYR2-domain jobs in boundary-review status.

## Fixing The Two Boundary-Needed Rows

| Design | Residues | Partner | Updated Status |
|---|---:|---|---|
| FKBP1B context | structure-derived RYR2 contact residues across ~633-750, 1309-1313, 1631-1781; broader design candidates include 599-809, 1300-1320, 1620-1790 or one continuous 599-1790 fragment | FKBP1B | evidence-constrained, but still needs AF3 design choice |
| CASQ2/TRDN luminal context | RYR2 first luminal loop candidate 4521-4573; TRDN KEKE/CASQ2 motif around 200-224 needs topology review | CASQ2/TRDN | candidate boundary identified, draft pair jobs possible |

Detailed resolution table:

```text
08_boundary_needed_resolution.tsv
```

### FKBP1B Context

The missing evidence was not whether FKBP1B interacts with RYR2. That is supported. The missing piece was the RYR2 sequence design.

I downloaded and parsed the RCSB 7U9T mmCIF file:

```text
structure_7u9t.cif
07_7U9T_RYR2_FKBP1B_contact_residues_5A.tsv
```

The structure contains RYR2 chains `B/G/H/I`, FKBP1B chains `A/C/D/E`, and calmodulin chains `F/J/K/L`. RYR2 residues within 5A of FKBP1B are not one clean peptide. They cluster across multiple RYR2 regions:

```text
633-750 region, near/within B30.2/SPRY1
1309-1313 region
1631-1781 region
```

Therefore, FKBP1B is not yet a simple JSON job. We need to choose one of three scientifically defensible routes:

| Route | Meaning | Pros | Risk |
|---|---|---|---|
| Large continuous fragment | RYR2 599-1790 + FKBP1B | Keeps native sequence continuity across the contact context | Large and may still be hard for AF3 Server |
| Separate fragment/domain entities | RYR2 599-809, 1300-1320, 1620-1790 + FKBP1B | Captures actual contact regions | Artificial disconnected fragments may be biologically awkward |
| PDB-derived features only | Use 7U9T contact mapping, no AF3 job | Most structurally grounded | Not an AF3-generated feature |

My recommendation is to defer FKBP1B JSON until we decide between continuous-fragment AF3 and PDB-derived interface features.

### CASQ2/TRDN Luminal Context

This row is now closer to JSON-ready.

The missing RYR2 boundary can be addressed using the reported RYR2 first luminal loop:

```text
RYR2 LL1: residues 4521-4573
```

The safer AF3 plan is not to jump immediately to a 3-chain complex. Instead:

| Candidate Job | Status | Why |
|---|---|---|
| RYR2 LL1 4521-4573 + CASQ2 | draftable after fragment extraction | Direct RYR2-CASQ2 luminal-loop hypothesis |
| CASQ2 + TRDN luminal/KEKE region | draftable after TRDN topology review | Tests CASQ2-triadin interaction separately |
| RYR2 LL1 + CASQ2 + TRDN | secondary | Only after pair jobs look plausible |

This is more reviewer-safe than forcing one large luminal complex immediately.

## Next Action For RYR2

Do not create final AF3 JSON yet. First, choose domain boundaries for:

1. Extract RYR2 CaMBD2 and CaMBD3 canonical fragments and create draft AF3 JSONs for CALM1/Ca2+ review.
2. Decide FKBP1B strategy: continuous 599-1790 fragment, separate fragment/domain job, or PDB-derived interface features only.
3. Extract RYR2 LL1 4521-4573 and review TRDN topology/KEKE motif before making CASQ2/TRDN draft jobs.

Once domain boundaries are selected, create server JSON(s) in this gene folder.

## AF3 Draft Job JSONs

Draft AlphaFold Server JSONs were generated for the RYR2 designs that are technically draftable now.

Combined upload candidate:

```text
af3_job_json/AF3_SYS_gene_001_RYR2_all_draft_jobs_alphafoldserver.json
```

Manifest:

```text
af3_job_json/AF3_SYS_gene_001_RYR2_json_manifest.tsv
```

Jobs intentionally not generated yet:

```text
af3_job_json/AF3_SYS_gene_001_RYR2_not_generated_jobs.tsv
```

| AF3 Job | Entities | Upload Status | Evidence IDs |
|---|---|---|---|
| `AF3_SYS_G001_RYR2_CaMBD1_1940_1965_CALM1_CA4_DRAFT` | RYR2 1940-1965 + CALM1 + 4 Ca2+ | draft_review_needed | `SRC_RYR2_CAMBD_REVIEW`, `SRC_UNIPROT_RYR2`, `SRC_UNIPROT_CALM1` |
| `AF3_SYS_G001_RYR2_CaMBD2_3580_3611_CALM1_CA4_DRAFT` | RYR2 3580-3611 + CALM1 + 4 Ca2+ | draft_review_needed | `SRC_UNIPROT_RYR2_CALM_REGION`, `SRC_RYR2_CAMBD_REVIEW` |
| `AF3_SYS_G001_RYR2_CaMBD3_4246_4275_CALM1_CA5_DRAFT` | RYR2 4246-4275 + CALM1 + 5 Ca2+ | draft_review_needed | `SRC_RCSB_7KL5`, `SRC_RYR2_CAMBD_REVIEW` |
| `AF3_SYS_G001_RYR2_LL1_4521_4573_CASQ2_DRAFT` | RYR2 4521-4573 + CASQ2 | draft_review_needed | `SRC_CSQ2_RYR2_LL1_PMID27609834` |
| `AF3_SYS_G001_CASQ2_TRDN_FULL_LENGTH_DRAFT` | CASQ2 + TRDN | draft_review_needed | `SRC_CSQ2_TRDN_PMID27609834` |

Not generated yet:

| Job | Reason |
|---|---|
| `AF3_SYS_G001_RYR2_FKBP1B_CONTEXT` | 7U9T shows a multi-region RYR2 interface; need choose continuous fragment, separate fragments, or PDB-derived features |
| `AF3_SYS_G001_RYR2_LL1_CASQ2_TRDN_THREE_CHAIN` | Should wait until RYR2_LL1-CASQ2 and CASQ2-TRDN pair jobs are reviewed |

## AF3 Server Results Folder

Dheyaa added the returned AlphaFold Server outputs here:

```text
af3_job_results/folds_2026_07_06_02_03/
```

Future gene folders should always include:

```text
af3_job_results/
```

This folder is reserved for raw downloaded AlphaFold Server outputs. Keep JSON inputs in `af3_job_json/` and keep parsed/processed summaries separate.

## AF3 Result QC Summary

The first RYR2 AF3 result batch was reviewed here:

```text
af3_result_qc/AF3_RESULT_REVIEW_gene_001_RYR2.md
```

Main result:

| Job | Main pair | Support | Feature decision |
|---|---|---|---|
| `RYR2_CaMBD2_3580_3611_CALM1_CA4` | RYR2 fragment-CALM1 | strong | promote to trusted AF3Draft interface feature |
| `RYR2_CaMBD1_1940_1965_CALM1_CA4` | RYR2 fragment-CALM1 | moderate | promote with sensitivity flag |
| `RYR2_CaMBD3_4246_4275_CALM1_CA5` | RYR2 fragment-CALM1 | borderline but below threshold | do not promote yet; review/redesign |
| `RYR2_LL1_4521_4573_CASQ2` | RYR2 LL1-CASQ2 | unsupported | do not promote |
| `CASQ2_TRDN_FULL_LENGTH` | CASQ2-TRDN | unsupported | do not promote |

Interpretation:

```text
RYR2-CALM biology is supported by the AF3 systematic run, especially CaMBD2.
The first CASQ2/TRDN luminal designs did not produce trusted interfaces and need redesign.
```
