# Gene 005: SCN5A AF3 Systematic Review

Date: 2026-07-06

## Current Decision

SCN5A is gene #5 in the systematic AF3 queue and is one of the strongest AF3 candidates so far. Unlike several large sarcomere genes, SCN5A already has prior full-length AF3 evidence in this project.

Two completed local AlphaFold Server jobs support usable interaction features:

| Prior job | Pair | Pair ipTM | Min PAE | Max contact probability | Interface residues | QC call |
|---|---|---:|---:|---:|---:|---|
| AF3_JOB_001_SCN5A_CALM | SCN5A-CALM1 | 0.77 | 3.6 | 0.86 | 98 | strong supported interface |
| AF3_PILOTV1_REVIEW_SCN5A_SCN1B | SCN5A-SCN1B | 0.82 | 1.8 | 0.98 | 110 | strong supported interface |
| AF3_JOB_001_SCN5A_CALM duplicate run | SCN5A-CALM1 | 0.74 | 4.1 | 0.84 | 98 | strong supported interface |

So this packet does two things:

1. It preserves those two supported contexts as top-priority AF3 feature families.
2. It adds SCN5A-SNTA1 as a third draft job because syntrophin/PDZ scaffold biology is cardiac-relevant, but it is marked weaker until the AF3 output is QC-reviewed.

## Local Dataset Context

| Metric | Value |
|---|---:|
| Included SCN5A rows from priority/ready table | 4,603 |
| Binary rows | 2,123 |
| VUS rows | 2,429 |
| Pathogenic rows | 530 |
| Benign rows | 1,593 |
| Missing label rows in ready model table | 51 |
| Rows with VEP protein position | 4,603 |
| Rows with dbNSFP protein position | 2,603 |
| Rows with older `protein_position` | 89 |

Important mapping decision: SCN5A AF3 feature mapping should use VEP/UniProt Q14524 reconciliation, not the older sparse `protein_position` column.

## Sequence Audit

| Field | Value |
|---|---|
| SCN5A UniProt accession | Q14524 |
| Reviewed entry | yes |
| Canonical protein length | 2,016 aa |
| Baseline sequence policy | reviewed UniProt canonical SCN5A amino-acid sequence |
| PTM/isoform policy | no PTMs, glycosylation states, or alternative isoforms in baseline; sensitivity designs only if justified |
| Full-length AF3 | acceptable here because prior full-length SCN5A jobs worked, but membrane/gating context remains missing |

## Candidate Partner Summary

| Partner | Priority | Why |
|---|---|---|
| CALM1 + 4 Ca2+ | top-3 | Direct UniProt/PubMed-supported SCN5A-CALM interaction and prior strong AF3 interface |
| SCN1B | top-3 | Regulatory beta-1 sodium-channel subunit with direct UniProt/PubMed support and prior strong AF3 interface |
| SNTA1 | top-3 | Sarcolemmal syntrophin scaffold context; useful but weaker because evidence is by-similarity-heavy and needs AF3 QC |
| FGF13 | hold | Strong regulatory partner, but sequence/domain setup not ready and top-3 quota is filled |
| ANK3 | hold | Important targeting/scaffold biology, but exact ANK3 domain boundary is needed |
| PKP2 | hold | Relevant desmosome/channel context, but better handled in PKP2/DSP packet or later sensitivity run |
| NEDD4/NEDD4L/WWP2/GPD1L | hold | Regulatory/trafficking context, not a clean first-pass interface feature |

Detailed evidence table:

```text
03_candidate_partner_evidence.tsv
```

## Proposed AF3 Designs

| Design | Residues | Partner | Status |
|---|---:|---|---|
| SCN5A-CALM1+Ca4 | 1-2016 | CALM1 + 4 Ca2+ | top-3 draft JSON ready; prior AF3 strong |
| SCN5A-SCN1B | 1-2016 | SCN1B | top-3 draft JSON ready; prior AF3 strong |
| SCN5A-SNTA1 | 1-2016 | SNTA1 | top-3 draft JSON ready; needs AF3 output QC |
| SCN5A-FGF13 | boundary not selected | FGF13 | held, not generated |
| SCN5A-ANK3 | boundary needed | ANK3 domain | held, not generated |
| SCN5A-PKP2 | boundary needed | PKP2 | held, not generated |

Detailed design table:

```text
04_af3_design_decisions.tsv
```

## Why Full-Length Is Acceptable Here

SCN5A is long and membrane-bound, but the prior full-length SCN5A-CALM1 and SCN5A-SCN1B AF3 runs were feasible and produced strong pair-level confidence. For this gene, full-length jobs preserve the cytosolic loops, C-terminus, extracellular beta-subunit context, and long-range channel context better than arbitrary fragments.

However, the feature interpretation still has limits. AF3 Server does not include the lipid bilayer, membrane potential, sodium ions through the pore, glycosylation, phosphorylation state, or channel open/inactivated state. Therefore, SCN5A AF3 features should be stored as interaction-context features, not as full channel-gating truth.

Recommended feature flags:

```text
af3_context_level = full_length_no_membrane
af3_full_length_not_modeled_reason = none; full length modeled, but membrane/gating state absent
af3_partner_boundary_status = full_length_partner
af3_promoted_feature = true only after pair QC supports interface
```

## Source Links For Manual Review

- UniProt SCN5A: https://www.uniprot.org/uniprotkb/Q14524/entry
- UniProt CALM1: https://www.uniprot.org/uniprotkb/P0DP23/entry
- UniProt SCN1B: https://www.uniprot.org/uniprotkb/Q07699/entry
- UniProt SNTA1: https://www.uniprot.org/uniprotkb/Q13424/entry
- SCN5A-CALM PubMed link: https://pubmed.ncbi.nlm.nih.gov/21167176/
- SCN5A-SCN1B PubMed link: https://pubmed.ncbi.nlm.nih.gov/21994374/
- Prior AF3 pair QC table: `results/af3/current_completed_af3_summary/completed_af3_pair_qc.tsv`

Specific Ctrl-F source table:

```text
09_specific_citation_sources.tsv
```

## AF3 Draft Job JSONs

Top-3 upload file:

```text
af3_job_json/AF3_SYS_gene_005_SCN5A_TOP3_upload_alphafoldserver.json
```

Individual JSONs:

```text
af3_job_json/AF3_SYS_G005_SCN5A_CALM1_CA4_FULL_LENGTH_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G005_SCN5A_SCN1B_FULL_LENGTH_DRAFT.alphafoldserver.json
af3_job_json/AF3_SYS_G005_SCN5A_SNTA1_FULL_LENGTH_DRAFT.alphafoldserver.json
```

Place AlphaFold Server outputs here:

```text
af3_job_results/
```

## What To Watch When AF3 Results Come Back

CALM1 and SCN1B already have strong prior local AF3 support, so their reruns are mainly for systematic record-keeping or updated output consistency. SNTA1 should not be promoted automatically. It should require pair-level support, low interface PAE, credible contact probability, and a biologically plausible SCN5A cytosolic/PDZ-region contact before becoming a model feature.

## Next Action For SCN5A

Upload only the top-3 JSON unless you explicitly want sensitivity jobs. After output is added, run pair-level ipTM/PAE/contact QC and map interface-distance/contact features to UniProt Q14524 positions using the VEP protein-position field.
