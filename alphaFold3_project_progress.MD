# AlphaFold 3 Project Progress

## Current Decision

This file starts a new project branch for an **AlphaFold 3 augmented version**
of the cardiogenetics CatBoost model. It does not replace or rewrite the
previous primary model. The current primary model remains:

```text
Binary CatBoost P/LP vs B/LB + 0.1/0.9 deferral zone
```

The AlphaFold 3 branch will create a new feature family:

```text
AF3 interaction-aware structural features
```

## Safety Boundary

The AF3 work is an additive branch only. It must not overwrite, delete, or
rename the existing modeling matrices, source registries, trained model outputs,
or previous protein-structure feature tables.

Current rule:

```text
Existing model/data branch = frozen unless explicitly changed for a separate task
AF3 branch = new files under datasets/feature_sources/alphafold3/ and results/af3/
```

AF3 input sequences are **not** taken from AlphaFold 2 structures. AF2/UniProt
files may only be used as local reference material for mapping and comparison.
For AF3 submissions, the intended input source is:

```text
Reviewed human canonical UniProt protein sequence + curated partner/ion/ligand entities
```

This keeps the new AF3 layer separate from the previous AlphaFold 2-derived
pLDDT/DSSP/FreeSASA feature work.

Important sequence/PTM note:

```text
Baseline AF3 jobs currently use reviewed UniProt canonical amino-acid sequences.
They do not automatically include all mature-chain processing, alternative isoforms,
or post-translational modifications.
```

This is intentional for the first reproducible AF3 feature layer. PTM-aware,
alternative-isoform, or mature-chain AF3 jobs should be separate targeted
sensitivity/mechanism jobs with explicit evidence and labels.

The goal is to test whether these features improve or explain variant
prioritization beyond the existing model, especially for missense variants,
VUS, and cases where stability predictors do not explain pathogenicity.

## Core Idea

Do **not** run AlphaFold 3 once per variant.

That would be too slow, too expensive in quota, and not easy to defend for an
85,677-row variant matrix.

Instead:

```text
1 AlphaFold 3 job = 1 biologically justified wild-type molecular complex
```

Then map many variants onto that complex.

Example:

```text
KCNQ1 + KCNE1 + CALM1 + Ca2+
```

One AF3 complex prediction can then help annotate many KCNQ1, KCNE1, and
calmodulin variants by asking:

- Is the variant residue at the KCNQ1-KCNE1 interface?
- Is the variant residue at the KCNQ1-calmodulin interface?
- How close is the residue to Ca2+?
- Is the local structure confident?
- Is the predicted interface reliable by PAE/ipTM?

This is different from FoldX. FoldX asks whether a mutation destabilizes the
folded protein. AF3 complex mapping asks whether a residue sits at a biologic
interaction surface, ligand pocket, ion-binding region, or partner interface.

## Why We Are Doing This

The current model is already strong, but the protein-structure feature family is
sparse and mostly stability-focused.

Current manuscript-facing matrix:

| Item | Count |
|---|---:|
| Analysis-ready variant rows | `85,677` |
| Primary binary training rows | `42,990` |
| Exploratory three-class rows | `84,721` |
| FoldX DDG coverage | `313 / 85,677` |
| FoldX DDG coverage fraction | `0.37%` |

This means FoldX is useful when available, but it cannot carry the full model.
AF3 complex mapping may cover more missense residues if we model wild-type
complexes and map variants by protein position.

The scientific reason is simple:

| Mechanism | Existing / New Feature |
|---|---|
| Population constraint and allele frequency | gnomAD/dbNSFP-derived features |
| Sequence/in silico damagingness | REVEL, SIFT, PolyPhen, CADD, AlphaMissense |
| Splice/consequence context | VEP, SpliceAI |
| Protein stability | FoldX DDG |
| Protein domains/local structure | UniProt, AlphaFold, DSSP, FreeSASA |
| Molecular interaction context | **New AF3 complex-mapping features** |

AF3 is most useful for variants that may disrupt binding without strongly
destabilizing the whole protein.

## Gene Panel Summary

The AF3 branch should recover the exact gene list from the local matrix, not
from pasted prose.

Current matrix checked:

```text
datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv
```

Observed gene counts:

| Gene set | Count | Notes |
|---|---:|---|
| Raw `primary_gene` strings in ready matrix | `61` | Includes case/alias artifacts |
| Normalized uppercase genes in ready matrix | `60` | `RyR2` and `RYR2` collapse |
| Raw included gene strings | `54` | `primary_model_inclusion == include` |
| Normalized included genes | `53` | This matches the expected ~52-53 gene panel |

Included normalized genes:

```text
ACTC1, BAG3, CACNA1C, CACNB2, CALM1, CALM2, CALM3, CASQ2, CSRP3,
DES, DSC2, DSG2, DSP, EMD, FHOD3, FLNC, GLA, HCN4, HRAS, JUP,
KCNE1, KCNE2, KCNH2, KCNJ2, KCNQ1, KRAS, LAMP2, LMNA, MAP2K2,
MIB1, MYBPC3, MYH7, MYL2, MYL3, NRAS, PKP2, PLN, PRKAG2,
PTPN11, RAF1, RBM20, RIT1, RYR2, SCN1B, SCN5A, SLC22A5, SNTA1,
TCAP, TMEM43, TNNI3, TNNT2, TPM1, TRDN
```

Scope-specific raw gene summary from the current matrix:

| Modeling scope | Raw gene count | Genes |
|---|---:|---|
| `primary` | `38` | ACTC1, BAG3, CACNA1C, CALM1, CALM2, CALM3, CASQ2, CSRP3, DES, DSC2, DSG2, DSP, FLNC, HCN4, JUP, KCNE1, KCNE2, KCNH2, KCNJ2, KCNQ1, LMNA, MYBPC3, MYH7, MYL2, MYL3, PKP2, PLN, PRKAG2, RBM20, RYR2, RyR2, SCN5A, TCAP, TMEM43, TNNI3, TNNT2, TPM1, TRDN |
| `expanded_scope` | `11` | EMD, GLA, HRAS, KRAS, LAMP2, MAP2K2, NRAS, PTPN11, RAF1, RIT1, SLC22A5 |
| `primary_monitor` | `5` | CACNB2, FHOD3, MIB1, SCN1B, SNTA1 |
| `sensitivity_only` | `3` | ANK2, SOS1, TRPM4 |
| `ttn_separate` | `1` | TTN |
| `exclude_or_quarantine` | `4` | AKAP9, DMD, FPGT-TNNI3K, KCNE1 |

Important QC notes:

- `RYR2` and `RyR2` both appear and must be normalized.
- The historical text mentioned `KNCH2`; the current ready matrix contains
  `KCNH2`, not `KNCH2`. Codex must explicitly check this in the AF3 audit.
- `KCNE1` appears in both included and exclude/quarantine scope in different
  rows. Decision: **include KCNE1 for AF3** because it is biologically central
  to KCNQ1/IKs complex modeling. Preserve the mixed-scope warning for audit
  provenance, but do not block KCNE1 from AF3 partner assignment.
- `TTN` is separate/special. Full-length TTN is not feasible for AF3 server
  input. Any TTN AF3 use must be domain-level or deferred.

## Desired AF3 Coverage Philosophy

The final AF3 plan is:

```text
Full-panel AF3 complex coverage
not
prioritized AF3 pilot only
```

But “full-panel” does **not** mean combining unrelated proteins into AF3 jobs.
It means every gene gets reviewed and assigned one AF3 coverage status.

Allowed statuses:

| AF3 coverage status | Meaning |
|---|---|
| `shared_complex` | Gene is covered by a multi-gene biological complex job |
| `direct_protein_complex` | Gene has one direct protein partner job |
| `protein_ligand_or_ion_context` | Gene is modeled with ligand/ion context |
| `domain_level_complex` | Full-length input is too large or biologically inappropriate, so a domain/fragment is used |
| `no_curated_partner` | No defensible partner found yet |
| `manual_review_needed` | Potential partner exists but needs human review |
| `excluded_or_quarantine` | Gene is outside current primary AF3 model scope |
| `ttn_domain_only_or_defer` | TTN requires special handling |

This prevents selection bias while staying biologically defensible.

## Update 2026-07-03: Sequence And Starter Partner Manifest Built

The next real step was completed:

```text
Build the first partner/sequence manifest, not final AF3 jobs.
```

New script:

```text
datasets/feature_sources/alphafold3/scripts/build_af3_sequence_partner_manifests.py
```

Inputs:

```text
results/af3/gene_panel_audit.tsv
datasets/feature_sources/protein_structure/raw/uniprot_human_reviewed_2026_02.fasta
```

Outputs:

```text
datasets/feature_sources/alphafold3/interim/uniprot_gene_to_sequence.tsv
datasets/feature_sources/alphafold3/interim/af3_all_entities.fasta
datasets/feature_sources/alphafold3/interim/af3_manual_cardiac_complex_whitelist.tsv
results/af3/af3_full_panel_partner_manifest.tsv
results/af3/af3_sequence_partner_manifest.summary.json
```

Important: these are **review manifests**, not final AF3 server jobs. No row is
ready for submission until its evidence and domain/sequence design are checked.

Current manifest summary:

| Item | Count |
|---|---:|
| AF3 target genes for review | `56` |
| Sequence entities in FASTA/TSV | `59` |
| Missing reviewed human UniProt sequences | `0` |
| Manual cardiac whitelist complexes | `10` |
| Partner manifest rows | `172` |
| Targets with starter whitelist context | `31` |
| Targets still needing partner review | `25` |

The AF3 review set is broader than the strict 53-gene primary included modeling
panel because it also keeps sensitivity/monitor genes visible for AF3 planning.
This does **not** mean all 56 will be included in the AF3 model. It means they
are available for review and will receive explicit coverage decisions.

External partner sequences added because they are needed for obvious cardiac
complex contexts:

| External partner | Why included |
|---|---|
| `ATP2A2` | PLN/SERCA2a calcium-handling context |
| `FKBP1B` | RYR2 regulatory complex context |
| `TNNC1` | Cardiac troponin/thin-filament calcium complex |

Starter manual complex whitelist:

| Complex group | Main genes covered | Current status |
|---|---|---|
| `AF3_COMPLEX_KCNQ1_KCNE1_CALM_CA` | KCNQ1, KCNE1, CALM1/2/3 | Review pending |
| `AF3_COMPLEX_CACNA1C_CACNB2_CALM_CA` | CACNA1C, CACNB2, CALM1/2/3 | Review pending |
| `AF3_COMPLEX_SCN5A_CALM` | SCN5A, CALM1/2/3 | Review pending |
| `AF3_COMPLEX_RYR2_CALM_FKBP1B_CA` | RYR2, CALM1/2/3 | Review pending |
| `AF3_COMPLEX_CASQ2_TRDN_RYR2_CA` | CASQ2, TRDN, RYR2 | Review pending |
| `AF3_COMPLEX_TROPONIN_THIN_FILAMENT_CA` | ACTC1, TPM1, TNNI3, TNNT2 | Review pending |
| `AF3_COMPLEX_DESMOSOME_CORE` | DSP, PKP2, DSG2, DSC2, JUP | Review pending |
| `AF3_COMPLEX_LMNA_EMD` | LMNA, EMD | Review pending |
| `AF3_COMPLEX_PLN_ATP2A2` | PLN | Review pending |
| `AF3_COMPLEX_RAS_MAPK_GTP` | HRAS, KRAS, NRAS, RAF1, MAP2K2, PTPN11, RIT1, SOS1 | Review pending |

Rows in `af3_full_panel_partner_manifest.tsv` intentionally use:

```text
final_include = review_pending
manual_review_status = needs_literature_api_confirmation_before_submission
```

for whitelist-derived rows. Genes without a starter complex are kept as:

```text
partner_type = none_found
manual_review_status = manual_review_needed
final_include = no
exclusion_reason = no_curated_partner_assigned_yet
```

This is deliberate. It prevents weak or unverified interaction assumptions from
becoming AF3 server submissions.

Targets still needing partner review:

```text
ANK2, BAG3, CSRP3, DES, FHOD3, FLNC, GLA, HCN4, KCNE2, KCNH2,
KCNJ2, LAMP2, MIB1, MYBPC3, MYH7, MYL2, MYL3, PRKAG2, RBM20,
SCN1B, SLC22A5, SNTA1, TCAP, TMEM43, TRPM4
```

Next evidence-confirmation step:

1. Query/curate UniProt, Reactome, IntAct, and RCSB/PDB for each starter
   whitelist row.
2. Fill `pubmed_ids`, `reactome_pathway_or_complex`, `uniprot_interaction_note`,
   `intact_or_string_id`, and `pdb_or_rcsb_id`.
3. Add additional curated partners for the 25 genes still marked
   `manual_review_needed`.
4. Only after that, create AF3 complex groups and manual server submission
   packets.

## Update 2026-07-03: UniProt Evidence Pass Started

The first evidence-confirmation pass was completed using the official UniProt
REST API.

New scripts:

```text
datasets/feature_sources/alphafold3/scripts/enrich_af3_partner_manifest_uniprot.py
datasets/feature_sources/alphafold3/scripts/summarize_af3_evidence_review.py
```

New outputs:

```text
datasets/feature_sources/alphafold3/raw/uniprot_evidence/
datasets/feature_sources/alphafold3/interim/evidence/uniprot_evidence_by_gene.tsv
results/af3/af3_full_panel_partner_manifest.uniprot_evidence.tsv
results/af3/af3_uniprot_partner_evidence.summary.json
results/af3/af3_complex_evidence_readiness.tsv
results/af3/af3_partner_evidence_review_queue.tsv
results/af3/af3_evidence_review.summary.json
```

The script cached raw UniProt JSON for all sequence entities, then extracted:

- SUBUNIT text.
- FUNCTION/activity context.
- PubMed IDs attached to UniProt evidence statements.
- UniProt interaction partners from the INTERACTION comment.
- IntAct IDs exposed through UniProt interaction records.
- Reactome, ComplexPortal, PDB, CORUM, BioGRID, and DIP cross-references.

Current UniProt evidence summary:

| Item | Count |
|---|---:|
| Sequence entities checked | `59` |
| UniProt JSON files cached | `59` |
| Gene evidence rows | `59` |
| Partner rows enriched | `172` |

Pair-support counts:

| UniProt pair support level | Rows |
|---|---:|
| `direct_uniprot_intact_interaction` | `18` |
| `uniprot_subunit_text_mention` | `39` |
| `uniprot_calcium_context` | `24` |
| `uniprot_gtp_context` | `6` |
| `calmodulin_family_proxy_review_needed` | `13` |
| `calmodulin_family_context_review_needed` | `2` |
| `nonprotein_entity_review_needed` | `4` |
| `not_confirmed_by_uniprot_yet` | `41` |
| `no_partner_assigned` | `25` |

Complex-level readiness after UniProt-only evidence:

| Readiness status | Complexes |
|---|---:|
| `uniprot_supported_ready_for_manual_final_check` | `3` |
| `partially_supported_needs_external_evidence` | `7` |

The three complexes with enough UniProt support for manual final checking are:

```text
AF3_COMPLEX_LMNA_EMD
AF3_COMPLEX_PLN_ATP2A2
AF3_COMPLEX_SCN5A_CALM
```

The seven partially supported complexes are biologically plausible but still
need Reactome, IntAct, RCSB/PDB, or literature confirmation before submission:

```text
AF3_COMPLEX_CACNA1C_CACNB2_CALM_CA
AF3_COMPLEX_CASQ2_TRDN_RYR2_CA
AF3_COMPLEX_DESMOSOME_CORE
AF3_COMPLEX_KCNQ1_KCNE1_CALM_CA
AF3_COMPLEX_RAS_MAPK_GTP
AF3_COMPLEX_RYR2_CALM_FKBP1B_CA
AF3_COMPLEX_TROPONIN_THIN_FILAMENT_CA
```

Important QC notes:

- The original starter manifest was preserved.
- The UniProt-enriched table is a new copy, not an overwrite.
- CALM1, CALM2, and CALM3 canonical local UniProt sequences are identical
  length-149 proteins. Therefore a single calmodulin sequence can be used as an
  AF3 entity, but the manifest still preserves CALM1/2/3 as separate target
  genes for variant mapping and review.
- UniProt is strong for some pairwise interactions, such as KCNQ1-KCNE1,
  CACNA1C-CACNB2, RYR2-FKBP1B, LMNA-EMD, PLN-ATP2A2, and SCN5A-calmodulin
  context.
- UniProt alone is not sufficient for several expected structural complexes,
  especially thin filament and CASQ2/TRDN/RYR2. These require external evidence
  review before AF3 jobs.

The next evidence step should fill the `not_confirmed_by_uniprot_yet` and
`none_found` rows using Reactome, IntAct, RCSB/PDB, and targeted literature
review. The review queue is:

```text
results/af3/af3_partner_evidence_review_queue.tsv
```

## Update 2026-07-03: Reactome, RCSB/PDB, And IntAct Evidence Added

The next structured evidence pass was completed using:

- Reactome ContentService.
- RCSB Search API.
- EBI IntAct PSICQUIC service.

New script:

```text
datasets/feature_sources/alphafold3/scripts/enrich_af3_partner_manifest_external_sources.py
datasets/feature_sources/alphafold3/scripts/summarize_af3_combined_evidence.py
```

New raw caches:

```text
datasets/feature_sources/alphafold3/raw/reactome_evidence/
datasets/feature_sources/alphafold3/raw/rcsb_evidence/
datasets/feature_sources/alphafold3/raw/intact_evidence/
```

New outputs:

```text
datasets/feature_sources/alphafold3/interim/evidence/external_evidence_by_gene.tsv
results/af3/af3_full_panel_partner_manifest.external_evidence.tsv
results/af3/af3_none_found_candidate_partners.tsv
results/af3/af3_external_partner_evidence.summary.json
results/af3/af3_complex_combined_evidence_readiness.tsv
results/af3/af3_partner_combined_evidence_review_queue.tsv
results/af3/af3_combined_evidence_readiness.summary.json
```

External structured support summary:

| External support level | Rows |
|---|---:|
| `intact_psicquic_pair_support` | `38` |
| `rcsb_shared_pdb_support` | `19` |
| `reactome_shared_pathway_context` | `41` |
| `reactome_target_context_for_nonprotein_entity` | `34` |
| `not_supported_by_external_structured_sources_yet` | `15` |
| `no_partner_assigned` | `25` |

Combined UniProt + Reactome + RCSB + IntAct pair support:

| Combined support level | Rows |
|---|---:|
| `uniprot_subunit_text_mention` | `39` |
| `direct_uniprot_intact_interaction` | `18` |
| `intact_psicquic_pair_support` | `13` |
| `rcsb_shared_pdb_support` | `11` |
| `reactome_shared_pathway_context` | `21` |
| `uniprot_calcium_context` | `24` |
| `uniprot_gtp_context` | `6` |
| `reactome_target_context_for_nonprotein_entity` | `4` |
| `calmodulin_family_proxy_review_needed` | `7` |
| `needs_manual_literature_or_domain_review` | `4` |
| `no_partner_assigned` | `25` |

Current complex readiness:

| Status | Complexes |
|---|---:|
| `structured_evidence_ready_for_manual_af3_design` | `6` |
| `mostly_supported_needs_domain_and_literature_review` | `4` |

Complexes now ready for manual AF3 design review:

```text
AF3_COMPLEX_CASQ2_TRDN_RYR2_CA
AF3_COMPLEX_DESMOSOME_CORE
AF3_COMPLEX_LMNA_EMD
AF3_COMPLEX_PLN_ATP2A2
AF3_COMPLEX_SCN5A_CALM
AF3_COMPLEX_TROPONIN_THIN_FILAMENT_CA
```

Complexes still needing targeted domain/literature review before an AF3 packet:

```text
AF3_COMPLEX_CACNA1C_CACNB2_CALM_CA
AF3_COMPLEX_KCNQ1_KCNE1_CALM_CA
AF3_COMPLEX_RAS_MAPK_GTP
AF3_COMPLEX_RYR2_CALM_FKBP1B_CA
```

Important examples of gaps filled by external structured sources:

- Thin filament gained RCSB shared-PDB support for TNNI3/TNNT2/TPM1/TNNC1
  combinations and IntAct support for ACTC1-TPM1 and TNNT2-TNNC1.
- Desmosome gained strong combined support through UniProt, IntAct, and Reactome
  context.
- CASQ2/TRDN/RYR2 gained Reactome pathway-context support, but still needs
  careful domain design because a single full-length RYR2 complex is not
  feasible.

The 25 genes still marked `no_partner_assigned` now have a raw candidate table:

```text
results/af3/af3_none_found_candidate_partners.tsv
```

This table is intentionally **not** final. IntAct can include aliases,
non-cardiac interactors, and high-throughput interactions. Candidate partners
must be manually reviewed before being promoted into the AF3 partner manifest.

Still no AF3 job packets have been generated. This remains an evidence and
review layer only.

## Update 2026-07-03: AF3 Design Review And First Draft Packets

The next safe step was completed:

```text
AF3 design review before server submission
```

New script:

```text
datasets/feature_sources/alphafold3/scripts/build_af3_complex_design_review.py
```

New outputs:

```text
results/af3/af3_complex_design_review.tsv
results/af3/af3_complex_groups.tsv
results/af3/af3_complex_design_review.summary.json
```

The design review checks:

- Combined evidence readiness.
- Protein entity lengths from reviewed UniProt sequences.
- Estimated AF3 token burden.
- Whether full-length input is feasible.
- Whether domain/split design is required.
- Protein-position coverage from the current modeling matrix.

Current design-review result:

| Packet/design status | Complexes |
|---|---:|
| `ready_for_manual_packet_draft` | `4` |
| `hold_for_manual_review` | `2` |
| `not_ready_for_packet` | `4` |

Ready for draft manual AF3 packets:

| Complex | Entities | Estimated tokens |
|---|---|---:|
| `AF3_COMPLEX_SCN5A_CALM` | SCN5A, CALM1, Ca2+ | `2166` |
| `AF3_COMPLEX_TROPONIN_THIN_FILAMENT_CA` | ACTC1, TPM1, TNNI3, TNNT2, TNNC1, Ca2+ | `1331` |
| `AF3_COMPLEX_LMNA_EMD` | LMNA, EMD | `918` |
| `AF3_COMPLEX_PLN_ATP2A2` | PLN, ATP2A2, Ca2+ | `1095` |

Hold for manual review before packet:

| Complex | Why |
|---|---|
| `AF3_COMPLEX_KCNQ1_KCNE1_CALM_CA` | Token-feasible, but membrane-span context and CALM isoform proxy should be explicitly reviewed. |
| `AF3_COMPLEX_CACNA1C_CACNB2_CALM_CA` | Token-feasible, but full-length CACNA1C is a large membrane protein; beta-subunit/CALM domain context should be reviewed. |

Not ready for packet:

| Complex | Why |
|---|---|
| `AF3_COMPLEX_RYR2_CALM_FKBP1B_CA` | Full RYR2 plus partners exceeds practical token limit; domain selection required. |
| `AF3_COMPLEX_CASQ2_TRDN_RYR2_CA` | Full-length RYR2 makes the complex too large; domain/split design required. |
| `AF3_COMPLEX_DESMOSOME_CORE` | Full desmosome core exceeds token limit; should be split into smaller extracellular/plaque contexts. |
| `AF3_COMPLEX_RAS_MAPK_GTP` | Pathway context, not one guaranteed physical complex; should be split into specific RAS-RAF/other domain jobs. |

Draft manual AF3 submission packets were then generated only for the four
`ready_for_manual_packet_draft` complexes.

New script:

```text
datasets/feature_sources/alphafold3/scripts/build_af3_manual_submission_packets.py
```

New outputs:

```text
datasets/feature_sources/alphafold3/jobs/af3_job_manifest.tsv
datasets/feature_sources/alphafold3/jobs/manual_submission_packets/
results/af3/af3_manual_submission_packets.summary.json
```

Draft packet folders:

```text
AF3_JOB_001_SCN5A_CALM
AF3_JOB_002_TROPONIN_THIN_FILAMENT_CA
AF3_JOB_003_LMNA_EMD
AF3_JOB_004_PLN_ATP2A2
```

Each packet contains:

```text
job_name.txt
entities.fasta
ligands_and_ions.txt
expected_chain_map.tsv
evidence_summary.txt
notes_for_server_ui.txt
```

Important status:

```text
manual_ready_for_submission = draft_needs_user_review
af3_status = not_submitted
```

Nothing has been submitted to AlphaFold Server. These packets are ready for
human review and manual copy/paste into the AF3 server only after approval.

## Update 2026-07-03: Gene-Level AF3 Coverage Review

The project was zoomed back to one row per gene so the AF3 setup can be checked
before moving forward.

New script:

```text
datasets/feature_sources/alphafold3/scripts/build_af3_gene_level_coverage_review.py
```

New outputs:

```text
results/af3/af3_gene_level_coverage_review.tsv
results/af3/af3_gene_level_coverage_review.summary.json
```

Current gene-level status:

| Category | Gene count | Meaning |
|---|---:|---|
| AF3 target genes reviewed | `56` | Included/sensitivity/monitor genes considered for AF3 |
| Reviewed UniProt protein sequence available | `56` | All AF3 target genes have canonical reviewed human protein sequence |
| Packet draft ready | `11` | Protein sequence + partner/complex + evidence + feasible draft packet |
| Partner assigned but held for manual review | `4` | Complex is token-feasible but needs domain/literature review |
| Partner assigned but domain/split required | `16` | Biology/evidence exists, but full AF3 job is too large or not a single clean complex |
| No partner assigned yet | `25` | Sequence exists, but no final defensible AF3 partner/complex has been selected |

Packet-draft-ready genes:

```text
ACTC1, CALM1, CALM2, CALM3, EMD, LMNA, PLN, SCN5A, TNNI3, TNNT2, TPM1
```

Assigned but not packet-ready genes:

```text
CACNA1C, CACNB2, CASQ2, DSC2, DSG2, DSP, HRAS, JUP, KCNE1, KCNQ1,
KRAS, MAP2K2, NRAS, PKP2, PTPN11, RAF1, RIT1, RYR2, SOS1, TRDN
```

Genes still without final AF3 partner assignment:

```text
ANK2, BAG3, CSRP3, DES, FHOD3, FLNC, GLA, HCN4, KCNE2, KCNH2,
KCNJ2, LAMP2, MIB1, MYBPC3, MYH7, MYL2, MYL3, PRKAG2, RBM20,
SCN1B, SLC22A5, SNTA1, TCAP, TMEM43, TRPM4
```

The gene-level review table includes these checkable evidence columns:

```text
pubmed_ids_from_uniprot
pubmed_ids_from_intact
intact_ids_from_uniprot
reactome_ids
pdb_ids
evidence_source_urls
fix_needed
```

Important interpretation:

- Protein sequences are complete for the AF3 target set.
- Partner/complex assignment is **not complete** yet.
- The next work should focus on the 25 `no_partner_assigned` genes and the 20
  assigned-but-not-packet-ready genes before claiming the AF3 setup is finished.
- Candidate partners for the 25 unresolved genes exist in
  `results/af3/af3_none_found_candidate_partners.tsv`, but these are raw
  IntAct candidates and must be manually checked before promotion.

## Update 2026-07-03: Second-Wave Partner Decisions For Previously Unassigned Genes

The 25 genes that previously had no final partner assignment were reviewed using
UniProt SUBUNIT text, Reactome pathway context, RCSB/PDB structure support, and
IntAct candidate evidence. A second-wave decision table was generated.

New script:

```text
datasets/feature_sources/alphafold3/scripts/build_af3_second_wave_partner_decisions.py
datasets/feature_sources/alphafold3/scripts/build_af3_gene_level_coverage_review_v2.py
```

New outputs:

```text
results/af3/af3_second_wave_curated_partner_decisions.tsv
datasets/feature_sources/alphafold3/interim/af3_second_wave_entities_to_sequence.tsv
datasets/feature_sources/alphafold3/interim/af3_second_wave_entities.fasta
results/af3/af3_second_wave_partner_decisions.summary.json
results/af3/af3_gene_level_coverage_review_v2.tsv
results/af3/af3_gene_level_coverage_review_v2.summary.json
```

New external partner sequences added:

```text
ACTN2, CRYAB, HSPA8, HSPB8, MYOT, PDZK1, PRKAA2, PRKAB2,
SPTBN1, SQSTM1, STUB1
```

All 25 previously unassigned genes now have an explicit path:

| Second-wave decision type | Gene count |
|---|---:|
| Full-length packet candidate after final evidence check | `12` |
| Ligand/homomer context needing CCD or stoichiometry review | `4` |
| Domain-range review required | `6` |
| Cardiac-specificity review required | `1` |
| Partner selection still needed | `1` |
| RNA/spliceosome design deferred | `1` |

Second-wave full-length packet candidates:

```text
BAG3, CSRP3, DES, KCNE2, KCNH2, KCNJ2, LAMP2, SCN1B, SNTA1,
TCAP, TMEM43, TRPM4
```

Second-wave ligand/homomer context genes:

```text
GLA, HCN4, PRKAG2, SLC22A5
```

Second-wave domain-review genes:

```text
ANK2, FLNC, MYBPC3, MYH7, MYL2, MYL3
```

Special unresolved/deferred genes:

```text
FHOD3 = needs cardiac-specificity review for FHOD3-SQSTM1
MIB1 = needs specific partner selection in Notch context
RBM20 = defer until RNA/protein spliceosome AF3 design is explicitly scoped
```

Updated gene-level AF3 status after second-wave decisions:

| Status | Gene count |
|---|---:|
| `packet_draft_ready` | `11` |
| `second_wave_packet_candidate_needs_final_evidence_check` | `12` |
| `second_wave_ligand_or_homomer_context_needs_ccd_review` | `4` |
| `assigned_hold_for_manual_review` | `4` |
| `assigned_domain_or_split_required` | `16` |
| `second_wave_domain_review_required` | `6` |
| `second_wave_needs_cardiac_specificity_review` | `1` |
| `second_wave_needs_partner_selection` | `1` |
| `defer_until_rna_complex_design` | `1` |
| `no_partner_assigned` | `0` |

Important interpretation:

- Every AF3 target gene now has either a protein partner/complex path,
  homomer/ligand context, domain/split path, or explicit defer/review reason.
- This is **not** the same as every gene being AF3-server-ready.
- The next step should be final evidence confirmation for the 12 second-wave
  packet candidates, then packet generation for those that pass.
- IntAct-derived fields can contain broad/high-throughput evidence. For final
  submission decisions, prioritize UniProt SUBUNIT PubMed IDs, Reactome IDs,
  PDB support, and manually verified IntAct rows.

## What AlphaFold 3 Gives Versus What We Calculate

AlphaFold 3 gives structure and confidence outputs. It does not directly give a
CatBoost-ready feature table.

| Feature/output | Direct from AF3 | Local engine calculates |
|---|---|---|
| Predicted 3D coordinates | Yes | No |
| mmCIF structure files | Yes | No |
| pLDDT confidence | Yes | Aggregate per residue/window |
| PAE matrix / confidence JSON | Yes | Extract residue-partner values |
| pTM / ipTM / chain-pair confidence | Yes | Store as QC and feature |
| Chain IDs | Yes | Map chain ID to gene/entity |
| Interface residue yes/no | No | Yes |
| Distance to partner chain | No | Yes |
| Distance to ligand/ion | No | Yes |
| Contact counts | No | Yes |
| Variant-level AF3 feature row | No | Yes |

Local feature extraction is therefore required.

## Full Workflow

### Phase 0 — Set Up AF3 Project Structure

Create the AF3 directories without touching the previous model outputs:

```text
datasets/feature_sources/alphafold3/
datasets/feature_sources/alphafold3/raw/
datasets/feature_sources/alphafold3/interim/
datasets/feature_sources/alphafold3/jobs/
datasets/feature_sources/alphafold3/jobs/manual_submission_packets/
datasets/feature_sources/alphafold3/processed/
datasets/feature_sources/alphafold3/scripts/
results/af3/
```

### Phase 1 — Recover and Audit the Full Gene Panel

Use the current ready matrix first:

```text
datasets/modeling/ready/modeling_table_with_splits_weights_hiro_emerge_rescued.tsv
```

Also keep the intermediate matrix as a reference:

```text
datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx_vep_hgvs_spliceai_hiro_agg.tsv
```

Generate:

```text
results/af3/gene_panel_audit.tsv
```

Required columns:

```text
gene_raw
gene_normalized
n_total_variants
n_missense
n_plp
n_blb
n_vus
n_hgvsp_nonmissing
n_protein_position_nonmissing
n_uniprot_mapped
modeling_scope_values
primary_model_inclusion_values
panel_bucket
ttn_separate_flag
excluded_or_quarantine_flag
case_alias_warning
af3_coverage_status
af3_complex_group
af3_manual_review_needed
```

Gene-symbol checks:

```text
KCNH2 present?
KNCH2 present?
RYR2/RyR2 duplicate?
KCNE1 present in both included and excluded scopes?
TTN separate?
```

### Phase 2 — Build the Gene-to-Partner Manifest

Goal:

```text
For every normalized gene, identify the best defensible biological context for AF3.
```

Create:

```text
results/af3/af3_full_panel_partner_manifest.tsv
```

Required columns:

```text
target_gene
target_uniprot
target_protein_name
target_length
target_panel_status
complex_group_id
complex_name
partner_gene_or_entity
partner_type
partner_uniprot_or_identifier
partner_length_or_token_count
interaction_source
interaction_evidence_type
pubmed_ids
reactome_pathway_or_complex
uniprot_interaction_note
intact_or_string_id
pdb_or_rcsb_id
ligand_or_ion_id
recommended_af3_input_mode
domain_start
domain_end
reason_for_domain_trimming
token_estimate
manual_review_status
final_include
exclusion_reason
```

Partner/entity types:

```text
protein
protein_domain
peptide
ion
ligand
DNA
RNA
PTM
none_found
```

Partner sources should be searched in this order:

1. Manual cardiac complex whitelist.
2. UniProt reviewed human entry and curated interaction notes.
3. Reactome complexes/pathways.
4. IntAct or curated physical-interaction evidence.
5. PDB/RCSB co-complex evidence.
6. STRING only as weak backup requiring manual review.

Do not let STRING alone justify a final AF3 job. STRING includes functional
associations, not only direct binding.

### Phase 3 — Retrieve Canonical Protein Sequences

Codex only has gene symbols, so each target and partner needs a sequence.

Preferred sequence source:

```text
Reviewed human Swiss-Prot UniProt entry
```

Save:

```text
datasets/feature_sources/alphafold3/interim/uniprot_gene_to_sequence.tsv
datasets/feature_sources/alphafold3/interim/af3_all_entities.fasta
```

Required columns:

```text
gene
uniprot_accession
entry_name
protein_name
organism
reviewed
canonical_sequence
sequence_length
isoform_note
retrieval_date
retrieval_source
manual_review_needed
```

Rules:

- Use reviewed Swiss-Prot human entries first.
- Prefer canonical isoform.
- Do not use unreviewed TrEMBL unless no reviewed human entry exists.
- If multiple reviewed entries exist, require manual review.
- If variant reference amino acid does not match the UniProt sequence later,
  mark the AF3 mapping as `ref_aa_mismatch`.

### Phase 4 — Group Genes Into Non-Redundant AF3 Complex Jobs

Create:

```text
results/af3/af3_complex_groups.tsv
```

Required columns:

```text
complex_group_id
complex_name
genes_covered
protein_entities
ligand_entities
ion_entities
dna_entities
rna_entities
total_tokens_estimated
full_length_feasible
domain_trim_required
manual_submission_name
expected_variant_coverage
notes
```

The key efficiency gain is shared complexes.

Examples:

| Complex group | Genes covered | Entities |
|---|---|---|
| `AF3_COMPLEX_KCNQ1_KCNE1_CALM_CA` | KCNQ1, KCNE1, CALM1, CALM2, CALM3 | KCNQ1, KCNE1, CALM1, Ca2+ |
| `AF3_COMPLEX_TROPONIN_THIN_FILAMENT_CA` | TNNI3, TNNT2, TNNC1 if added, TPM1, ACTC1 | TNNI3, TNNT2, TNNC1, TPM1, ACTC1, Ca2+ |
| `AF3_COMPLEX_LMNA_EMD` | LMNA, EMD | LMNA domain, EMD domain |
| `AF3_COMPLEX_RYR2_FKBP1B_OR_CALM` | RYR2, FKBP1B/CALM partner context | RYR2 domain, FKBP1B or CALM |
| `AF3_COMPLEX_SCN5A_CALM` | SCN5A, CALM1/2/3 | SCN5A IQ/C-terminal region, CALM |
| `AF3_COMPLEX_MYO_FILAMENT` | MYH7, MYBPC3, ACTC1 if feasible | MYH7/MYBPC3/actin fragments |

AF3 jobs must obey:

- No job exceeds the AF3 server token limit.
- Full length is preferred only when feasible and biologically appropriate.
- Domain-level jobs are allowed for very large proteins or specific domains.
- Do not submit unrelated proteins together to save quota.
- Every job must represent one plausible biological complex.

### Phase 5 — Build AF3 Server Submission Packets

Create:

```text
datasets/feature_sources/alphafold3/jobs/af3_job_manifest.tsv
datasets/feature_sources/alphafold3/jobs/manual_submission_packets/
```

Required job manifest columns:

```text
job_id
complex_group_id
job_name
genes_covered
entity_1_type
entity_1_name
entity_1_sequence_or_id
entity_1_range
entity_2_type
entity_2_name
entity_2_sequence_or_id
entity_2_range
entity_3_type
entity_3_name
entity_3_sequence_or_id
entity_3_range
ligands
ions
total_tokens_estimated
full_length_or_domain
reason_for_job
evidence_summary
manual_ready_for_submission
submission_day
submission_order
af3_status
af3_output_zip_path
```

Each manual submission packet should contain:

```text
job_name.txt
entities.fasta
ligands_and_ions.txt
notes_for_server_ui.txt
expected_chain_map.tsv
evidence_summary.txt
```

Example packet for KCNQ1:

```text
job_name: KCNQ1_KCNE1_CALM1_CA_WT

entities:
>KCNQ1_HUMAN
[canonical KCNQ1 sequence]

>KCNE1_HUMAN
[canonical KCNE1 sequence]

>CALM1_HUMAN
[canonical CALM1 sequence]

ion:
Ca2+

expected_chain_map:
A = KCNQ1
B = KCNE1
C = CALM1
D = Ca2+
```

Important: expected chain map is a prediction only. The actual chain map must
be verified from the AF3 output mmCIF and job JSON.

### Phase 6 — AF3 Server Production

Expected practical production:

| Step | Approximate scale |
|---|---:|
| Included normalized genes to review | `53` |
| AF3 jobs after shared-complex grouping | likely `30-55` |
| AF3 server quota assumption | `30 jobs/day` |
| Production time after manifests are ready | about `2 days` |

For each completed job, download the full output zip immediately.

Store:

```text
datasets/feature_sources/alphafold3/raw/{job_id}/
```

Expected files per job:

```text
job_request.json
fold_<job_name>_model_0.cif
fold_<job_name>_model_1.cif
fold_<job_name>_model_2.cif
fold_<job_name>_model_3.cif
fold_<job_name>_model_4.cif
fold_<job_name>_summary_confidences_0.json
fold_<job_name>_full_data_0.json
README.md
```

Also create:

```text
datasets/feature_sources/alphafold3/interim/af3_completed_jobs.tsv
```

### Phase 7 — Local AF3 Feature Extraction Engine

AlphaFold 3 gives the structure and confidence files. Our code creates the
features.

Create script:

```text
datasets/feature_sources/alphafold3/scripts/build_af3_residue_features.py
```

Parse:

```text
mmCIF structures
summary confidence JSON
full-data JSON / PAE
job request JSON
expected and verified chain maps
```

Recommended Python packages:

```text
gemmi
biopython
pandas
numpy
scipy
```

Residue-level output:

```text
datasets/feature_sources/alphafold3/processed/af3_residue_features.tsv
```

Required residue-level columns:

```text
af3_complex_group_id
af3_job_id
gene
uniprot_accession
chain_id
residue_number
residue_aa
local_plddt
min_distance_to_any_partner_A
min_distance_to_partner_gene_A
interface_any_partner_4A
interface_any_partner_5A
interface_any_partner_8A
contact_count_any_partner_5A
min_distance_to_ion_A
near_ion_5A
min_distance_to_ligand_A
near_ligand_5A
chain_pair_iptm
chain_pair_pae_min
local_interface_pae_mean
af3_model_rank
af3_complex_confidence_flag
```

Definitions:

```text
interface_any_partner_5A =
  any heavy atom in target residue within 5 Angstrom of any heavy atom
  in a partner chain

near_ion_5A =
  any heavy atom in target residue within 5 Angstrom of an ion atom

local_plddt =
  mean pLDDT over all atoms in the residue

local_interface_pae_mean =
  mean PAE between target residue token and contacted partner-chain residues
```

### Phase 8 — Map AF3 Features to Variants

Create script:

```text
datasets/feature_sources/alphafold3/scripts/build_af3_variant_features.py
```

Merge keys:

```text
gene
uniprot_accession
protein_position
reference_amino_acid
```

Reference amino-acid QC is mandatory:

```text
variant_ref_aa == uniprot_sequence[protein_position]
```

If this fails:

```text
af3_mapping_status = ref_aa_mismatch
```

Do not impute AF3 features.

Variant-level output:

```text
datasets/feature_sources/alphafold3/processed/af3_variant_features.tsv
```

Required columns:

```text
variant_id
gene
hgvsp
protein_position
ref_aa
alt_aa
af3_mapping_status
af3_complex_available
af3_complex_group_id
af3_best_job_id
af3_residue_in_any_interface_5A
af3_min_distance_to_any_partner_A
af3_min_distance_to_best_partner_A
af3_best_partner_gene
af3_near_ion_5A
af3_near_ligand_5A
af3_local_plddt
af3_local_interface_pae_mean
af3_chain_pair_iptm
af3_feature_missing_reason
```

Missing reasons:

```text
not_missense
no_protein_position
no_uniprot_mapping
ref_aa_mismatch
gene_not_in_af3_panel
no_curated_partner
domain_not_covering_residue
af3_job_failed
low_confidence_complex
excluded_or_quarantine_gene
ttn_not_modeled_full_length
```

### Phase 9 — Merge AF3 Features Into a New Model Version

Do not overwrite the current primary model.

Create a new matrix:

```text
datasets/modeling/ready/modeling_table_with_af3_features.tsv
```

Create a new model family:

```text
results/models/primary_binary_catboost_af3_v1/
results/model_performance/primary_binary_catboost_af3_v1/
```

Compare:

```text
Baseline primary binary CatBoost
vs
Primary binary CatBoost + AF3 features
```

Run comparisons in three layers:

| Analysis | Purpose |
|---|---|
| Full matrix | Does AF3 help overall, despite missingness? |
| AF3-covered subset | Does AF3 help where it actually applies? |
| Mechanistic case studies | Does AF3 explain pathogenic-like VUS or variants missed by ΔΔG? |

Do not claim global superiority unless the full held-out/source-held-out tests
support it.

### Phase 10 — Manuscript Framing

Best framing if AF3 works:

> We added an interaction-aware structural feature layer using AlphaFold 3
> wild-type complex predictions. Instead of predicting every variant structure,
> we modeled curated cardiac molecular complexes and mapped variants to
> interfaces, ion/ligand proximity, and confidence-weighted structural context.

Reviewer-safe language:

- AF3 features are exploratory interaction-aware annotations.
- AF3 features are not direct functional assays.
- A variant near an interface is not automatically pathogenic.
- Confidence metrics such as pLDDT, PAE, and ipTM must be reported.
- Low-confidence complexes should be flagged and excluded or sensitivity-tested.

## Initial Candidate Complex Groups

These are starting hypotheses, not final submissions. Every row needs evidence
review before AF3 submission.

| Candidate complex | Genes covered | Notes |
|---|---|---|
| KCNQ1-KCNE1-CALM-Ca2+ | KCNQ1, KCNE1, CALM1/2/3 | Strong LQTS/IKs/calmodulin context |
| KCNH2 auxiliary/context complex | KCNH2, KCNE2 or other reviewed partner | Needs careful partner evidence |
| SCN5A-CALM | SCN5A, CALM1/2/3 | C-terminal/IQ region likely domain-level |
| SCN5A-ANK3/PKP2/SNTA1 context | SCN5A, ANK3/PKP2/SNTA1 | Partner evidence and domain trimming needed |
| Troponin/thin-filament-Ca2+ | TNNI3, TNNT2, TNNC1, TPM1, ACTC1 | One shared job may cover several sarcomere genes |
| MYH7-MYBPC3/ACTC1 context | MYH7, MYBPC3, ACTC1 | Likely domain/fragment-level |
| RYR2-FKBP1B/CALM | RYR2, FKBP1B, CALM | RYR2 must be domain-level |
| CASQ2-TRDN | CASQ2, TRDN | SR calcium-release context |
| LMNA-EMD/SUN/TMEM43 | LMNA, EMD, TMEM43 | Nuclear envelope context, domain-level likely |
| PLN-ATP2A2/SERCA2a | PLN, ATP2A2 | ATP2A2 is a partner, not necessarily in current variant panel |
| Desmosome complex | DSP, PKP2, DSG2, DSC2, JUP | Complex design may need multiple smaller jobs |
| RASopathy/phenocopy complexes | HRAS, KRAS, NRAS, RAF1, RIT1, MAP2K2, PTPN11, SOS1 | May be pathway/ligand/GTP context rather than direct cardiac complex |

## Immediate Next Steps

1. Create AF3 folder structure.
2. Generate `results/af3/gene_panel_audit.tsv`.
3. Normalize gene symbols and resolve `RYR2/RyR2`, `KCNH2/KNCH2`, and KCNE1 scope inconsistency.
4. Build the first manual cardiac complex whitelist.
5. Retrieve reviewed UniProt sequences for all included genes.
6. Build `af3_full_panel_partner_manifest.tsv`.
7. Build `af3_complex_groups.tsv`.
8. Build AF3 manual submission packets.
9. After user submits jobs and downloads outputs, build local mmCIF/JSON parser.
10. Merge AF3 variant features into a new AF3 model matrix and train `primary_binary_catboost_af3_v1`.

## Current Status

| Task | Status |
|---|---|
| AF3 plan created | Done in this file |
| Existing primary model preserved | Yes |
| Current matrix gene count checked | Done |
| Included normalized gene count | `53` |
| AF3 folder structure | Done |
| Gene panel audit file | Done: `results/af3/gene_panel_audit.tsv` |
| Partner manifest | Not done yet |
| UniProt sequence retrieval | Not done yet |
| AF3 job manifest | Not done yet |
| AF3 server submissions | Not started |
| AF3 local feature engine | Not started |
| AF3-augmented model | Not started |

## Generated Files, 2026-07-03

AF3 folder structure:

```text
datasets/feature_sources/alphafold3/
datasets/feature_sources/alphafold3/raw/
datasets/feature_sources/alphafold3/interim/
datasets/feature_sources/alphafold3/jobs/
datasets/feature_sources/alphafold3/jobs/manual_submission_packets/
datasets/feature_sources/alphafold3/processed/
datasets/feature_sources/alphafold3/scripts/
results/af3/
```

Gene-panel audit outputs:

```text
results/af3/gene_panel_audit.tsv
results/af3/gene_panel_audit.summary.json
datasets/feature_sources/alphafold3/scripts/build_gene_panel_audit.py
datasets/feature_sources/alphafold3/README.md
```

Audit summary:

| Item | Count |
|---|---:|
| Normalized genes in ready matrix | `60` |
| Raw gene strings represented | `61` |
| Included normalized genes | `53` |
| Primary genes | `37` |
| Expanded-scope genes | `11` |
| Primary-monitor genes | `5` |
| Sensitivity-only genes | `3` |
| Exclude/quarantine genes | `3` |
| TTN-separate genes | `1` |

Important audit warnings:

| Warning | Status |
|---|---|
| `KCNH2` present | Yes |
| `KNCH2` present | No |
| `RYR2` / `RyR2` mixed case | Yes, collapsed to `RYR2` |
| `KCNE1` mixed include/quarantine scope | Yes, but user decision is to include KCNE1 for AF3 |

## Key Principle

Full-panel AF3 coverage means:

```text
Every gene gets a curated coverage decision.
```

It does **not** mean:

```text
Every gene is forced into a weak or unrelated AF3 complex.
```

That distinction is what will keep this feature layer scientifically defensible.



big issue to keep in mind : 

Yes. This is the exact right question.

The answer is:

> We get protein partners from curated biology databases and literature. Then we decide whether to submit **one multi-protein complex** or **separate pairwise complexes** based on whether those partners truly bind together at the same time.

AlphaFold 3 does **not** magically know which partners should be included for your gene. It predicts a structure for the molecules **you provide**. So the input design matters a lot.

## 1) How do we get the protein partners?

For each gene in your 52-ish gene panel, Codex should build a **partner evidence table** from multiple sources.

The hierarchy should be:

### Tier 1: Curated known complexes

Best source when available.

Examples:

```text
TNNI3 + TNNT2 + TNNC1 + TPM1 + ACTC1 + Ca2+
KCNQ1 + KCNE1 + calmodulin + Ca2+
MYBPC3 + MYH7 / actin fragment
```

These are ideal because they are real biological complexes, not random associations.

### Tier 2: UniProt reviewed entry

For each gene, get the reviewed human UniProt entry and read:

```text
Subunit structure
Interaction section
Cofactor / binding-site annotations
Domain annotations
Protein sequence
Isoform notes
```

UniProt interaction sections include binary protein interactions and can link to IntAct evidence, while the subunit section describes known quaternary structure and interaction context. ([UniProt][1])

### Tier 3: Reactome complexes/pathways

Reactome is useful for identifying whether several proteins participate in the **same biological complex or pathway**, not just whether they are loosely associated. Reactome is curated and peer reviewed, and it has a content service API that Codex can query programmatically. ([Reactome][2])

### Tier 4: IntAct / curated experimental interactions

IntAct is good for experimentally supported protein-protein interactions. UniProt imports a filtered subset of IntAct interactions and links back to the IntAct records. ([UniProt][3])

### Tier 5: PDB/RCSB co-complex evidence

If there is an experimentally solved structure of the protein with a partner, that is very strong evidence. It also tells you which **domain/residue range** to model if the full protein is too large.

### Tier 6: STRING as backup only

STRING is useful for broad discovery, but it combines **physical interactions and functional associations**. That means STRING can tell you proteins are functionally related without proving they directly bind. So STRING should be used as a backup, not as enough evidence by itself for a final AF3 complex. ([PMC][4])

## 2) Do we submit protein + all partners, or protein-partner1 separately?

Use this rule:

```text
If the partners are known to exist in the same biological complex at the same time:
    submit them together.

If the partners are alternative, context-dependent, or not known to bind simultaneously:
    submit separate pairwise/domain-level jobs.
```

### Case A: True multi-protein complex

Submit together.

Example:

```text
TNNI3 + TNNT2 + TNNC1 + TPM1 + ACTC1 + Ca2+
```

This makes sense because these proteins are part of the cardiac thin-filament/troponin regulatory system. A single AF3 job can annotate variants across multiple genes.

### Case B: One protein has multiple alternative partners

Submit separately.

Example:

```text
SCN5A + calmodulin
SCN5A + ANK3 domain
SCN5A + PKP2 domain
```

Do **not** necessarily submit:

```text
SCN5A + calmodulin + ANK3 + PKP2
```

unless there is evidence they bind the same SCN5A molecule at the same time in one complex.

### Case C: Regulatory ion/ligand needed for the biology

Include the ion/ligand.

Example:

```text
TNNI3 + TNNT2 + TNNC1 + Ca2+
CALM + Ca2+
KCNQ1 + KCNE1 + CALM + Ca2+
```

AlphaFold Server supports complexes with proteins, DNA, RNA, ligands, ions, and modifications. ([AlphaFold Server][5])

### Case D: Huge protein

Use a domain or fragment.

Example:

```text
RYR2 full-length + FKBP1B
```

may exceed the token limit or be structurally too huge. Better:

```text
RYR2 interaction domain + FKBP1B
```

AlphaFold Server jobs have token limits, so full-length mega-protein complexes are not always practical. ([AlphaFold Server][5])

## 3) Does AlphaFold 3 know whether they bind at the same time?

Not reliably in the way you need for a manuscript.

AF3 can predict the structure of the molecules you provide. But it does **not** replace biological curation of the complex.

If you give AF3:

```text
Protein A + Protein B + Protein C
```

it will attempt to model that system. If B and C do not actually bind A at the same time, AF3 might:

```text
predict weak/no interface,
place one partner poorly,
produce low inter-chain confidence,
or generate a biologically misleading arrangement.
```

That is why you need to classify partners before submission.

## 4) Does it matter if they bind at the same time?

Yes. It matters a lot.

If they truly bind simultaneously, a multi-protein AF3 job is valuable because it can show:

```text
shared interface
competition/cooperative positioning
ion/ligand placement
multi-chain confidence
```

If they do **not** bind simultaneously, forcing them into one job can create misleading features.

For CatBoost, a wrong complex could create false features like:

```text
variant is near Partner C
```

when Partner C should never have been in that complex.

So the safest rule is:

```text
Only put all partners together when there is curated evidence they are part of the same complex.

Otherwise run separate jobs:
target + partner1
target + partner2
target + partner3
```

## 5) What should Codex actually do?

For every gene, Codex should create a table like this:

```text
target_gene
partner_gene_or_entity
evidence_source
evidence_type
same_complex_evidence
simultaneous_binding_evidence
direct_physical_interaction
protein_domain_needed
include_in_multichain_job
include_as_pairwise_job
manual_review_needed
```

Then apply this decision logic:

```text
If Reactome/PDB/literature says same complex:
    group into one AF3 multichain job.

If UniProt/IntAct says direct physical interaction but no same-complex evidence:
    make pairwise job.

If only STRING functional association:
    do not make final production job unless manually approved.

If full-length job exceeds 5,000 tokens:
    use domain-level job.

If partner is an ion/ligand/cofactor with binding-site evidence:
    include it in the relevant job.
```

## 6) Example: KCNQ1

Possible partners/entities:

```text
KCNE1
CALM1/CALM2/CALM3
Ca2+
AKAP9
KCNE2
```

How to decide:

### Good multi-entity job

```text
KCNQ1 + KCNE1 + CALM + Ca2+
```

because this represents a coherent channel/regulatory complex.

### Separate jobs if needed

```text
KCNQ1 + KCNE1
KCNQ1 + CALM + Ca2+
KCNQ1 + AKAP9 domain
```

Do not automatically submit:

```text
KCNQ1 + KCNE1 + KCNE2 + KCNE3 + AKAP9 + CALM + Ca2+
```

That may be biologically messy unless there is evidence all entities bind simultaneously.

## 7) What AF3 gives you after submission

AF3 gives you:

```text
complex coordinates
chain positions
ligand/ion positions
pLDDT
PAE
pTM
ipTM / chain-pair confidence summaries
```

The AlphaFold Server output includes `.cif` structures plus confidence JSON files and full-data JSON files with detailed PAE information. ([EMBL-EBI][6])

Then your local engine calculates:

```text
distance from variant residue to partner
whether residue is at interface
distance to Ca2+/ligand
local pLDDT
local interface PAE
chain-pair ipTM
contact counts
```

So AF3 gives the raw structure/confidence. Your code turns it into CatBoost features.

## Bottom line

For your full-panel plan:

```text
Do not put every possible partner into one massive AF3 job.
Do not rely on AF3 to decide which proteins biologically belong together.
Do build a partner-evidence table first.
Do submit multi-protein jobs only for true known complexes.
Do submit separate pairwise/domain jobs for alternative or independent partners.
Do include ions/ligands only when biologically justified.
```

The clean manuscript framing is:

> Partner selection was evidence-guided. Proteins were modeled together only when curated sources supported participation in the same biological complex; otherwise, pairwise or domain-level complexes were modeled separately. AlphaFold 3 outputs were then converted locally into residue- and variant-level interface features.

[1]: https://www.uniprot.org/help/binary_interactions?utm_source=chatgpt.com "Binary interactions | UniProt help"
[2]: https://reactome.org/?utm_source=chatgpt.com "Reactome Pathway Database: Home"
[3]: https://www.uniprot.org/help/binary_interactions_import?utm_source=chatgpt.com "Why are only a subset of binary interactions from the IntAct ..."
[4]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11701646/?utm_source=chatgpt.com "The STRING database in 2025: protein networks with ... - PMC"
[5]: https://alphafoldserver.com/?utm_source=chatgpt.com "AlphaFold Server"
[6]: https://www.ebi.ac.uk/training/online/courses/alphafold/alphafold-3-and-alphafold-server/alphafold-server-your-gateway-to-alphafold-3/a-step-by-step-guide-to-generating-predictions-with-alphafold-server/?utm_source=chatgpt.com "AlphaFold"


*********

Exactly. That is the correct understanding.

This is **more complicated than “run AF3 for 52 proteins”**, but it is still manageable if we structure it properly.

The real workflow is:

```text
gene → canonical protein → known partners/entities → biological complex design → AF3 job → local feature extraction → variant-level CatBoost features
```

So yes, for each gene we need to understand:

```text
1. What protein does this gene encode?
2. What are the known binding partners?
3. Are those partners direct physical binders or just pathway-associated?
4. Do they bind at the same time or separately?
5. Should AF3 model them together or as separate pairwise/domain complexes?
6. Which variants in your matrix map to residues in those complexes?
7. How do we summarize all of that under the same gene/variant umbrella for CatBoost?
```

## The key idea

For each gene, we should not create one messy AF3 job with every possible partner.

Instead, we create **interaction contexts**.

Example:

```text
KCNQ1 interaction contexts:

1. KCNQ1 + KCNE1
2. KCNQ1 + calmodulin + Ca2+
3. KCNQ1 + KCNE1 + calmodulin + Ca2+ if biologically justified
4. KCNQ1 + AKAP9 domain if there is clear interaction-domain evidence
```

Then every KCNQ1 variant gets mapped across all KCNQ1-relevant contexts.

So the model does not only know:

```text
gene = KCNQ1
```

It can know:

```text
variant lies at KCNQ1-KCNE1 interface
variant lies near calmodulin interface
variant lies near Ca2+-related regulatory region
variant is not near any modeled interface
variant is in a low-confidence region
variant is outside modeled domain
```

That is the useful part.

## One-to-one vs one-to-many

You’re right: some genes will be simple, others will be complicated.

### One-to-one example

```text
Gene A protein + Partner B
```

One AF3 job may be enough.

### One-to-many example

```text
Gene A protein + Partner B
Gene A protein + Partner C
Gene A protein + Partner D
```

Here, we probably run **separate jobs**, unless there is evidence B, C, and D all bind together at the same time.

### Shared complex example

```text
TNNI3 + TNNT2 + TNNC1 + TPM1 + ACTC1 + Ca2+
```

One AF3 job can cover multiple genes at once.

That is efficient and biologically cleaner.

## How this should look in your data

You should not force everything into one flat table too early. Use a relational structure.

### Table 1: gene panel audit

```text
gene
n_total_variants
n_missense
n_plp
n_vus
n_uniprot_mapped
protein_name
uniprot_accession
protein_length
```

### Table 2: partner evidence table

```text
target_gene
target_protein
partner_gene_or_entity
partner_type
evidence_source
direct_physical_interaction
same_complex_evidence
simultaneous_binding_evidence
domain_or_full_length
manual_review_status
```

### Table 3: AF3 complex/job table

```text
complex_id
job_name
genes_covered
protein_entities
ion_entities
ligand_entities
full_length_or_domain
why_this_complex
af3_status
output_path
```

### Table 4: residue-level AF3 features

```text
complex_id
gene
chain_id
residue_position
residue_aa
local_plddt
min_distance_to_partner
interface_5A
contact_count
min_distance_to_ion
local_interface_pae
chain_pair_iptm
```

### Table 5: variant-level AF3 features

This is the one that merges into CatBoost.

```text
variant_id
gene
hgvsp
protein_position
ref_aa
alt_aa
af3_complex_available
af3_n_complexes_covering_residue
af3_any_interface_5A
af3_min_distance_to_any_partner
af3_best_partner_gene
af3_near_any_ion_5A
af3_local_plddt
af3_local_interface_pae
af3_feature_missing_reason
```

## How CatBoost should use it

Your final CatBoost table should still be **one row per variant**.

But behind that one row, the AF3 module can summarize multiple complexes.

For example, a KCNQ1 variant might have:

```text
af3_n_complexes_covering_residue = 3
af3_any_interface_5A = true
af3_best_partner_gene = KCNE1
af3_min_distance_to_any_partner = 3.6
af3_near_any_ion_5A = false
af3_local_plddt = 91.2
af3_local_interface_pae = 4.8
```

That is how you keep the training data comprehensive without making it chaotic.

## The most important rule

For each gene, Codex should classify partners into:

```text
1. same-complex partners
2. alternative pairwise partners
3. ligand/ion/cofactor entities
4. pathway-only associations
5. weak/uncertain partners
```

Only categories 1–3 should become AF3 jobs.

Category 4 and 5 should not become production AF3 jobs unless manually approved.

## Does AF3 know whether they all bind together?

No. Not in the way we need.

AF3 can predict a structure for the entities you give it. But it does not replace biological judgment.

So we decide:

```text
Do these proteins belong in one real biological complex?
```

Then AF3 predicts that complex.

If we give AF3 unrelated or mutually exclusive partners, it may still generate a structure, but the features could be biologically misleading.

## The clean final framing

What you are building is:

```text
A full-panel cardiogenetic protein-interaction annotation layer.
```

Not:

```text
An AF3 mutation-effect predictor.
```

The manuscript logic becomes:

> For each cardiogenetic gene, we identified curated protein, ligand, and ion interaction contexts. Wild-type complexes were modeled with AlphaFold 3. Variants were then mapped onto these complexes to generate interaction-aware features, including interface proximity, partner-specific contacts, ligand/ion distance, local confidence, and interface-confidence metrics. These features were evaluated as additions to the existing CatBoost pathogenicity-prioritization model.

That is comprehensive, defensible, and much stronger than simply saying “we ran AlphaFold 3.”

## AF3 Pilot-v1 Pathway Created

Date: 2026-07-03.

We decided not to wait for every difficult full-panel protein partner problem before testing whether AF3 features help the model. The incomplete/domain-heavy genes are being kept aside for later expansion. The first AF3 model pathway now focuses on genes that are either already draft-packet ready or close enough to create review packets.

This is additive only. No existing CatBoost matrix, previous AlphaFold 2 feature table, FoldX table, dbNSFP table, ClinVar/eMERGE/HiRO/CardioBoost data, or trained model output was overwritten.

### Pilot-v1 gene set

The first AF3 pilot set contains 23 target genes:

```text
ACTC1, BAG3, CALM1, CALM2, CALM3, CSRP3, DES, EMD, KCNE2, KCNH2,
KCNJ2, LAMP2, LMNA, PLN, SCN1B, SCN5A, SNTA1, TCAP, TMEM43,
TNNI3, TNNT2, TPM1, TRPM4
```

Status:

| Category | Genes | Meaning |
|---|---:|---|
| Already draft packet ready | 11 | First-wave AF3 packets already exist and need user review before server submission |
| Second-wave packet candidates | 12 | Review packets exist, but final citation/evidence check is still needed before submission |
| Total pilot-v1 target genes | 23 | These can support the first AF3-augmented feature experiment |

The phrase “11 + 12” means 23 target genes, not 32. Some genes share the same AF3 complex, so the number of jobs is smaller than the number of genes.

### Pilot-v1 AF3 complex plan

The pilot uses 14 complex/job rows:

| Stage | Jobs | Notes |
|---|---:|---|
| First-wave draft packets | 4 | Existing draft packets: SCN5A-CALM, troponin/thin filament, LMNA-EMD, PLN-ATP2A2 |
| Second-wave review packets | 10 | New review packets for BAG3, CSRP3/TCAP, DES, KCNH2/KCNE2, KCNJ2, LAMP2, SCN1B, SNTA1, TMEM43, TRPM4 |

First-wave draft packets:

```text
AF3_JOB_001_SCN5A_CALM
AF3_JOB_002_TROPONIN_THIN_FILAMENT_CA
AF3_JOB_003_LMNA_EMD
AF3_JOB_004_PLN_ATP2A2
```

Second-wave review packets:

```text
AF3_PILOTV1_REVIEW_BAG3_HSPB8_HSPA8_STUB1
AF3_PILOTV1_REVIEW_CSRP3_TCAP_ACTN2
AF3_PILOTV1_REVIEW_DES_CRYAB
AF3_PILOTV1_REVIEW_KCNH2_KCNE2
AF3_PILOTV1_REVIEW_KCNJ2_HOMOTETRAMER
AF3_PILOTV1_REVIEW_LAMP2_HSPA8
AF3_PILOTV1_REVIEW_SCN5A_SCN1B
AF3_PILOTV1_REVIEW_SCN5A_SNTA1
AF3_PILOTV1_REVIEW_TMEM43_EMD_LMNA
AF3_PILOTV1_REVIEW_TRPM4_HOMOTETRAMER
```

The second-wave folders are review packets, not approved submission packets. They contain sequences and evidence checklists so the user can inspect whether each complex is biologically acceptable before using AlphaFold Server.

### Files created for the pilot pathway

Scripts:

```text
datasets/feature_sources/alphafold3/scripts/build_af3_pilot_v1_pathway.py
datasets/feature_sources/alphafold3/scripts/build_af3_pilot_v1_variant_scaffold.py
```

Pilot result tables:

```text
results/af3/pilot_v1/af3_pilot_v1_gene_set.tsv
results/af3/pilot_v1/af3_pilot_v1_complex_plan.tsv
results/af3/pilot_v1/af3_pilot_v1_pathway.summary.json
results/af3/pilot_v1/af3_pilot_v1_variant_scaffold.summary.json
```

Submission and output-tracking files:

```text
datasets/feature_sources/alphafold3/jobs/pilot_v1/af3_pilot_v1_submission_manifest.tsv
datasets/feature_sources/alphafold3/jobs/pilot_v1/af3_pilot_v1_completed_job_manifest_template.tsv
datasets/feature_sources/alphafold3/jobs/pilot_v1/manual_submission_packets/
datasets/feature_sources/alphafold3/raw/pilot_v1_server_outputs/
```

Feature-design files:

```text
datasets/feature_sources/alphafold3/processed/pilot_v1/af3_pilot_v1_feature_schema.tsv
datasets/feature_sources/alphafold3/processed/pilot_v1/af3_pilot_v1_variant_feature_scaffold.tsv
```

### Variant coverage before AF3 output

The current frozen modeling table has 85,677 rows. The pilot-v1 AF3 scaffold covers:

| Quantity | Count |
|---|---:|
| Total rows in current modeling matrix | 85,677 |
| Rows in pilot-v1 AF3 genes | 22,119 |
| Covered rows with protein position available | 22,116 |
| Covered rows missing protein position | 3 |

Covered labels:

| Label | Rows |
|---|---:|
| Benign | 8,245 |
| Pathogenic | 2,709 |
| VUS | 11,053 |
| Missing/blank label | 112 |

This means the first AF3 pilot is large enough to test as a real feature family. It is not just a small demonstration layer. The most important limitation is that the actual AF3 structural features still require AlphaFold Server outputs.

### Planned AF3 features

The pilot feature schema includes:

```text
af3_pilot_v1_gene_covered
af3_pilot_v1_complex_count_covering_gene
af3_pilot_v1_complex_ids
af3_residue_mapped_to_model
af3_local_plddt
af3_min_distance_to_any_partner_A
af3_any_partner_interface_5A
af3_any_partner_interface_8A
af3_best_partner_entity
af3_min_distance_to_ion_or_ligand_A
af3_near_ion_or_ligand_5A
af3_local_interface_pae
af3_interface_confidence_class
af3_feature_missing_reason
```

### Next step

Step-by-step from here:

1. User reviews the 4 first-wave draft packets and decides which to submit to AlphaFold Server.
2. User reviews the 10 second-wave evidence packets. Only biologically defensible complexes should be submitted.
3. After each AF3 result is downloaded, place the zip or extracted output under:

```text
datasets/feature_sources/alphafold3/raw/pilot_v1_server_outputs/
```

4. Fill or update:

```text
datasets/feature_sources/alphafold3/jobs/pilot_v1/af3_pilot_v1_completed_job_manifest_template.tsv
```

5. Codex will parse the returned AF3 structures, confirm chain mapping, map variants to residues, compute interface/ligand/confidence features, and write a standalone AF3 feature table.
6. Only after QC will we create a new AF3-augmented model matrix and train a new model version. The current primary CatBoost model remains preserved as the baseline comparator.

## HPA Evidence And AlphaFold Server JSON Added

Date: 2026-07-03.

We added Human Protein Atlas as a localization/expression evidence layer for the AF3 branch. This does not replace UniProt, Reactome, IntAct, or RCSB/PDB. It adds biological context for where a protein is observed and whether the proposed interaction is at least compatible with known cellular localization/expression.

Important rule:

```text
HPA is context evidence, not direct protein-protein interaction evidence.
```

Direct/structural interaction support still comes from:

```text
UniProt SUBUNIT / curated interaction notes
IntAct molecular interaction records
RCSB/PDB co-structure support
Reactome pathway/complex context
```

HPA contributes:

```text
subcellular location
subcellular main location
subcellular additional location
protein tissue specificity/distribution
RNA tissue specificity/distribution
protein class
biological process
molecular function
secretome annotation when present
```

### HPA data added

Downloaded file:

```text
datasets/feature_sources/alphafold3/raw/hpa_evidence/proteinatlas.tsv.zip
```

Generated files:

```text
datasets/feature_sources/alphafold3/interim/evidence/hpa_evidence_by_gene.tsv
results/af3/af3_full_panel_partner_manifest.final_evidence_stack.tsv
results/af3/af3_hpa_evidence.summary.json
```

Summary:

| Quantity | Count |
|---|---:|
| HPA rows loaded | 20,151 |
| AF3 evidence genes checked | 59 |
| AF3 evidence genes with HPA row | 59 |
| AF3 evidence genes with HPA subcellular annotation | 40 |

Some genes have HPA expression/tissue data but no HPA subcellular IF call. Those remain useful as expression/localization context, but they should not be overinterpreted.

### Evidence/readiness tables regenerated

After adding HPA, the AF3 evidence stack was regenerated.

Updated generated files:

```text
results/af3/af3_partner_combined_evidence_review_queue.tsv
results/af3/af3_complex_combined_evidence_readiness.tsv
results/af3/af3_combined_evidence_readiness.summary.json
results/af3/af3_complex_design_review.tsv
results/af3/af3_complex_groups.tsv
results/af3/af3_gene_level_coverage_review.tsv
results/af3/af3_gene_level_coverage_review_v2.tsv
results/af3/pilot_v1/af3_pilot_v1_gene_set.tsv
results/af3/pilot_v1/af3_pilot_v1_complex_plan.tsv
datasets/feature_sources/alphafold3/processed/pilot_v1/af3_pilot_v1_variant_feature_scaffold.tsv
```

The HPA columns are now present in the final pair evidence stack and combined evidence review queue. HPA did not change the direct interaction classification rules.

### AlphaFold Server JSON generated

I verified the AlphaFold Server JSON format from DeepMind's official server JSON documentation/example. The Server format uses:

```text
top-level list of job dictionaries
name
modelSeeds
sequences
dialect = alphafoldserver
version = 1
proteinChain / ion / ligand entities
```

This is different from the standalone open-source AlphaFold 3 JSON dialect.

Generated combined JSON:

```text
datasets/feature_sources/alphafold3/jobs/pilot_v1/server_json/af3_pilot_v1_all_jobs_alphafoldserver.json
```

Generated individual JSON directory:

```text
datasets/feature_sources/alphafold3/jobs/pilot_v1/server_json/individual_jobs/
```

Generated individual JSON files:

```text
AF3_JOB_001_SCN5A_CALM.alphafoldserver.json
AF3_JOB_002_TROPONIN_THIN_FILAMENT_CA.alphafoldserver.json
AF3_JOB_003_LMNA_EMD.alphafoldserver.json
AF3_JOB_004_PLN_ATP2A2.alphafoldserver.json
AF3_PILOTV1_REVIEW_BAG3_HSPB8_HSPA8_STUB1.alphafoldserver.json
AF3_PILOTV1_REVIEW_CSRP3_TCAP_ACTN2.alphafoldserver.json
AF3_PILOTV1_REVIEW_DES_CRYAB.alphafoldserver.json
AF3_PILOTV1_REVIEW_KCNH2_KCNE2.alphafoldserver.json
AF3_PILOTV1_REVIEW_KCNJ2_HOMOTETRAMER.alphafoldserver.json
AF3_PILOTV1_REVIEW_LAMP2_HSPA8.alphafoldserver.json
AF3_PILOTV1_REVIEW_SCN5A_SCN1B.alphafoldserver.json
AF3_PILOTV1_REVIEW_SCN5A_SNTA1.alphafoldserver.json
AF3_PILOTV1_REVIEW_TMEM43_EMD_LMNA.alphafoldserver.json
AF3_PILOTV1_REVIEW_TRPM4_HOMOTETRAMER.alphafoldserver.json
```

The JSON files validated as parseable JSON locally. They are still marked as draft upload inputs because the second-wave jobs require final user evidence review before submission.

### Testable pilot genes now

Detailed table:

```text
results/af3/pilot_v1/af3_pilot_v1_testable_genes_and_interactions.tsv
```

This table has one row per pilot gene with:

```text
gene
pilot stage
complex group
job ID
protein entities
nonprotein entities
one-to-one vs one-to-many/multichain design
submission status
HPA localization/expression fields
rationale/review note
individual AlphaFold Server JSON file path
```

Current testable pilot genes remain 23:

```text
ACTC1, BAG3, CALM1, CALM2, CALM3, CSRP3, DES, EMD, KCNE2, KCNH2,
KCNJ2, LAMP2, LMNA, PLN, SCN1B, SCN5A, SNTA1, TCAP, TMEM43,
TNNI3, TNNT2, TPM1, TRPM4
```

### What is left before using AF3 in training

For the 23 pilot genes:

| Step | Status | Notes |
|---|---|---|
| HPA localization/expression added | Done | 59 AF3 evidence genes have HPA rows; 40 have HPA subcellular annotations |
| Final evidence stack regenerated | Done | HPA included as context evidence |
| AlphaFold Server JSON generated | Done | 14 individual JSONs plus 1 combined JSON |
| First-wave packet review | Needed | 4 jobs are closest to submit-ready |
| Second-wave evidence review | Needed | 10 jobs have JSONs, but should not be submitted blindly |
| AF3 server submission | Not done | User must upload JSON/jobs to AlphaFold Server |
| AF3 output download | Not done | Place results under `raw/pilot_v1_server_outputs/` |
| Chain mapping/QC parser | Next after outputs | Needed to map returned chains/residues correctly |
| Variant-level AF3 feature extraction | Waiting on outputs | Interface distance, ligand/ion proximity, pLDDT, PAE |
| AF3-augmented model training | Waiting on extracted features | New model version only; current baseline model preserved |

Recommended practical upload order:

```text
1. AF3_JOB_001_SCN5A_CALM
2. AF3_JOB_002_TROPONIN_THIN_FILAMENT_CA
3. AF3_JOB_003_LMNA_EMD
4. AF3_JOB_004_PLN_ATP2A2
```

Then review and submit the second-wave jobs if their evidence checklists look biologically acceptable.

## First AF3 Server Output QC: TMEM43-EMD-LMNA

Date: 2026-07-04.

The first returned AlphaFold Server result was added here:

```text
datasets/feature_sources/alphafold3/raw/pilot_v1_server_outputs/fold_af3_pilotv1_review_tmem43_emd_lmna/
```

I inspected the full output folder and generated a formal QC report:

```text
results/af3/pilot_v1/output_qc/fold_af3_pilotv1_review_tmem43_emd_lmna/AF3_OUTPUT_QC_REPORT.md
results/af3/pilot_v1/output_qc/fold_af3_pilotv1_review_tmem43_emd_lmna/af3_output_qc_summary.json
results/af3/pilot_v1/output_qc/fold_af3_pilotv1_review_tmem43_emd_lmna/chain_confidence_summary.tsv
results/af3/pilot_v1/output_qc/fold_af3_pilotv1_review_tmem43_emd_lmna/pair_interaction_confidence.tsv
```

### File completeness

The output is complete enough for the AF3 parsing pipeline.

| File type | Count |
|---|---:|
| Model CIF structures | 5 |
| Full-data JSON files | 5 |
| Summary-confidence JSON files | 5 |
| Job request JSON | 1 |

The result contains the pieces needed for feature extraction:

```text
CIF coordinates
per-atom/per-residue confidence through CIF B-factor/pLDDT fields
full-data JSON with PAE
full-data JSON with contact probabilities
summary confidence JSON with ranking score, pTM, ipTM, and chain-pair metrics
job request JSON for chain/entity mapping
```

### Chain mapping

The AlphaFold Server job request confirms:

| Chain | Gene | Submitted length |
|---|---|---:|
| A | TMEM43 | 400 |
| B | EMD | 254 |
| C | LMNA | 664 |

### Best model confidence

| Metric | Value |
|---|---:|
| ranking_score | 0.46 |
| ipTM | 0.14 |
| pTM | 0.32 |
| fraction_disordered | 0.56 |
| has_clash | 0.0 |

### Chain-local confidence

| Chain | Gene | Residues | Mean pLDDT | Median pLDDT | Residues pLDDT <50 | Residues pLDDT >=70 |
|---|---|---:|---:|---:|---:|---:|
| A | TMEM43 | 400 | 79.79 | 83.35 | 15 | 357 |
| B | EMD | 254 | 46.29 | 35.36 | 177 | 57 |
| C | LMNA | 664 | 69.69 | 80.26 | 166 | 458 |

Interpretation: TMEM43 is locally modeled fairly well. LMNA has many confident regions but also many low-confidence residues. EMD is weak in this full-length three-chain run.

### Inter-chain confidence

| Gene pair | Min PAE | Median PAE | Max contact probability | Mean contact probability |
|---|---:|---:|---:|---:|
| TMEM43-EMD | 24.8 | 31.4 | 0.12 | 0.000205 |
| TMEM43-LMNA | 28.5 | 31.6 | 0.13 | 0.000112 |
| EMD-LMNA | 27.9 | 31.7 | 0.10 | 0.000078 |

This is the important result: the model does not support a reliable three-chain interface. ipTM is low, inter-chain PAE is high, and contact probabilities are near zero.

### Variant coverage if usable

These three genes cover 3,688 current modeling rows:

| Gene | Rows | Benign | Pathogenic | VUS | With protein position |
|---|---:|---:|---:|---:|---:|
| EMD | 568 | 265 | 91 | 212 | 568 |
| LMNA | 2,073 | 624 | 529 | 912 | 2,073 |
| TMEM43 | 1,047 | 450 | 2 | 594 | 1,047 |

### Decision for this output

This output is useful as a pipeline smoke test and confirms that the AlphaFold Server download contains the exact files we need to parse future AF3 results.

However, this specific TMEM43-EMD-LMNA trimer should not be used as a high-confidence interface model. It should be marked:

```text
af3_output_parseable = true
af3_local_structure_context_possible = partial
af3_interface_features_trusted = false
af3_interface_confidence_class = low
af3_missing_or_low_confidence_reason = low_iptm_high_interchain_pae_near_zero_contact_probability
```

Recommended next steps:

1. Run the simpler first-wave `LMNA-EMD` job and compare confidence.
2. For TMEM43, consider smaller pairwise jobs: `TMEM43-EMD` and `TMEM43-LMNA`.
3. If pairwise full-length jobs remain weak, move to domain/fragment designs.
4. Do not let this low-confidence trimer inject strong interface-distance features into the model.

## AF3 Pilot V1 10-Job Batch QC And Feature Extraction

Date: 2026-07-04.

The first single TMEM43-EMD-LMNA output was superseded by a full 10-job pilot batch from AlphaFold Server. The batch is stored here:

```text
datasets/feature_sources/alphafold3/raw/pilot_v1_server_outputs/folds_2026_07_04_21_37/
```

I created a general batch analyzer:

```text
datasets/feature_sources/alphafold3/scripts/analyze_af3_pilot_v1_batch_outputs.py
```

This script does not modify the prior registry, feature matrix, or trained CatBoost models. It only reads AlphaFold Server output folders and writes AF3-specific QC and feature tables.

### Outputs generated

```text
results/af3/pilot_v1/batch_output_qc/folds_2026_07_04_21_37/AF3_BATCH_OUTPUT_QC_REPORT.md
results/af3/pilot_v1/batch_output_qc/folds_2026_07_04_21_37/af3_batch_output_qc_summary.json
results/af3/pilot_v1/batch_output_qc/folds_2026_07_04_21_37/af3_batch_job_qc_summary.tsv
results/af3/pilot_v1/batch_output_qc/folds_2026_07_04_21_37/af3_batch_chain_confidence.tsv
results/af3/pilot_v1/batch_output_qc/folds_2026_07_04_21_37/af3_batch_pair_interaction_qc.tsv
results/af3/pilot_v1/batch_output_qc/folds_2026_07_04_21_37/af3_batch_gene_feature_coverage.tsv
datasets/feature_sources/alphafold3/processed/pilot_v1/af3_pilot_v1_batch_variant_features_long.tsv
```

The long variant-feature table has one row per variant-position per AF3 chain/job. This is deliberate. For genes tested in more than one AF3 context, such as `SCN5A`, we preserve separate AF3 evidence rows and defer aggregation until the modeling step.

### Batch-level result

| Metric | Result |
|---|---:|
| AF3 jobs analyzed | 10 |
| Chain confidence rows | 27 |
| Pairwise interaction rows | 28 |
| Long variant-feature rows | 28,580 |
| Jobs with at least one usable interface | 4 |
| Jobs not interface-usable | 6 |

The 10 jobs were:

```text
AF3_JOB_001_SCN5A_CALM
AF3_JOB_002_TROPONIN_THIN_FILAMENT_CA
AF3_JOB_003_LMNA_EMD
AF3_JOB_004_PLN_ATP2A2
AF3_PILOTV1_REVIEW_CSRP3_TCAP_ACTN2
AF3_PILOTV1_REVIEW_DES_CRYAB
AF3_PILOTV1_REVIEW_KCNH2_KCNE2
AF3_PILOTV1_REVIEW_KCNJ2_HOMOTETRAMER
AF3_PILOTV1_REVIEW_SCN5A_SCN1B
AF3_PILOTV1_REVIEW_TMEM43_EMD_LMNA
```

### Job-level calls

| Job | Protein entities | Ranking | ipTM | pTM | Call |
|---|---|---:|---:|---:|---|
| AF3_JOB_001_SCN5A_CALM | SCN5A,CALM1 | 0.86 | 0.77 | 0.66 | usable_complex_features |
| AF3_JOB_002_TROPONIN_THIN_FILAMENT_CA | ACTC1,TPM1,TNNI3,TNNT2,TNNC1 | 0.51 | 0.25 | 0.36 | usable only for the TNNI3-TNNC1 local interface |
| AF3_JOB_003_LMNA_EMD | LMNA,EMD | 0.47 | 0.10 | 0.23 | not_interface_usable |
| AF3_JOB_004_PLN_ATP2A2 | PLN,ATP2A2 | 0.35 | 0.23 | 0.68 | not_interface_usable |
| AF3_PILOTV1_REVIEW_CSRP3_TCAP_ACTN2 | CSRP3,TCAP,ACTN2 | 0.34 | 0.22 | 0.39 | not_interface_usable |
| AF3_PILOTV1_REVIEW_DES_CRYAB | DES,CRYAB | 0.30 | 0.09 | 0.22 | not_interface_usable |
| AF3_PILOTV1_REVIEW_KCNH2_KCNE2 | KCNH2,KCNE2 | 0.58 | 0.31 | 0.40 | not_interface_usable |
| AF3_PILOTV1_REVIEW_KCNJ2_HOMOTETRAMER | KCNJ2,KCNJ2,KCNJ2,KCNJ2 | 0.88 | 0.75 | 0.76 | usable_complex_features |
| AF3_PILOTV1_REVIEW_SCN5A_SCN1B | SCN5A,SCN1B | 0.93 | 0.82 | 0.70 | usable_complex_features |
| AF3_PILOTV1_REVIEW_TMEM43_EMD_LMNA | TMEM43,EMD,LMNA | 0.46 | 0.14 | 0.32 | not_interface_usable |

### Supported interfaces

The supported or partially supported AF3 interfaces are:

| Job | Pair | Pair ipTM | Min PAE | Max contact probability | Interface residues within 5 A | Support |
|---|---|---:|---:|---:|---:|---|
| AF3_JOB_001_SCN5A_CALM | SCN5A-CALM1 | 0.77 | 3.6 | 0.86 | 98 | strong |
| AF3_JOB_002_TROPONIN_THIN_FILAMENT_CA | TNNI3-TNNC1 | 0.75 | 3.4 | 0.92 | 134 | strong |
| AF3_PILOTV1_REVIEW_KCNJ2_HOMOTETRAMER | KCNJ2-KCNJ2 | 0.68-0.71 | 1.5-2.3 | 0.88-1.00 | 6-190 | strong/moderate across tetramer interfaces |
| AF3_PILOTV1_REVIEW_SCN5A_SCN1B | SCN5A-SCN1B | 0.82 | 1.8 | 0.98 | 110 | strong |

Interpretation: these are the first AF3 complex outputs that look usable for real interface-aware features. The most useful immediate feature families are:

```text
AF3 residue pLDDT
distance from variant residue to nearest partner chain
binary within-interface flags at 5 A and 8 A
nearest partner gene
pair-level ipTM
pair-level min PAE
pair-level max contact probability
pair support class
```

### Low-confidence or failed interface jobs

The following jobs should not be used as high-confidence interface evidence in the current form:

```text
AF3_JOB_003_LMNA_EMD
AF3_JOB_004_PLN_ATP2A2
AF3_PILOTV1_REVIEW_CSRP3_TCAP_ACTN2
AF3_PILOTV1_REVIEW_DES_CRYAB
AF3_PILOTV1_REVIEW_KCNH2_KCNE2
AF3_PILOTV1_REVIEW_TMEM43_EMD_LMNA
```

These outputs can still be useful as negative/diagnostic AF3 runs, but they should be flagged:

```text
af3_interface_features_trusted = false
af3_interface_confidence_class = low
```

They may need smaller domain/fragment designs, membrane-aware interpretation, or better-defined stoichiometry before they become useful for interface features.

### Gene-level AF3 feature coverage from this batch

Rows with residue features are rows where the variant had a protein position that mapped onto the returned AF3 chain.

| Gene | AF3 job | Feature rows | Rows with residue feature | Rows with strong/moderate pair | Within 5 A | Within 8 A |
|---|---|---:|---:|---:|---:|---:|
| SCN5A | SCN5A-CALM | 4,603 | 3,827 | 4,603 | 86 | 147 |
| SCN5A | SCN5A-SCN1B | 4,603 | 3,827 | 4,603 | 92 | 206 |
| CALM1 | SCN5A-CALM | 179 | 92 | 179 | 35 | 58 |
| CALM2 | SCN5A-CALM | 222 | 93 | 222 | 33 | 64 |
| CALM3 | SCN5A-CALM | 200 | 106 | 200 | 42 | 81 |
| TNNI3 | Troponin thin filament | 769 | 540 | 769 | 271 | 430 |
| KCNJ2 | KCNJ2 homotetramer | 2,444 | 2,044 | 2,444 | 930 | 1,284 |
| SCN1B | SCN5A-SCN1B | 636 | 321 | 636 | 82 | 153 |

Other genes in this batch have extracted residue-local fields but should not yet be treated as high-confidence interface features because the pair-level AF3 confidence did not support the submitted complex.

### Important modeling decision

Do not blindly merge all AF3 distances as if every complex is equally real.

Recommended first AF3 modeling layer:

1. Use all `af3_residue_plddt` values as local AF3 confidence/context features where residue mapping is available.
2. Use interface-distance and partner-contact features only when `af3_best_pair_support_class` is `strong_supported_interface` or `moderate_supported_interface`.
3. Keep low-confidence job flags as missing or low-confidence, not as true negative biology.
4. Aggregate long AF3 rows into variant-level features only after deciding whether to use best-supported complex, nearest partner, or multiple-complex summaries.

## AF3 Pilot V2 Second-20 Server Jobs

Date: 2026-07-04.

After the first 10-job AF3 batch, I generated the next 20 jobs as a separate additive batch. The first 10 results directly informed the design:

- Do not repeat weak full-length interfaces unchanged.
- Keep successful design patterns, especially `SCN5A` partner jobs and channel homotetramers.
- Use pair-only or reduced-complex versions for large or ambiguous complexes.
- Keep ligand/CCD review jobs protein-only first unless the nonprotein entity is simple calcium.

Script:

```text
datasets/feature_sources/alphafold3/scripts/build_af3_second20_server_json.py
```

Outputs:

```text
datasets/feature_sources/alphafold3/jobs/pilot_v2_second20/server_json/af3_pilot_v2_second20_alphafoldserver.json
datasets/feature_sources/alphafold3/jobs/pilot_v2_second20/server_json/individual_jobs/
results/af3/pilot_v2_second20/af3_pilot_v2_second20_job_summary.tsv
results/af3/pilot_v2_second20/af3_pilot_v2_partner_source_summary.tsv
results/af3/pilot_v2_second20/af3_pilot_v2_second20.summary.json
```

The combined server JSON contains 20 jobs and passed JSON validation.

### Second-20 batch summary

| Category | Jobs | Why included |
|---|---:|---|
| Second-wave ready-for-evidence-check jobs | 4 | Already generated from curated partner decisions but not included in the first completed 10 outputs |
| Manual cardiac whitelist core complexes | 2 | Major missing cardiac channel complexes: IKs and L-type calcium channel |
| Pair-only simplifications | 6 | Safer reduced jobs for complex systems; helps determine whether a specific pair is clearer than a large multichain job |
| Protein-only ligand-review jobs | 4 | Avoids uncertain ligand CCD handling while still testing the protein complex |
| Cardiac-specificity review job | 1 | FHOD3-SQSTM1 has UniProt support but needs caution |
| Replacement after weak first-10 job | 2 | KCNH2-KCNE2 and DES-CRYAB were weak, so test simpler homomeric contexts |
| Reduced desmosome multichain job | 1 | Captures plaque context without submitting the full oversized desmosome core |

### Second-20 jobs

| Order | Job | Target genes | Protein entities | Nonprotein | Token estimate | Reason |
|---:|---|---|---|---|---:|---|
| 1 | AF3_PILOTV2_001_BAG3_HSPB8_HSPA8_STUB1 | BAG3 | BAG3,HSPB8,HSPA8,STUB1 |  | 1,720 | CASA/chaperone complex |
| 2 | AF3_PILOTV2_002_LAMP2_HSPA8 | LAMP2 | LAMP2,HSPA8 |  | 1,056 | Chaperone-mediated autophagy context |
| 3 | AF3_PILOTV2_003_SCN5A_SNTA1 | SNTA1,SCN5A | SCN5A,SNTA1 |  | 2,521 | SCN5A partner extension after SCN5A-CALM and SCN5A-SCN1B succeeded |
| 4 | AF3_PILOTV2_004_TRPM4_HOMOTETRAMER | TRPM4 | TRPM4,TRPM4,TRPM4,TRPM4 |  | 4,856 | Channel homotetramer, following KCNJ2 homotetramer success |
| 5 | AF3_PILOTV2_005_KCNQ1_KCNE1_CALM1_CA | KCNQ1,KCNE1,CALM1,CALM2,CALM3 | KCNQ1,KCNE1,CALM1 | CA | 955 | Core IKs/calmodulin/calcium context |
| 6 | AF3_PILOTV2_006_KCNQ1_KCNE1 | KCNQ1,KCNE1 | KCNQ1,KCNE1 |  | 805 | Pair-only IKs subunit test |
| 7 | AF3_PILOTV2_007_KCNQ1_CALM1_CA | KCNQ1,CALM1,CALM2,CALM3 | KCNQ1,CALM1 | CA | 826 | Pair-only KCNQ1-calmodulin test |
| 8 | AF3_PILOTV2_008_CACNA1C_CACNB2_CALM1_CA | CACNA1C,CACNB2,CALM1,CALM2,CALM3 | CACNA1C,CACNB2,CALM1 | CA | 3,031 | L-type calcium-channel context |
| 9 | AF3_PILOTV2_009_CACNA1C_CACNB2 | CACNA1C,CACNB2 | CACNA1C,CACNB2 |  | 2,881 | Pair-only calcium-channel subunit test |
| 10 | AF3_PILOTV2_010_CASQ2_TRDN | CASQ2,TRDN | CASQ2,TRDN |  | 1,128 | Compact calcium-release pair; avoids full-length RYR2 problem |
| 11 | AF3_PILOTV2_011_GLA_HOMODIMER | GLA | GLA,GLA |  | 858 | Protein-only homodimer before substrate/ligand review |
| 12 | AF3_PILOTV2_012_HCN4_HOMOTETRAMER | HCN4 | HCN4,HCN4,HCN4,HCN4 |  | 4,812 | Protein-only channel tetramer before cAMP CCD review |
| 13 | AF3_PILOTV2_013_PRKAA2_PRKAB2_PRKAG2 | PRKAG2 | PRKAA2,PRKAB2,PRKAG2 |  | 1,393 | AMPK heterotrimer, protein-only before AMP ligand review |
| 14 | AF3_PILOTV2_014_SLC22A5_PDZK1 | SLC22A5 | SLC22A5,PDZK1 |  | 1,076 | Transporter scaffold context, protein-only before L-carnitine review |
| 15 | AF3_PILOTV2_015_FHOD3_SQSTM1 | FHOD3 | FHOD3,SQSTM1 |  | 1,862 | Review-grade sarcomere/autophagy partner |
| 16 | AF3_PILOTV2_016_KCNH2_HOMOTETRAMER | KCNH2 | KCNH2,KCNH2,KCNH2,KCNH2 |  | 4,636 | Replacement after weak KCNH2-KCNE2 first-10 output |
| 17 | AF3_PILOTV2_017_DES_HOMODIMER | DES | DES,DES |  | 940 | Replacement after weak DES-CRYAB first-10 output |
| 18 | AF3_PILOTV2_018_DSP_JUP | DSP,JUP | DSP,JUP |  | 3,616 | Desmosome linker pair |
| 19 | AF3_PILOTV2_019_DSG2_DSC2 | DSG2,DSC2 | DSG2,DSC2 |  | 2,019 | Desmosome cadherin pair |
| 20 | AF3_PILOTV2_020_PKP2_DSP_JUP | PKP2,DSP,JUP | PKP2,DSP,JUP |  | 4,497 | Reduced desmosome plaque context |

### Target-gene coverage in second-20 batch

The second-20 batch covers 26 target genes:

```text
BAG3, CACNA1C, CACNB2, CALM1, CALM2, CALM3, CASQ2, DES, DSC2,
DSG2, DSP, FHOD3, GLA, HCN4, JUP, KCNE1, KCNH2, KCNQ1, LAMP2,
PKP2, PRKAG2, SCN5A, SLC22A5, SNTA1, TRDN, TRPM4
```

### Cautions before upload

The following jobs are larger and should be watched for server acceptance and output confidence:

```text
AF3_PILOTV2_004_TRPM4_HOMOTETRAMER
AF3_PILOTV2_012_HCN4_HOMOTETRAMER
AF3_PILOTV2_016_KCNH2_HOMOTETRAMER
AF3_PILOTV2_020_PKP2_DSP_JUP
```

The following are intentionally protein-only first passes because ligand/CCD handling still needs review:

```text
AF3_PILOTV2_011_GLA_HOMODIMER
AF3_PILOTV2_012_HCN4_HOMOTETRAMER
AF3_PILOTV2_013_PRKAA2_PRKAB2_PRKAG2
AF3_PILOTV2_014_SLC22A5_PDZK1
```

This second-20 batch is ready for AlphaFold Server upload as a test batch. After output returns, analyze it with the same pair-level QC logic used for the first 10.

## AF3 Additional Returned Jobs: 2026-07-04 22:36 Batch

Date analyzed: 2026-07-04.

New downloaded AF3 Server output folder:

```text
results/af3/folds_2026_07_04_22_36/
```

The same batch also exists in Downloads:

```text
/Users/dheyaa/Downloads/folds_2026_07_04_22_36/
```

This batch contains 5 returned jobs:

```text
AF3_JOB_001_SCN5A_CALM
AF3_PILOTV2_017_DES_HOMODIMER
AF3_PILOTV2_018_DSP_JUP
AF3_PILOTV2_019_DSG2_DSC2
AF3_PILOTV2_020_PKP2_DSP_JUP
```

I updated the AF3 batch parser so it can analyze arbitrary batch folders and read both pilot-v1 and pilot-v2 manifests:

```text
datasets/feature_sources/alphafold3/scripts/analyze_af3_pilot_v1_batch_outputs.py
```

Command used:

```text
python3 datasets/feature_sources/alphafold3/scripts/analyze_af3_pilot_v1_batch_outputs.py results/af3/folds_2026_07_04_22_36
```

Outputs:

```text
results/af3/batch_output_qc/folds_2026_07_04_22_36/AF3_BATCH_OUTPUT_QC_REPORT.md
results/af3/batch_output_qc/folds_2026_07_04_22_36/af3_batch_job_qc_summary.tsv
results/af3/batch_output_qc/folds_2026_07_04_22_36/af3_batch_pair_interaction_qc.tsv
results/af3/batch_output_qc/folds_2026_07_04_22_36/af3_batch_gene_feature_coverage.tsv
datasets/feature_sources/alphafold3/processed/folds_2026_07_04_22_36/af3_folds_2026_07_04_22_36_variant_features_long.tsv
```

### Batch summary

| Metric | Result |
|---|---:|
| Jobs analyzed | 5 |
| Chain rows | 11 |
| Pairwise interaction rows | 7 |
| Long variant-feature rows | 26,916 |
| Jobs with usable strong/moderate interface | 2 |
| Jobs not interface-usable | 3 |

### Job-level calls

| Job | Proteins | Ranking | ipTM | pTM | Call |
|---|---|---:|---:|---:|---|
| AF3_JOB_001_SCN5A_CALM | SCN5A,CALM1 | 0.83 | 0.74 | 0.66 | usable complex features |
| AF3_PILOTV2_017_DES_HOMODIMER | DES,DES | 0.64 | 0.38 | 0.38 | not interface-usable |
| AF3_PILOTV2_018_DSP_JUP | DSP,JUP | 0.34 | 0.28 | 0.34 | not interface-usable |
| AF3_PILOTV2_019_DSG2_DSC2 | DSG2,DSC2 | 0.43 | 0.30 | 0.38 | not interface-usable |
| AF3_PILOTV2_020_PKP2_DSP_JUP | PKP2,DSP,JUP | 0.54 | 0.47 | 0.40 | usable complex features |

### Supported interactions in this batch

| Job | Pair | Pair ipTM | Min PAE | Max contact | Support |
|---|---|---:|---:|---:|---|
| AF3_JOB_001_SCN5A_CALM | SCN5A-CALM1 | 0.74 | 4.1 | 0.84 | strong |
| AF3_PILOTV2_020_PKP2_DSP_JUP | DSP-JUP | 0.59 | 5.8 | 0.90 | moderate |

Interpretation:

- The repeated `SCN5A-CALM1` job again supports a strong interface, consistent with the earlier first-10 result.
- The `PKP2-DSP-JUP` reduced desmosome trimer supports the `DSP-JUP` pair at moderate confidence.
- The standalone `DSP-JUP` pair did not pass the interface-confidence gate, so the trimer context appears more useful than the pair-only design for this interaction.
- `DES-DES`, `DSP-JUP` pair-only, and `DSG2-DSC2` show some close atomic contacts, but pair confidence is too weak for trusted interface features under the current thresholds.

### Combined completed AF3 summary

Combined completed summaries across the first batch and this new batch:

```text
results/af3/current_completed_af3_summary/completed_af3_jobs.tsv
results/af3/current_completed_af3_summary/completed_af3_pair_qc.tsv
results/af3/current_completed_af3_summary/completed_af3_gene_feature_coverage.tsv
```

Current completed AF3 output status:

| Item | Count |
|---|---:|
| Completed output rows across batches | 15 |
| Unique completed job IDs | 14 |
| Duplicate completed job | SCN5A-CALM |
| Strong/moderate supported pair rows | 11 |

Trusted or potentially usable supported interactions so far:

```text
SCN5A-CALM1
SCN5A-SCN1B
TNNI3-TNNC1
KCNJ2-KCNJ2 homotetramer interfaces
DSP-JUP in PKP2-DSP-JUP trimer
```

Important caveat: `SCN5A-CALM1` appears in two completed batches. Treat it as replicated support for the same interaction, not as a new independent gene coverage category.

## Systematic AF3 Coverage Policy

Date: 2026-07-04.

The AF3 branch should not become an ad hoc list of attractive complexes. For the manuscript, we need a reproducible coverage policy that explains:

```text
which genes were considered
which partners were considered
which complex was submitted
which AF3 outputs were trusted
which features entered the model
which citations support each partner decision
```

### Current gene universe

The current included modeling table contains:

| Item | Count |
|---|---:|
| Included genes | 53 |
| Included variant rows | 84,982 |

Gene-priority audit table:

```text
results/af3/af3_gene_priority_for_systematic_coverage.tsv
```

This table ranks all included genes by current variant-row burden and records current AF3 status.

Top included genes by current variant rows:

| Rank | Gene | Included rows | Current AF3 status |
|---:|---|---:|---|
| 1 | RYR2 | 9,465 | not yet planned/completed; needs domain-level design |
| 2 | FLNC | 5,671 | not yet planned/completed; needs domain-level design |
| 3 | DSP | 5,333 | pilot-v2 desmosome jobs planned |
| 4 | MYH7 | 5,228 | not yet planned/completed; needs domain-level design |
| 5 | SCN5A | 4,603 | pilot-v1 completed; pilot-v2 SNTA1 planned |
| 6 | MYBPC3 | 4,236 | not yet planned/completed; needs domain-level design |
| 7 | CACNA1C | 3,557 | pilot-v2 jobs planned |
| 8 | KCNH2 | 3,478 | pilot-v1 weak KCNH2-KCNE2 completed; pilot-v2 homotetramer planned |
| 9 | KCNQ1 | 2,792 | pilot-v2 jobs planned |
| 10 | HCN4 | 2,170 | pilot-v2 homotetramer planned |

Rare/sparse included genes are also important:

```text
PLN, KCNE2, CALM1, CALM2, CALM3, KCNE1, TCAP, KRAS, RIT1, FHOD3
```

These genes should not be ignored just because they have fewer rows. They often represent exactly the cases where protein context may help, but they should not drive the first full engineering pass.

### Recommended strategy

Use all 53 included genes as the formal review universe, but prioritize execution in tiers.

| Tier | Genes | Goal | AF3 handling |
|---|---|---|---|
| Tier 1 | Top 30 by variant burden | Maximize feature coverage and manuscript impact | Run one defensible primary complex per gene where possible |
| Tier 2 | Rare/sparse but clinically important genes | Test whether AF3 helps sparse-gene interpretation | Run compact, high-confidence complexes if available |
| Tier 3 | Large/domain-problem genes | Avoid bad full-length AF3 jobs | Do domain-level designs only after residue/domain mapping |
| Tier 4 | Genes without curated partners | Avoid weak evidence | Keep in partner audit table, do not submit until evidence improves |

This means the practical target is:

```text
all 53 genes reviewed
top 30 genes actively covered first where biologically feasible
rare genes covered when the partner/complex is strong and compact
large genes handled by domain design rather than full-length AF3
```

This is better than choosing only top 30 or only rare 30. It gives us a manuscript-safe reason for every gene:

- high-burden genes improve feature coverage;
- rare genes test the sparse-gene hypothesis;
- skipped genes are documented with a reason, not forgotten.

### Partner universe

The current evidence stack contains:

| Item | Count |
|---|---:|
| Partner/evidence rows | 172 |
| Target genes represented | 56 |
| Broad complex groups | 10 |

Publication-style partner audit table:

```text
results/af3/af3_partner_manifest_publication_audit.tsv
```

This table is the start of the manuscript supplement table. It keeps:

```text
target gene
target UniProt accession
target protein length
candidate partner/entity
partner UniProt accession or identifier
interaction source
evidence type
UniProt pair support
IntAct support
Reactome shared pathway or complex IDs
RCSB/PDB support
PubMed IDs
recommended AF3 input mode
manual review status
final include/exclusion flag
source URLs
```

### Should we submit one partner or all partners?

For each gene, we should keep all plausible partners in the audit table, but submit only a small number of AF3 jobs.

Recommended rule:

| Level | Meaning | Use |
|---|---|---|
| Candidate partner list | All plausible partners from UniProt, IntAct, Reactome, RCSB/PDB, HPA/literature | Supplement/audit table |
| Primary AF3 complex | Best-supported, most cardiac-relevant complex for that gene | Main AF3 feature source |
| Secondary AF3 complex | Additional partner only if biologically distinct and likely useful | Sensitivity/extra features |
| Diagnostic/replacement job | Used after weak AF3 result or uncertain stoichiometry | Helps refine design, not automatic model input |

So the goal is not:

```text
every possible partner goes into AF3
```

The goal is:

```text
every gene has a documented partner review,
each gene gets one primary AF3 complex when feasible,
some genes get secondary complexes when there is a strong reason,
all omitted partners have a recorded reason.
```

### How to avoid missing important partners

For each gene, partner review should use this source order:

1. UniProt reviewed entry: canonical sequence, function, domains, subunit notes, curated interactions, PubMed links.
2. IntAct: curated molecular interaction evidence.
3. Reactome: pathway and complex context.
4. RCSB/PDB: structural support for the protein or interaction.
5. Human Protein Atlas: localization and tissue/subcellular context.
6. STRING or broad pathway association only as weak backup, not as primary evidence.
7. Targeted literature review for genes where the automated evidence stack is incomplete or clinically important.

The publication table should eventually have one row per gene-partner candidate, not only one row per submitted AF3 job.

Suggested supplement columns:

```text
gene
protein_name
uniprot_accession
candidate_partner
partner_type
partner_uniprot_or_identifier
complex_or_pathway_name
evidence_sources
pubmed_ids
reactome_ids
intact_ids
pdb_ids
hpa_localization
af3_submission_status
af3_job_id
af3_output_status
af3_pair_iptm
af3_pair_min_pae
af3_max_contact_probability
af3_interface_support_class
used_for_model_features
exclusion_or_low_confidence_reason
```

### How AF3 outputs become model features

AF3 output is not automatically accepted just because the job was submitted.

Use this rule:

| AF3 result | Feature decision |
|---|---|
| Strong or moderate pair confidence | Use interface-distance/contact features |
| Low pair confidence but good local pLDDT | Use local pLDDT only, mark interface features missing/low-confidence |
| Low pLDDT and low interface confidence | Do not use as trusted structural feature |
| Weak full-length job in important gene | Redesign as pair-only, homomer, or domain-level job |

Current pair-level acceptance criteria:

```text
strong_supported_interface:
  pair ipTM >= 0.70
  pair min PAE <= 5
  max contact probability >= 0.60

moderate_supported_interface:
  pair ipTM >= 0.45
  pair min PAE <= 10
  max contact probability >= 0.35
```

Only strong/moderate outputs should feed interface-distance features into the model.

### Current interpretation after first 10 jobs

Strong usable first-batch interactions:

```text
SCN5A-CALM1
SCN5A-SCN1B
KCNJ2 homotetramer
TNNI3-TNNC1 inside the thin-filament job
```

Weak or not interface-usable first-batch outputs:

```text
LMNA-EMD
PLN-ATP2A2
CSRP3-TCAP-ACTN2
DES-CRYAB
KCNH2-KCNE2
TMEM43-EMD-LMNA
```

The second-20 batch was therefore built to:

- add missing high-burden/channel genes;
- add already curated second-wave jobs;
- add pair-only simplifications;
- replace weak KCNH2-KCNE2 with KCNH2 homotetramer;
- replace weak DES-CRYAB with DES homodimer;
- keep ligand-review jobs protein-only first.

### Next systematic tasks

1. Review the top-30 gene list and mark each gene as `covered`, `planned`, `domain_design_needed`, or `no_curated_partner`.
2. Build a domain-design queue for `RYR2`, `FLNC`, `MYH7`, `MYBPC3`, and large desmosome/sarcomere proteins.
3. After the second-20 outputs return, run the same pair-level AF3 QC.
4. Promote only strong/moderate AF3 interfaces to model-ready features.
5. Keep all candidate partners and citations in the publication audit table, including partners not submitted.

## AF3Draft Smoke-Test Layer: July 4, 2026

Important naming decision: the current AlphaFold 3 work is now treated as **AF3Draft**. This keeps the prototype structural layer separate from the established CatBoost model and avoids implying that the AF3 partner set is final.

No previous model tables were overwritten. New draft outputs were written under:

```text
results/af3Draft/
datasets/feature_sources/af3Draft/processed/
results/model_performance/af3Draft_smoke_test/
```

The current AF3Draft features were built from the completed pilot outputs only:

```text
datasets/feature_sources/alphafold3/processed/pilot_v1/af3_pilot_v1_batch_variant_features_long.tsv
datasets/feature_sources/alphafold3/processed/folds_2026_07_04_22_36/af3_folds_2026_07_04_22_36_variant_features_long.tsv
```

Aggregated AF3Draft feature table:

```text
datasets/feature_sources/af3Draft/processed/af3Draft_variant_features_aggregated.tsv
```

Augmented smoke-test modeling table:

```text
datasets/modeling/interim/modeling_table_with_splits_weights_plus_af3Draft.tsv
```

### Current AF3Draft variant coverage

| Scope | Rows | Genes | Any AF3Draft output | Strong/moderate AF3Draft pair | Strong AF3Draft pair |
|---|---:|---:|---:|---:|---:|
| All modeling rows | 85,677 | 61 | 32,143 | 13,884 | 7,220 |
| Primary binary rows | 42,990 | 53 | 15,557 | 6,762 | 3,479 |
| VUS rows | 42,361 | 61 | 16,461 | 7,059 | 3,684 |

Interpretation: AF3Draft already touches many rows because some completed jobs are high-burden genes. However, only the strong/moderate rows should be considered plausible interface-feature coverage. The rest are either local-residue coverage or low-confidence interface attempts.

More precise structural coverage:

| Scope | Rows | Trusted AF3Draft pair context | Trusted pair plus residue feature | Trusted pair plus 5A interface | Trusted pair plus 8A interface |
|---|---:|---:|---:|---:|---:|
| All modeling rows | 85,677 | 13,884 | 11,048 | 1,723 | 3,075 |
| Primary binary rows | 42,990 | 6,762 | 4,650 | 732 | 1,276 |
| VUS rows | 42,361 | 7,059 | 6,339 | 986 | 1,792 |

The trusted-pair count means the variant's gene has a strong/moderate AF3Draft interaction context. The residue-feature and interface columns are stricter variant-level structural coverage.

### Genes currently driving trusted AF3Draft coverage

Strong or moderate AF3Draft pair coverage currently comes mainly from:

```text
DSP, SCN5A, JUP, TNNI3, SCN1B, KCNJ2, CALM1, CALM2, CALM3
```

This is not yet a systematic 53-gene AF3 layer. It is a trial run of the pathway.

### CatBoost smoke test

Two AF3Draft-informed CatBoost models were trained using the same primary binary CatBoost setup and the same `split_source_heldout` split:

| Model | AF3Draft feature policy | Purpose |
|---|---|---|
| `af3Draft_plus_primary_binary_fullDraft` | Uses structural features plus draft/provenance fields such as AF3 status/counts/nearest partner gene | Exploratory only |
| `af3Draft_plus_primary_binary_strictStructural` | Uses only structure-like AF3Draft features: pLDDT, interface support rank, distance/contact flags, ipTM, PAE, contact probability | Cleaner reviewer-safe smoke test |

AF3 job IDs and complex names were excluded from both models.

Main report:

```text
results/model_performance/af3Draft_smoke_test/AF3DRAFT_SMOKE_TEST_REPORT.md
```

### Smoke-test result

The strict structural AF3Draft model did not materially change the already strong primary binary model. That is scientifically acceptable at this stage because AF3Draft coverage is incomplete and external AF3-covered subsets are small.

Compared with the current baseline primary binary model, the strict structural AF3Draft model showed:

| External subset | AUROC change | AUPRC change | Threshold metrics |
|---|---:|---:|---|
| CardioBoost | +0.0027 | +0.0015 | Sensitivity, specificity, PPV unchanged |
| eMERGE | -0.0025 | -0.0046 | Sensitivity, specificity, PPV unchanged |
| HiRO | +0.0027 | +0.0030 | Sensitivity, specificity, PPV unchanged |

On AF3Draft trusted-pair external rows, sample sizes are still small:

```text
CardioBoost: 35 rows
eMERGE: 42 rows
HiRO: 21 rows
```

So these subset results are useful for debugging, but they are not yet publication-grade evidence.

### Interpretation

This smoke test supports continuing the AF3Draft pathway, but it is not yet the smoking-gun result. The honest conclusion is:

```text
The current AF3Draft layer is technically integrated, does not destabilize model performance,
and shows small nonzero signal. The next step is to scale AF3 coverage systematically across
the 53 included genes, then rerun the same strict structural-only comparison.
```

The next AF3 runs should prioritize high-burden genes and complexes likely to produce trusted interfaces, especially:

```text
RYR2, FLNC, MYH7, MYBPC3, KCNQ1, CACNA1C, HCN4, PKP2, DSG2, DSC2, BAG3, TRDN, PRKAG2, GLA
```
