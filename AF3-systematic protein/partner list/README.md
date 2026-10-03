# AF3 Systematic Protein-Partner List

Created: 2026-07-04

This folder is the clean systematic AF3 workspace requested by Dheyaa. It is separate from earlier AF3 pilot and AF3Draft model outputs.

Execution order:

1. Phase 1: top 30 genes by included variant rows.
2. Phase 2: top 10 rare-end genes by lowest included row count.
3. Phase 3: remaining middle-burden genes.

Current panel size is 53 genes, so the remaining group is 13 genes, not 20. The older phrase "rest 20" was approximate.

Priority table:

```text
00_priority/af3_systematic_gene_priority.tsv
```

Gene #1 packet:

```text
gene_001_RYR2/
```

## Required Per-Gene Folder Structure

Every gene must get its own folder:

```text
gene_###_GENE/
```

Each gene folder should contain:

| File/Folder | Required | Purpose |
|---|---|---|
| `GENE_systematic_review.md` | yes | Human-readable summary, decisions, caveats, and next action |
| `01_gene_context.tsv` | yes | Local dataset burden, labels, variant/source context |
| `02_sequence_audit.tsv` | yes | UniProt protein sequence, length, isoform/PTM/mature-chain caveats |
| `03_candidate_partner_evidence.tsv` | yes | All plausible protein/ion/ligand partners and evidence |
| `04_af3_design_decisions.tsv` | yes | Selected/deferred AF3 designs and why |
| `05_variant_mapping_audit.tsv` | yes | Whether variants can map to protein positions/domains |
| `09_specific_citation_sources.tsv` | yes | Specific citations/sources for each evidence claim |
| `af3_job_json/` | yes | AlphaFold Server JSONs for upload or draft review |
| `af3_job_results/` | yes | User-added AlphaFold Server outputs after jobs finish |

Additional files such as domain-boundary reviews, PDB contact tables, topology
audits, or unresolved-design notes should be added when needed.

## Mandatory Workflow For Every Gene

For each gene, do this in order:

1. Build the systematic table:
   - local variant burden;
   - binary/VUS/pathogenic/benign counts;
   - current AF3 status;
   - current protein-position mapping coverage.
2. Map gene to protein:
   - reviewed UniProt accession;
   - canonical sequence;
   - sequence length;
   - isoform, mature-chain, signal peptide, topology, and PTM caveats.
3. Build the partner list:
   - all plausible partners from UniProt, IntAct, Reactome, RCSB/PDB, Human Protein Atlas, and literature;
   - cardiac relevance class;
   - whether each partner is primary, secondary, deferred, or excluded.
4. Make the decision table:
   - accepted AF3 design;
   - deferred design;
   - reason for each decision;
   - specific missing evidence if not ready.
5. Add specific citations:
   - one row per source/evidence claim in `09_specific_citation_sources.tsv`;
   - include URL, PMID/PDB/UniProt/Reactome IDs, and Ctrl-F terms.
6. Generate AF3 job JSONs:
   - write combined AlphaFold Server JSON;
   - write individual job JSONs;
   - write a JSON manifest with upload status.
7. Create or preserve `af3_job_results/`:
   - this is where Dheyaa will add downloaded AlphaFold Server outputs;
   - future analysis scripts should read results from this folder.

## AF3 JSON Status Labels

Every JSON job must be labeled:

| Status | Meaning |
|---|---|
| `ready_to_upload` | Evidence, sequence, and design are clean enough for server upload |
| `draft_review_needed` | JSON generated, but Dheyaa/manual review is required before upload |
| `not_generated_yet` | Not scientifically safe yet; reason must be recorded |

If a design is unsafe, do not silently skip it. Add it to:

```text
af3_job_json/GENE_not_generated_jobs.tsv
```

## Citation Rule

Every major partner/design decision needs a specific source row.

The source row must include:

```text
source_id
evidence_claim
source_type
citation_or_id
url
ctrl_f_terms
used_for
```

The goal is that someone can open the source and Ctrl-F the exact terms used to
support the decision.

## AF3 Result Folder Rule

Every gene folder must include:

```text
af3_job_results/
```

Dheyaa will add downloaded AlphaFold Server outputs there. Use a dated batch
folder when possible:

```text
af3_job_results/folds_YYYY_MM_DD_HH_MM/
```

Do not mix raw AF3 server outputs with JSON inputs or processed feature tables.


## AF3 Job Limit Policy

For each gene, keep the systematic partner evidence table broad enough to be auditable, but limit actual AlphaFold Server submissions to the top 3 jobs whenever possible. A 4th job is allowed only when it is clearly justified by strong cardiac relevance, nonredundant biology, and available AF3 quota.

Default upload classes:

| Class | Meaning | Action |
|---|---|---|
| `top3_submit_first` | Best three jobs for the gene | Submit first |
| `optional_4th_only_if_quota` | Useful but not essential | Submit only if quota/time allows |
| `hold` | Sensitivity or deferred design | Keep documented, do not submit now |

Selection rules:

1. Prefer cardiac disease mechanism over generic interaction abundance.
2. Prefer curated direct interaction or mapped domain evidence over STRING-style functional association.
3. Prefer compact domain-level jobs over full-length proteins when the full protein is large or modular.
4. Avoid repeating prior AF3 designs that already failed unless the new design changes the biological context.
5. Preserve all non-submitted jobs in the evidence tables and manifests so reviewers can see what was considered.

Current top-priority upload summary:

```text
00_priority/af3_top3_upload_policy_current_genes.tsv
```

## Gene Run Notes
- `gene_004_MYH7`: MYH7 systematic packet created with top-3-only AF3 jobs for motor-ACTC1, neck-MYL2/MYL3, and S2-MYBPC3 contexts.
- `gene_003_DSP`: DSP systematic packet created with N-terminal plaque and C-terminal DES/intermediate-filament draft jobs plus `af3_job_results/` folder.
- `gene_002_FLNC`: FLNC systematic packet created with domain-level AF3 draft jobs and `af3_job_results/` folder for server outputs.


## Gene 005 SCN5A Packet

Created `gene_005_SCN5A` as a systematic AF3 packet. Default upload is the top-3 JSON only: SCN5A-CALM1+Ca4, SCN5A-SCN1B, and SCN5A-SNTA1. CALM1 and SCN1B already have prior strong local AF3 support; SNTA1 is included as a cardiac scaffold candidate but must pass output QC before feature promotion. Held partners are FGF13, ANK3, PKP2, and NEDD4/NEDD4L/WWP2/GPD1L.


## Gene 006 MYBPC3 Packet

Created `gene_006_MYBPC3` as a systematic domain-level AF3 packet. Default upload is the top-3 JSON only: MYBPC3 N-terminal 1-452 with MYH7 S2 838-963, MYBPC3 N-terminal 1-452 with ACTC1, and exploratory MYBPC3 C-terminal 971-1274 with MYH7 S2 838-963. Full-length MYBPC3 and expanded thin-filament jobs are held for later sensitivity analysis.


## Gene 007 CACNA1C Packet

Created `gene_007_CACNA1C` as a systematic AF3 packet. Default upload is the top-3 JSON only: CACNA1C full length with CACNB2/CALM1/Ca4, CACNA1C full length with CACNB2 as a pair control, and CACNA1C IQ region 1648-1705 with CALM1/Ca4. CACNA2D/gamma/STAC and compact AID-CACNB2 sensitivity jobs are documented but held until output QC or boundary review justifies them.


## Gene 008 KCNH2 Packet

Created `gene_008_KCNH2` as a systematic AF3 packet. Default upload is the top-3 JSON only: KCNH2 full-length homotetramer, compact KCNH2 PAS/PAC 1-150 with KCNH2 CNBHD 720-870, and KCNH2-KCNE2 full-length comparator. The prior KCNH2-KCNE2 AF3 output was not interface-usable, so KCNE2 is documented as biologically supported but not the lead feature source.


## Gene 009 KCNQ1 Packet

Created `gene_009_KCNQ1` as a systematic AF3 packet. Default upload is the top-3 JSON only: tetrameric KCNQ1-KCNE1-CALM1-Ca16 IKs regulatory complex, KCNQ1-KCNE1 pair control, and KCNQ1-CALM1-Ca4 pair control. KCNE1 is included for AF3 partner modeling by user decision despite mixed broader modeling scope. KCNE2, AKAP9, PIP2, and compact peptide sensitivity jobs are documented but held.

## Gene 010 HCN4 Packet

Created `gene_010_HCN4` as a systematic AF3 packet. Default upload is the top-3 JSON only: full-length HCN4 homotetramer, compact HCN4 C-linker/CNBD tetramer with cAMP ligand draft, and compact HCN4 channel-core/C-linker/CNBD tetramer. HCN1/2/3 heteromers, kinase/signaling associations, and by-similarity regulatory partners are documented but held.

## Gene 011 PKP2 Packet

Created `gene_011_PKP2` as a systematic AF3 packet. Default upload is the top-3 JSON only: PKP2-DSP N-terminal-JUP reduced desmosomal plaque context, PKP2-JUP pair control, and PKP2-DSP N-terminal pair control. Cadherin-tail, SCN5A/GJA1/ANK3, and generic proteomics partners are documented but held.
- `gene_012_LMNA/`: LMNA systematic packet ready; top-3 JSON covers mature homodimer, LMNA-tail/EMD compact rescue, and IF-rod homodimer. Prior full-length LMNA-EMD and TMEM43-EMD-LMNA outputs are documented as not interface-usable.
- `gene_013_DSG2/`: DSG2 systematic packet ready; top-3 JSON covers extracellular DSG2-DSC2 rescue, DSG2 ectodomain homodimer, and DSG2 cytoplasmic-tail/PKP2 plaque context. Prior full-length DSG2-DSC2 is documented as not interface-usable.
- `gene_014_RBM20/`: RBM20 systematic packet ready as review-grade domain/monomer only; no curated protein partner promoted. RNA/U1/U2 snRNP designs are explicitly held.
- `gene_015_DSC2/`: DSC2 systematic packet ready; top-3 JSON covers extracellular DSC2-DSG2 rescue, DSC2-tail/PKP2, and DSC2-tail/JUP plaque contexts. Prior full-length DSG2-DSC2 is documented as not interface-usable.
