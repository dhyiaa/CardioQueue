# Point E: Protein Structure, Folding, and Stability Features

This feature block covers AlphaFold structure context, residue confidence,
structural neighborhood features, and mutation stability features such as
delta-delta-G. These are feature-engineering sources only. They must not be used
as labels.

## Core Principle

Separate these concepts:

| Concept | What it measures | Example features | Main source/tool |
|---|---|---|---|
| AlphaFold confidence/context | Whether a residue sits in a confident predicted structure | pLDDT, local confidence bin, PAE availability, model version | AlphaFold DB |
| Structural location | Where the residue is in the folded protein | secondary structure, solvent accessibility, interface/proximity flags, domain membership | AlphaFold model plus DSSP/FreeSASA/domain annotations |
| Stability effect | Whether an amino-acid substitution destabilizes the folded state | predicted delta-delta-G, destabilizing flag, absolute ddG magnitude | FoldX/Rosetta/PremPS/DynaMut-style tools |
| Protein language/constraint context | Learned or evolutionary intolerance of substitution | ESM1b/ESM-variant score, AlphaMissense, conservation | ESM/AlphaMissense/dbNSFP |

AlphaFold alone does not provide delta-delta-G. Delta-delta-G requires a
separate mutation-stability method.

## Current Public AlphaFold Access

| Resource | Status | Notes |
|---|---|---|
| AlphaFold DB targeted prediction API | API-ready | Example endpoint: `https://alphafold.ebi.ac.uk/api/prediction/P12883` |
| Human AlphaFold DB proteome tar | Downloaded and inspected | `UP000005640_9606_HUMAN_v6.tar`, 5,177,506,304 bytes |
| AlphaFold model files | API-linked | PDB/CIF URLs returned per UniProt accession |
| AlphaFold confidence JSON | API-linked | Per-residue pLDDT available through `confidence_v6.json` URLs |
| AlphaFold PAE JSON | API-linked | Pairwise predicted aligned error available through `predicted_aligned_error_v6.json` URLs |
| AlphaMissense hg38 annotations linked from AlphaFold entries | Already handled under Point C | Keep pathogenicity scores in Point C; structure links may be cross-referenced here |
| UniProt reviewed human FASTA | Downloaded | 20,431 reviewed human proteins, UniProt release `2026_02` |
| UniProt reviewed human feature TSV | Downloaded | 20,431 reviewed human protein records plus header; includes domains/regions/motifs/binding/active-site fields |

## Feature Set for ACMG-Oriented Modeling

### Minimum Useful Features

These are the first features to add because they are practical and defensible:

| Feature | Level | Reason |
|---|---|---|
| `uniprot_accession` | variant/gene | Required for protein-level joins |
| `protein_position` | variant | Required for residue-level features |
| `alphafold_model_version` | gene/protein | Reproducibility |
| `alphafold_global_plddt` | protein | Overall model confidence |
| `alphafold_residue_plddt` | residue | Whether residue-level structural interpretation is trustworthy |
| `alphafold_residue_confidence_bin` | residue | Very low / low / confident / very high |
| `in_confident_af_region` | residue | Do not over-interpret low-confidence/disordered regions |
| `alphafold_has_pae` | protein | Whether PAE can support domain/interface context |
| `alphafold_status` and missing reason | row | Prevent silent nulls |

### High-Value Structural Features

Add after exact protein mapping is stable:

| Feature | Tool/source | Notes |
|---|---|---|
| Secondary structure at residue | DSSP/mkdssp on AlphaFold model | Helix/sheet/coil context |
| Relative solvent accessibility | FreeSASA or DSSP | Buried substitutions are often more structurally disruptive |
| Buried/core flag | FreeSASA threshold | Useful for missense impact |
| Domain membership | UniProt/Pfam/InterPro | Important for cardiomyopathy/arrhythmia genes with functional domains |
| Distance to active/functional site | UniProt features or curated domain annotations | Use only when source is explicit |
| Distance to known pathogenic residues | Internal ClinVar/ClinGen mapped residues | Useful feature, but must avoid leakage in training splits if derived from labels |
| Coiled-coil/repeat region | UniProt/InterPro | Important for myosin, desmosomal, and cytoskeletal proteins |
| Interface/proximity flag | AlphaFold multimer/known complex/PDB if available | Use carefully; monomer AlphaFold may not represent interfaces |

### Delta-Delta-G / Stability Features

Delta-delta-G features should be added for missense variants only after a
stability engine is selected and reproducibly installed.

| Feature | Meaning |
|---|---|
| `ddg_method` | FoldX, Rosetta, PremPS, DynaMut2, or other chosen method |
| `ddg_kcal_mol` | Predicted mutation-induced folding free-energy change |
| `ddg_abs` | Absolute predicted perturbation |
| `ddg_destabilizing_flag` | Method-specific threshold, documented before modeling |
| `ddg_strong_destabilizing_flag` | Higher threshold for severe destabilization |
| `ddg_status` | `ok`, `not_applicable`, `not_queryable`, `failed` |
| `ddg_missing_reason` | Exact reason no score is present |

Do not mix delta-delta-G from different engines without a `ddg_method` feature,
because sign conventions and calibration differ by method.

## Practical Tool Decision

| Tool/source | Role | Current decision |
|---|---|---|
| AlphaFold DB API | Protein and residue confidence metadata | Implement targeted workflow first |
| AlphaFold human proteome bulk tar | Full local structural archive | Downloaded and parsed into a manifest |
| FoldX/FoldXPro | Fast local ddG candidate | Academic access obtained; FoldXPro Mac Silicon runs locally and one-variant BuildModel smoke test passed |
| Rosetta cartesian_ddg | More computationally expensive ddG candidate | Recommended sensitivity/validation engine for a smaller high-priority subset |
| DynaMut2/PremPS/DDGun-style tools | Useful for spot checks and comparison | Avoid for large automated cohorts unless bulk/API terms and reproducibility are explicit |
| DSSP/mkdssp | Secondary structure and accessibility | Installed locally and smoke-tested on MYH7 AlphaFold PDB |
| FreeSASA | Solvent accessibility | Installed locally through Python API and smoke-tested on MYH7 AlphaFold PDB |
| UniProt/Pfam/InterPro | Domains and functional sites | Add after UniProt mapping is stable |

### Delta-Delta-G Tool Review and Recommendation

For this cardiogenetics CatBoost study, the best delta-delta-G choice is not
the most elaborate physics engine by default. The practical requirement is a
reproducible, local, scalable feature generator for thousands of missense
variants across AlphaFold structures, with clear missingness and confidence
flags.

| Tool | Strengths | Limitations for this study | Recommendation |
|---|---|---|---|
| FoldX/FoldXPro | Fast local executable; designed for protein stability and mutation effects; works directly from PDB-like structures; practical for thousands of variants | Licensed academic software; PDB-only input; AlphaFold structures should be filtered by pLDDT and mapping QC; record whether the engine is FoldX or FoldXPro | Use as primary first-pass ddG engine |
| Rosetta `cartesian_ddg` | Strong, widely used structure-modeling framework; allows local backbone/cartesian refinement around mutations; useful as a robustness check | More computationally expensive; Rosetta energy units and protocol choices require careful calibration; heavier setup | Use on a smaller sensitivity subset, not as the first whole-cohort engine |
| PremPS | Published ML method for single-point stability changes; designed for large numbers of mutations; has a public implementation/server | Requires structure input and dependency review; external web use may be less reproducible than local locked tooling unless standalone workflow is pinned | Consider as secondary comparator after FoldX/Rosetta decision |
| DynaMut2 | Web/API-oriented stability and dynamics prediction; useful for individual variants and spot checks | Server/API queueing and terms may limit high-throughput use; external service adds reproducibility risk | Spot-check or small validation set only unless bulk use is explicitly allowed |
| DDGun / sequence-based tools | Can work when structure mapping is weak; useful fallback for low-confidence or unresolved structures | Measures a different signal than structure-relaxed ddG; can overlap conceptually with sequence predictors already planned | Optional fallback/comparator, not primary structure ddG |

Recommended implementation order:

1. Finish DSSP/FreeSASA residue-level extraction and add those columns to
   `protein_features.tsv`.
2. Resolve multi-fragment AlphaFold mapping for large proteins such as `TTN`,
   `DSP`, and `RYR2`.
3. Add FoldX/FoldXPro columns as the primary `ddg_method` feature block for a small pilot batch.
4. Scale FoldXPro to eligible high-confidence missense variants after pilot QC.
5. Run Rosetta `cartesian_ddg` only on a prioritized sensitivity subset, such
   as high-confidence P/LP vs B/LB missense variants in key genes, to test
   whether FoldX-derived stability signals are directionally robust.

## Implemented Targeted Workflow

Script:

```bash
python3 datasets/feature_sources/protein_structure/scripts/alphafold_api.py doctor

python3 datasets/feature_sources/protein_structure/scripts/alphafold_api.py gene-panel \
  --gene-panel datasets/feature_sources/gnomad_v4/external/gene_panel.example.txt \
  --output datasets/feature_sources/protein_structure/interim/alphafold_gene_panel.example.tsv \
  --raw-json datasets/feature_sources/protein_structure/raw/alphafold/alphafold_gene_panel.example.json

python3 datasets/feature_sources/protein_structure/scripts/alphafold_api.py validate-gene-panel \
  --input datasets/feature_sources/protein_structure/interim/alphafold_gene_panel.example.tsv
```

The gene-panel workflow records reviewed human UniProt accession, AlphaFold
model version, global pLDDT, pLDDT confidence fractions, model URLs, confidence
JSON URL, PAE JSON URL, and AlphaMissense hg38 URL.

Residue-level workflow expects already mapped protein positions:

```bash
python3 datasets/feature_sources/protein_structure/scripts/alphafold_api.py residue-features \
  --input datasets/feature_sources/protein_structure/external/residue_positions.example.tsv \
  --output datasets/feature_sources/protein_structure/interim/alphafold_residue_features.example.tsv
```

Required input columns:

- `input_id`
- `uniprot_accession`
- `protein_position`

## Mapping Requirements

Protein-structure features are only trustworthy if the variant maps correctly.

Required mapping checks:

1. Genomic variant must map to a transcript.
2. Transcript must map to a protein sequence.
3. Protein change must match the UniProt sequence at the wild-type residue.
4. UniProt isoform/canonical accession must be recorded.
5. Failed or ambiguous mappings must be explicit `not_queryable`, not blank.

## Status

| Component | Status |
|---|---|
| AlphaFold API access | Verified |
| UniProt reviewed-human mapping | Verified for seed genes |
| AlphaFold gene-panel metadata script | Added and validated: 19 seed genes, 16 `ok`, 3 `not_found` |
| AlphaFold bulk tar | Downloaded and inspected: 47,172 compressed PDB/CIF model files |
| AlphaFold bulk manifest annotation | Added and validated: 19 seed genes, 19 `ok` |
| Residue pLDDT extraction | Added and validated on two example positions |
| Unified protein feature table | Added: 3,591 rows across HiRO, eMERGE, and CardioBoost |
| Structural features from PDB/CIF | Tooling installed; residue-level table merge pending |
| Delta-delta-G | FoldXPro smoke test passed; pilot batch/table integration pending. Rosetta remains a future sensitivity subset option |
| Domain/functional-site annotation | Added from UniProt features in unified protein table |

Validated example outputs:

| Output | Path | Result |
|---|---|---|
| AlphaFold gene-panel features | `datasets/feature_sources/protein_structure/interim/alphafold_gene_panel.example.tsv` | 19 rows; 16 `ok`, 3 `not_found` |
| AlphaFold raw metadata cache | `datasets/feature_sources/protein_structure/raw/alphafold/alphafold_gene_panel.example.json` | Saved; 24 KB |
| AlphaFold human bulk tar | `datasets/feature_sources/protein_structure/raw/alphafold/UP000005640_9606_HUMAN_v6.tar` | Downloaded; 5,177,506,304 bytes |
| AlphaFold bulk manifest | `datasets/feature_sources/protein_structure/raw/alphafold/UP000005640_9606_HUMAN_v6.manifest.txt` | 47,172 members |
| AlphaFold bulk summary | `datasets/feature_sources/protein_structure/raw/alphafold/UP000005640_9606_HUMAN_v6.summary.json` | 47,172 `.gz` model files |
| AlphaFold bulk seed-panel coverage | `datasets/feature_sources/protein_structure/interim/alphafold_bulk_manifest_gene_panel.example.tsv` | 19 rows; 19 `ok` |
| Residue-position input example | `datasets/feature_sources/protein_structure/external/residue_positions.example.tsv` | 2 rows |
| Residue pLDDT features | `datasets/feature_sources/protein_structure/interim/alphafold_residue_features.example.tsv` | 2 rows; both `ok` |
| UniProt reviewed human FASTA | `datasets/feature_sources/protein_structure/raw/uniprot_human_reviewed_2026_02.fasta` | 20,431 sequences |
| UniProt reviewed human feature TSV | `datasets/feature_sources/protein_structure/raw/uniprot_human_reviewed_features_2026_02.tsv` | 20,432 lines including header |
| SHA256 checksums | `datasets/feature_sources/protein_structure/raw/SHA256SUMS` | Local checksums recorded |
| Unified protein features | `datasets/feature_sources/protein_structure/interim/protein_features.tsv` | 3,591 rows; protein mapping, UniProt features, and AlphaFold pLDDT where resolvable |
| DSSP/FreeSASA environment | `.tools/envs/protein-structure` | `mkdssp 4.6.1`; FreeSASA Python API available |
| FoldXPro environment | `.tools/foldx/` | Native Mac Silicon binary runs; `BuildModel` command listed |
| FoldXPro smoke test | `datasets/feature_sources/protein_structure/interim/foldxpro_smoke_test.tsv` | ACTC1 M134T; `ddg_kcal_mol=-0.168657`; status `ok` |

The three seed genes without AlphaFold DB model coverage in the targeted API
example were `TTN`, `DSP`, and `RYR2`. Their UniProt mappings were found, and
the local bulk tar does contain segmented models for them:

| Gene | UniProt | Bulk fragments |
|---|---|---:|
| `TTN` | `Q8WZ42` | 166 |
| `DSP` | `P15924` | 9 |
| `RYR2` | `Q92736` | 19 |

For large proteins, the bulk tar is therefore better than relying only on the
targeted AlphaFold prediction API.

## QC Rules

- Do not compute structure features for synonymous, intronic, UTR, or splice-only
  variants unless there is a defined protein consequence.
- Missense variants are the primary target for residue and delta-delta-G
  features.
- Truncating variants may use gene/domain/region features, but not missense ddG.
- Low pLDDT residues should be flagged; do not interpret them as stable folded
  structural sites.
- Every output must include status and missing-reason fields.
