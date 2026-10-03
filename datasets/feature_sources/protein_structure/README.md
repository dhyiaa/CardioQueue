# Protein Structure Sources

Role: protein mapping, AlphaFold, UniProt, and structural feature injection.

Current status:

- AlphaFold human v6 bulk tar is local.
- UniProt reviewed human FASTA and feature TSV are local.
- Example gene-panel and residue-level AlphaFold outputs generated.
- Unified protein feature table generated from HiRO, eMERGE, and CardioBoost baseline inputs.
- DSSP/mkdssp and FreeSASA are installed in a project-local micromamba environment; batch residue-level features were generated from local AlphaFold structures.
- FoldX 5.1 is installed locally. The accepted development/external matrix now has an outcome-blind maximal protein-mapping audit, and corrected `RepairPDB` plus batched `BuildModel` scale-up is in progress.

Folders:

- `raw/`: AlphaFold, UniProt, and checksums.
- `external/`: small residue/gene input examples.
- `interim/`: generated example outputs.
- `scripts/`: AlphaFold/UniProt helper scripts.

Tool environment:

```bash
.tools/micromamba/micromamba run -p .tools/envs/protein-structure mkdssp --version
.tools/micromamba/micromamba run -p .tools/envs/protein-structure python -c "import freesasa; print(hasattr(freesasa, 'calc'))"
.tools/foldx/bin/foldxpro --help
```

Primary generated table:

- `interim/protein_features.tsv`: joinable protein feature table with parsed protein positions, UniProt domain/region/motif/site overlaps, and AlphaFold residue pLDDT where the AlphaFold fragment can be resolved.
- `interim/protein_features.summary.json`: row counts and status counts for the generated table.
- `interim/dssp_freesasa_features.tsv`: batch DSSP secondary-structure and FreeSASA solvent-accessibility features for mapped AlphaFold residue rows.
- `interim/dssp_freesasa_features.summary.json`: DSSP/FreeSASA batch coverage and status counts.
- `../../modeling/interim/protein_structure_selected_features.tsv`: selected one-row-per-variant protein/structure features for modeling.
- `../../modeling/interim/final_modeling_table_with_clingen_protein.tsv`: current model-facing table with ClinGen plus protein/structure features merged.
- `interim/foldxpro_smoke_test.tsv`: ACTC1 M134T FoldXPro smoke test; `ddg_kcal_mol=-0.168657`.
- `interim/foldxpro_ddg_pilot.tsv`: small FoldXPro DDG pilot on clean mapped substitutions.
- `interim/foldxpro_ddg_pilot.summary.json`: pilot status counts and selected label/gene mix.
- `interim/foldxpro_ddg_batch_preflight.tsv`: preflight-selected primary FoldX batch rows using `pLDDT >= 70`.
- `interim/foldxpro_ddg_batch.tsv`: resumable primary FoldX DDG batch output, written incrementally while the background run progresses.
- `interim/foldxpro_ddg_batch.summary.json`: live/finished FoldX batch status counts.

Known current limitations:

- Splice/non-protein HGVS variants are intentionally marked `not_mapped`.
- Large multi-fragment AlphaFold proteins are marked `fragment_mapping_unresolved` unless the residue can be safely assigned to the first fragment.
- DSSP/FreeSASA features remain a separate source table, but selected fields
  are merged into the current model-facing table.
- FoldX DDG is defined here only for exact single-amino-acid substitutions. Indels, truncating, splice, synonymous, and unresolved/ambiguous transcript mappings are intentionally excluded.
- The primary FoldX analysis uses reference-validated mappings with AlphaFold residue `pLDDT >= 70`; `pLDDT >= 50` and all mapped residues are sensitivity analyses.
- Run the complete mapping, calculation, matched model controls, and paired evaluation with `scripts/run_maximized_foldx_analysis.sh`. The FoldX cache is resumable and verifies exact mutation keys before reuse.

Maximized FoldX mapping outputs:

- `interim/maximized_foldx_mapping_audit.tsv`: one row per registry allele with mapping method, explicit exclusion reason, fragment position, pLDDT, and tier eligibility.
- `interim/maximized_foldx_structure_features.tsv`: one-row-per-allele structure bridge used by the FoldX runner.
- `interim/maximized_foldx_mapping.summary.json`: split-level coverage and method accounting.
- `interim/maximized_foldx51_ddg.tsv`: final corrected FoldX values after the full runner completes.
- `interim/maximized_foldx_augmentations/`: matched eligibility-only and numeric-DDG feature tables by pLDDT tier.
