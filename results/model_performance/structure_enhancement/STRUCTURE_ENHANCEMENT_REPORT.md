# Registry-wide protein-structure enhancement

## Coverage

Reference-validated AlphaFold v6, DSSP, and FreeSASA features were available for 2,889 training, 493 validation, and 271 source-held-out variants. Long proteins used overlapping 1,400-residue AlphaFold fragments at 200-residue offsets. The selected fragment maximized distance from its nearest edge. Reference amino acids were checked against both reviewed UniProt sequence and the selected PDB fragment.

## External performance

Across all 430 held-out variants, baseline AUROC was 0.977 and structure-enhanced AUROC was 0.979. The gene-clustered paired difference was 0.003 (95% CI -0.000 to 0.005). Brier score changed from 0.0529 to 0.0520.

Among 271 structure-covered external variants, AUROC changed from 0.984 to 0.988. The gene-clustered paired difference was 0.003 (95% CI 0.000 to 0.006).

## Interpretation

The enhanced model assigned nonzero importance to pLDDT, DSSP, and FreeSASA values. The availability indicator had zero importance. Effect estimates are small and their paired intervals should determine whether the analysis is described as improvement or feasibility. FoldX DDG remains unavailable in development and is not part of the demonstrated incremental signal.
