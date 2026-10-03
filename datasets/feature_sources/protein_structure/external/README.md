# Protein Structure External Inputs

This directory is reserved for protein-structure input manifests, such as
variant-to-UniProt mappings and protein residue-position tables.

Residue-level AlphaFold features require a table with:

- `input_id`
- `uniprot_accession`
- `protein_position`

Do not populate these fields by guesswork. The wild-type amino acid from the
variant must match the selected UniProt protein sequence before using structural
features for modeling.
