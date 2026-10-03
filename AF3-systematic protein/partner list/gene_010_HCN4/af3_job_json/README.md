# HCN4 AF3 Job JSONs

Default upload file:

```text
AF3_SYS_gene_010_HCN4_TOP3_upload_alphafoldserver.json
```

Submit the three jobs in this packet only by default:

1. Full-length HCN4 homotetramer.
2. HCN4 C-linker/CNBD tetramer with cAMP ligand draft (`CCD_CMP`).
3. HCN4 channel-core/C-linker/CNBD tetramer.

Important: the full-length tetramer is 4,812 protein residues, close to practical server limits. The cAMP job uses RCSB/PDB CCD code `CMP` as `CCD_CMP`; if AlphaFold Server rejects it, run the same protein-only CNBD tetramer and document ligand rejection.
