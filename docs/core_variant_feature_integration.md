# Core Variant and Feature Integration Map

The project should keep two layers separate:

1. Core variant/label sources.
2. Feature-injection sources used to annotate those variants.

This separation prevents accidental leakage and keeps source provenance visible.

## Core Variant Sources

| Source | Role | Local status | Notes |
|---|---|---|---|
| ClinVar | Public large-scale submitted variant labels | Local, validated | Main broad label source; confidence should use review status and conflicts. |
| eMERGE/Glazer | Public arrhythmia-gene variant labels | Local, profiled | Real training/evaluation source; VUS-heavy three-class data. |
| HiRO 480 | Internal cardiogenetics variants | Not found locally yet | Add once the local file/folder is provided. |
| CardioBoost public | Public inherited-cardiac-condition variant labels | Local, profiled | Coordinate-rich GitHub-derived source; processed cohort is binary B/P with no VUS. |
| CardioBoost V2 / DYNA Zenodo | Protein sequence-pair CardioBoost/ClinVar data | Local, downloaded from Zenodo, MD5-validated, processed | Useful for protein-sequence modeling/pretraining; lacks genomic coordinates/gene symbols. |

## Feature-Injection Sources

| Source | Injected features | Status |
|---|---|---|
| gnomAD v4 | AF, popmax AF, homozygote/hemizygote counts, gene constraint, benign-proxy sampling | API-ready; bulk deferred. |
| In silico predictors | AlphaMissense, CADD, REVEL, SIFT, PolyPhen-2, MetaLR, FATHMM-XF, etc. | Partial; dbNSFP still needed for full old-predictor parity. |
| ClinGen | Gene-disease validity, curated evidence context, variant evidence when available | Partial; gene validity API works, variant-level pipeline pending. |
| AlphaFold/UniProt/protein features | pLDDT, protein position, domains, structural context, stability features | Partial; raw AlphaFold/UniProt local, residue/stability extraction pending. |

## Integration Rules

- Every row in a unified modeling table must keep `source_dataset`.
- Exact variant keys should be normalized before merging:
  `genome_build`, `chrom`, `pos`, `ref`, `alt`.
- Keep label provenance separate from annotation provenance.
- Do not silently overwrite labels when the same variant appears in multiple
  sources; create conflict-resolution fields.
- Use source-aware splits so near-duplicate public variants do not leak between
  training and validation.
- Treat missing phenotype fields from public variant-only sources as unavailable,
  not as normal phenotype values.
- CardioBoost can contribute B/P signal, but its lack of VUS must be modeled
  explicitly in any three-class experiment.
