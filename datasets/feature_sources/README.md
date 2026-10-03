# Feature Sources

These sources do not primarily provide the study's variant labels. Instead, they
inject features into variants from ClinVar, eMERGE, HiRO, and CardioBoost.

| Folder | Source | Main features |
|---|---|---|
| `gnomad_v4/` | gnomAD v4 | Allele frequency, popmax AF, homozygote/hemizygote counts, gene constraint, benign-proxy sampling |
| `insilico/` | AlphaMissense, CADD, dbNSFP/VEP/ESM planning | In silico pathogenicity and conservation scores |
| `clingen/` | ClinGen | Gene-disease validity and curated evidence context |
| `protein_structure/` | AlphaFold, UniProt, structural tools | pLDDT, protein mapping, domains, residue/structure features |

Keep feature-source outputs provenance-rich. Every generated table should record
source status, missing reason, source version/date, and join key.

