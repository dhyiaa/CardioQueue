# Datasets

This folder contains the project dataset sources that provide variant rows,
labels, or reusable external data packages.

## Layout

| Folder | Dataset | Role | Status |
|---|---|---|---|
| `clinvar/` | ClinVar `variant_summary.txt.gz` | Large public variant-label source | Downloaded, checksum-validated |
| `emerge/full_arrhythmia_gene_dataset/` | eMERGE/Glazer arrhythmia-gene dataset | Core variant-level training/evaluation source | Local, profiled |
| `hiro/full_dataset/` | HiRO / CASPER WES / VERDICT | Internal variant + phenotype source | Local, profiled |
| `cardioboost/public_dataset/` | CardioBoost public coordinate-rich data | Core public variant-level source | Local, profiled; old model artifacts removed |
| `cardioboost/public_dataset/v2_dyna_zenodo_sequence_dataset/` | CardioBoost V2 / DYNA Zenodo sequence-pair data | Auxiliary protein-sequence source | Downloaded, MD5-validated, processed |
| `feature_sources/gnomad_v4/` | gnomAD v4 | Feature injection and benign-proxy source | API-ready; examples validated |
| `feature_sources/insilico/` | AlphaMissense, CADD, dbNSFP/VEP/ESM notes | Feature injection | Partial; AlphaMissense local, CADD API-ready |
| `feature_sources/clingen/` | ClinGen APIs/data | Curated evidence feature source | Gene-validity API validated |
| `feature_sources/protein_structure/` | AlphaFold/UniProt/protein features | Structural feature source | AlphaFold and UniProt local |

## Important Distinction

The coordinate-rich datasets are the ones that can join cleanly into the main
variant table by `chrom-pos-ref-alt` or related variant identifiers:

- ClinVar
- eMERGE
- CardioBoost public dataset
- HiRO 480 once added

CardioBoost V2 / DYNA Zenodo is different. It contains protein sequence pairs
and labels, but not genomic coordinates or gene symbols. Keep it as an auxiliary
sequence-level dataset; do not use it as a replacement for the original
CardioBoost public dataset.

Feature-injection sources such as gnomAD, AlphaMissense, CADD, ClinGen,
AlphaFold, and UniProt are tracked under `datasets/feature_sources/`.

Each feature-source folder may contain:

- `raw/`: downloaded raw files or source metadata.
- `external/`: small inputs such as gene panels or query files.
- `interim/`: generated example outputs and API responses.
- `scripts/`: source-specific acquisition or annotation scripts.
