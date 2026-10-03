# CardioBoost V2 / DYNA Zenodo Dataset

This folder is a fresh, reproducible download from Zenodo record `13397296`:

- Title: DYNA: Disease-Specific Language Model for Variant Pathogenicity
- DOI: `10.5281/zenodo.13397296`
- License: CC-BY-4.0
- API metadata snapshot: `provenance/zenodo_record_13397296.json`
- Download manifest and MD5 validation: `provenance/download_manifest.json`

This is separate from `public_dataset`, which contains
the coordinate-rich GitHub/CardioBoost-derived handoff.

## What V2 Contains

The Zenodo files are sequence-pair datasets:

- `wt_seq` / `mut_seq` protein sequence pairs for cardiomyopathy and arrhythmia.
- `seq_a` / `seq_b` splicing sequence pairs for DYNA splicing tasks.
- `labels` from the Zenodo files.

Important limitation: the protein files do not include chromosome, position,
reference allele, alternate allele, or gene symbol. They are useful for
protein-sequence modeling and pretraining-style experiments, but they are not
drop-in replacements for coordinate-based variant tables.

## Processed Protein Table

Main processed file:

- `processed/cardioboost_v2_protein_variant_sequences.csv`

Rows: 21,472

| Clinical label | Count |
|---|---:|
| VUS | 18,619 |
| Benign | 1,702 |
| Pathogenic | 1,151 |

| Source family | Count |
|---|---:|
| `clinvar_protein_2024` | 20,352 |
| `cardioboost_original` | 1,120 |

| Disease panel | Count |
|---|---:|
| Cardiomyopathy | 20,027 |
| Arrhythmia | 1,445 |

The original CardioBoost sequence rows are:

| File group | Rows | Labels |
|---|---:|---|
| CM train | 440 | 238 Pathogenic, 202 Benign |
| CM test | 218 | 118 Pathogenic, 100 Benign |
| ARM train | 308 | 168 Pathogenic, 140 Benign |
| ARM test | 154 | 84 Pathogenic, 70 Benign |

## Label Mapping

For binary protein files:

- `labels = 0` -> Benign
- `labels = 1` -> Pathogenic

For VUS protein files:

- `clinical_label = VUS`
- the original `labels` column is preserved as `original_label`

This avoids incorrectly treating VUS-file `labels` as pathogenic labels.

## Splicing Files

The non-cardiogenetics splicing/DYNA files are retained raw and summarized in:

- `processed/dyna_splicing_sequence_summary.csv`

They are not merged into the cardiogenetics protein-variant table.

## Rebuild

```bash
python3 datasets/cardioboost/public_dataset/v2_dyna_zenodo_sequence_dataset/scripts/process_cardioboost_v2_zenodo.py \
  --base datasets/cardioboost/public_dataset/v2_dyna_zenodo_sequence_dataset
```

