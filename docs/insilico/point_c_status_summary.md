# Point C: In Silico Annotation Status

In silico annotations are feature-engineering sources only. They are not variant
sources and should not be used as clinical labels.

This document tracks Point C from the project proposal: dbNSFP, Ensembl VEP,
CADD, AlphaMissense, and ESM1b/ESM-variant. It also tracks parity with the
previous study's in silico panel so the current project does not regress.

## Status Definitions

| Status | Meaning |
|---|---|
| Downloaded and validated | The official file is present locally and basic file integrity/schema checks have passed. |
| API-ready | A tested script can query the source for small or targeted variant sets and records explicit status/missing-reason fields. |
| Deferred | The source is valid, but bulk download or local installation is intentionally delayed because of size, setup, or strategy. |
| Pending official access | The source requires registration, institutional access, or a non-public access link before download. |
| Pending source confirmation | The source is conceptually required, but the exact public file and mapping strategy still need to be fixed. |

## Summary Table

| Source / tool | Main features expected | Current status | Local/API path | What is done | What remains |
|---|---|---:|---|---|---|
| AlphaMissense hg38 | `am_pathogenicity`, `am_class`, UniProt ID, transcript ID, protein variant | Downloaded and validated | Local file: `datasets/feature_sources/insilico/raw/alphamissense/AlphaMissense_hg38.tsv.gz`; join script: `datasets/feature_sources/insilico/scripts/join_alphamissense.py` | Official hg38 file downloaded, gzip checked, MD5 recorded, line count checked, example join validated | Scale join to the full final ClinVar+gnomAD variant list |
| CADD v1.7 GRCh38 | CADD raw score, CADD PHRED | API-ready for small sets; bulk deferred | API script: `datasets/feature_sources/insilico/scripts/cadd_api.py`; CADD MD5 file saved in `datasets/feature_sources/insilico/raw/cadd/` | CADD doctor query passed; example ClinVar SNV annotation validated; status/missing-reason columns implemented | For large-scale use, choose either targeted API for small panels, tabix remote lookup, or selected bulk download. Do not download 81G whole-genome SNV file blindly |
| dbNSFP v5.3.1 academic branch | Old-SI parity: REVEL, SIFT, PolyPhen-2 HDIV, MetaLR, FATHMM-XF. New proposal additions: MutationTaster, MutationAssessor, FATHMM, PROVEAN, VEST, M-CAP, MetaSVM, Eigen, GERP++, PhyloP, PhastCons, SiPhy, and newer dbNSFP predictors | Downloaded, validated, and joined | GRCh38 BGZF: `datasets/feature_sources/insilico/raw/dbnsfp/downloads/dbNSFP5.3.1a_grch38.gz`; feature table: `datasets/feature_sources/insilico/interim/dbnsfp_features_master_variant_registry.tsv` | Official academic GRCh38 file downloaded cleanly, MD5 matched `afff632f13325c687477d64ca0c63b1a`, gene table passed gzip test, header/schema smoke test passed, batched `pysam` join completed for the master registry | Merge selected `dbnsfp_*` columns into the final in silico/modeling feature table |
| Ensembl VEP 116 GRCh38 | Consequence annotation, transcript effects, optional plugin framework | Deferred large local cache | Official cache checked: `homo_sapiens_vep_116_GRCh38.tar.gz`; documented in `datasets/feature_sources/insilico/raw/vep/README.md` | Current human GRCh38 cache size checked: `27,644,657,162` bytes compressed | Decide whether local VEP is needed. If yes, download cache and install VEP/plugins. If no, rely on gnomAD consequences plus dbNSFP/AlphaMissense/CADD |
| SpliceAI through VEP plugin | SpliceAI delta scores and splice-impact annotations | Pending VEP/plugin decision | Not installed yet | Requirement documented | Need local VEP/plugin setup or another verified precomputed SpliceAI source |
| ESM1b / ESM-variant | Protein language model variant-effect scores | Pending mapping/source confirmation | Documented in `datasets/feature_sources/insilico/raw/esm1b/README.md` | Recognized as new-proposal addition | Confirm exact public precomputed human proteome source; define transcript/protein mapping QC; then download and join |

## Current Local Deliverables

| Deliverable | Path | Purpose |
|---|---|---|
| AlphaMissense raw file | `datasets/feature_sources/insilico/raw/alphamissense/AlphaMissense_hg38.tsv.gz` | Local hg38 precomputed missense pathogenicity scores |
| AlphaMissense join script | `datasets/feature_sources/insilico/scripts/join_alphamissense.py` | Exact `chrom-pos-ref-alt` join with explicit status fields |
| CADD API script | `datasets/feature_sources/insilico/scripts/cadd_api.py` | Targeted CADD raw/PHRED annotation for small QC/test sets |
| dbNSFP join script | `datasets/feature_sources/insilico/scripts/join_dbnsfp.py` | Batched `pysam`/tabix-indexed dbNSFP join for registry variants |
| dbNSFP registry feature table | `datasets/feature_sources/insilico/interim/dbnsfp_features_master_variant_registry.tsv` | Full master-registry dbNSFP feature layer; 161,108 rows |
| dbNSFP registry feature summary | `datasets/feature_sources/insilico/interim/dbnsfp_features_master_variant_registry.summary.json` | Match-rate and selected feature coverage summary |
| Merged in silico example | `datasets/feature_sources/insilico/interim/insilico_features_clinvar_example.tsv` | Current example feature table combining implemented sources |
| Requirement metadata | `metadata/insilico_feature_requirements.json` | Machine-readable feature coverage requirements |
| Requirement validator | `datasets/feature_sources/insilico/scripts/validate_insilico_requirements.py` | Detects missing old-SI required predictor groups |

## Old-SI Feature Parity

| Predictor from old SI | Current coverage | Current source | Status |
|---|---|---|---|
| REVEL | Present in dbNSFP registry feature table | dbNSFP preferred | Joined; 70,257 nonmissing registry rows |
| CADD PHRED | Present in example workflow | CADD API, `GRCh38-v1.7` | Implemented for targeted use |
| SIFT | Present in dbNSFP registry feature table | dbNSFP preferred | Joined; 69,739 nonmissing registry rows |
| PolyPhen-2 HDIV | Present in dbNSFP registry feature table | dbNSFP preferred | Joined; 66,702 nonmissing registry rows |
| MetaLR | Present in dbNSFP registry feature table | dbNSFP preferred | Joined; 70,419 nonmissing registry rows |
| FATHMM-XF | Present in dbNSFP registry feature table | dbNSFP preferred | Joined; 67,283 nonmissing registry rows |
| AlphaMissense | Present | AlphaMissense hg38 local file | Downloaded and validated |

Current old-SI parity status:

- Complete now: `AlphaMissense`, `CADD PHRED`
- Complete as a source-specific feature layer: `REVEL`, `SIFT`,
  `PolyPhen-2 HDIV`, `MetaLR`, `FATHMM-XF`
- Still required: merge/rename selected `dbnsfp_*` columns into the final
  modeling feature table and run strict validation there

## Null and Failure Handling Rules

Every Point C output must distinguish biological absence from pipeline failure.
Blank cells alone are not acceptable.

| Situation | Required handling |
|---|---|
| Exact annotation found | Set source status to `ok` and populate required feature fields. |
| Variant not represented in the source | Set source status to `not_found` or source-specific equivalent and explain in missing-reason field. |
| Variant is outside source scope | Set status to `not_applicable`; example: non-missense variant in AlphaMissense. |
| Coordinates are insufficient | Set status to `not_queryable` and preserve the reason. |
| API/network/schema failure | Fail the run or mark explicit error status; do not silently create null feature values. |
| Protein-level mapping ambiguity | Do not populate ESM1b/ESM-variant features until transcript/protein mapping is resolved. |

## Validation Commands

Current example in silico table:

```bash
python3 datasets/feature_sources/insilico/scripts/merge_insilico_features.py validate \
  --input datasets/feature_sources/insilico/interim/insilico_features_clinvar_example.tsv

python3 datasets/feature_sources/insilico/scripts/validate_insilico_requirements.py \
  --input datasets/feature_sources/insilico/interim/insilico_features_clinvar_example.tsv
```

Expected result before the dbNSFP feature layer is merged into the final table:

- The merged table validates.
- The requirement check reports five missing groups: `REVEL`, `SIFT`,
  `PolyPhen2_HDIV`, `MetaLR`, and `FATHMM_XF`.

After the dbNSFP columns are merged into the final feature table, rerun the
requirement check with `--strict`.

```bash
python3 datasets/feature_sources/insilico/scripts/validate_insilico_requirements.py \
  --input datasets/feature_sources/insilico/interim/insilico_features_clinvar_full.tsv \
  --strict
```

## Immediate Next Work

1. Merge selected `dbnsfp_*` columns into the final in silico/modeling feature
   table.
2. Rebuild the merged in silico feature table and run strict validation.
3. Decide whether VEP/SpliceAI and ESM1b should be added before first modeling
   or deferred to the next feature block.

## Modeling Readiness

The current Point C feature layer is much closer to modeling-ready: dbNSFP has
been downloaded, validated, queried, and materialized for the master registry.
The remaining step is to merge/standardize dbNSFP columns into the final
modeling feature table alongside AlphaMissense, CADD, protein features, and
source labels.

Minimum Point C readiness for modeling:

1. AlphaMissense and CADD columns validated.
2. dbNSFP columns added for `REVEL`, `SIFT`, `PolyPhen-2 HDIV`, `MetaLR`, and
   `FATHMM-XF`.
3. New-proposal dbNSFP predictors selected from the official v5.3.1 README.
4. VEP/SpliceAI decision documented.
5. ESM1b/ESM-variant source and mapping QC documented or explicitly deferred.
