#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
cd "$ROOT"

TABLE="datasets/external_validation_candidates/expanded_external_2026/modeling_table_expanded_external_annotated.tsv"
AUDIT="datasets/feature_sources/protein_structure/interim/maximized_foldx_mapping_audit.tsv"
STRUCTURE="datasets/feature_sources/protein_structure/interim/maximized_foldx_structure_features.tsv"
DDG="datasets/feature_sources/protein_structure/interim/maximized_foldx51_ddg.tsv"
SUMMARY="datasets/feature_sources/protein_structure/interim/maximized_foldx51_ddg.summary.json"
WORK="datasets/feature_sources/protein_structure/interim/registry_foldx51_repaired"
PDB="datasets/feature_sources/protein_structure/interim/registry_structure/pdb"
AUGMENT_DIR="datasets/feature_sources/protein_structure/interim/maximized_foldx_augmentations"
MODEL_ROOT="results/models/maximized_foldx_2026"
PREDICTION_FILE="primary_binary_predictions_split_source_heldout.tsv"

python3 datasets/feature_sources/protein_structure/scripts/build_maximized_foldx_mapping.py

python3 datasets/feature_sources/protein_structure/scripts/run_registry_foldx51_repaired_ddg.py \
  --table "$TABLE" \
  --audit "$AUDIT" \
  --structure "$STRUCTURE" \
  --pdb-dir "$PDB" \
  --work-dir "$WORK" \
  --output "$DDG" \
  --summary "$SUMMARY" \
  --min-plddt 0 \
  --repair-workers 12 \
  --build-workers 12 \
  --repair-timeout-seconds 14400 \
  --build-timeout-seconds 7200

python3 datasets/feature_sources/protein_structure/scripts/build_maximized_foldx_augmentations.py \
  --audit "$AUDIT" \
  --ddg "$DDG" \
  --out-dir "$AUGMENT_DIR"

python3 datasets/modeling/scripts/train_primary_binary_catboost.py \
  --table "$TABLE" \
  --split-col split_source_heldout \
  --slice-col training_slice_primary_binary \
  --exclude-prefix hiro_ \
  --exclude-prefix foldx_ddg_ \
  --out-dir "$MODEL_ROOT/baseline"

for tier in ge70 ge50 any; do
  if [[ "$tier" == "any" ]]; then
    eligibility_column="foldx_eligible_any_plddt"
  else
    eligibility_column="foldx_eligible_plddt_${tier}"
  fi
  python3 datasets/modeling/scripts/train_primary_binary_catboost.py \
    --table "$TABLE" \
    --split-col split_source_heldout \
    --slice-col training_slice_primary_binary \
    --exclude-prefix hiro_ \
    --exclude-prefix foldx_ddg_ \
    --augment-table "$AUGMENT_DIR/eligibility_only_${tier}.tsv" \
    --out-dir "$MODEL_ROOT/eligibility_${tier}"

  python3 datasets/modeling/scripts/train_primary_binary_catboost.py \
    --table "$TABLE" \
    --split-col split_source_heldout \
    --slice-col training_slice_primary_binary \
    --exclude-prefix hiro_ \
    --exclude-prefix foldx_ddg_ \
    --augment-table "$AUGMENT_DIR/ddg_features_${tier}.tsv" \
    --out-dir "$MODEL_ROOT/ddg_${tier}"

  python3 results/model_performance/maximized_foldx_enhancement/evaluate_maximized_foldx_models.py \
    --baseline-predictions "$MODEL_ROOT/baseline/$PREDICTION_FILE" \
    --eligibility-predictions "$MODEL_ROOT/eligibility_${tier}/$PREDICTION_FILE" \
    --ddg-predictions "$MODEL_ROOT/ddg_${tier}/$PREDICTION_FILE" \
    --eligibility-column "$eligibility_column" \
    --out-dir "results/model_performance/maximized_foldx_enhancement/${tier}" \
    --bootstrap-replicates 2000
done
