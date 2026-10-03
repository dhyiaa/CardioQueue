#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../../../.."

INPUT="datasets/feature_sources/gnomad_v4/external/gnomad_missing_from_dbnsfp.clean_combined.tsv"
OUTPUT="datasets/feature_sources/gnomad_v4/interim/gnomad_browser_v4_1_missing_clean_combined.tsv.bgz"
LOG="datasets/feature_sources/gnomad_v4/interim/gnomad_browser_v4_1_missing_clean_combined.run.log"
PROJECT="cardiogenetics-gnomad"

mkdir -p "$(dirname "$OUTPUT")"

{
  echo "[$(date)] Starting full gnomAD browser Hail extraction"
  echo "input=$INPUT"
  echo "output=$OUTPUT"
  echo "project=$PROJECT"
  echo

  source "$HOME/opt/miniforge3/etc/profile.d/conda.sh"
  conda activate hail-gnomad
  export JAVA_HOME="$HOME/opt/jdk-11.0.31+11/Contents/Home"
  export PATH="$JAVA_HOME/bin:$PATH"

  python datasets/feature_sources/gnomad_v4/scripts/extract_gnomad_browser_hail_missing.py \
    --input "$INPUT" \
    --output "$OUTPUT" \
    --gcp-project "$PROJECT"

  echo
  echo "[$(date)] Finished full gnomAD browser Hail extraction"
} 2>&1 | tee "$LOG"
