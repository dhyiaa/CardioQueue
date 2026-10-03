#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../../../.."

PROJECT="cardiogenetics-gnomad"
CHUNK_SIZE="${CHUNK_SIZE:-1000}"
INPUT="datasets/feature_sources/gnomad_v4/external/gnomad_missing_from_dbnsfp.clean_combined.tsv"
CHUNK_DIR="datasets/feature_sources/gnomad_v4/external/gnomad_missing_from_dbnsfp.clean_combined.chunks"
OUT_DIR="datasets/feature_sources/gnomad_v4/interim/gnomad_browser_v4_1_missing_chunks"
COMBINED="datasets/feature_sources/gnomad_v4/interim/gnomad_browser_v4_1_missing_clean_combined.tsv.gz"
LOG="datasets/feature_sources/gnomad_v4/interim/gnomad_browser_v4_1_missing_chunks.run.log"
LOCK_DIR="datasets/feature_sources/gnomad_v4/interim/gnomad_browser_v4_1_missing_chunks.lock"

mkdir -p "$CHUNK_DIR" "$OUT_DIR" "$(dirname "$COMBINED")"

if ! mkdir "$LOCK_DIR" 2>/dev/null; then
  if ! pgrep -f "extract_gnomad_browser_hail_missing.py.*gnomad_missing_from_dbnsfp.clean_combined.chunks" >/dev/null 2>&1; then
    echo "Removing stale gnomAD chunk lock: $LOCK_DIR" >&2
    rm -rf "$LOCK_DIR"
    mkdir "$LOCK_DIR"
  else
  echo "Another gnomAD chunk run appears to be active: $LOCK_DIR" >&2
  echo "If no run is active, remove the stale lock and rerun:" >&2
  echo "  rm -rf $LOCK_DIR" >&2
  exit 1
  fi
fi
trap 'rm -rf "$LOCK_DIR"' EXIT

source "$HOME/opt/miniforge3/etc/profile.d/conda.sh"
conda activate hail-gnomad
export JAVA_HOME="$HOME/opt/jdk-11.0.31+11/Contents/Home"
export PATH="$JAVA_HOME/bin:$PATH"
export INPUT CHUNK_DIR CHUNK_SIZE

if [ ! -f "$CHUNK_DIR/.complete" ]; then
  rm -f "$CHUNK_DIR"/chunk_*.tsv "$CHUNK_DIR/.complete"
  python - <<'PY'
import os
from pathlib import Path

input_path = Path(os.environ["INPUT"])
chunk_dir = Path(os.environ["CHUNK_DIR"])
chunk_size = int(os.environ["CHUNK_SIZE"])

with input_path.open() as handle:
    header = handle.readline()
    chunk_index = 0
    rows_in_chunk = 0
    out = None
    for line in handle:
        if out is None or rows_in_chunk >= chunk_size:
            if out is not None:
                out.close()
            out = (chunk_dir / f"chunk_{chunk_index:04d}.tsv").open("w")
            out.write(header)
            chunk_index += 1
            rows_in_chunk = 0
        out.write(line)
        rows_in_chunk += 1
    if out is not None:
        out.close()

(chunk_dir / ".complete").write_text(f"chunks={chunk_index}\nchunk_size={chunk_size}\n")
PY
fi

is_complete_output() {
  local input="$1"
  local output="$2"
  [ -s "$output" ] || return 1
  gzip -t "$output" >/dev/null 2>&1 || return 1
  local input_lines output_lines
  input_lines=$(wc -l < "$input" | tr -d ' ')
  output_lines=$(gzip -dc "$output" | wc -l | tr -d ' ')
  [ "$input_lines" = "$output_lines" ]
}

{
  echo "[$(date)] Starting resumable gnomAD browser Hail chunk run"
  echo "project=$PROJECT"
  echo "input=$INPUT"
  echo "chunk_dir=$CHUNK_DIR"
  echo "out_dir=$OUT_DIR"
  echo "combined=$COMBINED"
  echo "chunk_size=$CHUNK_SIZE"
  echo

  for chunk in "$CHUNK_DIR"/chunk_*.tsv; do
    chunk_base="$(basename "$chunk" .tsv)"
    output="$OUT_DIR/${chunk_base}.tsv.bgz"
    chunk_rows=$(( $(wc -l < "$chunk" | tr -d ' ') - 1 ))

    if is_complete_output "$chunk" "$output"; then
      echo "[$(date)] SKIP complete $chunk_base rows=$chunk_rows"
      continue
    fi

    rm -f "$output"
    echo "[$(date)] RUN $chunk_base rows=$chunk_rows"
    python datasets/feature_sources/gnomad_v4/scripts/extract_gnomad_browser_hail_missing.py \
      --input "$chunk" \
      --output "$output" \
      --gcp-project "$PROJECT"
    echo "[$(date)] DONE $chunk_base"
  done

  tmp_combined="${COMBINED}.tmp"
  rm -f "$tmp_combined"
  first=1
  for output in "$OUT_DIR"/chunk_*.tsv.bgz; do
    if [ "$first" = 1 ]; then
      gzip -dc "$output"
      first=0
    else
      gzip -dc "$output" | tail -n +2
    fi
  done | gzip > "$tmp_combined"
  mv "$tmp_combined" "$COMBINED"

  echo
  echo "[$(date)] Finished resumable gnomAD browser Hail chunk run"
  echo "combined_rows=$(gzip -dc "$COMBINED" | wc -l | tr -d ' ')"
} 2>&1 | tee -a "$LOG"
