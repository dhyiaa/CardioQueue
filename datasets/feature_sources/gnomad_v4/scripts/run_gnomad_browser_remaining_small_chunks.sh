#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../../../.."

PROJECT="cardiogenetics-gnomad"
CHUNK_SIZE="${CHUNK_SIZE:-100}"
INPUT="datasets/feature_sources/gnomad_v4/external/gnomad_missing_from_dbnsfp.clean_combined.tsv"
COMPLETED_OUT_DIR="datasets/feature_sources/gnomad_v4/interim/gnomad_browser_v4_1_missing_chunks"
REMAINING_INPUT="datasets/feature_sources/gnomad_v4/external/gnomad_missing_from_dbnsfp.clean_combined.remaining_after_completed.tsv"
CHUNK_DIR="datasets/feature_sources/gnomad_v4/external/gnomad_missing_from_dbnsfp.clean_combined.remaining_small_chunks"
OUT_DIR="datasets/feature_sources/gnomad_v4/interim/gnomad_browser_v4_1_missing_remaining_small_chunks"
COMBINED="datasets/feature_sources/gnomad_v4/interim/gnomad_browser_v4_1_missing_combined_from_mixed_chunks.tsv.gz"
LOG="datasets/feature_sources/gnomad_v4/interim/gnomad_browser_v4_1_missing_remaining_small_chunks.run.log"
LOCK_DIR="datasets/feature_sources/gnomad_v4/interim/gnomad_browser_v4_1_missing_remaining_small_chunks.lock"

mkdir -p "$CHUNK_DIR" "$OUT_DIR" "$(dirname "$COMBINED")"

if ! mkdir "$LOCK_DIR" 2>/dev/null; then
  if ! pgrep -f "extract_gnomad_browser_hail_missing.py.*remaining_small_chunks" >/dev/null 2>&1; then
    echo "Removing stale gnomAD small-chunk lock: $LOCK_DIR" >&2
    rm -rf "$LOCK_DIR"
    mkdir "$LOCK_DIR"
  else
    echo "Another gnomAD small-chunk run appears to be active: $LOCK_DIR" >&2
    exit 1
  fi
fi
trap 'rm -rf "$LOCK_DIR"' EXIT

source "$HOME/opt/miniforge3/etc/profile.d/conda.sh"
conda activate hail-gnomad
export JAVA_HOME="$HOME/opt/jdk-11.0.31+11/Contents/Home"
export PATH="$JAVA_HOME/bin:$PATH"
export INPUT COMPLETED_OUT_DIR REMAINING_INPUT CHUNK_DIR CHUNK_SIZE

build_remaining_and_chunks() {
  rm -f "$REMAINING_INPUT" "$CHUNK_DIR"/chunk_*.tsv "$CHUNK_DIR/.complete"
  python - <<'PY'
import csv
import gzip
import os
from pathlib import Path

input_path = Path(os.environ["INPUT"])
completed_out_dir = Path(os.environ["COMPLETED_OUT_DIR"])
remaining_input = Path(os.environ["REMAINING_INPUT"])
chunk_dir = Path(os.environ["CHUNK_DIR"])
chunk_size = int(os.environ["CHUNK_SIZE"])

completed_ids = set()
for output in sorted(completed_out_dir.glob("chunk_*.tsv.bgz")):
    try:
        with gzip.open(output, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            ids = [row.get("variant_id", "") for row in reader]
        if ids:
            completed_ids.update(v for v in ids if v)
    except Exception:
        pass

with input_path.open(newline="") as in_handle, remaining_input.open("w", newline="") as out_handle:
    reader = csv.DictReader(in_handle, delimiter="\t")
    writer = csv.DictWriter(out_handle, fieldnames=reader.fieldnames, delimiter="\t")
    writer.writeheader()
    remaining_rows = 0
    skipped_rows = 0
    for row in reader:
        if row.get("variant_id") in completed_ids:
            skipped_rows += 1
            continue
        writer.writerow(row)
        remaining_rows += 1

with remaining_input.open(newline="") as handle:
    header = handle.readline()
    chunk_index = 0
    rows_in_chunk = 0
    out = None
    for line in handle:
        if out is None or rows_in_chunk >= chunk_size:
            if out is not None:
                out.close()
            out = (chunk_dir / f"chunk_{chunk_index:05d}.tsv").open("w")
            out.write(header)
            chunk_index += 1
            rows_in_chunk = 0
        out.write(line)
        rows_in_chunk += 1
    if out is not None:
        out.close()

(chunk_dir / ".complete").write_text(
    f"chunks={chunk_index}\nchunk_size={chunk_size}\n"
    f"completed_ids={len(completed_ids)}\nskipped_rows={skipped_rows}\nremaining_rows={remaining_rows}\n"
)
print((chunk_dir / ".complete").read_text(), end="")
PY
}

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
  echo "[$(date)] Starting remaining-only small-chunk gnomAD browser Hail run"
  echo "project=$PROJECT"
  echo "input=$INPUT"
  echo "completed_out_dir=$COMPLETED_OUT_DIR"
  echo "remaining_input=$REMAINING_INPUT"
  echo "chunk_dir=$CHUNK_DIR"
  echo "out_dir=$OUT_DIR"
  echo "combined=$COMBINED"
  echo "chunk_size=$CHUNK_SIZE"
  echo

  build_remaining_and_chunks
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
      --gcp-project "$PROJECT" \
      --tmp-dir "tmp/hail/gnomad_browser_missing_${chunk_base}"
    echo "[$(date)] DONE $chunk_base"
  done

  tmp_combined="${COMBINED}.tmp"
  rm -f "$tmp_combined"
  first=1
  for output in "$COMPLETED_OUT_DIR"/chunk_*.tsv.bgz "$OUT_DIR"/chunk_*.tsv.bgz; do
    [ -s "$output" ] || continue
    if [ "$first" = 1 ]; then
      gzip -dc "$output"
      first=0
    else
      gzip -dc "$output" | tail -n +2
    fi
  done | gzip > "$tmp_combined"
  mv "$tmp_combined" "$COMBINED"

  echo
  echo "[$(date)] Finished remaining-only small-chunk gnomAD browser Hail run"
  echo "combined_rows=$(gzip -dc "$COMBINED" | wc -l | tr -d ' ')"
} 2>&1 | tee -a "$LOG"
