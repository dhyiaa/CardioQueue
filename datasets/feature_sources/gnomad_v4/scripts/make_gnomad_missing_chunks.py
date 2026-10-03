#!/usr/bin/env python3
"""Split the dbNSFP-missing gnomAD input into resumable chunks."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        default="datasets/feature_sources/gnomad_v4/external/gnomad_missing_from_dbnsfp.clean_combined.tsv",
        type=Path,
    )
    parser.add_argument(
        "--out-dir",
        default="datasets/feature_sources/gnomad_v4/external/chunks",
        type=Path,
    )
    parser.add_argument("--chunk-size", default=2500, type=int)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    with args.input.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        chunk_index = 0
        row_count = 0
        writer = None
        out_handle = None

        try:
            for row in reader:
                if row_count % args.chunk_size == 0:
                    if out_handle is not None:
                        out_handle.close()
                    chunk_index += 1
                    out_path = args.out_dir / f"gnomad_missing_chunk_{chunk_index:04d}.tsv"
                    out_handle = out_path.open("w", newline="")
                    writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(header)
                writer.writerow(row)
                row_count += 1
        finally:
            if out_handle is not None:
                out_handle.close()

    print(f"rows={row_count}")
    print(f"chunks={chunk_index}")
    print(f"chunk_size={args.chunk_size}")
    print(f"out_dir={args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
