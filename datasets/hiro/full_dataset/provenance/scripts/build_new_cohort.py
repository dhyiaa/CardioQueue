#!/usr/bin/env python3
"""Build normalized CSVs from the April 2026 CASPER WES and VERDICT files."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cardio_pipeline.new_data import build_new_data_cohort, write_new_data_outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default="new_data", help="Directory containing Brianna/Doug CSV files.")
    parser.add_argument("--output-dir", default="outputs/new_data", help="Directory for normalized outputs.")
    args = parser.parse_args()

    result = build_new_data_cohort(args.data_dir)
    paths = write_new_data_outputs(result, args.output_dir)

    print(json.dumps(result["report"], indent=2))
    print("\nWrote:")
    for name, path in paths.items():
        print(f"  {name}: {path}")


if __name__ == "__main__":
    main()
