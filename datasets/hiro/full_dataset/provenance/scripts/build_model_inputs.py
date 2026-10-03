#!/usr/bin/env python3
"""Build frozen ML and LLM input artifacts from the annotated cohort."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cardio_pipeline.features import build_model_inputs, read_csv_rows, write_model_input_outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="outputs/annotations/annotated_cohort.csv")
    parser.add_argument("--output-dir", default="outputs/model_inputs")
    args = parser.parse_args()

    rows, _ = read_csv_rows(args.input)
    result = build_model_inputs(rows)
    paths = write_model_input_outputs(result, args.output_dir)

    print(json.dumps(result["report"], indent=2, sort_keys=True))
    print("\nWrote:")
    for name, path in paths.items():
        print(f"  {name}: {path}")


if __name__ == "__main__":
    main()
