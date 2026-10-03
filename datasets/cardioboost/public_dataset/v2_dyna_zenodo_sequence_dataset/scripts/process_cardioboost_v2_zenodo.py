#!/usr/bin/env python3
"""Process the Zenodo DYNA/CardioBoost V2 sequence datasets."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


PROTEIN_FILES = {
    "cm_train_data_1024.csv": ("cardioboost_original", "cardiomyopathy", "train", "binary"),
    "cm_test_data_1024.csv": ("cardioboost_original", "cardiomyopathy", "test", "binary"),
    "arm_train_data_1024.csv": ("cardioboost_original", "arrhythmia", "train", "binary"),
    "arm_test_data_1024.csv": ("cardioboost_original", "arrhythmia", "test", "binary"),
    "Clinvar_CM_protein.csv": ("clinvar_protein_2024", "cardiomyopathy", "unsplit", "binary"),
    "Clinvar_ARM_protein.csv": ("clinvar_protein_2024", "arrhythmia", "unsplit", "binary"),
    "Clinvar_CM_VUS_protein.csv": ("clinvar_protein_2024", "cardiomyopathy", "unsplit", "vus"),
    "Clinvar_ARM_VUS_protein.csv": ("clinvar_protein_2024", "arrhythmia", "unsplit", "vus"),
}


SPLICING_FILES = {
    "mfass_train.csv": ("mfass", "train"),
    "mfass_val.csv": ("mfass", "validation"),
    "mfass_test.csv": ("mfass", "test"),
    "clinvar_test.csv": ("clinvar_splicing", "test"),
}


def clinical_label(label: str, label_type: str) -> str:
    if label_type == "vus":
        return "VUS"
    if label == "1":
        return "Pathogenic"
    if label == "0":
        return "Benign"
    return "Unknown"


def diff_features(wt: str, mut: str) -> dict[str, object]:
    length = min(len(wt), len(mut))
    diffs = [i for i in range(length) if wt[i] != mut[i]]
    insertion_or_deletion = len(wt) != len(mut)
    first = diffs[0] if diffs else None
    return {
        "wt_len": len(wt),
        "mut_len": len(mut),
        "diff_count": len(diffs) + abs(len(wt) - len(mut)),
        "has_indel_or_length_change": insertion_or_deletion,
        "first_diff_local_pos_1based": "" if first is None else first + 1,
        "wt_aa": "" if first is None else wt[first],
        "mut_aa": "" if first is None else mut[first],
    }


def process_protein(raw_dir: Path, output: Path) -> dict:
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "source_file",
        "source_family",
        "disease_panel",
        "split",
        "clinical_label",
        "original_label",
        "start_pos",
        "wt_len",
        "mut_len",
        "diff_count",
        "has_indel_or_length_change",
        "first_diff_local_pos_1based",
        "wt_aa",
        "mut_aa",
        "wt_seq",
        "mut_seq",
    ]
    summary = {}
    with output.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields)
        writer.writeheader()
        for filename, (source_family, disease_panel, split, label_type) in PROTEIN_FILES.items():
            path = raw_dir / filename
            counts = {"rows": 0, "labels": {}, "diff_counts": {}}
            with path.open(newline="", errors="replace") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    label = clinical_label(row["labels"], label_type)
                    features = diff_features(row["wt_seq"], row["mut_seq"])
                    out_row = {
                        "source_file": filename,
                        "source_family": source_family,
                        "disease_panel": disease_panel,
                        "split": split,
                        "clinical_label": label,
                        "original_label": row["labels"],
                        "start_pos": row.get("start_pos", ""),
                        "wt_seq": row["wt_seq"],
                        "mut_seq": row["mut_seq"],
                        **features,
                    }
                    writer.writerow(out_row)
                    counts["rows"] += 1
                    counts["labels"][label] = counts["labels"].get(label, 0) + 1
                    diff_key = str(features["diff_count"])
                    counts["diff_counts"][diff_key] = counts["diff_counts"].get(diff_key, 0) + 1
            summary[filename] = counts
    return summary


def process_splicing_summary(raw_dir: Path, output: Path) -> dict:
    rows = []
    for filename, (source_family, split) in SPLICING_FILES.items():
        path = raw_dir / filename
        counts = {"rows": 0, "labels": {}}
        with path.open(newline="", errors="replace") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                counts["rows"] += 1
                counts["labels"][row["labels"]] = counts["labels"].get(row["labels"], 0) + 1
        rows.append({"source_file": filename, "source_family": source_family, "split": split, **counts})
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=["source_file", "source_family", "split", "rows", "labels"])
        writer.writeheader()
        for row in rows:
            row = {**row, "labels": json.dumps(row["labels"], sort_keys=True)}
            writer.writerow(row)
    return {row["source_file"]: {"rows": row["rows"], "labels": row["labels"]} for row in rows}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", default=Path("datasets/cardioboost/public_dataset/v2_dyna_zenodo_sequence_dataset"), type=Path)
    args = parser.parse_args()

    raw_dir = args.base / "raw_zenodo"
    processed = args.base / "processed"
    provenance = args.base / "provenance"
    protein_output = processed / "cardioboost_v2_protein_variant_sequences.csv"
    splicing_output = processed / "dyna_splicing_sequence_summary.csv"
    summary = {
        "protein": process_protein(raw_dir, protein_output),
        "splicing": process_splicing_summary(raw_dir, splicing_output),
        "outputs": {
            "protein_variant_sequences": str(protein_output),
            "splicing_sequence_summary": str(splicing_output),
        },
        "label_mapping": {
            "binary_files": {"0": "Benign", "1": "Pathogenic"},
            "vus_files": {"clinical_label": "VUS", "original_label": "preserved"},
        },
        "important_limitation": "Zenodo protein files contain sequence pairs, not genomic coordinates or gene symbols.",
    }
    json.dump(summary, open(provenance / "processed_summary.json", "w"), indent=2)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

