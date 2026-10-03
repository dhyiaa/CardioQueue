#!/usr/bin/env python3
"""Create a VCF for VEP from the current modeling matrix."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd


def open_text(path: Path, mode: str = "wt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return path.open(mode)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--out-vcf", required=True, type=Path)
    args = parser.parse_args()

    cols = ["variant_id", "chrom", "pos", "ref", "alt"]
    df = pd.read_csv(args.matrix, sep="\t", dtype=str, usecols=cols)
    df = df.drop_duplicates("variant_id").copy()
    df["chrom"] = df["chrom"].str.removeprefix("chr")
    df["pos_int"] = pd.to_numeric(df["pos"], errors="coerce")
    bad = df[df[["variant_id", "chrom", "pos", "ref", "alt"]].isna().any(axis=1) | df["pos_int"].isna()]
    if len(bad):
        raise SystemExit(f"Cannot build VCF: {len(bad)} rows have incomplete coordinates")

    df["pos_int"] = df["pos_int"].astype(int)
    sort_chrom = df["chrom"].replace({"X": "23", "Y": "24", "MT": "25", "M": "25"})
    df["_chrom_sort"] = pd.to_numeric(sort_chrom, errors="coerce").fillna(99).astype(int)
    df = df.sort_values(["_chrom_sort", "pos_int", "ref", "alt"], kind="mergesort")

    args.out_vcf.parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.out_vcf, "wt") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##reference=GRCh38\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in df.itertuples(index=False):
            out.write(f"{row.chrom}\t{row.pos_int}\t{row.variant_id}\t{row.ref}\t{row.alt}\t.\t.\t.\n")

    print(f"Wrote {len(df):,} variants to {args.out_vcf}")


if __name__ == "__main__":
    main()
