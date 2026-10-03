#!/usr/bin/env python3
"""Inspect a downloaded AlphaFold tar without extracting all structures."""

from __future__ import annotations

import argparse
import json
import tarfile
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tar", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    tar_path = Path(args.tar)
    manifest_path = Path(args.manifest)
    summary_path = Path(args.summary)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    counts: dict[str, int] = {}
    total_members = 0
    total_member_bytes = 0
    first_members: list[str] = []

    with tarfile.open(tar_path, "r") as archive, manifest_path.open("w") as manifest:
        for member in archive:
            total_members += 1
            total_member_bytes += member.size
            name = member.name
            manifest.write(name + "\n")
            suffix = Path(name).suffix.lower() or "<no_suffix>"
            counts[suffix] = counts.get(suffix, 0) + 1
            if len(first_members) < 20:
                first_members.append(name)

    summary = {
        "tar_path": str(tar_path),
        "tar_bytes": tar_path.stat().st_size,
        "total_members": total_members,
        "total_member_bytes": total_member_bytes,
        "suffix_counts": dict(sorted(counts.items())),
        "first_members": first_members,
    }
    summary_path.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
