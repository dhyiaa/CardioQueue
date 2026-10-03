#!/usr/bin/env python3
"""Download a large HTTP file with resumable byte-range chunks."""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path


def head(url: str, timeout: int) -> tuple[int, bool]:
    request = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(request, timeout=timeout) as response:
        length = int(response.headers["Content-Length"])
        accept_ranges = response.headers.get("Accept-Ranges", "").lower() == "bytes"
    return length, accept_ranges


def download_chunk(
    url: str,
    part_path: Path,
    start: int,
    end: int,
    timeout: int,
    retries: int,
    retry_delay: int,
) -> tuple[Path, int]:
    expected = end - start + 1
    if part_path.exists() and part_path.stat().st_size == expected:
        return part_path, expected

    tmp_path = part_path.with_suffix(part_path.suffix + ".tmp")
    last_error: Exception | None = None
    for attempt in range(retries + 1):
        try:
            request = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
            with urllib.request.urlopen(request, timeout=timeout) as response, tmp_path.open("wb") as handle:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)

            actual = tmp_path.stat().st_size
            if actual != expected:
                raise RuntimeError(f"{part_path.name}: expected {expected} bytes, got {actual}")
            tmp_path.replace(part_path)
            return part_path, actual
        except (OSError, urllib.error.URLError, RuntimeError) as exc:
            last_error = exc
            try:
                tmp_path.unlink()
            except FileNotFoundError:
                pass
            if attempt < retries:
                time.sleep(retry_delay)
                continue
    raise RuntimeError(f"{part_path.name}: failed after {retries + 1} attempts: {last_error}")


def concatenate(parts: list[Path], output: Path, expected_size: int) -> None:
    tmp_output = output.with_suffix(output.suffix + ".tmp")
    with tmp_output.open("wb") as out:
        for part in parts:
            with part.open("rb") as handle:
                while True:
                    chunk = handle.read(1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
    actual = tmp_output.stat().st_size
    if actual != expected_size:
        raise RuntimeError(f"final size mismatch: expected {expected_size}, got {actual}")
    tmp_output.replace(output)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--chunk-mb", type=int, default=64)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--retry-delay", type=int, default=10)
    parser.add_argument("--keep-parts", action="store_true")
    args = parser.parse_args()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    expected_size, accept_ranges = head(args.url, timeout=args.timeout)
    if not accept_ranges:
        raise RuntimeError("server does not advertise byte-range support")

    if output.exists() and output.stat().st_size == expected_size:
        print(f"already complete: {output} bytes={expected_size}")
        return 0
    if output.exists() and output.stat().st_size != expected_size:
        backup = output.with_suffix(output.suffix + ".single_curl_partial")
        output.replace(backup)
        print(f"moved incomplete existing file to {backup}", file=sys.stderr)

    part_dir = output.with_suffix(output.suffix + ".parts")
    part_dir.mkdir(parents=True, exist_ok=True)
    chunk_size = args.chunk_mb * 1024 * 1024
    ranges = []
    for start in range(0, expected_size, chunk_size):
        end = min(start + chunk_size - 1, expected_size - 1)
        index = len(ranges)
        ranges.append((part_dir / f"part_{index:05d}", start, end))

    print(
        f"downloading {expected_size} bytes in {len(ranges)} parts "
        f"with {args.workers} workers"
    )
    completed = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_to_part = {
            executor.submit(
                download_chunk,
                args.url,
                part,
                start,
                end,
                args.timeout,
                args.retries,
                args.retry_delay,
            ): part
            for part, start, end in ranges
        }
        for future in concurrent.futures.as_completed(future_to_part):
            part, size = future.result()
            completed += size
            if completed // (512 * 1024 * 1024) != (completed - size) // (512 * 1024 * 1024):
                print(f"downloaded at least {completed / (1024 ** 3):.2f} GiB")

    parts = [part for part, _, _ in ranges]
    concatenate(parts, output, expected_size)
    print(f"complete: {output} bytes={os.path.getsize(output)}")

    if not args.keep_parts:
        for part in parts:
            part.unlink()
        part_dir.rmdir()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
