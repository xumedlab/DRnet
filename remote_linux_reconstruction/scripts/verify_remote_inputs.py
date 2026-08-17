#!/usr/bin/env python3
"""Verify every uploaded FASTQ, reference, and tool archive in one read pass."""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path


FIELDS = [
    "role",
    "path",
    "status",
    "expected_bytes",
    "actual_bytes",
    "expected_md5",
    "actual_md5",
    "expected_sha256",
    "actual_sha256",
    "note",
]


def hash_file(path: Path) -> tuple[str, str]:
    md5 = hashlib.md5(usedforsecurity=False)
    sha256 = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            md5.update(chunk)
            sha256.update(chunk)
    return md5.hexdigest(), sha256.hexdigest()


def verify_manifest(package_root: Path, manifest_path: Path) -> list[dict[str, object]]:
    with manifest_path.open("r", encoding="utf-8-sig", newline="") as handle:
        expected_rows = list(csv.DictReader(handle, delimiter="\t"))
    results: list[dict[str, object]] = []
    for index, expected in enumerate(expected_rows, start=1):
        relative_path = expected["path"]
        path = package_root / Path(relative_path)
        print(f"[{index:02d}/{len(expected_rows):02d}] verifying {relative_path}", flush=True)
        result: dict[str, object] = {
            "role": expected["role"],
            "path": relative_path,
            "status": "FAIL",
            "expected_bytes": expected["expected_bytes"],
            "actual_bytes": "",
            "expected_md5": expected["expected_md5"],
            "actual_md5": "",
            "expected_sha256": expected["expected_sha256"],
            "actual_sha256": "",
            "note": "",
        }
        if not path.is_file():
            result["note"] = "missing file"
            results.append(result)
            continue
        actual_bytes = path.stat().st_size
        result["actual_bytes"] = actual_bytes
        if actual_bytes != int(expected["expected_bytes"]):
            result["note"] = "byte-size mismatch"
            results.append(result)
            continue
        actual_md5, actual_sha256 = hash_file(path)
        result["actual_md5"] = actual_md5
        result["actual_sha256"] = actual_sha256
        errors = []
        if expected["expected_md5"] and actual_md5.lower() != expected["expected_md5"].lower():
            errors.append("MD5 mismatch")
        if expected["expected_sha256"] and actual_sha256.lower() != expected["expected_sha256"].lower():
            errors.append("SHA-256 mismatch")
        if errors:
            result["note"] = "; ".join(errors)
        else:
            result["status"] = "PASS"
            result["note"] = "size, MD5, and SHA-256 verified"
        results.append(result)
    return results


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    results = verify_manifest(args.package_root.resolve(), args.manifest.resolve())
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with args.report.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    failures = [row for row in results if row["status"] != "PASS"]
    print(f"Input verification: {len(results) - len(failures)}/{len(results)} PASS")
    if failures:
        for row in failures:
            print(f"FAIL: {row['path']}: {row['note']}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
