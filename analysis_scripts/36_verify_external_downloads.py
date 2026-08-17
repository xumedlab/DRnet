#!/usr/bin/env python3
"""Verify downloaded FASTQ, reference, tool, PDF, and spreadsheet artifacts.

The FASTQ files are read once while MD5 and SHA-256 are updated in parallel.
Expected FASTQ sizes and MD5 values come from ENA manifests; expected GENCODE
MD5 values come from the release 42 ``MD5SUMS`` file downloaded by the author.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import tarfile
import zipfile
from pathlib import Path
from typing import Any


FASTQ_MANIFESTS = (
    ("PRJNA1159345_ENA_run_manifest_full.tsv", "GSE276892"),
    ("GSE147657_ENA_run_manifest.tsv", "GSE147657"),
)

REFERENCE_FILES = {
    "GRCh38.p13.genome.fa.gz": "tools/reference/gencode_v42_all_regions/GRCh38.p13.genome.fa.gz",
    "gencode.v42.chr_patch_hapl_scaff.annotation.gtf.gz": (
        "tools/reference/gencode_v42_all_regions/"
        "gencode.v42.chr_patch_hapl_scaff.annotation.gtf.gz"
    ),
    "GRCh38.primary_assembly.genome.fa.gz": (
        "tools/reference/gencode_v42_primary_assembly/"
        "GRCh38.primary_assembly.genome.fa.gz"
    ),
    "gencode.v42.primary_assembly.annotation.gtf.gz": (
        "tools/reference/gencode_v42_primary_assembly/"
        "gencode.v42.primary_assembly.annotation.gtf.gz"
    ),
}


def hash_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> tuple[str, str]:
    """Return MD5 and SHA-256 from a single sequential read."""
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(chunk_size):
            md5.update(chunk)
            sha256.update(chunk)
    return md5.hexdigest(), sha256.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def read_gencode_md5(path: Path) -> dict[str, str]:
    expected: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        fields = line.strip().split()
        if len(fields) >= 2 and len(fields[0]) == 32:
            expected[Path(fields[-1].lstrip("*./")).name] = fields[0].lower()
    return expected


def make_record(
    *,
    category: str,
    path: Path,
    root: Path,
    status: str,
    expected_bytes: int | None = None,
    actual_bytes: int | None = None,
    expected_md5: str = "",
    actual_md5: str = "",
    sha256: str = "",
    note: str = "",
) -> dict[str, Any]:
    return {
        "category": category,
        "path": path.relative_to(root).as_posix(),
        "status": status,
        "expected_bytes": expected_bytes,
        "actual_bytes": actual_bytes,
        "expected_md5": expected_md5,
        "actual_md5": actual_md5,
        "sha256": sha256,
        "note": note,
    }


def verify_fastqs(root: Path) -> list[dict[str, Any]]:
    validation = root / "analysis_data" / "independent_validation"
    records: list[dict[str, Any]] = []
    for manifest_name, dataset in FASTQ_MANIFESTS:
        manifest = validation / manifest_name
        rows = read_tsv(manifest)
        runs = [row["run_accession"] for row in rows]
        if len(runs) != len(set(runs)):
            raise ValueError(f"Duplicate run_accession values in {manifest}")
        fastq_dir = validation / "raw_reads" / dataset
        expected_names = {f"{run}.fastq.gz" for run in runs}
        actual_names = {path.name for path in fastq_dir.glob("*.fastq.gz")}
        extras = sorted(actual_names - expected_names)
        missing = sorted(expected_names - actual_names)
        if extras or missing:
            raise ValueError(
                f"FASTQ set mismatch for {dataset}: missing={missing}, extras={extras}"
            )
        for row in sorted(rows, key=lambda item: item["run_accession"]):
            fastq = fastq_dir / f"{row['run_accession']}.fastq.gz"
            expected_bytes = int(row["fastq_bytes"])
            expected_md5 = row["fastq_md5"].lower()
            actual_bytes = fastq.stat().st_size
            actual_md5, sha256 = hash_file(fastq)
            status = (
                "PASS"
                if actual_bytes == expected_bytes and actual_md5 == expected_md5
                else "FAIL"
            )
            records.append(
                make_record(
                    category="FASTQ",
                    path=fastq,
                    root=root,
                    status=status,
                    expected_bytes=expected_bytes,
                    actual_bytes=actual_bytes,
                    expected_md5=expected_md5,
                    actual_md5=actual_md5,
                    sha256=sha256,
                    note=dataset,
                )
            )
    return records


def verify_references(root: Path) -> list[dict[str, Any]]:
    md5_path = (
        root
        / "tools"
        / "reference"
        / "gencode_v42_primary_assembly"
        / "release 42 MD5SUMS.txt"
    )
    expected = read_gencode_md5(md5_path)
    records: list[dict[str, Any]] = []
    for source_name, relative_path in REFERENCE_FILES.items():
        path = root / relative_path
        expected_md5 = expected.get(source_name, "")
        if not expected_md5:
            raise ValueError(f"No official MD5 found for {source_name}")
        actual_md5, sha256 = hash_file(path)
        records.append(
            make_record(
                category="GENCODE_v42",
                path=path,
                root=root,
                status="PASS" if actual_md5 == expected_md5 else "FAIL",
                actual_bytes=path.stat().st_size,
                expected_md5=expected_md5,
                actual_md5=actual_md5,
                sha256=sha256,
                note="official release 42 MD5SUMS",
            )
        )
    return records


def verify_archives(root: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    download_dir = root / "tools" / "downloads"
    for path in sorted(download_dir.iterdir()):
        if not path.is_file():
            continue
        actual_md5, sha256 = hash_file(path)
        status = "PASS"
        note = ""
        try:
            if path.name.endswith(".tar.gz"):
                with tarfile.open(path, "r:gz") as archive:
                    members = archive.getmembers()
                    if not members:
                        raise ValueError("empty tar archive")
                    note = f"tar members={len(members)}"
            elif path.suffix.lower() == ".zip":
                with zipfile.ZipFile(path) as archive:
                    bad = archive.testzip()
                    if bad:
                        raise ValueError(f"CRC failure: {bad}")
                    note = f"zip members={len(archive.infolist())}"
            else:
                status = "NOT_CHECKED"
                note = "unsupported archive suffix"
        except (OSError, tarfile.TarError, zipfile.BadZipFile, ValueError) as exc:
            status = "FAIL"
            note = str(exc)
        records.append(
            make_record(
                category="tool_archive",
                path=path,
                root=root,
                status=status,
                actual_bytes=path.stat().st_size,
                actual_md5=actual_md5,
                sha256=sha256,
                note=note,
            )
        )
    return records


def verify_documents(root: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    directory = root / "analysis_data" / "independent_validation" / "GSE179568"
    for path in sorted(directory.iterdir()):
        if not path.is_file() or path.suffix.lower() not in {".pdf", ".xlsx"}:
            continue
        actual_md5, sha256 = hash_file(path)
        status = "PASS"
        note = ""
        if path.suffix.lower() == ".pdf":
            header = path.read_bytes()[:5]
            if header != b"%PDF-":
                status = "FAIL"
                note = f"invalid PDF header: {header!r}"
            else:
                note = "valid %PDF- header"
        else:
            try:
                with zipfile.ZipFile(path) as archive:
                    bad = archive.testzip()
                    if bad:
                        raise ValueError(f"CRC failure: {bad}")
                    note = f"valid XLSX zip; members={len(archive.infolist())}"
            except (OSError, zipfile.BadZipFile, ValueError) as exc:
                status = "FAIL"
                note = str(exc)
        records.append(
            make_record(
                category="GSE179568_document",
                path=path,
                root=root,
                status=status,
                actual_bytes=path.stat().st_size,
                actual_md5=actual_md5,
                sha256=sha256,
                note=note,
            )
        )
    return records


def verify_prjna744210_manifest(root: Path) -> list[dict[str, Any]]:
    path = (
        root
        / "analysis_data"
        / "independent_validation"
        / "filereport_read_run_PRJNA744210.tsv"
    )
    rows = read_tsv(path)
    runs = [row["run_accession"] for row in rows]
    required = {
        "run_accession",
        "sample_accession",
        "fastq_ftp",
        "fastq_bytes",
        "fastq_md5",
        "library_layout",
    }
    columns = set(rows[0]) if rows else set()
    status = (
        "PASS"
        if rows and required <= columns and len(runs) == len(set(runs))
        else "FAIL"
    )
    actual_md5, sha256 = hash_file(path)
    return [
        make_record(
            category="ENA_manifest",
            path=path,
            root=root,
            status=status,
            actual_bytes=path.stat().st_size,
            actual_md5=actual_md5,
            sha256=sha256,
            note=f"rows={len(rows)}; unique_runs={len(set(runs))}",
        )
    ]


def write_reports(root: Path, records: list[dict[str, Any]]) -> None:
    output_dir = root / "analysis_results"
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "download_integrity_report.csv"
    json_path = output_dir / "download_integrity_report.json"
    fieldnames = list(records[0])
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
    summary = {
        "record_count": len(records),
        "pass_count": sum(record["status"] == "PASS" for record in records),
        "fail_count": sum(record["status"] == "FAIL" for record in records),
        "not_checked_count": sum(
            record["status"] == "NOT_CHECKED" for record in records
        ),
        "records": records,
    }
    json_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps({key: value for key, value in summary.items() if key != "records"}))
    print(csv_path)
    print(json_path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Submission package root",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    records = []
    records.extend(verify_fastqs(root))
    records.extend(verify_references(root))
    records.extend(verify_archives(root))
    records.extend(verify_documents(root))
    records.extend(verify_prjna744210_manifest(root))
    write_reports(root, records)
    return 1 if any(record["status"] == "FAIL" for record in records) else 0


if __name__ == "__main__":
    raise SystemExit(main())
