#!/usr/bin/env python3
"""Build the audited input manifests for the Linux raw-read reconstruction.

This script does not copy the large FASTQ/reference files.  It validates the
two ENA manifests against a fixed run-to-donor map and writes compact manifests
consumed by the remote Linux verification and reconstruction scripts.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from collections import Counter
from pathlib import Path


NEW_RUN_TO_SAMPLE = {
    "SRR30622751": "PDR_S5",
    "SRR30622750": "PDR_S6",
    "SRR30622746": "PDR_S7",
    "SRR30622745": "PDR_S8",
    "SRR30622747": "PDR_S9",
    "SRR30622744": "PDR_S10",
    "SRR30622748": "PDR_S11",
    "SRR30622743": "PDR_S12",
    "SRR30622749": "MP_S13",
    "SRR30622742": "MH_S17",
}

OLD_SAMPLE_TO_RUNS = {
    "MP_S14": ["SRR11435044", "SRR11435045", "SRR11435046"],
    "MP_S15": ["SRR11435047", "SRR11435048", "SRR11435049"],
    "MH_S18": ["SRR11435050", "SRR11435051", "SRR11435052"],
    "MH_S19": ["SRR11435053", "SRR11435054", "SRR11435055"],
    "MH_S20": ["SRR11435056", "SRR11435057", "SRR11435058"],
    "MP_S16": ["SRR11435062", "SRR11435063", "SRR11435064"],
    "MH_S21": ["SRR11435065", "SRR11435066", "SRR11435067"],
}

SAMPLE_ORDER = [
    "PDR_S5",
    "PDR_S6",
    "PDR_S7",
    "PDR_S8",
    "PDR_S9",
    "PDR_S10",
    "PDR_S11",
    "PDR_S12",
    "MP_S13",
    "MP_S14",
    "MP_S15",
    "MP_S16",
    "MH_S17",
    "MH_S18",
    "MH_S19",
    "MH_S20",
    "MH_S21",
]

SAMPLE_SHEET_FIELDS = [
    "run_accession",
    "sample_id",
    "diagnosis",
    "source_dataset",
    "is_reused_control",
    "lane_index",
    "lanes_per_sample",
    "fastq_relative_path",
    "expected_bytes",
    "expected_md5",
    "expected_sha256",
    "official_fastq_url",
    "sample_accession",
    "experiment_accession",
    "secondary_sample_accession",
    "library_layout",
    "instrument_model",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def md5_file(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _index_unique(rows: list[dict[str, str]], key: str, label: str) -> dict[str, dict[str, str]]:
    values = [row[key] for row in rows]
    duplicates = sorted(value for value, count in Counter(values).items() if count > 1)
    if duplicates:
        raise ValueError(f"Duplicate {label}: {duplicates}")
    return {row[key]: row for row in rows}


def build_sample_sheet(
    new_manifest: list[dict[str, str]],
    old_manifest: list[dict[str, str]],
    integrity_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    new_by_run = _index_unique(new_manifest, "run_accession", "new run accession")
    old_by_run = _index_unique(old_manifest, "run_accession", "old run accession")
    expected_old_runs = {
        run for runs in OLD_SAMPLE_TO_RUNS.values() for run in runs
    }
    if set(new_by_run) != set(NEW_RUN_TO_SAMPLE):
        raise ValueError(
            "PRJNA1159345 manifest run set differs from the fixed 10-run map: "
            f"observed={sorted(new_by_run)}, expected={sorted(NEW_RUN_TO_SAMPLE)}"
        )
    if set(old_by_run) != expected_old_runs:
        raise ValueError(
            "GSE147657 manifest must contain exactly the 21 required control lanes: "
            f"observed={sorted(old_by_run)}, expected={sorted(expected_old_runs)}"
        )

    integrity_by_path = _index_unique(integrity_rows, "path", "integrity-report path")
    run_to_old_sample = {
        run: sample for sample, runs in OLD_SAMPLE_TO_RUNS.items() for run in runs
    }
    run_order = sorted(
        NEW_RUN_TO_SAMPLE,
        key=lambda run: (SAMPLE_ORDER.index(NEW_RUN_TO_SAMPLE[run]), run),
    )
    run_order.extend(
        run
        for sample in SAMPLE_ORDER
        for run in OLD_SAMPLE_TO_RUNS.get(sample, [])
    )

    output: list[dict[str, object]] = []
    for run in run_order:
        is_new = run in new_by_run
        manifest_row = new_by_run[run] if is_new else old_by_run[run]
        sample = NEW_RUN_TO_SAMPLE[run] if is_new else run_to_old_sample[run]
        source = "GSE276892_new" if is_new else "GSE147657_reused"
        folder = "GSE276892" if is_new else "GSE147657"
        relative_path = (
            "analysis_data/independent_validation/raw_reads/"
            f"{folder}/{run}.fastq.gz"
        )
        integrity = integrity_by_path.get(relative_path)
        if integrity is None or integrity.get("status") != "PASS":
            raise ValueError(f"Missing PASS integrity record for {relative_path}")
        if integrity.get("actual_md5") != manifest_row["fastq_md5"]:
            raise ValueError(f"MD5 mismatch between manifest and integrity report for {run}")

        if is_new:
            lane_index = 1
            lanes_per_sample = 1
        else:
            lanes = OLD_SAMPLE_TO_RUNS[sample]
            lane_index = lanes.index(run) + 1
            lanes_per_sample = len(lanes)
        output.append(
            {
                "run_accession": run,
                "sample_id": sample,
                "diagnosis": "PDR" if sample.startswith("PDR_") else "control",
                "source_dataset": source,
                "is_reused_control": 0 if is_new else 1,
                "lane_index": lane_index,
                "lanes_per_sample": lanes_per_sample,
                "fastq_relative_path": relative_path,
                "expected_bytes": manifest_row["fastq_bytes"],
                "expected_md5": manifest_row["fastq_md5"],
                "expected_sha256": integrity["sha256"],
                "official_fastq_url": manifest_row["fastq_ftp"],
                "sample_accession": manifest_row["sample_accession"],
                "experiment_accession": manifest_row["experiment_accession"],
                "secondary_sample_accession": manifest_row[
                    "secondary_sample_accession"
                ],
                "library_layout": manifest_row["library_layout"],
                "instrument_model": manifest_row["instrument_model"],
            }
        )
    return output


def build_input_manifest(
    package_root: Path,
    sample_rows: list[dict[str, object]],
    integrity_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    integrity_by_path = {row["path"]: row for row in integrity_rows}
    output: list[dict[str, object]] = []
    for row in sample_rows:
        output.append(
            {
                "role": "FASTQ",
                "path": row["fastq_relative_path"],
                "expected_bytes": row["expected_bytes"],
                "expected_md5": row["expected_md5"],
                "expected_sha256": row["expected_sha256"],
            }
        )

    selected = [
        ("GENCODE_all_regions_FASTA", "tools/reference/gencode_v42_all_regions/GRCh38.p13.genome.fa.gz"),
        ("GENCODE_all_regions_GTF", "tools/reference/gencode_v42_all_regions/gencode.v42.chr_patch_hapl_scaff.annotation.gtf.gz"),
        ("GENCODE_primary_FASTA", "tools/reference/gencode_v42_primary_assembly/GRCh38.primary_assembly.genome.fa.gz"),
        ("GENCODE_primary_GTF", "tools/reference/gencode_v42_primary_assembly/gencode.v42.primary_assembly.annotation.gtf.gz"),
        ("STAR_2.7.8a_archive", "tools/downloads/STAR-2.7.8a.tar.gz"),
    ]
    for role, relative_path in selected:
        integrity = integrity_by_path.get(relative_path)
        if integrity is None or integrity.get("status") != "PASS":
            raise ValueError(f"Missing PASS integrity record for {relative_path}")
        output.append(
            {
                "role": role,
                "path": relative_path,
                "expected_bytes": integrity["actual_bytes"],
                "expected_md5": integrity["actual_md5"],
                "expected_sha256": integrity["sha256"],
            }
        )

    linux_subread = package_root / "tools/downloads/subread-2.0.1-Linux-x86_64.tar.gz"
    if not linux_subread.is_file():
        raise FileNotFoundError(linux_subread)
    output.append(
        {
            "role": "Subread_2.0.1_Linux_archive",
            "path": linux_subread.relative_to(package_root).as_posix(),
            "expected_bytes": linux_subread.stat().st_size,
            "expected_md5": md5_file(linux_subread),
            "expected_sha256": sha256_file(linux_subread),
        }
    )
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--package-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
    )
    args = parser.parse_args()
    package_root = args.package_root.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else package_root / "remote_linux_reconstruction/config"
    )

    validation_dir = package_root / "analysis_data/independent_validation"
    new_manifest = read_tsv(validation_dir / "PRJNA1159345_ENA_run_manifest_full.tsv")
    old_manifest = read_tsv(validation_dir / "GSE147657_ENA_run_manifest.tsv")
    integrity_rows = read_csv(
        package_root / "analysis_results/download_integrity_report.csv"
    )
    sample_rows = build_sample_sheet(new_manifest, old_manifest, integrity_rows)
    input_rows = build_input_manifest(package_root, sample_rows, integrity_rows)

    write_tsv(output_dir / "sample_sheet.tsv", sample_rows, SAMPLE_SHEET_FIELDS)
    write_tsv(
        output_dir / "input_files.tsv",
        input_rows,
        ["role", "path", "expected_bytes", "expected_md5", "expected_sha256"],
    )
    write_tsv(
        output_dir / "provenance_sources.tsv",
        [
            {
                "record": "PRJNA1159345 run-to-sample mapping",
                "source_url": "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1159345&result=read_run",
                "accessed_on": "2026-08-14",
                "note": "ENA sample_title supplied the fixed 10-run mapping.",
            },
            {
                "record": "GSE147657 reused-control lane mapping",
                "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276892",
                "accessed_on": "2026-08-14",
                "note": "Seven controls, each represented by three technical lanes, were retained.",
            },
            {
                "record": "GENCODE human release 42",
                "source_url": "https://www.gencodegenes.org/human/release_42.html",
                "accessed_on": "2026-08-14",
                "note": "GRCh38.p13 all-regions source reconstruction and primary-assembly sensitivity.",
            },
        ],
        ["record", "source_url", "accessed_on", "note"],
    )

    print(f"Wrote {len(sample_rows)} run rows for {len(set(row['sample_id'] for row in sample_rows))} biological samples")
    print(f"Wrote {len(input_rows)} required input-file records to {output_dir}")


if __name__ == "__main__":
    main()
