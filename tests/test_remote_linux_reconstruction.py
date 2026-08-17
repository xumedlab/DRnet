from __future__ import annotations

import csv
import importlib.util
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).parents[1]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PREPARE = load_module(
    "prepare_remote_linux_reconstruction",
    ROOT / "analysis_scripts" / "37_prepare_remote_linux_reconstruction.py",
)
VERIFY = load_module(
    "verify_remote_inputs",
    ROOT / "remote_linux_reconstruction" / "scripts" / "verify_remote_inputs.py",
)
AGGREGATE = load_module(
    "aggregate_featurecounts",
    ROOT / "remote_linux_reconstruction" / "scripts" / "aggregate_featurecounts.py",
)


def manifest_row(run: str, md5: str = "a" * 32) -> dict[str, str]:
    return {
        "run_accession": run,
        "sample_accession": f"SAMN_{run}",
        "experiment_accession": f"SRX_{run}",
        "secondary_sample_accession": f"SRS_{run}",
        "fastq_ftp": f"https://example.invalid/{run}.fastq.gz",
        "fastq_bytes": "100",
        "fastq_md5": md5,
        "library_layout": "SINGLE",
        "instrument_model": "test instrument",
    }


def integrity_row(run: str, is_new: bool, md5: str = "a" * 32) -> dict[str, str]:
    folder = "GSE276892" if is_new else "GSE147657"
    return {
        "path": (
            "analysis_data/independent_validation/raw_reads/"
            f"{folder}/{run}.fastq.gz"
        ),
        "status": "PASS",
        "actual_md5": md5,
        "sha256": "b" * 64,
    }


def build_sample_rows() -> list[dict[str, object]]:
    new_manifest = [manifest_row(run) for run in PREPARE.NEW_RUN_TO_SAMPLE]
    old_runs = [run for runs in PREPARE.OLD_SAMPLE_TO_RUNS.values() for run in runs]
    old_manifest = [manifest_row(run) for run in old_runs]
    integrity = [integrity_row(run, True) for run in PREPARE.NEW_RUN_TO_SAMPLE]
    integrity.extend(integrity_row(run, False) for run in old_runs)
    return PREPARE.build_sample_sheet(new_manifest, old_manifest, integrity)


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_fixed_mapping_has_31_lanes_and_17_biological_samples() -> None:
    rows = build_sample_rows()
    lane_counts = Counter(str(row["sample_id"]) for row in rows)

    assert len(rows) == 31
    assert len(lane_counts) == 17
    assert lane_counts["PDR_S5"] == 1
    assert lane_counts["MP_S14"] == 3
    assert lane_counts["MH_S21"] == 3
    assert sum(row["diagnosis"] == "PDR" for row in rows) == 8
    assert sum(row["is_reused_control"] == 1 for row in rows) == 21


def test_remote_verifier_detects_hash_and_size(tmp_path: Path) -> None:
    payload = tmp_path / "payload.bin"
    payload.write_bytes(b"DRnet remote reconstruction")
    md5, sha256 = VERIFY.hash_file(payload)
    manifest = tmp_path / "input_files.tsv"
    write_tsv(
        manifest,
        [
            {
                "role": "test",
                "path": "payload.bin",
                "expected_bytes": payload.stat().st_size,
                "expected_md5": md5,
                "expected_sha256": sha256,
            }
        ],
    )

    results = VERIFY.verify_manifest(tmp_path, manifest)

    assert results[0]["status"] == "PASS"


def test_lane_counts_are_summed_to_donor_level(tmp_path: Path) -> None:
    sample_rows = build_sample_rows()
    sample_sheet = tmp_path / "sample_sheet.tsv"
    write_tsv(sample_sheet, sample_rows)
    workflow = tmp_path / "workflow"
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text(
        '1\ttest\tgene\t1\t10\t.\t+\t.\tgene_id "ENSG_P2RX4.1"; '
        'gene_type "protein_coding"; gene_name "P2RX4";\n'
        '1\ttest\tgene\t20\t30\t.\t+\t.\tgene_id "ENSG_OTHER.1"; '
        'gene_type "protein_coding"; gene_name "OTHER";\n',
        encoding="utf-8",
    )
    for row in sample_rows:
        run = str(row["run_accession"])
        counts_path = workflow / "lane_counts" / f"{run}.featureCounts.txt"
        counts_path.parent.mkdir(parents=True, exist_ok=True)
        counts_path.write_text(
            "# Program:featureCounts v2.0.1\n"
            "Geneid\tChr\tStart\tEnd\tStrand\tLength\ttest.bam\n"
            "ENSG_P2RX4.1\t1\t1\t10\t+\t10\t1\n"
            "ENSG_OTHER.1\t1\t20\t30\t+\t11\t2\n",
            encoding="utf-8",
        )
        Path(f"{counts_path}.summary").write_text(
            "Status\ttest.bam\nAssigned\t3\nUnassigned_NoFeatures\t2\n",
            encoding="utf-8",
        )
        log_path = workflow / "align" / run / "Log.final.out"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(
            "Number of input reads | 100\n"
            "Uniquely mapped reads number | 80\n"
            "Uniquely mapped reads % | 80.00%\n"
            "% of reads mapped to multiple loci | 5.00%\n"
            "% of reads unmapped: too short | 10.00%\n",
            encoding="utf-8",
        )

    output = tmp_path / "output"
    AGGREGATE.aggregate(sample_sheet, workflow, gtf, output)
    rows = {row["sample_id"]: row for row in AGGREGATE.read_tsv(output / "sample_qc.tsv")}
    p2rx4 = AGGREGATE.read_tsv(output / "P2RX4_counts.tsv")[0]

    assert rows["MP_S14"]["number_of_lanes"] == "3"
    assert rows["MP_S14"]["input_reads"] == "300"
    assert p2rx4["PDR_S5"] == "1"
    assert p2rx4["MP_S14"] == "3"
