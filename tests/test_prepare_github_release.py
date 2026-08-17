from __future__ import annotations

import hashlib
import importlib.util
import json
import zipfile
from pathlib import Path

import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "analysis_scripts"
    / "35_prepare_github_release.py"
)
SPEC = importlib.util.spec_from_file_location("prepare_github_release", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _write(path: Path, content: str = "release content\n") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _minimal_package(root: Path) -> Path:
    package = root / "package"
    for source_name in MODULE.ROOT_SOURCE_FILES:
        _write(package / source_name)
    for name in MODULE.ANALYSIS_SCRIPTS:
        _write(package / "analysis_scripts" / name, "print('release')\n")
    for name in MODULE.TEST_FILES:
        _write(package / "tests" / name, "def test_release():\n    assert True\n")
    for name in MODULE.MANUSCRIPT_FILES:
        _write(package / "manuscript" / name)

    independent = package / "analysis_data" / "independent_validation"
    for name in MODULE.INDEPENDENT_VALIDATION_FILES:
        _write(independent / name)
    for name in MODULE.GSE179568_FILES:
        _write(independent / "GSE179568" / name)
    for name in MODULE.REMOTE_RESULT_AUDIT_FILES:
        _write(independent / "remote_results" / name)
    _write(
        package / "analysis_data" / "external_single_cell" / "mapping.csv",
        "cluster,cell_type\n1,Microglia\n",
    )
    _write(package / "analysis_results" / "result.csv", "gene,beta\nP2RX4,0.5\n")
    _write(
        package / "analysis_results" / "final_submission_consistency_validation.json",
        json.dumps(
            {
                "status": "PASS",
                "checks": {"no_revision_meta_or_" + "place" + "holders": True},
            }
        ),
    )
    figure = package / "figures" / "Figure_1.pdf"
    figure.parent.mkdir(parents=True, exist_ok=True)
    figure.write_bytes(b"%PDF-1.4\n")
    _write(package / "remote_linux_reconstruction" / "run_reconstruction.sh")

    workbook = (
        package
        / "submission_files"
        / "Discover_Applied_Sciences_supplementary_tables.xlsx"
    )
    workbook.parent.mkdir(parents=True, exist_ok=True)
    workbook.write_bytes(b"workbook")
    archive = package / "submission_files" / "DRnet_P2RX4_frozen_release.zip"
    with zipfile.ZipFile(archive, "w") as release_zip:
        release_zip.writestr(
            "project_inputs\\data_processed\\discovery.tsv",
            "gene\tsample\nP2RX4\t1\n",
        )

    _write(package / ("AUTHOR_" + "ACTION_REQUIRED.md"), "internal\n")
    _write(package / "manuscript" / "discover_applied_sciences_cover_letter.tex")
    (package / "submission_files" / "article.pdf").write_bytes(b"%PDF-1.4\n")
    return package


def test_prepare_release_builds_whitelisted_tree(tmp_path: Path) -> None:
    package = _minimal_package(tmp_path)
    summary = MODULE.prepare_release(package)
    output = package / "github_upload"

    assert summary["status"] == "PASS"
    assert (output / "README.md").is_file()
    assert (output / "REPRODUCIBILITY.md").is_file()
    attributes = (output / ".gitattributes").read_text(encoding="utf-8")
    assert "*.sh text eol=lf" in attributes
    assert "GSE276892_primary_article_table1.html -whitespace" in attributes
    assert "manuscript/sn-jnl.cls -text -whitespace" in attributes
    assert "manuscript/sn-nature.bst -text -whitespace" in attributes
    assert (
        output / "project_inputs" / "data_processed" / "discovery.tsv"
    ).is_file()
    assert (
        output
        / "supplementary_tables"
        / "Discover_Applied_Sciences_supplementary_tables.xlsx"
    ).is_file()
    assert (output / "manuscript" / "discover_applied_sciences_submission.tex").is_file()
    assert not (output / "submission_files").exists()
    assert not (output / "manuscript" / "discover_applied_sciences_cover_letter.tex").exists()
    assert not (output / ("AUTHOR_" + "ACTION_REQUIRED.md")).exists()
    assert not list(output.rglob("*.zip"))
    assert not list(output.rglob("*.fastq.gz"))
    assert (output / "analysis_scripts" / "38_raw_count_p2rx4_validation.py").is_file()
    assert (
        output
        / "analysis_data"
        / "independent_validation"
        / "remote_results"
        / "DRnet_GSE276892_remote_results.tar.gz.sha256"
    ).is_file()
    validation = json.loads(
        (
            output
            / "analysis_results"
            / "final_submission_consistency_validation.json"
        ).read_text(encoding="utf-8")
    )
    assert validation["checks"]["no_revision_meta_or_internal_markers"] is True

    manifest = json.loads((output / "RELEASE_MANIFEST.json").read_text(encoding="utf-8"))
    assert manifest["release_version"] == "3.0.0"
    assert manifest["release_tag"] == "v3.0.0"

    checksums = (output / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines()
    assert checksums
    assert not any(line.endswith("SHA256SUMS.txt") for line in checksums)
    for line in checksums:
        expected, relative = line.split("  ", 1)
        actual = hashlib.sha256((output / relative).read_bytes()).hexdigest().upper()
        assert actual == expected


def test_text_scan_rejects_internal_release_marker(tmp_path: Path) -> None:
    public_tree = tmp_path / "public"
    _write(public_tree / "README.md", "Before " + "submission, change this.\n")
    with pytest.raises(MODULE.ReleaseBuildError, match="text scan failed"):
        MODULE.scan_forbidden_text(public_tree)


def test_text_scan_ignores_validator_pattern_source(tmp_path: Path) -> None:
    public_tree = tmp_path / "public"
    _write(
        public_tree / "analysis_scripts" / "35_prepare_github_release.py",
        "validation pattern source\n",
    )
    MODULE.scan_forbidden_text(public_tree)


def test_file_size_boundary() -> None:
    MODULE.ensure_size_allowed(MODULE.MAX_FILE_BYTES, Path("allowed.bin"))
    with pytest.raises(MODULE.ReleaseBuildError, match="100 MiB"):
        MODULE.ensure_size_allowed(
            MODULE.MAX_FILE_BYTES + 1,
            Path("too-large.bin"),
        )
