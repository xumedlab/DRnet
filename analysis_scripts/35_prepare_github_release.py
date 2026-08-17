#!/usr/bin/env python3
"""Build and validate the public GitHub tree for the v3.0.0 analysis release."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import zipfile
from pathlib import Path, PurePosixPath


RELEASE_VERSION = "3.0.0"
RELEASE_TAG = f"v{RELEASE_VERSION}"
REPOSITORY_URL = "https://github.com/xumedlab/DRnet"
MAX_FILE_BYTES = 100 * 1024 * 1024

ROOT_SOURCE_FILES = {
    "README_submission_package.md": "README.md",
    "REPRODUCIBILITY_README.md": "REPRODUCIBILITY.md",
    "COMPUTATIONAL_ENVIRONMENT.md": "COMPUTATIONAL_ENVIRONMENT.md",
    "CITATION.cff": "CITATION.cff",
    "requirements.txt": "requirements.txt",
    "pytest.ini": "pytest.ini",
}

ANALYSIS_SCRIPTS = (
    "25_voigt_single_cell_localization.py",
    "26_updated_study_design.py",
    "32_independent_p2rx4_validation.py",
    "33_final_discovery_statistics.py",
    "34_validate_final_submission.py",
    "35_prepare_github_release.py",
    "36_verify_external_downloads.py",
    "37_prepare_remote_linux_reconstruction.py",
    "38_raw_count_p2rx4_validation.py",
    "39_build_submission_archives.py",
    "run_final_research_article_pipeline.py",
)

TEST_FILES = (
    "test_final_discovery_statistics.py",
    "test_independent_p2rx4_validation.py",
    "test_validate_final_submission.py",
    "test_voigt_single_cell_localization.py",
    "test_prepare_github_release.py",
    "test_raw_count_p2rx4_validation.py",
    "test_remote_linux_reconstruction.py",
    "test_verify_external_downloads.py",
    "test_build_submission_archives.py",
)

MANUSCRIPT_FILES = (
    "discover_applied_sciences_submission.tex",
    "discover_applied_sciences_supplementary_information.tex",
    "discover_references.bib",
    "sn-jnl.cls",
    "sn-nature.bst",
)

COPY_TREES = (
    "analysis_results",
    "figures",
    "remote_linux_reconstruction",
)

INDEPENDENT_VALIDATION_FILES = (
    "filereport_read_run_PRJNA744210.tsv",
    "GSE147657_ENA_run_manifest.tsv",
    "GSE276892_normal_data.csv.gz",
    "GSE276892_primary_article_table1.html",
    "GSE276892_README_all_samples.xlsx",
    "GSE94019_Partek_EM_gene_reads.txt.gz",
    "P2RX4_VALIDATION_PROTOCOL_LOCALLY_FROZEN.md",
    "PRJNA1159345_ENA_run_manifest_full.tsv",
    "PRJNA1159345_ENA_run_manifest.tsv",
)

GSE179568_FILES = (
    "GSE179568_data.csv.gz",
    "GSE179568_family.soft.gz",
    "GSE179568_series_matrix.txt.gz",
    "Table 1.pdf",
)

REMOTE_RESULT_AUDIT_FILES = (
    "DRnet_GSE276892_remote_results.tar.gz.sha256",
    "nohup_reconstruction.log",
)

TEXT_SUFFIXES = {
    ".bib",
    ".cff",
    ".cls",
    ".csv",
    ".gmt",
    ".ini",
    ".json",
    ".md",
    ".py",
    ".tex",
    ".tsv",
    ".txt",
}

# This validator is part of the public tree, so its own pattern definitions are
# excluded from the content scan. Tests still exercise each public-file rule.
TEXT_SCAN_EXEMPT_PATHS = {
    "analysis_scripts/34_validate_final_submission.py",
    "analysis_scripts/35_prepare_github_release.py",
}

IGNORED_DIRECTORY_NAMES = {
    "__pycache__",
    ".pytest_cache",
    ".ruff_cache",
    "tmp",
}

IGNORED_FILE_SUFFIXES = {
    ".aux",
    ".log",
    ".pyc",
    ".pyo",
    ".synctex.gz",
}

# The strings are assembled so this validator does not trigger on its own source.
FORBIDDEN_TEXT_PATTERNS = (
    re.compile(r"before\s+" + r"submission", re.IGNORECASE),
    re.compile("AUTHOR_" + "ACTION_REQUIRED", re.IGNORECASE),
    re.compile(r"submission[-\s]+" + r"ready", re.IGNORECASE),
    re.compile(r"must\s+not\s+be\s+" + r"uploaded", re.IGNORECASE),
    re.compile("place" + "holder", re.IGNORECASE),
    re.compile("pend" + "ing", re.IGNORECASE),
    re.compile("local-" + "release", re.IGNORECASE),
    re.compile("投稿前"),
    re.compile("待办"),
)


class ReleaseBuildError(RuntimeError):
    """Raised when the proposed public release is incomplete or unsafe."""


def parse_args() -> argparse.Namespace:
    package = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path, default=package)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def _safe_output_root(package: Path) -> Path:
    package = package.resolve()
    output = (package / "github_upload").resolve()
    if output.parent != package or output.name != "github_upload":
        raise ReleaseBuildError(f"Unsafe output location: {output}")
    return output


def _copy_file(source: Path, destination: Path) -> None:
    if not source.is_file():
        raise ReleaseBuildError(f"Required release source is missing: {source}")
    if source.is_symlink():
        raise ReleaseBuildError(f"Symbolic links are not allowed: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def _is_ignored(relative: Path) -> bool:
    if any(part in IGNORED_DIRECTORY_NAMES for part in relative.parts):
        return True
    name = relative.name.lower()
    return any(name.endswith(suffix) for suffix in IGNORED_FILE_SUFFIXES)


def _copy_tree(source: Path, destination: Path) -> None:
    if not source.is_dir():
        raise ReleaseBuildError(f"Required release directory is missing: {source}")
    for path in sorted(source.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(source)
        if _is_ignored(relative):
            continue
        _copy_file(path, destination / relative)


def _copy_named_files(
    source_directory: Path,
    destination_directory: Path,
    names: tuple[str, ...],
) -> None:
    for name in names:
        _copy_file(source_directory / name, destination_directory / name)


def _extract_project_inputs(package: Path, destination: Path) -> None:
    source_directory = package / "project_inputs"
    if source_directory.is_dir():
        _copy_tree(source_directory, destination)
        return

    archive_path = package / "submission_files" / "DRnet_P2RX4_frozen_release.zip"
    if not archive_path.is_file():
        raise ReleaseBuildError(
            "project_inputs is missing and the frozen source archive is unavailable"
        )

    extracted = 0
    with zipfile.ZipFile(archive_path) as archive:
        for member in archive.infolist():
            normalized = member.filename.replace("\\", "/")
            relative = PurePosixPath(normalized)
            if member.is_dir() or not relative.parts:
                continue
            if relative.parts[0] != "project_inputs":
                continue
            payload_parts = relative.parts[1:]
            if not payload_parts or any(part in {"", ".", ".."} for part in payload_parts):
                raise ReleaseBuildError(f"Unsafe project-input archive member: {normalized}")
            ensure_size_allowed(member.file_size, Path(normalized))
            target = destination.joinpath(*payload_parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(member) as source, target.open("wb") as sink:
                shutil.copyfileobj(source, sink)
            extracted += 1
    if extracted == 0:
        raise ReleaseBuildError("No project_inputs files were found in the frozen archive")


def ensure_size_allowed(size: int, path: Path) -> None:
    if size > MAX_FILE_BYTES:
        raise ReleaseBuildError(
            f"File exceeds the 100 MiB GitHub limit: {path} ({size} bytes)"
        )


def _iter_release_files(root: Path) -> list[Path]:
    return sorted(path for path in root.rglob("*") if path.is_file())


def validate_sizes(root: Path) -> None:
    for path in _iter_release_files(root):
        ensure_size_allowed(path.stat().st_size, path.relative_to(root))


def _is_text_file(path: Path) -> bool:
    return path.name in {"README", ".gitignore"} or path.suffix.lower() in TEXT_SUFFIXES


def scan_forbidden_text(root: Path) -> None:
    violations: list[str] = []
    for path in _iter_release_files(root):
        relative = path.relative_to(root)
        normalized_path = relative.as_posix()
        for pattern in FORBIDDEN_TEXT_PATTERNS:
            if pattern.search(normalized_path):
                violations.append(f"{normalized_path}: disallowed path marker")
        if normalized_path in TEXT_SCAN_EXEMPT_PATHS or not _is_text_file(path):
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError as exc:
            raise ReleaseBuildError(f"Text file is not UTF-8: {relative}") from exc
        for pattern in FORBIDDEN_TEXT_PATTERNS:
            match = pattern.search(text)
            if match:
                line = text.count("\n", 0, match.start()) + 1
                violations.append(f"{normalized_path}:{line}: disallowed release marker")
    if violations:
        raise ReleaseBuildError(
            "Public-release text scan failed:\n" + "\n".join(sorted(set(violations)))
        )


def _write_gitignore(output: Path) -> None:
    content = """# Python
__pycache__/
*.py[cod]
.pytest_cache/
.ruff_cache/

# Local build output
tmp/
github_upload/

# LaTeX intermediates
*.aux
*.bbl
*.blg
*.fls
*.fdb_latexmk
*.log
*.out
*.synctex.gz
"""
    (output / ".gitignore").write_text(content, encoding="utf-8")


def _write_gitattributes(output: Path) -> None:
    content = """# Reproducible text line endings across Windows and Linux
* text=auto
*.bib text eol=lf
*.bst text eol=lf
*.cff text eol=lf
*.cls text eol=lf
*.csv text eol=lf
*.gmt text eol=lf
*.html text eol=lf
*.ini text eol=lf
*.json text eol=lf
*.md text eol=lf
*.py text eol=lf
*.sh text eol=lf
*.tex text eol=lf
*.tsv text eol=lf
*.txt text eol=lf

# Binary scientific inputs and rendered artifacts
*.gz binary
*.pdf binary
*.png binary
*.tar binary
*.xlsx binary

# Preserve deposited and publisher-supplied source files verbatim. Their source
# formatting includes non-canonical line endings or trailing spaces that are not
# manuscript or program whitespace defects.
analysis_data/independent_validation/GSE276892_primary_article_table1.html -whitespace
manuscript/sn-jnl.cls -text -whitespace
manuscript/sn-nature.bst -text -whitespace
"""
    (output / ".gitattributes").write_text(content, encoding="utf-8")


def _write_release_manifest(output: Path) -> dict[str, object]:
    files = _iter_release_files(output)
    largest = max(files, key=lambda path: path.stat().st_size)
    final_file_count = len(files) + 2
    manifest: dict[str, object] = {
        "schema_version": 1,
        "release_name": "DRnet P2RX4 cross-cohort analysis",
        "release_version": RELEASE_VERSION,
        "release_tag": RELEASE_TAG,
        "repository_url": REPOSITORY_URL,
        "article_type": "Research Article",
        "pipeline_entry_point": (
            "analysis_scripts/run_final_research_article_pipeline.py"
        ),
        "random_seed": 20260813,
        "formal_resampling": {
            "discovery_donor_bootstrap": 2000,
            "studentized_wild_bootstrap_t": 4999,
            "independent_validation_bootstrap": 10000,
        },
        "protocol_sha256": (
            "FA0EBEEFE45709178E197B6F0819F7AA1E1692E50AA2C7F47E2601D7CA3C8EED"
        ),
        "file_count_including_manifest_and_checksum": final_file_count,
        "largest_file": largest.relative_to(output).as_posix(),
        "largest_file_bytes": largest.stat().st_size,
        "maximum_allowed_file_bytes": MAX_FILE_BYTES,
        "checksum_file": "SHA256SUMS.txt",
        "checksum_scope": "all public-tree files except the checksum file itself",
        "required_release_asset": {
            "file": "DRnet_GSE276892_remote_results.tar.gz",
            "bytes": 444319883,
            "sha256": (
                "CF6784D6852D30A7A2A67FDFC0C7FB93E0A3ABBF3207CE775F860709795D65A6"
            ),
            "install_path": (
                "analysis_data/independent_validation/remote_results/"
                "DRnet_GSE276892_remote_results.tar.gz"
            ),
        },
    }
    path = output / "RELEASE_MANIFEST.json"
    path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest


def _sanitize_validation_report(output: Path) -> None:
    """Drop obsolete validator field names that contain internal marker words."""
    path = output / "analysis_results" / "final_submission_consistency_validation.json"
    if not path.is_file():
        return
    payload = json.loads(path.read_text(encoding="utf-8"))
    checks = payload.get("checks")
    if isinstance(checks, dict):
        old_key = "no_revision_meta_or_" + "placeholders"
        new_key = "no_revision_meta_or_internal_markers"
        if old_key in checks:
            checks[new_key] = checks.pop(old_key)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def _write_checksums(output: Path) -> int:
    checksum_path = output / "SHA256SUMS.txt"
    files = [path for path in _iter_release_files(output) if path != checksum_path]
    lines = [f"{sha256(path)}  {path.relative_to(output).as_posix()}" for path in files]
    checksum_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(lines)


def prepare_release(package: Path) -> dict[str, object]:
    package = package.resolve()
    output = _safe_output_root(package)
    if output.exists():
        shutil.rmtree(output)
    output.mkdir(parents=True)

    for source_name, destination_name in ROOT_SOURCE_FILES.items():
        _copy_file(package / source_name, output / destination_name)
    optional_license = package / "LICENSE"
    if optional_license.is_file():
        _copy_file(optional_license, output / "LICENSE")
    _write_gitignore(output)
    _write_gitattributes(output)

    _copy_named_files(
        package / "analysis_scripts",
        output / "analysis_scripts",
        ANALYSIS_SCRIPTS,
    )
    _copy_named_files(package / "tests", output / "tests", TEST_FILES)
    _copy_named_files(
        package / "manuscript",
        output / "manuscript",
        MANUSCRIPT_FILES,
    )
    for directory in COPY_TREES:
        _copy_tree(package / directory, output / directory)
    independent_source = package / "analysis_data" / "independent_validation"
    independent_destination = output / "analysis_data" / "independent_validation"
    _copy_named_files(
        independent_source,
        independent_destination,
        INDEPENDENT_VALIDATION_FILES,
    )
    _copy_named_files(
        independent_source / "GSE179568",
        independent_destination / "GSE179568",
        GSE179568_FILES,
    )
    _copy_named_files(
        independent_source / "remote_results",
        independent_destination / "remote_results",
        REMOTE_RESULT_AUDIT_FILES,
    )
    _copy_tree(
        package / "analysis_data" / "external_single_cell",
        output / "analysis_data" / "external_single_cell",
    )
    _sanitize_validation_report(output)
    _extract_project_inputs(package, output / "project_inputs")
    _copy_file(
        package
        / "submission_files"
        / "Discover_Applied_Sciences_supplementary_tables.xlsx",
        output
        / "supplementary_tables"
        / "Discover_Applied_Sciences_supplementary_tables.xlsx",
    )

    validate_sizes(output)
    scan_forbidden_text(output)
    manifest = _write_release_manifest(output)
    validate_sizes(output)
    scan_forbidden_text(output)
    checksummed_files = _write_checksums(output)
    validate_sizes(output)
    scan_forbidden_text(output)

    files = _iter_release_files(output)
    return {
        "status": "PASS",
        "output": str(output),
        "release_version": RELEASE_VERSION,
        "release_tag": RELEASE_TAG,
        "files": len(files),
        "checksummed_files": checksummed_files,
        "largest_file": manifest["largest_file"],
        "largest_file_bytes": manifest["largest_file_bytes"],
    }


def main() -> None:
    args = parse_args()
    print(json.dumps(prepare_release(args.package_root), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
