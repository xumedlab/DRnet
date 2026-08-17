#!/usr/bin/env python3
"""Build deterministic, internally verified submission and reproducibility ZIPs."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import zipfile
from pathlib import Path, PurePosixPath


FIXED_ZIP_TIME = (2026, 8, 14, 0, 0, 0)
LATEX_FILES = (
    "discover_applied_sciences_submission.tex",
    "discover_applied_sciences_supplementary_information.tex",
    "discover_references.bib",
    "sn-jnl.cls",
    "sn-nature.bst",
)
PDF_FILES = (
    "discover_applied_sciences_submission.pdf",
    "discover_applied_sciences_supplementary_information.pdf",
    "Discover_Applied_Sciences_cover_letter.pdf",
)
WORKBOOK = "Discover_Applied_Sciences_supplementary_tables.xlsx"
LATEX_ARCHIVE = "Discover_Applied_Sciences_LaTeX_Source.zip"
REPRO_ARCHIVE = "Discover_Applied_Sciences_reproducibility_files.zip"
FROZEN_ARCHIVE = "DRnet_P2RX4_frozen_release.zip"
PAYLOAD_SUMS = "DRnet_P2RX4_frozen_payload_SHA256SUMS.txt"
REMOTE_ARCHIVE = (
    "analysis_data/independent_validation/remote_results/"
    "DRnet_GSE276892_remote_results.tar.gz"
)


class ArchiveBuildError(RuntimeError):
    """Raised when an expected release input or archive check fails."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def _validate_member_name(name: PurePosixPath) -> None:
    if name.is_absolute() or not name.parts or any(part in {"", ".", ".."} for part in name.parts):
        raise ArchiveBuildError(f"Unsafe ZIP member path: {name}")


def _write_zip(destination: Path, members: list[tuple[Path, PurePosixPath]]) -> dict[str, object]:
    normalized: dict[str, Path] = {}
    for source, member in members:
        _validate_member_name(member)
        if not source.is_file():
            raise ArchiveBuildError(f"Missing archive input: {source}")
        key = member.as_posix()
        if key in normalized:
            raise ArchiveBuildError(f"Duplicate ZIP member: {key}")
        normalized[key] = source

    destination.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(
        destination,
        mode="w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=6,
        allowZip64=True,
    ) as archive:
        for name in sorted(normalized):
            source = normalized[name]
            info = zipfile.ZipInfo(name, date_time=FIXED_ZIP_TIME)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = (0o100644 & 0xFFFF) << 16
            with source.open("rb") as input_stream, archive.open(info, "w") as output_stream:
                shutil.copyfileobj(input_stream, output_stream, length=1024 * 1024)

    if destination.read_bytes()[:4] != b"PK\x03\x04":
        raise ArchiveBuildError(f"Invalid ZIP header: {destination}")
    with zipfile.ZipFile(destination) as archive:
        corrupt = archive.testzip()
        if corrupt is not None:
            raise ArchiveBuildError(f"ZIP CRC failure in {destination}: {corrupt}")
        names = archive.namelist()
        if names != sorted(names) or len(names) != len(set(names)):
            raise ArchiveBuildError(f"ZIP member ordering/uniqueness failure: {destination}")
    return {
        "file": destination.name,
        "bytes": destination.stat().st_size,
        "sha256": sha256(destination),
        "members": len(normalized),
    }


def _tree_members(root: Path, prefix: PurePosixPath | None = None) -> list[tuple[Path, PurePosixPath]]:
    if not root.is_dir():
        raise ArchiveBuildError(f"Missing archive tree: {root}")
    members: list[tuple[Path, PurePosixPath]] = []
    for path in sorted(candidate for candidate in root.rglob("*") if candidate.is_file()):
        relative = PurePosixPath(path.relative_to(root).as_posix())
        member = relative if prefix is None else prefix / relative
        members.append((path, member))
    return members


def _write_checksum_file(destination: Path, entries: list[tuple[str, Path]]) -> None:
    lines = []
    for name, path in sorted(entries):
        if not path.is_file():
            raise ArchiveBuildError(f"Missing checksum input: {path}")
        lines.append(f"{sha256(path)}  {name}")
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def build(package: Path) -> dict[str, object]:
    package = package.resolve()
    manuscript = package / "manuscript"
    figures = package / "figures"
    public_tree = package / "github_upload"
    submission = package / "submission_files"
    submission.mkdir(parents=True, exist_ok=True)

    latex_members = [
        (manuscript / name, PurePosixPath("manuscript") / name) for name in LATEX_FILES
    ]
    latex_members.extend(_tree_members(figures, PurePosixPath("figures")))
    latex_result = _write_zip(submission / LATEX_ARCHIVE, latex_members)

    repro_result = _write_zip(submission / REPRO_ARCHIVE, _tree_members(public_tree))

    manifest_source = package / "RELEASE_MANIFEST.json"
    manifest_destination = submission / "RELEASE_MANIFEST.json"
    shutil.copy2(manifest_source, manifest_destination)

    payload_entries = [
        *[(name, submission / name) for name in PDF_FILES],
        (WORKBOOK, submission / WORKBOOK),
        (LATEX_ARCHIVE, submission / LATEX_ARCHIVE),
        (REPRO_ARCHIVE, submission / REPRO_ARCHIVE),
        ("RELEASE_MANIFEST.json", manifest_destination),
    ]
    payload_sums = submission / PAYLOAD_SUMS
    _write_checksum_file(payload_sums, payload_entries)

    frozen_members = [
        (path, PurePosixPath(name)) for name, path in payload_entries
    ]
    frozen_members.extend(
        [
            (payload_sums, PurePosixPath(PAYLOAD_SUMS)),
            (package / "README_submission_package.md", PurePosixPath("README_submission_package.md")),
            (package / "REPRODUCIBILITY_README.md", PurePosixPath("REPRODUCIBILITY_README.md")),
            (package / "GITHUB_UPLOAD_GUIDE.md", PurePosixPath("GITHUB_UPLOAD_GUIDE.md")),
        ]
    )
    frozen_result = _write_zip(submission / FROZEN_ARCHIVE, frozen_members)

    remote_archive = package / REMOTE_ARCHIVE
    final_entries = [
        (remote_archive.name, remote_archive),
        *payload_entries,
        (PAYLOAD_SUMS, payload_sums),
        (FROZEN_ARCHIVE, submission / FROZEN_ARCHIVE),
    ]
    final_sums = submission / "SHA256SUMS.txt"
    _write_checksum_file(final_sums, final_entries)
    shutil.copy2(final_sums, package / "SHA256SUMS.txt")

    return {
        "status": "PASS",
        "latex_source": latex_result,
        "reproducibility": repro_result,
        "frozen_release": frozen_result,
        "submission_checksums": str(final_sums),
        "remote_release_asset": {
            "file": remote_archive.name,
            "bytes": remote_archive.stat().st_size,
            "sha256": sha256(remote_archive),
        },
    }


def parse_args() -> argparse.Namespace:
    default_package = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path, default=default_package)
    return parser.parse_args()


def main() -> None:
    print(json.dumps(build(parse_args().package_root), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
