from __future__ import annotations

import hashlib
import importlib.util
import zipfile
from pathlib import Path, PurePosixPath


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "analysis_scripts"
    / "39_build_submission_archives.py"
)
SPEC = importlib.util.spec_from_file_location("build_submission_archives", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_deterministic_zip_is_sorted_unique_and_crc_clean(tmp_path: Path) -> None:
    first = tmp_path / "first.txt"
    second = tmp_path / "second.txt"
    first.write_text("alpha\n", encoding="utf-8")
    second.write_text("beta\n", encoding="utf-8")
    destination = tmp_path / "release.zip"
    members = [
        (second, PurePosixPath("z/second.txt")),
        (first, PurePosixPath("a/first.txt")),
    ]

    result_one = MODULE._write_zip(destination, members)
    digest_one = _digest(destination)
    result_two = MODULE._write_zip(destination, list(reversed(members)))

    assert destination.read_bytes()[:4] == b"PK\x03\x04"
    assert digest_one == _digest(destination)
    assert result_one == result_two
    with zipfile.ZipFile(destination) as archive:
        assert archive.namelist() == ["a/first.txt", "z/second.txt"]
        assert archive.testzip() is None


def test_unsafe_or_duplicate_member_is_rejected(tmp_path: Path) -> None:
    source = tmp_path / "input.txt"
    source.write_text("x", encoding="utf-8")
    destination = tmp_path / "release.zip"

    try:
        MODULE._write_zip(destination, [(source, PurePosixPath("../input.txt"))])
    except MODULE.ArchiveBuildError:
        pass
    else:
        raise AssertionError("unsafe member was not rejected")

    duplicate = [(source, PurePosixPath("same.txt"))] * 2
    try:
        MODULE._write_zip(destination, duplicate)
    except MODULE.ArchiveBuildError:
        pass
    else:
        raise AssertionError("duplicate member was not rejected")
