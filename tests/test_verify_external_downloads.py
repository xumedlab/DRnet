from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "analysis_scripts" / "36_verify_external_downloads.py"
SPEC = importlib.util.spec_from_file_location("verify_external_downloads", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_hash_file_updates_md5_and_sha256_in_one_result(tmp_path: Path) -> None:
    path = tmp_path / "payload.bin"
    payload = b"DRnet-download-integrity" * 1000
    path.write_bytes(payload)

    md5, sha256 = MODULE.hash_file(path, chunk_size=17)

    assert md5 == hashlib.md5(payload).hexdigest()
    assert sha256 == hashlib.sha256(payload).hexdigest()


def test_read_gencode_md5_accepts_common_md5sums_layouts(tmp_path: Path) -> None:
    path = tmp_path / "MD5SUMS"
    path.write_text(
        "0123456789abcdef0123456789abcdef  GRCh38.p13.genome.fa.gz\n"
        "fedcba9876543210fedcba9876543210 *./gencode.v42.annotation.gtf.gz\n",
        encoding="utf-8",
    )

    result = MODULE.read_gencode_md5(path)

    assert result["GRCh38.p13.genome.fa.gz"] == "0123456789abcdef0123456789abcdef"
    assert (
        result["gencode.v42.annotation.gtf.gz"]
        == "fedcba9876543210fedcba9876543210"
    )
