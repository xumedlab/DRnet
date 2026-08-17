from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "analysis_scripts" / "34_validate_final_submission.py"
SPEC = importlib.util.spec_from_file_location("validate_final_submission", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_abstract_word_count_handles_tex_commands() -> None:
    tex = r"\begin{abstract}Alpha \textit{beta gamma} $P=0.05$.\end{abstract}"
    assert MODULE.abstract_word_count(tex) < 10


def test_protocol_hash_is_frozen() -> None:
    assert MODULE.sha256(MODULE.PROTOCOL) == MODULE.EXPECTED_PROTOCOL_SHA256


def test_layout_resolver_uses_first_existing_file(tmp_path: Path) -> None:
    missing = tmp_path / "missing.txt"
    existing = tmp_path / "existing.txt"
    existing.write_text("ok", encoding="utf-8")
    assert MODULE.first_existing(missing, existing) == existing


def test_final_submission_consistency() -> None:
    assert MODULE.validate()["status"] == "PASS"
