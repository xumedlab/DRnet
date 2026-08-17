from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "analysis_scripts"
    / "run_final_research_article_pipeline.py"
)
SPEC = importlib.util.spec_from_file_location("final_pipeline", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_submission_validator_is_optional(tmp_path: Path, monkeypatch) -> None:
    scripts = tmp_path / "analysis_scripts"
    scripts.mkdir()
    calls: list[tuple[list[str], Path]] = []
    monkeypatch.setattr(MODULE, "run", lambda command, workdir: calls.append((command, workdir)))

    MODULE.run_submission_validator_if_present(scripts, "python", tmp_path)
    assert calls == []

    validator = scripts / "34_validate_final_submission.py"
    validator.write_text("print('validate')\n", encoding="utf-8")
    MODULE.run_submission_validator_if_present(scripts, "python", tmp_path)

    assert calls == [(["python", str(validator)], tmp_path)]
