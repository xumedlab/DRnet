from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = (
    Path(__file__).parents[1]
    / "analysis_scripts"
    / "38_raw_count_p2rx4_validation.py"
)
SPEC = importlib.util.spec_from_file_location("raw_count_p2rx4_validation", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_verify_output_manifest_detects_exact_sha256(tmp_path: Path) -> None:
    payload = tmp_path / "results" / "table.tsv"
    payload.parent.mkdir()
    payload.write_text("gene\tcount\nP2RX4\t10\n", encoding="utf-8")
    digest = hashlib.sha256(payload.read_bytes()).hexdigest()
    (tmp_path / "output_sha256.txt").write_text(
        f"{digest}  results/table.tsv\n", encoding="utf-8"
    )

    result = MODULE.verify_output_manifest(tmp_path)

    assert result["manifest_entries"] == 1
    assert result["checked_files"] == 1
    assert result["failure_count"] == 0


def test_json_safe_replaces_nonfinite_numpy_values() -> None:
    result = MODULE.json_safe(
        {"nan": np.float64(np.nan), "inf": float("inf"), "ok": np.int64(3)}
    )

    assert result == {"nan": None, "inf": None, "ok": 3}


def test_one_sided_greater_p_uses_upper_normal_tail() -> None:
    assert np.isclose(MODULE.one_sided_greater_p(0.0), 0.5)
    assert np.isclose(MODULE.one_sided_greater_p(1.96), 0.0249978951)
    assert MODULE.one_sided_greater_p(-1.0) > 0.5


def test_leave_one_out_preserves_biological_sample_grain() -> None:
    samples = [f"PDR_S{i}" for i in range(1, 9)] + [
        f"CTRL_S{i}" for i in range(1, 10)
    ]
    table = pd.DataFrame(
        {
            "diagnosis": ["PDR"] * 8 + ["control"] * 9,
            "p2rx4_deseq2_normalized_count": np.r_[
                np.arange(20, 28), np.arange(1, 10)
            ],
        },
        index=samples,
    )

    result = MODULE.leave_one_out(table)

    assert len(result) == 18
    assert result.iloc[0]["omitted_sample"] == "none"
    assert result.iloc[0]["n_pdr"] == 8
    assert result.iloc[0]["n_control"] == 9
    assert result["direction_positive"].all()


def test_reference_sensitivity_reports_one_read_difference() -> None:
    samples = ["PDR_S5", "MP_S13"]
    genes = [MODULE.TARGET_GENE_ID, "ENSG_OTHER.1"]
    all_counts = pd.DataFrame([[10, 100], [5, 80]], index=samples, columns=genes)
    primary_counts = pd.DataFrame([[11, 100], [5, 79]], index=samples, columns=genes)
    all_table = pd.DataFrame(
        {
            "diagnosis": ["PDR", "control"],
            "p2rx4_deseq2_normalized_count": [9.5, 5.2],
        },
        index=samples,
    )
    primary_table = pd.DataFrame(
        {"p2rx4_deseq2_normalized_count": [10.2, 5.1]}, index=samples
    )

    result = MODULE.reference_sensitivity(
        all_counts, primary_counts, all_table, primary_table
    )

    assert result.loc[0, "p2rx4_count_difference_primary_minus_all"] == 1
    assert result.loc[1, "p2rx4_count_difference_primary_minus_all"] == 0
    assert result["common_expressed_gene_ids"].eq(2).all()
