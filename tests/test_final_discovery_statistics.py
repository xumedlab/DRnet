"""Regression tests for the final Research Article discovery analysis."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t


SCRIPT = Path(__file__).resolve().parents[1] / "analysis_scripts" / "33_final_discovery_statistics.py"
SPEC = importlib.util.spec_from_file_location("final_discovery_statistics", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def synthetic_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    index = pd.Index([f"d{i}" for i in range(12)], name="donor")
    severity = np.arange(12, dtype=float)
    metadata = pd.DataFrame(
        {
            "detailed_group": np.repeat(MODULE.GROUP_ORDER, 3),
            "severity_mean": severity,
            "severity_worst_eye": severity,
            "severity_rank": pd.Series(severity).rank().to_numpy(),
            "stage_ordinal": np.repeat([0, 1, 1, 2], 3),
            "dme": np.repeat([0, 0, 1, 1], 3),
            "age": np.linspace(50, 80, 12),
            "sex_male": np.tile([0, 1], 6),
            "pmi": np.linspace(100, 300, 12),
            "rin": np.linspace(5, 9, 12),
        },
        index=index,
    )
    expression = pd.DataFrame(
        {
            "POS": severity + np.sin(severity),
            "NEG": -severity + np.cos(severity),
            "NOISE": np.sin(severity * 2.1),
        },
        index=index,
    )
    return expression, metadata


def test_hc3_uses_residual_df_t_inference() -> None:
    expression, metadata = synthetic_inputs()
    table = MODULE.fit_gene_models(
        expression,
        metadata,
        list(expression.columns),
        model_name="test",
        include_dme=False,
    )
    design = MODULE.design_matrix(metadata).to_numpy(dtype=float)
    expected_df = len(metadata) - np.linalg.matrix_rank(design)
    assert table["residual_df"].eq(expected_df).all()
    positive = table.set_index("gene_symbol").loc["POS"]
    expected_p = 2 * student_t.sf(abs(positive["severity_t"]), expected_df)
    assert np.isclose(positive["severity_pvalue_t"], expected_p)
    assert table.iloc[0]["gene_symbol"] == "POS"
    assert table.iloc[-1]["gene_symbol"] == "NEG"


def test_bootstrap_reports_rank_and_topk_uncertainty() -> None:
    expression, metadata = synthetic_inputs()
    table = MODULE.bootstrap_stability(
        expression,
        metadata,
        list(expression.columns),
        iterations=20,
        seed=7,
    )
    required = {
        "bootstrap_top5_frequency",
        "bootstrap_top10_frequency",
        "bootstrap_top20_frequency",
        "bootstrap_rank_median",
        "bootstrap_rank_ci_low",
        "bootstrap_rank_ci_high",
        "bootstrap_beta_ci_low",
        "bootstrap_beta_ci_high",
    }
    assert required.issubset(table.columns)
    assert table["bootstrap_top5_frequency"].between(0, 1).all()


def test_contrasts_include_within_and_global_multiplicity() -> None:
    expression, metadata = synthetic_inputs()
    table = MODULE.clinical_contrasts(expression, metadata, list(expression.columns))
    assert len(table) == 5 * len(expression.columns)
    assert "mannwhitney_padj_within_contrast_158" in table
    assert "mannwhitney_padj_global_790" in table
    assert set(table["multiplicity_family"]) == {
        "primary_contrast_family",
        "exploratory_contrast_family",
    }


def test_wild_bootstrap_is_reproducible_and_bounded() -> None:
    expression, metadata = synthetic_inputs()
    first = MODULE.wild_bootstrap_t(
        expression,
        metadata,
        list(expression.columns),
        include_dme=False,
        iterations=19,
        seed=11,
    )
    second = MODULE.wild_bootstrap_t(
        expression,
        metadata,
        list(expression.columns),
        include_dme=False,
        iterations=19,
        seed=11,
    )
    pd.testing.assert_frame_equal(first, second)
    assert first["wild_bootstrap_t_two_sided_p"].between(1 / 20, 1).all()
