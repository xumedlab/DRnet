"""Focused tests for author-mapped GSE130636 localisation."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "analysis_scripts"
    / "25_voigt_single_cell_localization.py"
)
SPEC = importlib.util.spec_from_file_location("voigt_localization", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_author_mapping_has_complete_original_labels() -> None:
    mapping = pd.read_csv(
        Path(__file__).resolve().parents[1]
        / "analysis_data"
        / "external_single_cell"
        / "voigt2019_author_cluster_mapping.csv",
        dtype=str,
    )
    MODULE.validate_author_mapping(mapping)
    assert mapping.loc[mapping["cluster_label"] == "7", "author_label"].item() == "Retinal ganglion cells"
    assert mapping.loc[mapping["cluster_label"] == "8A", "author_label"].item() == "Horizontal cells"
    assert mapping.loc[mapping["cluster_label"] == "8B", "author_label"].item() == "Amacrine cells"
    assert mapping.loc[mapping["cluster_label"] == "9", "author_label"].item() == "Unknown"


def test_library_filename_parser() -> None:
    result = MODULE.parse_library_name(
        Path("GSM3745992_fovea_donor_1_expression.tsv.gz")
    )
    assert result == ("GSM3745992", "fovea", "donor_1")


def test_localization_panel_uses_union_of_both_estimands() -> None:
    summary = {
        "total_top5": ["P2RX4", "TLR2", "CD82", "NLRP3", "FPR1"],
        "dme_conditioned_top5": ["P2RX4", "SLC31A1", "NLRP3", "CD82", "FPR1"]
    }
    assert MODULE.localization_panel_from_discovery_summary(summary) == [
        "P2RX4",
        "TLR2",
        "CD82",
        "NLRP3",
        "FPR1",
        "SLC31A1",
    ]


def test_localization_panel_rejects_retired_schema() -> None:
    with np.testing.assert_raises(ValueError):
        MODULE.localization_panel_from_discovery_summary(
            {
                "candidate_genes": ["P2RX4", "SLC31A1", "NLRP3", "CD82", "FPR1"],
                "dme_conditioned_top5": ["P2RX4", "SLC31A1", "NLRP3", "CD82", "FPR1"],
            }
        )


def test_linearized_expression_uses_max_index_not_compositional_share() -> None:
    pseudobulk = pd.DataFrame(
        {
            "cell_type": ["Microglia", "Pericytes", "Microglia", "Pericytes"],
            "gene_symbol": ["A", "A", "B", "B"],
            "n_cells": [10, 10, 10, 10],
            "mean_log_normalized_expression": [3.0, 1.0, 0.0, 2.0],
            "mean_expm1_log_normalized_expression": [
                np.expm1(3.0),
                np.expm1(1.0),
                np.expm1(0.0),
                np.expm1(2.0),
            ],
            "detection_fraction": [0.5, 0.2, 0.0, 0.4],
        }
    )
    result = MODULE.aggregate_cell_types(pseudobulk)
    assert "relative_expression_share" not in result
    maxima = result.groupby("gene_symbol")["expression_relative_to_gene_max"].max()
    np.testing.assert_allclose(maxima.to_numpy(), np.ones(len(maxima)))
    assert result["tau_cell_type_specificity"].between(0, 1).all()


def test_paired_region_summary_uses_donor_pairs() -> None:
    pseudobulk = pd.DataFrame(
        {
            "cell_type": ["Microglia"] * 6,
            "gene_symbol": ["A"] * 6,
            "donor": ["d1", "d1", "d2", "d2", "d3", "d3"],
            "region": ["fovea", "peripheral"] * 3,
            "mean_log_normalized_expression": [2.0, 1.0, 4.0, 2.0, 1.0, 1.0],
            "mean_expm1_log_normalized_expression": [2.0, 1.0, 4.0, 2.0, 1.0, 1.0],
        }
    )
    result = MODULE.paired_region_summary(pseudobulk).iloc[0]
    assert result["n_paired_donors"] == 3
    assert result["mean_paired_fovea_minus_peripheral"] == 1.0
