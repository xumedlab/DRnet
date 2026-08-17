"""Regression tests for independent P2RX4 validation safeguards."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "analysis_scripts"
    / "32_independent_p2rx4_validation.py"
)
SPEC = importlib.util.spec_from_file_location("independent_p2rx4_validation", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_gse276892_source_and_qc_metadata_are_explicit() -> None:
    data_dir = SCRIPT.parents[1] / "analysis_data" / "independent_validation"
    metadata, _, _ = MODULE.load_gse276892(data_dir)

    assert metadata["reused_control"].sum() == 7
    assert (
        metadata.loc[metadata["disease_group"].eq("control"), "reused_control"].eq(0).sum()
        == 2
    )
    assert metadata["cell_count"].notna().all()
    assert metadata["rna_concentration_pg_per_ul"].notna().all()


def test_gse276892_source_adjusted_results_and_qc_correlations() -> None:
    data_dir = SCRIPT.parents[1] / "analysis_data" / "independent_validation"
    metadata, _, _ = MODULE.load_gse276892(data_dir)
    adjusted = MODULE.adjusted_gse276892(metadata)
    qc = MODULE.gse276892_qc_correlations(metadata)

    assert "source_adjusted_ols_hc3_t_beta" in adjusted
    assert "source_age_sex_adjusted_ols_hc3_t_p_two_sided" in adjusted
    assert np.isfinite(list(adjusted.values())).all()
    assert len(qc) == 4
    assert qc["n_profiles"].eq(17).all()
    assert qc["spearman_p_two_sided"].between(0, 1).all()


def test_gse179568_clinical_mapping_matches_supplementary_table() -> None:
    data_dir = SCRIPT.parents[1] / "analysis_data" / "independent_validation"
    metadata, _, _ = MODULE.load_gse179568(data_dir)

    assert len(metadata) == 24
    assert metadata["age"].notna().all()
    assert metadata.loc[metadata["sample_id"].eq("PDR_S1"), "age"].item() == 40
    assert metadata.loc[metadata["sample_id"].eq("Gliose_S1"), "age"].item() == 71
    assert metadata.loc[metadata["sample_id"].eq("ILM_S7"), "dm_type"].item() == "II"
    assert metadata.loc[metadata["disease_group"].eq("PDR"), "previous_prp"].sum() == 5
    assert (
        metadata.loc[
            metadata["disease_group"].eq("PDR"),
            "previous_anti_vegf_over_3_months",
        ].sum()
        == 2
    )


def test_gse179568_clinical_sensitivity_reports_no_age_common_support() -> None:
    data_dir = SCRIPT.parents[1] / "analysis_data" / "independent_validation"
    metadata, _, _ = MODULE.load_gse179568(data_dir)
    results, balance = MODULE.gse179568_clinical_sensitivity(metadata)

    assert len(results) == 3
    assert not results["age_common_support"].any()
    assert results["age_overlap_interval"].eq("none").all()
    assert results["age_pdr_range"].eq("20-60").all()
    assert {
        "unadjusted_ols_hc3_t_beta",
        "age_sex_adjusted_ols_hc3_t_beta",
        "age_sex_treatment_adjusted_ols_hc3_t_beta",
        "age_sex_adjusted_huber_h1_beta",
        "mann_whitney_p_two_sided",
    }.issubset(results.columns)
    numeric = results.filter(regex="(_beta|_p_two_sided|_asymptotic)$")
    assert np.isfinite(numeric.to_numpy(dtype=float)).all()
    assert len(balance) == 21
    assert balance["comparison"].nunique() == 3


def test_post_protocol_holm_uses_three_two_sided_tests() -> None:
    table = pd.DataFrame(
        {
            "dataset": ["GSE276892", "GSE179568", "GSE179568", "GSE94019"],
            "comparison": [
                "PDR vs all surgical controls",
                "RNV vs macular-pucker membrane",
                "RNV vs macular-hole ILM",
                "PDR FVM CD31+ vs control retinal CD31+",
            ],
            "mann_whitney_p_two_sided": [0.09, 0.01, 0.20, 0.30],
        }
    )
    adjusted = MODULE.add_post_protocol_holm(table)
    family = adjusted["post_protocol_holm_family"]

    assert family.sum() == 3
    assert adjusted.loc[~family, "mann_whitney_p_two_sided_holm_post_protocol_3"].isna().all()
    assert np.allclose(
        adjusted.loc[family, "mann_whitney_p_two_sided_holm_post_protocol_3"],
        [0.03, 0.40, 0.40],
    )


def test_screening_table_records_included_and_excluded_datasets() -> None:
    table = MODULE.dataset_eligibility().set_index("dataset")

    assert {"GSE276892", "GSE179568", "GSE94019", "GSE102485", "GSE130636"}.issubset(
        table.index
    )
    assert table["search_date"].eq(MODULE.SEARCH_DATE).all()
    assert table.loc["GSE102485", "screening_decision"].startswith("excluded")
    assert "separate GEO" in table.loc["GSE179568", "analysis_status"]
    assert table["search_log_status"].str.contains("retrospectively").all()
    assert "unknown" in table.loc["GSE179568", "p2rx4_value_inspection_status"]


def test_patient_overlap_audit_does_not_overclaim_independence() -> None:
    audit = MODULE.patient_overlap_audit().iloc[0]

    assert bool(audit["shared_institution_or_team"])
    assert not bool(audit["author_confirmation_obtained"])
    assert "cannot be established" in audit["public_metadata_conclusion"]
    assert "separate GEO datasets" in audit["manuscript_language"]
