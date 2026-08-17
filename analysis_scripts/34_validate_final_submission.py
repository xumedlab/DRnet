#!/usr/bin/env python3
"""Validate final Research Article data, claims, protocol, figures, and workbook."""

from __future__ import annotations

import hashlib
import json
import math
import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd


PACKAGE = Path(__file__).resolve().parents[1]
RESULTS = PACKAGE / "analysis_results"
MANUSCRIPT = PACKAGE / "manuscript" / "discover_applied_sciences_submission.tex"
SUPPLEMENT = PACKAGE / "manuscript" / "discover_applied_sciences_supplementary_information.tex"
COVER_LETTER = PACKAGE / "manuscript" / "discover_applied_sciences_cover_letter.tex"


def first_existing(*candidates: Path) -> Path:
    """Resolve a file across the full submission and public-release layouts."""
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return candidates[0]


README = first_existing(PACKAGE / "README_submission_package.md", PACKAGE / "README.md")
REPRODUCIBILITY = first_existing(
    PACKAGE / "REPRODUCIBILITY_README.md",
    PACKAGE / "REPRODUCIBILITY.md",
)
WORKBOOK = first_existing(
    PACKAGE / "submission_files" / "Discover_Applied_Sciences_supplementary_tables.xlsx",
    PACKAGE / "supplementary_tables" / "Discover_Applied_Sciences_supplementary_tables.xlsx",
)
PROTOCOL = PACKAGE / "analysis_data" / "independent_validation" / "P2RX4_VALIDATION_PROTOCOL_LOCALLY_FROZEN.md"
REPORT = RESULTS / "final_submission_consistency_validation.json"
EXPECTED_PROTOCOL_SHA256 = "FA0EBEEFE45709178E197B6F0819F7AA1E1692E50AA2C7F47E2601D7CA3C8EED"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def abstract_word_count(tex: str) -> int:
    match = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", tex, re.S)
    if not match:
        return 10_000
    text = re.sub(r"\\[A-Za-z]+(?:\[[^]]*\])?\{([^}]*)\}", r"\1", match.group(1))
    text = re.sub(r"\\[A-Za-z]+|[{}$~]", " ", text)
    return len(re.findall(r"\b[\w'-]+\b", text))


def workbook_sheet_names(path: Path) -> list[str]:
    with zipfile.ZipFile(path) as archive:
        workbook_xml = archive.read("xl/workbook.xml")
    root = ET.fromstring(workbook_xml)
    namespace = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    return [sheet.attrib["name"] for sheet in root.findall("x:sheets/x:sheet", namespace)]


def isclose(value: object, expected: float, tolerance: float = 1e-10) -> bool:
    return math.isclose(float(value), expected, rel_tol=tolerance, abs_tol=tolerance)


def validate() -> dict[str, object]:
    tex = MANUSCRIPT.read_text(encoding="utf-8")
    supplement = SUPPLEMENT.read_text(encoding="utf-8")
    cover_letter = COVER_LETTER.read_text(encoding="utf-8") if COVER_LETTER.is_file() else ""
    public_text = "\n".join(
        [
            tex,
            supplement,
            cover_letter,
            README.read_text(encoding="utf-8"),
            REPRODUCIBILITY.read_text(encoding="utf-8"),
        ]
    )
    discovery = json.loads((RESULTS / "final_discovery_summary.json").read_text(encoding="utf-8"))
    primary = pd.read_csv(RESULTS / "final_primary_ranking_and_stability.csv")
    severity_sensitivity = pd.read_csv(RESULTS / "final_severity_model_sensitivities.csv")
    external = pd.read_csv(RESULTS / "Independent_validation_P2RX4_results.csv")
    raw_models = pd.read_csv(RESULTS / "raw_count_p2rx4_deseq2_results.csv")
    raw_qc = pd.read_csv(RESULTS / "raw_count_p2rx4_qc_correlations.csv")
    raw_reference = pd.read_csv(RESULTS / "raw_count_reference_sensitivity_by_sample.csv")
    raw_integrity = json.loads(
        (RESULTS / "raw_count_reconstruction_integrity.json").read_text(encoding="utf-8")
    )
    gse179_clinical = pd.read_csv(
        RESULTS / "Independent_validation_GSE179568_clinical_sensitivity.csv"
    )
    overlap_audit = pd.read_csv(RESULTS / "Independent_validation_patient_overlap_audit.csv")
    eligibility = pd.read_csv(RESULTS / "Independent_validation_dataset_eligibility.csv")
    lodo_normal = pd.read_csv(RESULTS / "normal_retina_leave_one_donor_dominance.csv")
    library_normal = pd.read_csv(RESULTS / "normal_retina_library_dominance.csv")
    author_mapping = pd.read_csv(RESULTS / "normal_retina_author_cluster_mapping.csv", dtype=str)
    p2rx4 = primary.loc[primary["gene_symbol"].eq("P2RX4")].iloc[0]
    gse179 = external.loc[
        external["dataset"].eq("GSE179568")
        & external["comparison"].eq("RNV vs macular-pucker membrane")
    ].iloc[0]
    gse179_main = gse179_clinical.loc[
        gse179_clinical["comparison"].eq("RNV vs macular-pucker membrane")
    ].iloc[0]
    rank_p2rx4 = severity_sensitivity.loc[
        severity_sensitivity["model"].eq("total_rank_transformed_severity")
        & severity_sensitivity["gene_symbol"].eq("P2RX4")
    ].iloc[0]
    raw_disease = raw_models.loc[
        raw_models["workflow"].eq("all_regions_source_reconstruction")
        & raw_models["model"].eq("disease_only")
    ].iloc[0]
    raw_primary = raw_models.loc[
        raw_models["workflow"].eq("primary_assembly_sensitivity")
        & raw_models["model"].eq("disease_only")
    ].iloc[0]
    raw_source = raw_models.loc[
        raw_models["workflow"].eq("all_regions_source_reconstruction")
        & raw_models["model"].eq("source_adjusted")
    ].iloc[0]
    raw_without_s10 = raw_models.loc[
        raw_models["workflow"].eq("all_regions_source_reconstruction")
        & raw_models["model"].eq("disease_only_without_PDR_S10")
    ].iloc[0]
    forbidden = re.compile(
        r"reviewer-driven|earlier five|retired five|supersede incompatible|"
        r"AUTHOR TO INSERT|TODO|PLACEHOLDER|before submission|submission-ready|"
        r"must not be uploaded|AUTHOR_ACTION_REQUIRED|matched-membrane",
        re.I,
    )
    expected_figures = [
        "Figure_1_final_study_design.pdf",
        "Figure_2_total_and_dme_conditioned_associations.pdf",
        "Figure_3_candidate_stability_and_uncertainty.pdf",
        "Figure_4_independent_P2RX4_validation.pdf",
        "Figure_5_cell_type_localization.pdf",
        "Supplementary_Figure_raw_count_reconstruction.pdf",
        "Supplementary_Figure_GSE179568_clinical_sensitivity.pdf",
    ]
    sheets = workbook_sheet_names(WORKBOOK)
    checks = {
        "research_article_title_and_scope": "Research Article" in supplement
        and "Brief Report" not in tex
        and "Donor-aware cross-cohort transcriptomics prioritizes P2RX4" in tex
        and "Target-locked" not in re.search(r"\\title\{.*?\}", tex, re.S).group(0),
        "abstract_under_250_words": abstract_word_count(tex) < 250,
        "discovery_donors_26": discovery["n_diabetic_donors"] == 26,
        "inflammatory_genes_158": discovery["n_inflammatory_genes"] == 158,
        "formal_resampling_counts": discovery["donor_bootstrap_iterations"] == 2000
        and discovery["wild_bootstrap_iterations"] == 4999,
        "p2rx4_primary_numbers": isclose(p2rx4["severity_beta"], 0.6879902439039749)
        and isclose(p2rx4["severity_pvalue_t"], 0.0014857541720058982)
        and isclose(p2rx4["bootstrap_top5_frequency"], 0.6165)
        and isclose(p2rx4["lodo_top5_frequency"], 1.0),
        "no_discovery_fdr_hits": discovery["fdr_significant_total_model"] == 0,
        "rank_sensitivity_fdr_is_reported": isclose(
            rank_p2rx4["severity_padj_bh_158"], 0.0458182, tolerance=1e-6
        )
        and "only tested severity representation" in tex,
        "gse179_unadjusted_numbers": int(gse179["n_case"]) == 7
        and int(gse179["n_control"]) == 10
        and isclose(gse179["mann_whitney_p_two_sided"], 0.009666803784450843)
        and isclose(gse179["cliff_delta"], 0.7428571428571429)
        and isclose(gse179_main["unadjusted_ols_hc3_t_p_two_sided"], 0.03145654660356175),
        "gse179_clinical_confounding_is_explicit": not bool(gse179_main["age_common_support"])
        and isclose(gse179_main["age_standardized_mean_difference"], -2.8265895085100694)
        and isclose(gse179_main["age_sex_adjusted_ols_hc3_t_p_two_sided"], 0.7726181437300494)
        and "cannot isolate disease from age or tissue composition" in str(
            gse179_main["identifiability_statement"]
        ),
        "raw_count_negative_binomial_results": int(raw_disease["n_samples"]) == 17
        and isclose(raw_disease["log2_fold_change_pdr_vs_control"], 0.29457265461702686)
        and isclose(raw_disease["wald_p_two_sided"], 0.6152642168527742)
        and isclose(raw_primary["log2_fold_change_pdr_vs_control"], 0.28944674336684617)
        and isclose(raw_source["wald_p_two_sided"], 0.40737182447961895)
        and isclose(raw_without_s10["wald_p_two_sided"], 0.10508246320350323),
        "raw_count_qc_and_reference_sensitivity": len(raw_qc) == 18
        and raw_qc["spearman_p_bh_18_tests"].notna().all()
        and raw_reference["p2rx4_count_difference_primary_minus_all"].abs().max() <= 1
        and raw_reference["spearman_log1p_gene_counts"].min() > 0.98,
        "remote_archive_and_manifests_verified": raw_integrity["archive_sha256"]
        == "CF6784D6852D30A7A2A67FDFC0C7FB93E0A3ABBF3207CE775F860709795D65A6"
        and raw_integrity["gzip_header"] == "1F8B08"
        and raw_integrity["unsafe_member_count"] == 0
        and all(
            workflow["failure_count"] == 0 and workflow["checked_files"] == 221
            for workflow in raw_integrity["workflow_checks"].values()
        ),
        "protocol_hash_unchanged": sha256(PROTOCOL) == EXPECTED_PROTOCOL_SHA256,
        "normal_retina_six_libraries": library_normal["geo_accession"].nunique() == 6,
        "normal_retina_candidate_union_and_lodo": set(lodo_normal["gene_symbol"])
        == set(discovery["candidate_union"])
        == {"P2RX4", "TLR2", "CD82", "NLRP3", "FPR1", "SLC31A1"}
        and lodo_normal.groupby("gene_symbol").size().eq(3).all()
        and int(
            lodo_normal.loc[
                lodo_normal["gene_symbol"].eq("TLR2"), "matches_full_cohort_dominance"
            ].sum()
        )
        == 3,
        "voigt_author_mapping_complete": (
            author_mapping["cluster_label"].nunique() == 18
            and not author_mapping["cluster_label"].duplicated().any()
            and author_mapping["source_evidence"].str.contains("Voigt et al. 2019").all()
        ),
        "workbook_s1_to_s28": sheets == ["README", *[f"S{i}" for i in range(1, 29)]],
        "final_figures_exist": all((PACKAGE / "figures" / name).is_file() for name in expected_figures),
        "manuscript_uses_figure_4_filename": "Figure_4_independent_P2RX4_validation" in tex,
        "supplementary_figures_prefixed": r"\renewcommand{\thefigure}{S\arabic{figure}}" in supplement,
        "no_revision_meta_or_internal_markers": not forbidden.search(public_text),
        "local_hash_not_timing_evidence": "does not prove its creation time" in tex
        and "no public preregistration or independently timestamped target lock is claimed" in tex,
        "separate_geo_not_independent_patients": overlap_audit["manuscript_language"].eq(
            "separate GEO datasets; no claim of independent patients"
        ).all()
        and "separate GEO datasets rather than independent patient cohorts" in tex,
        "retrospective_search_log_disclosed": len(eligibility) >= 5
        and eligibility["search_log_status"].str.contains("retrospectively reconstructed").all()
        and "retrospectively reconstructed GEO/PubMed search log" in tex,
        "github_code_location_is_factual": "https://github.com/xumedlab/DRnet" in tex
        and "unverified release tag" in tex,
        "auc_omitted_from_abstract": "AUC" not in re.search(
            r"\\begin\{abstract\}(.*?)\\end\{abstract\}", tex, re.S
        ).group(1),
    }
    checks = {key: bool(value) for key, value in checks.items()}
    failed = [key for key, value in checks.items() if not value]
    payload: dict[str, object] = {
        "status": "PASS" if not failed else "FAIL",
        "abstract_word_count": abstract_word_count(tex),
        "protocol_sha256": sha256(PROTOCOL),
        "workbook_sheets": sheets,
        "checks": checks,
        "failed_checks": failed,
    }
    REPORT.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if failed:
        raise AssertionError(f"Final submission consistency checks failed: {failed}")
    return payload


if __name__ == "__main__":
    print(json.dumps(validate(), ensure_ascii=False, indent=2))
