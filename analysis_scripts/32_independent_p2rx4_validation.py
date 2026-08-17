#!/usr/bin/env python3
"""Single-target disease-state checks for P2RX4 in separate GEO datasets.

The script evaluates exactly one gene (P2RX4) in deposited processed matrices.
It never treats normalized or abundance values as raw counts.  GSE276892 raw
reads were reconstructed separately by script 38; that post-protocol workflow
is intentionally kept distinct from these processed-value checks.  For
GSE179568, individual clinical covariates are restored from Supplementary
Table 1 and used for explicitly sensitivity-only robust analyses.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import mannwhitneyu, norm, spearmanr
from statsmodels.stats.multitest import multipletests


TARGET_GENE = "P2RX4"
DEFAULT_SEED = 20260813
DEFAULT_BOOTSTRAP = 10_000
SEARCH_DATE = "2026-08-13"
SEARCH_DATABASES = "NCBI GEO and PubMed"
SEARCH_STRING = (
    'human diabetic retinopathy AND (retina OR vitreous OR membrane) AND '
    '(RNA-seq OR transcriptomics)'
)

PDR_276892 = [
    ("PDR_S5", 77, "M", "PDR"),
    ("PDR_S6", 40, "M", "PDR"),
    ("PDR_S7", 56, "M", "PDR"),
    ("PDR_S8", 29, "F", "PDR"),
    ("PDR_S9", 66, "F", "PDR"),
    ("PDR_S10", 62, "M", "PDR"),
    ("PDR_S11", 63, "M", "PDR"),
    ("PDR_S12", 64, "F", "PDR"),
]
CONTROL_276892 = [
    ("MP_S13", 79, "F", "macular pucker; type II diabetes without PDR"),
    ("MP_S14", 72, "M", "macular pucker"),
    ("MP_S15", 73, "M", "macular pucker"),
    ("MP_S16", 73, "F", "macular pucker"),
    ("MH_S17", 68, "M", "macular hole"),
    ("MH_S18", 88, "M", "macular hole"),
    ("MH_S19", 78, "F", "macular hole"),
    ("MH_S20", 69, "M", "macular hole"),
    ("MH_S21", 72, "F", "macular hole"),
]
REUSED_CONTROL_IDS_276892 = {
    "MP_S14",
    "MP_S15",
    "MP_S16",
    "MH_S18",
    "MH_S19",
    "MH_S20",
    "MH_S21",
}
NEW_CONTROL_IDS_276892 = {"MP_S13", "MH_S17"}

# Audited transcription of the 24 RNA-seq rows in GSE179568 Supplementary
# Table 1 (Table 1.pdf, SHA256
# 375C2F9A341DB8AA90CF2B37847DE31B8FAE864B4D06C1BD2156E19B8D4D91B4).
# The row order is the deposited matrix order: PDR 1--7, MP 20--29, MH 35--41.
GSE179568_CLINICAL = [
    ("PDR_S1", 1, 40, "M", "I", "PDR", "phakic", 0, 0, 1),
    ("PDR_S2", 2, 35, "M", "I", "PDR", "phakic", 0, 0, 1),
    ("PDR_S3", 3, 29, "M", "I", "PDR", "phakic", 1, 1, 1),
    ("PDR_S4", 4, 29, "F", "I", "PDR", "phakic", 0, 0, 0),
    ("PDR_S5", 5, 56, "M", "II", "PDR", "phakic", 0, 0, 0),
    ("PDR_S6", 6, 20, "F", "I", "PDR", "phakic", 0, 0, 1),
    ("PDR_S7", 7, 60, "M", "I", "PDR", "pseudophakic", 0, 1, 1),
    ("Gliose_S1", 20, 71, "M", "none reported", "MP", "phakic", 0, 0, 0),
    ("Gliose_S2", 21, 83, "M", "none reported", "MP", "pseudophakic", 0, 0, 0),
    ("Gliose_S3", 22, 79, "M", "none reported", "MP", "phakic", 0, 0, 0),
    ("Gliose_S4", 23, 66, "F", "none reported", "MP", "phakic", 0, 0, 0),
    ("Gliose_S5", 24, 72, "F", "none reported", "MP", "phakic", 0, 0, 0),
    ("Gliose_S6", 25, 62, "M", "none reported", "MP", "pseudophakic", 0, 0, 0),
    ("Gliose_S7", 26, 62, "M", "none reported", "MP", "phakic", 0, 0, 0),
    ("Gliose_S8", 27, 75, "M", "none reported", "MP", "phakic", 0, 0, 0),
    ("Gliose_S9", 28, 70, "F", "none reported", "MP", "pseudophakic", 0, 0, 0),
    ("Gliose_S10", 29, 69, "M", "none reported", "MP", "phakic", 0, 0, 0),
    ("ILM_S1", 35, 65, "F", "none reported", "MH", "phakic", 0, 0, 0),
    ("ILM_S2", 36, 69, "M", "none reported", "MH", "phakic", 0, 0, 0),
    ("ILM_S3", 37, 79, "F", "none reported", "MH", "pseudophakic", 0, 0, 0),
    ("ILM_S4", 38, 67, "F", "none reported", "MH", "phakic", 0, 0, 0),
    ("ILM_S5", 39, 67, "M", "none reported", "MH", "phakic", 0, 0, 0),
    ("ILM_S6", 40, 61, "F", "none reported", "MH", "phakic", 0, 0, 0),
    ("ILM_S7", 41, 66, "F", "II", "MH", "pseudophakic", 0, 0, 0),
]


@dataclass(frozen=True)
class Comparison:
    dataset: str
    name: str
    case_label: str
    control_label: str
    case_ids: tuple[str, ...]
    control_ids: tuple[str, ...]
    role: str
    scale: str


def parse_args() -> argparse.Namespace:
    package = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path, default=package)
    parser.add_argument("--bootstrap", type=int, default=DEFAULT_BOOTSTRAP)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def cliff_delta(case: np.ndarray, control: np.ndarray) -> float:
    comparisons = np.sign(case[:, None] - control[None, :])
    return float(comparisons.mean())


def bootstrap_differences(
    case: np.ndarray,
    control: np.ndarray,
    iterations: int,
    rng: np.random.Generator,
) -> dict[str, float]:
    case_draws = rng.choice(case, size=(iterations, case.size), replace=True)
    control_draws = rng.choice(control, size=(iterations, control.size), replace=True)
    mean_diffs = case_draws.mean(axis=1) - control_draws.mean(axis=1)
    median_diffs = np.median(case_draws, axis=1) - np.median(
        control_draws, axis=1
    )
    log_case = np.log2(case_draws + 1.0)
    log_control = np.log2(control_draws + 1.0)
    log_mean_diffs = log_case.mean(axis=1) - log_control.mean(axis=1)
    return {
        "bootstrap_mean_difference_ci_low": float(np.quantile(mean_diffs, 0.025)),
        "bootstrap_mean_difference_ci_high": float(np.quantile(mean_diffs, 0.975)),
        "bootstrap_median_difference_ci_low": float(
            np.quantile(median_diffs, 0.025)
        ),
        "bootstrap_median_difference_ci_high": float(
            np.quantile(median_diffs, 0.975)
        ),
        "bootstrap_log2_mean_difference_ci_low": float(
            np.quantile(log_mean_diffs, 0.025)
        ),
        "bootstrap_log2_mean_difference_ci_high": float(
            np.quantile(log_mean_diffs, 0.975)
        ),
    }


def compare_groups(
    comparison: Comparison,
    values: dict[str, float],
    iterations: int,
    rng: np.random.Generator,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    case = np.asarray([values[sample] for sample in comparison.case_ids], dtype=float)
    control = np.asarray(
        [values[sample] for sample in comparison.control_ids], dtype=float
    )
    if not np.isfinite(np.r_[case, control]).all() or (np.r_[case, control] < 0).any():
        raise ValueError(f"Invalid deposited values for {comparison.name}")

    mw_one = mannwhitneyu(case, control, alternative="greater", method="auto")
    mw_two = mannwhitneyu(case, control, alternative="two-sided", method="auto")
    log_case = np.log2(case + 1.0)
    log_control = np.log2(control + 1.0)
    result: dict[str, object] = {
        "dataset": comparison.dataset,
        "comparison": comparison.name,
        "analysis_role": comparison.role,
        "case_label": comparison.case_label,
        "control_label": comparison.control_label,
        "deposited_value_scale": comparison.scale,
        "target_gene": TARGET_GENE,
        "n_case": case.size,
        "n_control": control.size,
        "case_mean": float(case.mean()),
        "control_mean": float(control.mean()),
        "mean_difference": float(case.mean() - control.mean()),
        "case_median": float(np.median(case)),
        "control_median": float(np.median(control)),
        "median_difference": float(np.median(case) - np.median(control)),
        "mean_log2_value_difference": float(log_case.mean() - log_control.mean()),
        "cliff_delta": cliff_delta(case, control),
        "mann_whitney_u": float(mw_one.statistic),
        "mann_whitney_p_one_sided_greater": float(mw_one.pvalue),
        "mann_whitney_p_two_sided": float(mw_two.pvalue),
    }
    result.update(bootstrap_differences(case, control, iterations, rng))

    loo_rows: list[dict[str, object]] = []
    all_ids = [*(('case', sample) for sample in comparison.case_ids),
               *(('control', sample) for sample in comparison.control_ids)]
    for omitted_group, omitted_id in all_ids:
        case_keep = np.asarray(
            [values[s] for s in comparison.case_ids if s != omitted_id], dtype=float
        )
        control_keep = np.asarray(
            [values[s] for s in comparison.control_ids if s != omitted_id], dtype=float
        )
        mean_diff = float(case_keep.mean() - control_keep.mean())
        median_diff = float(np.median(case_keep) - np.median(control_keep))
        loo_rows.append(
            {
                "dataset": comparison.dataset,
                "comparison": comparison.name,
                "omitted_group": omitted_group,
                "omitted_sample_id": omitted_id,
                "mean_difference_after_omission": mean_diff,
                "median_difference_after_omission": median_diff,
                "mean_log2_value_difference_after_omission": float(
                    np.log2(case_keep + 1).mean()
                    - np.log2(control_keep + 1).mean()
                ),
                "mean_direction_positive": mean_diff > 0,
                "median_direction_positive": median_diff > 0,
            }
        )
    result["loo_omissions"] = len(loo_rows)
    result["loo_mean_direction_positive_n"] = sum(
        bool(row["mean_direction_positive"]) for row in loo_rows
    )
    result["loo_median_direction_positive_n"] = sum(
        bool(row["median_direction_positive"]) for row in loo_rows
    )
    result["loo_mean_direction_positive_fraction"] = float(
        np.mean([row["mean_direction_positive"] for row in loo_rows])
    )
    result["loo_median_direction_positive_fraction"] = float(
        np.mean([row["median_direction_positive"] for row in loo_rows])
    )
    return result, loo_rows


def _table_number(value: object) -> float:
    """Parse a scalar Table 1 value while preserving unavailable entries."""
    parsed = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    return float(parsed) if pd.notna(parsed) else np.nan


def parse_table1_metadata(article_table: Path) -> pd.DataFrame:
    table = pd.read_html(article_table)[0]
    records: list[dict[str, object]] = []
    for _, row in table.iterrows():
        sample_text = str(row.iloc[0]).replace("*", "").strip()
        if not sample_text.isdigit() or int(sample_text) < 5:
            continue
        records.append(
            {
                "article_sample_number": int(sample_text),
                "age": int(row.iloc[1]),
                "sex": str(row.iloc[2]).strip(),
                "dm_type": str(row.iloc[3]).strip(),
                "ocular_diagnosis": str(row.iloc[4]).strip(),
                "cell_count": _table_number(row.iloc[12]),
                "rna_concentration_pg_per_ul": _table_number(row.iloc[13]),
            }
        )
    result = pd.DataFrame(records).drop_duplicates("article_sample_number")
    if result["article_sample_number"].tolist() != list(range(5, 22)):
        raise AssertionError("Could not reconstruct Table 1 samples 5-21")
    expected = pd.DataFrame(
        PDR_276892 + CONTROL_276892,
        columns=["sample_id", "age", "sex", "diagnosis_detail"],
    )
    expected["article_sample_number"] = expected["sample_id"].str.extract(
        r"S(\d+)$"
    )[0].astype(int)
    merged = expected.merge(result, on="article_sample_number", suffixes=("", "_table"))
    if not (merged["age"] == merged["age_table"]).all():
        raise AssertionError("Age mapping disagrees with downloaded article Table 1")
    if not (merged["sex"] == merged["sex_table"]).all():
        raise AssertionError("Sex mapping disagrees with downloaded article Table 1")
    if merged[["cell_count", "rna_concentration_pg_per_ul"]].isna().any().any():
        raise AssertionError("Table 1 cell-count or RNA-concentration mapping is incomplete")
    return merged


def load_gse276892(data_dir: Path) -> tuple[pd.DataFrame, dict[str, float], dict[str, float]]:
    matrix = pd.read_csv(
        data_dir / "GSE276892_normal_data.csv.gz", sep=";", decimal=","
    )
    target = matrix.loc[matrix["Symbol"].astype(str).str.upper().eq(TARGET_GENE)]
    if len(target) != 1:
        raise AssertionError("GSE276892 must contain exactly one P2RX4 row")
    target = target.iloc[0]
    metadata = parse_table1_metadata(
        data_dir / "GSE276892_primary_article_table1.html"
    )
    sample_columns = [row[0] for row in PDR_276892 + CONTROL_276892]
    missing = sorted(set(sample_columns) - set(matrix.columns))
    if missing:
        raise AssertionError(f"Missing GSE276892 columns: {missing}")
    values = {sample: float(target[sample]) for sample in sample_columns}
    metadata["dataset"] = "GSE276892"
    metadata["disease_group"] = np.where(
        metadata["sample_id"].str.startswith("PDR"), "PDR", "control"
    )
    metadata["control_subtype"] = np.select(
        [
            metadata["sample_id"].str.startswith("MP"),
            metadata["sample_id"].str.startswith("MH"),
        ],
        ["macular pucker", "macular hole"],
        default="not applicable",
    )
    metadata["tissue_or_cell_fraction"] = "vitreous hyalocytes"
    metadata["deposited_value_scale"] = "DESeq2-normalized reads (not raw counts)"
    metadata["p2rx4_value"] = metadata["sample_id"].map(values)
    metadata["log2_p2rx4_plus_1"] = np.log2(metadata["p2rx4_value"] + 1)
    metadata["sex_male"] = metadata["sex"].eq("M").astype(int)
    metadata["source_dataset"] = np.where(
        metadata["sample_id"].isin(REUSED_CONTROL_IDS_276892),
        "reused GSE147657 control",
        "new GSE276892 profile",
    )
    metadata["reused_control"] = metadata["sample_id"].isin(
        REUSED_CONTROL_IDS_276892
    ).astype(int)
    metadata["protocol_role"] = (
        "authors state locally specified before processed-matrix inspection; "
        "no independent timestamp"
    )
    source = {
        "source_reported_log2_fold_change": float(target["log2FC"]),
        "source_reported_adjusted_p": float(target["padj"]),
        "source_reported_mean_pdr": float(target["mean_Diabetes"]),
        "source_reported_mean_control": float(target["mean_Control"]),
    }
    pdr_values = np.asarray([values[row[0]] for row in PDR_276892])
    control_values = np.asarray([values[row[0]] for row in CONTROL_276892])
    if not np.isclose(pdr_values.mean(), source["source_reported_mean_pdr"], rtol=1e-8):
        raise AssertionError("GSE276892 PDR mean does not match deposited summary")
    if not np.isclose(
        control_values.mean(), source["source_reported_mean_control"], rtol=1e-8
    ):
        raise AssertionError("GSE276892 control mean does not match deposited summary")
    return metadata, values, source


def _ols_hc3_terms(
    outcome: pd.Series,
    design: pd.DataFrame,
    prefix: str,
) -> dict[str, float]:
    design = sm.add_constant(design, has_constant="add")
    fit = sm.OLS(outcome.astype(float), design).fit(cov_type="HC3", use_t=True)
    ci = fit.conf_int().loc["pdr"]
    beta = float(fit.params["pdr"])
    p_two = float(fit.pvalues["pdr"])
    p_one = p_two / 2 if beta >= 0 else 1 - p_two / 2
    return {
        f"{prefix}_beta": beta,
        f"{prefix}_ci_low": float(ci.iloc[0]),
        f"{prefix}_ci_high": float(ci.iloc[1]),
        f"{prefix}_p_two_sided": p_two,
        f"{prefix}_p_one_sided_greater": p_one,
        f"{prefix}_residual_df": float(fit.df_resid),
        f"{prefix}_design_condition_number": float(np.linalg.cond(design)),
    }


def adjusted_gse276892(metadata: pd.DataFrame) -> dict[str, float]:
    base_design = pd.DataFrame(
        {
            "pdr": metadata["disease_group"].eq("PDR").astype(float),
            "age_per_10y_centered": (metadata["age"] - metadata["age"].mean()) / 10,
            "sex_male": metadata["sex_male"].astype(float),
            "reused_control": metadata["reused_control"].astype(float),
        },
        index=metadata.index,
    )
    outcome = metadata["log2_p2rx4_plus_1"].astype(float)
    age_sex_design = base_design[["pdr", "age_per_10y_centered", "sex_male"]]
    source_design = base_design[["pdr", "reused_control"]]
    source_age_sex_design = base_design[
        ["pdr", "reused_control", "age_per_10y_centered", "sex_male"]
    ]

    result: dict[str, float] = {}
    result.update(
        _ols_hc3_terms(
            outcome,
            age_sex_design,
            "age_sex_adjusted_ols_hc3_t",
        )
    )
    result.update(
        _ols_hc3_terms(
            outcome,
            source_design,
            "source_adjusted_ols_hc3_t",
        )
    )
    result.update(
        _ols_hc3_terms(
            outcome,
            source_age_sex_design,
            "source_age_sex_adjusted_ols_hc3_t",
        )
    )

    robust_design = sm.add_constant(age_sex_design, has_constant="add")
    rlm = sm.RLM(outcome, robust_design, M=sm.robust.norms.HuberT()).fit(cov="H1")
    rlm_beta = float(rlm.params["pdr"])
    rlm_se = float(rlm.bse["pdr"])
    rlm_z = rlm_beta / rlm_se
    rlm_p_two = float(2 * norm.sf(abs(rlm_z)))
    result.update(
        {
        "robust_huber_beta": rlm_beta,
        "robust_huber_se_h1": rlm_se,
        "robust_huber_z": rlm_z,
        "robust_huber_p_two_sided_asymptotic": rlm_p_two,
        "robust_huber_ci_low_asymptotic": rlm_beta - 1.96 * rlm_se,
        "robust_huber_ci_high_asymptotic": rlm_beta + 1.96 * rlm_se,
        }
    )
    return result


def gse276892_qc_correlations(metadata: pd.DataFrame) -> pd.DataFrame:
    """Relate P2RX4 to the two per-patient QC proxies reported in Table 1."""
    rows: list[dict[str, object]] = []
    for qc_metric in ["cell_count", "rna_concentration_pg_per_ul"]:
        for expression_scale in ["p2rx4_value", "log2_p2rx4_plus_1"]:
            frame = metadata[[expression_scale, qc_metric]].dropna()
            test = spearmanr(
                frame[expression_scale].to_numpy(dtype=float),
                frame[qc_metric].to_numpy(dtype=float),
            )
            rows.append(
                {
                    "dataset": "GSE276892",
                    "expression_scale": expression_scale,
                    "qc_metric": qc_metric,
                    "n_profiles": len(frame),
                    "spearman_rho": float(test.statistic),
                    "spearman_p_two_sided": float(test.pvalue),
                    "scope": (
                        "descriptive QC correlation; mapping rate, library size, and "
                        "detected-gene count were not deposited per patient"
                    ),
                }
            )
    return pd.DataFrame(rows)


def parse_gse179568_sample_map(matrix_file: Path) -> pd.DataFrame:
    lines: dict[str, list[str]] = {}
    import gzip

    with gzip.open(matrix_file, "rt", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            if line.startswith("!Sample_title") or line.startswith("!Sample_geo_accession"):
                fields = [value.strip('"') for value in line.rstrip().split("\t")[1:]]
                lines[line.split("\t", 1)[0]] = fields
    titles = lines["!Sample_title"]
    accessions = lines["!Sample_geo_accession"]
    if len(titles) != 24 or len(accessions) != 24:
        raise AssertionError("GSE179568 sample map must contain 24 patients")
    matrix_ids = [
        *(f"PDR_S{i}" for i in range(1, 8)),
        *(f"Gliose_S{i}" for i in range(1, 11)),
        *(f"ILM_S{i}" for i in range(1, 8)),
    ]
    return pd.DataFrame(
        {"sample_id": matrix_ids, "sample_title": titles, "geo_accession": accessions}
    )


def load_gse179568(
    data_dir: Path,
) -> tuple[pd.DataFrame, dict[str, float], dict[str, float]]:
    cohort_dir = data_dir / "GSE179568"
    matrix = pd.read_csv(
        cohort_dir / "GSE179568_data.csv.gz", sep=";", decimal=","
    )
    target = matrix.loc[matrix["Symbol"].astype(str).str.upper().eq(TARGET_GENE)]
    if len(target) != 1:
        raise AssertionError("GSE179568 must contain exactly one P2RX4 row")
    target = target.iloc[0]
    sample_map = parse_gse179568_sample_map(
        cohort_dir / "GSE179568_series_matrix.txt.gz"
    )
    values = {sample: float(target[sample]) for sample in sample_map["sample_id"]}
    sample_map["dataset"] = "GSE179568"
    sample_map["disease_group"] = np.select(
        [
            sample_map["sample_id"].str.startswith("PDR"),
            sample_map["sample_id"].str.startswith("Gliose"),
        ],
        ["PDR", "control"],
        default="control",
    )
    sample_map["control_subtype"] = np.select(
        [
            sample_map["sample_id"].str.startswith("Gliose"),
            sample_map["sample_id"].str.startswith("ILM"),
        ],
        ["macular pucker", "macular hole"],
        default="not applicable",
    )
    sample_map["tissue_or_cell_fraction"] = np.select(
        [
            sample_map["sample_id"].str.startswith("PDR"),
            sample_map["sample_id"].str.startswith("Gliose"),
        ],
        ["retinal neovascularization membrane", "epiretinal membrane"],
        default="inner limiting membrane",
    )
    sample_map["deposited_value_scale"] = "DESeq2-normalized reads (not raw counts)"
    sample_map["p2rx4_value"] = sample_map["sample_id"].map(values)
    sample_map["log2_p2rx4_plus_1"] = np.log2(sample_map["p2rx4_value"] + 1)
    sample_map["protocol_role"] = "post-protocol secondary cohort"
    clinical_columns = [
        "sample_id",
        "supplementary_table_row",
        "age",
        "sex",
        "dm_type",
        "ocular_diagnosis",
        "lens_status",
        "previous_vitrectomy",
        "previous_anti_vegf_over_3_months",
        "previous_prp",
    ]
    clinical = pd.DataFrame(GSE179568_CLINICAL, columns=clinical_columns)
    sample_map = sample_map.merge(
        clinical,
        on="sample_id",
        how="left",
        validate="one_to_one",
    )
    if sample_map["age"].isna().any():
        raise AssertionError("GSE179568 clinical mapping is incomplete")
    sample_map["sex_male"] = sample_map["sex"].eq("M").astype(int)
    sample_map["pseudophakic"] = sample_map["lens_status"].eq(
        "pseudophakic"
    ).astype(int)
    sample_map["clinical_source"] = (
        "Source article Supplementary Table 1; audited fixed transcription"
    )
    source = {
        "source_reported_log2_fold_change_pdr_vs_ilm": float(
            target["log2FC_PDR_ILM"]
        ),
        "source_reported_adjusted_p_pdr_vs_ilm": float(target["padj_PDR_ILM"]),
        "source_reported_mean_pdr": float(target["mean_PDR"]),
        "source_reported_mean_macular_pucker": float(target["mean_Gliosis"]),
        "source_reported_mean_ilm": float(target["mean_ILM"]),
    }
    for key, ids in {
        "source_reported_mean_pdr": [f"PDR_S{i}" for i in range(1, 8)],
        "source_reported_mean_macular_pucker": [
            f"Gliose_S{i}" for i in range(1, 11)
        ],
        "source_reported_mean_ilm": [f"ILM_S{i}" for i in range(1, 8)],
    }.items():
        observed = float(np.mean([values[sample] for sample in ids]))
        if not np.isclose(observed, source[key], rtol=1e-8):
            raise AssertionError(f"GSE179568 deposited mean mismatch: {key}")
    return sample_map, values, source


def _huber_terms(
    outcome: pd.Series,
    design: pd.DataFrame,
    prefix: str,
) -> dict[str, float]:
    """Huber M-estimator sensitivity with asymptotic H1 covariance."""
    design = sm.add_constant(design, has_constant="add")
    fit = sm.RLM(
        outcome.astype(float),
        design.astype(float),
        M=sm.robust.norms.HuberT(),
    ).fit(cov="H1")
    beta = float(fit.params["pdr"])
    se = float(fit.bse["pdr"])
    z_stat = beta / se
    return {
        f"{prefix}_beta": beta,
        f"{prefix}_se_h1": se,
        f"{prefix}_z": z_stat,
        f"{prefix}_p_two_sided_asymptotic": float(2 * norm.sf(abs(z_stat))),
        f"{prefix}_ci_low_asymptotic": beta - 1.96 * se,
        f"{prefix}_ci_high_asymptotic": beta + 1.96 * se,
    }


def _standardized_mean_difference(case: pd.Series, control: pd.Series) -> float:
    pooled_sd = math.sqrt((float(case.var(ddof=1)) + float(control.var(ddof=1))) / 2)
    if pooled_sd == 0:
        return 0.0 if float(case.mean()) == float(control.mean()) else math.inf
    return float((case.mean() - control.mean()) / pooled_sd)


def gse179568_clinical_sensitivity(
    metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Quantify clinical confounding without claiming covariate exchangeability.

    Age has no common support between PDR and either control compartment.  The
    adjusted regressions are therefore reported as extrapolative sensitivities,
    not as deconfounded estimates.
    """
    comparisons = {
        "RNV vs macular-pucker membrane": "macular pucker",
        "RNV vs macular-hole ILM": "macular hole",
        "RNV vs pooled membrane controls": "pooled",
    }
    result_rows: list[dict[str, object]] = []
    balance_rows: list[dict[str, object]] = []
    for name, control_subtype in comparisons.items():
        control_mask = metadata["disease_group"].eq("control")
        if control_subtype != "pooled":
            control_mask &= metadata["control_subtype"].eq(control_subtype)
        frame = metadata.loc[metadata["disease_group"].eq("PDR") | control_mask].copy()
        frame["pdr"] = frame["disease_group"].eq("PDR").astype(float)
        frame["age_per_10y_centered"] = (frame["age"] - frame["age"].mean()) / 10
        case = frame.loc[frame["pdr"].eq(1)]
        control = frame.loc[frame["pdr"].eq(0)]
        outcome = frame["log2_p2rx4_plus_1"].astype(float)

        case_age_min = float(case["age"].min())
        case_age_max = float(case["age"].max())
        control_age_min = float(control["age"].min())
        control_age_max = float(control["age"].max())
        overlap_low = max(case_age_min, control_age_min)
        overlap_high = min(case_age_max, control_age_max)
        age_common_support = overlap_low <= overlap_high

        row: dict[str, object] = {
            "dataset": "GSE179568",
            "comparison": name,
            "n_pdr": len(case),
            "n_control": len(control),
            "age_pdr_mean": float(case["age"].mean()),
            "age_control_mean": float(control["age"].mean()),
            "age_pdr_range": f"{case_age_min:.0f}-{case_age_max:.0f}",
            "age_control_range": f"{control_age_min:.0f}-{control_age_max:.0f}",
            "age_common_support": bool(age_common_support),
            "age_overlap_interval": (
                f"{overlap_low:.0f}-{overlap_high:.0f}"
                if age_common_support
                else "none"
            ),
            "age_standardized_mean_difference": _standardized_mean_difference(
                case["age"], control["age"]
            ),
            "identifiability_statement": (
                "No age common support; adjusted disease coefficients require "
                "extrapolation and cannot isolate disease from age or tissue composition."
                if not age_common_support
                else "Age common support is present."
            ),
        }
        row.update(
            _ols_hc3_terms(outcome, frame[["pdr"]], "unadjusted_ols_hc3_t")
        )
        row.update(
            _ols_hc3_terms(
                outcome,
                frame[["pdr", "age_per_10y_centered"]],
                "age_adjusted_ols_hc3_t",
            )
        )
        row.update(
            _ols_hc3_terms(
                outcome,
                frame[["pdr", "age_per_10y_centered", "sex_male"]],
                "age_sex_adjusted_ols_hc3_t",
            )
        )
        row.update(
            _ols_hc3_terms(
                outcome,
                frame[
                    [
                        "pdr",
                        "age_per_10y_centered",
                        "sex_male",
                        "previous_prp",
                        "previous_anti_vegf_over_3_months",
                        "previous_vitrectomy",
                    ]
                ],
                "age_sex_treatment_adjusted_ols_hc3_t",
            )
        )
        row.update(
            _huber_terms(
                outcome,
                frame[["pdr", "age_per_10y_centered", "sex_male"]],
                "age_sex_adjusted_huber_h1",
            )
        )
        mw = mannwhitneyu(
            case["p2rx4_value"],
            control["p2rx4_value"],
            alternative="two-sided",
            method="auto",
        )
        row.update(
            {
                "mann_whitney_u": float(mw.statistic),
                "mann_whitney_p_two_sided": float(mw.pvalue),
                "cliff_delta": cliff_delta(
                    case["p2rx4_value"].to_numpy(dtype=float),
                    control["p2rx4_value"].to_numpy(dtype=float),
                ),
                "propensity_or_weighting_attempted": False,
                "propensity_or_weighting_reason": (
                    "Not estimable: age distributions do not overlap."
                ),
                "interpretation_role": (
                    "post-protocol, tissue-composition- and age-confounded sensitivity"
                ),
            }
        )
        result_rows.append(row)

        for covariate in [
            "age",
            "sex_male",
            "dm_type_I",
            "previous_prp",
            "previous_anti_vegf_over_3_months",
            "previous_vitrectomy",
            "pseudophakic",
        ]:
            if covariate == "dm_type_I":
                case_values = case["dm_type"].eq("I").astype(float)
                control_values = control["dm_type"].eq("I").astype(float)
            else:
                case_values = case[covariate].astype(float)
                control_values = control[covariate].astype(float)
            balance_rows.append(
                {
                    "dataset": "GSE179568",
                    "comparison": name,
                    "covariate": covariate,
                    "pdr_mean_or_proportion": float(case_values.mean()),
                    "control_mean_or_proportion": float(control_values.mean()),
                    "standardized_mean_difference": _standardized_mean_difference(
                        case_values, control_values
                    ),
                    "n_pdr": len(case_values),
                    "n_control": len(control_values),
                }
            )
    return pd.DataFrame(result_rows), pd.DataFrame(balance_rows)


def load_gse94019(data_dir: Path) -> tuple[pd.DataFrame, dict[str, float], int]:
    matrix = pd.read_csv(
        data_dir / "GSE94019_Partek_EM_gene_reads.txt.gz", sep="\t"
    )
    target = matrix.loc[
        matrix["Gene Symbol"].fillna("").astype(str).str.upper().eq(TARGET_GENE)
    ]
    if len(target) != 6:
        raise AssertionError("Expected six deposited P2RX4 transcript rows in GSE94019")
    metadata_columns = {"Gene Symbol", "Transcript ID", "Length"}
    sample_columns = [column for column in matrix if column not in metadata_columns]
    values = {sample: float(target[sample].sum()) for sample in sample_columns}
    case_ids = [sample for sample in sample_columns if sample.startswith("MEEI")]
    control_ids = [sample for sample in sample_columns if sample not in case_ids]
    if len(case_ids) != 9 or len(control_ids) != 4:
        raise AssertionError("GSE94019 deposited matrix must map to 9 PDR and 4 controls")
    metadata = pd.DataFrame({"sample_id": sample_columns})
    metadata["dataset"] = "GSE94019"
    metadata["disease_group"] = np.where(
        metadata["sample_id"].isin(case_ids), "PDR", "control"
    )
    metadata["control_subtype"] = np.where(
        metadata["sample_id"].isin(control_ids),
        "non-diabetic post-mortem retina",
        "not applicable",
    )
    metadata["tissue_or_cell_fraction"] = np.where(
        metadata["sample_id"].isin(case_ids),
        "fibrovascular-membrane CD31+ cells",
        "retinal CD31+ cells",
    )
    metadata["deposited_value_scale"] = (
        "sum of six Partek E/M transcript abundance measurements; unit not "
        "further specified by GEO (not labelled CPM or raw counts)"
    )
    metadata["p2rx4_value"] = metadata["sample_id"].map(values)
    metadata["log2_p2rx4_plus_1"] = np.log2(metadata["p2rx4_value"] + 1)
    metadata["protocol_role"] = "post-protocol secondary cohort"
    metadata["geo_accession"] = np.nan
    for column in ["age", "sex", "sex_male", "dm_type", "ocular_diagnosis"]:
        metadata[column] = np.nan
    return metadata, values, len(target)


def dataset_eligibility() -> pd.DataFrame:
    rows = [
        {
            "dataset": "GSE276892",
            "search_date": SEARCH_DATE,
            "search_databases": SEARCH_DATABASES,
            "search_string": SEARCH_STRING,
            "search_log_status": "retrospectively reconstructed on the stated date",
            "screening_decision": "included as the locally specified primary external dataset",
            "exclusion_reason": "not excluded",
            "p2rx4_value_inspection_status": (
                "authors state that processed P2RX4 values were not inspected before "
                "local specification; no independent timestamp proves the sequence"
            ),
            "primary_comparison": "PDR versus all surgical-control vitreous hyalocytes",
            "candidate_selection_independent": True,
            "human": True,
            "disease": "PDR",
            "case_material": "fluorescence-sorted vitreous hyalocytes",
            "control_material": "vitreous hyalocytes from macular-pucker/hole surgery",
            "n_case": 8,
            "n_control": 9,
            "analysis_status": "locally specified separate GEO ocular-compartment dataset",
            "eligibility_judgment": "eligible for an external single-target check",
            "not_a_claim": "not whole-retina replication or proof of patient non-overlap",
            "principal_limitation": (
                "controls are older; seven of nine controls were reused from GSE147657; "
                "one control has type II diabetes; disease and source are strongly "
                "associated; the exact source count matrix was not deposited"
            ),
            "geo_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276892",
            "article_identifier": "PMID 39543723; DOI 10.1186/s12974-024-03291-5",
        },
        {
            "dataset": "GSE179568",
            "search_date": SEARCH_DATE,
            "search_databases": SEARCH_DATABASES,
            "search_string": SEARCH_STRING,
            "search_log_status": "retrospectively reconstructed on the stated date",
            "screening_decision": "included post-protocol as a secondary membrane-compartment dataset",
            "exclusion_reason": "not excluded",
            "p2rx4_value_inspection_status": (
                "unknown from contemporaneous records; dataset and comparison are "
                "therefore treated as post-protocol exploratory"
            ),
            "primary_comparison": "PDR RNV versus macular-pucker membrane compartment",
            "candidate_selection_independent": True,
            "human": True,
            "disease": "PDR",
            "case_material": "retinal neovascularization membranes",
            "control_material": "macular-pucker membranes; macular-hole ILM sensitivity",
            "n_case": 7,
            "n_control": 17,
            "analysis_status": "post-protocol secondary separate GEO dataset",
            "eligibility_judgment": "eligible for secondary directional support",
            "not_a_claim": (
                "not matched tissue, whole-retina replication, or proof of patient "
                "non-overlap with other Freiburg datasets"
            ),
            "principal_limitation": (
                "small groups, severe age and cellular-composition differences, "
                "normalized reads only, and comparison rule set after the protocol"
            ),
            "geo_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179568",
            "article_identifier": "PMID 34795670 and PMID 38153746",
        },
        {
            "dataset": "GSE94019",
            "search_date": SEARCH_DATE,
            "search_databases": SEARCH_DATABASES,
            "search_string": SEARCH_STRING,
            "search_log_status": "retrospectively reconstructed on the stated date",
            "screening_decision": "included post-protocol as a tissue-mismatched direction check",
            "exclusion_reason": "not excluded",
            "p2rx4_value_inspection_status": (
                "unknown from contemporaneous records; dataset and comparison are "
                "therefore treated as post-protocol exploratory"
            ),
            "primary_comparison": "PDR FVM CD31+ versus control-retina CD31+ profiles",
            "candidate_selection_independent": True,
            "human": True,
            "disease": "PDR",
            "case_material": "fibrovascular-membrane CD31+ cells",
            "control_material": "post-mortem non-diabetic retinal CD31+ cells",
            "n_case": 9,
            "n_control": 4,
            "analysis_status": "post-protocol secondary separate GEO dataset",
            "eligibility_judgment": "eligible for secondary directional support",
            "not_a_claim": "not matched-tissue or whole-retina replication",
            "principal_limitation": (
                "case/control tissue context differs; deposited abundance unit is not "
                "fully specified; source article figure used 8 PDR although GEO has 9"
            ),
            "geo_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94019",
            "article_identifier": "PMID 28400392; DOI 10.2337/db16-1035",
        },
        {
            "dataset": "GSE102485",
            "search_date": SEARCH_DATE,
            "search_databases": SEARCH_DATABASES,
            "search_string": SEARCH_STRING,
            "search_log_status": "retrospectively reconstructed on the stated date",
            "screening_decision": "excluded from the external validation chain",
            "exclusion_reason": (
                "no non-PDR control group and mixed retina/neovascular-membrane "
                "material; retained only as post-hoc membrane context"
            ),
            "p2rx4_value_inspection_status": (
                "P2RX4 values had been inspected before final exclusion; retained only "
                "as post-hoc context"
            ),
            "primary_comparison": "none",
            "candidate_selection_independent": False,
            "human": True,
            "disease": "diabetic and vein-occlusion membrane material",
            "case_material": "mixed retinal and neovascular proliferative membranes",
            "control_material": "none suitable",
            "n_case": np.nan,
            "n_control": 0,
            "analysis_status": "excluded; post-hoc context only",
            "eligibility_judgment": "not eligible for disease-state validation",
            "not_a_claim": "not external validation",
            "principal_limitation": "no valid disease-versus-control estimand",
            "geo_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102485",
            "article_identifier": "GSE102485",
        },
        {
            "dataset": "GSE130636",
            "search_date": SEARCH_DATE,
            "search_databases": SEARCH_DATABASES,
            "search_string": SEARCH_STRING,
            "search_log_status": "retrospectively reconstructed on the stated date",
            "screening_decision": "excluded from disease-state validation",
            "exclusion_reason": "normal-retina atlas with no DR disease contrast",
            "p2rx4_value_inspection_status": (
                "not applicable to disease-state eligibility; normal-retina atlas only"
            ),
            "primary_comparison": "none",
            "candidate_selection_independent": True,
            "human": True,
            "disease": "none",
            "case_material": "normal retinal single cells",
            "control_material": "not applicable",
            "n_case": 0,
            "n_control": 3,
            "analysis_status": "localization context only",
            "eligibility_judgment": "not eligible for disease-state validation",
            "not_a_claim": "not disease-state replication",
            "principal_limitation": "normal tissue only",
            "geo_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130636",
            "article_identifier": "PMID 31141684; DOI 10.1016/j.exer.2019.05.001",
        },
    ]
    return pd.DataFrame(rows)


def patient_overlap_audit() -> pd.DataFrame:
    """Record what public metadata can and cannot establish about patient overlap."""
    return pd.DataFrame(
        [
            {
                "dataset_pair": "GSE276892 versus GSE179568",
                "shared_institution_or_team": True,
                "shared_corresponding_unit": True,
                "exact_deposited_sample_id_overlap": False,
                "exact_pdr_age_sex_matches": "40/M; 56/M; 29/F",
                "collection_date_available_per_patient_in_both": False,
                "eye_or_surgery_identifier_available_in_both": False,
                "author_confirmation_obtained": False,
                "public_metadata_conclusion": (
                    "Patient non-overlap cannot be established or excluded from public "
                    "metadata; different GEO accessions and sample labels are insufficient."
                ),
                "manuscript_language": (
                    "separate GEO datasets; no claim of independent patients"
                ),
            }
        ]
    )


def protocol_deviations(protocol_hash: str, gse94019_transcripts: int) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "item": "single target",
                "protocol_specification": "P2RX4 only",
                "implemented_analysis": "P2RX4 only in every dataset",
                "deviation": False,
                "reason_and_consequence": (
                    "Single-target scope preserved; the local hash authenticates content "
                    "but does not independently establish creation time."
                ),
            },
            {
                "item": "GSE276892 primary count model",
                "protocol_specification": (
                    "negative-binomial disease-only model and one-sided Wald test"
                ),
                "implemented_analysis": (
                    "source-reported DESeq2 log2FC/padj, patient-level analyses of "
                    "deposited DESeq2-normalized values, and a separate post-protocol "
                    "SRA reconstruction in script 38"
                ),
                "deviation": True,
                "reason_and_consequence": (
                    "The exact source gene-count matrix and unadjusted P2RX4 Wald P were "
                    "not deposited. Reconstructed counts depend on post-protocol "
                    "alignment and quantification choices, so the resulting negative-"
                    "binomial estimate remains a protocol-deviation sensitivity rather "
                    "than the originally specified confirmatory test."
                ),
            },
            {
                "item": "GSE276892 source-dataset sensitivity",
                "protocol_specification": "not specified",
                "implemented_analysis": (
                    "seven reused GSE147657 controls and two new controls identified; "
                    "source-adjusted and source+age+sex HC3 models plus an 8-versus-2 "
                    "new-profile-only comparison"
                ),
                "deviation": True,
                "reason_and_consequence": (
                    "Disease and data source are strongly associated, so source-adjusted "
                    "estimates are low-precision sensitivities rather than confirmatory "
                    "deconfounding."
                ),
            },
            {
                "item": "GSE276892 covariate sensitivity",
                "protocol_specification": (
                    "age/sex-adjusted negative-binomial model if covariates complete"
                ),
                "implemented_analysis": (
                    "age/sex-adjusted OLS on log2(normalized value + 1), HC3 covariance "
                    "with finite-sample t inference; Huber robust regression sensitivity"
                ),
                "deviation": True,
                "reason_and_consequence": (
                    "Age and sex were complete, but only normalized values were "
                    "available. These models are labelled sensitivity analyses."
                ),
            },
            {
                "item": "additional cohorts",
                "protocol_specification": "GSE276892 only",
                "implemented_analysis": "GSE179568 and GSE94019 added after the local protocol",
                "deviation": True,
                "reason_and_consequence": (
                    "Both are post-protocol secondary directional checks and cannot "
                    "rescue or redefine the locally specified criterion."
                ),
            },
            {
                "item": "GSE94019 transcript aggregation",
                "protocol_specification": "not specified",
                "implemented_analysis": (
                    f"sum of {gse94019_transcripts} P2RX4 Partek E/M transcript "
                    "abundance rows per sample"
                ),
                "deviation": True,
                "reason_and_consequence": (
                    "GEO calls the values abundance measurements but does not identify "
                    "them as CPM; results retain that neutral unit description."
                ),
            },
            {
                "item": "GSE94019 sample count",
                "protocol_specification": "not specified",
                "implemented_analysis": "all 9 PDR and 4 control columns in GEO matrix",
                "deviation": True,
                "reason_and_consequence": (
                    "The source article figure caption reports 8 PDR and 4 controls, "
                    "but GEO provides 9 PDR columns without an exclusion flag. The "
                    "deposited-matrix analysis is reported transparently."
                ),
            },
            {
                "item": "local protocol integrity and timing",
                "protocol_specification": (
                    "authors state that P2RX4 was locally specified before processed "
                    "GSE276892 gene-value inspection"
                ),
                "implemented_analysis": f"SHA256 {protocol_hash}",
                "deviation": False,
                "reason_and_consequence": (
                    "The hash verifies the current file content, not when it was created. "
                    "No public registration or independent pre-inspection timestamp is "
                    "claimed."
                ),
            },
        ]
    )


def annotate_result_status(results: pd.DataFrame) -> pd.DataFrame:
    results = results.copy()
    results["interpretation"] = "post-protocol descriptive comparison"
    main_276 = results["comparison"].eq("PDR vs all surgical controls")
    results.loc[main_276, "interpretation"] = (
        "inconclusive processed-value check; source adjustment and log-scale results "
        "do not support a separable disease effect"
    )
    main_179 = results["comparison"].eq("RNV vs macular-pucker membrane")
    results.loc[main_179, "interpretation"] = (
        "post-protocol secondary membrane-compartment association with age and "
        "cell-composition confounding"
    )
    main_940 = results["comparison"].eq("PDR FVM CD31+ vs control retinal CD31+")
    results.loc[main_940, "interpretation"] = (
        "post-protocol secondary direction-only support with tissue-context mismatch"
    )
    return results


def add_post_protocol_holm(results: pd.DataFrame) -> pd.DataFrame:
    """Adjust the three predesignated post-protocol dataset-level comparisons."""
    results = results.copy()
    results["post_protocol_holm_family"] = False
    results["mann_whitney_p_two_sided_holm_post_protocol_3"] = np.nan
    family = (
        (results["dataset"].eq("GSE179568") & results["comparison"].eq(
            "RNV vs macular-pucker membrane"
        ))
        | (results["dataset"].eq("GSE179568") & results["comparison"].eq(
            "RNV vs macular-hole ILM"
        ))
        | (results["dataset"].eq("GSE94019") & results["comparison"].eq(
            "PDR FVM CD31+ vs control retinal CD31+"
        ))
    )
    if int(family.sum()) != 3:
        raise AssertionError("Expected exactly three post-protocol Holm comparisons")
    adjusted = multipletests(
        results.loc[family, "mann_whitney_p_two_sided"].to_numpy(dtype=float),
        method="holm",
    )[1]
    results.loc[family, "post_protocol_holm_family"] = True
    results.loc[family, "mann_whitney_p_two_sided_holm_post_protocol_3"] = adjusted
    return results


def plot_validation(sample_table: pd.DataFrame, results: pd.DataFrame, output: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10.5,
            "axes.titlesize": 12.5,
            "axes.labelsize": 11,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    colors = {"PDR": "#B33A3A", "control": "#3B6EA8"}
    rng = np.random.default_rng(4421)
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.5), constrained_layout=True)

    panels = [
        (
            "GSE276892",
            ["control", "PDR"],
            {"control": "Surgical controls\n(n=9)", "PDR": "PDR\n(n=8)"},
            "A  Vitreous hyalocytes",
            "PDR vs all surgical controls",
        ),
        (
            "GSE179568",
            ["macular pucker", "macular hole", "PDR"],
            {
                "macular pucker": "Macular pucker\n(n=10)",
                "macular hole": "Macular hole\n(n=7)",
                "PDR": "PDR RNV\n(n=7)",
            },
            "B  Ocular membranes",
            "RNV vs macular-pucker membrane",
        ),
        (
            "GSE94019",
            ["control", "PDR"],
            {"control": "Control retina\n(n=4)", "PDR": "PDR FVM\n(n=9)"},
            "C  CD31+ fractions",
            "PDR FVM CD31+ vs control retinal CD31+",
        ),
    ]
    for ax, (dataset, order, labels, title, comparison) in zip(axes, panels, strict=True):
        frame = sample_table.loc[sample_table["dataset"].eq(dataset)].copy()
        if dataset == "GSE179568":
            frame["plot_group"] = np.where(
                frame["disease_group"].eq("PDR"), "PDR", frame["control_subtype"]
            )
        else:
            frame["plot_group"] = frame["disease_group"]
        positions = np.arange(len(order), dtype=float)
        for x, group in zip(positions, order, strict=True):
            group_frame = frame.loc[frame["plot_group"].eq(group)].copy()
            values = group_frame["log2_p2rx4_plus_1"].to_numpy()
            color = colors["PDR"] if group == "PDR" else colors["control"]
            jitter = rng.uniform(-0.09, 0.09, size=len(values))
            ax.scatter(
                np.full(len(values), x) + jitter,
                values,
                s=30,
                color=color,
                edgecolor="white",
                linewidth=0.6,
                zorder=3,
            )
            if dataset == "GSE276892" and group == "PDR":
                lowest_index = int(np.argmin(values))
                lowest = group_frame.iloc[lowest_index]
                ax.annotate(
                    str(lowest["sample_id"]),
                    (x + jitter[lowest_index], values[lowest_index]),
                    xytext=(6, 5),
                    textcoords="offset points",
                    fontsize=9,
                    fontweight="bold",
                    color="#6E1F1F",
                )
            median = float(np.median(values))
            q1, q3 = np.quantile(values, [0.25, 0.75])
            ax.vlines(x, q1, q3, color=color, linewidth=5, alpha=0.34, zorder=2)
            ax.hlines(median, x - 0.17, x + 0.17, color="#222222", linewidth=1.4)
        row = results.loc[
            results["dataset"].eq(dataset) & results["comparison"].eq(comparison)
        ].iloc[0]
        ax.set_xticks(positions, [labels[group] for group in order])
        ax.set_title(title, loc="left", fontweight="bold")
        ax.set_ylabel("log2(deposited P2RX4 value + 1)")
        ax.text(
            0.02,
            0.98,
            (
                f"MW two-sided P={row['mann_whitney_p_two_sided']:.3g}\n"
                f"Cliff's delta={row['cliff_delta']:.2f}"
            ),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9.5,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.9},
        )
        ax.grid(axis="y", color="#D9D9D9", linewidth=0.6, alpha=0.7)
    fig.suptitle(
        "P2RX4 across separate human PDR-associated GEO datasets",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.5,
        -0.02,
        (
            "Points are deposited patient profiles; bars show medians and interquartile "
            "ranges. Panels represent different ocular compartments. All displayed "
            "Mann-Whitney P values are two-sided. GSE276892 is inconclusive after "
            "source/covariate sensitivities; GSE179568 and GSE94019 are post-protocol."
        ),
        ha="center",
        va="top",
        fontsize=9.2,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_gse179568_clinical_sensitivity(
    metadata: pd.DataFrame,
    clinical_results: pd.DataFrame,
    output: Path,
) -> None:
    """Show age separation and model dependence in the primary membrane contrast."""
    frame = metadata.loc[
        metadata["sample_id"].str.startswith(("PDR_", "Gliose_"))
    ].copy()
    frame["plot_group"] = np.where(
        frame["disease_group"].eq("PDR"), "PDR RNV", "Macular-pucker membrane"
    )
    palette = {"PDR RNV": "#B93B36", "Macular-pucker membrane": "#3D74AE"}
    row = clinical_results.loc[
        clinical_results["comparison"].eq("RNV vs macular-pucker membrane")
    ].iloc[0]

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.8), constrained_layout=True)
    ax = axes[0]
    for group in ["Macular-pucker membrane", "PDR RNV"]:
        subset = frame.loc[frame["plot_group"].eq(group)]
        ax.scatter(
            subset["age"],
            subset["log2_p2rx4_plus_1"],
            s=48,
            color=palette[group],
            edgecolor="white",
            linewidth=0.7,
            label=f"{group} (n={len(subset)})",
            zorder=3,
        )
    ax.axvspan(60, 62, color="#F1C75B", alpha=0.24, zorder=0)
    ax.text(
        61,
        ax.get_ylim()[1],
        "no observed age overlap",
        ha="center",
        va="top",
        fontsize=9,
        color="#765C14",
    )
    ax.set_xlabel("Age (years)")
    ax.set_ylabel("log2(deposited P2RX4 value + 1)")
    ax.set_title("A  Age and expression", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(color="#D9D9D9", linewidth=0.6, alpha=0.7)

    model_specs = [
        ("Unadjusted HC3", "unadjusted_ols_hc3_t"),
        ("Age-adjusted HC3", "age_adjusted_ols_hc3_t"),
        ("Age + sex HC3", "age_sex_adjusted_ols_hc3_t"),
        (
            "Age + sex + treatment HC3",
            "age_sex_treatment_adjusted_ols_hc3_t",
        ),
        ("Age + sex Huber", "age_sex_adjusted_huber_h1"),
    ]
    estimates = []
    lows = []
    highs = []
    labels = []
    for label, prefix in model_specs:
        estimates.append(float(row[f"{prefix}_beta"]))
        if prefix.endswith("huber_h1"):
            lows.append(float(row[f"{prefix}_ci_low_asymptotic"]))
            highs.append(float(row[f"{prefix}_ci_high_asymptotic"]))
        else:
            lows.append(float(row[f"{prefix}_ci_low"]))
            highs.append(float(row[f"{prefix}_ci_high"]))
        labels.append(label)
    positions = np.arange(len(labels))[::-1]
    ax = axes[1]
    estimates_array = np.asarray(estimates)
    ax.errorbar(
        estimates_array,
        positions,
        xerr=np.vstack(
            [estimates_array - np.asarray(lows), np.asarray(highs) - estimates_array]
        ),
        fmt="o",
        color="#363636",
        ecolor="#6A6A6A",
        capsize=3,
        markersize=6,
    )
    ax.axvline(0, color="#B93B36", linewidth=1.1, linestyle="--")
    ax.set_yticks(positions, labels)
    ax.set_xlabel("PDR minus control difference on log2(value + 1) scale")
    ax.set_title("B  Model-dependence sensitivity", loc="left", fontweight="bold")
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.6, alpha=0.7)
    fig.suptitle(
        "GSE179568 clinical-confounding sensitivity for P2RX4",
        fontsize=14,
        fontweight="bold",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if args.bootstrap < 100:
        raise ValueError("Use at least 100 bootstrap iterations")
    package = args.package_root.resolve()
    data_dir = package / "analysis_data" / "independent_validation"
    results_dir = package / "analysis_results"
    figures_dir = package / "figures"
    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    protocol = data_dir / "P2RX4_VALIDATION_PROTOCOL_LOCALLY_FROZEN.md"
    protocol_hash = sha256(protocol)
    gse276_meta, gse276_values, gse276_source = load_gse276892(data_dir)
    gse179_meta, gse179_values, gse179_source = load_gse179568(data_dir)
    gse940_meta, gse940_values, gse940_transcripts = load_gse94019(data_dir)
    gse179_clinical, gse179_balance = gse179568_clinical_sensitivity(gse179_meta)

    comparisons = [
        Comparison(
            "GSE276892",
            "PDR vs all surgical controls",
            "PDR vitreous hyalocytes",
            "macular-pucker/hole vitreous hyalocytes",
            tuple(row[0] for row in PDR_276892),
            tuple(row[0] for row in CONTROL_276892),
            "locally pre-specified patient-level sensitivity",
            "DESeq2-normalized reads; not raw counts",
        ),
        Comparison(
            "GSE276892",
            "PDR vs macular-pucker controls",
            "PDR vitreous hyalocytes",
            "macular-pucker vitreous hyalocytes",
            tuple(row[0] for row in PDR_276892),
            tuple(row[0] for row in CONTROL_276892 if row[0].startswith("MP")),
            "locally pre-specified descriptive subgroup sensitivity",
            "DESeq2-normalized reads; not raw counts",
        ),
        Comparison(
            "GSE276892",
            "PDR vs macular-hole controls",
            "PDR vitreous hyalocytes",
            "macular-hole vitreous hyalocytes",
            tuple(row[0] for row in PDR_276892),
            tuple(row[0] for row in CONTROL_276892 if row[0].startswith("MH")),
            "locally pre-specified descriptive subgroup sensitivity",
            "DESeq2-normalized reads; not raw counts",
        ),
        Comparison(
            "GSE276892",
            "PDR vs non-diabetic surgical controls",
            "PDR vitreous hyalocytes",
            "macular-pucker/hole vitreous hyalocytes excluding MP_S13",
            tuple(row[0] for row in PDR_276892),
            tuple(row[0] for row in CONTROL_276892 if row[0] != "MP_S13"),
            "additional diabetes-status sensitivity",
            "DESeq2-normalized reads; not raw counts",
        ),
        Comparison(
            "GSE179568",
            "RNV vs macular-pucker membrane",
            "PDR retinal neovascularization membrane",
            "macular-pucker epiretinal membrane",
            tuple(f"PDR_S{i}" for i in range(1, 8)),
            tuple(f"Gliose_S{i}" for i in range(1, 11)),
            "post-protocol membrane-compartment secondary contrast",
            "DESeq2-normalized reads; not raw counts",
        ),
        Comparison(
            "GSE179568",
            "RNV vs macular-hole ILM",
            "PDR retinal neovascularization membrane",
            "macular-hole inner limiting membrane",
            tuple(f"PDR_S{i}" for i in range(1, 8)),
            tuple(f"ILM_S{i}" for i in range(1, 8)),
            "post-protocol secondary sensitivity",
            "DESeq2-normalized reads; not raw counts",
        ),
        Comparison(
            "GSE179568",
            "RNV vs pooled membrane controls",
            "PDR retinal neovascularization membrane",
            "macular-pucker membrane and macular-hole ILM",
            tuple(f"PDR_S{i}" for i in range(1, 8)),
            tuple(
                [*(f"Gliose_S{i}" for i in range(1, 11)),
                 *(f"ILM_S{i}" for i in range(1, 8))]
            ),
            "post-protocol pooled-control sensitivity",
            "DESeq2-normalized reads; not raw counts",
        ),
        Comparison(
            "GSE94019",
            "PDR FVM CD31+ vs control retinal CD31+",
            "PDR fibrovascular-membrane CD31+ cells",
            "non-diabetic post-mortem retinal CD31+ cells",
            tuple(sample for sample in gse940_values if sample.startswith("MEEI")),
            tuple(sample for sample in gse940_values if not sample.startswith("MEEI")),
            "post-protocol secondary tissue-mismatched contrast",
            (
                "summed Partek E/M transcript abundance measurements; unit not "
                "specified as CPM or raw counts"
            ),
        ),
        Comparison(
            "GSE276892",
            "PDR vs newly generated surgical controls",
            "new GSE276892 PDR vitreous-hyalocyte profiles",
            "new GSE276892 macular-pucker/hole control profiles",
            tuple(row[0] for row in PDR_276892),
            tuple(sorted(NEW_CONTROL_IDS_276892)),
            "source-restriction sensitivity with only two new controls",
            "DESeq2-normalized reads; not raw counts",
        ),
    ]
    value_sets = {
        "GSE276892": gse276_values,
        "GSE179568": gse179_values,
        "GSE94019": gse940_values,
    }
    rng = np.random.default_rng(args.seed)
    result_rows: list[dict[str, object]] = []
    loo_rows: list[dict[str, object]] = []
    for comparison in comparisons:
        result, loo = compare_groups(
            comparison, value_sets[comparison.dataset], args.bootstrap, rng
        )
        result_rows.append(result)
        loo_rows.extend(loo)
    results = pd.DataFrame(result_rows)

    main_276 = results["comparison"].eq("PDR vs all surgical controls")
    for key, value in gse276_source.items():
        results.loc[main_276, key] = value
    for key, value in adjusted_gse276892(gse276_meta).items():
        results.loc[main_276, key] = value
    results.loc[main_276, "locally_specified_validation_status"] = (
        "inconclusive: processed values and the post-protocol raw-read reconstruction "
        "do not provide a precise disease effect; preprocessing was not prespecified"
    )

    pdr_ilm = results["comparison"].eq("RNV vs macular-hole ILM")
    results.loc[pdr_ilm, "source_reported_log2_fold_change"] = gse179_source[
        "source_reported_log2_fold_change_pdr_vs_ilm"
    ]
    results.loc[pdr_ilm, "source_reported_adjusted_p"] = gse179_source[
        "source_reported_adjusted_p_pdr_vs_ilm"
    ]
    clinical_columns = [
        column
        for column in gse179_clinical.columns
        if column not in {"dataset", "comparison"}
    ]
    results = results.merge(
        gse179_clinical[["dataset", "comparison", *clinical_columns]],
        on=["dataset", "comparison"],
        how="left",
        suffixes=("", "_clinical"),
        validate="one_to_one",
    )
    results = annotate_result_status(results)
    results = add_post_protocol_holm(results)

    sample_table = pd.concat(
        [gse276_meta, gse179_meta, gse940_meta], ignore_index=True, sort=False
    )
    preferred_columns = [
        "dataset",
        "sample_id",
        "sample_title",
        "geo_accession",
        "disease_group",
        "control_subtype",
        "diagnosis_detail",
        "age",
        "sex",
        "dm_type",
        "ocular_diagnosis",
        "supplementary_table_row",
        "lens_status",
        "previous_vitrectomy",
        "previous_anti_vegf_over_3_months",
        "previous_prp",
        "pseudophakic",
        "clinical_source",
        "source_dataset",
        "reused_control",
        "cell_count",
        "rna_concentration_pg_per_ul",
        "tissue_or_cell_fraction",
        "deposited_value_scale",
        "p2rx4_value",
        "log2_p2rx4_plus_1",
        "protocol_role",
    ]
    sample_table = sample_table.reindex(columns=preferred_columns)
    eligibility = dataset_eligibility()
    overlap_audit = patient_overlap_audit()
    deviations = protocol_deviations(protocol_hash, gse940_transcripts)
    qc_correlations = gse276892_qc_correlations(gse276_meta)
    loo_table = pd.DataFrame(loo_rows)

    eligibility.to_csv(
        results_dir / "Independent_validation_dataset_eligibility.csv", index=False
    )
    overlap_audit.to_csv(
        results_dir / "Independent_validation_patient_overlap_audit.csv", index=False
    )
    sample_table.to_csv(
        results_dir / "Independent_validation_P2RX4_sample_level.csv", index=False
    )
    results.to_csv(
        results_dir / "Independent_validation_P2RX4_results.csv", index=False
    )
    deviations.to_csv(
        results_dir / "Independent_validation_protocol_deviations.csv", index=False
    )
    qc_correlations.to_csv(
        results_dir / "Independent_validation_GSE276892_qc_correlations.csv",
        index=False,
    )
    loo_table.to_csv(
        results_dir / "Independent_validation_P2RX4_leave_one_out.csv", index=False
    )
    gse179_clinical.to_csv(
        results_dir / "Independent_validation_GSE179568_clinical_sensitivity.csv",
        index=False,
    )
    gse179_balance.to_csv(
        results_dir / "Independent_validation_GSE179568_covariate_balance.csv",
        index=False,
    )

    summary = {
        "target_gene": TARGET_GENE,
        "protocol_sha256": protocol_hash,
        "bootstrap_iterations": args.bootstrap,
        "random_seed": args.seed,
        "locally_specified_primary_status": (
            "Inconclusive. The exact source count matrix was not deposited; the "
            "post-protocol raw-read reconstruction is available but depends on "
            "preprocessing choices absent from the local protocol."
        ),
        "source_reported_gse276892": gse276_source,
        "gse276892_patient_level_direction": results.loc[
            main_276,
            [
                "mean_difference",
                "mean_log2_value_difference",
                "cliff_delta",
                "mann_whitney_p_one_sided_greater",
                "mann_whitney_p_two_sided",
                "loo_mean_direction_positive_fraction",
            ],
        ].iloc[0].to_dict(),
        "source_reported_gse179568": gse179_source,
        "gse179568_clinical_sensitivity": {
            "table": "Independent_validation_GSE179568_clinical_sensitivity.csv",
            "covariate_balance": (
                "Independent_validation_GSE179568_covariate_balance.csv"
            ),
            "age_common_support_in_all_three_comparisons": bool(
                gse179_clinical["age_common_support"].all()
            ),
            "interpretation": (
                "No age common support exists between PDR and either control "
                "compartment. HC3 and Huber models are extrapolative sensitivities, "
                "not deconfounded disease-effect estimates."
            ),
        },
        "external_dataset_search": {
            "date": SEARCH_DATE,
            "databases": SEARCH_DATABASES,
            "search_string": SEARCH_STRING,
            "screening_table": "Independent_validation_dataset_eligibility.csv",
            "status": "retrospectively reconstructed; not a contemporaneous registry",
        },
        "patient_overlap_audit": {
            "table": "Independent_validation_patient_overlap_audit.csv",
            "conclusion": (
                "Public metadata do not prove patient non-overlap between the two "
                "Freiburg datasets; the manuscript therefore says separate GEO "
                "datasets rather than independent patients."
            ),
        },
        "scope_statement": (
            "The results describe P2RX4 in separate human PDR-associated GEO datasets. "
            "They do not prove patient non-overlap, whole-retina replication, "
            "causality, or cell-type specificity. GSE276892 is inconclusive after "
            "source/covariate sensitivity analyses."
        ),
        "input_sha256": {
            str(path.relative_to(package)): sha256(path)
            for path in [
                protocol,
                data_dir / "GSE276892_normal_data.csv.gz",
                data_dir / "GSE276892_README_all_samples.xlsx",
                data_dir / "GSE276892_primary_article_table1.html",
                data_dir / "GSE94019_Partek_EM_gene_reads.txt.gz",
                data_dir / "GSE179568" / "GSE179568_data.csv.gz",
                data_dir / "GSE179568" / "GSE179568_family.soft.gz",
                data_dir / "GSE179568" / "GSE179568_series_matrix.txt.gz",
                data_dir / "GSE179568" / "Table 1.pdf",
            ]
        },
    }
    (results_dir / "independent_validation_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    plot_validation(
        sample_table,
        results,
        figures_dir / "Figure_4_independent_P2RX4_validation",
    )
    plot_gse179568_clinical_sensitivity(
        gse179_meta,
        gse179_clinical,
        figures_dir / "Supplementary_Figure_GSE179568_clinical_sensitivity",
    )

    if not math.isclose(
        float(results.loc[main_276, "mann_whitney_p_two_sided"].iloc[0]),
        2
        * float(
            results.loc[
                main_276, "mann_whitney_p_one_sided_greater"
            ].iloc[0]
        ),
        rel_tol=1e-12,
    ):
        raise AssertionError("Unexpected one-/two-sided MW relation in GSE276892")
    print(
        json.dumps(
            {
                "target": TARGET_GENE,
                "samples": len(sample_table),
                "comparisons": len(results),
                "protocol_sha256": protocol_hash,
                "figure_png": str(
                    figures_dir / "Figure_4_independent_P2RX4_validation.png"
                ),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
