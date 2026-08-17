#!/usr/bin/env python3
"""Final donor-level discovery and robustness analyses for the Research Article.

The primary estimand is the total association between donor-level DR severity
and standardized macular expression among donors with diabetes.  A separate
DME-conditioned model is reported because the causal position of DME cannot be
identified from this cross-sectional public dataset.  All small-sample HC3
tests use residual degrees of freedom, and a studentized wild bootstrap is
reported as a heteroskedasticity-aware sensitivity analysis.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import mannwhitneyu
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, log_loss, roc_auc_score
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests


DEFAULT_SEED = 20260813
GROUP_ORDER = ["Diabetic", "NPDR", "NPDR + DME", "PDR + DME"]
BASE_COVARIATES = ["age", "sex_male", "pmi", "rin"]
STAGE_ORDER = {
    "Diabetic": 0,
    "NPDR": 1,
    "NPDR + DME": 1,
    "PDR + DME": 2,
}


def parse_args() -> argparse.Namespace:
    package = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=package.parents[1])
    parser.add_argument("--package-root", type=Path, default=package)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--wild-bootstrap", type=int, default=4999)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    return parser.parse_args()


def zscore(values: pd.Series | np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    scale = float(array.std(ddof=0))
    if not np.isfinite(scale) or scale == 0:
        return np.zeros_like(array)
    return (array - float(array.mean())) / scale


def bh(values: pd.Series | np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    output = np.full(array.shape, np.nan, dtype=float)
    valid = np.isfinite(array)
    if valid.any():
        output[valid] = multipletests(array[valid], method="fdr_bh")[1]
    return output


def cliff_delta(case: np.ndarray, control: np.ndarray) -> float:
    return float(np.sign(case[:, None] - control[None, :]).mean())


def read_hallmark(path: Path) -> list[str]:
    fields = path.read_text(encoding="utf-8").strip().split("\t")
    return list(dict.fromkeys(value.upper() for value in fields[2:] if value))


def load_donor_data(project: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    processed = project / "data_processed"
    matrix = pd.read_csv(processed / "log2cpm_macula_4groups.tsv", sep="\t")
    matrix["ensemblID"] = matrix["ensemblID"].astype(str).str.split(".").str[0]
    mapping = pd.read_csv(processed / "ensembl_to_symbol_mapping.csv")
    mapping["ensembl_id"] = mapping["ensembl_id"].astype(str).str.split(".").str[0]
    mapping["gene_symbol"] = mapping["gene_symbol"].astype(str).str.upper()
    matrix = matrix.merge(mapping, left_on="ensemblID", right_on="ensembl_id", how="inner")
    sample_columns = [column for column in matrix if column.startswith("sample_")]
    matrix["row_mean"] = matrix[sample_columns].mean(axis=1)
    matrix = matrix.sort_values("row_mean", ascending=False).drop_duplicates("gene_symbol")
    by_sample = matrix.set_index("gene_symbol")[sample_columns].T

    manifest = pd.read_csv(processed / "manifest_macula_4groups.csv")
    joined = manifest.set_index("sample_id").join(by_sample, how="inner")
    genes = by_sample.columns.tolist()
    donor_expression = joined.groupby("donor", sort=False)[genes].mean()
    donor_metadata = (
        manifest.groupby("donor", sort=False)
        .agg(
            disease_group=("disease_group", "first"),
            detailed_group=("disease_group_detailed", "first"),
            severity_mean=("dr_severity_score", "mean"),
            severity_worst_eye=("dr_severity_score", "max"),
            severity_best_eye=("dr_severity_score", "min"),
            severity_eye_values=("dr_severity_score", lambda x: ";".join(map(str, x))),
            dme=("dme", "max"),
            age=("age", "mean"),
            sex=("sex", "first"),
            pmi=("post_mortem_interval_min", "mean"),
            rin=("rin", "mean"),
            n_eyes=("sample_id", "size"),
            sample_ids=("sample_id", lambda x: ";".join(map(str, x))),
        )
        .loc[donor_expression.index]
    )
    donor_metadata["sex_male"] = donor_metadata["sex"].eq("male").astype(int)
    donor_metadata["stage_ordinal"] = donor_metadata["detailed_group"].map(STAGE_ORDER)
    donor_metadata["severity_rank"] = donor_metadata["severity_mean"].rank(method="average")
    return donor_expression, donor_metadata


def design_matrix(
    metadata: pd.DataFrame,
    *,
    include_severity: bool = True,
    include_dme: bool = False,
    severity_column: str = "severity_mean",
) -> pd.DataFrame:
    data: dict[str, np.ndarray] = {}
    if include_severity:
        data["severity"] = zscore(metadata[severity_column])
    if include_dme:
        data["dme"] = metadata["dme"].astype(float).to_numpy()
    data.update(
        {
            "age": zscore(metadata["age"]),
            "sex_male": metadata["sex_male"].astype(float).to_numpy(),
            "pmi": zscore(metadata["pmi"]),
            "rin": zscore(metadata["rin"]),
        }
    )
    return sm.add_constant(pd.DataFrame(data, index=metadata.index), has_constant="add")


def standardized_outcomes(expression: pd.DataFrame, genes: list[str]) -> np.ndarray:
    return np.column_stack([zscore(expression[gene]) for gene in genes])


def fit_gene_models(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    genes: list[str],
    *,
    model_name: str,
    include_dme: bool = False,
    severity_column: str = "severity_mean",
) -> pd.DataFrame:
    design = design_matrix(
        metadata,
        include_dme=include_dme,
        severity_column=severity_column,
    )
    rows: list[dict[str, float | int | str]] = []
    for gene in genes:
        outcome = pd.Series(zscore(expression[gene]), index=metadata.index)
        fit = sm.OLS(outcome, design).fit(cov_type="HC3", use_t=True)
        interval = fit.conf_int().loc["severity"]
        rows.append(
            {
                "model": model_name,
                "gene_symbol": gene,
                "severity_beta": float(fit.params["severity"]),
                "severity_se_hc3": float(fit.bse["severity"]),
                "severity_ci_low_t": float(interval.iloc[0]),
                "severity_ci_high_t": float(interval.iloc[1]),
                "severity_t": float(fit.tvalues["severity"]),
                "severity_pvalue_t": float(fit.pvalues["severity"]),
                "residual_df": int(fit.df_resid),
                "n_donors": len(metadata),
                "dme_beta": float(fit.params["dme"]) if include_dme else np.nan,
                "dme_pvalue_t": float(fit.pvalues["dme"]) if include_dme else np.nan,
            }
        )
    table = pd.DataFrame(rows)
    table["severity_padj_bh_158"] = bh(table["severity_pvalue_t"])
    table["rank_descending_beta"] = table["severity_beta"].rank(
        method="first", ascending=False
    ).astype(int)
    return table.sort_values("rank_descending_beta").reset_index(drop=True)


def hc3_t_matrix(x: np.ndarray, y: np.ndarray, coefficient: int = 1) -> tuple[np.ndarray, np.ndarray]:
    inverse = np.linalg.pinv(x)
    beta = inverse @ y
    residual = y - x @ beta
    leverage = np.sum(x * inverse.T, axis=1)
    adjusted = residual / np.clip(1.0 - leverage[:, None], 1e-8, None)
    weights = inverse[coefficient, :, None]
    variance = np.sum((weights * adjusted) ** 2, axis=0)
    standard_error = np.sqrt(np.maximum(variance, 0.0))
    statistic = np.divide(
        beta[coefficient],
        standard_error,
        out=np.zeros_like(beta[coefficient]),
        where=standard_error > 0,
    )
    return beta[coefficient], statistic


def wild_bootstrap_t(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    genes: list[str],
    *,
    include_dme: bool,
    iterations: int,
    seed: int,
) -> pd.DataFrame:
    full = design_matrix(metadata, include_dme=include_dme).to_numpy(dtype=float)
    reduced = design_matrix(
        metadata, include_severity=False, include_dme=include_dme
    ).to_numpy(dtype=float)
    outcomes = standardized_outcomes(expression, genes)
    observed_beta, observed_t = hc3_t_matrix(full, outcomes)

    reduced_inverse = np.linalg.pinv(reduced)
    fitted = reduced @ (reduced_inverse @ outcomes)
    residual = outcomes - fitted
    leverage = np.sum(reduced * reduced_inverse.T, axis=1)
    adjusted_residual = residual / np.clip(1.0 - leverage[:, None], 1e-8, None)
    rng = np.random.default_rng(seed)
    exceed = np.zeros(len(genes), dtype=int)
    for _ in range(iterations):
        multipliers = rng.choice(np.array([-1.0, 1.0]), size=(len(metadata), 1))
        bootstrap_outcome = fitted + adjusted_residual * multipliers
        _, statistic = hc3_t_matrix(full, bootstrap_outcome)
        exceed += np.abs(statistic) >= np.abs(observed_t)
    pvalues = (exceed + 1) / (iterations + 1)
    return pd.DataFrame(
        {
            "gene_symbol": genes,
            "observed_beta": observed_beta,
            "observed_hc3_t": observed_t,
            "wild_bootstrap_t_two_sided_p": pvalues,
            "wild_bootstrap_t_padj_bh_158": bh(pvalues),
            "wild_bootstrap_iterations": iterations,
            "wild_bootstrap_weight": "Rademacher",
            "null_residual_scaling": "HC3 (residual/(1-leverage))",
        }
    )


def coefficient_vector(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    genes: list[str],
    *,
    include_dme: bool = False,
    severity_column: str = "severity_mean",
) -> np.ndarray:
    design = design_matrix(
        metadata,
        include_dme=include_dme,
        severity_column=severity_column,
    ).to_numpy(dtype=float)
    return (np.linalg.pinv(design) @ standardized_outcomes(expression, genes))[1]


def rank_array(coefficients: np.ndarray) -> np.ndarray:
    order = np.argsort(-coefficients, kind="stable")
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, len(order) + 1)
    return ranks


def bootstrap_stability(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    genes: list[str],
    iterations: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    groups = [
        np.flatnonzero(metadata["detailed_group"].to_numpy() == group)
        for group in GROUP_ORDER
    ]
    coefficients = np.empty((iterations, len(genes)), dtype=np.float32)
    ranks = np.empty((iterations, len(genes)), dtype=np.int16)
    for iteration in range(iterations):
        sampled = np.concatenate(
            [rng.choice(index, len(index), replace=True) for index in groups]
        )
        sampled_expression = expression.iloc[sampled].reset_index(drop=True)
        sampled_metadata = metadata.iloc[sampled].reset_index(drop=True)
        beta = coefficient_vector(sampled_expression, sampled_metadata, genes)
        coefficients[iteration] = beta
        ranks[iteration] = rank_array(beta)
    return pd.DataFrame(
        {
            "gene_symbol": genes,
            "bootstrap_top5_frequency": (ranks <= 5).mean(axis=0),
            "bootstrap_top10_frequency": (ranks <= 10).mean(axis=0),
            "bootstrap_top20_frequency": (ranks <= 20).mean(axis=0),
            "bootstrap_positive_frequency": (coefficients > 0).mean(axis=0),
            "bootstrap_beta_median": np.median(coefficients, axis=0),
            "bootstrap_beta_ci_low": np.quantile(coefficients, 0.025, axis=0),
            "bootstrap_beta_ci_high": np.quantile(coefficients, 0.975, axis=0),
            "bootstrap_rank_median": np.median(ranks, axis=0),
            "bootstrap_rank_ci_low": np.quantile(ranks, 0.025, axis=0),
            "bootstrap_rank_ci_high": np.quantile(ranks, 0.975, axis=0),
            "bootstrap_iterations": iterations,
        }
    )


def lodo_stability(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    genes: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    coefficients = np.empty((len(metadata), len(genes)), dtype=float)
    ranks = np.empty((len(metadata), len(genes)), dtype=int)
    detailed_rows: list[dict[str, float | str]] = []
    for donor_index, held_out in enumerate(metadata.index):
        keep = metadata.index.difference([held_out], sort=False)
        beta = coefficient_vector(expression.loc[keep], metadata.loc[keep], genes)
        donor_ranks = rank_array(beta)
        coefficients[donor_index] = beta
        ranks[donor_index] = donor_ranks
        for gene_index, gene in enumerate(genes):
            detailed_rows.append(
                {
                    "held_out_donor": held_out,
                    "gene_symbol": gene,
                    "severity_beta": float(beta[gene_index]),
                    "rank_descending_beta": int(donor_ranks[gene_index]),
                }
            )
    summary = pd.DataFrame(
        {
            "gene_symbol": genes,
            "lodo_top5_frequency": (ranks <= 5).mean(axis=0),
            "lodo_top10_frequency": (ranks <= 10).mean(axis=0),
            "lodo_top20_frequency": (ranks <= 20).mean(axis=0),
            "lodo_positive_frequency": (coefficients > 0).mean(axis=0),
            "lodo_beta_min": coefficients.min(axis=0),
            "lodo_beta_max": coefficients.max(axis=0),
            "lodo_rank_min": ranks.min(axis=0),
            "lodo_rank_median": np.median(ranks, axis=0),
            "lodo_rank_max": ranks.max(axis=0),
        }
    )
    return summary, pd.DataFrame(detailed_rows)


def influence_diagnostics(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    genes: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    design = design_matrix(metadata)
    leverage = sm.OLS(np.zeros(len(metadata)), design).fit().get_influence().hat_matrix_diag
    threshold = 2.0 * design.shape[1] / len(metadata)
    detail: list[dict[str, float | str | bool]] = []
    summary: list[dict[str, float | str]] = []
    for gene in genes:
        outcome = pd.Series(zscore(expression[gene]), index=metadata.index)
        fit = sm.OLS(outcome, design).fit()
        influence = fit.get_influence()
        cooks = influence.cooks_distance[0]
        dfbeta = influence.dfbetas[:, list(design.columns).index("severity")]
        max_index = int(np.nanargmax(cooks))
        summary.append(
            {
                "gene_symbol": gene,
                "maximum_cooks_distance": float(cooks[max_index]),
                "maximum_cooks_donor": str(metadata.index[max_index]),
                "maximum_absolute_severity_dfbeta": float(np.nanmax(np.abs(dfbeta))),
            }
        )
        for index, donor in enumerate(metadata.index):
            detail.append(
                {
                    "gene_symbol": gene,
                    "donor": donor,
                    "leverage": float(leverage[index]),
                    "high_leverage_threshold_2p_over_n": float(threshold),
                    "is_high_leverage": bool(leverage[index] > threshold),
                    "cooks_distance": float(cooks[index]),
                    "severity_dfbeta": float(dfbeta[index]),
                    "pmi_minutes": float(metadata.loc[donor, "pmi"]),
                }
            )
    return pd.DataFrame(summary), pd.DataFrame(detail)


def model_sensitivities(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    genes: list[str],
    high_leverage_donors: list[str],
) -> pd.DataFrame:
    specifications = [
        ("total_mean_severity", False, "severity_mean", metadata.index),
        ("DME_conditioned_mean_severity", True, "severity_mean", metadata.index),
        ("total_worst_eye_severity", False, "severity_worst_eye", metadata.index),
        ("total_rank_transformed_severity", False, "severity_rank", metadata.index),
        ("total_clinical_stage_ordinal", False, "stage_ordinal", metadata.index),
        (
            "total_excluding_high_leverage",
            False,
            "severity_mean",
            metadata.index.difference(high_leverage_donors, sort=False),
        ),
    ]
    frames = []
    for name, include_dme, severity_column, donors in specifications:
        frames.append(
            fit_gene_models(
                expression.loc[donors],
                metadata.loc[donors],
                genes,
                model_name=name,
                include_dme=include_dme,
                severity_column=severity_column,
            )
        )

    robust_design = design_matrix(metadata)
    robust_rows = []
    for gene in genes:
        outcome = pd.Series(zscore(expression[gene]), index=metadata.index)
        fit = sm.RLM(outcome, robust_design, M=sm.robust.norms.HuberT()).fit()
        robust_rows.append(
            {
                "model": "total_Huber_RLM",
                "gene_symbol": gene,
                "severity_beta": float(fit.params["severity"]),
                "severity_se_hc3": np.nan,
                "severity_ci_low_t": np.nan,
                "severity_ci_high_t": np.nan,
                "severity_t": np.nan,
                "severity_pvalue_t": np.nan,
                "residual_df": np.nan,
                "n_donors": len(metadata),
                "dme_beta": np.nan,
                "dme_pvalue_t": np.nan,
                "severity_padj_bh_158": np.nan,
            }
        )
    robust = pd.DataFrame(robust_rows)
    robust["rank_descending_beta"] = robust["severity_beta"].rank(
        method="first", ascending=False
    ).astype(int)
    frames.append(robust.sort_values("rank_descending_beta"))
    return pd.concat(frames, ignore_index=True)


def dme_interaction_p2rx4(
    expression: pd.DataFrame, metadata: pd.DataFrame
) -> pd.DataFrame:
    design = design_matrix(metadata, include_dme=True)
    design.insert(2, "severity_by_dme", design["severity"] * design["dme"])
    outcome = pd.Series(zscore(expression["P2RX4"]), index=metadata.index)
    fit = sm.OLS(outcome, design).fit(cov_type="HC3", use_t=True)
    rows = []
    for term in ["severity", "dme", "severity_by_dme"]:
        interval = fit.conf_int().loc[term]
        rows.append(
            {
                "gene_symbol": "P2RX4",
                "model": "exploratory_DME_interaction",
                "term": term,
                "beta": float(fit.params[term]),
                "se_hc3": float(fit.bse[term]),
                "ci_low_t": float(interval.iloc[0]),
                "ci_high_t": float(interval.iloc[1]),
                "pvalue_t": float(fit.pvalues[term]),
                "residual_df": int(fit.df_resid),
            }
        )
    return pd.DataFrame(rows)


def clinical_contrasts(
    expression: pd.DataFrame, metadata: pd.DataFrame, genes: list[str]
) -> pd.DataFrame:
    definitions = [
        ("NPDR_vs_diabetes", ["NPDR"], ["Diabetic"], "primary_contrast_family"),
        (
            "advanced_DME_vs_diabetes",
            ["NPDR + DME", "PDR + DME"],
            ["Diabetic"],
            "exploratory_contrast_family",
        ),
        (
            "NPDR_DME_vs_diabetes",
            ["NPDR + DME"],
            ["Diabetic"],
            "exploratory_contrast_family",
        ),
        (
            "PDR_DME_vs_diabetes",
            ["PDR + DME"],
            ["Diabetic"],
            "exploratory_contrast_family",
        ),
        (
            "PDR_DME_vs_NPDR_DME",
            ["PDR + DME"],
            ["NPDR + DME"],
            "exploratory_contrast_family",
        ),
    ]
    rows: list[dict[str, float | int | str | bool]] = []
    for name, case_groups, control_groups, family in definitions:
        case_mask = metadata["detailed_group"].isin(case_groups)
        control_mask = metadata["detailed_group"].isin(control_groups)
        subset = case_mask | control_mask
        subset_metadata = metadata.loc[subset].copy()
        subset_metadata["case"] = case_mask.loc[subset].astype(int)
        estimable = int(case_mask.sum()) >= 3 and int(control_mask.sum()) >= 3 and int(subset.sum()) > 7
        if estimable:
            adjusted_design = pd.DataFrame(
                {
                    "case": subset_metadata["case"].astype(float),
                    "age": zscore(subset_metadata["age"]),
                    "sex_male": subset_metadata["sex_male"].astype(float),
                    "pmi": zscore(subset_metadata["pmi"]),
                    "rin": zscore(subset_metadata["rin"]),
                },
                index=subset_metadata.index,
            )
            adjusted_design = sm.add_constant(adjusted_design, has_constant="add")
        for gene in genes:
            case = expression.loc[case_mask, gene].to_numpy(dtype=float)
            control = expression.loc[control_mask, gene].to_numpy(dtype=float)
            test = mannwhitneyu(case, control, alternative="two-sided")
            adjusted_beta = np.nan
            adjusted_pvalue = np.nan
            adjusted_df = np.nan
            if estimable:
                outcome = pd.Series(zscore(expression.loc[subset, gene]), index=subset_metadata.index)
                fit = sm.OLS(outcome, adjusted_design).fit(cov_type="HC3", use_t=True)
                adjusted_beta = float(fit.params["case"])
                adjusted_pvalue = float(fit.pvalues["case"])
                adjusted_df = int(fit.df_resid)
            rows.append(
                {
                    "contrast": name,
                    "multiplicity_family": family,
                    "gene_symbol": gene,
                    "case_group": ";".join(case_groups),
                    "control_group": ";".join(control_groups),
                    "n_case_donors": len(case),
                    "n_control_donors": len(control),
                    "mean_difference": float(case.mean() - control.mean()),
                    "cliff_delta": cliff_delta(case, control),
                    "mannwhitney_pvalue": float(test.pvalue),
                    "covariate_adjusted_case_beta": adjusted_beta,
                    "covariate_adjusted_case_pvalue_t": adjusted_pvalue,
                    "covariate_adjusted_residual_df": adjusted_df,
                    "covariate_adjusted_estimable": estimable,
                }
            )
    table = pd.DataFrame(rows)
    table["mannwhitney_padj_within_contrast_158"] = table.groupby("contrast")[
        "mannwhitney_pvalue"
    ].transform(bh)
    table["mannwhitney_padj_global_790"] = bh(table["mannwhitney_pvalue"])
    table["adjusted_padj_within_contrast_158"] = table.groupby("contrast")[
        "covariate_adjusted_case_pvalue_t"
    ].transform(bh)
    table["adjusted_padj_global_estimable_contrasts"] = bh(
        table["covariate_adjusted_case_pvalue_t"]
    )
    return table


def parse_geo_metadata(path: Path) -> pd.DataFrame:
    titles: list[str] = []
    accessions: list[str] = []
    characteristics: list[list[str]] = []
    sources: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = re.findall(r'"([^"]*)"', line)
            elif line.startswith("!Sample_geo_accession"):
                accessions = re.findall(r'"([^"]*)"', line)
            elif line.startswith("!Sample_source_name_ch1"):
                sources = re.findall(r'"([^"]*)"', line)
            elif line.startswith("!Sample_characteristics_ch1"):
                characteristics.append(re.findall(r'"([^"]*)"', line))
            elif line.startswith("!series_matrix_table_begin"):
                break
    rows = []
    for index, title in enumerate(titles):
        attributes: dict[str, str] = {}
        for item in characteristics:
            if index < len(item) and ":" in item[index]:
                key, value = item[index].split(":", 1)
                attributes[key.strip().lower()] = value.strip()
        rows.append(
            {
                "sample_title": title,
                "geo_accession": accessions[index],
                "source_name": sources[index],
                **attributes,
            }
        )
    table = pd.DataFrame(rows)
    table["analysis_group"] = np.select(
        [
            table["disease"].str.contains("diabetes", case=False),
            table["disease"].str.contains("branch retinal vein", case=False),
            table["disease"].str.contains("periphlebitis", case=False),
            table["disease"].str.contains("normal", case=False),
        ],
        ["diabetic_samples", "BRVO_samples", "retinal_periphlebitis_samples", "normal_retina"],
        default="other",
    )
    table["is_membrane_source"] = table["source_name"].str.contains("membrane", case=False, na=False)
    table["patient_identifier_available"] = False
    return table


def gse102485_membrane_sensitivity(
    project: Path,
    genes: list[str],
    metadata: pd.DataFrame,
    seed: int,
) -> pd.DataFrame:
    matrix = pd.read_csv(
        project / "GSE102485_expressed_gene_FPKM.txt.gz", sep="\t", compression="gzip"
    ).rename(columns={"Symbol": "gene_symbol"})
    matrix["gene_symbol"] = matrix["gene_symbol"].astype(str).str.upper()
    samples = [sample for sample in metadata["sample_title"] if sample in matrix]
    matrix = matrix[["gene_symbol", *samples]].copy()
    matrix["row_mean"] = matrix[samples].mean(axis=1)
    matrix = matrix.sort_values("row_mean", ascending=False).drop_duplicates("gene_symbol").set_index("gene_symbol")
    cases = metadata.loc[
        metadata["analysis_group"].eq("diabetic_samples") & metadata["is_membrane_source"],
        "sample_title",
    ].tolist()
    controls = metadata.loc[
        metadata["analysis_group"].isin(["BRVO_samples", "retinal_periphlebitis_samples"])
        & metadata["is_membrane_source"],
        "sample_title",
    ].tolist()
    rng = np.random.default_rng(seed)
    rows = []
    for gene in genes:
        case = np.log2(matrix.loc[gene, cases].astype(float).to_numpy() + 1)
        control = np.log2(matrix.loc[gene, controls].astype(float).to_numpy() + 1)
        bootstrap = np.empty(2000)
        for index in range(len(bootstrap)):
            bootstrap[index] = rng.choice(case, len(case), replace=True).mean() - rng.choice(
                control, len(control), replace=True
            ).mean()
        test = mannwhitneyu(case, control, alternative="two-sided")
        rows.append(
            {
                "gene_symbol": gene,
                "n_diabetic_membrane_profiles": len(case),
                "n_disease_control_membrane_profiles": len(control),
                "control_composition": "3 BRVO; 2 retinal periphlebitis",
                "patient_mapping_available": False,
                "mean_log2_fpkm_difference": float(case.mean() - control.mean()),
                "bootstrap_ci_low": float(np.quantile(bootstrap, 0.025)),
                "bootstrap_ci_high": float(np.quantile(bootstrap, 0.975)),
                "mannwhitney_pvalue": float(test.pvalue),
                "role": "supplementary post-hoc membrane-only context; not independent validation",
            }
        )
    table = pd.DataFrame(rows)
    table["mannwhitney_padj_bh"] = bh(table["mannwhitney_pvalue"])
    return table


def nested_predictive_audit(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    genes: list[str],
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, float | int | str]] = []
    covariates = BASE_COVARIATES
    for fold, held_out in enumerate(metadata.index):
        train = metadata.index.difference([held_out], sort=False)
        selected_beta = coefficient_vector(expression.loc[train], metadata.loc[train], genes)
        selected = [genes[index] for index in np.argsort(-selected_beta, kind="stable")[:5]]

        train_cov = metadata.loc[train, covariates].astype(float).copy()
        test_cov = metadata.loc[[held_out], covariates].astype(float).copy()
        for column in ["age", "pmi", "rin"]:
            mean = float(train_cov[column].mean())
            scale = float(train_cov[column].std(ddof=0)) or 1.0
            train_cov[column] = (train_cov[column] - mean) / scale
            test_cov[column] = (test_cov[column] - mean) / scale
        cov_train = sm.add_constant(train_cov, has_constant="add").to_numpy()
        cov_test = sm.add_constant(test_cov, has_constant="add").to_numpy()
        train_values = expression.loc[train, selected].to_numpy(dtype=float)
        residual_fit = np.linalg.pinv(cov_train) @ train_values
        x_train = train_values - cov_train @ residual_fit
        x_test = expression.loc[[held_out], selected].to_numpy(dtype=float) - cov_test @ residual_fit
        scaler = StandardScaler().fit(x_train)
        target = metadata.loc[train, "detailed_group"].ne("Diabetic").astype(int)
        classifier = LogisticRegression(C=1.0, solver="liblinear", random_state=seed + fold, max_iter=5000)
        classifier.fit(scaler.transform(x_train), target)
        gene_probability = float(classifier.predict_proba(scaler.transform(x_test))[:, 1][0])

        cov_scaler = StandardScaler().fit(metadata.loc[train, covariates].astype(float))
        cov_classifier = LogisticRegression(C=1.0, solver="liblinear", random_state=seed + 1000 + fold, max_iter=5000)
        cov_classifier.fit(cov_scaler.transform(metadata.loc[train, covariates].astype(float)), target)
        cov_probability = float(cov_classifier.predict_proba(cov_scaler.transform(metadata.loc[[held_out], covariates].astype(float)))[:, 1][0])
        label = int(metadata.loc[held_out, "detailed_group"] != "Diabetic")
        for method, probability, selected_genes in [
            ("nested_total_association_top5", gene_probability, ";".join(selected)),
            ("covariate_only", cov_probability, ""),
        ]:
            rows.append(
                {
                    "method": method,
                    "held_out_donor": held_out,
                    "true_label_any_DR": label,
                    "predicted_probability": probability,
                    "selected_genes": selected_genes,
                }
            )
    predictions = pd.DataFrame(rows)
    summaries = []
    for method, group in predictions.groupby("method", sort=False):
        y = group["true_label_any_DR"].to_numpy()
        p = group["predicted_probability"].to_numpy()
        summaries.append(
            {
                "method": method,
                "n_held_out_donors": len(group),
                "auc": float(roc_auc_score(y, p)),
                "brier_score": float(brier_score_loss(y, p)),
                "log_loss": float(log_loss(y, p, labels=[0, 1])),
                "interpretation": "internal falsification audit; not a clinical performance estimate",
            }
        )
    return predictions, pd.DataFrame(summaries)


def draw_discovery_figure(
    ranking: pd.DataFrame,
    candidate_genes: list[str],
    output: Path,
) -> None:
    total = ranking.loc[ranking["model"].eq("total_association")].set_index("gene_symbol")
    conditioned = ranking.loc[ranking["model"].eq("DME_conditioned")].set_index("gene_symbol")
    ordered = total.head(15).sort_values("severity_beta")
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.2), constrained_layout=True)
    colors = ["#b84235" if gene == "P2RX4" else "#6f98ad" for gene in ordered.index]
    axes[0].barh(ordered.index, ordered["severity_beta"], color=colors)
    axes[0].errorbar(
        ordered["severity_beta"],
        np.arange(len(ordered)),
        xerr=[
            ordered["severity_beta"] - ordered["severity_ci_low_t"],
            ordered["severity_ci_high_t"] - ordered["severity_beta"],
        ],
        fmt="none",
        ecolor="#333333",
        capsize=2,
    )
    axes[0].axvline(0, color="black", linewidth=0.8)
    axes[0].set_xlabel("Total-association severity coefficient (HC3+t 95% CI)")
    axes[0].set_title("a  Primary donor-level ranking", loc="left", fontweight="bold")
    axes[0].grid(axis="x", alpha=0.2)

    y = np.arange(len(candidate_genes))
    axes[1].scatter(
        total.loc[candidate_genes, "severity_beta"], y - 0.12, color="#2f6f8f", label="Total association"
    )
    axes[1].scatter(
        conditioned.loc[candidate_genes, "severity_beta"], y + 0.12, color="#d67b5d", label="DME-conditioned"
    )
    for index, gene in enumerate(candidate_genes):
        axes[1].plot(
            [total.loc[gene, "severity_beta"], conditioned.loc[gene, "severity_beta"]],
            [index - 0.12, index + 0.12],
            color="#aaaaaa",
            linewidth=0.8,
        )
    axes[1].axvline(0, color="black", linewidth=0.8)
    axes[1].set_yticks(y, candidate_genes)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Standardized severity coefficient")
    axes[1].set_title("b  Estimand sensitivity", loc="left", fontweight="bold")
    axes[1].legend(frameon=False)
    axes[1].grid(axis="x", alpha=0.2)
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def draw_stability_figure(
    final: pd.DataFrame,
    candidate_genes: list[str],
    output: Path,
) -> None:
    subset = final.set_index("gene_symbol").loc[candidate_genes]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7), constrained_layout=True)
    y = np.arange(len(candidate_genes))
    axes[0].barh(y - 0.18, subset["lodo_top5_frequency"], 0.36, color="#2f6f8f", label="LODO top 5")
    axes[0].barh(y + 0.18, subset["bootstrap_top5_frequency"], 0.36, color="#d67b5d", label="Bootstrap top 5")
    axes[0].set_yticks(y, candidate_genes)
    axes[0].invert_yaxis()
    axes[0].set_xlim(0, 1.02)
    axes[0].set_xlabel("Selection frequency")
    axes[0].set_title("a  Identity stability", loc="left", fontweight="bold")
    axes[0].legend(frameon=False)
    axes[0].grid(axis="x", alpha=0.2)

    axes[1].errorbar(
        subset["bootstrap_beta_median"],
        y,
        xerr=[
            subset["bootstrap_beta_median"] - subset["bootstrap_beta_ci_low"],
            subset["bootstrap_beta_ci_high"] - subset["bootstrap_beta_median"],
        ],
        fmt="o",
        color="#2f6f8f",
        ecolor="#6f98ad",
        capsize=3,
    )
    axes[1].axvline(0, color="black", linewidth=0.8)
    axes[1].set_yticks(y, candidate_genes)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Bootstrap severity coefficient (median, 95% percentile interval)")
    axes[1].set_title("b  Resampling uncertainty", loc="left", fontweight="bold")
    axes[1].grid(axis="x", alpha=0.2)
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    project = args.project_root.resolve()
    package = args.package_root.resolve()
    bundled_project = package / "project_inputs"
    if bundled_project.is_dir() and not (project / "data_processed").is_dir():
        project = bundled_project
    results = package / "analysis_results"
    figures = package / "figures"
    results.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    expression, metadata = load_donor_data(project)
    diabetic = metadata["disease_group"].ne("healthy control")
    expression = expression.loc[diabetic]
    metadata = metadata.loc[diabetic]
    hallmark = read_hallmark(project / "data_raw" / "hallmark_inflammatory_response.gmt")
    genes = [gene for gene in hallmark if gene in expression.columns]

    total = fit_gene_models(
        expression, metadata, genes, model_name="total_association", include_dme=False
    )
    conditioned = fit_gene_models(
        expression, metadata, genes, model_name="DME_conditioned", include_dme=True
    )
    total_wild = wild_bootstrap_t(
        expression,
        metadata,
        genes,
        include_dme=False,
        iterations=args.wild_bootstrap,
        seed=args.seed + 100,
    )
    conditioned_wild = wild_bootstrap_t(
        expression,
        metadata,
        genes,
        include_dme=True,
        iterations=args.wild_bootstrap,
        seed=args.seed + 200,
    )
    total = total.merge(total_wild, on="gene_symbol", how="left")
    conditioned = conditioned.merge(conditioned_wild, on="gene_symbol", how="left")
    ranking = pd.concat([total, conditioned], ignore_index=True)

    stability = bootstrap_stability(
        expression, metadata, genes, args.bootstrap, args.seed + 300
    )
    lodo, lodo_detail = lodo_stability(expression, metadata, genes)
    primary = total.merge(stability, on="gene_symbol").merge(lodo, on="gene_symbol")
    primary = primary.sort_values("rank_descending_beta").reset_index(drop=True)

    total_top5 = total.head(5)["gene_symbol"].tolist()
    conditioned_top5 = conditioned.head(5)["gene_symbol"].tolist()
    candidate_genes = list(dict.fromkeys(["P2RX4", *total_top5, *conditioned_top5]))
    overlap = set(total_top5) & set(conditioned_top5)
    primary["candidate_tier"] = np.select(
        [
            primary["gene_symbol"].eq("P2RX4"),
            primary["gene_symbol"].isin(overlap - {"P2RX4"}),
            primary["gene_symbol"].isin((set(total_top5) | set(conditioned_top5)) - overlap),
        ],
        [
            "discovery-selected; externally target-locked",
            "cross-estimand secondary candidate",
            "estimand-dependent candidate",
        ],
        default="not prioritized",
    )

    influence_summary, influence_detail_all = influence_diagnostics(
        expression, metadata, candidate_genes
    )
    influence_detail = influence_detail_all.loc[
        influence_detail_all["gene_symbol"].isin(candidate_genes)
    ].copy()
    primary = primary.merge(influence_summary, on="gene_symbol", how="left")
    leverage_table = influence_detail[["donor", "leverage", "high_leverage_threshold_2p_over_n", "is_high_leverage", "pmi_minutes"]].drop_duplicates()
    high_leverage_donors = leverage_table.loc[leverage_table["is_high_leverage"], "donor"].tolist()

    sensitivities = model_sensitivities(
        expression, metadata, genes, high_leverage_donors
    )
    interaction = dme_interaction_p2rx4(expression, metadata)
    contrasts = clinical_contrasts(expression, metadata, genes)
    geo_metadata = parse_geo_metadata(project / "GSE102485_series_matrix.txt.gz")
    external = gse102485_membrane_sensitivity(
        project, candidate_genes, geo_metadata, args.seed + 400
    )
    predictions, predictive_summary = nested_predictive_audit(
        expression, metadata, genes, args.seed + 500
    )

    donor_table = metadata.reset_index().merge(expression[genes].reset_index(), on="donor")
    donor_table.to_csv(results / "final_donor_aggregated_data.csv", index=False)
    ranking.to_csv(results / "final_total_and_dme_conditioned_models.csv", index=False)
    primary.to_csv(results / "final_primary_ranking_and_stability.csv", index=False)
    stability.to_csv(results / "final_bootstrap_stability.csv", index=False)
    lodo.to_csv(results / "final_lodo_stability.csv", index=False)
    lodo_detail.loc[lodo_detail["gene_symbol"].isin(candidate_genes)].to_csv(
        results / "final_candidate_lodo_detail.csv", index=False
    )
    influence_detail.to_csv(results / "final_candidate_influence_diagnostics.csv", index=False)
    leverage_table.to_csv(results / "final_design_leverage.csv", index=False)
    sensitivities.to_csv(results / "final_severity_model_sensitivities.csv", index=False)
    interaction.to_csv(results / "final_p2rx4_dme_interaction.csv", index=False)
    contrasts.to_csv(results / "final_clinical_contrasts.csv", index=False)
    external.to_csv(results / "final_gse102485_membrane_only_context.csv", index=False)
    geo_metadata.to_csv(results / "final_gse102485_sample_metadata.csv", index=False)
    predictions.to_csv(results / "final_nested_predictive_predictions.csv", index=False)
    predictive_summary.to_csv(results / "final_nested_predictive_summary.csv", index=False)

    draw_discovery_figure(
        ranking,
        candidate_genes,
        figures / "Figure_2_total_and_dme_conditioned_associations",
    )
    draw_stability_figure(
        primary,
        candidate_genes,
        figures / "Figure_3_candidate_stability_and_uncertainty",
    )

    primary_p2rx4 = primary.set_index("gene_symbol").loc["P2RX4"]
    summary = {
        "analysis_version": "discover-research-article-final-v3",
        "seed": args.seed,
        "n_diabetic_donors": len(metadata),
        "group_counts": metadata["detailed_group"].value_counts().to_dict(),
        "n_inflammatory_genes": len(genes),
        "primary_estimand": "total association: standardized expression ~ donor-mean supplied DR severity + age + sex + PMI + RIN",
        "secondary_estimand": "DME-conditioned association: primary model + DME",
        "small_sample_inference": "HC3 covariance with residual-df t inference plus studentized Rademacher wild bootstrap-t",
        "total_top5": total_top5,
        "dme_conditioned_top5": conditioned_top5,
        "candidate_union": candidate_genes,
        "high_leverage_donors": high_leverage_donors,
        "p2rx4_total_beta": float(primary_p2rx4["severity_beta"]),
        "p2rx4_total_p_t": float(primary_p2rx4["severity_pvalue_t"]),
        "p2rx4_total_wild_bootstrap_p": float(primary_p2rx4["wild_bootstrap_t_two_sided_p"]),
        "p2rx4_bootstrap_top5_frequency": float(primary_p2rx4["bootstrap_top5_frequency"]),
        "p2rx4_lodo_top5_frequency": float(primary_p2rx4["lodo_top5_frequency"]),
        "fdr_significant_total_model": int((total["severity_padj_bh_158"] < 0.05).sum()),
        "wild_bootstrap_iterations": args.wild_bootstrap,
        "donor_bootstrap_iterations": args.bootstrap,
        "gse102485_role": "supplementary membrane-only post-hoc context; excluded from selection and validation",
    }
    (results / "final_discovery_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
