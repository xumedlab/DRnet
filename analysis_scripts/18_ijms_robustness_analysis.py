from __future__ import annotations

import runpy
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

from pipeline_utils import (
    bh_adjust,
    ensure_dirs,
    load_ensembl_symbol_mapping,
    log_message,
    matrix_ensembl_to_symbol,
    normalize_ensembl_id,
    read_gmt,
    spearman,
)


cfg = runpy.run_path("00_config.py")
RAW_DIR: Path = cfg["RAW_DIR"]
PROC_DIR: Path = cfg["PROC_DIR"]
RESULT_DIR: Path = cfg["RESULT_DIR"]
TABLE_DIR: Path = cfg["TABLE_DIR"]
FIG_DIR: Path = cfg["FIG_DIR"]

SELECTED_GENES = ["MSR1", "NMI", "FZD5", "TIMP1", "CMKLR1", "LYN", "TLR3"]
PRIMARY_CTRL = cfg["PRIMARY_CTRL"]
PRIMARY_CASE = cfg["PRIMARY_CASE"]
SEED = cfg.get("RANDOM_SEED", 202501)


def load_matrix_tsv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    first_col = df.columns[0]
    if first_col not in {"gene", "ensemblID"}:
        df = df.rename(columns={first_col: "gene"})
    else:
        df = df.rename(columns={first_col: "gene"})
    df["gene"] = df["gene"].map(normalize_ensembl_id)
    return df.set_index("gene").T


def build_symbol_expression_from_ensembl(expr: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    matrix_ensembl = {
        normalize_ensembl_id(gene): expr[gene].astype(float).tolist()
        for gene in expr.columns
    }
    matrix_symbol, _ = matrix_ensembl_to_symbol(matrix_ensembl, mapping)
    return pd.DataFrame(matrix_symbol, index=expr.index).astype(float)


def read_selected_coefficients() -> dict[str, float]:
    path = TABLE_DIR / "lasso_coefficients.csv"
    rows = pd.read_csv(path)
    return dict(zip(rows["gene_symbol"], rows["coefficient"]))


def load_inflammatory_symbols() -> set[str]:
    sets = read_gmt("hallmark_inflammatory_response.gmt")
    if "HALLMARK_INFLAMMATORY_RESPONSE" in sets:
        return {gene.upper() for gene in sets["HALLMARK_INFLAMMATORY_RESPONSE"]}
    first = next(iter(sets.values()))
    return {gene.upper() for gene in first}


def threshold_sensitivity() -> pd.DataFrame:
    deg = pd.read_csv(TABLE_DIR / "deg_primary_healthy_vs_npdr_pdr_dme.csv")
    trend = pd.read_csv(TABLE_DIR / "severity_trend_all_genes.csv")
    merged = deg.merge(
        trend[["gene", "spearman_rho", "padj"]].rename(
            columns={"padj": "severity_fdr", "spearman_rho": "severity_rho"}
        ),
        on="gene",
        how="left",
    )
    inflammatory = load_inflammatory_symbols()
    rows = []
    for padj_cutoff in [0.01, 0.05, 0.10]:
        for lfc_cutoff in [0.25, 0.50, 1.00]:
            sig = merged[
                (merged["padj"] < padj_cutoff)
                & (merged["log2FC"].abs() >= lfc_cutoff)
                & (merged["gene_symbol"].fillna("").str.upper().isin(inflammatory))
            ].copy()
            progressive = sig[(sig["severity_rho"] > 0) & (sig["severity_fdr"] < 0.1)]
            sig_symbols = set(sig["gene_symbol"].dropna().astype(str))
            progressive_symbols = set(progressive["gene_symbol"].dropna().astype(str))
            selected_in_sig = [gene for gene in SELECTED_GENES if gene in sig_symbols]
            selected_in_progressive = [gene for gene in SELECTED_GENES if gene in progressive_symbols]
            rows.append(
                {
                    "de_fdr_cutoff": padj_cutoff,
                    "abs_log2fc_cutoff": lfc_cutoff,
                    "inflammatory_candidate_count": len(sig),
                    "progression_related_candidate_count": len(progressive),
                    "selected_gene_recovery_count": len(selected_in_sig),
                    "selected_gene_recovery_fraction": len(selected_in_sig) / len(SELECTED_GENES),
                    "selected_progressive_recovery_count": len(selected_in_progressive),
                    "selected_progressive_recovery_fraction": len(selected_in_progressive) / len(SELECTED_GENES),
                    "selected_genes_recovered": ";".join(selected_in_sig),
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "ijms_robustness_threshold_sensitivity.csv", index=False)
    return out


def median_ratio_normalize(counts: pd.DataFrame) -> pd.DataFrame:
    positive = counts.replace(0, np.nan)
    geometric_means = np.exp(np.log(positive).mean(axis=0, skipna=True))
    valid = np.isfinite(geometric_means) & (geometric_means > 0)
    ratios = counts.loc[:, valid] / geometric_means[valid]
    size_factors = ratios.replace([np.inf, -np.inf], np.nan).median(axis=1, skipna=True)
    size_factors = size_factors.replace(0, np.nan).fillna(size_factors.median())
    return counts.div(size_factors, axis=0)


def upper_quartile_normalize(counts: pd.DataFrame) -> pd.DataFrame:
    uq = counts.quantile(0.75, axis=1)
    uq = uq.replace(0, np.nan).fillna(uq.median())
    scale = uq / uq.median()
    return counts.div(scale, axis=0)


def selected_gene_score_auc(expr: pd.DataFrame, pheno: pd.DataFrame, coefficients: dict[str, float]) -> float:
    primary_ids = pheno[pheno["disease_group"].isin([PRIMARY_CTRL, PRIMARY_CASE])].index
    genes = [gene for gene in SELECTED_GENES if gene in expr.columns and gene in coefficients]
    if not genes:
        return float("nan")
    z = pd.DataFrame(
        StandardScaler().fit_transform(expr.loc[primary_ids, genes]),
        index=primary_ids,
        columns=genes,
    )
    score = sum(z[gene] * coefficients[gene] for gene in genes)
    y = (pheno.loc[primary_ids, "disease_group"] == PRIMARY_CASE).astype(int)
    return float(roc_auc_score(y, score))


def normalization_sensitivity() -> pd.DataFrame:
    counts = load_matrix_tsv(PROC_DIR / "counts_macula_4groups.tsv").astype(float)
    current_log2cpm = load_matrix_tsv(PROC_DIR / "log2cpm_macula_4groups.tsv").astype(float)
    pheno = pd.read_csv(PROC_DIR / "pheno_macula_4groups.csv").set_index("sample_id")
    mapping = load_ensembl_symbol_mapping(RAW_DIR, PROC_DIR)
    coefficients = read_selected_coefficients()

    transforms = {
        "current_log2cpm": current_log2cpm,
        "median_ratio_log2norm": np.log2(median_ratio_normalize(counts) + 1),
        "upper_quartile_log2norm": np.log2(upper_quartile_normalize(counts) + 1),
    }

    rows = []
    severity = pheno.loc[counts.index, "severity_code"].astype(int).tolist()
    primary_ids = pheno[pheno["disease_group"].isin([PRIMARY_CTRL, PRIMARY_CASE])].index
    for method, ensembl_expr in transforms.items():
        symbol_expr = build_symbol_expression_from_ensembl(ensembl_expr, mapping)
        score_auc = selected_gene_score_auc(symbol_expr, pheno, coefficients)
        for gene in SELECTED_GENES:
            if gene not in symbol_expr.columns:
                rows.append(
                    {
                        "normalization_method": method,
                        "gene_symbol": gene,
                        "severity_rho": np.nan,
                        "severity_pvalue": np.nan,
                        "primary_delta_case_minus_control": np.nan,
                        "direction_positive": 0,
                        "apparent_weighted_score_auc": score_auc,
                    }
                )
                continue
            values = symbol_expr.loc[counts.index, gene].tolist()
            rho, pvalue = spearman(severity, values)
            case_ids = pheno.loc[primary_ids].index[
                pheno.loc[primary_ids, "disease_group"] == PRIMARY_CASE
            ]
            ctrl_ids = pheno.loc[primary_ids].index[
                pheno.loc[primary_ids, "disease_group"] == PRIMARY_CTRL
            ]
            case_mean = symbol_expr.loc[case_ids, gene].mean()
            ctrl_mean = symbol_expr.loc[ctrl_ids, gene].mean()
            rows.append(
                {
                    "normalization_method": method,
                    "gene_symbol": gene,
                    "severity_rho": rho,
                    "severity_pvalue": pvalue,
                    "primary_delta_case_minus_control": float(case_mean - ctrl_mean),
                    "direction_positive": int(rho > 0 and case_mean > ctrl_mean),
                    "apparent_weighted_score_auc": score_auc,
                }
            )

    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "ijms_robustness_normalization_sensitivity.csv", index=False)
    return out


def prepare_metadata_for_design(pheno: pd.DataFrame, covariates: list[str]) -> tuple[pd.DataFrame, str]:
    subset = pheno[pheno["disease_group"].isin([PRIMARY_CTRL, PRIMARY_CASE])].copy()
    subset["condition"] = pd.Categorical(
        subset["disease_group"],
        categories=[PRIMARY_CTRL, PRIMARY_CASE],
        ordered=True,
    )
    included = []
    for covariate in covariates:
        if covariate not in subset.columns:
            continue
        if covariate in {"age", "post_mortem_interval_min", "rin"}:
            values = pd.to_numeric(subset[covariate].replace("", np.nan), errors="coerce")
            if values.isna().any() or values.nunique(dropna=True) < 2 or values.std(ddof=0) == 0:
                continue
            subset[covariate] = (values - values.mean()) / values.std(ddof=0)
            included.append(covariate)
        elif covariate == "sex":
            values = subset[covariate].astype(str).str.strip()
            if values.isna().any() or values.nunique(dropna=True) < 2:
                continue
            subset[covariate] = pd.Categorical(values)
            included.append(covariate)
    design_terms = [*included, "condition"]
    design = "~" + " + ".join(design_terms)
    return subset[[*included, "condition"]], design


def run_primary_deseq_for_design(
    counts: pd.DataFrame,
    pheno: pd.DataFrame,
    mapping: dict[str, str],
    design_name: str,
    covariates: list[str],
) -> pd.DataFrame:
    meta, design = prepare_metadata_for_design(pheno, covariates)
    subset_counts = counts.loc[meta.index].astype(int)
    dds = DeseqDataSet(
        counts=subset_counts,
        metadata=meta,
        design=design,
        n_cpus=1,
        quiet=True,
    )
    dds.deseq2()
    stat = DeseqStats(dds, contrast=["condition", PRIMARY_CASE, PRIMARY_CTRL], quiet=True, n_cpus=1)
    stat.summary()
    result = stat.results_df.reset_index().rename(columns={"index": "gene", "log2FoldChange": "log2FC"})
    result["gene"] = result["gene"].map(normalize_ensembl_id)
    result["gene_symbol"] = result["gene"].map(mapping)
    for column in ["baseMean", "log2FC", "lfcSE", "stat", "pvalue", "padj"]:
        result[column] = pd.to_numeric(result[column], errors="coerce")
    result["pvalue"] = result["pvalue"].fillna(1.0)
    result["padj"] = result["padj"].fillna(1.0)
    result["design_name"] = design_name
    result["design_formula"] = design
    result.to_csv(TABLE_DIR / f"ijms_robustness_deseq_{design_name}.csv", index=False)
    return result


def covariate_design_sensitivity() -> pd.DataFrame:
    counts = load_matrix_tsv(PROC_DIR / "counts_macula_4groups.tsv")
    pheno = pd.read_csv(PROC_DIR / "pheno_macula_4groups.csv").set_index("sample_id")
    mapping = load_ensembl_symbol_mapping(RAW_DIR, PROC_DIR)
    inflammatory = load_inflammatory_symbols()
    designs = {
        "unadjusted": [],
        "age_sex": ["age", "sex"],
        "full_available": ["age", "sex", "post_mortem_interval_min", "rin"],
    }

    rows = []
    for design_name, covariates in designs.items():
        result = run_primary_deseq_for_design(counts, pheno, mapping, design_name, covariates)
        sig = result[
            (result["padj"] < 0.05)
            & (result["log2FC"].abs() >= 0.5)
            & (result["gene_symbol"].fillna("").str.upper().isin(inflammatory))
        ]
        sig_symbols = set(sig["gene_symbol"].dropna().astype(str))
        selected = result[result["gene_symbol"].isin(SELECTED_GENES)].copy()
        for gene in SELECTED_GENES:
            gene_row = selected[selected["gene_symbol"] == gene]
            if gene_row.empty:
                rows.append(
                    {
                        "design_name": design_name,
                        "design_formula": result["design_formula"].iloc[0],
                        "gene_symbol": gene,
                        "log2FC": np.nan,
                        "padj": np.nan,
                        "passes_primary_rule": 0,
                        "inflammatory_candidate_count": len(sig),
                    }
                )
            else:
                rec = gene_row.iloc[0]
                rows.append(
                    {
                        "design_name": design_name,
                        "design_formula": rec["design_formula"],
                        "gene_symbol": gene,
                        "log2FC": rec["log2FC"],
                        "padj": rec["padj"],
                        "passes_primary_rule": int(gene in sig_symbols),
                        "inflammatory_candidate_count": len(sig),
                    }
                )

    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "ijms_robustness_covariate_design_sensitivity.csv", index=False)
    return out


def leave_one_stage_stability() -> pd.DataFrame:
    log2cpm = load_matrix_tsv(PROC_DIR / "log2cpm_macula_4groups.tsv").astype(float)
    pheno = pd.read_csv(PROC_DIR / "pheno_macula_4groups.csv").set_index("sample_id")
    mapping = load_ensembl_symbol_mapping(RAW_DIR, PROC_DIR)
    symbol_expr = build_symbol_expression_from_ensembl(log2cpm, mapping)

    rows = []
    groups = list(cfg["GROUPS"])
    for omitted_group in groups:
        keep_ids = pheno[pheno["disease_group"] != omitted_group].index.tolist()
        severity = pheno.loc[keep_ids, "severity_code"].astype(int).tolist()
        pvalues = []
        gene_records = []
        for gene in SELECTED_GENES:
            values = symbol_expr.loc[keep_ids, gene].tolist()
            rho, pvalue = spearman(severity, values)
            pvalues.append(pvalue)
            gene_records.append((gene, rho, pvalue))
        adjusted = bh_adjust(pvalues)
        for (gene, rho, pvalue), padj in zip(gene_records, adjusted):
            rows.append(
                {
                    "omitted_stage": omitted_group,
                    "n_samples_remaining": len(keep_ids),
                    "gene_symbol": gene,
                    "severity_rho": rho,
                    "pvalue": pvalue,
                    "padj_within_selected_genes": padj,
                    "direction_positive": int(rho > 0),
                    "positive_and_fdr_lt_0_1": int(rho > 0 and padj < 0.1),
                }
            )

    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "ijms_robustness_leave_one_stage_stability.csv", index=False)
    return out


def make_summary(
    threshold_df: pd.DataFrame,
    norm_df: pd.DataFrame,
    covar_df: pd.DataFrame,
    leave_df: pd.DataFrame,
) -> pd.DataFrame:
    canonical_threshold = threshold_df[
        (threshold_df["de_fdr_cutoff"] == 0.05)
        & (threshold_df["abs_log2fc_cutoff"] == 0.50)
    ].iloc[0]
    threshold_min = threshold_df["selected_gene_recovery_count"].min()
    norm_positive = (
        norm_df.groupby("normalization_method")["direction_positive"].sum().min()
    )
    covar_recovered = covar_df.groupby("design_name")["passes_primary_rule"].sum()
    leave_positive = leave_df.groupby("omitted_stage")["direction_positive"].sum()
    leave_fdr = leave_df.groupby("omitted_stage")["positive_and_fdr_lt_0_1"].sum()

    rows = [
        {
            "analysis": "Threshold sensitivity",
            "summary": (
                f"Canonical thresholds recovered {int(canonical_threshold['selected_gene_recovery_count'])}/"
                f"{len(SELECTED_GENES)} selected genes; across all tested thresholds the minimum recovery was "
                f"{int(threshold_min)}/{len(SELECTED_GENES)}."
            ),
            "support_level": "Stable at the canonical and less stringent thresholds; stricter log2FC cutoffs reduce recovery.",
        },
        {
            "analysis": "Normalization sensitivity",
            "summary": (
                f"All normalization methods retained positive severity/primary directions for at least "
                f"{int(norm_positive)}/{len(SELECTED_GENES)} selected genes."
            ),
            "support_level": "Directionally stable across log2CPM, median-ratio, and upper-quartile transforms.",
        },
        {
            "analysis": "Covariate-design sensitivity",
            "summary": (
                "Selected genes passing the primary DE rule by design: "
                + "; ".join(f"{name}={int(value)}/{len(SELECTED_GENES)}" for name, value in covar_recovered.items())
                + "."
            ),
            "support_level": "Most candidates remain directionally positive, but adjusted significance is design-sensitive.",
        },
        {
            "analysis": "Leave-one-stage stability",
            "summary": (
                "Positive severity direction by omitted stage: "
                + "; ".join(f"{name}={int(value)}/{len(SELECTED_GENES)}" for name, value in leave_positive.items())
                + ". Positive and FDR<0.1: "
                + "; ".join(f"{name}={int(value)}/{len(SELECTED_GENES)}" for name, value in leave_fdr.items())
                + "."
            ),
            "support_level": "Direction is stable; significance weakens when endpoint stages are omitted.",
        },
    ]
    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "ijms_robustness_summary.csv", index=False)
    write_summary_md(out)
    return out


def write_summary_md(summary: pd.DataFrame) -> None:
    lines = ["# IJMS Computational Robustness Summary", ""]
    for _, row in summary.iterrows():
        lines.extend(
            [
                f"## {row['analysis']}",
                "",
                str(row["summary"]),
                "",
                f"Interpretation: {row['support_level']}",
                "",
            ]
        )
    (TABLE_DIR / "ijms_robustness_summary.md").write_text("\n".join(lines), encoding="utf-8")


def make_robustness_figure(summary: pd.DataFrame) -> None:
    threshold = pd.read_csv(TABLE_DIR / "ijms_robustness_threshold_sensitivity.csv")
    norm = pd.read_csv(TABLE_DIR / "ijms_robustness_normalization_sensitivity.csv")
    covar = pd.read_csv(TABLE_DIR / "ijms_robustness_covariate_design_sensitivity.csv")
    leave = pd.read_csv(TABLE_DIR / "ijms_robustness_leave_one_stage_stability.csv")

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), dpi=200)
    axes = axes.ravel()

    pivot = threshold.pivot(
        index="abs_log2fc_cutoff",
        columns="de_fdr_cutoff",
        values="selected_gene_recovery_count",
    ).sort_index(ascending=False)
    im = axes[0].imshow(pivot.values, vmin=0, vmax=len(SELECTED_GENES), cmap="YlGnBu")
    axes[0].set_xticks(range(len(pivot.columns)), [str(x) for x in pivot.columns])
    axes[0].set_yticks(range(len(pivot.index)), [str(x) for x in pivot.index])
    axes[0].set_xlabel("DE FDR cutoff")
    axes[0].set_ylabel("|log2FC| cutoff")
    axes[0].set_title("A. Threshold recovery")
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            axes[0].text(j, i, int(pivot.values[i, j]), ha="center", va="center")
    fig.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04)

    norm_counts = norm.groupby("normalization_method")["direction_positive"].sum()
    norm_labels = {
        "current_log2cpm": "log2CPM",
        "median_ratio_log2norm": "median ratio",
        "upper_quartile_log2norm": "upper quartile",
    }
    axes[1].bar([norm_labels.get(x, x) for x in norm_counts.index], norm_counts.values, color="#2a9d8f")
    axes[1].set_ylim(0, len(SELECTED_GENES))
    axes[1].set_ylabel("Positive genes")
    axes[1].set_title("B. Normalization direction")
    axes[1].tick_params(axis="x", rotation=25)

    covar_counts = covar.groupby("design_name")["passes_primary_rule"].sum()
    covar_labels = {
        "age_sex": "age + sex",
        "full_available": "full covariates",
        "unadjusted": "unadjusted",
    }
    axes[2].bar([covar_labels.get(x, x) for x in covar_counts.index], covar_counts.values, color="#5271a3")
    axes[2].set_ylim(0, len(SELECTED_GENES))
    axes[2].set_ylabel("Genes passing primary rule")
    axes[2].set_title("C. Covariate design")
    axes[2].tick_params(axis="x", rotation=15)

    leave_counts = leave.groupby("omitted_stage")["direction_positive"].sum().reindex(cfg["GROUPS"])
    axes[3].bar(leave_counts.index, leave_counts.values, color="#c77d1a")
    axes[3].set_ylim(0, len(SELECTED_GENES))
    axes[3].set_ylabel("Positive genes")
    axes[3].set_title("D. Leave-one-stage direction")
    axes[3].tick_params(axis="x", rotation=20)

    fig.tight_layout()
    for ext in cfg.get("FIG_FORMATS", ("png", "pdf")):
        fig.savefig(FIG_DIR / f"Supplementary_Figure_S6_robustness.{ext}", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ensure_dirs(TABLE_DIR, FIG_DIR, RESULT_DIR / "logs")
    threshold_df = threshold_sensitivity()
    norm_df = normalization_sensitivity()
    covar_df = covariate_design_sensitivity()
    leave_df = leave_one_stage_stability()
    summary = make_summary(threshold_df, norm_df, covar_df, leave_df)
    make_robustness_figure(summary)
    log_message("18_ijms_robustness_analysis", "IJMS robustness analyses completed.")


if __name__ == "__main__":
    main()
