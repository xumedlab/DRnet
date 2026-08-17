#!/usr/bin/env python3
"""Validate remote RNA-seq reconstruction and fit raw-count P2RX4 models.

The analysis is deliberately target focused.  It verifies the transferred
archive and every remotely checksummed output, aggregates at biological-sample
grain, fits PyDESeq2 negative-binomial models, and quantifies source, QC,
outlier, and reference-assembly sensitivity.  The raw-read reconstruction is a
post-protocol deviation and is not presented as a preregistered confirmation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import tarfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.stats import mannwhitneyu, norm, spearmanr
from statsmodels.stats.multitest import multipletests


TARGET_GENE = "P2RX4"
TARGET_GENE_ID = "ENSG00000135124.16"
EXPECTED_ARCHIVE_SHA256 = (
    "CF6784D6852D30A7A2A67FDFC0C7FB93E0A3ABBF3207CE775F860709795D65A6"
)
WORKFLOWS = {
    "all_regions_source_reconstruction": "all-regions source reconstruction",
    "primary_assembly_sensitivity": "primary-assembly sensitivity",
}


def parse_args() -> argparse.Namespace:
    package = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path, default=package)
    parser.add_argument("--n-cpus", type=int, default=4)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def json_safe(value: object) -> object:
    """Convert pandas/NumPy scalars and non-finite values to strict JSON."""
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def verify_output_manifest(workflow_root: Path) -> dict[str, object]:
    manifest = workflow_root / "output_sha256.txt"
    lines = manifest.read_text(encoding="utf-8").splitlines()
    pattern = re.compile(r"^([0-9a-fA-F]{64})\s+(.+)$")
    failures: list[dict[str, str]] = []
    checked = 0
    for line in lines:
        match = pattern.match(line)
        if match is None:
            failures.append({"path": "", "reason": f"malformed line: {line}"})
            continue
        expected, relative = match.groups()
        path = workflow_root / relative
        if not path.is_file():
            failures.append({"path": relative, "reason": "missing"})
            continue
        actual = sha256(path)
        if actual != expected.upper():
            failures.append(
                {"path": relative, "reason": f"SHA-256 {actual} != {expected}"}
            )
        checked += 1
    return {
        "manifest_entries": len(lines),
        "checked_files": checked,
        "failure_count": len(failures),
        "failures": failures,
    }


def validate_remote_bundle(remote_dir: Path) -> dict[str, object]:
    archive = remote_dir / "DRnet_GSE276892_remote_results.tar.gz"
    checksum_file = remote_dir / "DRnet_GSE276892_remote_results.tar.gz.sha256"
    extracted = remote_dir / "extracted"
    checksum_text = checksum_file.read_text(encoding="utf-8").strip()
    checksum_match = re.match(r"^([0-9a-fA-F]{64})\s+", checksum_text)
    if checksum_match is None:
        raise AssertionError("Malformed transferred archive SHA-256 file")
    recorded = checksum_match.group(1).upper()
    actual = sha256(archive)
    with archive.open("rb") as stream:
        gzip_header = stream.read(3).hex().upper()
    with tarfile.open(archive, "r:gz") as tar:
        members = tar.getmembers()
        unsafe = [
            member.name
            for member in members
            if Path(member.name).is_absolute() or ".." in Path(member.name).parts
        ]
    if actual != EXPECTED_ARCHIVE_SHA256 or recorded != actual:
        raise AssertionError(
            f"Archive SHA-256 mismatch: expected={EXPECTED_ARCHIVE_SHA256}, "
            f"recorded={recorded}, actual={actual}"
        )
    if gzip_header != "1F8B08":
        raise AssertionError(f"Unexpected gzip header: {gzip_header}")
    if unsafe:
        raise AssertionError(f"Unsafe archive members: {unsafe[:5]}")
    if not extracted.is_dir():
        raise AssertionError(f"Extracted remote result directory missing: {extracted}")

    workflow_checks = {
        workflow: verify_output_manifest(extracted / workflow)
        for workflow in WORKFLOWS
    }
    if any(check["failure_count"] for check in workflow_checks.values()):
        raise AssertionError(f"Remote output checksum failure: {workflow_checks}")
    return {
        "archive_path": str(archive),
        "archive_bytes": archive.stat().st_size,
        "archive_sha256": actual,
        "gzip_header": gzip_header,
        "tar_member_count": len(members),
        "unsafe_member_count": len(unsafe),
        "extracted_file_count": sum(1 for path in extracted.rglob("*") if path.is_file()),
        "workflow_checks": workflow_checks,
    }


def load_clinical_metadata(package: Path) -> pd.DataFrame:
    path = package / "analysis_results/Independent_validation_P2RX4_sample_level.csv"
    metadata = pd.read_csv(path)
    metadata = metadata.loc[metadata["dataset"].eq("GSE276892")].copy()
    if len(metadata) != 17 or metadata["sample_id"].duplicated().any():
        raise AssertionError("Expected 17 unique GSE276892 sample metadata rows")
    metadata["diagnosis"] = metadata["disease_group"].replace(
        {"PDR": "PDR", "control": "control"}
    )
    metadata["source_reused"] = metadata["reused_control"].astype(int)
    metadata["age_per_10y_centered"] = (
        metadata["age"] - metadata["age"].mean()
    ) / 10.0
    metadata["sex"] = metadata["sex"].astype(str)
    return metadata.set_index("sample_id", drop=False)


def load_workflow(
    extracted: Path,
    workflow: str,
    clinical: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    root = extracted / workflow / "results"
    matrix = pd.read_csv(root / "sample_gene_counts.tsv", sep="\t")
    sample_qc = pd.read_csv(root / "sample_qc.tsv", sep="\t")
    lane_qc = pd.read_csv(root / "lane_qc.tsv", sep="\t")
    p2rx4_file = pd.read_csv(root / "P2RX4_counts.tsv", sep="\t")
    if matrix["Geneid"].duplicated().any():
        raise AssertionError(f"Duplicate gene IDs in {workflow}")
    if len(sample_qc) != 17 or sample_qc["sample_id"].duplicated().any():
        raise AssertionError(f"Expected 17 unique sample-QC rows in {workflow}")
    if len(lane_qc) != 31 or lane_qc["run_accession"].nunique() != 31:
        raise AssertionError(f"Expected 31 unique lane-QC rows in {workflow}")
    if lane_qc["sample_id"].nunique() != 17:
        raise AssertionError(f"Expected 17 lane-QC sample IDs in {workflow}")
    if len(p2rx4_file) != 1 or p2rx4_file.iloc[0]["Geneid"] != TARGET_GENE_ID:
        raise AssertionError(f"Expected one release-42 P2RX4 row in {workflow}")
    sample_ids = sample_qc["sample_id"].tolist()
    if set(sample_ids) != set(clinical.index):
        raise AssertionError(f"Clinical/count sample mismatch in {workflow}")
    if any(sample not in matrix.columns for sample in sample_ids):
        raise AssertionError(f"Missing sample count column in {workflow}")
    counts = matrix.set_index("Geneid")[sample_ids].T
    if (counts < 0).any().any() or not np.equal(counts, np.floor(counts)).all().all():
        raise AssertionError(f"Non-integer or negative counts in {workflow}")
    counts = counts.loc[:, counts.sum(axis=0) > 0].astype(int)
    annotations = matrix.set_index("Geneid")[["gene_name", "gene_type"]]
    sample_qc = sample_qc.set_index("sample_id").loc[sample_ids]
    return counts, annotations, sample_qc, lane_qc


def _gene_value(dds: DeseqDataSet, column: str, gene_id: str) -> object:
    if column not in dds.var.columns:
        return np.nan
    value = dds.var.loc[gene_id, column]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return float(value) if pd.notna(value) else np.nan


def one_sided_greater_p(wald_statistic: float) -> float:
    """Return the standard normal upper-tail P for a directional Wald test."""
    return float(norm.sf(wald_statistic))


def fit_deseq2_model(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str,
    model_name: str,
    workflow: str,
    n_cpus: int,
    refit_cooks: bool = True,
    cooks_filter: bool = True,
) -> tuple[dict[str, object], pd.Series, pd.Series, pd.DataFrame]:
    metadata = metadata.loc[counts.index].copy()
    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design=design,
        fit_type="mean",
        refit_cooks=refit_cooks,
        n_cpus=n_cpus,
        quiet=True,
    )
    dds.deseq2()
    two_sided = DeseqStats(
        dds,
        contrast=["diagnosis", "PDR", "control"],
        cooks_filter=cooks_filter,
        independent_filter=False,
        n_cpus=n_cpus,
        quiet=True,
    )
    two_sided.summary()
    target = two_sided.results_df.loc[TARGET_GENE_ID]
    gene_index = counts.columns.get_loc(TARGET_GENE_ID)
    normalized = pd.Series(
        dds.layers["normed_counts"][:, gene_index],
        index=counts.index,
        name="p2rx4_deseq2_normalized_count",
    )
    size_factors = dds.obs["size_factors"].astype(float).rename("size_factor")
    design_matrix = np.asarray(dds.obsm["design_matrix"], dtype=float)
    cooks = np.asarray(dds.layers["cooks"][:, gene_index], dtype=float)
    result: dict[str, object] = {
        "workflow": workflow,
        "model": model_name,
        "design_formula": design,
        "n_samples": len(counts),
        "n_pdr": int(metadata["diagnosis"].eq("PDR").sum()),
        "n_control": int(metadata["diagnosis"].eq("control").sum()),
        "n_expressed_genes": counts.shape[1],
        "dispersion_fit_type": "mean",
        "base_mean": float(target["baseMean"]),
        "log2_fold_change_pdr_vs_control": float(target["log2FoldChange"]),
        "lfc_se": float(target["lfcSE"]),
        "wald_statistic": float(target["stat"]),
        "wald_p_two_sided": float(target["pvalue"]),
        "wald_p_one_sided_greater": one_sided_greater_p(float(target["stat"])),
        "bh_adjusted_p_all_expressed_genes": float(target["padj"]),
        "lfc_ci_low_95": float(target["log2FoldChange"] - 1.96 * target["lfcSE"]),
        "lfc_ci_high_95": float(target["log2FoldChange"] + 1.96 * target["lfcSE"]),
        "design_rank": int(np.linalg.matrix_rank(design_matrix)),
        "design_columns": int(design_matrix.shape[1]),
        "design_condition_number": float(np.linalg.cond(design_matrix)),
        "refit_cooks": refit_cooks,
        "cooks_filter": cooks_filter,
        "p2rx4_max_cooks_distance": float(np.nanmax(cooks)),
        "p2rx4_replaced": _gene_value(dds, "replaced", TARGET_GENE_ID),
        "p2rx4_refitted": _gene_value(dds, "refitted", TARGET_GENE_ID),
        "p2rx4_dispersion": _gene_value(dds, "dispersions", TARGET_GENE_ID),
        "analysis_role": "post-protocol raw-read reconstruction",
    }
    full_results = two_sided.results_df.reset_index(names="Geneid")
    return result, normalized, size_factors, full_results


def qc_correlations(sample_table: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    metrics = [
        "input_reads",
        "uniquely_mapped_percent",
        "featurecounts_assigned_percent",
        "total_gene_counts",
        "detected_genes",
        "size_factor",
        "cell_count",
        "rna_concentration_pg_per_ul",
        "age",
    ]
    for expression in ["p2rx4_raw_count", "p2rx4_deseq2_normalized_count"]:
        for metric in metrics:
            frame = sample_table[[expression, metric]].dropna()
            test = spearmanr(frame[expression], frame[metric])
            rows.append(
                {
                    "expression_scale": expression,
                    "qc_or_clinical_metric": metric,
                    "n_samples": len(frame),
                    "spearman_rho": float(test.statistic),
                    "spearman_p_two_sided": float(test.pvalue),
                }
            )
    result = pd.DataFrame(rows)
    result["spearman_p_bh_18_tests"] = multipletests(
        result["spearman_p_two_sided"], method="fdr_bh"
    )[1]
    result["scope"] = "descriptive QC/clinical association; not causal adjustment"
    return result


def leave_one_out(sample_table: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for omitted in [None, *sample_table.index.tolist()]:
        frame = sample_table if omitted is None else sample_table.drop(index=omitted)
        pdr = frame.loc[
            frame["diagnosis"].eq("PDR"), "p2rx4_deseq2_normalized_count"
        ].to_numpy(dtype=float)
        control = frame.loc[
            frame["diagnosis"].eq("control"),
            "p2rx4_deseq2_normalized_count",
        ].to_numpy(dtype=float)
        test = mannwhitneyu(pdr, control, alternative="two-sided", method="auto")
        log_difference = float(np.log2(pdr + 1).mean() - np.log2(control + 1).mean())
        rows.append(
            {
                "omitted_sample": "none" if omitted is None else omitted,
                "omitted_group": (
                    "none" if omitted is None else sample_table.loc[omitted, "diagnosis"]
                ),
                "n_pdr": len(pdr),
                "n_control": len(control),
                "log2_normalized_mean_difference": log_difference,
                "direction_positive": log_difference > 0,
                "mann_whitney_p_two_sided": float(test.pvalue),
            }
        )
    return pd.DataFrame(rows)


def reference_sensitivity(
    all_counts: pd.DataFrame,
    primary_counts: pd.DataFrame,
    all_sample_table: pd.DataFrame,
    primary_sample_table: pd.DataFrame,
) -> pd.DataFrame:
    common = all_counts.columns.intersection(primary_counts.columns)
    rows: list[dict[str, object]] = []
    for sample in all_counts.index:
        test = spearmanr(
            np.log1p(all_counts.loc[sample, common].to_numpy(dtype=float)),
            np.log1p(primary_counts.loc[sample, common].to_numpy(dtype=float)),
        )
        all_p2rx4 = int(all_counts.loc[sample, TARGET_GENE_ID])
        primary_p2rx4 = int(primary_counts.loc[sample, TARGET_GENE_ID])
        rows.append(
            {
                "sample_id": sample,
                "diagnosis": all_sample_table.loc[sample, "diagnosis"],
                "common_expressed_gene_ids": len(common),
                "spearman_log1p_gene_counts": float(test.statistic),
                "all_regions_total_gene_counts": int(all_counts.loc[sample].sum()),
                "primary_assembly_total_gene_counts": int(
                    primary_counts.loc[sample].sum()
                ),
                "total_count_ratio_primary_over_all": float(
                    primary_counts.loc[sample].sum() / all_counts.loc[sample].sum()
                ),
                "all_regions_p2rx4_count": all_p2rx4,
                "primary_assembly_p2rx4_count": primary_p2rx4,
                "p2rx4_count_difference_primary_minus_all": (
                    primary_p2rx4 - all_p2rx4
                ),
                "all_regions_normalized_p2rx4": float(
                    all_sample_table.loc[sample, "p2rx4_deseq2_normalized_count"]
                ),
                "primary_assembly_normalized_p2rx4": float(
                    primary_sample_table.loc[
                        sample, "p2rx4_deseq2_normalized_count"
                    ]
                ),
            }
        )
    return pd.DataFrame(rows)


def draw_figure(
    sample_table: pd.DataFrame,
    reference_table: pd.DataFrame,
    model_table: pd.DataFrame,
    output: Path,
) -> None:
    colors = {"PDR": "#b33c3c", "control": "#2f6f8f"}
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.3))

    for x, group in enumerate(["control", "PDR"]):
        subset = sample_table.loc[sample_table["diagnosis"].eq(group)]
        jitter = np.linspace(-0.12, 0.12, len(subset))
        axes[0, 0].scatter(
            np.full(len(subset), x) + jitter,
            subset["p2rx4_deseq2_normalized_count"],
            c=[colors[group]] * len(subset),
            edgecolor="white",
            linewidth=0.6,
            s=55,
        )
        for sample, row in subset.iterrows():
            if sample == "PDR_S10":
                axes[0, 0].annotate(
                    sample,
                    (x, row["p2rx4_deseq2_normalized_count"]),
                    xytext=(6, 4),
                    textcoords="offset points",
                    fontsize=8,
                )
    axes[0, 0].set_xticks([0, 1], ["Control", "PDR"])
    axes[0, 0].set_ylabel("DESeq2-normalized P2RX4 count")
    axes[0, 0].set_title("a  Raw-read reconstruction", loc="left", fontweight="bold")

    axes[0, 1].scatter(
        sample_table["detected_genes"],
        sample_table["p2rx4_deseq2_normalized_count"],
        c=sample_table["diagnosis"].map(colors),
        s=50,
        edgecolor="white",
        linewidth=0.6,
    )
    axes[0, 1].set_xlabel("Detected genes")
    axes[0, 1].set_ylabel("DESeq2-normalized P2RX4 count")
    axes[0, 1].set_title("b  QC relationship", loc="left", fontweight="bold")

    axes[1, 0].scatter(
        reference_table["all_regions_p2rx4_count"],
        reference_table["primary_assembly_p2rx4_count"],
        c=reference_table["diagnosis"].map(colors),
        s=50,
        edgecolor="white",
        linewidth=0.6,
    )
    lower = min(
        reference_table["all_regions_p2rx4_count"].min(),
        reference_table["primary_assembly_p2rx4_count"].min(),
    )
    upper = max(
        reference_table["all_regions_p2rx4_count"].max(),
        reference_table["primary_assembly_p2rx4_count"].max(),
    )
    axes[1, 0].plot([lower, upper], [lower, upper], color="black", linewidth=0.8)
    axes[1, 0].set_xlabel("All-regions P2RX4 count")
    axes[1, 0].set_ylabel("Primary-assembly P2RX4 count")
    axes[1, 0].set_title("c  Reference sensitivity", loc="left", fontweight="bold")

    display = model_table.loc[
        model_table["model"].isin(
            [
                "disease_only",
                "age_sex_adjusted",
                "source_adjusted",
                "source_age_sex_adjusted",
                "disease_only_without_PDR_S10",
                "new_profiles_only_8_vs_2",
            ]
        )
        & model_table["workflow"].eq("all_regions_source_reconstruction")
    ].copy()
    y = np.arange(len(display))
    axes[1, 1].errorbar(
        display["log2_fold_change_pdr_vs_control"],
        y,
        xerr=[
            display["log2_fold_change_pdr_vs_control"] - display["lfc_ci_low_95"],
            display["lfc_ci_high_95"] - display["log2_fold_change_pdr_vs_control"],
        ],
        fmt="o",
        color="#6b4c9a",
        capsize=3,
    )
    axes[1, 1].axvline(0, color="black", linewidth=0.8)
    axes[1, 1].set_yticks(y, display["model"].str.replace("_", " "))
    axes[1, 1].set_xlabel("PDR vs control log2 fold change (95% Wald CI)")
    axes[1, 1].set_title("d  Model sensitivity", loc="left", fontweight="bold")
    axes[1, 1].invert_yaxis()

    for axis in axes.flat:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(alpha=0.15)
    fig.tight_layout()
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    package = args.package_root.resolve()
    remote_dir = package / "analysis_data/independent_validation/remote_results"
    extracted = remote_dir / "extracted"
    results_dir = package / "analysis_results"
    figures_dir = package / "figures"
    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    integrity = validate_remote_bundle(remote_dir)
    clinical = load_clinical_metadata(package)
    workflow_data: dict[str, dict[str, object]] = {}
    for workflow in WORKFLOWS:
        counts, annotations, sample_qc, lane_qc = load_workflow(
            extracted, workflow, clinical
        )
        workflow_data[workflow] = {
            "counts": counts,
            "annotations": annotations,
            "sample_qc": sample_qc,
            "lane_qc": lane_qc,
        }

    model_rows: list[dict[str, object]] = []
    full_results_to_save: dict[str, pd.DataFrame] = {}
    primary_sample_tables: dict[str, pd.DataFrame] = {}
    for workflow in WORKFLOWS:
        data = workflow_data[workflow]
        counts = data["counts"]
        sample_qc = data["sample_qc"]
        metadata = clinical.loc[counts.index].copy()
        model, normalized, size_factors, full_results = fit_deseq2_model(
            counts,
            metadata,
            "~diagnosis",
            "disease_only",
            workflow,
            args.n_cpus,
        )
        model_rows.append(model)
        full_results_to_save[workflow] = full_results
        target_raw = counts[TARGET_GENE_ID].rename("p2rx4_raw_count")
        primary_sample_tables[workflow] = clinical.loc[counts.index].join(
            sample_qc[
                [
                    "input_reads",
                    "uniquely_mapped_reads",
                    "uniquely_mapped_percent",
                    "featurecounts_assigned_reads",
                    "featurecounts_assigned_percent",
                    "total_gene_counts",
                    "detected_genes",
                ]
            ]
        ).join(target_raw).join(normalized).join(size_factors)

    all_workflow = "all_regions_source_reconstruction"
    all_counts = workflow_data[all_workflow]["counts"]
    all_metadata = clinical.loc[all_counts.index].copy()
    additional_models = [
        ("~age_per_10y_centered + sex + diagnosis", "age_sex_adjusted", all_counts, all_metadata, True, True),
        ("~source_reused + diagnosis", "source_adjusted", all_counts, all_metadata, True, True),
        ("~source_reused + age_per_10y_centered + sex + diagnosis", "source_age_sex_adjusted", all_counts, all_metadata, True, True),
        ("~diagnosis", "disease_only_without_PDR_S10", all_counts.drop(index="PDR_S10"), all_metadata.drop(index="PDR_S10"), True, True),
        (
            "~diagnosis",
            "new_profiles_only_8_vs_2",
            all_counts.loc[all_metadata["source_reused"].eq(0)],
            all_metadata.loc[all_metadata["source_reused"].eq(0)],
            False,
            False,
        ),
    ]
    for design, name, counts, metadata, refit_cooks, cooks_filter in additional_models:
        model, _, _, _ = fit_deseq2_model(
            counts,
            metadata,
            design,
            name,
            all_workflow,
            args.n_cpus,
            refit_cooks=refit_cooks,
            cooks_filter=cooks_filter,
        )
        model_rows.append(model)

    model_table = pd.DataFrame(model_rows)
    all_sample_table = primary_sample_tables[all_workflow]
    primary_workflow = "primary_assembly_sensitivity"
    primary_sample_table = primary_sample_tables[primary_workflow]
    qc_table = qc_correlations(all_sample_table)
    loo_table = leave_one_out(all_sample_table)
    reference_table = reference_sensitivity(
        workflow_data[all_workflow]["counts"],
        workflow_data[primary_workflow]["counts"],
        all_sample_table,
        primary_sample_table,
    )

    model_table.to_csv(results_dir / "raw_count_p2rx4_deseq2_results.csv", index=False)
    all_sample_table.reset_index(drop=True).to_csv(
        results_dir / "raw_count_p2rx4_sample_level.csv", index=False
    )
    qc_table.to_csv(results_dir / "raw_count_p2rx4_qc_correlations.csv", index=False)
    loo_table.to_csv(results_dir / "raw_count_p2rx4_leave_one_out.csv", index=False)
    reference_table.to_csv(
        results_dir / "raw_count_reference_sensitivity_by_sample.csv", index=False
    )
    for workflow, full_results in full_results_to_save.items():
        annotations = workflow_data[workflow]["annotations"].reset_index()
        full_results.merge(annotations, on="Geneid", how="left").to_csv(
            results_dir / f"raw_count_deseq2_all_genes_{workflow}.csv", index=False
        )

    primary_model = model_table.loc[
        model_table["workflow"].eq(all_workflow)
        & model_table["model"].eq("disease_only")
    ].iloc[0]
    without_outlier = model_table.loc[
        model_table["model"].eq("disease_only_without_PDR_S10")
    ].iloc[0]
    source_adjusted = model_table.loc[
        model_table["model"].eq("source_adjusted")
    ].iloc[0]
    new_only = model_table.loc[
        model_table["model"].eq("new_profiles_only_8_vs_2")
    ].iloc[0]
    summary = {
        "analysis_version": "raw-count-reconstruction-v1",
        "analysis_role": "post-protocol raw-read reconstruction sensitivity",
        "integrity": integrity,
        "biological_sample_grain": {
            "n_samples": 17,
            "n_pdr": 8,
            "n_control": 9,
            "n_new_controls": 2,
            "n_reused_controls": 7,
            "technical_lanes": 31,
        },
        "primary_all_regions_disease_only": primary_model.to_dict(),
        "without_PDR_S10": without_outlier.to_dict(),
        "source_adjusted": source_adjusted.to_dict(),
        "new_profiles_only_8_vs_2": new_only.to_dict(),
        "reference_sensitivity": {
            "max_absolute_p2rx4_raw_count_difference": int(
                reference_table["p2rx4_count_difference_primary_minus_all"]
                .abs()
                .max()
            ),
            "minimum_sample_gene_count_spearman": float(
                reference_table["spearman_log1p_gene_counts"].min()
            ),
        },
        "interpretation_guardrail": (
            "Direction and inferential strength are judged across disease-only, "
            "source-adjusted, new-control-only, outlier, QC, and reference "
            "sensitivities; a single P value is not treated as validation."
        ),
    }
    summary = json_safe(summary)
    (results_dir / "raw_count_p2rx4_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    (results_dir / "raw_count_reconstruction_integrity.json").write_text(
        json.dumps(integrity, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    draw_figure(
        all_sample_table,
        reference_table,
        model_table,
        figures_dir / "Supplementary_Figure_raw_count_reconstruction",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
