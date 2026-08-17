#!/usr/bin/env python3
"""Describe candidate localization context in the normal GSE130636 retina atlas.

The biological unit is the donor-region library, not the cell. Expression and
detection are first summarised within donor x region x author-defined cell type,
then averaged across the six libraries. Library-, donor-, and leave-one-donor
summaries expose stability across the three normal-retina donors. This atlas is
not a diabetic-retinopathy replication cohort and cannot establish disease-state
cellular origin.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


CELL_TYPE_ORDER = [
    "Rod photoreceptors",
    "Cone photoreceptors",
    "Bipolar cells",
    "Retinal ganglion cells",
    "Horizontal cells",
    "Amacrine cells",
    "Pericytes",
    "Endothelial cells",
    "Microglia",
    "Müller-enriched glia",
]
GLIAL_MARKERS = ["RLBP1", "GFAP", "ALDH1L1"]
CELL_TYPE_SORT_ORDER = {cell_type: index for index, cell_type in enumerate(CELL_TYPE_ORDER)}


def parse_args() -> argparse.Namespace:
    script = Path(__file__).resolve()
    package = script.parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path, default=package)
    return parser.parse_args()


def parse_library_name(path: Path) -> tuple[str, str, str]:
    match = re.match(
        r"(?P<gsm>GSM\d+)_(?P<region>fovea|peripheral)_donor_(?P<donor>\d+)_expression\.tsv\.gz",
        path.name,
    )
    if match is None:
        raise ValueError(f"Unexpected GSE130636 filename: {path.name}")
    return match.group("gsm"), match.group("region"), f"donor_{match.group('donor')}"


def validate_author_mapping(mapping: pd.DataFrame) -> None:
    expected = {str(i) for i in range(1, 18)} | {"8A", "8B"}
    expected.remove("8")
    observed = set(mapping["cluster_label"].astype(str))
    if observed != expected:
        raise ValueError(
            f"Author mapping mismatch. Missing={sorted(expected-observed)}, "
            f"unexpected={sorted(observed-expected)}"
        )
    if mapping["cluster_label"].duplicated().any():
        raise ValueError("Author mapping contains duplicate cluster labels")
    if not mapping["source_evidence"].str.contains("Voigt et al. 2019").all():
        raise ValueError("Every cluster mapping must retain author-source evidence")


def localization_panel_from_discovery_summary(summary: dict[str, object]) -> list[str]:
    """Return the descriptive atlas panel from the current discovery analysis.

    The panel is the ordered union of the total-association and DME-conditioned
    top-five sets. Atlas values do not contribute to ranking; deriving the panel
    from the current summary keeps the descriptive display aligned with both
    prespecified estimands.
    """

    total = summary.get("total_top5")
    conditioned = summary.get("dme_conditioned_top5")
    if not isinstance(total, list) or len(total) != 5:
        raise ValueError(
            "final_discovery_summary.json must contain a five-gene "
            "total_top5 list"
        )
    if not isinstance(conditioned, list) or len(conditioned) != 5:
        raise ValueError(
            "final_discovery_summary.json must contain a five-gene "
            "dme_conditioned_top5 list"
        )
    genes = list(
        dict.fromkeys(
            [str(gene).upper() for gene in [*total, *conditioned]]
        )
    )
    if len(genes) != 6 or "P2RX4" not in genes:
        raise ValueError(
            "The normal-retina panel must be the six-gene union including P2RX4"
        )
    return genes


def summarize_library(
    path: Path,
    genes: list[str],
    mapping: pd.DataFrame,
) -> pd.DataFrame:
    gsm, region, donor = parse_library_name(path)
    header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    missing = [gene for gene in genes if gene not in header]
    if missing:
        raise ValueError(f"{path.name} is missing candidate genes: {missing}")
    data = pd.read_csv(
        path,
        sep="\t",
        usecols=["barcode", "cluster_label", *genes],
        dtype={"barcode": str, "cluster_label": str},
    )
    map_by_cluster = mapping.set_index("cluster_label")["analysis_cell_type"]
    data["cell_type"] = data["cluster_label"].map(map_by_cluster)
    if data["cell_type"].isna().any():
        unknown = sorted(data.loc[data["cell_type"].isna(), "cluster_label"].unique())
        raise ValueError(f"Unmapped author cluster labels in {path.name}: {unknown}")
    data = data.loc[data["cell_type"] != "Unknown - excluded"].copy()

    rows = []
    for cell_type, group in data.groupby("cell_type", sort=False):
        for gene in genes:
            values = pd.to_numeric(group[gene], errors="coerce").fillna(0.0).to_numpy()
            linearized = np.expm1(values)
            rows.append(
                {
                    "geo_accession": gsm,
                    "donor": donor,
                    "region": region,
                    "cell_type": cell_type,
                    "gene_symbol": gene,
                    "n_cells": int(len(values)),
                    "mean_log_normalized_expression": float(values.mean()),
                    "median_log_normalized_expression": float(np.median(values)),
                    "mean_expm1_log_normalized_expression": float(linearized.mean()),
                    "median_expm1_log_normalized_expression": float(
                        np.median(linearized)
                    ),
                    "detection_fraction": float(np.mean(values > 0)),
                }
            )
    return pd.DataFrame(rows)


def aggregate_cell_types(pseudobulk: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (cell_type, gene), group in pseudobulk.groupby(["cell_type", "gene_symbol"]):
        weights = group["n_cells"].to_numpy(dtype=float)
        detection = group["detection_fraction"].to_numpy(dtype=float)
        rows.append(
            {
                "cell_type": cell_type,
                "gene_symbol": gene,
                "n_libraries": int(len(group)),
                "total_cells": int(group["n_cells"].sum()),
                "donor_balanced_mean_expm1_expression": float(
                    group["mean_expm1_log_normalized_expression"].mean()
                ),
                "donor_balanced_median_expm1_expression": float(
                    group["mean_expm1_log_normalized_expression"].median()
                ),
                "pooled_detection_fraction": float(np.average(detection, weights=weights)),
                "min_cells_per_library": int(group["n_cells"].min()),
            }
        )
    out = pd.DataFrame(rows)
    out["expression_relative_to_gene_max"] = out.groupby("gene_symbol")[
        "donor_balanced_mean_expm1_expression"
    ].transform(lambda x: x / x.max() if x.max() else 0.0)
    specificity = {}
    for gene, group in out.groupby("gene_symbol"):
        values = group["donor_balanced_mean_expm1_expression"].to_numpy(dtype=float)
        maximum = float(values.max())
        specificity[gene] = (
            float(np.sum(1 - values / maximum) / (len(values) - 1))
            if maximum > 0 and len(values) > 1
            else 0.0
        )
    out["tau_cell_type_specificity"] = out["gene_symbol"].map(specificity)
    return out


def _dominance_summary(
    data: pd.DataFrame,
    group_columns: list[str],
    expression_column: str,
    detection_column: str,
) -> pd.DataFrame:
    """Return one deterministic dominant cell type per analytical unit and gene."""
    rows = []
    for keys, group in data.groupby(group_columns, sort=True, observed=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        ranked = group.assign(
            _cell_type_order=group["cell_type"].map(CELL_TYPE_SORT_ORDER).fillna(999)
        ).sort_values(
            [expression_column, detection_column, "_cell_type_order"],
            ascending=[False, False, True],
        )
        top = ranked.iloc[0]
        second_expression = (
            float(ranked.iloc[1][expression_column]) if len(ranked) > 1 else 0.0
        )
        top_expression = float(top[expression_column])
        row = dict(zip(group_columns, keys, strict=True))
        row.update(
            {
                "dominant_cell_type": str(top["cell_type"]),
                "dominant_mean_expm1_expression": top_expression,
                "dominant_detection_fraction": float(top[detection_column]),
                "second_mean_expm1_expression": second_expression,
                "dominance_margin": top_expression - second_expression,
                "n_cell_types_compared": int(len(ranked)),
                "dominance_identifiable": bool(top_expression > 0),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def library_dominance_summary(pseudobulk: pd.DataFrame) -> pd.DataFrame:
    """Summarise the dominant candidate localization in each donor-region library."""
    return _dominance_summary(
        pseudobulk,
        ["geo_accession", "donor", "region", "gene_symbol"],
        "mean_expm1_log_normalized_expression",
        "detection_fraction",
    )


def donor_aggregated_localization(pseudobulk: pd.DataFrame) -> pd.DataFrame:
    """Average regions within donors before comparing normal-retina cell types."""
    donor = (
        pseudobulk.groupby(
            ["donor", "cell_type", "gene_symbol"],
            as_index=False,
            sort=True,
            observed=True,
        )
        .agg(
            n_regions_observed=("region", "nunique"),
            n_libraries=("geo_accession", "nunique"),
            total_cells=("n_cells", "sum"),
            region_balanced_mean_expm1_expression=(
                "mean_expm1_log_normalized_expression",
                "mean",
            ),
            region_balanced_mean_detection_fraction=("detection_fraction", "mean"),
        )
    )
    donor["_cell_type_order"] = donor["cell_type"].map(CELL_TYPE_SORT_ORDER).fillna(999)
    donor = donor.sort_values(
        [
            "donor",
            "gene_symbol",
            "region_balanced_mean_expm1_expression",
            "region_balanced_mean_detection_fraction",
            "_cell_type_order",
        ],
        ascending=[True, True, False, False, True],
    )
    donor["cell_type_rank_within_donor"] = (
        donor.groupby(["donor", "gene_symbol"], observed=True).cumcount() + 1
    )
    maxima = donor.groupby(["donor", "gene_symbol"], observed=True)[
        "region_balanced_mean_expm1_expression"
    ].transform("max")
    donor["expression_relative_to_donor_gene_max"] = np.where(
        maxima > 0,
        donor["region_balanced_mean_expm1_expression"] / maxima,
        0.0,
    )
    donor["is_dominant_cell_type"] = donor["cell_type_rank_within_donor"] == 1
    return donor.drop(columns="_cell_type_order").reset_index(drop=True)


def leave_one_donor_dominance(donor_localization: pd.DataFrame) -> pd.DataFrame:
    """Recompute dominant cell types after omitting each normal-retina donor."""
    expression = "region_balanced_mean_expm1_expression"
    detection = "region_balanced_mean_detection_fraction"
    full = (
        donor_localization.groupby(["cell_type", "gene_symbol"], as_index=False)
        .agg(
            donor_balanced_mean_expm1_expression=(expression, "mean"),
            donor_balanced_mean_detection_fraction=(detection, "mean"),
        )
    )
    full_dominance = _dominance_summary(
        full,
        ["gene_symbol"],
        "donor_balanced_mean_expm1_expression",
        "donor_balanced_mean_detection_fraction",
    )[["gene_symbol", "dominant_cell_type"]].rename(
        columns={"dominant_cell_type": "full_cohort_dominant_cell_type"}
    )

    folds = []
    donors = sorted(donor_localization["donor"].unique())
    for omitted_donor in donors:
        retained = donor_localization.loc[donor_localization["donor"] != omitted_donor]
        pooled = (
            retained.groupby(["cell_type", "gene_symbol"], as_index=False)
            .agg(
                retained_donor_mean_expm1_expression=(expression, "mean"),
                retained_donor_mean_detection_fraction=(detection, "mean"),
            )
        )
        dominance = _dominance_summary(
            pooled,
            ["gene_symbol"],
            "retained_donor_mean_expm1_expression",
            "retained_donor_mean_detection_fraction",
        )
        dominance.insert(0, "omitted_donor", omitted_donor)
        dominance.insert(1, "n_donors_retained", int(len(donors) - 1))
        folds.append(dominance)

    result = pd.concat(folds, ignore_index=True).merge(
        full_dominance, on="gene_symbol", how="left", validate="many_to_one"
    )
    result["matches_full_cohort_dominance"] = (
        result["dominant_cell_type"] == result["full_cohort_dominant_cell_type"]
    )
    return result.sort_values(["gene_symbol", "omitted_donor"]).reset_index(drop=True)


def donor_interval_summary(donor_localization: pd.DataFrame) -> pd.DataFrame:
    """Provide descriptive between-donor ranges without treating three donors as cells."""
    expression = "region_balanced_mean_expm1_expression"
    detection = "region_balanced_mean_detection_fraction"
    return (
        donor_localization.groupby(["cell_type", "gene_symbol"], as_index=False)
        .agg(
            n_donors=("donor", "nunique"),
            donor_mean_expm1_expression=(expression, "mean"),
            donor_sd_expm1_expression=(expression, "std"),
            donor_min_expm1_expression=(expression, "min"),
            donor_max_expm1_expression=(expression, "max"),
            donor_mean_detection_fraction=(detection, "mean"),
            donor_sd_detection_fraction=(detection, "std"),
            donor_min_detection_fraction=(detection, "min"),
            donor_max_detection_fraction=(detection, "max"),
        )
        .sort_values(["gene_symbol", "donor_mean_expm1_expression"], ascending=[True, False])
        .reset_index(drop=True)
    )


def paired_region_summary(pseudobulk: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (cell_type, gene), group in pseudobulk.groupby(["cell_type", "gene_symbol"]):
        pivot = group.pivot_table(
            index="donor",
            columns="region",
            values="mean_expm1_log_normalized_expression",
            aggfunc="first",
        )
        if not {"fovea", "peripheral"}.issubset(pivot.columns):
            continue
        paired = pivot.dropna(subset=["fovea", "peripheral"]).copy()
        differences = paired["fovea"] - paired["peripheral"]
        rows.append(
            {
                "cell_type": cell_type,
                "gene_symbol": gene,
                "n_paired_donors": int(len(paired)),
                "mean_fovea_expression": float(paired["fovea"].mean()),
                "mean_peripheral_expression": float(paired["peripheral"].mean()),
                "mean_paired_fovea_minus_peripheral": float(differences.mean()),
                "median_paired_fovea_minus_peripheral": float(differences.median()),
                "positive_paired_differences": int((differences > 0).sum()),
            }
        )
    return pd.DataFrame(rows)


def top_localization_table(aggregated: pd.DataFrame) -> pd.DataFrame:
    ranked = aggregated.sort_values(
        ["gene_symbol", "donor_balanced_mean_expm1_expression"],
        ascending=[True, False],
    ).copy()
    ranked["cell_type_rank"] = ranked.groupby("gene_symbol").cumcount() + 1
    return ranked.loc[ranked["cell_type_rank"] <= 3].reset_index(drop=True)


def glial_cluster_markers(
    files: list[Path], mapping: pd.DataFrame, genes: list[str]
) -> pd.DataFrame:
    marker_genes = list(dict.fromkeys([*GLIAL_MARKERS, *genes]))
    rows = []
    for path in files:
        gsm, region, donor = parse_library_name(path)
        header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
        available = [gene for gene in marker_genes if gene in header]
        data = pd.read_csv(
            path,
            sep="\t",
            usecols=["cluster_label", *available],
            dtype={"cluster_label": str},
        )
        data = data.loc[data["cluster_label"].isin(["13", "14", "15", "16", "17"])]
        for cluster, group in data.groupby("cluster_label"):
            for gene in available:
                values = pd.to_numeric(group[gene], errors="coerce").fillna(0).to_numpy()
                rows.append(
                    {
                        "geo_accession": gsm,
                        "donor": donor,
                        "region": region,
                        "cluster_label": cluster,
                        "author_label": "Glial cells",
                        "analysis_label": "Müller-enriched glia",
                        "label_basis": (
                            "author label retained; analysis label supported by "
                            "RLBP1/GFAP/ALDH1L1 marker profile"
                        ),
                        "gene_symbol": gene,
                        "gene_role": (
                            "glial_identity_marker" if gene in GLIAL_MARKERS else "candidate"
                        ),
                        "n_cells": len(values),
                        "mean_expm1_log_normalized_expression": float(
                            np.expm1(values).mean()
                        ),
                        "detection_fraction": float(np.mean(values > 0)),
                    }
                )
    table = pd.DataFrame(rows)
    return (
        table.groupby(
            [
                "cluster_label",
                "author_label",
                "analysis_label",
                "label_basis",
                "gene_symbol",
                "gene_role",
            ],
            as_index=False,
        )
        .agg(
            n_libraries=("geo_accession", "size"),
            total_cells=("n_cells", "sum"),
            donor_region_balanced_mean_expm1_expression=(
                "mean_expm1_log_normalized_expression",
                "mean",
            ),
            pooled_detection_fraction=(
                "detection_fraction",
                "mean",
            ),
        )
    )


def draw_figure(
    aggregated: pd.DataFrame,
    regional: pd.DataFrame,
    lodo_dominance: pd.DataFrame,
    genes: list[str],
    output: Path,
) -> None:
    data = aggregated.copy()
    data["cell_type"] = pd.Categorical(
        data["cell_type"], categories=CELL_TYPE_ORDER, ordered=True
    )
    data = data.sort_values(["cell_type", "gene_symbol"])

    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10})
    fig = plt.figure(figsize=(13.4, 9.2), constrained_layout=True)
    fig.suptitle(
        "Normal-retina localization context and donor stability (GSE130636)",
        fontsize=13,
        fontweight="bold",
    )
    grid = fig.add_gridspec(2, 2, width_ratios=[1.25, 1.0], height_ratios=[1.0, 0.85])

    ax = fig.add_subplot(grid[0, 0])
    for y, cell_type in enumerate(CELL_TYPE_ORDER):
        for x, gene in enumerate(genes):
            row = data.loc[(data["cell_type"] == cell_type) & (data["gene_symbol"] == gene)]
            if row.empty:
                continue
            size = 30 + 520 * float(row["pooled_detection_fraction"].iloc[0])
            color = float(row["expression_relative_to_gene_max"].iloc[0])
            ax.scatter(x, y, s=size, c=[color], vmin=0, vmax=1, cmap="viridis", edgecolor="white")
    ax.set_xticks(range(len(genes)), genes, rotation=35, ha="right")
    ax.set_yticks(range(len(CELL_TYPE_ORDER)), CELL_TYPE_ORDER)
    ax.invert_yaxis()
    ax.set_xlim(-0.6, len(genes) - 0.4)
    ax.grid(color="#dddddd", linewidth=0.5)
    ax.text(-0.12, 1.04, "a", transform=ax.transAxes, fontsize=14, fontweight="bold")
    for fraction, x_position in zip([0.1, 0.3, 0.5], [0.15, 0.38, 0.64]):
        ax.scatter(
            x_position,
            -0.18,
            s=30 + 520 * fraction,
            color="#6a9d79",
            edgecolor="white",
            transform=ax.transAxes,
            clip_on=False,
        )
        ax.text(
            x_position + 0.045,
            -0.18,
            f"{fraction:.0%}",
            transform=ax.transAxes,
            va="center",
            fontsize=8,
        )
    ax.text(0.15, -0.28, "cell detection fraction", transform=ax.transAxes, fontsize=8)

    ax = fig.add_subplot(grid[0, 1])
    gene_summary = (
        data.sort_values(
            ["gene_symbol", "donor_balanced_mean_expm1_expression"],
            ascending=[True, False],
        )
        .groupby("gene_symbol", as_index=False)
        .first()
        .set_index("gene_symbol")
        .reindex(genes)
    )
    ax.bar(
        genes,
        gene_summary["tau_cell_type_specificity"],
        color="#4f7cac",
        width=0.68,
    )
    for index, row in enumerate(gene_summary.itertuples()):
        ax.text(
            index,
            row.tau_cell_type_specificity + 0.025,
            str(row.cell_type).replace(" cells", ""),
            ha="center",
            va="bottom",
            fontsize=8.5,
            rotation=25,
        )
    ax.set_ylim(0, 1.14)
    ax.set_ylabel("Tau cell-type specificity")
    ax.tick_params(axis="x", rotation=35)
    ax.set_title("Pooled normal-retina localization context", fontsize=10, pad=14)
    ax.text(-0.12, 1.04, "b", transform=ax.transAxes, fontsize=14, fontweight="bold")

    ax = fig.add_subplot(grid[1, 0])
    regional_matrix = regional.pivot(
        index="cell_type", columns="gene_symbol", values="mean_paired_fovea_minus_peripheral"
    ).reindex(index=CELL_TYPE_ORDER, columns=genes)
    bound = float(np.nanmax(np.abs(regional_matrix.to_numpy())))
    bound = bound if bound > 0 else 1.0
    image = ax.imshow(regional_matrix, cmap="coolwarm", vmin=-bound, vmax=bound, aspect="auto")
    ax.set_xticks(range(len(genes)), genes, rotation=35, ha="right")
    ax.set_yticks(range(len(CELL_TYPE_ORDER)), CELL_TYPE_ORDER)
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title("Paired fovea - peripheral expression", fontsize=10)
    ax.text(-0.12, 1.04, "c", transform=ax.transAxes, fontsize=14, fontweight="bold")

    ax = fig.add_subplot(grid[1, 1])
    stability = (
        lodo_dominance.groupby("gene_symbol", as_index=False)
        .agg(
            matching_folds=("matches_full_cohort_dominance", "sum"),
            n_lodo_folds=("omitted_donor", "nunique"),
            full_cohort_dominant_cell_type=(
                "full_cohort_dominant_cell_type",
                "first",
            ),
        )
        .set_index("gene_symbol")
        .reindex(genes)
    )
    stability["agreement_fraction"] = (
        stability["matching_folds"] / stability["n_lodo_folds"]
    )
    colors = [
        "#4f7cac" if value == 1.0 else "#c9973b"
        for value in stability["agreement_fraction"]
    ]
    ax.bar(genes, stability["agreement_fraction"], color=colors, width=0.68)
    for index, row in enumerate(stability.itertuples()):
        ax.text(
            index,
            row.agreement_fraction + 0.035,
            f"{int(row.matching_folds)}/{int(row.n_lodo_folds)}",
            ha="center",
            va="bottom",
            fontsize=8,
        )
        ax.text(
            index,
            0.04,
            str(row.full_cohort_dominant_cell_type).replace(" cells", ""),
            ha="center",
            va="bottom",
            fontsize=8.5,
            rotation=90,
            color="#333333",
        )
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("LODO agreement with pooled dominance")
    ax.tick_params(axis="x", rotation=35)
    ax.axhline(1.0, color="#555555", linewidth=0.8, linestyle="--")
    ax.set_title("Normal-donor stability (three leave-one-donor folds)", fontsize=10)
    ax.text(-0.12, 1.04, "d", transform=ax.transAxes, fontsize=14, fontweight="bold")

    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    package = args.package_root.resolve()
    source = package / "analysis_data" / "external_single_cell"
    results = package / "analysis_results"
    figures = package / "figures"
    results.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    summary = json.loads(
        (results / "final_discovery_summary.json").read_text(encoding="utf-8")
    )
    genes = localization_panel_from_discovery_summary(summary)
    mapping = pd.read_csv(source / "voigt2019_author_cluster_mapping.csv", dtype=str)
    validate_author_mapping(mapping)

    files = sorted(source.glob("GSM*_expression.tsv.gz"))
    if len(files) != 6:
        raise ValueError(f"Expected six GSE130636 matrices, found {len(files)}")
    pseudobulk = pd.concat(
        [summarize_library(path, genes, mapping) for path in files], ignore_index=True
    )
    aggregated = aggregate_cell_types(pseudobulk)
    regional = paired_region_summary(pseudobulk)
    top = top_localization_table(aggregated)
    library_dominance = library_dominance_summary(pseudobulk)
    donor_localization = donor_aggregated_localization(pseudobulk)
    lodo_dominance = leave_one_donor_dominance(donor_localization)
    donor_intervals = donor_interval_summary(donor_localization)

    mapping.to_csv(results / "normal_retina_author_cluster_mapping.csv", index=False)
    library_dominance.to_csv(
        results / "normal_retina_library_dominance.csv", index=False
    )
    donor_localization.to_csv(
        results / "normal_retina_donor_aggregated_localization.csv", index=False
    )
    lodo_dominance.to_csv(
        results / "normal_retina_leave_one_donor_dominance.csv", index=False
    )
    donor_intervals.to_csv(
        results / "normal_retina_donor_interval_summary.csv", index=False
    )

    draw_figure(
        aggregated,
        regional,
        lodo_dominance,
        genes,
        figures / "Figure_5_cell_type_localization",
    )

    dominant = (
        top.loc[top["cell_type_rank"] == 1]
        .set_index("gene_symbol")
        .reindex(genes)
        .reset_index()
    )
    output_summary = {
        "dataset": "GSE130636",
        "source_article": "Voigt et al., Experimental Eye Research 2019",
        "source_mapping_provenance": (
            "analysis_data/external_single_cell/AUTHOR_MAPPING_PROVENANCE.md"
        ),
        "source_mapping_page": 38,
        "n_cells_total": 8217,
        "n_donors": 3,
        "n_regions_per_donor": 2,
        "n_libraries": 6,
        "author_clusters": 17,
        "analysis_cell_types": len(CELL_TYPE_ORDER),
        "excluded_author_cluster": "9 (Unknown)",
        "candidate_genes": genes,
        "evidence_scope": "normal-retina localization context",
        "dominant_cell_types": {
            row.gene_symbol: {
                "cell_type": row.cell_type,
                "mean_expm1_log_normalized_expression": (
                    row.donor_balanced_mean_expm1_expression
                ),
                "expression_relative_to_gene_max": row.expression_relative_to_gene_max,
                "tau_cell_type_specificity": row.tau_cell_type_specificity,
                "detection_fraction": row.pooled_detection_fraction,
            }
            for row in dominant.itertuples()
        },
        "leave_one_donor_stability": {
            gene: {
                "full_cohort_dominant_cell_type": str(
                    lodo_dominance.loc[
                        lodo_dominance["gene_symbol"] == gene,
                        "full_cohort_dominant_cell_type",
                    ].iloc[0]
                ),
                "matching_folds": int(
                    lodo_dominance.loc[
                        lodo_dominance["gene_symbol"] == gene,
                        "matches_full_cohort_dominance",
                    ].sum()
                ),
                "n_folds": int(
                    lodo_dominance.loc[
                        lodo_dominance["gene_symbol"] == gene,
                        "omitted_donor",
                    ].nunique()
                ),
            }
            for gene in genes
        },
        "figure_title": (
            "Normal-retina localization context and donor stability (GSE130636)"
        ),
        "figure_caption": (
            "Author-mapped normal-retina donor-region libraries provide descriptive "
            "localization context. Expression and detection are summarized within "
            "library and cell type, aggregated across regions within donor, and "
            "checked in three leave-one-donor folds. These data do not establish "
            "disease-state cellular origin or independent DR replication."
        ),
        "interpretation": (
            "Normal-retina localization context only; not disease-state cellular "
            "origin or disease replication. Donor-region libraries, not individual "
            "cells, were the aggregation units. Between-donor intervals are "
            "descriptive ranges across three donors, not population confidence "
            "intervals. No log-expression percentage or compositional share was "
            "calculated."
        ),
    }
    (results / "single_cell_localization_summary.json").write_text(
        json.dumps(output_summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(output_summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
