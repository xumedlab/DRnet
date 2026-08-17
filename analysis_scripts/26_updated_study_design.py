#!/usr/bin/env python3
"""Draw the final discovery, validation, and interpretation study design."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


PACKAGE = Path(__file__).resolve().parents[1]
OUTPUT = PACKAGE / "figures" / "Figure_1_final_study_design"


def box(ax, xy, width, height, label, color, fontsize=10):
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.018,rounding_size=0.025",
        linewidth=1.4,
        edgecolor="#40515b",
        facecolor=color,
    )
    ax.add_patch(patch)
    ax.text(
        xy[0] + width / 2,
        xy[1] + height / 2,
        label,
        ha="center",
        va="center",
        fontsize=fontsize,
        linespacing=1.25,
    )


def arrow(ax, start, end):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=14,
            linewidth=1.4,
            color="#40515b",
        )
    )


def main():
    fig, ax = plt.subplots(figsize=(12.2, 7.4))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.02, 0.96, "a", fontsize=18, fontweight="bold", va="top")
    ax.text(
        0.07,
        0.96,
        "Discovery cohort and analysis unit",
        fontsize=14,
        fontweight="bold",
        va="top",
    )
    box(
        ax,
        (0.05, 0.70),
        0.25,
        0.17,
        "GSE160306 macular RNA-seq\n39 eyes from 36 donors\n26 donors with diabetes",
        "#d9ebf5",
    )
    box(
        ax,
        (0.37, 0.70),
        0.25,
        0.17,
        "Donor aggregation\npaired-eye expression averaged\nseverity retained at donor level",
        "#d9ebf5",
    )
    box(
        ax,
        (0.69, 0.70),
        0.25,
        0.17,
        "Predefined inflammatory space\n158 expressed Hallmark genes\nno outcome-driven expansion",
        "#d9ebf5",
    )
    arrow(ax, (0.30, 0.785), (0.37, 0.785))
    arrow(ax, (0.62, 0.785), (0.69, 0.785))

    ax.text(0.02, 0.61, "b", fontsize=18, fontweight="bold", va="top")
    ax.text(
        0.07,
        0.61,
        "Target selection and independent disease-state checks",
        fontsize=14,
        fontweight="bold",
        va="top",
    )
    box(
        ax,
        (0.05, 0.40),
        0.25,
        0.14,
        "Primary total-association model\nseverity + age + sex + PMI + RIN\nHC3+t and wild bootstrap-t",
        "#dcefe4",
        fontsize=9.3,
    )
    box(
        ax,
        (0.37, 0.40),
        0.25,
        0.14,
        "DME-conditioned model\n2,000 donor bootstraps\nLODO and influence analyses",
        "#dcefe4",
        fontsize=9.3,
    )
    box(
        ax,
        (0.69, 0.40),
        0.25,
        0.14,
        "Single-target external analysis\nP2RX4 only\nlocal protocol SHA-256 recorded",
        "#f7d9d3",
        fontsize=9.3,
    )
    arrow(ax, (0.30, 0.47), (0.37, 0.47))
    arrow(ax, (0.62, 0.47), (0.69, 0.47))

    ax.text(0.02, 0.30, "c", fontsize=18, fontweight="bold", va="top")
    ax.text(
        0.07,
        0.30,
        "Independent human ocular cohorts and localization context",
        fontsize=14,
        fontweight="bold",
        va="top",
    )
    box(
        ax,
        (0.05, 0.075),
        0.26,
        0.14,
        "GSE276892\n8 PDR vs 9 surgical controls\nvitreous hyalocytes\nraw-read reconstruction + QC",
        "#eee6f5",
        fontsize=9.0,
    )
    box(
        ax,
        (0.37, 0.075),
        0.26,
        0.14,
        "GSE179568\n7 PDR RNV vs 10 membranes\nseparate GEO dataset\nage/treatment sensitivity",
        "#eee6f5",
        fontsize=9.0,
    )
    box(
        ax,
        (0.69, 0.075),
        0.26,
        0.14,
        "GSE130636\n3 normal-retina donors\nauthor-mapped localization context\nnot disease-state replication",
        "#e0eddc",
        fontsize=9.0,
    )

    fig.savefig(OUTPUT.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(OUTPUT.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
