#!/usr/bin/env python3
"""Generate the current CardioQueue study workflow and claim boundary."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[2]
OUTPUTS = [
    ROOT / "results/manuscript_figures_tables/figures/cardioqueue_workflow.png",
    ROOT / "paper_MS_SI/figures/Figure_1_cardioqueue_workflow.png",
]
INK = "#17262C"
MUTED = "#526168"


def stage(ax, x, y, width, height, number, title, lines, face, edge):
    ax.add_patch(FancyBboxPatch((x, y), width, height,
        boxstyle="round,pad=0.010,rounding_size=0.012", linewidth=1.5,
        edgecolor=edge, facecolor=face))
    ax.text(x + 0.018, y + height - 0.030, number, fontsize=9.5, fontweight="bold",
            color="white", ha="center", va="center",
            bbox=dict(boxstyle="circle,pad=0.28", facecolor=edge, edgecolor="none"))
    ax.text(x + 0.047, y + height - 0.030, title, fontsize=11.2,
            fontweight="bold", color=INK, va="center")
    ax.text(x + 0.018, y + height - 0.073, "\n".join(lines), fontsize=8.7,
            color="#2D3A40", va="top", linespacing=1.30)


def arrow(ax, start, end):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=13,
                                linewidth=1.5, color="#5B686E"))


def main() -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    fig, ax = plt.subplots(figsize=(13.4, 7.5))
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")

    ax.text(0.025, 0.963, "CardioQueue: from multi-source evidence to a bounded review queue",
            fontsize=18, fontweight="bold", color=INK, va="top")
    ax.text(0.025, 0.917,
            "Exact-allele quarantine, multi-source evaluation, and retrospective VUS prioritization were kept distinct.",
            fontsize=10.5, color=MUTED, va="top")

    y_top, y_bottom, height, width = 0.615, 0.315, 0.225, 0.292
    xs = [0.025, 0.354, 0.683]
    stage(ax, xs[0], y_top, width, height, "1", "Assemble and harmonize",
          ["85,968 GRCh38 registry rows", "ClinVar + five held-out sources", "Variant, gene, protein and structure data", "Source assertions retained separately"],
          "#E8F1F2", "#26747C")
    stage(ax, xs[1], y_top, width, height, "2", "Quarantine before fitting",
          ["1,233 candidate benchmark alleles", "Zero exact-allele development overlap", "35,796 train; 6,317 validation", "Outcome-blind source and overlap rules"],
          "#EEF2E7", "#657B46")
    stage(ax, xs[2], y_top, width, height, "3", "Fit primary model",
          ["CatBoost with 300 public fields", "Population, consequence and effect scores", "Protein context and structure availability", "CardioQueue score; no private HiRO fields"],
          "#F8EFE7", "#AD642E")
    arrow(ax, (xs[0] + width, y_top + height / 2), (xs[1], y_top + height / 2))
    arrow(ax, (xs[1] + width, y_top + height / 2), (xs[2], y_top + height / 2))

    stage(ax, xs[0], y_bottom, width, height, "6", "Translate score to review order",
          ["0.1-0.9 explicit deferral interval", "25,481 historical VUS ranked", "Review-budget recovery measured", "Prospective clinical utility remains untested"],
          "#F7ECE8", "#A94B3A")
    stage(ax, xs[1], y_bottom, width, height, "5", "Test robustness",
          ["CardioBoost head-to-head benchmark", "Calibration, leakage and subgroup audits", "4,137 FoldX stability calculations", "Numeric DDG did not add robust benefit"],
          "#EAF0F6", "#486F92")
    stage(ax, xs[2], y_bottom, width, height, "4", "Evaluate held-out alleles",
          ["775 unique variants across 40 genes", "HiRO, eMERGE, CardioBoost, PLOS, SHaRe", "436 P/LP; 339 B/LB", "Source results + cross-source consensus"],
          "#F0ECF4", "#705B7C")
    arrow(ax, (xs[2] + width / 2, y_top), (xs[2] + width / 2, y_bottom + height))
    arrow(ax, (xs[2], y_bottom + height / 2), (xs[1] + width, y_bottom + height / 2))
    arrow(ax, (xs[1], y_bottom + height / 2), (xs[0] + width, y_bottom + height / 2))

    ax.text(0.025, 0.257, "Clinical interpretation boundary", fontsize=12.5,
            fontweight="bold", color=INK)
    ax.add_patch(FancyBboxPatch((0.025, 0.065), 0.950, 0.155,
        boxstyle="round,pad=0.010,rounding_size=0.012", facecolor="#F4F6F6",
        edgecolor="#8C989D", linewidth=1.3))
    columns = [
        (0.047, "MODEL OUTPUT", "Score and evidence-review priority", "#26747C"),
        (0.365, "EXPERT INTERPRETATION", "Phenotype, inheritance, segregation,\nmechanism and ACMG/AMP evidence", "#657B46"),
        (0.690, "ONLY IF INDEPENDENTLY SUPPORTED", "Amended laboratory classification\nand any resulting management", "#A94B3A"),
    ]
    for x, heading, body, color in columns:
        ax.text(x, 0.184, heading, fontsize=8.4, fontweight="bold", color=color, va="top")
        ax.text(x, 0.143, body, fontsize=9.4, color=INK, va="top", linespacing=1.28)
    ax.plot([0.338, 0.338], [0.085, 0.195], color="#CAD0D2", lw=1)
    ax.plot([0.665, 0.665], [0.085, 0.195], color="#CAD0D2", lw=1)
    arrow(ax, (0.310, 0.135), (0.355, 0.135))
    arrow(ax, (0.637, 0.135), (0.680, 0.135))
    ax.text(0.365, 0.083, "Clinician override: urgent or management-sensitive cases bypass score order.",
            fontsize=7.8, color="#A94B3A", va="center")

    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    for output in OUTPUTS:
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
        fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
