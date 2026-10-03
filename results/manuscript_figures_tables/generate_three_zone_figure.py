#!/usr/bin/env python3
"""Generate the CardioBoost-style 0.1/0.9 CardioQueue VUS triage figure."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TRIAGE = ROOT / "results/model_performance/final_cardioqueue_three_zone_all_rows"
HELDOUT = ROOT / "results/model_performance/final_model_heldout_2026/combined_unique_predictions.tsv"
OUTPUTS = [
    ROOT / "results/manuscript_figures_tables/figures/vus_three_zone_triage_distribution.png",
    ROOT / "paper_MS_SI/figures/Figure_3_vus_three_zone_triage.png",
]
COLORS = {"Lower-priority review": "#3977A8", "Deferred review": "#E2B83F", "Accelerated review": "#B84A4A"}
INK = "#18262C"


def read_vus() -> pd.DataFrame:
    return pd.read_csv(TRIAGE / "all_registry_vus_scores.tsv", sep="\t", low_memory=False)


def main() -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    all_vus = read_vus()
    heldout = pd.read_csv(HELDOUT, sep="\t", low_memory=False)
    sources = [
        ("All registry", all_vus),
        ("ClinVar", all_vus.loc[all_vus["in_clinvar"].astype(bool)]),
        ("eMERGE", all_vus.loc[all_vus["in_emerge"].astype(bool)]),
        ("HiRO", all_vus.loc[all_vus["in_hiro"].astype(bool)]),
    ]

    fig = plt.figure(figsize=(13.2, 8.4), facecolor="white")
    grid = fig.add_gridspec(2, 2, height_ratios=[1.04, 1], hspace=0.43, wspace=0.28)
    ax1 = fig.add_subplot(grid[0, :])
    ax2 = fig.add_subplot(grid[1, 0])
    ax3 = fig.add_subplot(grid[1, 1])

    ax1.axvspan(0, 0.1, color=COLORS["Lower-priority review"], alpha=0.13, lw=0)
    ax1.axvspan(0.1, 0.9, color=COLORS["Deferred review"], alpha=0.15, lw=0)
    ax1.axvspan(0.9, 1, color=COLORS["Accelerated review"], alpha=0.13, lw=0)
    bins = np.linspace(0, 1, 31)
    for label, color in [("Benign", COLORS["Lower-priority review"]), ("Pathogenic", COLORS["Accelerated review"])]:
        values = heldout.loc[heldout["combined_label"].eq(label), "pathogenic_probability"]
        weights = np.full(len(values), 100 / len(values))
        ax1.hist(values, bins=bins, weights=weights, histtype="step", linewidth=2.0,
                 color=color, label=f"{'B/LB' if label == 'Benign' else 'P/LP'} reference alleles ({len(values):,})")
    for threshold in (0.1, 0.9):
        ax1.axvline(threshold, color=INK, linestyle=(0, (4, 3)), lw=1.3)
    ymax = ax1.get_ylim()[1]
    ax1.text(0.05, ymax * 0.94, "Lower priority\n(benign-like) ≤0.1", ha="center", va="top",
             color="#245E88", fontweight="bold")
    ax1.text(0.50, ymax * 0.94, "Deferred / standard review\n0.1-0.9", ha="center", va="top",
             color="#80600E", fontweight="bold")
    ax1.text(0.95, ymax * 0.94, "Accelerated review\n(pathogenic-like) ≥0.9", ha="center", va="top",
             color="#8F3535", fontweight="bold")
    ax1.set(xlim=(0, 1), xlabel="CardioQueue score (higher = greater similarity to development P/LP variants)",
            ylabel="Within-class variants per bin (%)")
    ax1.set_title(f"A  Score separation in {len(heldout):,} held-out alleles",
                  loc="left", fontsize=12.5, fontweight="bold", color=INK)
    ax1.legend(frameon=False, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.82))
    ax1.grid(axis="y", color="#D9DEDF", lw=0.7, alpha=0.8)
    ax1.set_axisbelow(True)

    heldout_groups = [("True B/LB", heldout.loc[heldout["y_true"].eq(0)]),
                      ("True P/LP", heldout.loc[heldout["y_true"].eq(1)])]
    y = np.arange(len(heldout_groups))
    left = np.zeros(len(heldout_groups))
    zone_specs_prob = [("Lower-priority review", lambda p: p <= 0.1),
                       ("Deferred review", lambda p: (p > 0.1) & (p < 0.9)),
                       ("Accelerated review", lambda p: p >= 0.9)]
    for display, rule in zone_specs_prob:
        values = np.array([rule(frame["pathogenic_probability"]).mean() for _, frame in heldout_groups])
        bars = ax2.barh(y, values, left=left, color=COLORS[display], height=0.52,
                        edgecolor="white", linewidth=0.6)
        for bar, value, (_, frame), offset in zip(bars, values, heldout_groups, left):
            count = int(rule(frame["pathogenic_probability"]).sum())
            if value >= 0.075:
                ax2.text(offset + value / 2, bar.get_y() + bar.get_height() / 2,
                         f"{count}\n{value:.0%}", ha="center", va="center", fontsize=8.4,
                         color="white" if display != "Deferred review" else INK, fontweight="bold")
            elif value >= 0.04:
                ax2.text(offset + value / 2, bar.get_y() + bar.get_height() / 2,
                         f"{count}\n{value:.1%}", ha="center", va="center", fontsize=6.8,
                         color="white" if display != "Deferred review" else INK, fontweight="bold")
        left += values
    low_plp = heldout_groups[1][1]
    low_count = int((low_plp["pathogenic_probability"] <= 0.1).sum())
    low_fraction = low_count / len(low_plp)
    ax2.annotate(f"{low_count} ({low_fraction:.1%})", xy=(low_fraction / 2, 1), xytext=(0.09, 0.72),
                 fontsize=7.2, fontweight="bold", color="#245E88", ha="center",
                 arrowprops=dict(arrowstyle="-", color="#245E88", lw=0.8))
    ax2.set_yticks(y, [f"{name}\n(n={len(frame):,})" for name, frame in heldout_groups])
    ax2.invert_yaxis()
    ax2.set(xlim=(0, 1), xlabel="Fraction of held-out class")
    ticks = np.linspace(0, 1, 6)
    ax2.set_xticks(ticks, [f"{int(v * 100)}%" for v in ticks])
    ax2.set_title("B  Held-out outcomes by review zone", loc="left",
                  fontsize=12.5, fontweight="bold", color=INK)
    ax2.grid(axis="x", color="#D9DEDF", lw=0.7, alpha=0.8)
    ax2.set_axisbelow(True)

    y = np.arange(len(sources))
    left = np.zeros(len(sources))
    zone_specs = [("Lower-priority review", "lower_priority"), ("Deferred review", "deferred"), ("Accelerated review", "accelerated")]
    for display, raw in zone_specs:
        values = np.array([(frame["review_zone"].eq(raw)).mean() for _, frame in sources])
        bars = ax3.barh(y, values, left=left, color=COLORS[display], height=0.58,
                        edgecolor="white", linewidth=0.6)
        for bar, value, (_, frame), offset in zip(bars, values, sources, left):
            count = int((frame["review_zone"].eq(raw)).sum())
            if value >= 0.10:
                ax3.text(offset + value / 2, bar.get_y() + bar.get_height() / 2,
                         f"{count:,}\n{value:.1%}", ha="center", va="center", fontsize=8.2,
                         color="white" if display != "Deferred review" else INK, fontweight="bold")
        left += values
    ax3.set_yticks(y, [f"{name}\n(n={len(frame):,})" for name, frame in sources])
    ax3.invert_yaxis()
    ax3.set(xlim=(0, 1), xlabel="Fraction of source VUS")
    ticks = np.linspace(0, 1, 6)
    ax3.set_xticks(ticks, [f"{int(v * 100)}%" for v in ticks])
    ax3.set_title("C  VUS review zones: overall and overlapping source subsets", loc="left",
                  fontsize=12.5, fontweight="bold", color=INK)
    ax3.grid(axis="x", color="#D9DEDF", lw=0.7, alpha=0.8)
    ax3.set_axisbelow(True)

    for ax in (ax1, ax2, ax3):
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color("#7A878C")

    handles = [Patch(facecolor=COLORS[label], label=label) for label in COLORS]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 0.052))
    fig.text(0.5, 0.022,
             "These zones set review order only. They are not ACMG/AMP classifications and cannot independently change patient care.",
             ha="center", color="#4D5C62", fontsize=9.4)
    fig.text(0.5, 0.006, "Exploratory CardioBoost-derived thresholds; not calibrated clinical cutoffs.",
             ha="center", color="#7A4B31", fontsize=8.7)
    fig.subplots_adjust(left=0.085, right=0.985, top=0.955, bottom=0.145)

    for output in OUTPUTS:
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
        fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
