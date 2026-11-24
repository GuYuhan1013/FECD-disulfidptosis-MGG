# -*- coding: utf-8 -*-
"""S07 plot: Fig. S4 for the GSE74123 direction check.

Reads:
  TableS2_GSE74123_validation.csv
  GSE74123_19gene_matrix.csv
Writes:
  FigS4_GSE74123_validation.png / .pdf
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
TAB = pd.read_csv(HERE / "TableS2_GSE74123_validation.csv")
MAT = pd.read_csv(HERE / "GSE74123_19gene_matrix.csv", index_col=0)

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7,
})

BLUE = "#245A8D"
ORANGE = "#C4581E"
DARK = "#333333"
GREY = "#9A9A9A"

# Fixed display order: energy first, then disease-up.
order = ["OXSM", "GYS1", "LRPPRC", "RPN1", "NCKAP1", "IQGAP1", "FLNB", "SLC7A11", "FLNA"]
d = TAB.set_index("gene").loc[order]

fig = plt.figure(figsize=(7.4, 8.6), dpi=600)
gs = fig.add_gridspec(3, 2, height_ratios=[2.15, 1.45, 1.35], hspace=0.42, wspace=0.42)

# Panel A: discovery vs validation log2FC
ax = fig.add_subplot(gs[0, :])
y = np.arange(len(order))[::-1]
h = 0.34
ax.barh(y + h / 2, d["discovery_log2FC"], height=h, color=BLUE, edgecolor="none", label="Discovery (GSE142538)")
ax.barh(y - h / 2, d["validation_log2FC"], height=h, color=ORANGE, edgecolor="none", label="Validation (GSE74123)")
for yi, (g, row) in zip(y, d.iterrows()):
    p = row["welch_P"]
    if p < 0.05:
        x = row["validation_log2FC"]
        ax.text(x + (0.04 if x >= 0 else -0.04), yi - h / 2, "*", ha="left" if x >= 0 else "right",
                va="center", fontsize=8, color=DARK)
ax.axvline(0, color=DARK, lw=0.8)
ax.set_yticks(y)
ax.set_yticklabels(order, fontstyle="italic")
ax.set_xlabel("log2 fold change (FECD - control)")
ax.set_xlim(-2.25, 1.18)
ax.set_title("A  Nine candidate genes", loc="left", fontweight="bold")
ax.text(-2.22, len(order) - 0.05, "energy / anabolic", color=BLUE, fontsize=7.5, va="top")
ax.text(1.15, len(order) - 0.05, "disease-up: cytoskeletal\n+ cystine transport", color=ORANGE,
        fontsize=7.2, va="top", ha="right")
ax.legend(loc="lower left", bbox_to_anchor=(0.0, 0.0), ncol=1, frameon=False,
          handles=[Patch(color=BLUE, label="Discovery (GSE142538)"),
                   Patch(color=ORANGE, label="Validation (GSE74123)")])
ax.spines[["top", "right"]].set_visible(False)

# Panel B: four core genes, direction-oriented values
axb = fig.add_subplot(gs[1, 0])
core = ["OXSM", "GYS1", "FLNB", "IQGAP1"]
ctrl = ["GSM1910841", "GSM1910843", "GSM1910845", "GSM1910847"]
fecd = ["GSM1910840", "GSM1910842", "GSM1910844", "GSM1910846"]
positions = np.arange(len(core))
for i, g in enumerate(core):
    vals_c = MAT.loc[ctrl, g].values
    vals_f = MAT.loc[fecd, g].values
    bp = axb.boxplot([vals_c, vals_f], positions=[i - 0.16, i + 0.16], widths=0.26,
                     patch_artist=True, showfliers=False,
                     medianprops=dict(color="black", lw=1),
                     boxprops=dict(linewidth=0.8),
                     whiskerprops=dict(linewidth=0.8),
                     capprops=dict(linewidth=0.8))
    for patch, color in zip(bp["boxes"], [BLUE, "#9E2A2B"]):
        patch.set_facecolor(color)
        patch.set_alpha(0.85)
    rng = np.random.default_rng(20261124)
    axb.scatter(np.full(4, i - 0.16) + rng.uniform(-0.045, 0.045, 4), vals_c,
                s=11, facecolor="white", edgecolor=BLUE, linewidth=0.7, zorder=3)
    axb.scatter(np.full(4, i + 0.16) + rng.uniform(-0.045, 0.045, 4), vals_f,
                s=11, facecolor="white", edgecolor="#9E2A2B", linewidth=0.7, zorder=3)
    p = TAB.loc[TAB["gene"] == g, "welch_P"].iloc[0]
axb.set_xticks(positions)
axb.set_xticklabels([f"{g}\nP = {TAB.loc[TAB['gene'] == g, 'welch_P'].iloc[0]:.3f}" for g in core],
                    fontstyle="italic")
axb.set_ylim(-8.5, 12.5)
axb.set_ylabel("Direction-oriented probe score (log2)")
axb.set_title("B  Four core genes", loc="left", fontweight="bold")
axb.legend(handles=[Patch(color=BLUE, label="Control (n = 4)"),
                    Patch(color="#9E2A2B", label="FECD (n = 4)")],
           loc="lower right", frameon=False, ncol=1)
axb.spines[["top", "right"]].set_visible(False)

# Panel C: discovery vs validation for 19 genes
axc = fig.add_subplot(gs[1, 1])
colors = np.where(TAB["direction"].eq("concordant"), BLUE, ORANGE)
axc.scatter(TAB["discovery_log2FC"], TAB["validation_log2FC"], s=18, c=colors,
            edgecolor="white", linewidth=0.3, zorder=3)
axc.axhline(0, color=GREY, lw=0.8, ls="--")
axc.axvline(0, color=GREY, lw=0.8, ls="--")
for g in ["OXSM", "GYS1", "FLNB", "IQGAP1"]:
    row = TAB.loc[TAB["gene"] == g].iloc[0]
    axc.annotate(g, (row["discovery_log2FC"], row["validation_log2FC"]),
                 xytext=(4, -7), textcoords="offset points", fontsize=6.3, fontstyle="italic")
r = np.corrcoef(TAB["discovery_log2FC"], TAB["validation_log2FC"])[0, 1]
axc.text(0.98, 0.03, f"r = {r:.2f}; P = 0.014", transform=axc.transAxes,
         ha="right", va="bottom", fontsize=7)
axc.set_xlabel("Discovery log2FC")
axc.set_ylabel("Validation log2FC")
axc.set_title("C  Direction across 19 genes", loc="left", fontweight="bold")
axc.text(0.03, 0.97, "Concordant 11/19\nDiscordant 8/19", transform=axc.transAxes,
         ha="left", va="top", fontsize=7)
axc.spines[["top", "right"]].set_visible(False)

# Panel D: sign concordance
axd = fig.add_subplot(gs[2, :])
labels = ["Core four", "Nine candidates", "All 19"]
vals = [4 / 4, 6 / 9, 11 / 19]
bars = axd.bar(np.arange(3), vals, width=0.55, color=["#8E3B46", ORANGE, BLUE],
               edgecolor="white", linewidth=0.8)
for i, (b, v, lab) in enumerate(zip(bars, vals, ["4/4", "6/9", "11/19"])):
    axd.text(b.get_x() + b.get_width() / 2, v + 0.03, lab, ha="center", va="bottom", fontsize=8)
axd.axhline(0.5, color=GREY, lw=0.8, ls="--")
axd.set_xticks(np.arange(3))
axd.set_xticklabels(labels)
axd.set_ylim(0, 1.18)
axd.set_ylabel("Direction concordance")
axd.set_title("D  Sign concordance with discovery", loc="left", fontweight="bold")
axd.text(0.99, 0.92, "dashed line = chance", transform=axd.transAxes,
         ha="right", va="top", fontsize=7, color=GREY)
axd.spines[["top", "right"]].set_visible(False)

fig.savefig(HERE / "FigS4_GSE74123_validation.png", dpi=600, bbox_inches="tight")
fig.savefig(HERE / "FigS4_GSE74123_validation.pdf", bbox_inches="tight")
print("saved FigS4 in", HERE)
