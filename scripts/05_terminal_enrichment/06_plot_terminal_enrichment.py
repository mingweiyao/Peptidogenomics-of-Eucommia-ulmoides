# -*- coding: utf-8 -*-
"""Supplementary Fig. 8: terminal enrichment of NCP candidates vs permutation nulls.
(A) baseline chromosome-matched permutation; (B) expression-restricted null model."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os as _os
import argparse

_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--out", default=r"./output", help="output root directory")
args = _ap.parse_args()
OUTROOT = os.path.abspath(os.path.expanduser(args.out))

BASE = OUTROOT + r"/SupFig5_terminal_enrichment"
OUTD = OUTROOT + r"/SupFig5_terminal_enrichment"
os.makedirs(OUTD, exist_ok=True)

base_df = pd.read_csv(os.path.join(BASE, "Fig3C_terminal_enrichment_permutation.csv"))
expr_df = pd.read_csv(os.path.join(BASE, "Fig3C_terminal_enrichment_ExprRestricted.csv"))

chrs = base_df["Chr"].astype(int).tolist()
obs_b = base_df["Observed_in_terminal10pct"].astype(float).values
null_b = base_df["Null_mean"].astype(float).values
p_b = base_df["Empirical_p_gt"].astype(float).values

obs_e = expr_df["Observed"].astype(float).values
null_e = expr_df["Null_mean_ExprRestricted"].astype(float).values
p_e = expr_df["Empirical_p_ExprRestricted"].astype(float).values

sig_e = p_e < 0.05
print("ExprRestricted significant chromosomes:", [chrs[i] for i in range(17) if sig_e[i]])
print("baseline observed range:", obs_b.min(), "-", obs_b.max(), "| null range:", null_b.min(), "-", null_b.max())
print("baseline all P<=0.005:", bool((p_b <= 0.005).all()))
print("genome-wide baseline: observed", obs_b.sum(), "null", round(null_b.sum(), 1))
print("genome-wide expr: observed", obs_e.sum(), "null", round(null_e.sum(), 1))

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 11,
    "axes.linewidth": 1.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.1), dpi=300)

panels = [
    (axes[0], obs_b, null_b, None,
     "A", "Baseline permutation model"),
    (axes[1], obs_e, null_e, sig_e,
     "B", "Expression-restricted null model"),
]
for ax, obs, nl, sig, letter, title in panels:
    x = np.arange(1, 18)
    ax.plot(x, nl, "--", color="#4d4d4d", lw=1.4, label="Permutation null (mean)", zorder=2)
    if sig is None:
        ax.scatter(x, obs, s=46, color="#128f8f", edgecolor="white", lw=0.7,
                   label="Observed", zorder=3)
        handles = [
            Line2D([], [], marker="o", ls="", mfc="#128f8f", mec="white", ms=8, label="Observed"),
            Line2D([], [], color="#4d4d4d", ls="--", lw=1.4, label="Permutation null (mean)"),
        ]
    else:
        for xi, yi, s in zip(x, obs, sig):
            ax.scatter(xi, yi, s=52,
                       color="#E8850C" if s else "#9a9a9a",
                       edgecolor="white", lw=0.7, zorder=3)
        handles = [
            Line2D([], [], marker="o", ls="", mfc="#E8850C", mec="white", ms=8,
                   label="Significant (empirical P < 0.05)"),
            Line2D([], [], marker="o", ls="", mfc="#9a9a9a", mec="white", ms=8,
                   label="Not significant"),
            Line2D([], [], color="#4d4d4d", ls="--", lw=1.4, label="Permutation null (mean)"),
        ]
    ax.set_ylim(top=max(obs.max(), nl.max()) + 38)
    ax.set_xticks(x)
    ax.set_xticklabels([str(c) for c in x])
    ax.set_xlabel("Chromosome")
    ax.set_xlim(0.3, 17.7)
    ax.set_title(title, fontsize=12, pad=26, loc="left")
    ax.text(-0.075, 1.13, letter, transform=ax.transAxes,
            fontsize=17, fontweight="bold", va="top", ha="left")
    ax.legend(handles=handles, frameon=False, fontsize=9, loc="upper left",
              bbox_to_anchor=(0.01, 1.02), handletextpad=0.2, borderaxespad=0.0)

axes[0].set_ylabel("NCP candidates in terminal 10%")

fig.tight_layout(rect=(0, 0, 1, 0.97))
png = os.path.join(OUTD, "SuppFig5_terminal_enrichment.png")
pdf = os.path.join(OUTD, "SuppFig5_terminal_enrichment.pdf")
fig.savefig(png, dpi=300, facecolor="white")
fig.savefig(pdf, facecolor="white")
print("saved:", png)
print("saved:", pdf)
