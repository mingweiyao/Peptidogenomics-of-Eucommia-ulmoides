from pathlib import Path
import sys, os
import pandas as pd, numpy as np
from scipy.stats import gaussian_kde as stats_gaussian_kde
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts_common"))
from plot_style import *  # noqa

B = str(Path(__file__).resolve().parents[2])
A = f"{B}/output/04_article"
OUT = f"{A}/R1_identity/figures"
MAT = ["Root", "Xylem", "Seed", "Leaf"]
TC = [TISSUE[m + "_MS"] if m + "_MS" in TISSUE else TISSUE[m] for m in MAT]

CP = pd.read_excel(f"{B}/output/01location/CP_peptides.xlsx", sheet_name="CP",
                   usecols=["material", "length", "molecular_weight"])
NCP = pd.read_excel(f"{B}/output/01location/NCP_peptides.xlsx", sheet_name="NCP",
                    usecols=["material", "length", "molecular_weight"])
cnts = pd.read_csv(f"{A}/R1_identity/fig1A_counts.tsv", sep="\t")
wl = pd.read_csv(f"{A}/R1_identity/fig1B_length_wilcoxon.tsv", sep="\t").set_index("tissue")

fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0))
ax = axes[0, 0]
x = np.arange(4); wdt = 0.38
ax.bar(x - wdt/2, cnts["CP_rows"][:4], wdt, color=CP_C, label="CP")
ax.bar(x + wdt/2, cnts["NCP_rows"][:4], wdt, color=NCP_C, label="NCP")
for i, (a, b) in enumerate(zip(cnts["CP_rows"][:4], cnts["NCP_rows"][:4])):
    ax.text(i - wdt/2, a, f"{a:,}", ha="center", va="bottom", fontsize=6.5)
    ax.text(i + wdt/2, b, f"{b:,}", ha="center", va="bottom", fontsize=6.5)
ax.set_xticks(x, MAT); ax.set_ylabel("Identified peptides")
ax.legend(frameon=False); panel(ax, "A")

ax = axes[0, 1]
data, pos_, cols = [], [], []
for i, m in enumerate(MAT):
    data.append(CP.loc[CP["material"] == m, "length"].astype(float)); pos_.append(i - 0.2); cols.append(CP_C)
    data.append(NCP.loc[NCP["material"] == m, "length"].astype(float)); pos_.append(i + 0.2); cols.append(NCP_C)
bp = ax.boxplot(data, positions=pos_, widths=0.34, patch_artist=True, showfliers=False,
                medianprops=dict(color="black", lw=0.8))
for p, c in zip(bp["boxes"], cols):
    p.set_facecolor(c); p.set_edgecolor("black"); p.set_linewidth(0.6)
ax.set_xticks(range(4), MAT); ax.set_ylabel("Peptide length (aa)")
for i, m in enumerate(MAT):
    significance_mark(ax, i - 0.2, i + 0.2, 64.5, 1.2, wl.loc[m, "p_value"])
ax.set_ylim(0, 70); panel(ax, "B")

for ax_i, src, lab, letter in [(axes[1, 0], CP, "CP", "C"), (axes[1, 1], NCP, "NCP", "D")]:
    for m in MAT:
        v = src.loc[src["material"] == m, "molecular_weight"].astype(float)
        if len(v) < 10: continue
        kde = stats_gaussian_kde(v)
        xs = np.linspace(0, 6000, 300)
        ax_i.plot(xs, kde(xs), color=TISSUE.get(m + "_MS", TISSUE.get(m)), lw=1.2, label=m)
    mean_all = src["molecular_weight"].astype(float).mean()
    ax_i.axvline(mean_all, color="grey", ls="--", lw=0.8)
    ax_i.text(mean_all + 60, ax_i.get_ylim()[1]*0.9 if ax_i.get_ylim()[1] > 0 else 0,
              f"mean {mean_all:,.0f} Da", fontsize=7, color="grey", va="top")
    ax_i.set_xlabel("Molecular weight (Da)"); ax_i.set_ylabel("Density")
    ax_i.set_xlim(0, 6000); panel(ax_i, letter)
axes[1, 0].legend(frameon=False, fontsize=7)
fig.tight_layout()
save(fig, OUT, "Fig1")
print("Fig1 done")
