from pathlib import Path
import sys
import pandas as pd, numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts_common"))
from plot_style import *  # noqa
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from scipy.stats import gaussian_kde

B = str(Path(__file__).resolve().parents[2])
A = f"{B}/output/04_article"
D = f"{A}/R2_peptide_properties/data"
OUT = f"{A}/R2_peptide_properties/figures"
props = pd.read_csv(f"{D}/peptide_properties_matched.tsv", sep="\t")
cat = pd.read_csv(f"{D}/fig2B_source_categories_matched.tsv", sep="\t")

SHORT = {"exon_diff(CDS out-of-frame)": "CDS-OF", "intergenic": "Intergenic",
         "intron": "Intron", "junction": "Junction",
         "3'UTR": "3'UTR", "5'UTR": "5'UTR"}
props["cat_short"] = props["source_cat"].map(SHORT)
order = cat.sort_values("n", ascending=False)["source_cat"].tolist()
cat_order = [SHORT[c] for c in order]
CCOL = dict(zip(cat_order, CLASS8))
for p in ["length", "molecular_weight", "isoelectric_point", "gravy"]:
    props[p] = pd.to_numeric(props[p], errors="coerce")

fig = plt.figure(figsize=(11.0, 7.0))
gs = fig.add_gridspec(2, 3, hspace=0.52, wspace=0.34,
                      left=0.05, right=0.97, top=0.95, bottom=0.08)


ax = fig.add_subplot(gs[0, 0])
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
def node(ax, x, y, w, h, big, unit, edge, big_fs=11):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012",
                 fc="white", ec=edge, lw=1.1))
    ax.add_patch(plt.Rectangle((x, y), 0.015, h, fc=edge, ec="none"))
    ax.text(x + 0.05, y + h*0.66, big, fontsize=big_fs, fontweight="bold", va="center")
    ax.text(x + 0.05, y + h*0.25, unit, fontsize=6.0, color="#555555", va="center")
# Level 1
node(ax, 0.04, 0.86, 0.60, 0.115, f"{103284:,}", "NCP mapping records (genomic position)", "#9DC3E6")
ax.annotate("", xy=(0.34, 0.72), xytext=(0.34, 0.855),
            arrowprops=dict(arrowstyle="-|>", lw=1.0, color="grey"))
# Level 2
node(ax, 0.04, 0.585, 0.60, 0.135, f"{7619:,}", "transcript-supported NCPs (unique sequences)", "#F4B183")
ax.text(0.665, 0.72, f"{29660:,} pairs", fontsize=6.6, color="#555555", va="center")
ax.text(0.665, 0.655, f"{16595:,} transcripts", fontsize=6.6, color="#555555", va="center")
ax.text(0.665, 0.59, f"{4466:,} loci", fontsize=6.6, color="#555555", va="center")
ax.plot([0.64, 0.655], [0.652, 0.652], color="#B0B0B0", lw=0.6)

ax.annotate("", xy=(0.26, 0.445), xytext=(0.26, 0.58),
            arrowprops=dict(arrowstyle="-|>", lw=1.0, color="#A9D18E"))
ax.annotate("", xy=(0.775, 0.445), xytext=(0.40, 0.58),
            arrowprops=dict(arrowstyle="-|>", lw=1.0, color="#B0B0B0",
                            connectionstyle="arc3,rad=-0.25"))
# Level 3 (single)
node(ax, 0.04, 0.31, 0.44, 0.135, f"{2519:,}",
     "single-transcript NCPs\n(1,741 transcripts | 1,715 loci)", "#A9D18E", big_fs=10)
# multi branch node
ax.add_patch(FancyBboxPatch((0.58, 0.31), 0.40, 0.135, boxstyle="round,pad=0.012",
             fc="#F5F5F5", ec="#B0B0B0", lw=0.9))
ax.text(0.78, 0.392, f"{7619-2519:,}", fontsize=10, fontweight="bold", ha="center", va="center")
ax.text(0.78, 0.322, "multi-transcript NCPs", fontsize=5.8, color="#555555", ha="center", va="center")

ax.text(0.04, 0.215, "Units differ across levels: mapping records vs peptide sequences;\n"
        "transcript/loci counts describe the matched relations, not NCP counts.",
        fontsize=6.0, color="#777777", va="top")
panel(ax, "A", dx=-0.02)


ax = fig.add_subplot(gs[0, 1])
cs = cat.sort_values("n")
ax.barh(range(len(cs)), cs["n"], color=[CCOL[SHORT[c]] for c in cs["source_cat"]])
for y, (v, p) in enumerate(zip(cs["n"], cs["pct"])):
    ax.text(v + 60, y, f"{v:,} ({p}%)", va="center", fontsize=6.8)
ax.set_yticks(range(len(cs)), [SHORT[c] for c in cs["source_cat"]], fontsize=7.5)
ax.set_xlabel("Number of NCPs")
ax.set_xlim(0, cs["n"].max() * 1.30)
panel(ax, "B")


ax = fig.add_subplot(gs[0, 2])
data = [props.loc[props["cat_short"] == c, "length"].dropna().values for c in cat_order]
vp = ax.violinplot(data, positions=range(len(cat_order)), widths=0.82,
                   showmedians=False, showextrema=False)
for body, c in zip(vp["bodies"], cat_order):
    body.set_facecolor(CCOL[c]); body.set_edgecolor("black"); body.set_linewidth(0.4); body.set_alpha(0.92)
for i, v in enumerate(data):
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    ax.plot([i, i], [q1, q3], color="black", lw=1.0)
    ax.plot([i - 0.14, i + 0.14], [med, med], color="white", lw=2.4, solid_capstyle="butt")
ax.set_xticks(range(len(cat_order)), cat_order, fontsize=6.3, rotation=30)
for t in ax.get_xticklabels(): t.set_ha("right")
ax.set_ylabel("Peptide length (aa)")
ax.text(0.99, 1.03, "Transcript-supported NCPs | n = 7,619", transform=ax.transAxes,
        ha="right", fontsize=6.5, color="#555555")
panel(ax, "C")


ax = fig.add_subplot(gs[1, 0])
xs = np.linspace(500, 4200, 320)
for i, c in enumerate(cat_order):
    v = props.loc[props["cat_short"] == c, "molecular_weight"].dropna().values
    k = gaussian_kde(v)
    base = len(cat_order) - 1 - i
    dens = k(xs) / k(xs).max() * 0.88
    ax.plot(xs, base + dens, color=CCOL[c], lw=1.1)
    ax.fill_between(xs, base, base + dens, color=CCOL[c], alpha=0.85, lw=0)
    ax.axhline(base, color="white", lw=0.6)
med_all = props["molecular_weight"].median()
ax.axvline(med_all, color="grey", ls="--", lw=0.8)
ax.set_yticks([len(cat_order)-1-i for i in range(len(cat_order))], cat_order, fontsize=7)
ax.set_xlabel("Molecular weight (Da)")
ax.set_xlim(500, 4200); ax.set_ylim(-0.3, len(cat_order) - 0.1)
panel(ax, "D")


ax = fig.add_subplot(gs[1, 1])
data = [props.loc[props["cat_short"] == c, "isoelectric_point"].dropna().values for c in cat_order]
bp = ax.boxplot(data, tick_labels=cat_order, patch_artist=True, showfliers=False, widths=0.62,
                medianprops=dict(color="black", lw=1.0),
                whiskerprops=dict(lw=0.8), capprops=dict(lw=0.8))
for body, c in zip(bp["boxes"], cat_order):
    body.set_facecolor(CCOL[c]); body.set_edgecolor("black"); body.set_linewidth(0.6); body.set_alpha(0.75)
ax.set_xticklabels(cat_order, fontsize=6.3, rotation=30)
for t in ax.get_xticklabels(): t.set_ha("right")
ax.set_ylabel("Theoretical pI")
panel(ax, "E")


ax = fig.add_subplot(gs[1, 2])
xs = np.linspace(-3, 3.5, 300)
for c in cat_order:
    v = props.loc[props["cat_short"] == c, "gravy"].dropna().values
    k = gaussian_kde(v)
    ax.plot(xs, k(xs), color=CCOL[c], lw=1.2, label=c)
    ax.fill_between(xs, k(xs), color=CCOL[c], alpha=0.10)
ax.axvline(0, color="grey", ls="--", lw=0.8)
ax.set_xlabel("GRAVY"); ax.set_ylabel("Density")
ax.legend(frameon=False, fontsize=5.8, ncol=2, loc="upper right")
panel(ax, "F")

fig.savefig(f"{OUT}/Fig2.pdf")
fig.savefig(f"{OUT}/Fig2.svg")
plt.close(fig)
print("Figure 2 rendered (A branching | B-F properties)")
